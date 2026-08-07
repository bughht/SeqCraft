"""
The MRzero bridge's gradient-area arithmetic, which decides where a sample lands in k-space.

``examples/lib`` is not part of the package and its simulator needs MRzeroCore and torch -- but the
part that was wrong is pure numpy, and it was wrong in a way that produced a Nyquist ghost on EPI
while leaving a spiral merely a little blurred.  So it is tested here, against
:func:`seqcraft.core.events.pwl_moment`, which is an independent exact integrator.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.core import events as ev

BRIDGE = Path(__file__).resolve().parents[2] / 'examples' / 'lib' / 'mr0_bridge.py'


@pytest.fixture(scope='module')
def bridge() -> object:
    """Load ``examples/lib/mr0_bridge.py`` by path; it is not an installed module."""
    if not BRIDGE.exists():  # pragma: no cover - only if the examples are stripped
        pytest.skip('examples/lib/mr0_bridge.py is not present')
    spec = importlib.util.spec_from_file_location('_mr0_bridge_under_test', BRIDGE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def system() -> sc.System:
    return sc.System.preset('cima_x').derate('epi', grad=0.85, slew=1.0)


def _block(*events: object) -> object:
    """A stand-in for a pypulseq block: the bridge only reads ``gx``/``gy``/``gz`` off it."""
    fields = {f'g{axis}': None for axis in 'xyz'}
    for event in events:
        fields[f'g{event.channel}'] = event
    return SimpleNamespace(**fields)


def test_cumulative_area_is_exact_between_knots(bridge, system) -> None:
    """
    Against ``pwl_moment``, at query times deliberately **off** the raster.

    On the raster any reasonable scheme agrees; between raster points a cumulative integral is a
    quadratic, and interpolating it linearly -- which the previous version did -- under-reads on a
    ramp.  An ADC's sample centres fall exactly there.
    """
    opts = system.limits('epi')
    grad = pp.make_extended_trapezoid(
        'x', amplitudes=np.array([0.0, 2.0787e6, 0.0, -2.0787e6]),
        times=np.array([0.0, 260e-6, 520e-6, 780e-6]), system=opts)
    times, amps = ev.knots_of(grad)

    query = np.array([0.0, 1e-6, 43.4e-6, 131e-6, 259.9e-6, 260e-6, 391e-6, 519.5e-6, 780e-6])
    got = bridge._cumulative_area(times, amps, query)
    for index, when in enumerate(query):
        piece = ev.knots_of(grad)
        keep = piece[0] <= when
        edge_t = np.concatenate([piece[0][keep], [when]])
        edge_a = np.concatenate([piece[1][keep], [np.interp(when, piece[0], piece[1])]])
        wanted = ev.pwl_moment(edge_t, edge_a, 0)
        assert got[index] == pytest.approx(wanted, rel=1e-12, abs=1e-15), f'at {when * 1e6:.1f} us'


def test_areas_preserve_the_total_however_the_edges_fall(bridge, system) -> None:
    """
    The property the whole function exists for, and the one a trapezoidal rule loses.

    A trapezoid starts and ends at zero, so sampling it *at* interval edges and applying the
    trapezoidal rule integrates one interval spanning a whole prephaser to exactly zero -- k-space
    then starts at the origin instead of its corner.
    """
    opts = system.limits('epi')
    prephaser = pp.make_trapezoid('x', area=-270.0, system=opts)
    duration = float(pp.calc_duration(prephaser))
    block = _block(prephaser)

    for edges in (
        np.array([0.0, duration]),                                # one interval: the hard case
        np.linspace(0.0, duration, 4),
        np.array([0.0, 13e-6, 197e-6, duration]),                 # off-raster, uneven
    ):
        areas = bridge._areas(block, duration, system.grad_raster.dt, edges)
        assert float(areas[:, 0].sum()) == pytest.approx(-270.0, rel=1e-9)


def test_areas_at_sample_centres_place_k_where_pulseq_does(bridge, system) -> None:
    """
    An EPI echo, summed at ADC sample centres, against ``calculate_kspacePP``.

    This is the check that catches the two faults the bridge had: a half-dwell offset in where a
    sample sits, and the lead and tail areas -- which on an EPI train are precisely where the
    phase-encode blips live -- being left out of the running total.  Both are invisible on a spiral,
    whose lead is 10 us of almost nothing and whose readout never reverses; on EPI the gradient
    reverses every echo, so a systematic offset becomes an **alternating** one, which is a Nyquist
    ghost rather than a blur.
    """
    ro = sc.modules.EPIReadout(
        system, fov_ro_mm=240, matrix_ro=64, fov_pe_mm=240, matrix_pe=64,
        partial_fourier_pe=0.75, regime='epi')
    out = sc.compile(sc.LogicBlock('epi').add(0.0, ro.build()), system, regime='epi')

    running = np.zeros(3)
    placed: list[np.ndarray] = []
    for index in sorted(out.seq.block_events):
        block = out.seq.get_block(index)
        duration = float(out.seq.block_durations[index])
        adc = getattr(block, 'adc', None)
        if adc is None:
            areas = bridge._areas(block, duration, system.grad_raster.dt,
                                  np.array([0.0, duration]))
            running = running + areas[0]
            continue
        n, dwell, delay = int(adc.num_samples), float(adc.dwell), float(adc.delay)
        centres = delay + (np.arange(n) + 0.5) * dwell
        edges = np.concatenate(([0.0], centres, [duration]))
        areas = bridge._areas(block, duration, system.grad_raster.dt, edges)
        for sample in range(n):
            running = running + areas[sample]
            placed.append(running.copy())
        running = running + areas[n]

    mine = np.stack(placed, axis=1)[:2]
    theirs = out.kspace()['k_adc'][:2, : mine.shape[1]]
    worst = float(np.abs(mine - theirs).max()) / ro.dk_pe_per_m
    assert worst < 1e-3, f'sample placement disagrees with pypulseq by {worst:.4f} dk'

    # And specifically: no *alternating* component, which is the one that makes stripes.
    drift = (mine[0] - theirs[0]).reshape(ro.n_echoes, ro.samples_per_echo).mean(axis=1)
    alternating = abs(float(drift[::2].mean() - drift[1::2].mean())) / ro.dk_pe_per_m
    assert alternating < 1e-3, f'alternating kx offset of {alternating:.4f} dk -- a Nyquist ghost'
