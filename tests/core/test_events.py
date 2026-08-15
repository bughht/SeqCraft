"""
Limit checking on gradients whose knots are neither uniform nor on the raster.

Every case here was silent before ``check_limits`` measured on the union of knots, and every one
is reachable only through an *extended trapezoid* -- which is what an EPI readout train has to be
if it is to be one event per axis rather than one per echo.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.core import events as ev


@pytest.fixture
def opts() -> pp.Opts:
    """A Cima.X derated for EPI: 200 mT/m at 85 %, full slew, which is where a train wants to run."""
    return sc.opts.derate(
        pp.Opts(
            max_grad=200, grad_unit='mT/m', max_slew=200, slew_unit='T/m/s', B0=2.8936,
            rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
            adc_samples_limit=8192,
        ),
        grad=0.85,
    )


def _epi_junction(opts: object) -> tuple[object, object]:
    """
    One echo junction of the measured 128 x 128 EPI train: gx through zero, a gy blip across it.

    Amplitudes and knots are the probe's own: a 48.83 mT/m lobe on a 260 us ramp, and a 3.26 mT/m
    blip centred on the junction at 520 us.
    """
    gx = pp.make_extended_trapezoid(
        'x',
        amplitudes=np.array([0.0, 2.0787e6, 0.0, -2.0787e6 * 10.0 / 260.0]),
        times=np.array([0.0, 260e-6, 520e-6, 530e-6]),
        system=opts,
    )
    gy = pp.make_extended_trapezoid(
        'y',
        amplitudes=np.array([0.0, 0.0, 1.3889e5, 1.3889e5 * 20.0 / 30.0]),
        times=np.array([0.0, 490e-6, 520e-6, 530e-6]),
        system=opts,
    )
    return gx, gy


def _exact_norm_slew(events: tuple[object, ...]) -> float:
    """The vector-norm slew, computed independently on the union of the events' knots."""
    pieces = [ev.knots_of(e) for e in events]
    grid = np.unique(np.concatenate([t for t, _ in pieces]))
    stack = np.stack([np.interp(grid, t, a, left=0.0, right=0.0) for t, a in pieces])
    norm = np.linalg.norm(stack, axis=0)
    return float(np.max(np.abs(np.diff(norm) / np.diff(grid))))


def test_norm_slew_uses_knot_spacing_not_the_raster(opts) -> None:
    """
    Dividing ``diff(norm)`` by the raster overstated an EPI junction's slew by 26x.

    The knots are 260 us apart and the raster is 10 us, so the reported figure was the peak
    amplitude divided by one raster -- 4882.9 T/m/s, or 2441 % of a 200 T/m/s limit, where the
    truth is 187.8 T/m/s.  It was a warning rather than an error, so it survived: a compile of a
    96-echo train produced 97 of them and reported ``ok=True``.
    """
    gx, gy = _epi_junction(opts)
    truth = _exact_norm_slew((gx, gy))
    assert truth == pytest.approx(
        sc.convert(187.8, 'T/m/s', 'Hz/m/s', gamma=opts.gamma), rel=1e-3
    )

    reported = {kind: got for kind, _, got, _ in ev.check_limits([gx, gy], opts)}
    assert 'slew_norm' not in reported, (
        f'a legal junction at {truth / opts.max_slew * 100:.0f} % of the limit was reported as a '
        f'violation at {reported.get("slew_norm", 0.0) / opts.max_slew * 100:.0f} %'
    )
    # And the number itself, when a violation is genuine, is the knot-spacing one.
    tight = sc.opts.derate(opts, slew=0.9 * truth / float(opts.max_slew))
    got = {kind: value for kind, _, value, _ in ev.check_limits([gx, gy], tight)}
    assert got['slew_norm'] == pytest.approx(truth, rel=1e-9)


def test_norm_aligns_axes_in_time_not_by_sample_index(opts) -> None:
    """
    Stacking by index combined gx at 260 us with gy at 490 us.

    Demonstrated here in the direction that **hides** a fault, which is the dangerous one.  Both
    lobes peak at 700 us, so the vector norm genuinely reaches ``sqrt(2)`` times a legal per-axis
    amplitude -- 106 % of the limit.  But gy is described with twice as many knots, so its peak
    sits at index 2 while gx's sits at index 1; aligning by index pairs each peak with the other
    axis's half-height point and reports 84 %, comfortably inside.
    """
    amplitude = 0.75 * float(opts.max_grad)          # 182 T/m/s over a 700 us ramp: legal
    coarse = pp.make_extended_trapezoid(
        'x', amplitudes=np.array([0.0, amplitude, 0.0]),
        times=np.array([0.0, 700e-6, 1400e-6]), system=opts)
    fine = pp.make_extended_trapezoid(
        'y', amplitudes=np.array([0.0, 0.5 * amplitude, amplitude, 0.5 * amplitude, 0.0]),
        times=np.array([0.0, 350e-6, 700e-6, 1050e-6, 1400e-6]), system=opts)

    assert not [k for k, *_ in ev.check_limits([coarse], opts)]
    assert not [k for k, *_ in ev.check_limits([fine], opts)]

    peaks = {where: got for kind, where, got, _ in ev.check_limits([coarse, fine], opts)
             if kind == 'grad_norm'}
    assert 'norm' in peaks, (
        'two axes peaking at the same instant were not combined; aligning by index would give '
        f'{np.hypot(amplitude, 0.5 * amplitude) / float(opts.max_grad) * 100:.0f} % of the limit '
        f'instead of {amplitude * np.sqrt(2) / float(opts.max_grad) * 100:.0f} %'
    )
    assert peaks['norm'] == pytest.approx(amplitude * np.sqrt(2.0), rel=1e-9)


def test_simultaneous_axes_still_reach_root_two(opts) -> None:
    """The genuine case must still be caught: two axes ramping together, in vector norm."""
    strong = {'amplitude': 0.9 * float(opts.max_grad), 'duration': 1e-3, 'system': opts}
    kinds = [kind for kind, *_ in ev.check_limits(
        [pp.make_trapezoid('x', **strong), pp.make_trapezoid('y', **strong)], opts)]
    assert 'grad_norm' in kinds


def test_two_events_on_one_axis_are_both_seen(opts) -> None:
    """
    A second event on an axis used to overwrite the first in the norm's dictionary.

    Given their node times they superpose; given none they are taken to play at their own delay,
    which is what makes ``starts`` necessary rather than decorative.
    """
    half = pp.make_trapezoid('z', amplitude=0.6 * float(opts.max_grad), duration=4e-4, system=opts)
    other = pp.make_trapezoid('y', amplitude=0.1 * float(opts.max_grad), duration=4e-4, system=opts)

    # Placed back to back, the z axis peaks at 0.6 -- legal -- and the norm sees each in its turn.
    apart = ev.check_limits([half, half, other], opts, starts=[0.0, 4e-4, 0.0])
    assert not [k for k, *_ in apart if k == 'grad']

    # Placed on top of each other, z reaches 1.2 of the limit and must be an error.
    together = ev.check_limits([half, half, other], opts, starts=[0.0, 0.0, 0.0])
    peaks = {where: got for kind, where, got, _ in together if kind == 'grad'}
    assert 'z' in peaks
    assert peaks['z'] == pytest.approx(1.2 * float(opts.max_grad), rel=1e-9)


def test_waveform_of_returns_a_grad_events_own_times(opts) -> None:
    """
    ``waveform_of`` samples a trapezoid onto the raster and leaves a ``grad`` event alone.

    Stated as a test because the docstring used to promise a uniform raster for both, and
    ``check_limits`` believed it.
    """
    raster = sc.Raster(opts.grad_raster_time).dt
    t_trap, _ = ev.waveform_of(pp.make_trapezoid('x', area=100.0, system=opts), raster)
    assert np.allclose(np.diff(t_trap), raster)

    gx, _ = _epi_junction(opts)
    t_grad, _ = ev.waveform_of(gx, raster)
    assert [round(v * 1e6) for v in t_grad] == [0, 260, 520, 530]
    assert not np.allclose(np.diff(t_grad), raster)


def test_an_epi_train_compiles_without_limit_warnings(opts) -> None:
    """
    The whole point: a train built as one extended trapezoid per axis reports nothing.

    Before the fix this produced one bogus ``slew_norm_limit`` warning per echo -- 97 of them for
    96 echoes -- which is the state in which a warning stops meaning anything.
    """
    n_echo, spacing, ramp = 8, 520e-6, 260e-6
    times_x: list[float] = [0.0]
    amps_x: list[float] = [0.0]
    for echo in range(n_echo):
        sign = 1.0 if echo % 2 == 0 else -1.0
        times_x += [echo * spacing + ramp, (echo + 1) * spacing]
        amps_x += [sign * 2.0787e6, 0.0]
    times_y: list[float] = [0.0]
    amps_y: list[float] = [0.0]
    for junction in range(1, n_echo):
        start = junction * spacing - 30e-6
        times_y += [start, start + 30e-6, start + 60e-6]
        amps_y += [0.0, 1.3889e5, 0.0]
    times_y.append(n_echo * spacing)
    amps_y.append(0.0)

    train = sc.LogicBlock('epi_train')
    train.add(0.0, pp.make_extended_trapezoid(
        'x', amplitudes=np.asarray(amps_x), times=np.asarray(times_x), system=opts))
    train.add(0.0, pp.make_extended_trapezoid(
        'y', amplitudes=np.asarray(amps_y), times=np.asarray(times_y), system=opts))
    adc = pp.make_adc(num_samples=228, dwell=2000e-9, delay=30e-6, system=opts)
    for echo in range(n_echo):
        train.add(echo * spacing, adc)

    out = sc.compile(train, opts, name='epi_train')
    limits = [i for i in out.check().issues if i.kind.endswith('_limit')]
    assert not limits, f'{len(limits)} limit issues on a legal train: {limits[:2]}'
    assert not [i for i in out.check().issues if i.kind == 'grad_merge']
