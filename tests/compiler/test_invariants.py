"""
The invariants that run on every compile -- and proof that each one can fail.

`verify_against_tree` is the last line of defence: it runs for every user, on every sequence,
including the ones no test thought of.  So each check gets two tests: it passes on a good compile,
and it *fires* on a deliberately broken one.  A check that cannot fail is worse than no check,
because it reads like coverage.  Every broken case here is checked by the same function, on the
same inputs, that ``sc.compile`` used (``recheck.py``).
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pypulseq as pp
import pytest
from recheck import compile_capturing, recheck

import seqcraft as sc
from seqcraft.design.events import knots_of as _grad_knots
from seqcraft.design.events import pwl_moment, trapz


def pwl_moment_m1(times, amps):
    """m1 of a PWL waveform: the shape the compiler-private helper had before it moved."""
    return pwl_moment(times, amps, 1)


# ------------------------------------------------------------------- the exact PWL first moment
def test_pwl_m1_matches_the_closed_form_for_a_triangle() -> None:
    """
    A symmetric triangle has centroid at its apex, so ``m1 == m0 * h == A * h^2``.

    Worth pinning independently: the whole point of the closed form is that the trapezoidal rule
    is *wrong* here, since ``g(t) * t`` is quadratic between knots.
    """
    amp, h = 3.0, 2.0
    times = np.array([0.0, h, 2 * h])
    amps = np.array([0.0, amp, 0.0])
    assert pwl_moment_m1(times, amps) == pytest.approx(amp * h**2, rel=1e-14)


def test_pwl_m1_is_exact_where_trapz_is_not() -> None:
    """
    A bare ramp, where trapz is genuinely wrong and refining it converges to the closed form.

    A *symmetric* triangle is the wrong test: its two halves have equal and opposite second
    derivatives, so the trapezoidal errors cancel exactly and trapz looks perfect.  A single ramp
    has no such cancellation, which is what exposes the difference.
    """
    amp, h = 3.0, 2.0
    exact = pwl_moment_m1(np.array([0.0, h]), np.array([0.0, amp]))
    assert exact == pytest.approx(amp * h**2 / 3.0, rel=1e-14), 'closed form for a ramp'

    errors = []
    for n in (10, 100, 1000):
        t = np.linspace(0.0, h, n + 1)
        g = np.interp(t, [0.0, h], [0.0, amp])
        errors.append(abs(float(trapz(g * t, t)) - exact))
    assert errors[0] > 0.0, 'trapz must actually be wrong here, or the test proves nothing'
    assert errors[0] > errors[1] > errors[2], 'and refining it must approach the closed form'
    assert errors[-1] < 1e-4 * exact


def test_pwl_m1_shifts_by_area_times_offset() -> None:
    """The property the invariant relies on: displacing a waveform changes m1 by ``m0 * dt``."""
    times = np.array([0.0, 100e-6, 300e-6, 400e-6])
    amps = np.array([0.0, 1e4, 1e4, 0.0])
    m0 = float(trapz(amps, times))
    dt = 10e-6
    before = pwl_moment_m1(times, amps)
    after = pwl_moment_m1(times + dt, amps)
    assert after - before == pytest.approx(m0 * dt, rel=1e-12)


def test_grad_knots_keeps_the_edge_values_of_an_arbitrary_gradient(opts) -> None:
    """
    ``first`` and ``last`` are real knots half a raster outside the samples.

    Dropping them truncates the waveform at both ends, which is the mistake that made the
    original merge lose peak amplitude.
    """
    wave = np.concatenate([
        np.linspace(0, 1, 21), np.linspace(1, 0, 21)[1:],
    ]) * 0.3 * opts.max_grad
    g = pp.make_arbitrary_grad('x', waveform=wave, first=0.0, last=0.0, system=opts)
    times, amps = _grad_knots(g, 0.0)
    assert times[1] == pytest.approx(0.5 * float(opts.grad_raster_time)), (
        'the first sample sits at a raster centre, which is why first/last are extra knots'
    )
    assert len(times) == len(wave) + 2, 'first and last must appear as knots'
    assert times[0] == pytest.approx(0.0)
    assert times[-1] == pytest.approx(float(g.shape_dur))
    assert amps[0] == pytest.approx(0.0), 'the edge knot carries `first`'
    assert amps[-1] == pytest.approx(0.0), 'the edge knot carries `last`'


def test_grad_knots_does_not_duplicate_an_extended_trapezoids_edges(opts) -> None:
    """Its ``tt`` already reaches both edges, so first/last coincide with existing knots."""
    g = pp.make_extended_trapezoid(
        'x', times=np.array([0.0, 100e-6, 300e-6, 400e-6]),
        amplitudes=np.array([0.0, 1e4, 1e4, 0.0]), system=opts)
    times, _ = _grad_knots(g, 0.0)
    assert len(times) == 4, f'expected 4 knots, got {len(times)}'


# ------------------------------------------------------------------------ the checks can fail
def _late(placed, axis: str, dt: float) -> list:
    """The placed events with every gradient on `axis` moved `dt` later: a mistimed expectation."""
    return [
        dataclasses.replace(p, node_t=p.node_t + dt, start=p.start + dt, end=p.end + dt,
                            res_start=p.res_start + dt, res_end=p.res_end + dt)
        if p.kind in ('trap', 'grad') and p.event.channel == axis else p
        for p in placed
    ]


def test_m1_catches_a_gradient_that_plays_at_the_wrong_time(opts) -> None:
    """
    The failure m0 is blind to.

    Area is exactly the quantity a time shift preserves, so a lobe playing a raster late leaves
    m0 untouched -- as long as no area crosses a block edge, which a slow y gradient spanning the
    whole block guarantees here.  The *tree* is moved after compiling, which is the same
    discrepancy a compiler that mis-placed the event would produce.
    """
    g = pp.make_trapezoid('x', area=500.0, duration=1e-3, system=opts)
    hold = pp.make_trapezoid('y', area=20.0, duration=4e-3, system=opts)
    tree = sc.LogicBlock('t').add(1e-3, g).add(0.0, hold)
    _, seen = compile_capturing(tree, opts)        # clean to begin with, or this would have raised

    with pytest.raises(sc.CompilerContractError) as err:
        recheck(seen, placed=_late(seen['placed'], 'x', 10e-6))
    text = str(err.value)
    assert 'compiled m1' in text, f'm1 must notice a 10 us displacement; got {text}'
    assert 'compiled m0' not in text, (
        'and m0 must not -- area is exactly what a time shift preserves, which is the whole '
        'reason m1 was added'
    )


def test_m0_still_catches_a_lost_lobe(opts) -> None:
    """The original invariant must keep working: dropping area is what it is for."""
    g = pp.make_trapezoid('x', area=500.0, duration=1e-3, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g)
    _, seen = compile_capturing(tree, opts)

    doubled = [*seen['placed'], *[p for p in seen['placed'] if p.kind == 'trap']]
    with pytest.raises(sc.CompilerContractError, match='compiled m0'):
        recheck(seen, placed=doubled)


def test_registration_catches_a_waveform_pypulseq_stores_wrong(opts) -> None:
    """
    The second stage: what PyPulseq stored, against what emission handed it.

    The stored trapezoid's amplitude is changed after the compile, 1e-3 relative -- a hundred
    times what PyPulseq's six-digit key can do.  The refusal names registration, not the compile.
    """
    g = pp.make_trapezoid('x', area=500.0, duration=1e-3, system=opts)
    seq, seen = compile_capturing(sc.LogicBlock('t').add(0.0, g), opts)

    grad_id = int(seq.block_events[1][2])
    amplitude, *rest = seq.grad_library.data[grad_id]
    seq.grad_library.data[grad_id] = (amplitude * (1.0 + 1e-3), *rest)
    seq.block_cache.clear()
    with pytest.raises(sc.CompilerContractError) as err:
        recheck(seen)
    text = str(err.value)
    assert 'registered m0' in text and 'registration tolerance' in text
    assert 'compiled m' not in text, 'the compile handed PyPulseq the right waveform'


def test_the_address_check_catches_a_label_on_the_wrong_readout(opts) -> None:
    """
    The check the duplicate-address test cannot do.

    Retargeting the label to the *first* ADC is exactly what containment used to do, and the
    resulting addressing is still unique -- so `_label_issues` sees nothing.
    """
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = (
        sc.LogicBlock('t')
        .add(0.0, adc)
        .add(2000e-6, pp.make_label('LIN', 'SET', 7))
        .add(5000e-6, adc)
    )
    _, seen = compile_capturing(tree, opts)    # clean to begin with, or this would have raised

    placed = seen['placed']
    adc_starts = sorted(p.res_start for p in placed if p.kind == 'adc')
    wrong = {
        i: adc_starts[0]                       # claim the label belongs to the first readout
        for i, p in enumerate(placed)
        if p.kind in ('labelset', 'labelinc')
    }
    with pytest.raises(sc.CompilerContractError) as err:
        recheck(seen, targets=wrong)
    text = str(err.value)
    assert 'wrong readout' in text, (
        f'the address check must notice a label attributed to the wrong readout; got {text}'
    )
    assert 'LIN' in text


def test_the_duplicate_address_check_would_not_have_caught_it(opts) -> None:
    """
    Why the address invariant had to be added.

    A labelling shifted by one readout stays unique, so the duplicate check passes -- which is
    how the label corruption escaped notice.
    """
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = sc.LogicBlock('t')
    for i in range(3):
        tree.add(i * 5e-3, adc, pp.make_label('LIN', 'SET', i + 1))
    # The compile returns, so check_label_addresses saw nothing: the addresses are unique.
    seq = sc.compile(tree, opts)
    labels = seq.evaluate_labels(evolution='adc')
    seen = [int(v) for v in np.atleast_1d(np.asarray(labels['LIN']))]
    assert len(set(seen)) == len(seen), 'unique, so the duplicate check is silent'
    assert seen == [1, 2, 3], 'and here they are also *right*, which is the other half of it'


def test_duration_is_still_checked(opts) -> None:
    """The oldest invariant, and the one that fences boundary merging (W6)."""
    tree = sc.LogicBlock('t').add(0.0, pp.make_trapezoid('x', area=100.0, system=opts))
    seq, seen = compile_capturing(tree, opts)
    with pytest.raises(sc.CompilerContractError, match='differs from the tree total'):
        recheck(seen, tree_duration_s=float(seq.duration()[0]) + 1e-3)
