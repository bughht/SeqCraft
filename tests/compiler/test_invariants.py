"""
The invariants that run on every compile -- and proof that each one can fail.

`verify_against_tree` is the last line of defence: it runs for every user, on every sequence,
including the ones no test thought of.  So each check gets two tests: it passes on a good compile,
and it *fires* on a deliberately broken one.  A check that cannot fail is worse than no check,
because it reads like coverage.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.compiler.boundaries import label_targets
from seqcraft.compiler.placement import place_events
from seqcraft.compiler.verification import _sequence_moments, verify_against_tree
from seqcraft.design.events import knots_of as _grad_knots
from seqcraft.design.events import pwl_moment


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
        errors.append(abs(float(np.trapezoid(g * t, t)) - exact))
    assert errors[0] > 0.0, 'trapz must actually be wrong here, or the test proves nothing'
    assert errors[0] > errors[1] > errors[2], 'and refining it must approach the closed form'
    assert errors[-1] < 1e-4 * exact


def test_pwl_m1_shifts_by_area_times_offset() -> None:
    """The property the invariant relies on: displacing a waveform changes m1 by ``m0 * dt``."""
    times = np.array([0.0, 100e-6, 300e-6, 400e-6])
    amps = np.array([0.0, 1e4, 1e4, 0.0])
    m0 = float(np.trapezoid(amps, times))
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
def _recheck(out, placed, targets) -> None:
    """
    Re-run the invariants against a (possibly tampered) tree.

    Returns nothing: the invariants raise now, so ``sc.compile`` having returned at all is the
    "clean to begin with" assertion each test below used to make explicitly.
    """
    verify_against_tree(
        placed,
        targets,
        duration_s=out.duration_s,
        tree_duration_s=out.tree_duration_s,
        moments=lambda order: _sequence_moments(out.seq, order),
        label_states=lambda: out.seq.evaluate_labels(evolution='adc'),
    )


def test_m1_catches_a_gradient_that_plays_at_the_wrong_time(opts) -> None:
    """
    The failure m0 is blind to.

    Area is exactly the quantity a time shift preserves, so a lobe playing a raster early leaves
    m0 untouched.  Here the *tree* is moved after compiling, which is the same discrepancy a
    compiler that mis-placed the event would produce.
    """
    g = pp.make_trapezoid('x', area=500.0, duration=1e-3, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g).add(3e-3, pp.make_delay(1e-3))
    out = sc.compile(tree, opts)              # clean to begin with, or this would have raised

    placed = place_events(tree, opts)
    shifted = [
        dataclasses.replace(p, node_t=p.node_t + 10e-6, start=p.start + 10e-6,
                            end=p.end + 10e-6, res_start=p.res_start + 10e-6,
                            res_end=p.res_end + 10e-6)
        if p.kind == 'trap' else p
        for p in placed
    ]
    with pytest.raises(sc.CompilerContractError) as err:
        _recheck(out, shifted, label_targets(shifted))
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
    out = sc.compile(tree, opts)

    placed = place_events(tree, opts)
    doubled = [*placed, *[p for p in placed if p.kind == 'trap']]
    with pytest.raises(sc.CompilerContractError, match='compiled m0'):
        _recheck(out, doubled, label_targets(doubled))


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
    out = sc.compile(tree, opts)              # clean to begin with, or this would have raised
    assert out.check().ok

    placed = place_events(tree, opts)
    adc_starts = sorted(p.res_start for p in placed if p.kind == 'adc')
    wrong = {
        i: adc_starts[0]                       # claim the label belongs to the first readout
        for i, p in enumerate(placed)
        if p.kind in ('labelset', 'labelinc')
    }
    with pytest.raises(sc.CompilerContractError) as err:
        _recheck(out, placed, wrong)
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
    out = sc.compile(tree, opts)
    labels = out.seq.evaluate_labels(evolution='adc')
    seen = [int(v) for v in np.atleast_1d(np.asarray(labels['LIN']))]
    assert len(set(seen)) == len(seen), 'unique, so the duplicate check is silent'
    assert not [i for i in out.check().issues if i.kind == 'label']


def test_duration_is_still_checked(opts) -> None:
    """The oldest invariant, and the one that fences boundary merging (W6)."""
    tree = sc.LogicBlock('t').add(0.0, pp.make_trapezoid('x', area=100.0, system=opts))
    out = sc.compile(tree, opts)
    object.__setattr__(out, 'tree_duration_s', out.tree_duration_s + 1e-3)
    with pytest.raises(sc.CompilerContractError, match='differs from the tree total'):
        _recheck(out, place_events(tree, opts), {})
