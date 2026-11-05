"""
A raster resample is an approximation with a derived bound, and anything beyond it is refused.

Legalization resamples a block-axis onto the gradient raster when no pulseq gradient can hold its
knots: a trapezoid, whose knots are on raster edges, summed with a gradient sampled at raster
centres.  Keeping only the raster-edge values replaces each raster interval's centre-sample tent by
its chord.  Over a stretch of samples the area changes telescope, so the change is two end terms
(ADR-008):

    dm0 = (r/8) (D_end - D_start)
    dm1 = (r/8) (t_b D_end - t_a D_start) - (r^2/8) (g_b - g_a)

``D`` is the sample step across each end of the stretch.  At a gradient's own start or end it uses
the virtual sample its ``first`` or ``last`` implies, and at a cut it uses the next sample.  ``t_a``
and ``t_b`` are the stretch's ends and ``g_a``, ``g_b`` the values there.  Within the slew limit
``|D| <= max_slew * r``, so a stretch moves k by at most ``r^2 max_slew / 4`` at every instant,
whatever its duration, and m1 by at most twice the block's length times that.

A gradient sampled every half raster has no such bound, and what it would move is measured and
refused.  ADR-007's compile stage still compares the emitted block with the tree resampled the same
way: these tests cover the contract on that resample, not the faithfulness of its emission.
"""

from __future__ import annotations

import re
import warnings

import numpy as np
import pypulseq as pp
import pytest
from recheck import compile_capturing

import seqcraft as sc
from seqcraft.compiler.verification import _tree_side
from seqcraft.design.events import knots_of, pwl_moment


def _opts() -> pp.Opts:
    """40 mT/m, 150 T/m/s, 10 us raster: the bound r^2 max_slew / 4 is 0.160 1/m."""
    return pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)


def _bound(opts: pp.Opts) -> float:
    return float(opts.grad_raster_time) ** 2 * float(opts.max_slew) / 4.0


def _centre(axis: str, samples: np.ndarray, opts: pp.Opts):
    return pp.make_arbitrary_grad(axis, np.asarray(samples, dtype=float), first=0.0, last=0.0,
                                  system=opts)


def _triangle(opts: pp.Opts, sign: float = 1.0) -> np.ndarray:
    """Raster-centre samples climbing and falling at 0.999 of the slew limit: 37.5 mT/m peak."""
    step = 0.999 * float(opts.max_slew) * float(opts.grad_raster_time)
    up = (np.arange(25) + 0.5) * step
    return sign * np.concatenate([up, up[::-1]])


def _holder(axis: str, opts: pp.Opts, duration: float = 3e-3):
    """A trapezoid on `axis` spanning `duration`.  On the same axis, a centre-sampled piece placed on
    its flat top forces the resample without adding slew; on another axis, it keeps the block open."""
    return pp.make_trapezoid(axis, area=50.0, duration=duration, system=opts)


#: Where a triangle sits on a same-axis holder: inside the flat top, so their sum stays legal.
ON_TOP = 1e-3


def _compile(tree: sc.LogicBlock, opts: pp.Opts):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        seq, seen = compile_capturing(tree, opts)
    raised = [str(w.message) for w in caught if issubclass(w.category, sc.SeqCraftWarning)]
    return seq, seen, raised


def _change(seen: dict, block: int, axis: str, opts: pp.Opts) -> tuple[float, float]:
    """What the compile did to one block-axis: the ready block's m0 and m1 minus the tree's."""
    handed = seen['handed']
    raster = round(float(opts.grad_raster_time) * 1e12)
    tree0, tree1, *_ = _tree_side(seen['placed'], handed.edges, (), raster, 0.0)
    a = 'xyz'.index(axis)
    return handed.m0[block, a] - tree0[block, a], handed.m1[block, a] - tree1[block, a]


def _end_terms(samples: np.ndarray, start: int, stop: int, r: float, offset: float):
    """The closed form for samples[start:stop] of a raster-centre gradient with first = last = 0.

    `offset` is where the stretch's first raster edge lies from its block's start."""
    c = np.concatenate([[-samples[0]], samples, [-samples[-1]]])     # virtual end samples
    i, j = start + 1, stop + 1                                        # into c
    d_start, d_end = c[i] - c[i - 1], c[j] - c[j - 1]
    g_a, g_b = (c[i - 1] + c[i]) / 2, (c[j - 1] + c[j]) / 2
    t_a, t_b = offset, offset + (stop - start) * r
    dm0 = r / 8 * (d_end - d_start)
    dm1 = r / 8 * (t_b * d_end - t_a * d_start) - r * r / 8 * (g_b - g_a)
    return dm0, dm1


# ------------------------------------------------------------------------- the derived contract
def test_a_resampled_raster_centre_stretch_moves_m0_and_m1_by_its_end_terms() -> None:
    """
    A raster-centre sine summed with a trapezoid and cut on a slope: two resampled blocks.

    In each, the compiled block's m0 and m1 minus the tree's equal the closed form's end terms,
    and no more.
    """
    opts = _opts()
    r = float(opts.grad_raster_time)
    period, n = 320e-6, 256
    samples = 3e-3 * float(opts.gamma) * np.sin(2 * np.pi * (np.arange(n) + 0.5) * r / period)
    tree = (sc.LogicBlock('sine').add(0.0, _centre('x', samples, opts))
            .add(0.5e-3, pp.make_trapezoid('x', area=20.0, duration=1e-3, system=opts))
            .add(1.0e-3, sc.barrier('on a slope')))
    _, seen, raised = _compile(tree, opts)
    assert any('resampled onto the raster' in text for text in raised), raised
    assert sorted(seen['resampled']) == [(0, 'x'), (1, 'x')]

    cut = round(1.0e-3 / r)
    for block, (start, stop) in enumerate([(0, cut), (cut, n)]):
        dm0, dm1 = _change(seen, block, 'x', opts)
        want0, want1 = _end_terms(samples, start, stop, r, 0.0)
        assert dm0 == pytest.approx(want0, rel=1e-9, abs=1e-12)
        assert dm1 == pytest.approx(want1, rel=1e-9, abs=1e-15)
        assert abs(dm0) <= _bound(opts)


def test_a_max_slew_triangle_cut_on_its_way_down_reaches_the_bound() -> None:
    """The largest change one stretch can hold: slopes of opposite sign at the limit at its ends."""
    opts = _opts()
    r = float(opts.grad_raster_time)
    tri = _triangle(opts)
    cut = 35                                       # ten samples down the falling side
    tree = (sc.LogicBlock('tri').add(0.0, _holder('z', opts)).add(ON_TOP, _centre('z', tri, opts))
            .add(ON_TOP + cut * r, sc.barrier('down-slope')))
    _, seen, raised = _compile(tree, opts)
    dm0, _ = _change(seen, 0, 'z', opts)
    assert dm0 == pytest.approx(-0.999 * _bound(opts), rel=1e-9)
    assert dm0 == pytest.approx(_end_terms(tri, 0, cut, r, ON_TOP)[0], rel=1e-9)
    assert any('resampled onto the raster' in text for text in raised)


def test_the_resample_warning_states_how_far_k_and_m1_moved() -> None:
    opts = _opts()
    r = float(opts.grad_raster_time)
    tree = (sc.LogicBlock('tri').add(0.0, _holder('z', opts))
            .add(ON_TOP, _centre('z', _triangle(opts), opts))
            .add(ON_TOP + 35 * r, sc.barrier('down-slope')))
    _, seen, raised = _compile(tree, opts)
    text = next(t for t in raised if 'resampled onto the raster' in t)
    k = float(re.search(r'k by up to ([-+0-9.e]+) 1/m of a ([-+0-9.e]+) 1/m bound', text).group(1))
    dm0, dm1 = _change(seen, 0, 'z', opts)
    assert k == pytest.approx(0.999 * _bound(opts), rel=1e-2)
    assert f'block m0 {dm0:+.3g} 1/m and m1 {dm1:+.3g} s/m' in text


def test_repeated_resamples_report_the_net_and_the_summed_magnitude() -> None:
    """Two resamples that cancel in m0: the net hides them, the summed magnitude does not."""
    opts = _opts()
    r = float(opts.grad_raster_time)
    tree = sc.LogicBlock('pair')
    for j, sign in enumerate((1.0, -1.0)):
        t0 = j * 4e-3
        tree.add(t0, _holder('z', opts)).add(t0 + ON_TOP, _centre('z', _triangle(opts, sign), opts))
        tree.add(t0 + ON_TOP + 35 * r, sc.barrier(f'down-slope {j}'))
    _, _, raised = _compile(tree, opts)
    total = next(t for t in raised if 'resampled in more than one block' in t)
    # Each triangle's cut leaves two resampled blocks: the climb and ten samples down, which
    # carries the whole end term, and the rest of the way down, which carries none.
    found = re.search(r'axis z over (\d+) blocks: m0 ([-+0-9.e]+) 1/m net; k ([-+0-9.e]+) 1/m', total)
    assert int(found.group(1)) == 4
    assert abs(float(found.group(2))) < 1e-9
    assert float(found.group(3)) == pytest.approx(2 * 0.999 * _bound(opts), rel=1e-2)


# ------------------------------------------------------------------- beyond the derivation
def _half_raster(opts: pp.Opts, signs: np.ndarray):
    """Sampled every half raster: zero on the centres, +-slew*r/2 on the edges, each step legal."""
    samples = np.zeros(2 * signs.size - 1)
    samples[1::2] = signs[:-1] * 0.999 * float(opts.max_slew) * float(opts.grad_raster_time) / 2
    return pp.make_arbitrary_grad('y', samples, first=0.0, last=0.0, oversampling=True, system=opts)


def test_a_resample_beyond_the_raster_bound_is_refused() -> None:
    """A half-raster waveform keeps its edge values, so resampling it adds area every interval."""
    opts = _opts()
    tree = (sc.LogicBlock('half').add(0.0, _holder('y', opts))
            .add(0.5e-3, _half_raster(opts, np.ones(200))))
    with pytest.raises(sc.CompileError, match=r'(?s)within the resample bound.*every half raster'):
        sc.compile(tree, opts)


def test_a_block_whose_m0_cancels_is_refused_when_k_strays() -> None:
    """
    The bound is on the running k, not on the block's m0, because m1 follows the running k.

    Half the edges up and half down: the block's m0 change cancels to within the bound, while k
    climbs by tens of times the bound and comes back.
    """
    opts = _opts()
    r = float(opts.grad_raster_time)
    over = _half_raster(opts, np.where(np.arange(200) < 100, 1.0, -1.0))
    t, g = knots_of(over, 0.0)
    edges = np.arange(int(round(t[-1] / r)) + 1) * r
    union = np.union1d(t, edges)
    d = np.interp(union, edges, np.interp(edges, t, g)) - np.interp(union, t, g)
    running = np.concatenate([[0.0], np.cumsum(np.diff(union) * (d[1:] + d[:-1]) / 2)])
    assert abs(running[-1]) <= _bound(opts) < 0.05 * np.max(np.abs(running))
    assert abs(pwl_moment(union, d, 1)) > 1e-3        # m1 moved: 0.016 s/m

    tree = sc.LogicBlock('cancel').add(0.0, _holder('y', opts)).add(0.5e-3, over)
    with pytest.raises(sc.CompileError, match='within the resample bound'):
        sc.compile(tree, opts)


def test_a_smooth_half_raster_sine_is_resampled_within_the_bound() -> None:
    """The guard measures; it does not refuse half-raster input as such."""
    opts = _opts()
    r = float(opts.grad_raster_time)
    sine = 3e-3 * float(opts.gamma) * np.sin(2 * np.pi * (np.arange(383) + 1) * r / 2 / 320e-6)
    over = pp.make_arbitrary_grad('y', sine, first=0.0, last=0.0, oversampling=True, system=opts)
    tree = sc.LogicBlock('smooth').add(0.0, _holder('y', opts)).add(0.5e-3, over)
    _, seen, raised = _compile(tree, opts)
    assert seen['resampled'] == ((0, 'y'),)
    assert any('resampled onto the raster' in text for text in raised)


# ---------------------------------------------------------------- held exactly, not resampled
def test_a_raster_centre_gradient_starting_inside_its_block_is_split_exactly() -> None:
    """
    A centre-sampled gradient that starts inside its block and is then cut is the same event with
    a delay, so it is held exactly rather than resampled.  So is one that ends before its block
    does.
    """
    opts = _opts()
    r = float(opts.grad_raster_time)
    tri = _triangle(opts)
    tree = (sc.LogicBlock('delayed').add(0.0, _holder('x', opts))
            .add(250e-6, _centre('z', tri, opts)).add(250e-6 + 35 * r, sc.barrier('cut'))
            .add(250e-6 + 35 * r, _holder('x', opts)))
    _, seen, raised = _compile(tree, opts)
    assert not any('resampled' in text for text in raised), raised
    assert seen['resampled'] == ()
    for block in (0, 1):
        dm0, dm1 = _change(seen, block, 'z', opts)
        assert abs(dm0) < 1e-12 and abs(dm1) < 1e-15


def test_a_delayed_run_of_four_samples_is_resampled_not_stored_short() -> None:
    """
    PyPulseq stores a centre-sampled shape of four samples or fewer without its last half raster,
    and registration refuses that.  Four samples before a cut stay a bounded resample.
    """
    opts = _opts()
    r = float(opts.grad_raster_time)
    tree = (sc.LogicBlock('short').add(0.0, _holder('x', opts))
            .add(250e-6, _centre('z', _triangle(opts), opts))
            .add(250e-6 + 4 * r, sc.barrier('cut')).add(250e-6 + 4 * r, _holder('x', opts)))
    _, seen, raised = _compile(tree, opts)
    assert seen['resampled'] == ((0, 'z'),)
    dm0, _ = _change(seen, 0, 'z', opts)
    assert abs(dm0) <= _bound(opts)
    assert any('resampled onto the raster' in text for text in raised)
