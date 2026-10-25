"""
The m0 tolerance scales with the area each event traverses, not with each event's net.

A gradient event that changes sign and nets zero -- a closed trajectory's axis, or a bipolar pair
held in one shape -- used to collapse the tolerance to its absolute floor of 1e-6 1/m, because the
scale summed each event's |net area|.  PyPulseq's registration rounds a compressible shape's
peak-normalised samples to 1e-7, which moves an ordinary waveform's area by a few 1e-6 to a few
1e-5 1/m.  So legal input was refused whenever compression engaged and the shape had no symmetry
to cancel the rounding.

Which of two identical sequences was refused also depended on packaging.  The same waveform held
as two overlapping events, whose nets do not vanish one by one, compiled to the identical emitted
shape with the identical rounding, and passed.

Traversed area is additive over a division in time and never smaller than the waveform's own, so
one waveform now gets one verdict however the tree divides it.  The check still refuses area that
is actually lost, and m1 -- deliberately unchanged -- still refuses a gradient played at the wrong
time.  One gap remains, pinned at the bottom of this file.

Nothing here builds a module.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pypulseq as pp
import pytest
from pypulseq.compress_shape import compress_shape

import seqcraft as sc
import seqcraft.compiler as compiler_package
from seqcraft.compiler.boundaries import label_targets
from seqcraft.compiler.placement import place_events
from seqcraft.compiler.verification import _sequence_moments, _traversed_area, verify_against_tree
from seqcraft.design.events import knots_of, pwl_moment

#: PyPulseq's quantum for a peak-normalised shape, applied when run-length compression engages.
QUANTUM = 1e-7

#: Legal zero-net compressible gradients the collapsed floor refused: (seed, smooth samples, peak).
REFUSED_BEFORE = [(0, 800, 0.5), (1, 2500, 0.9), (4, 2500, 0.5), (7, 800, 0.9)]


def _opts() -> pp.Opts:
    return pp.Opts(max_grad=32, grad_unit='mT/m', max_slew=120, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)


def _net(values: np.ndarray) -> float:
    return float(np.sum(0.5 * (values[1:] + values[:-1])))


def _zero_net_section(seed: int, n_smooth: int) -> np.ndarray:
    """
    Exact-slope ramps and flats, which give compression runs to shorten, around a random-phase
    smooth section that breaks every symmetry which could cancel the rounding.  A raised-cosine
    bump makes the net exactly zero.
    """
    rng = np.random.default_rng(seed)
    ramp = np.arange(65) / 64.0
    s = np.arange(1, n_smooth + 1) / n_smooth
    smooth = np.cos(np.pi * s) + np.sin(np.pi * s) ** 2 * (
        0.3 * np.sin(6 * np.pi * s + rng.uniform(0, 6))
        + 0.2 * np.sin(14 * np.pi * s + rng.uniform(0, 6)))
    body = np.concatenate([ramp, np.ones(40), smooth, -np.ones(40), -ramp[::-1]])
    bump = np.zeros_like(body)
    k = np.arange(ramp.size + 40, ramp.size + 40 + n_smooth)
    bump[k] = np.sin(np.pi * (k - k[0] + 1) / (n_smooth + 1)) ** 2
    return body - _net(body) / _net(bump) * bump


def _lobe_then_section(seed: int) -> tuple[np.ndarray, int]:
    """A positive lobe, an exact zero, then a section netting the opposite area; and that zero's index."""
    ramp = np.arange(65) / 64.0
    lobe = np.concatenate([ramp, np.ones(200), ramp[::-1]])
    section = _zero_net_section(seed, 1500)
    bump = np.sin(np.pi * np.arange(section.size) / (section.size - 1)) ** 2
    section = section - _net(lobe) / _net(bump) * bump
    return np.concatenate([lobe, section[1:]]), lobe.size - 1


def _extended(times: np.ndarray, amps: np.ndarray, opts: pp.Opts):
    return pp.make_extended_trapezoid('y', times=times, amplitudes=amps, system=opts)


def _gradient(shape: np.ndarray, fraction: float, opts: pp.Opts):
    amps = fraction * opts.max_grad * shape / np.max(np.abs(shape))
    return _extended(np.arange(amps.size) * opts.grad_raster_time, amps, opts)


def _two_events(opts: pp.Opts) -> sc.LogicBlock:
    whole, cut = _lobe_then_section(1)
    amps = 0.8 * opts.max_grad * whole / np.max(np.abs(whole))
    t = np.arange(amps.size) * opts.grad_raster_time
    return (sc.LogicBlock('two')
            .add(0.0, _extended(t[:cut + 1], amps[:cut + 1], opts))
            .add(float(t[cut]), _extended(t[cut:] - t[cut], amps[cut:], opts)))


def _compile_unchecked(tree: sc.LogicBlock, opts: pp.Opts):
    """Compile with the against-tree check switched off, to measure what it measures."""
    original = compiler_package.verify_against_tree
    compiler_package.verify_against_tree = lambda *a, **k: None       # the check is under test
    try:
        return sc.compile(tree, opts, name='m0')
    finally:
        compiler_package.verify_against_tree = original


def _recheck(seq, placed) -> None:
    """Re-run the invariants on a compiled sequence against a (possibly tampered) tree."""
    duration = float(seq.duration()[0])
    verify_against_tree(placed, label_targets(placed), duration_s=duration,
                        tree_duration_s=duration,
                        moments=lambda order: _sequence_moments(seq, order),
                        label_states=lambda: seq.evaluate_labels(evolution='adc'))


# ------------------------------------------------------------------------------------ regression
@pytest.mark.parametrize(('seed', 'n_smooth', 'fraction'), REFUSED_BEFORE)
def test_a_legal_zero_net_compressible_gradient_compiles(seed: int, n_smooth: int,
                                                         fraction: float) -> None:
    """The regression, stated as the thing that was broken: all of these are legal."""
    opts = _opts()
    grad = _gradient(_zero_net_section(seed, n_smooth), fraction, opts)
    shape = np.asarray(grad.waveform) / np.max(np.abs(grad.waveform))
    assert len(compress_shape(shape).data) < shape.size, 'the case must exercise compression'

    sc.compile(sc.LogicBlock('zero_net').add(0.0, grad), opts, name='legal')


def test_the_area_that_moved_is_the_rounding_of_the_stored_shape() -> None:
    """
    What the tolerance has to cover, measured rather than assumed.

    The compiled area is the area of the peak-normalised samples rounded half away from zero to
    1e-7 -- reproduced here independently of PyPulseq's code.  It exceeds the old floor, which is
    why these were refused.  It is a few percent of the worst case, half a quantum per sample, as
    independent rounding errors should be.  And it is two orders below the new tolerance.
    """
    opts = _opts()
    grad = _gradient(_zero_net_section(1, 2500), 0.9, opts)
    times, amps = knots_of(grad, 0.0)
    peak = float(np.max(np.abs(amps)))
    scaled = amps / peak / QUANTUM
    rounded = np.sign(scaled) * np.floor(np.abs(scaled) + 0.5) * QUANTUM * peak

    seq = _compile_unchecked(sc.LogicBlock('m').add(0.0, grad), opts)
    moved = _sequence_moments(seq, 0)['y'] - pwl_moment(times, amps, 0)

    # Each side sums ~1e4 1/m of partial areas, so float alone is a few 1e-12; 1e-10 is five
    # orders below the effect.
    assert moved == pytest.approx(pwl_moment(times, rounded, 0) - pwl_moment(times, amps, 0),
                                  abs=1e-10)
    assert abs(moved) > 1e-6, 'it must exceed the old floor, or this is not the regression'
    worst = 0.5 * QUANTUM * peak * float(times[-1] - times[0])
    assert abs(moved) < 0.1 * worst
    assert abs(moved) < 1e-2 * 1e-6 * _traversed_area(times, amps)


# ------------------------------------------------------------------------------------- packaging
def test_traversed_area_is_exact_through_a_zero_and_additive_at_knots() -> None:
    times, amps = np.array([0.0, 1.0, 3.0]), np.array([2.0, -2.0, 4.0])
    # 0..1 crosses zero at 0.5: two triangles of 0.5 each.  1..3 crosses at 1 + 2/3: 2/3 + 8/3.
    assert _traversed_area(times, amps) == pytest.approx(1.0 + 10.0 / 3.0, rel=1e-15)
    assert (_traversed_area(times[:2], amps[:2]) + _traversed_area(times[1:], amps[1:])
            == pytest.approx(_traversed_area(times, amps), rel=1e-15))


@pytest.mark.parametrize('seed', [1, 4])
def test_one_waveform_gets_one_verdict_however_the_tree_divides_it(seed: int) -> None:
    """
    The packaging dependence, and what removes it.

    One event, two events cut at an exact zero, two overlapping pieces, one event cut by a block
    barrier: the emitted waveform is the same, so the verdict must be.  Under the sum of each
    event's |net| the first and last were refused and the middle two passed.
    """
    opts = _opts()
    whole, cut = _lobe_then_section(seed)
    amps = 0.8 * opts.max_grad * whole / np.max(np.abs(whole))
    t = np.arange(amps.size) * opts.grad_raster_time
    s = 0.5 * (1.0 - np.cos(np.pi * np.arange(amps.size) / (amps.size - 1)))
    trees = {
        'one event': sc.LogicBlock('w').add(0.0, _extended(t, amps, opts)),
        'cut at a zero': (sc.LogicBlock('w')
                          .add(0.0, _extended(t[:cut + 1], amps[:cut + 1], opts))
                          .add(float(t[cut]), _extended(t[cut:] - t[cut], amps[cut:], opts))),
        'overlapping': (sc.LogicBlock('w')
                        .add(0.0, _extended(t, amps * s, opts))
                        .add(0.0, _extended(t, amps * (1.0 - s), opts))),
        'barrier': (sc.LogicBlock('w').add(0.0, _extended(t, amps, opts))
                    .add(float(t[amps.size // 2]), sc.barrier())),
    }
    for name, tree in trees.items():
        if name == 'overlapping':
            with pytest.warns(sc.SeqCraftWarning, match='same-axis gradient merge'):
                sc.compile(tree, opts, name='overlapping')
        else:
            sc.compile(tree, opts, name=name.replace(' ', '_'))


# ----------------------------------------------------------------------------- negative controls
def test_area_actually_lost_on_a_sign_changing_axis_is_still_refused() -> None:
    """
    The wider tolerance must not hide what m0 exists to catch.

    The compiled sequence is checked against a tree asking for more than it emitted: a whole lobe
    more, and then only a sliver -- 0.05 1/m, one percent of a k-space step at 200 mm, and still
    a few times this axis's tolerance.  The message reports the difference and the tolerance, and
    leaves the cause to the reader rather than guessing a split or a merge.
    """
    opts = _opts()
    tree = _two_events(opts)
    seq = sc.compile(tree, opts, name='clean')                # clean to begin with
    placed = place_events(tree, opts)

    lobe = min((p for p in placed if p.kind in ('trap', 'grad')), key=lambda p: p.node_t)
    with pytest.raises(sc.CompilerContractError) as err:
        _recheck(seq, [*placed, lobe])
    text = str(err.value)
    assert 'compiled m0' in text and 'beyond the tolerance' in text and 'traverses' in text
    assert 'split or a merge' not in text

    sliver = pp.make_trapezoid('y', area=0.05, duration=200e-6, system=opts)
    with pytest.raises(sc.CompilerContractError, match='compiled m0'):
        _recheck(seq, place_events(tree.add(0.0, sliver), opts))


def test_m1_still_refuses_a_piece_played_one_raster_late() -> None:
    """m1 is unchanged by this fix, and on the waveforms it concerns it still sees a time shift."""
    opts = _opts()
    tree = _two_events(opts)
    seq = sc.compile(tree, opts, name='clean')
    placed = place_events(tree, opts)
    first = min(p.node_t for p in placed if p.kind in ('trap', 'grad'))
    dt = float(opts.grad_raster_time)
    shifted = [dataclasses.replace(p, node_t=p.node_t + dt, start=p.start + dt, end=p.end + dt,
                                   res_start=p.res_start + dt, res_end=p.res_end + dt)
               if p.kind in ('trap', 'grad') and p.node_t == first else p
               for p in placed]

    with pytest.raises(sc.CompilerContractError) as err:
        _recheck(seq, shifted)
    assert 'compiled m1' in str(err.value)
    assert 'compiled m0' not in str(err.value), 'a time shift preserves area exactly'


# ------------------------------------------------------------------------------- the gap left open
@pytest.mark.xfail(strict=True, raises=sc.CompilerContractError, reason=(
    'Known: 1e-6 of the traversed area covers worst-case storage rounding only while a shape '
    "spends at least 5 % of its peak on average.  Below that, a long constant section rounds the "
    'same way at every sample and the error accumulates.  Remove this mark when a storage-quantum '
    'floor covers it.'))
def test_a_dense_low_duty_shape_with_systematic_rounding() -> None:
    """
    A 0.5 ms spike cancelled by a long weak constant tail, in one densely sampled shape.

    Mean |g| is 2 % of the peak.  The tail's level is chosen so its samples round almost half a
    quantum each way the same way, so the stored area moves by nearly the worst case.
    """
    opts = _opts()
    n = max(range(3000, 3040), key=lambda m: abs(30.0 / m / QUANTUM - round(30.0 / m / QUANTUM)))
    ramp = np.arange(21) / 20.0
    shape = np.concatenate([ramp, np.ones(9), ramp[::-1], -np.full(n, 30.0 / n), [0.0]])
    grad = _extended(np.arange(shape.size) * opts.grad_raster_time,
                     20e-3 * opts.gamma * shape, opts)
    sc.compile(sc.LogicBlock('low_duty').add(0.0, grad), opts, name='low_duty')
