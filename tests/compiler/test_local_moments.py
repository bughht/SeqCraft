"""
m0 and m1 are compared per compiled block, each about the block's own start, in two stages.

Referencing m1 to the start of the sequence made a legal waveform's verdict depend on where it
sat: registration moves a stored shape's area by ``dm0``, which moves the absolute m1 by
``t0 * dm0``, and a compressible shape accepted at 0 s was refused at 1 s.  No tolerance alone
fixes that without hiding a blip played one raster late far into a sequence.  About each block's
own start the lever arm is at most a block long, so the same waveform compares the same anywhere.

The two stages keep two error sources apart.  **Compile**: the tree's events, restricted exactly
to each block, against the ready blocks emission hands PyPulseq; its tolerance is SeqCraft's own
arithmetic.  **Registration**: those ready blocks against what PyPulseq stored; its tolerance is
PyPulseq's storage precision.  A refusal says which one moved the waveform.

Every negative control tampers the expected side and re-runs the check ``sc.compile`` ran, on the
inputs it ran with (``recheck.py``).  Nothing here builds a module.
"""

from __future__ import annotations

import dataclasses
import warnings

import numpy as np
import pypulseq as pp
import pytest
from recheck import compile_capturing, recheck

import seqcraft as sc
from seqcraft.compiler.placement import place_events
from seqcraft.compiler.verification import _sequence_moments, _tree_side
from seqcraft.design.events import knots_of, pwl_moment

#: Where the controls are placed: from the start of the sequence to 100 s into it.
PLACEMENTS = [0.0, 1.0, 10.0, 100.0]


def _opts() -> pp.Opts:
    """32 mT/m and 120 T/m/s, as in ``test_m0_tolerance.py``."""
    return pp.Opts(max_grad=32, grad_unit='mT/m', max_slew=120, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)


def _net(values: np.ndarray) -> float:
    return float(np.sum(0.5 * (values[1:] + values[:-1])))


def _compressible(seed: int, n_smooth: int, *, zero_net: bool = True) -> np.ndarray:
    """Exact-slope ramps and flats around a random-phase smooth section; optionally zero net."""
    rng = np.random.default_rng(seed)
    ramp = np.arange(65) / 64.0
    s = np.arange(1, n_smooth + 1) / n_smooth
    smooth = np.cos(np.pi * s) + np.sin(np.pi * s) ** 2 * (
        0.3 * np.sin(6 * np.pi * s + rng.uniform(0, 6))
        + 0.2 * np.sin(14 * np.pi * s + rng.uniform(0, 6)))
    if not zero_net:
        return np.concatenate([ramp, np.ones(400), smooth, -np.ones(40), -ramp[::-1]])
    body = np.concatenate([ramp, np.ones(40), smooth, -np.ones(40), -ramp[::-1]])
    bump = np.zeros_like(body)
    k = np.arange(ramp.size + 40, ramp.size + 40 + n_smooth)
    bump[k] = np.sin(np.pi * (k - k[0] + 1) / (n_smooth + 1)) ** 2
    return body - _net(body) / _net(bump) * bump


def _extended(opts: pp.Opts, amps: np.ndarray, axis: str = 'y'):
    return pp.make_extended_trapezoid(axis, times=np.arange(amps.size) * opts.grad_raster_time,
                                      amplitudes=amps, system=opts)


def _shape(opts: pp.Opts, shape: np.ndarray, fraction: float = 0.9, axis: str = 'y'):
    return _extended(opts, fraction * opts.max_grad * shape / np.max(np.abs(shape)), axis)


def _blip(opts: pp.Opts):
    """4.5 1/m: one raster late moves m1 by 4.5e-5 s/m."""
    return pp.make_trapezoid('y', area=4.5, system=opts)


def _holder(opts: pp.Opts, duration: float):
    """A gentle x gradient spanning `duration`, so no block boundary can fall inside it."""
    return pp.make_trapezoid('x', area=50.0, duration=duration, system=opts)


def _end(event) -> float:
    return float(knots_of(event, 0.0)[0][-1])


def _late(placed, match, dt: float) -> list:
    """The placed events, with those `match` selects moved `dt` later."""
    return [dataclasses.replace(p, node_t=p.node_t + dt, start=p.start + dt, end=p.end + dt,
                                res_start=p.res_start + dt, res_end=p.res_end + dt)
            if match(p) else p for p in placed]


def _gradients(tree: sc.LogicBlock, opts: pp.Opts) -> list:
    return [p for p in place_events(tree, opts) if p.kind in ('trap', 'grad')]


def _is_blip(p) -> bool:
    return p.kind == 'trap' and p.event.channel == 'y' and abs(p.event.area - 4.5) < 1e-9


def _compile_warned(tree: sc.LogicBlock, opts: pp.Opts, *expected: str):
    """Compile capturing the check's inputs, and assert each expected warning was raised."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        result = compile_capturing(tree, opts)
    raised = ' '.join(str(w.message) for w in caught if issubclass(w.category, sc.SeqCraftWarning))
    assert all(text in raised for text in expected), raised
    return result


# ------------------------------------------------------------------- correct global translation
@pytest.mark.parametrize('t0', PLACEMENTS)
@pytest.mark.parametrize('zero_net', [True, False], ids=['zero net', 'non-zero net'])
def test_a_compressible_shape_compiles_wherever_it_sits(t0: float, zero_net: bool) -> None:
    """
    The translation reproducers: the global m1 refused both shapes at 1, 10 and 100 s.

    Compression engages, so registration moves the stored area by about 3e-5 1/m.
    """
    opts = _opts()
    grad = _shape(opts, _compressible(1 if zero_net else 5, 2500, zero_net=zero_net))
    sc.compile(sc.LogicBlock('shape').add(t0, grad), opts, name='translated')


def test_a_shape_compares_the_same_wherever_it_sits() -> None:
    """
    Both stages measure the same difference for the same waveform at 0 s and at 100 s.

    Registration stores the same shape at both placements, and about each block's start nothing
    else enters, so stored minus handed is identical.  The compile stage agrees to float.
    """
    opts = _opts()
    grad = _shape(opts, _compressible(1, 2500))
    measured = []
    for t0 in (0.0, 100.0):
        _, seen = compile_capturing(sc.LogicBlock('shape').add(t0, grad), opts)
        handed, stored = seen['handed'], seen['stored']()
        row = int(np.flatnonzero(handed.error0[:, 1])[0])
        measured.append((stored.m0[row, 1] - handed.m0[row, 1],
                         stored.m1[row, 1] - handed.m1[row, 1]))
    assert abs(measured[0][0]) > 1e-6, 'registration must actually move the shape here'
    assert measured[0] == measured[1]


@pytest.mark.parametrize('t0', [0.0, 10.0])
def test_block_moments_sum_to_the_whole_sequence_moments(t0: float) -> None:
    """
    Why the global m0 and m1 need no check of their own.

    For any partition, ``M0 = sum m0_k`` and ``M1 = sum (m1_k + a_k m0_k)``, so a sequence whose
    every block agrees cannot disagree as a whole by more than the summed tolerances.  The stored
    side is checked against the whole-sequence integrator here, which integrates at absolute
    times and so carries a little float of its own at 10 s.
    """
    opts = _opts()
    tree = (sc.LogicBlock('sum').add(t0, _shape(opts, _compressible(1, 2500)))
            .add(t0 + 0.02, sc.barrier()).add(t0 + 0.03, _blip(opts)))
    seq, seen = compile_capturing(tree, opts)
    stored = seen['stored']()
    start = (stored.edges[:-1] / 1e12)[:, None]
    whole0, whole1 = _sequence_moments(seq, 0), _sequence_moments(seq, 1)
    assert np.sum(stored.m0[:, 1]) == pytest.approx(whole0['y'], abs=1e-9)
    assert np.sum(stored.m1 + start * stored.m0, axis=0)[1] == pytest.approx(
        whole1['y'], abs=1e-7, rel=1e-12)


# ------------------------------------------------------------------------------------ lost area
@pytest.mark.parametrize('t0', PLACEMENTS)
def test_a_lost_sliver_is_refused_wherever_it_sits(t0: float) -> None:
    """0.05 1/m inside the zero-net shape: one percent of a k-space step at 200 mm."""
    opts = _opts()
    tree = sc.LogicBlock('shape').add(t0, _shape(opts, _compressible(1, 2500)))
    _, seen = compile_capturing(tree, opts)
    sliver = pp.make_trapezoid('y', area=0.05, duration=200e-6, system=opts)
    extra = _gradients(sc.LogicBlock('x').add(t0 + 0.01, sliver), opts)
    with pytest.raises(sc.CompilerContractError, match='compiled m0'):
        recheck(seen, placed=[*seen['placed'], *extra])


# --------------------------------------------------------------- a mistimed piece of non-zero area
@pytest.mark.parametrize('t0', PLACEMENTS)
def test_a_blip_one_raster_late_is_refused_wherever_it_sits(t0: float) -> None:
    """
    4.5 1/m inside a 3 ms block: 4.5e-5 s/m of m1, the size of a 1e-4 s/m tolerance at 100 s.

    The block is held open, so no area crosses an edge and m0 cannot see it: only m1 can.
    """
    opts = _opts()
    tree = sc.LogicBlock('b').add(t0, _holder(opts, 3e-3)).add(t0 + 1e-3, _blip(opts))
    _, seen = compile_capturing(tree, opts)
    with pytest.raises(sc.CompilerContractError) as err:
        recheck(seen, placed=_late(seen['placed'], _is_blip, float(opts.grad_raster_time)))
    assert 'axis y m1' in str(err.value) and 'compiled m0' not in str(err.value)


# --------------------------------------------------------------------- errors across a boundary
@pytest.mark.parametrize('t0', [0.0, 1.0, 10.0])
def test_area_moved_across_a_block_boundary_is_refused(t0: float) -> None:
    """
    0.05 1/m taken from just after a block edge and added just before it.

    The sequence's whole m0 does not change at all, which is why a global check missed it.
    """
    opts = _opts()
    r = float(opts.grad_raster_time)
    crossing = pp.make_trapezoid('y', area=400.0, duration=2e-3, system=opts)
    tree = sc.LogicBlock('c').add(t0, crossing).add(t0 + 1e-3, sc.barrier())
    _, seen = compile_capturing(tree, opts)

    def triangle(area: float):
        return pp.make_extended_trapezoid('y', times=np.array([0.0, r, 2 * r]),
                                          amplitudes=np.array([0.0, area / r, 0.0]), system=opts)

    moved = _gradients(sc.LogicBlock('m').add(t0 + 1e-3 - 2 * r, triangle(0.05))
                       .add(t0 + 1e-3, triangle(-0.05)), opts)
    assert sum(pwl_moment(*knots_of(p.event, p.node_t), 0) for p in moved) == pytest.approx(
        0.0, abs=1e-12)
    with pytest.raises(sc.CompilerContractError, match='compiled m0'):
        recheck(seen, placed=[*seen['placed'], *moved])


# ---------------------------------------------------------------------------------- long blocks
@pytest.mark.parametrize('t0', [0.0, 10.0])
def test_a_long_block_compiles_and_still_sees_a_late_blip(t0: float) -> None:
    """
    A 97 ms shape and a blip after it, held in one block.

    Registration error times a 97 ms lever arm would swamp the blip's 4.5e-5 s/m.  It does not
    enter the compile stage at all, which compares the tree with what SeqCraft emitted.
    """
    opts = _opts()
    long_shape = _shape(opts, _compressible(1, 9500))
    tree = (sc.LogicBlock('L').add(t0, long_shape).add(t0 + _end(long_shape) + 1e-3, _blip(opts))
            .add(t0, _holder(opts, _end(long_shape) + 2e-3)))
    _, seen = _compile_warned(tree, opts, 'same-axis gradient merge')
    assert len(seen['handed'].edges) <= 4, 'the shape and the blip must share one long block'
    with pytest.raises(sc.CompilerContractError, match='axis y m1'):
        recheck(seen, placed=_late(seen['placed'], _is_blip, float(opts.grad_raster_time)))


@pytest.mark.parametrize(('top', 'up', 'down', 'cut', 't0'), [
    (1229870.7682268561, 18.98e-3, 2.09e-3, 0.81e-3, 0.0),
    (1180614.015904176, 41.6e-3, 49.56e-3, 1.03e-3, 0.0),
    (422673.65083461127, 26.11e-3, 51.14e-3, 0.7e-3, 10.0),
])
def test_a_ramp_cut_far_from_its_knots_compiles(top: float, up: float, down: float, cut: float,
                                                t0: float) -> None:
    """
    A long ramp up and down, cut into short blocks.

    Near the foot of the ramp down a value at a cut is small, but it is interpolated from a knot
    tens of milliseconds away and up to a million Hz/m: far more than the block itself holds.  Its
    float error scales with that knot.  Scaled by the block's own peak instead, the compile
    tolerance refused each of these legal trees.
    """
    opts = _suite()
    ramp = pp.make_extended_trapezoid('x', times=np.array([0.0, up, up + down]),
                                      amplitudes=np.array([0.0, top, 0.0]), system=opts)
    tree = sc.LogicBlock('ramp').add(t0, ramp)
    k = 1
    while k * cut < up + down:
        tree.add(t0 + k * cut, sc.barrier())
        k += 1
    sc.compile(tree, opts, name='ramp')


@pytest.mark.parametrize(('t0', 'off'), [(100.0, 0.4999), (300.0, 0.49), (4000.0, 0.3)])
def test_a_start_near_half_a_tick_off_compiles(t0: float, off: float) -> None:
    """
    The raster check accepts a gradient start within half a picosecond tick of its tick.

    Superposing turns each knot's absolute float time into ticks itself, and for a start that close
    to half a tick it can round a knot the other way; far into a sequence a block edge can too.  The
    compile tolerance allows a tick for exactly those pieces, and without it each of these legal
    trees was refused.
    """
    opts = _suite()
    t = (round(t0 * 1e12) + off) / 1e12
    s = np.linspace(0.0, 1.0, 400)
    curve = np.sin(np.pi * s) * np.cos(2.0 * np.pi * s) * 0.3 * opts.max_grad
    trap = pp.make_trapezoid('x', area=-35.0, duration=900e-6, system=opts)
    tree = (sc.LogicBlock('near').add(t, pp.make_arbitrary_grad('x', waveform=curve, first=0.0,
                                                                 last=0.0, system=opts))
            .add(t + 1e-3, pp.make_trapezoid('x', area=20.0, duration=600e-6, system=opts))
            .add(t + 0.7e-3, sc.barrier()).add(t + 1.3e-3, sc.barrier())
            .add(t + 5e-3, trap).add(t + 5.3e-3, trap).add(t + 5.5e-3, sc.barrier()))
    for k in range(int(t0 // 50)):                  # PyPulseq refuses a block of 150 s
        tree.add(k * 50.0, sc.barrier())
    _compile_warned(tree, opts, 'resampled')


# ------------------------------------------------------------------------ same-axis superposition
@pytest.mark.parametrize('t0', [0.0, 10.0])
def test_superposed_pieces_compile_and_a_late_one_is_refused(t0: float) -> None:
    """The zero-net waveform held as two overlapping pieces, summed; then one of them late."""
    opts = _opts()
    shape = _compressible(1, 2500)
    amps = 0.9 * opts.max_grad * shape / np.max(np.abs(shape))
    weight = 0.5 * (1.0 - np.cos(np.pi * np.arange(amps.size) / (amps.size - 1)))
    first, second = _extended(opts, amps * weight), _extended(opts, amps * (1.0 - weight))
    tree = (sc.LogicBlock('o').add(t0, first).add(t0, second)
            .add(t0, _holder(opts, _end(first) + 2e-3)))
    _, seen = _compile_warned(tree, opts, 'same-axis gradient merge')
    with pytest.raises(sc.CompilerContractError, match='axis y m1'):
        recheck(seen, placed=_late(seen['placed'], lambda p: p.event is first,
                                   float(opts.grad_raster_time)))


@pytest.mark.parametrize('t0', [0.0, 10.0])
def test_a_blip_summed_onto_a_shape_compiles(t0: float) -> None:
    opts = _opts()
    shape = _shape(opts, _compressible(1, 2500))
    tree = sc.LogicBlock('s').add(t0, shape).add(t0 + _end(shape) - 2e-3, _blip(opts))
    with pytest.warns(sc.SeqCraftWarning, match='same-axis gradient merge'):
        sc.compile(tree, opts, name='summed')


# ---------------------------------------------------------------------------- resampled blocks
def _suite() -> pp.Opts:
    return pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)


def _resampled_stretch(t0: float, opts: pp.Opts) -> sc.LogicBlock:
    """
    A curved arbitrary gradient with a trapezoid summed onto it, cut through the curve.

    The two blocks the trapezoid overlaps cannot hold the sum exactly and are resampled.  The
    stretch they form starts and ends where the curve is non-zero and sloped, so the resample's
    error does not cancel across the stretch.
    """
    s = np.linspace(0.0, 1.0, 400)
    curve = np.sin(np.pi * s) * np.cos(2.0 * np.pi * s) * 0.3 * opts.max_grad
    spiral = pp.make_arbitrary_grad('x', waveform=curve, first=0.0, last=0.0, system=opts)
    trap = pp.make_trapezoid('x', area=20.0, duration=600e-6, system=opts)
    return (sc.LogicBlock('rs').add(t0, spiral).add(t0 + 1e-3, trap)
            .add(t0 + 0.7e-3, sc.barrier()).add(t0 + 1.3e-3, sc.barrier())
            .add(t0 + 2.2e-3, sc.barrier()))


@pytest.mark.parametrize('t0', [0.0, 10.0])
def test_a_resampled_stretch_with_sloped_ends_compiles(t0: float) -> None:
    """
    Each resampled block is compared with the tree resampled the same way, on its own axis.

    Against the tree as given, each block of the stretch differs by about 1e-2 1/m, and the two
    do not cancel: merging the stretch into one interval would have refused it.
    """
    opts = _suite()
    _, seen = _compile_warned(_resampled_stretch(t0, opts), opts, 'resampled')
    resampled = seen['resampled']
    assert len(resampled) == 2 and {axis for _, axis in resampled} == {'x'}

    handed = seen['handed']
    raster = round(float(opts.grad_raster_time) * 1e12)
    plain, _, _, _, _, _ = _tree_side(seen['placed'], handed.edges, (), raster)
    rows = [k for k, _ in resampled]
    per_block = handed.m0[rows, 0] - plain[rows, 0]
    assert np.all(np.abs(per_block) > 5e-3)
    assert abs(np.sum(per_block)) > 5e-3, 'the sloped ends keep the stretch from cancelling'


@pytest.mark.parametrize('t0', [0.0, 10.0])
def test_a_trapezoid_late_inside_a_resampled_block_is_refused(t0: float) -> None:
    opts = _suite()
    _, seen = _compile_warned(_resampled_stretch(t0, opts), opts, 'resampled')
    with pytest.raises(sc.CompilerContractError, match='axis x m'):
        recheck(seen, placed=_late(seen['placed'], lambda p: p.kind == 'trap',
                                   float(opts.grad_raster_time)))


def test_a_resample_on_one_axis_leaves_the_others_checked_block_by_block() -> None:
    """
    Resampling x changes nothing about how y is compared.

    A y blip inside the resampled blocks, played one raster late, is refused on y; so is y area
    moved across the edge between the two resampled blocks.
    """
    opts = _suite()
    r = float(opts.grad_raster_time)
    blip = pp.make_trapezoid('y', area=4.5, system=opts)
    tree = _resampled_stretch(0.0, opts).add(0.8e-3, blip)
    _, seen = _compile_warned(tree, opts, 'resampled')
    assert {axis for _, axis in seen['resampled']} == {'x'}

    with pytest.raises(sc.CompilerContractError, match='axis y m1'):
        recheck(seen, placed=_late(seen['placed'], lambda p: p.event is blip, r))

    def triangle(area: float):
        return pp.make_extended_trapezoid('y', times=np.array([0.0, r, 2 * r]),
                                          amplitudes=np.array([0.0, area / r, 0.0]), system=opts)

    moved = _gradients(sc.LogicBlock('m').add(1.3e-3 - 2 * r, triangle(0.05))
                       .add(1.3e-3, triangle(-0.05)), opts)
    with pytest.raises(sc.CompilerContractError, match='axis y m0'):
        recheck(seen, placed=[*seen['placed'], *moved])


# ------------------------------------------------------------------------ registration precision
@pytest.mark.parametrize('t0', [0.0, 10.0])
def test_six_digit_deduplication_is_within_the_registration_tolerance(t0: float) -> None:
    """
    Two trapezoids at 100000.0 and 100000.4 Hz/m share PyPulseq's ``'%.6g'`` key.

    The second is stored as the first, which moves its area by 4.4e-4 1/m in memory.  That is
    PyPulseq's precision, so it is accepted -- and it is still what the sequence plays.
    """
    opts = _opts()
    first = pp.make_trapezoid('y', amplitude=100000.0, flat_time=1e-3, rise_time=1e-4, system=opts)
    second = pp.make_trapezoid('y', amplitude=100000.4, flat_time=1e-3, rise_time=1e-4,
                               system=opts)
    tree = sc.LogicBlock('d').add(t0, first).add(t0 + 2e-3, sc.barrier()).add(t0 + 3e-3, second)
    seq = sc.compile(tree, opts, name='dedup')
    stored = [seq.get_block(i).gy for i in sorted(seq.block_events)
              if getattr(seq.get_block(i), 'gy', None) is not None]
    assert [float(g.amplitude) for g in stored] == [100000.0, 100000.0]


def test_shapes_sharing_a_six_digit_key_are_within_the_registration_tolerance() -> None:
    """
    Two compressible shapes that print alike under ``'%.6g'``, one with a tail the other lacks.

    The tail is 3e-7 of the peak.  It shifts two of the stored differences by less than their sixth
    significant digit, so PyPulseq stores the second shape as the first and plays it without its
    100 ms tail: 7e-3 1/m of area, in memory.  Half a quantum per sample and the key's pointwise
    allowance together do not cover that; the key's allowance for a coarse shape difference,
    carried into every later sample, does.
    """
    opts = _suite()
    times = np.array([0.0, 0.1, 0.2, 20.0, 40.0, 60.0, 80.0, 99.8, 99.9, 100.0]) * 1e-3
    eps, c = 3e-7, 0.5                              # eps - 1 prints as -1, and c - eps as 0.5
    first = pp.make_extended_trapezoid('x', times=times, amplitudes=234567.1 * np.array(
        [0, 1, 0, 0, 0, 0, 0, 0, c, 0]), system=opts)
    second = pp.make_extended_trapezoid('x', times=times, amplitudes=234567.4 * np.array(
        [0, 1, eps, eps, eps, eps, eps, eps, c, 0]), system=opts)
    tree = sc.LogicBlock('key').add(0.0, first).add(0.101, sc.barrier()).add(0.2, second)
    seq, seen = compile_capturing(tree, opts)
    handed, stored = seen['handed'], seen['stored']()
    row = int(np.flatnonzero(handed.error0[:, 0])[-1])
    moved = float(stored.m0[row, 0] - handed.m0[row, 0])
    assert moved == pytest.approx(-7.06e-3, rel=0.01)

    played = knots_of(seq.get_block(row + 1).gx, 0.0)[1]
    assert played[2] == 0.0, 'the second shape plays the first one\'s tail'
    times, amps = knots_of(second, 0.0)
    traversed = float(np.sum(0.5 * np.diff(times) * (np.abs(amps[1:]) + np.abs(amps[:-1]))))
    assert abs(moved) > 1.5e-7 * float(np.max(amps)) * 0.1 + 1e-5 * (1 + 1e-5) * traversed


def test_the_compression_edge_case_is_inside_the_registration_tolerance() -> None:
    """
    A shape whose first stored sample lies exactly on a half quantum.

    Rounding half away from zero then carries one more quantum into every later sample, so the
    stored shape sits a quantum above the handed one wherever it is zero.  On this low-duty shape
    that is beyond what half a quantum per sample plus the key allowance would cover, and inside
    the stated bound of a quantum and a half.
    """
    opts = _suite()
    peak = 2.0**20                                    # a power of two keeps g / peak exact
    ramp = np.arange(1, 21) / 20.0
    samples = np.concatenate([[2.5 * 1e-7], np.zeros(5000), ramp, np.ones(10), ramp[::-1],
                              np.zeros(5000)]) * peak
    assert samples[0] / peak / 1e-7 == 2.5, 'the first sample must sit exactly on a half quantum'
    grad = pp.make_arbitrary_grad('x', waveform=samples, first=0.0, last=0.0, system=opts)
    _, seen = compile_capturing(sc.LogicBlock('edge').add(0.0, grad), opts)
    handed, stored = seen['handed'], seen['stored']()
    span = samples.size * float(opts.grad_raster_time)
    moved = float(stored.m0[0, 0] - handed.m0[0, 0])
    assert moved == pytest.approx(1e-7 * peak * span, rel=0.02)

    times, amps = knots_of(grad, 0.0)
    traversed = float(np.sum(0.5 * np.diff(times) * (np.abs(amps[1:]) + np.abs(amps[:-1]))))
    half_quantum = 0.5e-7 * peak * span + 1e-5 * (1 + 1e-5) * traversed
    assert moved > half_quantum, 'the case must need more than half a quantum per sample'
    assert moved < float(handed.storage0[0, 0])


# ----------------------------------------------------------------------------- the refusal itself
def test_a_refusal_names_the_stage_block_time_difference_and_tolerance() -> None:
    opts = _opts()
    tree = sc.LogicBlock('b').add(10.0, _holder(opts, 3e-3)).add(10.001, _blip(opts))
    _, seen = compile_capturing(tree, opts)
    with pytest.raises(sc.CompilerContractError) as err:
        recheck(seen, placed=_late(seen['placed'], _is_blip, float(opts.grad_raster_time)))
    line = next(line for line in str(err.value).splitlines() if 'axis y m1' in line)
    assert 'PyPulseq block 2 at 10000.000 ms' in line
    assert 'compiled m1' in line and 'about the block start' in line
    assert 'by -4.5e-05' in line and 'beyond the compile tolerance' in line


def test_warnings_are_unchanged_by_the_check() -> None:
    """The check adds no warning of its own: a clean, unmerged compile stays silent."""
    opts = _opts()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        sc.compile(sc.LogicBlock('q').add(5.0, _shape(opts, _compressible(1, 2500))), opts)
    assert not [w for w in caught if issubclass(w.category, sc.SeqCraftWarning)]
