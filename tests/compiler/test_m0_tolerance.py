"""
Registration rounding of a compressible shape is accepted; area actually lost is not.

PyPulseq rounds a compressible shape's peak-normalised samples to 1e-7 when it registers them,
which moves an ordinary waveform's area by a few 1e-6 to a few 1e-5 1/m.  The check once compared
the whole sequence's m0 against a tolerance scaled by each event's ``|net area|``.  That collapsed
to 1e-6 1/m for a gradient that changes sign and nets zero, so legal input was refused -- and which
of two identical sequences was refused depended on how the tree packaged the waveform.

m0 and m1 are now compared per compiled block, in two stages (``test_local_moments.py``).
Registration is compared with what emission handed PyPulseq, within PyPulseq's storage precision:
at most 1.5 quanta of the peak per stored sample over the shape's span, plus its six-digit key.
That bound has no duty condition, so the dense low-duty shape at the bottom of this file, which
exceeded the earlier traversed-area tolerance, compiles as well.  The compile stage still refuses
area that is actually lost, and a piece played one raster late, at the start of a sequence or
ten seconds into it.

Nothing here builds a module.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pypulseq as pp
import pytest
from pypulseq.compress_shape import compress_shape
from recheck import compile_capturing, recheck

import seqcraft as sc
from seqcraft.compiler.placement import place_events
from seqcraft.compiler.verification import _sequence_moments
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


def _late(placed, first: float, dt: float) -> list:
    """The placed events with the gradient that starts at `first` moved `dt` later."""
    return [dataclasses.replace(p, node_t=p.node_t + dt, start=p.start + dt, end=p.end + dt,
                                res_start=p.res_start + dt, res_end=p.res_end + dt)
            if p.kind in ('trap', 'grad') and p.node_t == first else p
            for p in placed]


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
    What the registration tolerance has to cover, measured rather than assumed.

    The stored area is the area of the peak-normalised samples rounded half away from zero to
    1e-7 -- reproduced here independently of PyPulseq's code.  It exceeds the old floor, which is
    why these were refused.  It is a few percent of the worst case, half a quantum per sample, as
    independent rounding errors should be, and inside the registration tolerance by two orders.
    """
    opts = _opts()
    grad = _gradient(_zero_net_section(1, 2500), 0.9, opts)
    times, amps = knots_of(grad, 0.0)
    peak = float(np.max(np.abs(amps)))
    scaled = amps / peak / QUANTUM
    rounded = np.sign(scaled) * np.floor(np.abs(scaled) + 0.5) * QUANTUM * peak

    seq, seen = compile_capturing(sc.LogicBlock('m').add(0.0, grad), opts, check=False)
    moved = _sequence_moments(seq, 0)['y'] - pwl_moment(times, amps, 0)

    # Each side sums ~1e4 1/m of partial areas, so float alone is a few 1e-12; 1e-10 is five
    # orders below the effect.
    assert moved == pytest.approx(pwl_moment(times, rounded, 0) - pwl_moment(times, amps, 0),
                                  abs=1e-10)
    assert abs(moved) > 1e-6, 'it must exceed the old floor, or this is not the regression'
    worst = 0.5 * QUANTUM * peak * float(times[-1] - times[0])
    assert abs(moved) < 0.1 * worst
    assert abs(moved) < 1e-2 * float(seen['handed'].storage0[0, 1])


# ------------------------------------------------------------------------------------- packaging
@pytest.mark.parametrize('seed', [1, 4])
def test_four_divisions_of_one_waveform_get_one_verdict(seed: int) -> None:
    """
    The packaging dependence these four cases exposed, and what removes it for them.

    One event; two events cut at an exact zero, a division in time that leaves the scale
    unchanged; two overlapping pieces, whose scale can only be larger; and one event cut by a
    block barrier.  The emitted waveform is the same each time, and so is the verdict.  Under the
    sum of each event's |net|, the first and last were refused and the middle two passed.  This
    is not a claim about every possible decomposition.
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
    Accepting registration rounding must not hide what m0 exists to catch.

    The compiled sequence is checked against a tree asking for more than it emitted: a whole lobe
    more, and then only a sliver -- 0.05 1/m, one percent of a k-space step at 200 mm.  The
    message reports the stage, the difference and the tolerance, and leaves the cause to the
    reader rather than guessing a split or a merge.
    """
    opts = _opts()
    tree = _two_events(opts)
    _, seen = compile_capturing(tree, opts, name='clean')     # clean to begin with
    placed = seen['placed']

    lobe = min((p for p in placed if p.kind in ('trap', 'grad')), key=lambda p: p.node_t)
    with pytest.raises(sc.CompilerContractError) as err:
        recheck(seen, placed=[*placed, lobe])
    text = str(err.value)
    assert 'compiled m0' in text and 'beyond the compile tolerance' in text
    assert 'split or a merge' not in text

    sliver = pp.make_trapezoid('y', area=0.05, duration=200e-6, system=opts)
    with pytest.raises(sc.CompilerContractError, match='compiled m0'):
        recheck(seen, placed=place_events(tree.add(0.0, sliver), opts))


def test_m1_still_refuses_a_piece_played_one_raster_late() -> None:
    """
    Accepting registration rounding must not hide a time shift either.

    The lobe ends where a block does, so one raster late also carries some of its ramp into the
    next block, which m0 may report too; m1 is the moment that sees the shift itself.
    """
    opts = _opts()
    tree = _two_events(opts)
    _, seen = compile_capturing(tree, opts, name='clean')
    first = min(p.node_t for p in seen['placed'] if p.kind in ('trap', 'grad'))

    with pytest.raises(sc.CompilerContractError, match='compiled m1'):
        recheck(seen, placed=_late(seen['placed'], first, float(opts.grad_raster_time)))


# ------------------------------------------------- the same controls, far from the sequence start
def _two_lobes_at(t0: float, opts: pp.Opts) -> sc.LogicBlock:
    """
    Two trapezoids of non-zero area on one axis, starting at `t0`.

    Trapezoids are stored as parameters, so registration rounds nothing, and only what a control
    puts there changes with `t0`.  A translated zero-net waveform would be the wrong probe: it moves
    m1 by `t0` times its net area, which is zero.
    """
    first = pp.make_trapezoid('y', area=500.0, duration=2e-3, system=opts)
    second = pp.make_trapezoid('y', area=-300.0, duration=2e-3, system=opts)
    return sc.LogicBlock('late').add(t0, first).add(t0 + 3e-3, second)


@pytest.mark.parametrize('t0', [0.0, 1.0, 10.0])
def test_lost_area_is_refused_wherever_the_axis_sits(t0: float) -> None:
    opts = _opts()
    _, seen = compile_capturing(_two_lobes_at(t0, opts), opts, name='clean')
    placed = seen['placed']
    lobe = min((p for p in placed if p.kind in ('trap', 'grad')), key=lambda p: p.node_t)

    with pytest.raises(sc.CompilerContractError, match='compiled m0'):
        recheck(seen, placed=[*placed, lobe])


@pytest.mark.parametrize('t0', [0.0, 1.0, 10.0])
def test_a_lobe_one_raster_late_is_refused_wherever_the_axis_sits(t0: float) -> None:
    """
    Only one side of the comparison moves, and the piece that moves has area.

    The 500 1/m lobe played one raster late changes its block's m1 by about 500 1/m * 10 us,
    however far into the sequence it sits.
    """
    opts = _opts()
    _, seen = compile_capturing(_two_lobes_at(t0, opts), opts, name='clean')
    first = min(p.node_t for p in seen['placed'] if p.kind in ('trap', 'grad'))

    with pytest.raises(sc.CompilerContractError, match='compiled m1'):
        recheck(seen, placed=_late(seen['placed'], first, float(opts.grad_raster_time)))


# ------------------------------------------------------------------------- low duty, now covered
def test_a_dense_low_duty_shape_with_systematic_rounding_compiles() -> None:
    """
    A 0.5 ms spike cancelled by a long weak constant tail, in one densely sampled shape.

    Mean |g| is 2 % of the peak.  The tail's level is chosen so its samples round almost half a
    quantum each way the same way, so the stored area moves by nearly the worst case.  1e-6 of
    the traversed area did not cover that; a bound per stored sample over the shape's span does,
    at any duty.
    """
    opts = _opts()
    n = max(range(3000, 3040), key=lambda m: abs(30.0 / m / QUANTUM - round(30.0 / m / QUANTUM)))
    ramp = np.arange(21) / 20.0
    shape = np.concatenate([ramp, np.ones(9), ramp[::-1], -np.full(n, 30.0 / n), [0.0]])
    grad = _extended(np.arange(shape.size) * opts.grad_raster_time,
                     20e-3 * opts.gamma * shape, opts)
    sc.compile(sc.LogicBlock('low_duty').add(0.0, grad), opts, name='low_duty')
