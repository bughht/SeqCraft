"""
Boundary selection: the contract, and the cost of computing it.

The functional contract is short -- at most one RF and one ADC per block, nothing cut that
cannot be cut -- but the *cost* of arriving at it used to be quadratic on the one case the
fallback path exists for, so that gets a test too.
"""

from __future__ import annotations

import statistics
import time

import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.compiler.boundaries import find_boundaries
from seqcraft.compiler.placement import place_events


def _epi_tree(opts, n_echo: int) -> sc.LogicBlock:
    """
    One continuous readout gradient spanning `n_echo` ADC windows.

    The shape that forces the gap-midpoint fallback on every echo: the gradient covers both
    natural candidates (the end of one ADC's reservation and the start of the next), so neither
    is acceptable and a midpoint has to be invented each time.
    """
    lb = sc.LogicBlock('epi')
    adc = pp.make_adc(num_samples=64, dwell=4e-6, system=opts)
    period = 500e-6
    lb.add(0.0, pp.make_trapezoid('x', amplitude=0.2 * opts.max_grad,
                                  duration=n_echo * period, system=opts))
    for i in range(n_echo):
        lb.add(100e-6 + i * period, adc)
    return lb


def _time_boundaries(opts, n_echo: int) -> float:
    placed = place_events(_epi_tree(opts, n_echo), opts)
    raster = sc.Raster(opts.block_duration_raster)
    total = raster.ceil(max(p.res_end for p in placed))
    max_block = float(opts.block_duration_raster) * 2**24
    best = float('inf')
    for _ in range(3):                          # best of three; we care about the floor, not noise
        t0 = time.perf_counter()
        find_boundaries(placed, total, raster, max_block)
        best = min(best, time.perf_counter() - t0)
    return best


def _scaling_ratio(opts) -> tuple[float, float, float]:
    """One paired measurement of what quadrupling the echo count costs."""
    small = _time_boundaries(opts, 400)
    large = _time_boundaries(opts, 1600)
    return large / max(small, 1e-9), small, large


def test_boundary_selection_is_not_quadratic(opts) -> None:
    """
    Quadrupling the echo count must not multiply the work by ~16.

    Measured across a 4x step rather than 2x so the linear and quadratic predictions (4x versus
    16x) are far enough apart to assert on without flaking.  Before the fix this path re-sorted
    the whole mark set once per echo: 253 ms at 3200 echoes against 17 ms after, and the ratio
    per doubling was 3.7 rather than 2.

    **The median of five paired ratios, not one.**  Each size is already a best of three, but that
    only removes interference that hits *some* of the three repetitions: the quantities here are
    about 1.3 ms and 5.5 ms on the wall clock, so a single ratio is still dominated by the machine
    descheduling the large measurement and not the small one.  On a quiet machine that is
    invisible -- the ratio sits at 4.2 with no spread worth mentioning -- but CI runs several
    xdist workers on a shared runner, and under deliberate CPU contention a single ratio reached
    12.3 and crossed this threshold in 3 of 20 attempts.  The signal survives that; the *tail* is
    what does not.  Taking the median of five held the worst of 20 attempts at 4.96 under the same
    load, which was the smallest aggregation tested that removed the failures.

    The obvious alternative, a CPU-time clock, is not available here: on the Windows runner
    ``time.process_time`` reports 100 ns resolution and in fact advances only once per 15.625 ms
    scheduler quantum, so a 1.3 ms measurement reads as either 0 or 15.625 ms.  ``thread_time``
    is the same.  Measuring for longer does not rescue it, and measuring for longer on the wall
    clock makes things worse rather than better -- an 800/3200 pair failed 7 of 20 under the load
    above, because a longer window accumulates more interference, not less.

    What the threshold still buys: against the quadratic implementation this replaced, the same
    statistic measured at least 22.1 under that contention, so the two live on opposite sides of
    8.0 with a factor of four between them.
    """
    measured = [_scaling_ratio(opts) for _ in range(5)]
    ratio = statistics.median(r for r, _, _ in measured)
    assert ratio < 8.0, (
        f'boundary selection scaled by {ratio:.1f}x (median of five) for a 4x larger EPI train; '
        f'linear predicts ~4x, quadratic ~16x.  Measured '
        + ', '.join(f'{r:.1f}x ({lo * 1e3:.1f} -> {hi * 1e3:.1f} ms)' for r, lo, hi in measured)
    )


def test_every_reservation_gap_gets_a_boundary(opts) -> None:
    """
    The invariant the midpoint fallback exists to guarantee.

    With it broken, two ADCs share a block and pypulseq raises 'Multiple ADC events' on a block
    index far from anything that looks wrong.
    """
    n = 40
    tree = _epi_tree(opts, n)
    out = sc.compile(tree, opts)
    per_block = []
    for index in sorted(out.block_events):
        block = out.get_block(index)
        per_block.append(
            (1 if getattr(block, 'adc', None) is not None else 0,
             1 if getattr(block, 'rf', None) is not None else 0)
        )
    assert sum(a for a, _ in per_block) == n, 'every ADC must survive'
    assert all(a <= 1 for a, _ in per_block), 'no block may hold two ADCs'
    assert all(r <= 1 for _, r in per_block), 'no block may hold two RF pulses'


def test_a_long_empty_stretch_is_subdivided_to_fit_the_duration_field(opts) -> None:
    """pulseq stores a block duration in a fixed-width field, so one long delay must split."""
    out = sc.compile(sc.LogicBlock('t').add(0.0, pp.make_delay(1.0)), opts)
    assert out.duration()[0] == pytest.approx(1.0)
