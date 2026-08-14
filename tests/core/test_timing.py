"""
Exact raster arithmetic.

Every test here is a float trap that was observed rather than imagined: an exact multiple pushed up
a whole raster, ``n * raster`` landing a femtosecond low, a duration that drifted off the raster
after being summed.  The point of :class:`~seqcraft.Raster` is that none of them can happen at a
call site, because the call site never sees the arithmetic.
"""

from __future__ import annotations

import pytest

import seqcraft as sc
from seqcraft.core.timing import (
    TICKS_PER_SECOND,
    Raster,
    exact_diff,
    exact_sum,
    from_ticks,
    to_ticks,
)

#: Rasters from four vendors plus two invented ones.  Nothing in timing.py may assume 10 us.
RASTERS = [
    10e-6,      # Siemens gradient / block
    1e-6,       # Siemens RF
    100e-9,     # Siemens ADC
    4e-6,       # GE gradient
    2e-6,       # GE RF
    6.4e-6,     # Philips
    2.5e-6,     # invented: a finer future raster
    12.5e-6,    # invented: a coarser one, not a whole number of microseconds
]


# ------------------------------------------------------------------------------ the four ops
@pytest.mark.parametrize('dt', RASTERS)
def test_exact_multiples_are_left_alone(dt: float) -> None:
    """
    The trap that motivates the whole module: ``ceil(1.5e-3 / 1e-5) * 1e-5`` moves 1.5 ms to 1.51.

    ``1.5e-3 / 1e-5`` is 149.99999999999997 in float64, so a naive ceiling adds a whole raster to a
    duration that was already legal.
    """
    r = Raster(dt)
    for n in (1, 2, 7, 150, 1000, 99_991):
        t = r.at(n)
        assert r.ceil(t) == t
        assert r.floor(t) == t
        assert r.nearest(t) == t
        assert r.holds(t)
        assert r.count(t) == n


@pytest.mark.parametrize('dt', RASTERS)
def test_ceil_floor_bracket_the_value(dt: float) -> None:
    r = Raster(dt)
    for t in (1e-6, 3.33e-4, 1.2345e-3, 0.5):
        lo, hi = r.floor(t), r.ceil(t)
        assert lo <= t <= hi
        assert hi - lo == pytest.approx(dt) or lo == hi
        assert r.holds(lo)
        assert r.holds(hi)


def test_at_is_exact_where_multiplication_is_not() -> None:
    """``250 * 1e-7`` is 2.4999999999999998e-05, which is how noise reached ADC dwell times."""
    assert 250 * 1e-7 != 2.5e-05
    assert Raster(1e-7).at(250) == 2.5e-05
    assert Raster(1e-5).at(150) == 0.0015


def test_nearest_rounds_half_up_and_is_symmetric_about_zero() -> None:
    r = Raster(1e-5)
    assert r.nearest(1.4e-5) == 1e-05
    assert r.nearest(1.6e-5) == 2e-05
    assert r.nearest(-1.4e-5) == -1e-05      # nearest, not "away from zero"
    assert r.nearest(-1.6e-5) == -2e-05
    assert r.ceil(-1.23e-5) == -1e-05
    assert r.floor(-1.23e-5) == -2e-05


def test_holds_tolerates_noise_finer_than_a_tick() -> None:
    r = Raster(1e-5)
    assert r.holds(1.5e-3 * (1 + 1e-15))
    assert not r.holds(1.5e-3 + 1e-6)


# --------------------------------------------------------------------------------- exactness
def test_exact_sum_stays_on_the_raster() -> None:
    """
    Plain addition of raster-aligned floats does not stay aligned, and the drift is not academic:
    a solved TE then compares unequal to the value it came from, and the next ceil adds 10 us.
    """
    terms = [1e-5] * 121
    r = Raster(1e-5)
    assert sum(terms) != 0.00121
    assert exact_sum(terms) == 0.00121
    assert r.holds(exact_sum(terms))


def test_exact_diff_matches_the_decimal_answer() -> None:
    """3.7 ms minus 1.9 ms is 1.8000000000000002 ms in float64: close, but not equal."""
    assert 3.7e-3 - 1.9e-3 != 1.8e-3
    assert exact_diff(3.7e-3, 1.9e-3) == 1.8e-3
    assert 0.05 - 0.0123 != 0.0377
    assert exact_diff(0.05, 0.0123) == 0.0377


def test_accumulating_a_tr_leaves_the_raster_entirely() -> None:
    """
    The measured consequence of ``t += tr`` in a loop, and the reason a start time is a closed form.

    Accumulating an exactly-legal 1.5 ms TR drifts off the 10 us raster after 9813 repetitions --
    14.7 s into the sequence -- and then more than half of all later start times are illegal.
    pypulseq rejects them individually, thousands of blocks after anything that looks wrong.
    """
    r = Raster(1e-5)
    running = 0.0
    off = 0
    for _ in range(20_000):
        running += 1.5e-3
        off += not r.holds(running)
    assert off > 10_000

    total = exact_sum([1.5e-3] * 20_000)
    assert total == 30.0
    assert r.holds(total)


def test_a_start_time_late_in_a_sequence_stays_on_its_own_raster() -> None:
    """
    40 TRs of 1.997 ms starting 3 s in: the shape that produced an RF delay of 129.9999999986 us.

    Accumulation lands 8e-15 s low -- inside a tick, so still *on* the raster, but no longer equal
    to the value any other route computes, which is how two supposedly-identical times diverge.
    """
    rf = Raster(1e-6)
    accumulated = 3.0
    for _ in range(40):
        accumulated += 1.997e-3
    exact = exact_sum([3.0, *[1.997e-3] * 40])

    assert accumulated != exact
    assert exact == 3.07988
    assert rf.holds(exact)
    assert rf.count(exact) == 3_079_880


def test_ticks_round_trip() -> None:
    assert to_ticks(2.5e-05) == 25_000_000
    assert from_ticks(25_000_000) == 2.5e-05
    assert TICKS_PER_SECOND == 10**12
    for t in (1e-9, 1e-6, 1.5e-3, 3.0, 3600.0):
        assert from_ticks(to_ticks(t)) == pytest.approx(t, rel=1e-15)


# ------------------------------------------------------------------------------- the errors
def test_count_refuses_to_round_silently() -> None:
    """
    An integer count is what makes durations summable, so an off-raster value is a hard error.

    Rounding here would hide the caller's mistake in an integer that then looks authoritative.
    """
    r = Raster(1e-5, 'block')
    with pytest.raises(sc.RasterError, match='is not a multiple'):
        r.count(1.505e-3)


def test_require_names_the_two_nearest_legal_values() -> None:
    r = Raster(1e-5, 'block')
    with pytest.raises(sc.RasterError) as excinfo:
        r.require(1.505e-3)
    message = str(excinfo.value)
    assert '1.500000 ms' in message
    assert '1.510000 ms' in message
    assert 'block raster' in message


@pytest.mark.parametrize('dt', [0.0, -1e-5, 1e-15])
def test_an_impossible_raster_is_rejected_at_construction(dt: float) -> None:
    with pytest.raises(sc.RasterError):
        Raster(dt)


# ------------------------------------------------------------------------------ value semantics
def test_rasters_compare_and_hash_by_interval() -> None:
    """A raster is a value, so two reads of the same scanner give the same raster."""
    assert Raster(1e-5, 'block') == Raster(1e-5, 'gradient')
    assert Raster(1e-5) != Raster(1e-6)
    assert len({Raster(1e-5), Raster(1e-5, 'x'), Raster(1e-6)}) == 2
    assert float(Raster(1e-5)) == 1e-5


def test_a_raster_cannot_be_mutated() -> None:
    r = Raster(1e-5, 'block')
    with pytest.raises(AttributeError):
        r.dt = 1e-6


def test_repr_reads_as_a_duration() -> None:
    assert repr(Raster(1e-5, 'block')) == 'Raster(block, 10 us)'
    assert repr(Raster(1e-7, 'ADC')) == 'Raster(ADC, 100 ns)'
    assert repr(Raster(6.4e-6)) == 'Raster(6.4 us)'


# ---------------------------------------------------------------------------- System wiring
def test_system_exposes_its_rasters_as_objects() -> None:
    system = sc.System.preset('prisma')
    assert system.block_raster == Raster(10e-6)
    assert system.grad_raster == Raster(10e-6)
    assert system.rf_raster == Raster(1e-6)
    assert system.adc_raster == Raster(100e-9)
    assert system.grad_raster.dt == 1e-5


def test_a_scanner_with_unusual_rasters_works_end_to_end() -> None:
    """
    The user's point: 10 us is today's value, not an assumption.  A 4 us / 2 us system compiles.
    """
    system = sc.System.from_limits(
        sc.Limits(50.0, 180.0), name='four_us', grad_raster_us=4.0, rf_raster_us=2.0,
        block_raster_us=4.0,
    )
    assert system.grad_raster == Raster(4e-6)
    exc = sc.modules.HardExcitation(system, flip_deg=90, duration_us=500)
    out = sc.compile(sc.LogicBlock('t').add(0.0, exc.build()), system)
    assert out.check().ok
    for duration in out.seq.block_durations.values():
        assert system.block_raster.holds(float(duration))
