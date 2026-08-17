"""
``PhaseEncode`` -- the k mapping, the parity, and the invariant the module exists for.

The invariant is *equal duration for every line*.  It is not a property of the base class and
nothing else in the package checks it, so it is asserted here: the caller's placement arithmetic
depends on it, and a module that designed each blip at its own minimum would be shorter on
average and would move the echo line by line.
"""

from __future__ import annotations

import pytest

import seqcraft as sc


@pytest.mark.parametrize('matrix', [64, 65])
def test_the_area_follows_the_k_mapping_for_either_parity(opts, matrix) -> None:
    """
    ``k = (line - matrix // 2) * 1000 / fov_mm``, with no special case for an odd matrix.

    64 spans −32 … +31 and 65 spans −32 … +32, which is what zero-based indexing gives for free
    and what a signed ``-matrix/2 … matrix/2`` convention has to special-case.
    """
    pe = sc.modules.PhaseEncode(opts=opts, fov_mm=250.0, matrix=matrix)
    dk = 1e3 / 250.0

    for line in (0, 1, matrix // 2, matrix - 1):
        area = float(pe(line=line).nodes[0].item.area)
        assert area == pytest.approx((line - matrix // 2) * dk, abs=1e-9)


def test_the_centre_line_encodes_nothing(opts) -> None:
    """``k(center_line) == 0``, and ``center_line`` is ``matrix // 2`` for either parity."""
    for matrix in (64, 65):
        pe = sc.modules.PhaseEncode(opts=opts, fov_mm=250.0, matrix=matrix)
        assert pe.center_line == matrix // 2
        assert pe.k_per_m(pe.center_line) == 0.0
        assert float(pe(line=pe.center_line).nodes[0].item.area) == 0.0


def test_every_line_takes_the_same_time(opts) -> None:
    """The invariant this module exists for."""
    pe = sc.modules.PhaseEncode(opts=opts, fov_mm=250.0, matrix=64)

    durations = {pe(line=line).duration for line in range(64)}

    assert len(durations) == 1, f'{len(durations)} different blip lengths: {sorted(durations)}'


def test_a_rewinder_is_the_same_blip_backwards(opts) -> None:
    """Same line, same duration, opposite area -- and both go through ``__call__``."""
    pe = sc.modules.PhaseEncode(opts=opts, fov_mm=250.0, matrix=64)

    forward, back = pe(line=17), pe(line=17, rewind=True)

    assert float(back.nodes[0].item.area) == pytest.approx(-float(forward.nodes[0].item.area))
    assert back.duration == forward.duration
    assert back.tag == forward.tag == 'PhaseEncode', 'both were finalized, so both are tagged'


@pytest.mark.parametrize('line', [-1, 64, 100])
def test_a_line_outside_the_matrix_raises(opts, line) -> None:
    """No clamping and no wrap-around: an index nobody meant is a bug, not a request."""
    pe = sc.modules.PhaseEncode(opts=opts, fov_mm=250.0, matrix=64)

    with pytest.raises(sc.ConfigurationError, match='0 ... 63'):
        pe(line=line)


def test_line_has_no_default(opts) -> None:
    """There is no sensible default line, so there is none."""
    pe = sc.modules.PhaseEncode(opts=opts, fov_mm=250.0, matrix=64)

    with pytest.raises(TypeError, match='line'):
        pe()


def test_a_requested_duration_stretches_every_line(opts) -> None:
    """The winder override: longer blips, still all equal, still the same areas."""
    short = sc.modules.PhaseEncode(opts=opts, fov_mm=250.0, matrix=64)
    stretched = sc.modules.PhaseEncode(
        opts=opts, fov_mm=250.0, matrix=64, duration_s=short.min_duration_s + 200e-6,
    )

    assert stretched(line=0).duration == pytest.approx(short.min_duration_s + 200e-6)
    assert {stretched(line=line).duration for line in range(64)} == {stretched(line=0).duration}
    assert float(stretched(line=40).nodes[0].item.area) == pytest.approx(
        float(short(line=40).nodes[0].item.area),
    )


def test_min_duration_s_is_the_minimum_and_stays_it(opts) -> None:
    """
    ``min_duration_s`` reports the shortest legal blip whether or not one was requested.

    Which is what lets the caller of the winder coupling ask both leaves for their minimum and
    hand the maximum back, without the answer moving underneath it.
    """
    short = sc.modules.PhaseEncode(opts=opts, fov_mm=250.0, matrix=64)
    stretched = sc.modules.PhaseEncode(
        opts=opts, fov_mm=250.0, matrix=64, duration_s=short.min_duration_s + 200e-6,
    )

    assert stretched.min_duration_s == short.min_duration_s


def test_a_duration_below_the_minimum_raises_naming_it(opts) -> None:
    pe = sc.modules.PhaseEncode(opts=opts, fov_mm=250.0, matrix=64)

    with pytest.raises(sc.ConfigurationError, match='min_duration_s'):
        sc.modules.PhaseEncode(
            opts=opts, fov_mm=250.0, matrix=64, duration_s=pe.min_duration_s / 2,
        )


def test_handing_back_its_own_minimum_is_accepted(opts) -> None:
    """
    The winder coupling's exact call, which float arithmetic makes non-trivial.

    ``max(ro.prephaser_duration_s, pe.min_duration_s)`` can land one ulp below where it came
    from once it is snapped onto the raster, so the comparison is against a nanosecond tolerance
    rather than exact.  Getting this wrong makes the composite refuse its own arithmetic.
    """
    pe = sc.modules.PhaseEncode(opts=opts, fov_mm=250.0, matrix=64)

    again = sc.modules.PhaseEncode(
        opts=opts, fov_mm=250.0, matrix=64, duration_s=pe.min_duration_s,
    )

    assert again(line=0).duration == pytest.approx(pe.min_duration_s)


def test_it_is_pure_and_compiles_alone(opts, component_checks) -> None:
    """
    ``line=17`` rather than the centre: a zero-amplitude blip passes a limit check for the
    wrong reason, and would not notice a scaling bug either.
    """
    component_checks.all(sc.modules.PhaseEncode(opts=opts, fov_mm=250.0, matrix=64), line=17)
