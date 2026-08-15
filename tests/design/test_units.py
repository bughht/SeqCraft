"""
Unit conversion: exactness, symmetry, and the errors.

``convert`` is the one place a factor is written down, so these tests check the properties that
make that safe -- every pair round-trips, decimal conversions are exact rather than nearly exact,
and a wrong unit is a loud error rather than a plausible number.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import seqcraft as sc
from seqcraft.design.units import GAMMA_1H, convert, dimension_of, dimensions, known_units

# ------------------------------------------------------------------------------ known values
# Each is an independently quoted equivalence, not a number the code produced.
KNOWN = [
    (1.0, 's', 'ms', 1e3),
    (1.0, 's', 'us', 1e6),
    (1.0, 's', '10us', 1e5),
    (1.0, 'min', 's', 60.0),
    (220.0, 'mm', 'm', 0.22),
    (1.0, 'cm', 'mm', 10.0),
    (180.0, 'deg', 'rad', math.pi),
    (1.0, 'turn', 'deg', 360.0),
    (1.0, 'G/cm', 'mT/m', 10.0),                    # 1 G/cm == 10 mT/m
    (1.0, 'T/m', 'mT/m', 1e3),
    (1.0, 'mT/m/ms', 'T/m/s', 1.0),                 # the same slew rate, spelled two ways
    (1.0, 'G/cm/ms', 'T/m/s', 10.0),
    (1000.0, 's/mm^2', 's/m^2', 1e9),
    (1.0, 'ms/um^2', 's/mm^2', 1e3),
    (1.0, '1/mm', '1/m', 1e3),
    (1.0, 'ppm', '1', 1e-6),
    (1.0, '%', 'ppm', 1e4),
    (1.0, 'rad/s', 'Hz', 1.0 / (2.0 * math.pi)),
]


@pytest.mark.parametrize(('value', 'src', 'dst', 'want'), KNOWN)
def test_known_equivalences(value: float, src: str, dst: str, want: float) -> None:
    assert convert(value, src, dst) == pytest.approx(want, rel=1e-12)


def test_gamma_dependent_known_values() -> None:
    """The three conversions that need a gyromagnetic ratio, at the proton value."""
    assert convert(40.0, 'mT/m', 'Hz/m') == pytest.approx(40e-3 * GAMMA_1H)
    assert convert(200.0, 'T/m/s', 'Hz/m/s') == pytest.approx(200.0 * GAMMA_1H)
    assert convert(1.0, 'T', 'Hz') == pytest.approx(GAMMA_1H)
    # 1 T = 42.576 MHz for 1H, the number every MR textbook quotes.
    assert convert(1.0, 'T', 'MHz') == pytest.approx(42.576, rel=1e-9)


def test_gamma_cancels_within_the_tesla_family() -> None:
    """
    mT/m -> G/cm needs no gamma, so passing an absurd one must not change the answer.

    The implementation tracks the power of gamma per unit rather than routing through Hz/m, which
    is what keeps this exact instead of merely close.
    """
    assert convert(40.0, 'mT/m', 'G/cm', gamma=1.0) == convert(40.0, 'mT/m', 'G/cm', gamma=1e9)
    assert convert(40.0, 'mT/m', 'G/cm') == 4.0


# ---------------------------------------------------------------------------------- exactness
@pytest.mark.parametrize('unit', ['ms', 'us', 'ns', '10us'])
@pytest.mark.parametrize('value', [1.0, 4.2, 1500.0, 3200.0, 4200.0, 0.0001])
def test_time_round_trips_bit_exactly(unit: str, value: float) -> None:
    """A parameter converted to seconds and back must be the number the user typed."""
    assert convert(convert(value, unit, 's'), 's', unit) == value


def test_decimal_conversions_are_exact_not_approximate() -> None:
    """
    The reason the scale is a Fraction: a single correctly-rounded operation, not two.

    ``4200 * 1e-6`` is 0.004200000000000001; ``4200 / 1e6`` is 0.0042.  Both are "close", but only
    the second compares equal to the raster value a module then quantises against.
    """
    assert convert(4200, 'us', 's') == 0.0042
    assert convert(1.5, 'ms', 'us') == 1500.0
    assert convert(0.0121, 's', '10us') == 1210.0
    assert convert(3.2, 'ms', 's') == 0.0032


@pytest.mark.parametrize(('src', 'dst'), [(row[1], row[2]) for row in KNOWN])
def test_every_known_pair_round_trips(src: str, dst: str) -> None:
    back = convert(convert(7.0, src, dst), dst, src)
    assert back == pytest.approx(7.0, rel=1e-12)


def test_all_pairs_within_a_dimension_round_trip() -> None:
    """Exhaustive over the table: no unit is a one-way street."""
    for dimension in dimensions():
        units = known_units(dimension)
        for src in units:
            for dst in units:
                back = convert(convert(3.0, src, dst), dst, src)
                assert back == pytest.approx(3.0, rel=1e-9), f'{src} -> {dst} -> {src}'


def test_canonical_default_matches_the_explicit_call() -> None:
    for unit in known_units():
        assert convert(2.0, unit) == convert(2.0, unit, None)


# ------------------------------------------------------------------------------- ppm and f0
def test_ppm_to_hz_needs_a_larmor_frequency() -> None:
    with pytest.raises(sc.ConfigurationError, match='Larmor'):
        convert(-3.4, 'ppm', 'Hz')


def test_ppm_to_hz_matches_the_hand_calculation() -> None:
    f0 = convert(3.0, 'T', 'Hz')
    assert convert(-3.4, 'ppm', 'Hz', f0=f0) == pytest.approx(-3.4e-6 * GAMMA_1H * 3.0)
    assert convert(-434.0, 'Hz', 'ppm', f0=f0) == pytest.approx(-434.0 / f0 * 1e6)


def test_the_scanners_own_gamma_is_passed_not_assumed(opts) -> None:
    """
    ``gamma`` and ``f0`` are arguments rather than globals, so a non-proton scanner cannot
    silently pick up the proton values.  There is no scanner object to fill them in -- the
    ``Opts`` carries them, and the call site passes them.
    """
    assert convert(-3.4, 'ppm', 'Hz', f0=opts.gamma * opts.B0) == pytest.approx(
        -3.4e-6 * opts.gamma * opts.B0
    )
    assert convert(opts.max_grad, 'Hz/m', 'mT/m', gamma=opts.gamma) == pytest.approx(40.0)


# ---------------------------------------------------------------------------------- arrays
def test_arrays_convert_elementwise() -> None:
    """Waveforms are one call, which is why there is no Quantity type in the way."""
    got = convert(np.array([10.0, 20.0, 30.0]), 'mT/m', 'G/cm')
    assert np.allclose(got, [1.0, 2.0, 3.0])
    assert isinstance(got, np.ndarray)


# ----------------------------------------------------------------------------------- errors
def test_unknown_unit_names_the_closest_match() -> None:
    with pytest.raises(sc.ConfigurationError, match="did you mean.*'mT/m'"):
        convert(1.0, 'mt/m', 'Hz/m')


def test_crossing_unrelated_dimensions_is_an_error_that_lists_the_alternatives() -> None:
    with pytest.raises(sc.ConfigurationError, match='cannot convert') as excinfo:
        convert(1.0, 'ms', 'mm')
    assert 'time units' in str(excinfo.value)


def test_dimension_lookup() -> None:
    assert dimension_of('mT/m') == 'gradient'
    assert dimension_of('uT') == 'frequency'          # B1 is carried in hertz
    assert dimension_of('s/mm^2') == 'bvalue'
    with pytest.raises(sc.ConfigurationError, match='unknown dimension'):
        known_units('nonsense')
