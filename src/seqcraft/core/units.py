"""
Unit conversion: one function, both directions, every pair.

:func:`convert` takes the value, the unit it is in, and the unit you want -- the shape
``pypulseq.convert`` uses, generalised from gradients and slew rates to the nine dimensions a
sequence actually deals in.  Modules therefore never spell out a factor:

>>> convert(220, 'mm', 'm')
0.22
>>> convert(1.5, 'ms', 'us')
1500.0
>>> convert(40, 'mT/m', 'Hz/m')                        # at the proton gamma
1703040.0
>>> convert(1, 'G/cm', 'mT/m')                         # no gamma needed: same family
10.0
>>> round(convert(3.0, 'T', 'MHz'), 3)                 # Larmor frequency
127.728

The two-tier naming rule
------------------------
seqcraft converts **once**, at the boundary where a user-facing parameter enters a module:

===============================  ===============================================
tier                             convention
===============================  ===============================================
public parameters                researcher-natural units, unit **in the name**:
                                 ``fov_mm``, ``te_ms``, ``flip_deg``,
                                 ``max_grad_mT_m``, ``readout_duration_us``
internal / derived values        strict SI or pulseq-native, SI suffix: ``fov_m``,
                                 ``te_s``, ``flip_rad``, ``max_grad_Hz_per_m``
===============================  ===============================================

There is no ``Quantity`` type and no ``pint``.  A unit type would poison numpy interop for
waveform arrays and would still not prevent the actual observed bug class, which is a *plausible
number in the wrong unit* -- ``fov=220`` meaning millimetres reaching code that reads metres.
That is caught by the range bands in :mod:`seqcraft.core.validate`, not by a type.

Precision
---------
Each unit's scale is an exact :class:`~fractions.Fraction` wherever it is a decimal multiple, so
the ratio between two units is computed exactly and collapsed to a float **once**.  The result
carries a single rounding, which is the best float64 can do:

>>> convert(convert(3.2e-3 / 128, 's', 'ns'), 'ns', 's')       # exact round trip
2.5e-05

Exact *raster* arithmetic -- summing durations, quantising onto a 10 us grid -- is a different
job and lives in :mod:`seqcraft.core.timing`.

Gamma, and other system parameters
----------------------------------
Pairs that cross the tesla/hertz divide (``mT/m`` <-> ``Hz/m``, ``uT`` <-> ``Hz``,
``T/m/s`` <-> ``Hz/m/s``) need a gyromagnetic ratio, and ``ppm`` <-> ``Hz`` needs a Larmor
frequency.  Both are keyword arguments rather than globals, so a nucleus other than ``1H`` cannot
silently pick up the proton value -- take them from the ``Opts`` the sequence is built against:

>>> import pypulseq as pp
>>> import seqcraft as sc
>>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', B0=3.0)
>>> round(sc.convert(-3.4, 'ppm', 'Hz', f0=opts.gamma * opts.B0))   # fat/water shift at 3 T
-434
>>> round(sc.convert(opts.max_grad, 'Hz/m', 'mT/m', gamma=opts.gamma), 1)
40.0
"""

from __future__ import annotations

import difflib
import math
from fractions import Fraction
from typing import NamedTuple

from .errors import ConfigurationError, format_error

__all__ = ['GAMMA_1H', 'convert', 'dimension_of', 'dimensions', 'known_units']

#: Proton gyromagnetic ratio in Hz/T, matching ``pypulseq.Opts`` default ``gamma``.
GAMMA_1H = 42_576_000.0

_TWO_PI = 2.0 * math.pi


class _Unit(NamedTuple):
    """
    One row of the unit table.

    ``value * scale * gamma**gamma_power`` is the value in the dimension's canonical unit.
    `scale` is a :class:`~fractions.Fraction` unless pi enters it.
    """

    dimension: str
    scale: Fraction | float
    gamma_power: int = 0


def _f(numerator: int, denominator: int = 1) -> Fraction:
    return Fraction(numerator, denominator)


#: The unit table.  A unit is legal exactly when it is a key here.
#:
#: Field and frequency share one dimension on purpose: in pulseq they *are* one dimension -- B1 is
#: carried in hertz -- so a B1 limit in uT, a Larmor frequency in MHz and an off-resonance in Hz are
#: all the same conversion.
_UNITS: dict[str, _Unit] = {
    # ------------------------------------------------------------------------------- time
    's': _Unit('time', _f(1)),
    'ms': _Unit('time', _f(1, 10**3)),
    'us': _Unit('time', _f(1, 10**6)),
    'ns': _Unit('time', _f(1, 10**9)),
    'ps': _Unit('time', _f(1, 10**12)),
    # Ten microseconds: the block-duration quantum a .seq file counts in.  This is a plain unit,
    # *not* a synonym for "the block raster" -- for counts of this scanner's raster use
    # Raster(opts.block_duration_raster).count(t), which follows the scanner if it ever changes.
    '10us': _Unit('time', _f(1, 10**5)),
    'min': _Unit('time', _f(60)),
    # ----------------------------------------------------------------------------- length
    'm': _Unit('length', _f(1)),
    'cm': _Unit('length', _f(1, 10**2)),
    'mm': _Unit('length', _f(1, 10**3)),
    'um': _Unit('length', _f(1, 10**6)),
    'nm': _Unit('length', _f(1, 10**9)),
    # ------------------------------------------------------------------------------ angle
    'rad': _Unit('angle', _f(1)),
    'deg': _Unit('angle', math.pi / 180.0),
    'turn': _Unit('angle', _TWO_PI),
    # --------------------------------------------- frequency and field, related by gamma
    'Hz': _Unit('frequency', _f(1)),
    'kHz': _Unit('frequency', _f(10**3)),
    'MHz': _Unit('frequency', _f(10**6)),
    'rad/s': _Unit('frequency', 1.0 / _TWO_PI),
    'T': _Unit('frequency', _f(1), 1),
    'mT': _Unit('frequency', _f(1, 10**3), 1),
    'uT': _Unit('frequency', _f(1, 10**6), 1),
    'G': _Unit('frequency', _f(1, 10**4), 1),
    # ------------------------------------------- gradient amplitude, tesla family via gamma
    'Hz/m': _Unit('gradient', _f(1)),
    'rad/ms/mm': _Unit('gradient', 1e6 / _TWO_PI),
    'T/m': _Unit('gradient', _f(1), 1),
    'mT/m': _Unit('gradient', _f(1, 10**3), 1),
    'G/cm': _Unit('gradient', _f(1, 10**2), 1),
    # --------------------------------------------------- slew rate, tesla family via gamma
    'Hz/m/s': _Unit('slew', _f(1)),
    'rad/ms/mm/ms': _Unit('slew', 1e9 / _TWO_PI),
    'T/m/s': _Unit('slew', _f(1), 1),
    'mT/m/ms': _Unit('slew', _f(1), 1),
    'mT/m/s': _Unit('slew', _f(1, 10**3), 1),
    'G/cm/ms': _Unit('slew', _f(10), 1),
    # ---------------------------------- spatial frequency: k-space positions, gradient areas
    '1/m': _Unit('kspace', _f(1)),
    '1/mm': _Unit('kspace', _f(10**3)),
    '1/cm': _Unit('kspace', _f(10**2)),
    'rad/m': _Unit('kspace', 1.0 / _TWO_PI),
    # Higher gradient moments, named as validate.DEFAULT_RANGES names them, so a report and an
    # error message can quote a unit the reader can pass straight back into convert.
    '1/m/s': _Unit('kspace_rate', _f(1)),
    '1/mm/s': _Unit('kspace_rate', _f(10**3)),
    '1/cm/s': _Unit('kspace_rate', _f(10**2)),
    '1/m^2': _Unit('kspace_area', _f(1)),
    '1/mm^2': _Unit('kspace_area', _f(10**6)),
    '1/cm^2': _Unit('kspace_area', _f(10**4)),
    # --------------------------------------------------------------------------- b-value
    's/m^2': _Unit('bvalue', _f(1)),
    's/mm^2': _Unit('bvalue', _f(10**6)),
    'ms/um^2': _Unit('bvalue', _f(10**9)),
    # ----------------------------------------------------------------------------- ratio
    '1': _Unit('ratio', _f(1)),
    '%': _Unit('ratio', _f(1, 10**2)),
    'ppm': _Unit('ratio', _f(1, 10**6)),
    'ppb': _Unit('ratio', _f(1, 10**9)),
}

#: The unit a dimension converts to when `to_unit` is omitted.
_CANONICAL: dict[str, str] = {
    'time': 's',
    'length': 'm',
    'angle': 'rad',
    'frequency': 'Hz',
    'gradient': 'Hz/m',
    'slew': 'Hz/m/s',
    'kspace': '1/m',
    'kspace_rate': '1/m/s',
    'kspace_area': '1/m^2',
    'bvalue': 's/m^2',
    'ratio': '1',
}

#: Ratios are the hot path -- one per module parameter per build -- so the exact Fraction division
#: happens once per unit pair and the result is reused.
_RATIO: dict[tuple[str, str], tuple[float, bool]] = {}


def _resolve(unit: str) -> _Unit:
    """Look a unit up, naming the closest known spelling if it is not one."""
    try:
        return _UNITS[unit]
    except (KeyError, TypeError):
        fields: dict[str, object] = {}
        close = difflib.get_close_matches(str(unit), list(_UNITS), n=3)
        if close:
            fields['did you mean'] = ', '.join(repr(c) for c in close)
        fields['known units'] = ', '.join(sorted(_UNITS))
        raise ConfigurationError(
            format_error(f'unknown unit {unit!r}.', fields)
        ) from None


def _ratio(from_unit: str, to_unit: str) -> tuple[float, bool]:
    """
    Return ``(factor, divide)`` for ``scale_from / scale_to``.

    The exact Fraction is computed first, then collapsed to **one** float operation.  When the
    ratio is a reciprocal integer -- ``us`` to ``ms``, ``mm`` to ``m``, the common cases -- that
    operation is a division by an exactly-representable integer rather than a multiplication by an
    inexact fraction, so ``4200 us`` is ``4.2 ms`` and not ``4.200000000000001 ms``.
    """
    key = (from_unit, to_unit)
    cached = _RATIO.get(key)
    if cached is None:
        exact = _UNITS[from_unit].scale / _UNITS[to_unit].scale
        if isinstance(exact, Fraction) and exact.numerator == 1 and exact.denominator != 1:
            cached = (float(exact.denominator), True)
        else:
            cached = (float(exact), False)
        _RATIO[key] = cached
    return cached


def convert(
    value: float,
    from_unit: str,
    to_unit: str | None = None,
    *,
    gamma: float = GAMMA_1H,
    f0: float | None = None,
) -> float:
    """
    Convert `value` from one unit to another.

    Parameters
    ----------
    value
        Scalar or numpy array.  Arrays convert elementwise, so a waveform is one call.
    from_unit
        The unit `value` is in.  See :func:`known_units`.
    to_unit
        The unit wanted.  Defaults to the dimension's canonical unit (``'s'``, ``'m'``, ``'rad'``,
        ``'Hz'``, ``'Hz/m'``, ``'Hz/m/s'``, ``'1/m'``, ``'s/m^2'``, ``'1'``).  Passing it
        explicitly reads better and is what module code should do.
    gamma
        Gyromagnetic ratio, Hz/T.  Used only by pairs that cross the tesla/hertz divide.  Defaults
        to the proton value; pass ``gamma=opts.gamma`` to use the scanner's own.
    f0
        Larmor frequency, Hz.  Required only for ``ppm``/``%`` <-> frequency.

    Returns
    -------
    float
        `value` in `to_unit`.  Same shape as `value`.

    Raises
    ------
    ConfigurationError
        If a unit is unknown, if the two units belong to unrelated dimensions, or if a
        conversion needing `f0` did not get one.

    Examples
    --------
    Time, exactly, in both directions:

    >>> convert(4200, 'us', 'ms')
    4.2
    >>> convert(0.0121, 's', '10us')
    1210.0

    Gradients and slew rates, with and without gamma:

    >>> round(convert(1703040.0, 'Hz/m', 'mT/m'), 6)
    40.0
    >>> convert(200, 'T/m/s', 'mT/m/ms')                  # the same rate, spelled twice
    200.0
    >>> round(convert(180, 'mT/m/ms', 'G/cm/ms'), 1)
    18.0

    Field, B1 and frequency are one dimension:

    >>> round(convert(12, 'uT', 'Hz'), 1)                  # a B1 limit
    510.9
    >>> round(convert(500.0, 'Hz', 'uT'), 3)
    11.744

    Chemical shift needs a Larmor frequency:

    >>> round(convert(-3.4, 'ppm', 'Hz', f0=convert(3.0, 'T', 'Hz')))
    -434

    b-values, k-space and angles:

    >>> convert(1000, 's/mm^2', 's/m^2')
    1000000000.0
    >>> convert(2.0, '1/mm', '1/m')
    2000.0
    >>> round(convert(90, 'deg', 'rad'), 6)
    1.570796

    Arrays convert elementwise:

    >>> import numpy as np
    >>> convert(np.array([10.0, 20.0]), 'mT/m', 'G/cm')
    array([1., 2.])
    """
    src = _resolve(from_unit)
    if to_unit is None:
        to_unit = _CANONICAL[src.dimension]
    dst = _resolve(to_unit)

    if src.dimension == dst.dimension:
        factor, divide = _ratio(from_unit, to_unit)
        out = value / factor if divide else value * factor
        power = src.gamma_power - dst.gamma_power
        return out * gamma**power if power else out

    # The one bridge between dimensions: a dimensionless shift becomes a frequency at f0, which is
    # the fat/water offset, a CEST offset, and every "ppm" a spectroscopist quotes.
    if {src.dimension, dst.dimension} == {'ratio', 'frequency'}:
        if f0 is None:
            raise ConfigurationError(
                format_error(
                    f'converting {from_unit!r} to {to_unit!r} needs a Larmor frequency.',
                    {'why': 'a chemical shift is a fraction of f0, so f0 sets the hertz per ppm'},
                    [
                        'pass f0=<Hz>, e.g. f0=convert(3.0, "T", "Hz")',
                        'or use system.convert(...), which supplies gamma * B0',
                    ],
                )
            )
        if src.dimension == 'ratio':
            as_hz = value * float(src.scale) * f0
            return convert(as_hz, 'Hz', to_unit, gamma=gamma)
        return convert(value, from_unit, 'Hz', gamma=gamma) / f0 / float(dst.scale)

    raise ConfigurationError(
        format_error(
            f'cannot convert {from_unit!r} to {to_unit!r}.',
            {
                from_unit: f'a {src.dimension}',
                to_unit: f'a {dst.dimension}',
                f'{src.dimension} units': ', '.join(known_units(src.dimension)),
            },
        )
    )


def dimension_of(unit: str) -> str:
    """
    Return the dimension `unit` measures.

    Examples
    --------
    >>> dimension_of('mT/m')
    'gradient'
    >>> dimension_of('uT')                # B1 is carried in hertz, so field is frequency
    'frequency'
    """
    return _resolve(unit).dimension


def known_units(dimension: str | None = None) -> tuple[str, ...]:
    """
    Return the units :func:`convert` understands, sorted, optionally one dimension's.

    Examples
    --------
    >>> known_units('bvalue')
    ('ms/um^2', 's/m^2', 's/mm^2')
    >>> 'G/cm' in known_units()
    True
    """
    if dimension is None:
        return tuple(sorted(_UNITS))
    if dimension not in _CANONICAL:
        raise ConfigurationError(
            format_error(
                f'unknown dimension {dimension!r}.',
                {'known': ', '.join(sorted(_CANONICAL))},
            )
        )
    return tuple(sorted(u for u, spec in _UNITS.items() if spec.dimension == dimension))


def dimensions() -> tuple[str, ...]:
    """
    Return every dimension name.

    Examples
    --------
    >>> dimensions()[:3]
    ('angle', 'bvalue', 'frequency')
    """
    return tuple(sorted(_CANONICAL))
