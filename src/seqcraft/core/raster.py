"""
Exact raster arithmetic.

Pulseq quantises time on four rasters (typical Siemens values in brackets):

=====================  ==========  ==========================================
raster                 value       what must land on it
=====================  ==========  ==========================================
``block_duration``     10 us       every block duration
``grad_raster_time``   10 us       gradient waveform samples
``rf_raster_time``     1 us        RF waveform samples
``adc_raster_time``    100 ns      ADC dwell time
=====================  ==========  ==========================================

Two float traps have to be avoided, and both bite in practice.

**Rounding at the boundary.**  Naive ``math.ceil(t / r) * r`` is wrong for an exact
multiple: ``1.5e-3 / 1e-5`` evaluates to ``149.99999999999997``, so a value already on the
raster gets pushed up by a whole raster.

**Reconstructing the result.**  Even with the right integer count, ``n * raster`` is
inexact, because a decimal raster such as ``1e-7`` has no exact binary representation:
``250 * 1e-7`` is ``2.4999999999999998e-05``, not ``2.5e-05``.  That noise then propagates
into ADC dwell times, block durations and, eventually, into ``.seq`` text that differs
between builds.

The fix for both is to do the arithmetic in **integer picoseconds** and divide once at the
end.  A picosecond is 10^5 times finer than the finest pulseq raster (100 ns) and float64
represents integers exactly up to 2^53 ps, which is 2.5 hours -- far longer than any
sequence.  Dividing an exact integer by ``1e12`` lands on the nearest float to the true
decimal value, which is what makes results reproducible and doctestable.

Examples
--------
>>> ceil_to(1.5e-3, 1e-5) == 1.5e-3
True
>>> ceil_to(1.5001e-3, 1e-5)
0.00151
>>> floor_to(3.2e-3 / 128, 1e-7)          # the exact dwell, not 2.4999999999999998e-05
2.5e-05
>>> on_raster(1.5e-3, 1e-5)
True
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

from .errors import RasterError, format_error

if TYPE_CHECKING:
    from collections.abc import Iterable

__all__ = [
    'EPS',
    'PS',
    'ceil_to',
    'floor_to',
    'on_raster',
    'picoseconds',
    'quantise',
    'require_on_raster',
    'round_to',
    'seconds',
    'sub_exact',
    'sum_exact',
]

#: Absolute time tolerance, in seconds.  One nanosecond: 100x finer than the finest
#: pulseq raster and ~1e6x coarser than float64 resolution at millisecond magnitudes.
EPS = 1e-9

#: Picoseconds per second.  The internal integer time unit.
PS = 1_000_000_000_000


def picoseconds(t: float) -> int:
    """
    Convert seconds to integer picoseconds.

    Examples
    --------
    >>> picoseconds(2.5e-05)
    25000000
    """
    return int(round(t * PS))


def seconds(t_ps: int) -> float:
    """
    Convert integer picoseconds back to seconds.

    Examples
    --------
    >>> seconds(25000000)
    2.5e-05
    """
    return t_ps / PS


def sum_exact(values: Iterable[float]) -> float:
    """
    Sum durations without accumulating float error.

    Adding raster-aligned floats does not stay raster-aligned.  A solved TE would then
    compare unequal to the raster value it was derived from, and a downstream ``ceil_to``
    would push it up by a whole raster.  Summing in integer picoseconds removes the problem.

    Examples
    --------
    >>> terms = [1.3e-4, 5.0e-4, 3.7e-4, 2.1e-4]      # all exact multiples of 10 us
    >>> sum(terms)                                     # plain addition drifts
    0.0012100000000000001
    >>> sum_exact(terms)
    0.00121
    >>> on_raster(sum(terms), 1e-5), on_raster(sum_exact(terms), 1e-5)
    (True, True)
    """
    return seconds(sum(picoseconds(v) for v in values))


def sub_exact(a: float, b: float) -> float:
    """
    Return ``a - b`` computed in integer picoseconds.

    Examples
    --------
    >>> sub_exact(0.008, 0.00324)
    0.00476
    """
    return seconds(picoseconds(a) - picoseconds(b))


def _raster_ps(raster: float) -> int:
    """Validate a raster and return it in integer picoseconds."""
    if raster <= 0:
        msg = format_error(
            f'raster must be positive, got {raster!r}.',
            {'hint': 'pass e.g. system.block_raster_s (10 us), not 0'},
        )
        raise RasterError(msg)
    r_ps = int(round(raster * PS))
    if r_ps <= 0:
        msg = format_error(
            f'raster {raster!r} s is finer than one picosecond.',
            {'hint': 'seqcraft quantises time in picoseconds; pulseq rasters are >= 100 ns'},
        )
        raise RasterError(msg)
    return r_ps


def ceil_to(t: float, raster: float) -> float:
    """
    Round `t` up to a multiple of `raster`, leaving exact multiples untouched.

    Parameters
    ----------
    t
        Time in seconds. May be negative.
    raster
        Raster in seconds; must be positive.

    Returns
    -------
    float
        The smallest multiple of `raster` that is ``>= t``.

    Examples
    --------
    >>> ceil_to(1.5e-3, 1e-5) == 1.5e-3        # exact multiple: unchanged
    True
    >>> ceil_to(1.23e-5, 1e-5)
    2e-05
    >>> ceil_to(-1.23e-5, 1e-5)
    -1e-05
    """
    r = _raster_ps(raster)
    return seconds(-(-picoseconds(t) // r) * r)


def floor_to(t: float, raster: float) -> float:
    """
    Round `t` down to a multiple of `raster`, leaving exact multiples untouched.

    Used for ADC dwell times, which Siemens truncates rather than rounds up so the
    readout never outgrows its gradient.

    Examples
    --------
    >>> floor_to(1.5e-3, 1e-5) == 1.5e-3
    True
    >>> floor_to(1.29e-5, 1e-5)
    1e-05
    """
    r = _raster_ps(raster)
    return seconds((picoseconds(t) // r) * r)


def round_to(t: float, raster: float) -> float:
    """
    Round `t` to the nearest multiple of `raster`.

    Examples
    --------
    >>> round_to(1.4e-5, 1e-5)
    1e-05
    >>> round_to(1.6e-5, 1e-5)
    2e-05
    """
    r = _raster_ps(raster)
    t_ps = picoseconds(t)
    return seconds(((t_ps + (r // 2 if t_ps >= 0 else -(r // 2))) // r) * r)


def on_raster(t: float, raster: float) -> bool:
    """
    Report whether `t` is an exact multiple of `raster` to within :data:`EPS`.

    Examples
    --------
    >>> on_raster(1.5e-3, 1e-5)
    True
    >>> on_raster(1.5e-3 + 1e-6, 1e-5)
    False

    Float noise below one picosecond is tolerated, since that is finer than any raster:

    >>> on_raster(1.5e-3 * (1 + 1e-15), 1e-5)
    True
    """
    return picoseconds(t) % _raster_ps(raster) == 0


def quantise(t: float, raster: float) -> int:
    """
    Return `t` as an exact integer count of `raster` units.

    Raises
    ------
    RasterError
        If `t` is not a multiple of `raster`.  Integer counts are what let durations be
        summed and compared without accumulating float error, so this is a hard error
        rather than a silent round.

    Examples
    --------
    >>> quantise(1.5e-3, 1e-5)
    150
    """
    require_on_raster(t, raster, what='time')
    return picoseconds(t) // _raster_ps(raster)


def require_on_raster(t: float, raster: float, *, what: str = 'duration') -> None:
    """
    Raise :class:`~seqcraft.core.errors.RasterError` unless `t` lands on `raster`.

    The message names the two nearest valid values so the fix is mechanical.

    Examples
    --------
    >>> require_on_raster(1.5e-3, 1e-5)
    >>> require_on_raster(1.505e-3, 1e-5)
    Traceback (most recent call last):
        ...
    seqcraft.core.errors.RasterError: duration 1.505000 ms is not a multiple of the 10.0 us raster.
      nearest below:  1.500000 ms
      nearest above:  1.510000 ms
      fix
        round with seqcraft.core.raster.ceil_to(t, raster)
    """
    if on_raster(t, raster):
        return
    msg = format_error(
        f'{what} {t * 1e3:.6f} ms is not a multiple of the {raster * 1e6:.1f} us raster.',
        {
            'nearest below': f'{floor_to(t, raster) * 1e3:.6f} ms',
            'nearest above': f'{ceil_to(t, raster) * 1e3:.6f} ms',
        },
        ['round with seqcraft.core.raster.ceil_to(t, raster)'],
    )
    raise RasterError(msg)
