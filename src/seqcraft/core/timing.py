"""
Exact time arithmetic: rasters as objects, addition without drift.

Pulseq quantises time on four rasters.  Every one of them is a **value the scanner supplies**, not
a constant in this file -- Siemens' 10 us gradient raster, GE's 4 us and Philips' 6.4 us are all
just numbers passed in, and a raster invented next year needs no change here.

=====================  ==================  ==========================================
raster                 ``System``          what must land on it
=====================  ==================  ==========================================
block duration         ``block_raster``    every block duration
gradient               ``grad_raster``     gradient waveform samples
RF                     ``rf_raster``       RF waveform samples
ADC                    ``adc_raster``      ADC dwell time
=====================  ==================  ==========================================

A raster is a :class:`Raster`, so the quantising operations are methods on the thing being
quantised onto and no call site carries a bare float:

>>> block = Raster(10e-6, 'block')
>>> block.ceil(1.5001e-3)
0.00151
>>> block.ceil(1.5e-3) == 1.5e-3                    # already on the raster: untouched
True
>>> block.count(1.5e-3)                             # exact integer count
150
>>> block.at(150)                                   # and back, exactly
0.0015

Two float traps, and why the arithmetic is integer
--------------------------------------------------
**Rounding at the boundary.**  Naive ``math.ceil(t / r) * r`` is wrong for an exact multiple:
``1.5e-3 / 1e-5`` evaluates to ``149.99999999999997``, so a value already on the raster gets pushed
up by a whole raster.

**Reconstructing the result.**  Even with the right integer count, ``n * raster`` is inexact,
because a decimal raster such as ``1e-7`` has no exact binary representation: ``250 * 1e-7`` is
``2.4999999999999998e-05``, not ``2.5e-05``.  That noise propagates into ADC dwell times, block
durations and eventually into ``.seq`` text that differs between builds.

Both are fixed by doing the arithmetic in **integer ticks** and dividing once at the end.

>>> adc = Raster(100e-9, 'adc')
>>> adc.floor(3.2e-3 / 128)                         # the exact dwell
2.5e-05
>>> Raster(6.4e-6, 'grad').ceil(1e-3)               # a non-Siemens raster works the same
0.0010048

Examples
--------
Durations are summed exactly, so a solved TE still compares equal to the raster value it came from:

>>> terms = [1e-5] * 121                             # all exact multiples of 10 us
>>> sum(terms) == 0.00121                           # plain addition drifts
False
>>> exact_sum(terms)
0.00121
>>> exact_diff(0.008, 0.00324)
0.00476
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .errors import RasterError, format_error

if TYPE_CHECKING:
    from collections.abc import Iterable

__all__ = [
    'EPS',
    'TICKS_PER_SECOND',
    'Raster',
    'exact_diff',
    'exact_sum',
    'from_ticks',
    'to_ticks',
]

#: Absolute time tolerance, in seconds.  One nanosecond: 100x finer than the finest pulseq raster
#: and ~1e6x coarser than float64 resolution at millisecond magnitudes.
EPS = 1e-9

#: Ticks per second.  One tick is a picosecond: 10^5 finer than the finest pulseq raster (100 ns),
#: and float64 holds integers exactly to 2^53 ticks, which is 2.5 hours -- far longer than any
#: sequence.  Dividing an exact integer count of ticks by this lands on the nearest float to the
#: true decimal value, which is what makes results reproducible and doctestable.
TICKS_PER_SECOND = 10**12


def to_ticks(t: float) -> int:
    """
    Return `t` seconds as an integer number of ticks.

    This is **not** a unit conversion -- :func:`seqcraft.units.convert` does that, and
    ``convert(t, 's', 'ps')`` is the conversion this resembles.  Ticks are the exact-integer time
    domain that lets durations be added, compared and quantised without float drift, and they are
    an implementation detail of this module and the compiler.  Module code should not need them.

    Examples
    --------
    >>> to_ticks(2.5e-05)
    25000000
    """
    return int(round(t * TICKS_PER_SECOND))


def from_ticks(n: int) -> float:
    """
    Return `n` ticks as seconds, correctly rounded.

    Examples
    --------
    >>> from_ticks(25000000)
    2.5e-05
    """
    return n / TICKS_PER_SECOND


def exact_sum(values: Iterable[float]) -> float:
    """
    Sum durations without accumulating float error.

    Adding raster-aligned floats does not stay raster-aligned.  A solved TE would then compare
    unequal to the raster value it was derived from, and a downstream :meth:`Raster.ceil` would
    push it up by a whole raster.  Summing in ticks removes the problem.

    Examples
    --------
    >>> exact_sum([1.3e-4, 5.0e-4, 3.7e-4, 2.1e-4])
    0.00121
    >>> Raster(1e-5).holds(exact_sum([1.3e-4, 5.0e-4, 3.7e-4, 2.1e-4]))
    True
    """
    return from_ticks(sum(to_ticks(v) for v in values))


def exact_diff(a: float, b: float) -> float:
    """
    Return ``a - b`` computed in ticks.

    Examples
    --------
    >>> exact_diff(0.008, 0.00324)
    0.00476
    """
    return from_ticks(to_ticks(a) - to_ticks(b))


class Raster:
    """
    A time raster: exact arithmetic on integer multiples of one interval.

    Parameters
    ----------
    dt
        The raster interval in seconds.  Must be positive and at least one tick (1 ps).
    name
        Optional label -- ``'block'``, ``'gradient'`` -- used in error messages.

    Attributes
    ----------
    dt : float
        The interval in seconds.  This is the number to use in array arithmetic, e.g.
        ``np.diff(g) / system.grad_raster.dt`` for a slew rate.
    name : str

    Notes
    -----
    Hashable and comparable by value, so two rasters read from the same ``Opts`` are the same
    raster.  Get them from :class:`~seqcraft.System` rather than constructing them: the scanner is
    the authority on what its rasters are.

    Examples
    --------
    >>> import seqcraft as sc
    >>> system = sc.System.preset('prisma')
    >>> system.grad_raster
    Raster(gradient, 10 us)
    >>> system.adc_raster.dt
    1e-07
    >>> system.block_raster == Raster(10e-6)
    True
    """

    _ticks: int
    dt: float
    name: str

    __slots__ = ('_ticks', 'dt', 'name')

    def __init__(self, dt: float, name: str = '') -> None:
        if not dt > 0:
            msg = format_error(
                f'raster must be positive, got {dt!r}.',
                {'hint': 'pass e.g. system.block_raster.dt (10 us on Siemens), not 0'},
            )
            raise RasterError(msg)
        ticks = int(round(dt * TICKS_PER_SECOND))
        if ticks <= 0:
            msg = format_error(
                f'raster {dt!r} s is finer than one tick.',
                {
                    'one tick': f'{1 / TICKS_PER_SECOND:g} s',
                    'why': 'seqcraft quantises time in ticks; pulseq rasters are >= 100 ns',
                },
            )
            raise RasterError(msg)
        object.__setattr__(self, '_ticks', ticks)
        object.__setattr__(self, 'dt', from_ticks(ticks))
        object.__setattr__(self, 'name', name)

    # ------------------------------------------------------------------------ quantising
    def holds(self, t: float) -> bool:
        """
        Report whether `t` is an exact multiple of this raster, to within :data:`EPS`.

        Float noise below one tick is tolerated, since that is finer than any raster.

        Examples
        --------
        >>> r = Raster(1e-5)
        >>> r.holds(1.5e-3), r.holds(1.5e-3 + 1e-6)
        (True, False)
        >>> r.holds(1.5e-3 * (1 + 1e-15))
        True
        """
        return to_ticks(t) % self._ticks == 0

    def ceil(self, t: float) -> float:
        """
        Return the smallest multiple of this raster that is ``>= t``, exactly.

        Examples
        --------
        >>> r = Raster(1e-5)
        >>> r.ceil(1.5e-3) == 1.5e-3        # exact multiple: unchanged
        True
        >>> r.ceil(1.23e-5)
        2e-05
        >>> r.ceil(-1.23e-5)
        -1e-05
        """
        n = self._ticks
        return from_ticks(-(-to_ticks(t) // n) * n)

    def floor(self, t: float) -> float:
        """
        Return the largest multiple of this raster that is ``<= t``, exactly.

        Used for ADC dwell times, which Siemens truncates rather than rounds up so the readout
        never outgrows its gradient.

        Examples
        --------
        >>> r = Raster(1e-5)
        >>> r.floor(1.5e-3) == 1.5e-3
        True
        >>> r.floor(1.29e-5)
        1e-05
        """
        n = self._ticks
        return from_ticks((to_ticks(t) // n) * n)

    def nearest(self, t: float) -> float:
        """
        Return the multiple of this raster closest to `t`, exactly.

        Ties go to the larger value, on both sides of zero -- ``+half`` then floor, with no branch
        on the sign.  Branching on it (and subtracting half for negatives) rounds ``-1.4`` rasters
        to ``-2``, which is not the nearest multiple.

        Examples
        --------
        >>> r = Raster(1e-5)
        >>> r.nearest(1.4e-5), r.nearest(1.6e-5)
        (1e-05, 2e-05)
        >>> r.nearest(-1.4e-5), r.nearest(-1.6e-5)
        (-1e-05, -2e-05)
        """
        n = self._ticks
        return from_ticks(((to_ticks(t) + n // 2) // n) * n)

    def count(self, t: float) -> int:
        """
        Return `t` as an exact integer count of this raster.

        Raises
        ------
        RasterError
            If `t` is not a multiple.  Integer counts are what let durations be summed and
            compared without accumulating float error, so this is a hard error rather than a
            silent round -- use :meth:`ceil` or :meth:`nearest` first if that is what you meant.

        Examples
        --------
        >>> Raster(1e-5).count(1.5e-3)
        150
        """
        self.require(t, what='time')
        return to_ticks(t) // self._ticks

    def at(self, n: int) -> float:
        """
        Return ``n`` rasters as seconds, exactly.

        Examples
        --------
        >>> Raster(1e-7).at(250)                  # not 2.4999999999999998e-05
        2.5e-05
        >>> Raster(1e-5).at(150)
        0.0015
        """
        return from_ticks(int(n) * self._ticks)

    def require(self, t: float, *, what: str = 'duration') -> None:
        """
        Raise :class:`~seqcraft.core.errors.RasterError` unless `t` lands on this raster.

        The message names the two nearest valid values, so the fix is mechanical.

        Examples
        --------
        >>> block = Raster(1e-5, 'block')
        >>> block.require(1.5e-3)
        >>> block.require(1.505e-3)
        Traceback (most recent call last):
            ...
        seqcraft.core.errors.RasterError: duration 1.505000 ms is not a multiple of the 10.0 us block raster.
          nearest below:  1.500000 ms
          nearest above:  1.510000 ms
          fix
            round it first with <raster>.ceil(t)
        """
        if self.holds(t):
            return
        label = f'{self.dt * 1e6:.1f} us{" " + self.name if self.name else ""} raster'
        msg = format_error(
            f'{what} {t * 1e3:.6f} ms is not a multiple of the {label}.',
            {
                'nearest below': f'{self.floor(t) * 1e3:.6f} ms',
                'nearest above': f'{self.ceil(t) * 1e3:.6f} ms',
            },
            ['round it first with <raster>.ceil(t)'],
        )
        raise RasterError(msg)

    # -------------------------------------------------------------------------- protocol
    def __setattr__(self, key: str, value: object) -> None:
        """Reject mutation: a raster is a property of the scanner, not a variable."""
        msg = f'{type(self).__name__} is immutable'
        raise AttributeError(msg)

    def __float__(self) -> float:
        """The interval in seconds, for ``float(raster)``.  Prefer ``.dt`` in arithmetic."""
        return self.dt

    def __eq__(self, other: object) -> bool:
        """Equal when the interval is, to the tick.  The name is a label and does not count."""
        return isinstance(other, Raster) and other._ticks == self._ticks

    def __hash__(self) -> int:
        """Hash on the interval, so a raster can key a cache."""
        return hash((Raster, self._ticks))

    def __repr__(self) -> str:
        """``Raster(gradient, 10 us)``."""
        us = self.dt * 1e6
        pretty = f'{us:g} us' if us >= 1 else f'{self.dt * 1e9:g} ns'
        return f'Raster({self.name}, {pretty})' if self.name else f'Raster({pretty})'
