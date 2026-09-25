"""
Value-domain checks: refuse a bad parameter before anything is built from it.

Separate from :mod:`seqcraft.design.units`, which converts between units, and separate from
``modules/_support.py``, which needs pypulseq and knows about waveforms.  These five need neither,
and they sit here so that anything in the package can use them **without importing
:mod:`seqcraft.modules`**.

That is not tidiness.  The public augmentation intents live above ``modules`` and are imported
*by* the kernels, so an intent reaching into ``modules._support`` for `require_axis` would close
the loop ``augmentation -> modules -> kernel -> augmentation``.  Import cycles of that shape work
until someone changes an import order.

**Internal.**  ``modules/_support.py`` re-exports every name here, so existing call sites are
unchanged and the public surface is exactly what it was.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..errors import ConfigurationError, format_error
from .events import AXES

if TYPE_CHECKING:
    from collections.abc import Iterable

__all__ = [
    'require_axis',
    'require_count',
    'require_pair',
    'require_positive',
    'require_range',
]


def require_positive(value: float, name: str, *, fixes: Iterable[str] = ()) -> float:
    """
    Return `value` as a float, having checked that it is finite and above zero.

    Examples
    --------
    >>> require_positive(5.0, 'thickness_mm')
    5.0
    >>> require_positive(0.0, 'thickness_mm')
    Traceback (most recent call last):
        ...
    seqcraft.errors.ConfigurationError: thickness_mm must be positive, got 0.
    """
    number = float(value)
    if not np.isfinite(number) or number <= 0.0:
        msg = format_error(f'{name} must be positive, got {number:g}.', {name: number}, fixes)
        raise ConfigurationError(msg)
    return number


def require_range(value: float, name: str, *, low: float, high: float) -> float:
    """
    Return `value` as a float, having checked ``low < value <= high``.

    Half-open on purpose: every caller here bounds a *fraction* whose lower end is degenerate
    (a partial Fourier factor of zero samples nothing) and whose upper end is the ordinary case.

    Examples
    --------
    >>> require_range(0.75, 'partial_fourier', low=0.0, high=1.0)
    0.75
    >>> require_range(1.4, 'partial_fourier', low=0.0, high=1.0)
    Traceback (most recent call last):
        ...
    seqcraft.errors.ConfigurationError: partial_fourier must be in (0, 1], got 1.4.
    """
    number = float(value)
    if not np.isfinite(number) or not low < number <= high:
        msg = format_error(
            f'{name} must be in ({low:g}, {high:g}], got {number:g}.',
            {name: number},
        )
        raise ConfigurationError(msg)
    return number


def require_pair(value: float | tuple[float, float], name: str) -> tuple[float, float]:
    """
    Return `value` as an ``(x, y)`` pair, accepting a scalar as "the same on both axes".

    Here rather than in either caller because it has two, and they are in **different
    folders** -- :class:`~seqcraft.modules.GRE2DTR` in ``kernel/`` and
    :class:`~seqcraft.modules.EPI2D` in ``readout/``.  That is the same argument
    :func:`shift_slice` makes below it: both are modules that know the scan is
    two-dimensional, so both have to split a pair the single-axis leaves beneath them never
    see, and a cross-folder import of another folder's private name is worse than one shared
    file whose whole purpose is being shared.

    Examples
    --------
    >>> require_pair(250.0, 'fov_mm')
    (250.0, 250.0)
    >>> require_pair((250.0, 180.0), 'fov_mm')
    (250.0, 180.0)
    >>> require_pair((256, 128, 64), 'matrix')
    Traceback (most recent call last):
        ...
    seqcraft.errors.ConfigurationError: matrix must be a pair (x, y), got 3 values.
    """
    if isinstance(value, (int, float)):
        return (float(value), float(value))
    try:
        first, second = value
    except (TypeError, ValueError):
        count = len(value) if hasattr(value, '__len__') else '?'
        msg = format_error(
            f'{name} must be a pair (x, y), got {count} values.',
            {name: value},
            [f'{name}=250.0 means square', f'{name}=(250.0, 180.0) is readout then phase'],
        )
        raise ConfigurationError(msg) from None
    return (float(first), float(second))


def require_count(value: int, name: str, *, low: int = 1, hint: str = '') -> int:
    """
    Return `value` as an int at or above `low`.

    Five arguments across two modules are counts with a floor -- `oversampling`, `blip_lines` and
    `navigator_echoes` on :class:`~seqcraft.modules.EPI2D`, `echoes` on
    :class:`~seqcraft.modules.CartesianLine` -- and one function is better than five that differ
    only in the noun.  It was private to ``epi_2d.py`` while those were three, on the rule that
    keeps this file the shared things and not the leftovers; the fourth caller, in a module
    ``epi_2d`` must not be imported *by*, is what moved it.

    Examples
    --------
    >>> require_count(4, 'blip_lines')
    4
    >>> require_count(-1, 'navigator_echoes', low=0)
    Traceback (most recent call last):
        ...
    seqcraft.errors.ConfigurationError: navigator_echoes must not be negative, got -1.
    """
    count = int(value)
    if count < low:
        floor = 'not be negative' if low == 0 else f'be at least {low}'
        msg = format_error(f'{name} must {floor}, got {count}.', {name: count},
                           [hint] if hint else ())
        raise ConfigurationError(msg)
    return count


def require_axis(axis: str, name: str = 'axis') -> str:
    """
    Return `axis` having checked it names a logical gradient channel.

    Examples
    --------
    >>> require_axis('y')
    'y'
    >>> require_axis('ky')
    Traceback (most recent call last):
        ...
    seqcraft.errors.ConfigurationError: axis must be one of 'x', 'y', 'z', got 'ky'.
    """
    if axis not in AXES:
        listed = ', '.join(repr(a) for a in AXES)
        msg = format_error(
            f'{name} must be one of {listed}, got {axis!r}.',
            {name: axis},
            ['the logical axes are the pulseq channel names, not k-space names'],
        )
        raise ConfigurationError(msg)
    return axis

