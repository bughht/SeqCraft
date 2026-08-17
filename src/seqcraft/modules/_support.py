"""
Argument checks and exact gradient arithmetic shared by the modules.

Nothing here is a :class:`~seqcraft.design.module.Module`, and nothing here is MR knowledge.  It
is the handful of lines that would otherwise be copied into five files: the same range check
written five ways drifts, and the version that drifts is the one whose message stops naming the
argument that fixes it.

:func:`area_until` is the one that is arithmetic rather than validation, and it is here rather
than in :mod:`seqcraft.design.events` because only the module library has ever needed it -- the
compiler integrates whole events, and a *partial* integral is a question about physics (where is
the echo?) rather than about legality.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..design.events import AXES, knots_of, pwl_moment
from ..design.timing import Raster
from ..errors import ConfigurationError, format_error

if TYPE_CHECKING:
    from collections.abc import Iterable

    from ..design.events import Event

__all__ = ['area_until', 'ceil_raster', 'require_axis', 'require_positive', 'require_range']


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


def ceil_raster(value: float, raster: float) -> float:
    """
    Return `value` rounded up onto `raster`.

    Thin, and named, because the alternative spelling -- ``ceil(v / dt) * dt`` -- accumulates
    float error at the millisecond magnitudes a TR is written in, and :class:`Raster` already
    does the arithmetic in integer ticks.

    Examples
    --------
    >>> ceil_raster(4.992e-3, 1e-5)
    0.005
    >>> ceil_raster(5e-3, 1e-5)                 # already on the raster: unchanged
    0.005
    """
    return float(Raster(raster).ceil(value))


def area_until(event: Event, t_end: float) -> float:
    """
    Return the exact area of a gradient event from its node up to `t_end`, in 1/m.

    Exact because a pulseq gradient *is* piecewise linear: :func:`~seqcraft.events.knots_of`
    gives its knots and :func:`~seqcraft.events.pwl_moment` integrates them in closed form, so
    the answer carries no discretisation error however coarse the raster.  The event's own
    ``delay`` is included, and `t_end` is measured from the event's node -- the same origin
    :func:`knots_of` uses.

    This is what a readout prephaser must cancel, and computing it rather than a closed form is
    what makes ramp sampling, partial echo and the half-dwell sample offset work without a
    derivation each.  ``-gx.area / 2`` is right for a symmetric full echo over the flat top up to
    half a ``dk``; ``-gx.flat_area / 2``, one letter away in spelling, is right for nothing at
    all and drops the whole ramp-up.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s')
    >>> g = pp.make_trapezoid(channel='x', flat_area=1000.0, flat_time=2e-3, system=o)
    >>> round(area_until(g, pp.calc_duration(g)), 6) == round(float(g.area), 6)
    True
    >>> round(area_until(g, 0.0), 6)
    0.0
    >>> half = area_until(g, g.rise_time + g.flat_time / 2)     # ramp-up plus half the flat top
    >>> abs(half - float(g.area) / 2) < 1e-9                    # exactly half, ramps included
    True
    """
    times, amps = knots_of(event)
    if times.size < 2 or t_end <= times[0]:
        return 0.0
    if t_end >= times[-1]:
        return float(pwl_moment(times, amps))
    cut = int(np.searchsorted(times, t_end))
    edge = float(np.interp(t_end, times, amps))
    return float(pwl_moment(np.append(times[:cut], t_end), np.append(amps[:cut], edge)))
