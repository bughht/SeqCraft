"""
Argument checks and exact gradient arithmetic shared by the modules.

Nothing here is a :class:`~seqcraft.design.module.Module`, and nothing here is MR knowledge.  It
is the handful of lines that would otherwise be copied into five files: the same range check
written five ways drifts, and the version that drifts is the one whose message stops naming the
argument that fixes it.

:func:`area_until` is the one that is arithmetic rather than validation, and it is here rather
than in :mod:`seqcraft.design.events` because only the module library has ever needed it -- the
compiler integrates whole events, and a *partial* integral is a question about physics (where is
the echo?) rather than about legality.  :func:`shift_slice` is the second of those, and it is here
for a sharper reason: it has two consumers, :class:`~seqcraft.modules.Excitation` in ``rf/`` and
:class:`~seqcraft.modules.IRPrep` in ``preparation/``, and they are in **different folders**.  A
cross-folder import of another folder's private name is worse than one shared file whose whole
purpose is being shared.
"""

from __future__ import annotations

from math import pi
from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ..design.events import AXES, derive, knots_of, pwl_moment
from ..design.timing import Raster
from ..errors import ConfigurationError, format_error

if TYPE_CHECKING:
    from collections.abc import Iterable

    from ..design.events import Event

__all__ = [
    'area_until', 'ceil_raster', 'require_axis', 'require_pair', 'require_positive',
    'require_range', 'shift_slice',
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


def shift_slice(rf: Event, gz: Event, *, position_m: float) -> Event:
    """
    Return `rf` retuned to excite the slice `position_m` from isocentre along `gz`.

    Two lines, and the second is the one that gets omitted::

        freq_offset  = gz.amplitude * position_m
        phase_offset = rf.phase_offset - 2*pi * freq_offset * calc_rf_center(rf)[0]

    Without the phase reference every off-centre slice carries a phase that depends on where the
    pulse's effective centre falls, so multi-slice is silently wrong -- and for a minimum-phase
    pulse, whose centre is not its midpoint, it is wrong by a different amount than for a sinc.
    A 10 ms hyperbolic secant is the extreme case: its centre is 5 ms in, so the term it drops is
    five milliseconds' worth of a frequency that can be tens of kilohertz.

    This is the one real gap in pypulseq's coverage: ``make_sinc_pulse`` takes ``slice_thickness``
    but no position, because it builds one event at a time from explicit numbers, and this needs
    two events plus a geometric input.  It lives in the module library rather than in ``design/``
    because it is **MR imaging knowledge** -- it knows what a slice is -- and the core holds none.
    Returning an event rather than a block makes it *not a module*; it does not make it core.

    Not public API: the spelling a caller writes is ``exc(position_mm=20.0)`` or
    ``irprep(position_mm=20.0)``.  It is unprefixed inside this already-private file so that both
    of them import a name rather than a name with an underscore in it.

    Examples
    --------
    >>> import numpy as np, pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s')
    >>> rf, gz, _ = pp.make_sinc_pulse(flip_angle=np.pi / 2, duration=3e-3,
    ...                                slice_thickness=5e-3, return_gz=True, system=o)
    >>> shifted = shift_slice(rf, gz, position_m=0.02)
    >>> round(shifted.freq_offset)                      # gz.amplitude * 20 mm
    5333
    >>> abs(shift_slice(rf, gz, position_m=0.0).freq_offset)
    0.0
    """
    freq_offset = float(gz.amplitude) * position_m
    phase_offset = float(rf.phase_offset) - 2 * pi * freq_offset * float(
        pp.calc_rf_center(rf)[0]
    )
    return derive(rf, freq_offset=freq_offset, phase_offset=phase_offset)
