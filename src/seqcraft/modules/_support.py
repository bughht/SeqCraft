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

:func:`halve_onto` and :func:`dwell_quantum` are here on a *related* argument rather than that
one, and the difference is worth stating because it is the first time this file has taken
something for a reason other than the folders.  Both came out of
:class:`~seqcraft.modules.EPI2D`, and both are now also called by
:class:`~seqcraft.modules.CartesianLine` -- which is in the **same** folder, ``readout/``.  What
makes them shared rather than imported across is the *direction*: ``CartesianLine`` is the older,
smaller module that ``GRE2DTR``, ``MPRAGE2D`` and three example notebooks already stand on, and
``EPI2D`` is the newer one that stands on nothing.  ``from .epi_2d import _halve`` would point the
dependency the wrong way round -- every existing GRE would import the EPI train to find out how to
split a guard -- and it would do it through a private name.  One shared file, and neither module
knows the other exists.
"""

from __future__ import annotations

from math import ceil, gcd, hypot, lcm, pi
from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ..design.events import AXES, derive, knots_of, pwl_moment
from ..design.timing import Raster, from_ticks, to_ticks
from ..errors import ConfigurationError, format_error

if TYPE_CHECKING:
    from collections.abc import Iterable

    from pypulseq.opts import Opts

    from ..design.events import Event

__all__ = [
    'area_until', 'ceil_raster', 'dwell_quantum', 'halve_onto', 'oblique_trapezoid',
    'require_axis', 'require_count', 'require_pair', 'require_positive', 'require_range',
    'sample_quantum', 'segment_samples', 'shift_slice',
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


def halve_onto(total_s: float, sample_s: float, raster_s: float, *, opts: Opts) -> float:
    """
    Return ``(total_s - sample_s) / 2`` **in integer ticks**, refusing a half that is not legal.

    Not a float subtraction.  ``(470 us - 409.6 us) / 2`` evaluates to 30.200000000000003 us,
    which is 302.00000000000006 ADC rasters, and pypulseq's timing check rejects it -- naming the
    block rather than the femtosecond.  This is the drift :mod:`seqcraft.design.timing` exists to
    remove, and it is the one place in either caller where the arithmetic has to leave floating
    point.

    The refusal is Rule 1's parity condition made concrete: the guard is half of a gap that has to
    land on `raster_s`, so the sampling duration must be an **even** number of them.

    **`raster_s` is an argument because the two callers pass different rasters, and which one is
    not guessable.**  :class:`~seqcraft.modules.EPI2D` halves onto the *ADC* raster, because its
    lobe duration and its sample count already force the guard onto the RF raster between them.
    :class:`~seqcraft.modules.CartesianLine`'s bipolar lobe halves onto the *RF* raster, because
    that is what pypulseq's ``_check_timing_block`` divides an ``adc`` event's ``delay`` by -- see
    :func:`dwell_quantum`.  `opts` is here to name the raster in the refusal rather than quote it
    in nanoseconds, which is the difference between a message that says what to change and one
    that says what happened.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s')
    >>> round(halve_onto(4080e-6, 4032e-6, o.rf_raster_time, opts=o) * 1e6, 3)
    24.0

    The worked example, which is why this is not a subtraction.  In floating point
    ``(470e-6 - 409.6e-6) / 2`` is 302.00000000000006 ADC rasters and pypulseq refuses it; in
    ticks it is 302 of them exactly:

    >>> (470e-6 - 409.6e-6) / 2 / o.adc_raster_time
    302.00000000000006
    >>> halve_onto(470e-6, 409.6e-6, o.adc_raster_time, opts=o) / o.adc_raster_time
    302.0

    An **odd** gap is the one this refuses, because half of it is not a delay at all:

    >>> halve_onto(470e-6, 409.7e-6, o.adc_raster_time, opts=o)
    Traceback (most recent call last):
        ...
    seqcraft.errors.ConfigurationError: the guard would be 30.1500 us, which is not a whole adc_raster_time (0.100 us).
    """
    raster_ticks = to_ticks(float(raster_s))
    gap = to_ticks(float(total_s)) - to_ticks(float(sample_s))
    if gap % (2 * raster_ticks):
        name = _raster_name(raster_s, opts)
        msg = format_error(
            f'the guard would be {from_ticks(gap) / 2 * 1e6:.4f} us, which is not a whole '
            f'{name} ({float(raster_s) * 1e6:.3f} us).',
            {'total_s': float(total_s), 'sample_s': float(sample_s), name: float(raster_s)},
            [f'the sampling duration must be an even number of {name}'],
        )
        raise ConfigurationError(msg)
    return from_ticks(gap // 2)


def dwell_quantum(num_samples: int, *, opts: Opts) -> float:
    """
    Return the coarsest step a dwell must land on for ``num_samples * dwell`` to be a legal guard.

    The fact this exists for is **not in pulseq's specification** and is two lines of pypulseq
    source: ``_check_timing_block`` divides an ``adc``, ``rf`` or ``output`` event's ``delay`` by
    ``rf_raster_time`` -- 1 us on Siemens -- and everything else by ``grad_raster_time``.  The
    guard *is* the ADC's delay, so it lives on the RF raster and not on the ten-times-finer ADC
    raster the dwell itself lives on.  Since the lobe duration is on the gradient raster, the
    whole burden falls on ``num_samples * dwell``, which must be an even number of RF rasters.

    Solved in integer ticks rather than searched.  With ``a`` the ADC raster and
    ``u = 2 * rf_raster``, a dwell of ``k*a`` works exactly when ``N*k*a`` is a multiple of ``u``,
    i.e. when ``k`` is a multiple of ``u / gcd(N*a, u)``.  Coarse, visible in the protocol table,
    and much better found here than in a timing-check failure that names a block and not a
    femtosecond.

    The **quantum** rather than the dwell, so each caller keeps its own "exactly one of bandwidth
    and dwell" refusal and its own rounding direction.  Both round *up*: rounding down raises the
    bandwidth, which shortens the lobe, which is the direction that turns a feasible protocol into
    a refusal.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s')
    >>> round(dwell_quantum(128, opts=o) * 1e9)                 # 500 ns, five ADC rasters
    500
    >>> round(dwell_quantum(320, opts=o) * 1e9)                 # 320 samples: the ADC raster
    100
    """
    adc_ticks = to_ticks(float(opts.adc_raster_time))
    pair_ticks = 2 * to_ticks(float(opts.rf_raster_time))
    return from_ticks(adc_ticks * (pair_ticks // gcd(int(num_samples) * adc_ticks, pair_ticks)))


def _raster_name(raster_s: float, opts: Opts) -> str:
    """Return the ``Opts`` field `raster_s` came from, or a spelled-out fallback."""
    for name in ('adc_raster_time', 'rf_raster_time', 'grad_raster_time', 'block_duration_raster'):
        if to_ticks(float(getattr(opts, name, -1.0))) == to_ticks(float(raster_s)):
            return name
    return 'raster'


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


def oblique_trapezoid(
    *,
    area_per_m: tuple[float, float],
    opts: Opts,
    axes: tuple[str, str] = ('x', 'y'),
    duration_s: float | None = None,
) -> tuple[Event, Event]:
    """
    Return one trapezoid of area ``hypot(ax, ay)``, split across `axes` by direction cosines.

    The amplifier's limit on a 2D gradient is on the **vector** slew, and two axes designed
    independently and stretched to a common duration each reach ``max_slew`` and together ask for
    ``sqrt(2)`` times it.  :func:`~seqcraft.events.check_limits` measures the summed waveform per
    axis and does not see that -- it reports a ``slew_norm`` warning, which is the right report
    for a phase-encode blip beside a prewinder and the wrong one for a dephaser that *is* one
    oblique vector.  Measured on the reference protocol, a naive per-axis pair asks **250.3 T/m/s
    of a 180 T/m/s amplifier at 45 degrees and 395.3 at 57.3 degrees**, compiling cleanly at both;
    this asks 177.93 at every angle, with the timing identical at every angle too.  The worst case
    is not at 45 degrees, which is the second reason to solve this rather than bound it: two
    independent trapezoids of unequal area are stretched to a common duration, and it is the
    *shorter* one that then slews hardest.

    Here rather than in :class:`~seqcraft.modules.Spiral2D` because it is what a caller builds
    when the readout does *not*: ``prephase=False`` hands back
    :meth:`~seqcraft.modules.Spiral2D.k_start_per_m` and gets out of the way, and what the caller
    then writes is a separate event -- so the vector limit becomes theirs to respect.  A spiral's
    own joins live inside its arbitrary gradient, where they are triangles rather than trapezoids
    and there is no sum for the compiler to take.

    Parameters
    ----------
    area_per_m
        The 2D area to deliver, ``(on axes[0], on axes[1])``, in 1/m.
    opts
        The scanner.  The trapezoid is designed against its limits as a single vector.
    axes
        The two logical channels, which must differ.
    duration_s
        Stretch both to this length.  ``None`` is the shortest legal one.

    Examples
    --------
    >>> import numpy as np, pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s')
    >>> gx, gy = oblique_trapezoid(area_per_m=(200.0, 0.0), opts=o)
    >>> round(float(gx.area), 6), round(float(gy.area), 6)
    (200.0, 0.0)

    Timing and vector peak do not depend on the angle, which is the property it exists for:

    >>> shapes = set()
    >>> for angle in (0.0, np.pi / 4, 1.0):
    ...     pair = oblique_trapezoid(
    ...         area_per_m=(200.0 * np.cos(angle), 200.0 * np.sin(angle)), opts=o)
    ...     shapes.add((round(pp.calc_duration(*pair) * 1e6),
    ...                 round(float(np.hypot(pair[0].amplitude, pair[1].amplitude)))))
    >>> len(shapes)
    1
    """
    if axes[0] == axes[1]:
        msg = format_error(f'axes must differ, got {axes!r}.', {'axes': axes})
        raise ConfigurationError(msg)
    components = (float(area_per_m[0]), float(area_per_m[1]))
    total = require_positive(hypot(*components), 'hypot(*area_per_m)',
                             fixes=['a zero area needs no event at all'])
    kwargs: dict[str, object] = {'channel': require_axis(axes[0]), 'area': total, 'system': opts}
    if duration_s is not None:
        kwargs['duration'] = require_positive(duration_s, 'duration_s')
    duration = float(pp.calc_duration(pp.make_trapezoid(**kwargs)))
    return tuple(
        derive(pp.scale_grad(
            pp.make_trapezoid(channel=require_axis(axis), area=total, duration=duration,
                              system=opts),
            component / total))
        for axis, component in zip(axes, components)
    )


def sample_quantum(*, dwell_s: float, opts: Opts) -> int:
    """
    Return the multiple every ADC event's sample count must be, when the count comes from a time.

    Two constraints, and only one of them is seqcraft's.  A sampling window's duration is
    ``n * dwell``, and for the block boundary that ends it to land on the **gradient raster**,
    ``n`` must be a multiple of ``grad_raster / gcd(grad_raster, dwell)``.  And pypulseq's own
    timing check enforces ``adc_samples_divisor`` -- 4 on Siemens -- which nothing in this library
    had ever met by accident, because every other readout here picks a power-of-two sample count
    for its own reasons.  A spiral picks its count from a *duration*, so it meets neither unless
    it is made to.

    The companion of :func:`dwell_quantum` and its mirror image: there the sample count is known
    and the *dwell* is snapped, here the dwell is the caller's and the *count* is snapped.

    Examples
    --------
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s')
    >>> sample_quantum(dwell_s=2.5e-6, opts=o)              # four dwells to a raster, divisor 4
    4
    >>> sample_quantum(dwell_s=3e-6, opts=o)                # 10 us / gcd(10, 3) us = 10, then 4
    20
    """
    raster_ticks = to_ticks(float(opts.grad_raster_time))
    dwell_ticks = to_ticks(require_positive(dwell_s, 'dwell_s'))
    by_raster = raster_ticks // gcd(raster_ticks, dwell_ticks)
    return lcm(by_raster, int(getattr(opts, 'adc_samples_divisor', 1) or 1))


def segment_samples(total: int, *, dwell_s: float, limit: int, opts: Opts) -> tuple[int, ...]:
    """
    Split `total` samples into the fewest ADC events under `limit`, each a whole quantum long.

    Solved rather than searched: a segment's duration is ``n * dwell``, so `n` being a multiple of
    :func:`sample_quantum` is exactly what puts every seam on the gradient raster and every count
    on ``adc_samples_divisor``.  `total` is rounded **down** onto the quantum, because rounding up
    is what would put the last sample past the gradient that positions it.

    Filled greedily to the limit with the short segment **last**, which is not cosmetic.  A seam
    is where samples are lost -- one ``adc_dead_time`` at each side, 3.28 ``dk`` at the reference
    protocol -- and a greedy fill puts it as far as possible from the start of the region.  For a
    spiral in/out pair that start is where ``|k|`` is largest, so the seam falls as far as it can
    from the echo, which is the sample the readout exists for.

    Parameters
    ----------
    total
        Samples that fit in the acquisition, before quantisation.
    dwell_s, opts
        What :func:`sample_quantum` needs.
    limit
        ``opts.adc_samples_limit``, or any interpreter's ceiling.  Zero or below is pypulseq's
        "no limit" convention and returns one segment.

    Examples
    --------
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s')
    >>> segment_samples(5286, dwell_s=2.5e-6, limit=8192, opts=o)       # fits in one
    (5284,)
    >>> segment_samples(20956, dwell_s=2.5e-6, limit=8192, opts=o)      # the short one last
    (8192, 8192, 4572)
    >>> segment_samples(20956, dwell_s=2.5e-6, limit=0, opts=o)         # pypulseq's "no limit"
    (20956,)
    """
    quantum = sample_quantum(dwell_s=dwell_s, opts=opts)
    whole = (int(total) // quantum) * quantum
    if whole <= 0:
        return ()
    if int(limit) <= 0 or whole <= int(limit):
        return (whole,)
    usable = (int(limit) // quantum) * quantum
    count = ceil(whole / usable)
    return (*(usable,) * (count - 1), whole - usable * (count - 1))
