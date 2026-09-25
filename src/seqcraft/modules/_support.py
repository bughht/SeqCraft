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

from math import gcd, pi
from typing import TYPE_CHECKING, Any

import numpy as np
import pypulseq as pp

from ..design.events import derive, knots_of, pwl_moment
from ..design.timing import Raster, from_ticks, to_ticks
from ..design.validation import (
    require_axis,
    require_count,
    require_pair,
    require_positive,
    require_range,
)
from ..errors import ConfigurationError, format_error

if TYPE_CHECKING:
    from collections.abc import Iterable

    from pypulseq.opts import Opts

    from ..design.events import Event

__all__ = [
    'area_until', 'ceil_raster', 'check_peak_b1', 'duration_for_peak_b1', 'duration_remedy', 'dwell_quantum',
    'halve_onto', 'peak_b1_hz', 'peak_scales_with_duration', 'require_axis',
    'require_usable_max_b1', 'require_count', 'require_pair',
    'require_positive', 'require_range', 'shift_slice',
]


def peak_b1_hz(rf: Event) -> float:
    """Peak :math:`|B_1|` of an RF event, in hertz -- ``max(abs(rf.signal))`` and nothing else."""
    return float(np.abs(np.asarray(rf.signal)).max())


def check_peak_b1(rf: Event, opts: Opts, *, described: str, remedies: Iterable[str]) -> None:
    """
    Refuse an RF event whose peak amplitude exceeds ``opts.max_b1``.

    The **shared invariant** every RF-producing module holds, and it is literal::

        peak_b1_hz(rf) <= opts.max_b1

    measured on the waveform that will be emitted rather than predicted from the pulse's
    parameters.  There is no sentinel: ``+inf`` passes the comparison on its own, which is how a
    caller designs without a transmit limit, and **it is the only value that does**.

    ============  ==================================================
    ``max_b1``    behaviour
    ============  ==================================================
    ``+inf``      enforcement disabled
    finite > 0    the literal comparison above
    0, negative   refused as an unusable limit
    ``NaN``       refused -- every comparison against it is False
    ``None``      refused; it is not a spelling of "unlimited"
    ============  ==================================================

    Zero is not pypulseq's convention for this field: unlike ``adc_samples_limit``, whose ``0`` is
    documented as no limit, ``make_sinc_pulse`` compares ``rf_amplitude > system.max_b1``
    unconditionally, so a zero limit makes it warn at ``inf %``.  And ``pp.Opts(max_b1=None)``
    falls back to the 20 uT default rather than to no limit, so ``None`` cannot mean one here.

    Note that ``pp.Opts()`` **defaults ``max_b1`` to 851.52 Hz (20 uT)**, so the limit is live
    unless a caller has deliberately raised it.

    What differs between modules is not the measurement but the **remedy**, so that is the
    caller's to supply: lengthening is right for a sinc whose envelope merely stretches, and
    meaningless for an adiabatic inversion whose peak is set by its frequency sweep.

    Parameters
    ----------
    rf
        The RF event to measure.
    opts
        Supplies ``max_b1`` in hertz, and ``gamma`` for reporting the same numbers in microtesla.
    described
        How to name the pulse in the error, e.g. ``"a 2 ms 180 degree sinc pulse"``.
    remedies
        Module-specific fixes, in the order a caller should try them.

    Raises
    ------
    ConfigurationError
        Naming the measured peak as a percentage of the limit, and the remedies given.
    """
    configured = getattr(opts, 'max_b1', None)
    peak = peak_b1_hz(rf)
    # One comparison decides both questions.  ``limit > 0.0`` is True for a finite positive limit
    # and for ``+inf``, and False for zero, negatives and NaN alike; ``peak <= limit`` is then
    # trivially true for ``+inf``, which is how "design without a transmit limit" is spelled.
    # An absent or ``None`` limit is **not** a third spelling of that: ``pp.Opts(max_b1=None)``
    # falls back to the 20 uT default rather than to no limit, so nothing here may invent one.
    if configured is not None:
        limit = float(configured)
        if limit > 0.0 and peak <= limit:
            return
    else:
        limit = float('nan')
    # In hertz throughout; microtesla is for the reader, so it follows this scanner's nucleus
    # rather than a hard-coded proton value.
    gamma = abs(float(getattr(opts, 'gamma', 0.0) or 0.0))
    detail: dict[str, float] = {'peak_b1_hz': round(peak, 2), 'max_b1_hz': limit}
    if gamma > 0.0:
        detail['peak_b1_uT'] = round(peak / gamma * 1e6, 2)
        detail['max_b1_uT'] = round(limit / gamma * 1e6, 2) if np.isfinite(limit) else limit
    if not limit > 0.0:
        _raise_unusable_max_b1(configured, limit, described, detail)
    raise ConfigurationError(format_error(
        f'{described} peaks at {peak / limit * 100:.0f} % of max_b1.',
        detail,
        [*remedies, 'or raise max_b1, if the scanner really does deliver it'],
    ))


def _raise_unusable_max_b1(configured: object, limit: float, described: str,
                           detail: dict[str, float] | None = None) -> None:
    """Refuse a ``max_b1`` that is not a transmit limit: zero, negative, NaN or absent."""
    stated = 'max_b1 is not set' if configured is None else f'max_b1 is {limit:g}'
    raise ConfigurationError(format_error(
        f'{stated}, which is not a transmit limit {described} can meet.',
        detail if detail is not None else {'max_b1': repr(configured)},
        ['set a real max_b1 -- sc.convert(20, "uT", "Hz") is a common value',
         'or use max_b1=inf to design without a transmit limit'],
    ))


def require_usable_max_b1(opts: Opts, *, described: str) -> None:
    """
    Refuse an unusable ``opts.max_b1`` **before** a pulse is built from it.

    :func:`check_peak_b1` reaches the same verdict, but it needs a waveform to measure -- and
    pypulseq's shaped factories compare ``rf_amplitude > system.max_b1`` while designing one, so
    an absent limit raises ``TypeError`` from inside pypulseq before this package ever sees it.
    Calling this first turns that into the same :class:`~seqcraft.errors.ConfigurationError` every
    other unusable limit produces.

    ``+inf`` and any finite positive value pass.
    """
    configured = getattr(opts, 'max_b1', None)
    if configured is None:
        _raise_unusable_max_b1(None, float('nan'), described)
    limit = float(configured)
    if not limit > 0.0:
        _raise_unusable_max_b1(configured, limit, described)


def duration_for_peak_b1(duration_s: float, peak_hz: float, max_b1_hz: float) -> float:
    """
    Return the duration at which a pulse of **unchanged envelope** would fit ``max_b1``, seconds.

    Peak :math:`B_1` scales as ``1 / duration`` for a waveform that is merely stretched -- the
    shape and the flip angle held fixed -- so this inverts that relation.  Rounded up onto a tenth
    of a millisecond, because a floor quoted to the nanosecond is refused again by its own last
    digit.

    **It is only a floor where that premise holds.**  Changing the duration re-runs the pulse
    design, and several designs are not a simple stretch: an SLR filter is recomputed, and a
    module that derives ``time_bw_product`` from a fixed bandwidth changes the design by changing
    the duration.  It does not hold at all for adiabatic pulses, whose peak is set by the
    frequency sweep -- a hyperbolic secant's peak does not move with duration.

    Use :func:`duration_remedy` rather than this directly, so that the wording matches the
    strength of the claim.
    """
    return float(np.ceil(duration_s * peak_hz / max_b1_hz * 1e4) / 10.0) / 1e3


#: Shapes whose envelope is merely stretched by a longer duration *when the time--bandwidth
#: product is held fixed*, so peak B1 then scales as ``1 / duration`` exactly.  An SLR filter is
#: recomputed instead, so it never qualifies.
_STRETCHED_SHAPES = frozenset({'sinc', 'gauss'})


def peak_scales_with_duration(pulse: str, design_opts: dict[str, Any]) -> bool:
    """
    Whether lengthening this pulse *stretches* it, so that peak B1 scales as ``1 / duration``.

    Two things have to hold, and the second is why the pulse's name is not enough.  The shape has
    to be one that stretches at all -- ``'sinc'`` or ``'gauss'``, not an SLR filter, which is
    recomputed.  And the design has to be pinned by a **time--bandwidth product** rather than by a
    **bandwidth**: ``make_gauss_pulse`` takes an explicit ``bandwidth``, and with one supplied a
    longer duration buys almost nothing.  Measured, at a 4 kHz bandwidth and a 90 degree flip, the
    peak is 1000.0 Hz at 1 ms and 1000.0 Hz at 2 ms.

    So a caller who passes ``pulse_opts={'bandwidth': ...}`` must not be told a duration floor,
    because rebuilding there would still be over the limit.
    """
    return pulse in _STRETCHED_SHAPES and 'bandwidth' not in design_opts


def duration_remedy(duration_s: float, peak_hz: float, max_b1_hz: float, *,
                    exact: bool) -> str:
    """
    Return the "make it longer" remedy line, worded to match what is actually guaranteed.

    `exact` says whether the envelope is unchanged by the duration, so that the inverted
    ``1 / duration`` relation is a **floor** rather than a starting point.  It is true for a sinc
    or gauss at a fixed time--bandwidth product, and false wherever changing the duration
    redesigns the pulse.
    """
    floor_ms = duration_for_peak_b1(duration_s, peak_hz, max_b1_hz) * 1e3
    if exact:
        return f'pass duration_s >= {floor_ms:.1f} ms, which is where this shape fits'
    return (f'try duration_s >= {floor_ms:.1f} ms and check again -- this pulse is redesigned '
            f'when the duration changes, so that is a starting point rather than a floor')



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
