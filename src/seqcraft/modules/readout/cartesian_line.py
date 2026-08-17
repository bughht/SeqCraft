"""
:class:`CartesianLine` -- prephaser, readout gradient and ADC designed as one thing.

``readout/`` because it contains an ADC.  The line index is deliberately *not* an argument: a
Cartesian line is identical every TR, and which line of k-space it lands on is the phase
encoding's business.  What is this module's business is where the echo falls inside it, which
nothing can measure from the tree -- a block knows times, not meanings.

The one piece of arithmetic worth reading twice
-----------------------------------------------
The readout gradient accumulates moment from the instant it starts *ramping*, so k at the echo is
the prephaser plus the ramp-up area plus whatever flat top precedes the echo::

    A_pre = -( 1/2 * G * t_rise  +  G * t_flat_before_echo )

For the symmetric case -- ADC over the flat top, equal rise and fall -- that collapses to a
formula which differs from the classic error by a single word::

    prephaser_area = -gx.flat_area / 2      # WRONG: the ramp-up is missing
    prephaser_area = -gx.area / 2           # nearly right: total area, ramps included

``pp.make_trapezoid`` supplies both fields, one letter apart in spelling and a whole ramp apart in
meaning.  Getting it wrong offsets the whole readout in k by the ramp area: the echo stops
falling at the ADC centre, and the image picks up a linear phase ramp -- or, if the shift is
large, ringing and signal loss that read as a hardware fault.

**The second one is still not exact**, and the residue is worth naming rather than absorbing:
k = 0 belongs at sample index ``matrix // 2``, whose centre time is *half a dwell* after the
midpoint of the ADC window.  That is the convention a centred inverse FFT assumes, and it is what
puts the sampled positions on ``-N/2 ... N/2-1`` times ``dk`` instead of on a half-integer grid.
So the exact answer is ``-gx.area / 2 - dk / 2``, and a formula derived from the trapezoid's
geometry alone cannot know that.

This module computes no closed form at all.  It integrates the gradient's own knots up to the
echo (:func:`~seqcraft.modules._support.area_until`), which produces the half-dwell term for
free and generalises without a second derivation to the two cases where even the corrected
closed form stops holding: **ramp sampling**, where the ADC opens during the ramp, and **partial
echo**, where the echo is not at the centre of the readout.  Hand-deriving each of those
separately is how the ramp gets dropped from one of them.
"""

from __future__ import annotations

from math import pi
from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ...design.events import derive
from ...design.logic import LogicBlock
from ...design.module import Module
from ...design.timing import EPS
from ...errors import ConfigurationError, format_error
from .._support import area_until, ceil_raster, require_axis, require_positive, require_range

if TYPE_CHECKING:
    from pypulseq.opts import Opts

    from ...design.events import Event

__all__ = ['CartesianLine']


class CartesianLine(Module):
    """
    One Cartesian readout: prephaser, readout gradient, and the ADC that samples it.

    Parameters
    ----------
    opts
        The scanner.
    fov_mm
        Field of view along `axis`.  Sets the k-space step, ``dk = 1000 / fov_mm`` in 1/m.
    matrix
        Samples across the **full** k-space extent.  With `partial_fourier` below 1 the ADC takes
        fewer than this; the extent one sample spans is unchanged, which is what keeps the image
        resolution attached to `matrix` rather than to how much of k-space was sampled.
    bandwidth_hz_px, dwell_s
        The sampling rate, either way round.  **Exactly one is required**; both or neither raises.
        ``dwell = 1 / (bandwidth_hz_px * matrix)``, then snapped onto the ADC raster, so the
        achieved values are readable back as :attr:`dwell_s` and :attr:`bandwidth_hz_px`.
    partial_fourier
        Fraction of the full k-space extent that is sampled, with the missing part taken off the
        **pre-echo** side.  ``1.0`` is a symmetric full echo.  ``0 < partial_fourier <= 1``.
    prephaser_duration_s
        Lengthen the prephaser beyond its own minimum.  ``None`` is the minimum.  A composite
        passes the winder maximum here; see :class:`~seqcraft.modules.GRE2DTR`.
    axis
        Logical gradient channel.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
    ...          adc_dead_time=10e-6)
    >>> ro = CartesianLine(opts=o, fov_mm=250.0, matrix=128, bandwidth_hz_px=250.0)
    >>> ro.num_samples, round(ro.dwell_s * 1e9)
    (128, 31200)

    The prephaser cancels everything the readout has accumulated by the echo, ramp included:

    >>> gx = ro.gx
    >>> from seqcraft.modules._support import area_until
    >>> accumulated = area_until(gx, ro.time_to_echo() - ro.prephaser_duration_s)
    >>> abs(ro.prephaser_area_per_m + accumulated) < 1e-9
    True

    Partial echo moves the echo forward and shortens the ADC, in one factor with no new branch:

    >>> pf = CartesianLine(opts=o, fov_mm=250.0, matrix=128, bandwidth_hz_px=250.0,
    ...                    partial_fourier=0.75)
    >>> pf.num_samples
    96
    >>> round(pf.pre_echo_samples / pf.num_samples, 3)          # (2*pf - 1) / (2*pf) = 1/3
    0.333

    Notes
    -----
    **`acquire=False` plays the same gradients with no ADC**, for dummy repetitions.  The
    duration is unchanged, which is the point: a dummy has to load the gradients exactly as a real
    TR does, or the steady state it is establishing is not the one that will be acquired.

    **`phase_deg` is the receiver phase, and it is not optional in a spoiled sequence.**  The
    receiver is phase-locked to the transmitter: whatever carrier phase the excitation was given,
    the ADC has to demodulate at.  Leaving it at zero while the RF runs a spoiling schedule
    writes the schedule's quadratic phase straight into ``ky`` -- measured on a single off-centre
    voxel, that scattered a point source over the whole phase-encode direction and moved its peak
    by thirteen pixels, with the readout direction perfectly correct beside it.  The composite
    that owns the schedule passes the same number to both; this module only applies it.

    **`offset_mm` shifts the FOV along the readout axis** by demodulating at a shifted frequency,
    ``adc.freq_offset = gx.amplitude * offset_m``, with the phase referenced to the echo so the
    shift contributes no phase at k = 0.  It is the receive-side twin of the slice offset an
    :class:`~seqcraft.modules.Excitation` applies to its RF, and the sign is a trap in both:
    against a symmetric phantom a mirrored FOV looks entirely correct.  Note that MRzero's
    ``Sequence.import_file`` **parses ``freq`` and then ignores it**, so a simulation through a
    written ``.seq`` cannot see this offset at all -- checking it needs a scanner or a simulator
    that demodulates.
    """

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: float,
        matrix: int,
        bandwidth_hz_px: float | None = None,
        dwell_s: float | None = None,
        partial_fourier: float = 1.0,
        prephaser_duration_s: float | None = None,
        axis: str = 'x',
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.fov_mm = require_positive(fov_mm, 'fov_mm')
        self.matrix = int(matrix)
        self.axis = require_axis(axis)
        self.partial_fourier = require_range(partial_fourier, 'partial_fourier', low=0.0, high=1.0)
        if self.matrix < 1:
            msg = format_error(f'matrix must be at least 1, got {self.matrix}.',
                               {'matrix': self.matrix})
            raise ConfigurationError(msg)

        #: k-space step between adjacent samples, 1/m.
        self.dk_per_m = 1e3 / self.fov_mm
        self.dwell_s = self._resolve_dwell(bandwidth_hz_px, dwell_s)
        #: Achieved receive bandwidth per pixel, Hz.  Reads back the snapped dwell.
        self.bandwidth_hz_px = 1.0 / (self.dwell_s * self.matrix)

        #: Samples the ADC actually takes: the full matrix, less the missing pre-echo part.
        self.num_samples = max(1, int(round(self.partial_fourier * self.matrix)))
        #: Samples before the echo.  ``(2*pf - 1) / (2*pf)`` of the way into the readout.
        self.pre_echo_samples = self.num_samples - (self.matrix - self.matrix // 2)
        if self.pre_echo_samples < 0:
            msg = format_error(
                f'partial_fourier = {self.partial_fourier:g} leaves no pre-echo samples at '
                f'matrix = {self.matrix}.',
                {'partial_fourier': self.partial_fourier, 'matrix': self.matrix,
                 'num_samples': self.num_samples},
                ['partial_fourier below ~0.55 has too little pre-echo data for a phase estimate'],
            )
            raise ConfigurationError(msg)

        self.gx, self.adc = self._design_readout()
        #: Seconds from the readout gradient's own start to k = 0.
        self._echo_in_gx = float(self.adc.delay) + (self.pre_echo_samples + 0.5) * self.dwell_s
        self._prephaser = self._design_prephaser(prephaser_duration_s)

    # ------------------------------------------------------------------ what it knows
    @property
    def prephaser_duration_s(self) -> float:
        """Seconds the prephaser occupies -- its own minimum unless one was requested."""
        return float(pp.calc_duration(self._prephaser))

    @property
    def prephaser_area_per_m(self) -> float:
        """
        The prephaser's area, 1/m -- negative, and exactly minus what the readout accumulates.

        Exposed so a test can assert that identity without reaching into the block.
        """
        return float(self._prephaser.area)

    def time_to_echo(self) -> float:
        """
        Seconds from the start of this module's block to k = 0.

        Not measurable from the tree: the block knows when its events play, not which instant
        among them is the echo.
        """
        return self.prephaser_duration_s + self._echo_in_gx

    # ----------------------------------------------------------------------- assembly
    def build(
        self, *, acquire: bool = True, phase_deg: float = 0.0, offset_mm: float = 0.0,
    ) -> LogicBlock:
        """
        Return prephaser, readout gradient and -- unless `acquire` is false -- the ADC.

        Parameters
        ----------
        acquire
            ``False`` plays identical gradients with no ADC, for a dummy repetition.
        phase_deg
            Receiver phase for this repetition, degrees.  **Give it the same value the
            excitation got**: the receiver is phase-locked to the transmitter, so an RF-spoiling
            schedule that moves one has to move the other.  See the note below.
        offset_mm
            Shift the FOV along `axis`, by demodulating at a shifted frequency.
        """
        start = self.prephaser_duration_s
        out = LogicBlock().add(0.0, self._prephaser).add(start, self.gx)
        if acquire:
            out.add(start, self._adc_for(float(offset_mm) / 1e3, float(np.deg2rad(phase_deg))))
        return out

    def _adc_for(self, offset_m: float, phase_rad: float) -> Event:
        """Return the ADC with the receiver phase and the FOV shift this call asked for."""
        if offset_m == 0.0 and phase_rad == 0.0:
            return self.adc
        freq_offset = float(self.gx.amplitude) * offset_m
        # Referenced to the echo rather than to the first sample, so the shift adds no phase at
        # k = 0.  Omitting this term is a linear phase across the image that reads as a gradient
        # delay -- the same correction _shift_slice applies at an RF pulse's effective centre.
        phase_offset = phase_rad - 2 * pi * freq_offset * (
            self._echo_in_gx - float(self.adc.delay)
        )
        return derive(self.adc, freq_offset=freq_offset, phase_offset=phase_offset)

    # -------------------------------------------------------------------------- design
    def _resolve_dwell(self, bandwidth_hz_px: float | None, dwell_s: float | None) -> float:
        """Return the dwell time on the ADC raster, from whichever of the two was given."""
        if (bandwidth_hz_px is None) == (dwell_s is None):
            msg = format_error(
                'exactly one of bandwidth_hz_px and dwell_s is required.',
                {'bandwidth_hz_px': bandwidth_hz_px, 'dwell_s': dwell_s},
                [
                    'bandwidth_hz_px=200 is how a protocol quotes it',
                    'dwell_s=5e-6 is how the receiver is actually programmed',
                    'they are one number: dwell = 1 / (bandwidth_hz_px * matrix)',
                ],
            )
            raise ConfigurationError(msg)
        wanted = (
            float(dwell_s) if dwell_s is not None
            else 1.0 / (require_positive(bandwidth_hz_px, 'bandwidth_hz_px') * self.matrix)
        )
        raster = float(self.opts.adc_raster_time)
        snapped = round(require_positive(wanted, 'dwell_s') / raster) * raster
        if snapped <= 0.0:
            msg = format_error(
                f'a dwell of {wanted * 1e9:.3f} ns rounds to zero on the '
                f'{raster * 1e9:.0f} ns ADC raster.',
                {'dwell_s': wanted, 'adc_raster_time': raster},
                ['lower bandwidth_hz_px, or reduce matrix'],
            )
            raise ConfigurationError(msg)
        return snapped

    def _design_readout(self) -> tuple[Event, Event]:
        """Return the readout gradient and its ADC, sampling the flat top."""
        amplitude = self.dk_per_m / self.dwell_s
        if abs(amplitude) > self.opts.max_grad:
            need = self.dk_per_m / self.opts.max_grad
            msg = format_error(
                f'the readout needs {amplitude:.0f} Hz/m, above the '
                f'{self.opts.max_grad:.0f} Hz/m limit.',
                {'fov_mm': self.fov_mm, 'dwell_s': self.dwell_s,
                 'bandwidth_hz_px': self.bandwidth_hz_px},
                [
                    f'lower the bandwidth to at most {1.0 / (need * self.matrix):.0f} Hz/px',
                    'or widen fov_mm, which lowers dk and so lowers the amplitude',
                ],
            )
            raise ConfigurationError(msg)
        # The flat top is the ADC window rounded up onto the gradient raster; the few spare
        # microseconds sit after the last sample, where they change no k-space position.
        flat_time = ceil_raster(self.num_samples * self.dwell_s, self.opts.grad_raster_time)
        gx = pp.make_trapezoid(
            channel=self.axis, amplitude=amplitude, flat_time=flat_time, system=self.opts,
        )
        # The sampling offset goes in the ADC's own delay, and the node goes where the gradient
        # starts: pp.make_adc raises a delay below adc_dead_time up to it, and seqcraft preserves
        # an event's delay rather than folding it away, so a node offset would add to it twice.
        adc = pp.make_adc(
            num_samples=self.num_samples, dwell=self.dwell_s, delay=gx.rise_time, system=self.opts,
        )
        # make_adc raises a delay below adc_dead_time up to it, so on a scanner whose receive
        # dead time exceeds the ramp the window starts inside the flat top and would otherwise
        # run off the end of it -- sampling the fall ramp, which is not what the k arithmetic
        # here assumes.  Lengthen the flat top instead; the ramp time is unaffected.
        overrun = (
            float(adc.delay) - float(gx.rise_time) + self.num_samples * self.dwell_s - flat_time
        )
        if overrun > 0.0:
            gx = pp.make_trapezoid(
                channel=self.axis,
                amplitude=amplitude,
                flat_time=ceil_raster(flat_time + overrun, self.opts.grad_raster_time),
                system=self.opts,
            )
        return gx, adc

    def _design_prephaser(self, requested_s: float | None) -> Event:
        """Return the prephaser that puts k = 0 at the echo, stretched if one was requested."""
        area = -area_until(self.gx, self._echo_in_gx)
        shortest = pp.make_trapezoid(channel=self.axis, area=area, system=self.opts)
        self._min_prephaser_duration_s = float(pp.calc_duration(shortest))
        if requested_s is None:
            return shortest
        wanted = ceil_raster(
            require_positive(requested_s, 'prephaser_duration_s'), self.opts.grad_raster_time,
        )
        # EPS, not exact: a caller passing `max(ro.prephaser_duration_s, pe.min_duration_s)`
        # back in hands us its own minimum, and snapping that onto the raster can land one ulp
        # below where it came from.  Times a nanosecond apart are the same time.
        if wanted < self._min_prephaser_duration_s - EPS:
            msg = format_error(
                f'prephaser_duration_s = {requested_s * 1e6:.1f} us is shorter than the '
                f'minimum this readout needs.',
                {'prephaser_duration_s': requested_s,
                 'min_prephaser_duration_s': self._min_prephaser_duration_s,
                 'prephaser_area_per_m': area},
                [
                    f'pass prephaser_duration_s >= {self._min_prephaser_duration_s:.6g}',
                    'or pass None for the shortest legal prephaser',
                ],
            )
            raise ConfigurationError(msg)
        return pp.make_trapezoid(
            channel=self.axis, area=area, duration=wanted, system=self.opts,
        )
