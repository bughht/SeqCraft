"""
RF pulse modules: excitation, refocusing, inversion and saturation.

Every pulse here shares one piece of arithmetic that is worth stating once, because it is the
most commonly omitted line in hand-written pulseq code.  Selecting an off-centre slice needs
**two** terms:

.. code-block:: text

    freq_offset  = G_slice * z                       # which slice
    phase_offset = base + extra - 2*pi*f*rf.center   # keeps that slice's phase consistent

The second compensates for the RF centre not coinciding with the start of its block.  Omit it and
every slice acquires a different phase, which shows up in the reconstruction as slice-dependent
phase and usually gets blamed on the scanner.  :meth:`RFPulse.offset_rf` holds it, so every
module in the package gets it right or none do -- and it *composes* with the existing phase
rather than overwriting it, so a CPMG pi/2 on a refocusing pulse survives.

The other shared piece is :attr:`RFPulse.isodelay`: the time from the start of the built block to
the pulse's effective centre.  That is the quantity TE is measured from, and reading it from
``rf.center`` rather than recomputing a centre of mass means an asymmetric or minimum-phase pulse
needs no special-casing.

Examples
--------
>>> import seqcraft as sc
>>> system = sc.System.preset('generic_3t')
>>> exc = SincExcitation(system, flip_deg=15, duration_us=1000, slice_thickness_mm=5)
>>> exc.build()
LogicBlock(exc, 3 nodes, 1.80 ms)
>>> round(exc.isodelay * 1e6)          # 100 us dead time + 500 us to the sinc centre
630
>>> round(exc.bandwidth_hz)            # time-bandwidth 4 over 1 ms
4000
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any

import numpy as np
import pypulseq as pp

from ...core import events as ev
from ...core.errors import ConfigurationError, format_error
from ...core.logic import LogicBlock
from ...core.module import Module
from ...core.raster import ceil_to
from ...core.registry import register
from ...core.validate import Range, require_in, require_in_range, require_positive

if TYPE_CHECKING:
    from types import SimpleNamespace

    from ...core.system import System

__all__ = [
    'AdiabaticInversion',
    'GaussSaturation',
    'HardExcitation',
    'HardRefocusing',
    'RFPulse',
    'RefocusingPulse',
    'SLRExcitation',
    'SLRRefocusing',
    'SincExcitation',
    'SincRefocusing',
    'SlabExcitation',
]

_FLIP_RANGE = Range(-720.0, 720.0, 'deg', ((57.29578, 'rad'),))
_DURATION_RANGE = Range(1.0, 100_000.0, 'us', ((1e6, 's'), (1e3, 'ms')))
_THICKNESS_RANGE = Range(0.05, 500.0, 'mm', ((1e3, 'm'),))


class RFPulse(Module):
    """
    Shared behaviour for RF pulse modules.

    Parameters
    ----------
    system
        The scanner.
    flip_deg
        Flip angle, degrees.
    duration_us
        Pulse duration, microseconds.
    slice_thickness_mm
        Slice or slab thickness.  ``None`` makes the pulse non-selective, in which case there
        is no slice-select gradient and ``slice_offset_m`` has no effect.
    rephase
        Include the slice rephasing lobe.  Excitation pulses need it; a refocusing pulse does
        not, because its slice-select gradient is symmetric about the pulse and crushers take
        the rephaser's place.
    time_bw_product
        Time-bandwidth product; sets slice profile sharpness for a given duration.
    apodization
        Hamming apodization fraction for windowed pulses.
    phase_offset_rad
        Design-time carrier phase, *composed* with ``build``'s ``rf_phase_rad``.  A CPMG
        refocusing pulse sets ``pi/2`` here.
    freq_offset_hz
        Design-time carrier frequency offset, for chemical-shift-selective pulses such as fat
        saturation.
    slice_axis
        Which logical axis carries the slice-select gradient.
    passband_ripple, stopband_ripple
        SLR design tolerances -- how much the profile may deviate inside the slice and how much
        signal may leak outside it.  Ignored by the sinc pulses, which have no design freedom: a
        sinc's ripple is whatever the truncation and apodization leave.
    filter_type
        SLR filter design: ``'ms'`` (maximally flat), ``'pm'`` (Parks--McClellan equiripple),
        ``'min'`` or ``'max'`` phase.  A minimum-phase pulse puts its energy late, which shortens
        the effective TE -- and is why :attr:`isodelay` reads ``rf.center`` rather than assuming
        the middle.
    use
        The pypulseq RF ``use`` string.  Fixed by each subclass, so ``'undefined'`` -- and the
        ``detect_rf_use`` crutch that exists because people leave it undefined -- is
        unreachable.

    Properties
    ----------
    isodelay
        Seconds from the start of the built block to the pulse's effective centre.  **This is
        what TE is measured from.**
    pulse_duration
        Seconds occupied by the pulse and its slice-select gradient, excluding any rephaser.
    duration
        Seconds occupied by the whole built block, rephaser included.
    bandwidth_hz, peak_b1_uT, energy_hz2_s
        Reported properties of the designed pulse.

    Build arguments
    ---------------
    slice_offset_m : float, default 0.0
        Slice position along `slice_axis`, in metres.
    rf_phase_rad : float, default 0.0
        Extra carrier phase, e.g. from :func:`~seqcraft.core.ordering.rf_spoil_phase`.  Added
        to `phase_offset_rad`.

    `rephase` is deliberately **not** a build argument: dropping the rephaser changes the block's
    duration, which would silently invalidate every timing property the caller places the block
    by.  Construct with ``rephase=False`` and place :meth:`rephaser` yourself instead -- then the
    properties and the block always agree.  That is also the arrangement you want when the
    rephaser is to overlap a phase encode and a readout prephaser.

    Notes
    -----
    Peak B1 for a hard pulse is ``flip / (2*pi*gamma*duration)``, so a short high-flip pulse can
    exceed ``max_b1`` on real hardware -- a 2 ms 180 degree slice-selective sinc needs about
    23 uT, above the 20 uT typical limit.  When pypulseq rejects a pulse it reports only
    ``Amplitude violation``; this module raises first, naming the three parameters that fix it.
    """

    #: The pypulseq RF ``use`` string.  Set by each subclass.
    use = 'undefined'

    def __init__(
        self,
        system: System,
        *,
        flip_deg: float,
        duration_us: float,
        slice_thickness_mm: float | None = None,
        rephase: bool = True,
        time_bw_product: float = 4.0,
        apodization: float = 0.5,
        phase_offset_rad: float = 0.0,
        freq_offset_hz: float = 0.0,
        slice_axis: str = 'z',
        passband_ripple: float = 0.01,
        stopband_ripple: float = 0.01,
        filter_type: str = 'ms',
        regime: str = 'default',
    ) -> None:
        super().__init__(system, regime=regime)
        self.passband_ripple = float(passband_ripple)
        self.stopband_ripple = float(stopband_ripple)
        self.filter_type = str(filter_type)
        self.flip_deg = float(flip_deg)
        self.duration_us = float(duration_us)
        self.slice_thickness_mm = None if slice_thickness_mm is None else float(slice_thickness_mm)
        self.rephase = bool(rephase)
        self.time_bw_product = float(time_bw_product)
        self.apodization = float(apodization)
        self.phase_offset_rad = float(phase_offset_rad)
        self.freq_offset_hz = float(freq_offset_hz)
        self.slice_axis = str(slice_axis)

        require_in_range(self, 'flip_deg', _FLIP_RANGE.lo, _FLIP_RANGE.hi, unit='deg')
        require_in_range(self, 'duration_us', _DURATION_RANGE.lo, _DURATION_RANGE.hi, unit='us')
        require_positive(self, 'duration_us')
        require_in(self, 'slice_axis', ('x', 'y', 'z'))
        if self.slice_thickness_mm is not None:
            require_in_range(
                self, 'slice_thickness_mm', _THICKNESS_RANGE.lo, _THICKNESS_RANGE.hi, unit='mm'
            )
        if self.slice_thickness_mm is None and self.rephase:
            msg = format_error(
                f'{type(self).__name__}: rephase=True needs a slice_thickness_mm.',
                {'reason': 'there is no slice-select gradient to rephase'},
                ['pass slice_thickness_mm=..., or rephase=False for a non-selective pulse'],
            )
            raise ConfigurationError(msg)

        self._check_b1()
        self.rf, self.gz, self.gzr = self._design()
        self._check_designed_b1()

    # ------------------------------------------------------------------ what a subclass does
    def _design(self) -> tuple[SimpleNamespace, SimpleNamespace | None, SimpleNamespace | None]:
        """
        Build the RF event and, if slice-selective, the slice-select and rephasing lobes.

        Returns
        -------
        rf, gz, gzr
            `gz` and `gzr` are ``None`` for a non-selective pulse; `gzr` is ``None`` when
            ``rephase`` is false.
        """
        raise NotImplementedError

    def _check_designed_b1(self) -> None:
        """
        Raise if the **designed** waveform's peak exceeds ``max_b1``.

        :meth:`_check_b1` runs before designing and can only use the hard-pulse floor, which is a
        lower bound -- for a 180 in 4 ms it is 2.9 uT against a 20 uT limit, so it always passes.
        That makes it a necessary condition and nothing more.  The peak of the shaped pulse that
        actually gets played is two to five times the floor, and until this check existed nothing
        compared the two: a pulse that no amplifier could produce compiled, wrote and got as far as
        the scanner, which refused it with ``RF amplitude could not be realized``.

        ``max_b1`` has to come from the site.  It is not a property of the pulse or of the gradient
        system but of the transmit chain and the loading, so the default of 20 uT is a placeholder --
        set ``max_b1_uT`` on the :class:`~seqcraft.core.system.System` from the reference voltage the
        scanner reports, and keep a margin.
        """
        limit_t = float(self.opts.max_b1) / self.system.gamma
        peak_t = self.peak_b1_uT / 1e6
        if peak_t > limit_t * 1.001:
            # Peak B1 falls roughly as 1/duration at fixed flip and time-bandwidth product.
            suggested = self.duration_us * peak_t / limit_t
            msg = format_error(
                f'the designed {self.flip_deg:g} degree pulse peaks at {peak_t * 1e6:.2f} uT, '
                f'above the {limit_t * 1e6:.2f} uT limit.',
                {
                    'duration_us': self.duration_us,
                    'time_bw_product': self.time_bw_product,
                    'peak_b1_uT': round(peak_t * 1e6, 2),
                    'max_b1_uT': round(limit_t * 1e6, 2),
                },
                [
                    f'lengthen the pulse to about {suggested:.0f} us -- peak B1 falls as 1/duration',
                    'reduce time_bw_product, which lowers the peak at the cost of profile sharpness',
                    'raise max_b1_uT if the scanner really does allow it',
                ],
            )
            raise ConfigurationError(msg)

    def _check_b1(self) -> None:
        """
        Raise before designing if the requested flip cannot be reached within ``max_b1``.

        pypulseq's own message is ``Amplitude violation (117%)``, which does not say which
        parameter to change.  This one names all three.
        """
        duration_s = self.duration_us / 1e6
        gamma = self.system.gamma
        # A shaped pulse needs more peak B1 than a hard one of the same flip; the sinc factor
        # is ~2.2 for tbw 4 with Hamming apodization.  Use the hard-pulse floor, which is a
        # lower bound, so this never rejects something that would in fact fit.
        floor_t = abs(math.radians(self.flip_deg)) / (2.0 * math.pi * gamma * duration_s)
        limit_t = float(self.opts.max_b1) / gamma
        if floor_t > limit_t:
            needed_us = abs(math.radians(self.flip_deg)) / (2.0 * math.pi * gamma * limit_t) * 1e6
            msg = format_error(
                f'a {self.flip_deg:g} degree pulse in {self.duration_us:g} us needs at least '
                f'{floor_t * 1e6:.1f} uT, above the {limit_t * 1e6:.1f} uT limit.',
                {
                    'flip_deg': self.flip_deg,
                    'duration_us': self.duration_us,
                    'max_b1_uT': limit_t * 1e6,
                },
                [
                    f'lengthen the pulse to at least {needed_us:.0f} us '
                    '(a shaped pulse needs roughly twice that)',
                    'reduce flip_deg',
                    'use an adiabatic pulse, which reaches a full inversion at lower peak B1',
                ],
            )
            raise ConfigurationError(msg)

    # -------------------------------------------------------------------------- properties
    @property
    def flip_rad(self) -> float:
        """Flip angle in radians."""
        return math.radians(self.flip_deg)

    @property
    def isodelay(self) -> float:
        """
        Seconds from the start of the built block to the pulse's effective centre.

        **The quantity TE is measured from.**  Taken from ``rf.center``, which pypulseq 1.5
        carries as a field, rather than recomputing a centre of mass -- the two disagree for
        asymmetric and minimum-phase pulses, and the designed value is the correct one for
        timing.
        """
        return float(self.rf.delay) + _rf_center(self.rf)

    @property
    def pulse_duration(self) -> float:
        """Seconds occupied by the pulse and its slice-select gradient, excluding the rephaser."""
        events = (self.rf,) if self.gz is None else (self.rf, self.gz)
        return ceil_to(float(pp.calc_duration(*events)), self.system.block_raster_s)

    @property
    def rephase_duration(self) -> float:
        """Seconds occupied by the slice rephasing lobe, or zero when there is none."""
        if self.gzr is None:
            return 0.0
        return ceil_to(float(pp.calc_duration(self.gzr)), self.system.block_raster_s)

    @property
    def duration(self) -> float:
        """Seconds occupied by the whole built block, including the rephaser if played."""
        return self.pulse_duration + (self.rephase_duration if self.rephase else 0.0)

    @property
    def slice_select_amplitude(self) -> float:
        """Slice-select gradient amplitude in Hz/m, or zero for a non-selective pulse."""
        return 0.0 if self.gz is None else float(self.gz.amplitude)

    @property
    def bandwidth_hz(self) -> float:
        """Excitation bandwidth in hertz, from ``time_bw_product / duration``."""
        return self.time_bw_product / (self.duration_us / 1e6)

    @property
    def peak_b1_uT(self) -> float:
        """Peak B1 amplitude, microtesla."""
        signal = np.asarray(self.rf.signal)
        peak_hz = float(np.max(np.abs(signal))) if signal.size else 0.0
        return peak_hz / self.system.gamma * 1e6

    @property
    def energy_hz2_s(self) -> float:
        """
        Integrated ``|B1|^2 dt`` in Hz^2 s -- the quantity SAR is proportional to.

        Reported rather than limited: an absolute SAR figure needs a coil and a body model,
        which is outside seqcraft's scope.
        """
        signal = np.asarray(self.rf.signal)
        return float(np.sum(np.abs(signal) ** 2) * self.system.rf_raster_s)

    # -------------------------------------------------------------------------------- build
    def build(self, *, slice_offset_m: float = 0.0, rf_phase_rad: float = 0.0) -> LogicBlock:
        """
        Return the pulse, its slice-select gradient, and the rephasing lobe if it owns one.

        Parameters
        ----------
        slice_offset_m
            Slice position along the slice axis, metres.
        rf_phase_rad
            Extra carrier phase, added to `phase_offset_rad`.

        Examples
        --------
        >>> import seqcraft as sc
        >>> system = sc.System.preset('generic_3t')
        >>> exc = SincExcitation(system, flip_deg=15, duration_us=1000, slice_thickness_mm=5)
        >>> exc.build()
        LogicBlock(exc, 3 nodes, 1.80 ms)
        >>> bare = SincExcitation(system, flip_deg=15, duration_us=1000,
        ...                       slice_thickness_mm=5, rephase=False)
        >>> bare.build(), round(bare.duration * 1e6)
        (LogicBlock(exc, 2 nodes, 1.26 ms), 1260)
        """
        rf = self.rf
        if slice_offset_m or rf_phase_rad:
            rf = self.offset_rf(rf, slice_offset_m, rf_phase_rad)

        out = LogicBlock(self._tag())
        out.add(0.0, rf) if self.gz is None else out.add(0.0, rf, self.gz)
        if self.rephase and self.gzr is not None:
            out.add(self.pulse_duration, self.gzr)
        return out

    def rephaser(self) -> LogicBlock:
        """
        Return just the slice rephasing lobe, for placing it separately.

        Raises
        ------
        ConfigurationError
            If the pulse has no rephasing lobe.

        Examples
        --------
        >>> import seqcraft as sc
        >>> exc = SincExcitation(sc.System.preset('generic_3t'), flip_deg=15,
        ...                      duration_us=1000, slice_thickness_mm=5)
        >>> exc.rephaser()
        LogicBlock(exc_rephaser, 1 node, 0.54 ms)
        """
        if self.gzr is None:
            msg = format_error(
                f'{type(self).__name__} has no rephasing lobe.',
                {'slice_thickness_mm': self.slice_thickness_mm, 'rephase': self.rephase},
                ['construct it with a slice_thickness_mm to get a slice-select gradient'],
            )
            raise ConfigurationError(msg)
        return LogicBlock(f'{self._tag()}_rephaser').add(0.0, self.gzr)

    @property
    def rephaser_area_per_m(self) -> float:
        """Area of the slice rephasing lobe in 1/m, or zero when there is none."""
        return 0.0 if self.gzr is None else float(self.gzr.area)

    def _tag(self) -> str:
        """The tag given to built blocks.  Short, and derived from the pulse's role."""
        return {
            'excitation': 'exc',
            'refocusing': 'refoc',
            'inversion': 'inv',
            'saturation': 'sat',
            'preparation': 'prep',
        }.get(self.use, 'rf')

    def offset_rf(
        self,
        rf: SimpleNamespace,
        offset_m: float,
        extra_phase_rad: float,
    ) -> SimpleNamespace:
        """
        Return `rf` shifted to an off-centre slice, with the phase compensation applied.

        Parameters
        ----------
        rf
            The designed RF event.  Never mutated.
        offset_m
            Slice position along the slice axis, metres.
        extra_phase_rad
            Carrier phase to add: RF spoiling, phase cycling, a CPMG offset.

        Notes
        -----
        Both terms are needed::

            freq_offset  = G_slice * z                        # selects the slice
            phase_offset = base + extra - 2*pi*f*rf.center     # keeps its phase consistent

        The second term compensates for the RF centre not coinciding with the start of the
        block; without it each slice has a different phase.  The phase *adds* to whatever the
        event already carries, so a design-time offset survives.
        """
        freq = self.slice_select_amplitude * float(offset_m)
        phase = (
            float(rf.phase_offset)
            + float(extra_phase_rad)
            - 2.0 * math.pi * freq * _rf_center(rf)
        )
        return ev.derive(
            rf,
            freq_offset=float(rf.freq_offset) + freq,
            phase_offset=phase % (2.0 * math.pi),
        )

    def _slr(self, *, use: str) -> tuple[Any, Any | None, Any | None]:
        """
        Design a Shinnar--Le Roux pulse with its slice-select and rephasing lobes.

        pypulseq picks the SLR pulse type from `use`: ``'ex'`` for an excitation above 30 degrees
        (``'st'`` below, where small-tip is valid anyway), ``'se'`` for refocusing, ``'inv'`` for
        inversion.  Each solves a different design problem -- an ``'se'`` pulse is designed so that
        its *refocusing* profile is flat, which is not the same as a flat excitation profile, and it
        is why a 180 designed as an excitation pulse has a poor refocusing profile even when its
        excitation profile looks fine.
        """
        duration_s = ceil_to(self.duration_us / 1e6, self.system.rf_raster_s)
        common: dict[str, Any] = {
            'flip_angle': self.flip_rad,
            'duration': duration_s,
            'system': self.opts,
            'time_bw_product': self.time_bw_product,
            'passband_ripple': self.passband_ripple,
            'stopband_ripple': self.stopband_ripple,
            'filter_type': self.filter_type,
            'phase_offset': self.phase_offset_rad,
            'freq_offset': self.freq_offset_hz,
            'use': use,
        }
        if self.slice_thickness_mm is None:
            return pp.make_slr_pulse(**common, return_gz=False), None, None
        rf, gz, gzr = pp.make_slr_pulse(
            **common, slice_thickness=self.slice_thickness_mm / 1e3, return_gz=True
        )
        return rf, _on_axis(gz, self.slice_axis), _on_axis(gzr, self.slice_axis)

    def _sinc(self, *, use: str) -> tuple[Any, Any | None, Any | None]:
        """Design a windowed-sinc pulse with its slice-select and rephasing lobes."""
        duration_s = ceil_to(self.duration_us / 1e6, self.system.rf_raster_s)
        common: dict[str, Any] = {
            'flip_angle': self.flip_rad,
            'duration': duration_s,
            'system': self.opts,
            'time_bw_product': self.time_bw_product,
            'apodization': self.apodization,
            'phase_offset': self.phase_offset_rad,
            'freq_offset': self.freq_offset_hz,
            'use': use,
        }
        if self.slice_thickness_mm is None:
            rf = pp.make_sinc_pulse(**common, return_gz=False)
            return rf, None, None
        rf, gz, gzr = pp.make_sinc_pulse(
            **common,
            slice_thickness=self.slice_thickness_mm / 1e3,
            return_gz=True,
        )
        gz = _on_axis(gz, self.slice_axis)
        gzr = _on_axis(gzr, self.slice_axis)
        return rf, gz, gzr


def _rf_center(rf: SimpleNamespace) -> float:
    """
    Return the RF centre offset within the pulse, preferring the pypulseq 1.5 ``center`` field.

    ``calc_rf_center`` recomputes a centre of mass; for an asymmetric or minimum-phase pulse
    that disagrees with the designed centre, and the designed value is the one timing needs.
    """
    centre = getattr(rf, 'center', None)
    if centre is not None:
        return float(centre)
    return float(pp.calc_rf_center(rf)[0])


def _on_axis(grad: SimpleNamespace | None, axis: str) -> SimpleNamespace | None:
    """Move a gradient event onto `axis`, leaving it alone if it is already there."""
    if grad is None or grad.channel == axis:
        return grad
    return ev.derive(grad, channel=axis)


# --------------------------------------------------------------------------------- concrete
@register()
class SincExcitation(RFPulse):
    """
    Slice-selective windowed-sinc excitation.

    Parameters
    ----------
    See :class:`RFPulse`.  `rephase` defaults to ``True``, because an excitation's slice-select
    gradient leaves through-slice dephasing that must be undone before the readout.

    Properties
    ----------
    See :class:`RFPulse`.

    Build arguments
    ---------------
    See :meth:`RFPulse.build`.

    Notes
    -----
    Slice-select amplitude is ``bandwidth / thickness`` in Hz/m, so a thinner slice or a higher
    time-bandwidth product needs a stronger gradient -- which is what eventually limits how thin
    a slice a given system can excite.

    Examples
    --------
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> exc = SincExcitation(system, flip_deg=15, duration_us=1000, slice_thickness_mm=5)
    >>> round(exc.slice_select_amplitude / (4000 / 5e-3), 6)     # bandwidth / thickness
    1.0
    >>> round(exc.isodelay * 1e6)
    630
    """

    use = 'excitation'

    def _design(self) -> tuple[Any, Any | None, Any | None]:
        """Design the sinc pulse and its slice gradients."""
        return self._sinc(use='excitation')


@register()
class SlabExcitation(SincExcitation):
    """
    Slab-selective excitation for 3D imaging.

    Identical physics to :class:`SincExcitation`; the separate name exists because a slab is
    thick, so the plausibility range for its thickness is different and a 5 mm value passed by
    mistake is worth catching.

    Parameters
    ----------
    See :class:`RFPulse`.  `slice_thickness_mm` is the slab thickness and is required.

    Notes
    -----
    A thick slab needs a weak slice-select gradient, which means a long rephaser is not the
    concern; the concern is the low bandwidth-per-millimetre making the slab profile soft at
    its edges.  Raise `time_bw_product` for a sharper slab, at the cost of peak B1.

    Examples
    --------
    >>> import seqcraft as sc
    >>> slab = SlabExcitation(sc.System.preset('generic_3t'), flip_deg=10,
    ...                       duration_us=1500, slice_thickness_mm=80, time_bw_product=8)
    >>> round(slab.slice_select_amplitude / 1e3, 2)      # weak: 8/1.5ms over 80 mm, kHz/m
    66.67
    """

    def __init__(self, system: System, *, slice_thickness_mm: float, **kwargs: Any) -> None:
        if slice_thickness_mm < 5.0:
            msg = format_error(
                f'SlabExcitation thickness {slice_thickness_mm} mm is thinner than 5 mm.',
                {'slice_thickness_mm': slice_thickness_mm},
                [
                    'a slab covers the whole 3D volume -- 40 to 200 mm is typical',
                    'use SincExcitation for a single thin slice',
                ],
            )
            raise ConfigurationError(msg)
        super().__init__(system, slice_thickness_mm=slice_thickness_mm, **kwargs)


class RefocusingPulse(RFPulse):
    """
    Slice-selective windowed-sinc refocusing pulse, optionally with crushers.

    Parameters
    ----------
    See :class:`RFPulse`, plus:
    crusher_twists
        Crusher strength in dephasing cycles per voxel.  ``None`` for no crushers.
    crusher_voxel_mm
        Voxel size the twists are counted across.  Defaults to the slice thickness.
    crusher_axes
        Which axes carry the crushers.

    Properties
    ----------
    crusher_duration
        Seconds occupied by one crusher lobe.
    crusher_area_per_m
        Crusher area in 1/m.

    Build arguments
    ---------------
    See :meth:`RFPulse.build`.  `rephase` is ignored: a refocusing pulse's slice-select gradient
    is symmetric about the pulse, so there is nothing to rephase.

    Notes
    -----
    **The carrier phase defaults to pi/2**, a quarter turn from an excitation at phase zero -- the
    CPMG condition, and the reason a spin echo is written 90x-180y rather than 90x-180x.  It makes
    the echo amplitude first-order insensitive to flip-angle error on the refocusing pulse, and it
    keeps the stimulated echoes an imperfect pulse creates from interfering with the primary one.

    Measured on a simulated spin echo, the difference is not small: at 20 % low B1 the echo loses
    **17 %** of its amplitude with the refocusing pulse in phase and **2 %** with it at pi/2.  Since
    transmit fields are never uniform across a slice, that is signal lost everywhere off-centre for
    no reason.

    The default assumes the excitation sits at phase zero, which is the usual case -- a spin echo is
    not RF-spoiled.  If you phase-cycle the excitation, pass ``phase_offset_rad`` yourself so the
    **difference** stays a quarter turn; :meth:`RFPulse.offset_rf` composes with it rather than
    overwriting, so ``build(rf_phase_rad=...)`` adds on top.

    The two crushers are **equal**, not opposite.  The refocusing pulse inverts the phase in
    between, so an equal pair leaves the wanted echo untouched while dephasing the FID that an
    imperfect pulse creates.  Making them opposite would crush the echo instead.

    Examples
    --------
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> refoc = SincRefocusing(system, flip_deg=180, duration_us=3000,
    ...                        slice_thickness_mm=5, crusher_twists=4)
    >>> refoc.build()
    LogicBlock(refoc, 4 nodes, 4.63 ms)
    >>> round(refoc.phase_offset_rad, 4)          # a quarter turn: the CPMG condition
    1.5708
    >>> round(refoc.crusher_area_per_m)              # 4 twists across 5 mm
    800
    >>> round(refoc.isodelay * 1e6)                  # includes the leading crusher
    2340
    """

    use = 'refocusing'

    def _design(self) -> tuple[Any, Any | None, Any | None]:
        """Subclasses choose the pulse shape."""
        raise NotImplementedError

    def __init__(
        self,
        system: System,
        *,
        crusher_twists: float | None = None,
        crusher_voxel_mm: float | None = None,
        crusher_axes: tuple[str, ...] = ('z',),
        rephase: bool = False,
        phase_offset_rad: float = math.pi / 2,
        **kwargs: Any,
    ) -> None:
        self.crusher_twists = None if crusher_twists is None else float(crusher_twists)
        self.crusher_voxel_mm = None if crusher_voxel_mm is None else float(crusher_voxel_mm)
        self.crusher_axes = tuple(crusher_axes)
        super().__init__(system, rephase=rephase, phase_offset_rad=phase_offset_rad, **kwargs)
        if self.crusher_twists is not None:
            require_positive(self, 'crusher_twists')
            if self.crusher_voxel_mm is None and self.slice_thickness_mm is None:
                msg = format_error(
                    'crusher_twists needs a voxel size to count twists across.',
                    {'crusher_voxel_mm': None, 'slice_thickness_mm': None},
                    ['pass crusher_voxel_mm=..., or give the pulse a slice_thickness_mm'],
                )
                raise ConfigurationError(msg)
        self.crushers = self._design_crushers()

    def _design_crushers(self) -> tuple[SimpleNamespace, ...]:
        """Build the crusher gradients, or an empty tuple when none were requested."""
        if self.crusher_twists is None:
            return ()
        voxel_mm = self.crusher_voxel_mm or self.slice_thickness_mm
        area = self.crusher_twists * 1e3 / float(voxel_mm)  # twists / voxel, in 1/m
        return tuple(
            pp.make_trapezoid(channel=axis, area=area, system=self.opts)
            for axis in self.crusher_axes
        )

    @property
    def crusher_area_per_m(self) -> float:
        """Crusher area in 1/m, or zero when there are no crushers."""
        return 0.0 if not self.crushers else float(self.crushers[0].area)

    @property
    def crusher_duration(self) -> float:
        """Seconds occupied by one crusher lobe."""
        if not self.crushers:
            return 0.0
        return ceil_to(float(pp.calc_duration(*self.crushers)), self.system.block_raster_s)

    @property
    def isodelay(self) -> float:
        """Seconds from the start of the built block to the pulse centre, crusher included."""
        return self.crusher_duration + super().isodelay

    @property
    def duration(self) -> float:
        """Seconds occupied by the whole built block: crusher, pulse, crusher."""
        return 2.0 * self.crusher_duration + self.pulse_duration

    def build(self, *, slice_offset_m: float = 0.0, rf_phase_rad: float = 0.0) -> LogicBlock:
        """Return crusher, pulse and crusher as one block."""
        rf = self.rf
        if slice_offset_m or rf_phase_rad:
            rf = self.offset_rf(rf, slice_offset_m, rf_phase_rad)

        out = LogicBlock('refoc')
        crusher_dur = self.crusher_duration
        if self.crushers:
            out.add(0.0, *self.crushers)
        out.add(crusher_dur, rf) if self.gz is None else out.add(crusher_dur, rf, self.gz)
        if self.crushers:
            out.add(crusher_dur + self.pulse_duration, *self.crushers)
        return out


@register()
class SincRefocusing(RefocusingPulse):
    """
    Slice-selective windowed-sinc refocusing pulse, optionally with crushers.

    See :class:`RefocusingPulse` for everything except the pulse shape.

    Notes
    -----
    A sinc is a small-tip design, and a 180 is not a small tip -- so its refocusing profile is
    noticeably worse than its excitation profile looks.  Prefer :class:`SLRRefocusing` unless you
    have a reason not to; it costs nothing but a sigpy dependency.

    Examples
    --------
    >>> import seqcraft as sc
    >>> refoc = SincRefocusing(sc.System.preset('generic_3t'), flip_deg=180, duration_us=3000,
    ...                        slice_thickness_mm=5, crusher_twists=4)
    >>> refoc.build()
    LogicBlock(refoc, 4 nodes, 4.63 ms)
    """

    def _design(self) -> tuple[Any, Any | None, Any | None]:
        """Design the sinc pulse; a refocusing pulse needs no rephaser."""
        rf, gz, _gzr = self._sinc(use='refocusing')
        return rf, gz, None


@register()
class SLRExcitation(RFPulse):
    """
    Slice-selective Shinnar--Le Roux excitation.

    Parameters
    ----------
    See :class:`RFPulse`, including `passband_ripple`, `stopband_ripple` and `filter_type`.

    Notes
    -----
    A windowed sinc is derived from the small-tip approximation, where the slice profile is the
    Fourier transform of the envelope.  At 90 degrees that approximation is already poor and at 180
    it is wrong: the real profile has rounded edges and a broadened transition, so the slice is not
    the thickness you asked for and its edges are partially excited.

    SLR solves the actual Bloch problem instead, as a polynomial design with explicit tolerances --
    so you state how flat the passband must be and how much leaks outside, and get a pulse that
    delivers it at the flip angle you asked for.

    Needs ``sigpy``, which pypulseq uses for the design.

    Examples
    --------
    >>> import seqcraft as sc
    >>> exc = SLRExcitation(sc.System.preset('generic_3t'), flip_deg=90, duration_us=2000,
    ...                     slice_thickness_mm=4)
    >>> exc.use
    'excitation'
    >>> exc.build()
    LogicBlock(exc, 3 nodes, 2.76 ms)
    """

    use = 'excitation'

    def _design(self) -> tuple[Any, Any | None, Any | None]:
        """Design the SLR pulse and its slice gradients."""
        return self._slr(use='excitation')


@register()
class SLRRefocusing(RefocusingPulse):
    """
    Slice-selective Shinnar--Le Roux refocusing pulse, optionally with crushers.

    Parameters
    ----------
    See :class:`RefocusingPulse`.

    Notes
    -----
    pypulseq designs this with SLR type ``'se'``, which optimises the **refocusing** profile rather
    than the excitation profile.  Those are different problems, and the distinction is the whole
    point: a 180 designed as though it were an excitation pulse has a poor refocusing profile even
    when its excitation profile looks acceptable, which shows up as signal loss at the slice edges
    and as a slice thinner than requested.

    Carries the same ``pi/2`` default carrier phase as :class:`SincRefocusing` -- see its notes.

    Examples
    --------
    >>> import seqcraft as sc
    >>> refoc = SLRRefocusing(sc.System.preset('generic_3t'), flip_deg=180, duration_us=4000,
    ...                       slice_thickness_mm=4, crusher_twists=4)
    >>> refoc.use
    'refocusing'
    >>> round(refoc.phase_offset_rad, 4)
    1.5708
    """

    def _design(self) -> tuple[Any, Any | None, Any | None]:
        """Design the SLR pulse; a refocusing pulse needs no rephaser."""
        rf, gz, _gzr = self._slr(use='refocusing')
        return rf, gz, None


@register()
class HardExcitation(RFPulse):
    """
    Non-selective rectangular (hard) excitation pulse.

    Parameters
    ----------
    See :class:`RFPulse`.  `slice_thickness_mm` must be ``None`` and `rephase` false: a hard
    pulse excites everything the coil reaches.

    Notes
    -----
    Peak B1 is exactly ``flip / (2*pi*gamma*duration)``, so this is the pulse where the B1 limit
    bites first -- and the one to reach for when TE must be as short as possible, because there
    is no slice-select gradient and no rephaser to fit in.

    Examples
    --------
    >>> import seqcraft as sc
    >>> hard = HardExcitation(sc.System.preset('generic_3t'), flip_deg=90, duration_us=500)
    >>> round(hard.peak_b1_uT, 2)
    11.74
    >>> round(hard.isodelay * 1e6)          # 100 us dead time + half of 500 us
    350
    """

    use = 'excitation'

    def __init__(self, system: System, **kwargs: Any) -> None:
        kwargs.setdefault('rephase', False)
        super().__init__(system, **kwargs)

    def _design(self) -> tuple[Any, Any | None, Any | None]:
        """Design a rectangular pulse."""
        rf = pp.make_block_pulse(
            flip_angle=self.flip_rad,
            duration=ceil_to(self.duration_us / 1e6, self.system.rf_raster_s),
            system=self.opts,
            phase_offset=self.phase_offset_rad,
            freq_offset=self.freq_offset_hz,
            use=self.use,
        )
        return rf, None, None


@register()
class HardRefocusing(HardExcitation):
    """
    Non-selective rectangular refocusing pulse.

    Parameters
    ----------
    See :class:`RFPulse`.  `slice_thickness_mm` must be ``None``: a hard pulse refocuses everything
    the coil reaches, which is what makes it the right choice when a slice-selective 180 would not
    fit in the available TE.  `phase_offset_rad` defaults to ``pi/2`` for the same reason as
    :class:`SincRefocusing` -- see its notes.

    Examples
    --------
    >>> import math
    >>> import seqcraft as sc
    >>> refoc = HardRefocusing(sc.System.preset('generic_3t'), flip_deg=180, duration_us=1000)
    >>> refoc.use
    'refocusing'
    >>> round(refoc.phase_offset_rad, 6) == round(math.pi / 2, 6)
    True
    """

    use = 'refocusing'

    def __init__(
        self, system: System, *, phase_offset_rad: float = math.pi / 2, **kwargs: Any
    ) -> None:
        super().__init__(system, phase_offset_rad=phase_offset_rad, **kwargs)


@register()
class GaussSaturation(RFPulse):
    """
    Gaussian saturation pulse, for fat saturation and saturation bands.

    Parameters
    ----------
    See :class:`RFPulse`.  Pass `freq_offset_hz` for chemical-shift selection; the fat--water
    shift is 3.4 ppm, which is ``-3.4e-6 * gamma * B0`` hertz.

    Notes
    -----
    Gaussian rather than sinc because a saturation pulse wants a smooth, sidelobe-free spectral
    profile: a sinc's sidelobes would partially saturate water.

    Examples
    --------
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> shift = -3.4e-6 * system.gamma * system.b0_T
    >>> fat = GaussSaturation(system, flip_deg=90, duration_us=8000, freq_offset_hz=shift)
    >>> round(fat.freq_offset_hz)
    -434
    """

    use = 'saturation'

    def __init__(self, system: System, **kwargs: Any) -> None:
        kwargs.setdefault('rephase', False)
        super().__init__(system, **kwargs)

    def _design(self) -> tuple[Any, Any | None, Any | None]:
        """Design a Gaussian pulse, slice-selective only if a thickness was given."""
        duration_s = ceil_to(self.duration_us / 1e6, self.system.rf_raster_s)
        common: dict[str, Any] = {
            'flip_angle': self.flip_rad,
            'duration': duration_s,
            'system': self.opts,
            'time_bw_product': self.time_bw_product,
            'apodization': self.apodization,
            'phase_offset': self.phase_offset_rad,
            'freq_offset': self.freq_offset_hz,
            'use': self.use,
        }
        if self.slice_thickness_mm is None:
            return pp.make_gauss_pulse(**common, return_gz=False), None, None
        rf, gz, gzr = pp.make_gauss_pulse(
            **common,
            slice_thickness=self.slice_thickness_mm / 1e3,
            return_gz=True,
        )
        return rf, _on_axis(gz, self.slice_axis), _on_axis(gzr, self.slice_axis)


@register()
class AdiabaticInversion(RFPulse):
    """
    Hyperbolic-secant adiabatic inversion pulse.

    Parameters
    ----------
    See :class:`RFPulse`, plus:
    beta
        Hyperbolic-secant modulation parameter, 1/s.
    mu
        Dimensionless sweep parameter.  ``mu * beta`` is the frequency sweep half-width.

    Notes
    -----
    Adiabatic inversion is insensitive to B1 amplitude above a threshold, which is why it is
    used for inversion recovery at depth where the transmit field is uneven.  The `flip_deg`
    parameter is nominal -- an adiabatic pulse inverts once the adiabatic condition is met,
    regardless of the requested angle -- so it is recorded for provenance and not used to scale
    the amplitude.

    Examples
    --------
    >>> import seqcraft as sc
    >>> inv = AdiabaticInversion(sc.System.preset('generic_3t'), duration_us=10000)
    >>> inv.use
    'inversion'
    >>> round(inv.duration * 1e3, 1)
    10.1
    """

    use = 'inversion'

    def __init__(
        self,
        system: System,
        *,
        duration_us: float = 10_000.0,
        beta: float = 800.0,
        mu: float = 4.9,
        **kwargs: Any,
    ) -> None:
        self.beta = float(beta)
        self.mu = float(mu)
        kwargs.setdefault('flip_deg', 180.0)
        kwargs.setdefault('rephase', False)
        super().__init__(system, duration_us=duration_us, **kwargs)

    def _check_b1(self) -> None:
        """An adiabatic pulse does not follow the hard-pulse B1 floor, so skip that check."""

    def _design(self) -> tuple[Any, Any | None, Any | None]:
        """Design a hyperbolic-secant pulse."""
        rf = pp.make_adiabatic_pulse(
            pulse_type='hypsec',
            duration=ceil_to(self.duration_us / 1e6, self.system.rf_raster_s),
            system=self.opts,
            beta=self.beta,
            mu=self.mu,
            phase_offset=self.phase_offset_rad,
            freq_offset=self.freq_offset_hz,
            use=self.use,
        )
        return rf, None, None
