"""
Cartesian frequency-encoded readout: one k-space line per ADC.

The module owns its own prephaser, because the prephaser's area is exactly minus half the
readout's own and nothing outside the readout can know that.  It can still be placed separately
via :meth:`CartesianLine.prephaser_block`, which is what lets it overlap a slice rewinder and a
phase-encode blip -- three gradients on three axes that the compiler then puts in one block.

Examples
--------
>>> import seqcraft as sc
>>> system = sc.System.preset('generic_3t')
>>> ro = CartesianLine(system, fov_ro_mm=250, matrix_ro=64, readout_duration_us=3200)
>>> ro.build()
LogicBlock(readout, 3 nodes, 3.54 ms)
>>> round(ro.bandwidth_per_pixel_hz, 1)
312.5
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pypulseq as pp

from ...core import events as ev
from ...core.errors import ConfigurationError, format_error
from ...core.logic import LogicBlock
from ...core.module import Module
from ...core.units import convert
from ...core.validate import Range, require_in, require_in_range, require_int_in, require_positive

if TYPE_CHECKING:
    from types import SimpleNamespace

    from ...core.system import System

__all__ = ['CartesianLine', 'NoiseAcquisition']

_FOV_RANGE = Range(0.5, 2000.0, 'mm', ((1e3, 'm'), (10.0, 'cm')))


class CartesianLine(Module):
    """
    One frequency-encoded k-space line, with its prephaser.

    Parameters
    ----------
    system
        The scanner.
    fov_ro_mm
        Field of view along the readout direction, millimetres.
    matrix_ro
        Number of samples across the FOV, before oversampling.
    readout_duration_us
        Duration of the flat top of the readout gradient.  Longer means finer bandwidth per
        pixel, hence more SNR, but more off-resonance distortion and a longer minimum TE.
    oversampling
        Readout oversampling factor.  Two is the Siemens default and suppresses aliasing along
        the readout direction.
    axis
        Which logical axis carries the readout.
    ramp_sampling
        Sample during the gradient ramps as well as the flat top, which shortens the readout for
        a given resolution.  The trajectory is then non-uniform, so the reconstruction must grid
        it; ``calculate_kspacePP`` gives the correct sample positions either way.
    partial_echo
        Fraction of the **leading** half of k-space to acquire, from 0.5 to 1.0.  ``1.0`` samples
        symmetrically about the echo; ``0.75`` starts at ``-0.75 k_max`` and still ends at
        ``+k_max``, so a quarter of the leading half is dropped and TE shortens by that much.

        The resolution does not change -- ``dk`` per sample is fixed by ``fov_ro_mm`` -- so what is
        given up is the conjugate-symmetric part of k-space, which the reconstruction must fill in
        (zero-filling blurs, a homodyne or POCS reconstruction does better).  The sample count is
        snapped to the ADC divisor, and the attribute is then updated to the fraction actually
        achieved rather than the one requested.
    prephase
        Include the prephaser in :meth:`build`.

    Properties
    ----------
    time_to_echo
        Seconds from the start of the built block to k=0.  **What TE is measured from.**
    duration_after_echo
        Seconds from k=0 to the end of the block.  The pair, rather than one "echo time", is
        what lets symmetric and asymmetric readouts be placed by the same arithmetic.
    bandwidth_per_pixel_hz
        ``1 / (dwell * n_samples)``.
    dk_per_m, k_max_per_m, resolution_mm, duration, prephase_duration

    Build arguments
    ---------------
    polarity : {1, -1}, default 1
        Flip the readout direction, for an alternating or bipolar readout.

        `prephase` is deliberately **not** a build argument: dropping the prephaser changes the
        block's duration, which would silently invalidate :attr:`time_to_echo` and move the echo.
        Construct with ``prephase=False`` and place :meth:`prephaser_block` yourself instead --
        then the property and the block always agree.

    Notes
    -----
    With an even sample count there is **no sample exactly at k=0**: the two central samples
    straddle it by half a dwell.  So a test asserting "a sample lands exactly at TE" would be wrong
    rather than the code; what has to be true is that ``time_to_echo`` points at k=0, which is
    checked against pypulseq's own trajectory.

    **The prephaser has to cancel the ramp's area as well as half the flat top.**  Minus ``k_max``
    is the tempting value and it is wrong by ``amplitude * rise_time / 2`` -- one ``dk`` on a wide
    readout and nearly four on a short low-bandwidth one.  The whole line then sits off-centre in
    k-space, and since a k-space offset is a linear phase ramp across the image, a magnitude image
    looks completely normal while anything built on the phase is wrong.  Both the prephaser and
    ``time_to_echo`` come from :meth:`_moment_to`, so they cannot disagree.

    Examples
    --------
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> ro = CartesianLine(system, fov_ro_mm=250, matrix_ro=64, readout_duration_us=3200)
    >>> round(ro.k_max_per_m, 1)                 # matrix / (2 * FOV)
    128.0
    >>> round(ro.resolution_mm, 4)
    3.9062
    >>> abs(ro.duration - ro.time_to_echo - ro.duration_after_echo) < 1e-12
    True

    The prephaser is larger than ``k_max`` by exactly the ramp's area:

    >>> ramp = float(ro.gx.amplitude) * float(ro.gx.rise_time) / 2.0
    >>> abs(float(ro.pre.area) + ro.k_max_per_m + ramp) < 1e-6
    True

    A partial echo drops leading samples and brings the echo forward:

    >>> short = CartesianLine(system, fov_ro_mm=250, matrix_ro=64, readout_duration_us=3200,
    ...                       partial_echo=0.75)
    >>> short.n_samples, ro.n_samples
    (56, 64)
    >>> short.time_to_echo < ro.time_to_echo
    True
    >>> round(short.resolution_mm, 4)            # unchanged: dk per sample is fixed
    3.9062
    """

    def __init__(
        self,
        system: System,
        *,
        fov_ro_mm: float,
        matrix_ro: int,
        readout_duration_us: float,
        oversampling: float = 1.0,
        axis: str = 'x',
        ramp_sampling: bool = False,
        partial_echo: float = 1.0,
        prephase: bool = True,
        regime: str = 'default',
    ) -> None:
        super().__init__(system, regime=regime)
        self.fov_ro_mm = float(fov_ro_mm)
        self.matrix_ro = int(matrix_ro)
        self.readout_duration_us = float(readout_duration_us)
        self.oversampling = float(oversampling)
        self.axis = str(axis)
        self.ramp_sampling = bool(ramp_sampling)
        self.partial_echo = float(partial_echo)
        self.prephase = bool(prephase)

        require_positive(self, 'fov_ro_mm', 'readout_duration_us', 'oversampling')
        require_in_range(self, 'fov_ro_mm', _FOV_RANGE.lo, _FOV_RANGE.hi, unit='mm')
        require_in_range(self, 'partial_echo', 0.5, 1.0, unit='')
        require_int_in(self, 'matrix_ro', lo=2, hi=8192)
        require_in(self, 'axis', ('x', 'y', 'z'))

        divisor = self.system.adc_samples_divisor
        n_full = self.matrix_ro * self.oversampling
        if abs(int(round(n_full / divisor)) * divisor - n_full) > 1e-9:
            msg = format_error(
                f'{self.matrix_ro} x {self.oversampling:g} = {n_full:g} ADC samples is not a '
                f'multiple of the required divisor {divisor}.',
                {'nearest multiple': int(round(n_full / divisor)) * divisor},
                [
                    f'use a matrix_ro that is a multiple of {divisor}',
                    f'or an oversampling that makes the product one',
                ],
            )
            raise ConfigurationError(msg)

        # A partial echo drops samples from the *leading* half only, so k runs from
        # -partial_echo * k_max up to +k_max.  Snap to the divisor and report back what that
        # achieved, rather than pretending the request was met exactly.
        requested = n_full * (1.0 + self.partial_echo) / 2.0
        self._n_samples = max(divisor, int(round(requested / divisor)) * divisor)
        self.partial_echo = 2.0 * self._n_samples / n_full - 1.0

        flat_time = self.system.grad_raster.ceil(
            convert(self.readout_duration_us, 'us', 's')
        )
        self.gx = pp.make_trapezoid(
            channel=self.axis,
            # Delta-k per sample is fixed by the resolution, so a shorter readout covers less of
            # k-space rather than covering the same extent more coarsely.
            flat_area=2.0 * self.k_max_per_m * self._n_samples / n_full,
            flat_time=flat_time,
            system=self.opts,
        )
        window = flat_time
        if self.ramp_sampling:
            window += float(self.gx.rise_time) + float(self.gx.fall_time)
        dwell = self.system.adc_raster.nearest(window / self._n_samples)
        self.adc = pp.make_adc(
            num_samples=self._n_samples,
            dwell=dwell,
            delay=0.0 if self.ramp_sampling else float(self.gx.rise_time),
            system=self.opts,
        )
        self.pre = pp.make_trapezoid(
            channel=self.axis,
            area=-(self.partial_echo * self.k_max_per_m + self._moment_to(float(self.adc.delay))),
            system=self.opts,
        )

    def _moment_to(self, when: float) -> float:
        """
        Area under the readout gradient from its own start to `when`, in 1/m.

        One expression for every question about where k is, which is what keeps the prephaser and
        :attr:`time_to_echo` from disagreeing.  The ramp is the part that gets forgotten: a
        prephaser of exactly minus ``k_max`` leaves the whole line offset by the ramp's own area,
        which is about one ``dk`` here -- and a k-space offset is a linear phase ramp across the
        image, so it corrupts a phase measurement while leaving a magnitude image looking fine.
        """
        amplitude = float(self.gx.amplitude)
        rise, flat = float(self.gx.rise_time), float(self.gx.flat_time)
        fall = float(self.gx.fall_time)
        if when <= 0.0:
            return 0.0
        if when < rise:
            return amplitude * when * when / (2.0 * rise)
        if when <= rise + flat:
            return amplitude * (rise / 2.0 + (when - rise))
        past = min(when - rise - flat, fall)
        return amplitude * (rise / 2.0 + flat + past * (1.0 - past / (2.0 * fall)))

    # --------------------------------------------------------------------------- properties
    @property
    def dk_per_m(self) -> float:
        """k-space step between adjacent lines, ``1 / FOV``, in 1/m."""
        return 1e3 / self.fov_ro_mm

    @property
    def k_max_per_m(self) -> float:
        """Outer k-space radius along the readout, ``matrix / (2 * FOV)``, in 1/m."""
        return self.matrix_ro / (2.0 * convert(self.fov_ro_mm, 'mm', 'm'))

    @property
    def resolution_mm(self) -> float:
        """Readout resolution, ``FOV / matrix``, in millimetres."""
        return self.fov_ro_mm / self.matrix_ro

    @property
    def n_samples(self) -> int:
        """ADC samples per line, oversampling included."""
        return self._n_samples

    @property
    def dwell_s(self) -> float:
        """ADC dwell time, seconds."""
        return float(self.adc.dwell)

    @property
    def adc_duration(self) -> float:
        """Length of the ADC window, seconds.  An ADC event has no ``duration`` attribute."""
        return float(self.adc.num_samples) * float(self.adc.dwell)

    @property
    def bandwidth_per_pixel_hz(self) -> float:
        """Receive bandwidth per pixel, ``1 / (dwell * n_samples)``."""
        return 1.0 / (self.dwell_s * self._n_samples)

    @property
    def prephase_duration(self) -> float:
        """Seconds occupied by the prephaser."""
        return self.system.block_raster.ceil(float(pp.calc_duration(self.pre)))

    @property
    def readout_block_duration(self) -> float:
        """Seconds occupied by the readout gradient and its ADC."""
        return self.system.block_raster.ceil(float(pp.calc_duration(self.gx, self.adc)))

    @property
    def duration(self) -> float:
        """Seconds occupied by the whole built block."""
        lead = self.prephase_duration if self.prephase else 0.0
        return lead + self.readout_block_duration

    @property
    def time_to_echo(self) -> float:
        """
        Seconds from the start of the built block to **k = 0**.

        **What TE is measured from**: placing this readout at ``te - time_to_echo`` *is* the
        definition of TE for the sequence.  With an even sample count no sample sits exactly at
        k=0; the two central ones straddle it.

        Found from the moment rather than assumed to be the middle of the ADC window.  For a full
        echo the two coincide, and for a partial echo they do not: k=0 arrives earlier, which is the
        entire reason to shorten the leading half.  Taking the window centre would put TE later than
        requested by half the samples that were dropped, silently.
        """
        lead = self.prephase_duration if self.prephase else 0.0
        # k(t) = pre.area + moment(t), and the moment is monotonic, so invert it on the flat top.
        amplitude = float(self.gx.amplitude)
        rise = float(self.gx.rise_time)
        wanted = -float(self.pre.area) - amplitude * rise / 2.0
        return lead + rise + wanted / amplitude

    @property
    def duration_after_echo(self) -> float:
        """Seconds from k=0 to the end of the built block."""
        return self.duration - self.time_to_echo

    def definitions(self) -> dict[str, float]:
        """Return the ``.seq`` definitions this readout is responsible for."""
        return {
            'BandwidthPerPixelHz': self.bandwidth_per_pixel_hz,
            'ReadoutOversamplingFactor': self.oversampling,
        }

    # -------------------------------------------------------------------------------- build
    def build(self, *, polarity: int = 1, rf_phase_rad: float = 0.0) -> LogicBlock:
        """
        Return the prephaser (if this readout owns one), the readout gradient and the ADC.

        Parameters
        ----------
        polarity
            ``-1`` reverses the readout direction, for an alternating or bipolar readout.
        rf_phase_rad
            Receiver phase, which must be **the same value passed to the excitation**.

            Whatever phase the transmitter used is carried by the magnetisation, so the receiver has
            to demodulate with it or the phase of every shot differs. For a magnitude image that is
            invisible; for anything that reads the phase -- a field map, a flow measurement, any
            complex averaging -- it is fatal, and it is invisible in simulation if the simulator
            aligns the two for you.

        Examples
        --------
        >>> import seqcraft as sc
        >>> system = sc.System.preset('generic_3t')
        >>> ro = CartesianLine(system, fov_ro_mm=250, matrix_ro=64, readout_duration_us=3200)
        >>> ro.build()
        LogicBlock(readout, 3 nodes, 3.54 ms)
        >>> bare = CartesianLine(system, fov_ro_mm=250, matrix_ro=64,
        ...                      readout_duration_us=3200, prephase=False)
        >>> bare.build()
        LogicBlock(readout, 2 nodes, 3.24 ms)
        >>> round(bare.time_to_echo * 1e6)            # the property follows the constructor
        1620
        >>> spoiled = bare.build(rf_phase_rad=1.25)
        >>> round(float(spoiled.nodes[-1].item.phase_offset), 3)
        1.25
        """
        if polarity not in (1, -1):
            msg = format_error(
                f'polarity must be 1 or -1, got {polarity!r}.',
                {},
                ['use -1 for a reversed readout, as in alternating EPI lines'],
            )
            raise ConfigurationError(msg)

        out = LogicBlock('readout')
        lead = 0.0
        if self.prephase:
            out.add(0.0, self.pre if polarity == 1 else _flip(self.pre))
            lead = self.prephase_duration
        adc = (self.adc if rf_phase_rad == 0.0
               else ev.derive(self.adc, phase_offset=float(rf_phase_rad)))
        out.add(lead, self.gx if polarity == 1 else _flip(self.gx), adc)
        return out

    def prephaser_block(self, *, polarity: int = 1) -> LogicBlock:
        """
        Return just the prephaser, for placing it separately.

        Examples
        --------
        >>> import seqcraft as sc
        >>> ro = CartesianLine(sc.System.preset('generic_3t'), fov_ro_mm=250, matrix_ro=64,
        ...                    readout_duration_us=3200)
        >>> ro.prephaser_block()
        LogicBlock(readout_pre, 1 node, 0.30 ms)
        """
        return LogicBlock('readout_pre').add(
            0.0, self.pre if polarity == 1 else _flip(self.pre)
        )


def _flip(grad: SimpleNamespace) -> SimpleNamespace:
    """Return a trapezoid with its sign reversed."""
    return ev.derive(
        grad,
        amplitude=-float(grad.amplitude),
        area=-float(grad.area),
        flat_area=-float(getattr(grad, 'flat_area', 0.0)),
    )


class NoiseAcquisition(Module):
    """
    An ADC with no RF and no gradients, for a noise-covariance measurement.

    Parameters
    ----------
    system
        The scanner.
    n_samples
        Number of samples.  Must be a multiple of the ADC sample divisor.
    dwell_ns
        Dwell time, nanoseconds.

    Notes
    -----
    A noise scan must be labelled ``NOISE`` so the reconstruction excludes it from k-space -- and
    so seqcraft's own duplicate-address check skips it.  Add the label alongside this block::

        seq.add(t, noise.build(), pp.make_label('NOISE', 'SET', True))

    Examples
    --------
    >>> import seqcraft as sc
    >>> noise = NoiseAcquisition(sc.System.preset('generic_3t'), n_samples=256)
    >>> noise.build()
    LogicBlock(noise, 1 node, 0.53 ms)
    """

    def __init__(
        self,
        system: System,
        *,
        n_samples: int = 256,
        dwell_ns: float = 2000.0,
        regime: str = 'default',
    ) -> None:
        super().__init__(system, regime=regime)
        self.n_samples = int(n_samples)
        self.dwell_ns = float(dwell_ns)
        require_int_in(self, 'n_samples', lo=4, hi=1_000_000)
        self.adc = pp.make_adc(
            num_samples=self.n_samples,
            dwell=self.system.adc_raster.nearest(convert(self.dwell_ns, 'ns', 's')),
            system=self.opts,
        )

    @property
    def duration(self) -> float:
        """Seconds occupied."""
        return self.system.block_raster.ceil(float(pp.calc_duration(self.adc)))

    def build(self) -> LogicBlock:
        """Return the bare ADC."""
        return LogicBlock('noise').add(0.0, self.adc)
