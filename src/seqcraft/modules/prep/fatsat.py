"""
Magnetisation preparation: fat saturation and inversion recovery.

Each of these is a nested module -- an RF pulse plus a spoiler, or a pulse plus a delay -- which
is the case that shows why nesting needs no special mechanism: a preparation's ``build()``
returns a block whose children are other modules' blocks.

Examples
--------
>>> import seqcraft as sc
>>> system = sc.System.preset('generic_3t')
>>> fat = FatSat(system, voxel_mm=5)
>>> fat.build()
LogicBlock(fatsat, 2 nodes, ... ms)
>>> round(fat.freq_offset_hz)              # -3.4 ppm at 3 T
-434
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ...core.errors import ConfigurationError, format_error
from ...core.logic import LogicBlock
from ...core.module import Module
from ...core.raster import ceil_to
from ...core.registry import register
from ...core.validate import require_positive
from ..encoding.cartesian import Spoiler
from ..rf.pulses import AdiabaticInversion, GaussSaturation

if TYPE_CHECKING:
    from ...core.system import System

__all__ = ['FatSat', 'InversionRecovery']

#: Fat--water chemical shift, parts per million.  Negative: fat resonates below water.
FAT_SHIFT_PPM = -3.4


@register()
class FatSat(Module):
    """
    Chemical-shift-selective fat saturation: a spectrally selective pulse, then a spoiler.

    Parameters
    ----------
    system
        The scanner.
    voxel_mm
        Voxel size the spoiler dephases across.
    flip_deg
        Saturation flip angle.  90 degrees nulls the fat magnetisation; slightly more
        over-saturates, which is more robust to B1 inhomogeneity at the cost of some water
        saturation through the pulse's spectral tails.
    duration_us
        Pulse duration.  Must be long enough that its bandwidth is comfortably narrower than the
        fat--water separation, or water is saturated too -- see Notes.
    time_bw_product
        Time-bandwidth product.  Default 1.6, which gives a 200 Hz bandwidth over 8 ms: narrow
        enough to miss water at 1.5 T as well as 3 T.  The 4.0 that suits a slice-selective pulse
        is far too wide here.
    shift_ppm
        Chemical shift to target.  Defaults to the fat--water shift.
    twists
        Spoiler strength, dephasing cycles per voxel.
    spoil_axes
        Which axes carry the spoiler.
    max_bandwidth_fraction
        The largest bandwidth accepted, as a fraction of the fat--water separation.  Default 0.7:
        a Gaussian's spectral profile has tails beyond its nominal width, so matching the
        separation exactly still saturates some water.

    Properties
    ----------
    freq_offset_hz
        The carrier offset, ``shift_ppm * 1e-6 * gamma * B0``.
    bandwidth_hz
        The pulse's spectral width.
    duration

    Notes
    -----
    The pulse must be **spectrally narrow enough to miss water**.  At 3 T the fat--water
    separation is only about 434 Hz -- and at 1.5 T it is 217 Hz -- so a pulse whose bandwidth
    approaches that saturates water too, and the image loses signal everywhere rather than just in
    fat.  The constructor checks this and refuses, because the failure looks like poor SNR rather
    than like a bug.  That is also why the default `time_bw_product` is 1.6 and not the 4.0 that
    suits a slice-selective pulse.

    The spoiler is not optional and has no default: forgetting to spoil after a saturation pulse
    leaves the saturated magnetisation coherent, so it refocuses into the image as a ghost -- a
    classic silent artifact.

    Examples
    --------
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> fat = FatSat(system, voxel_mm=5, duration_us=8000)
    >>> round(fat.bandwidth_hz)
    200
    >>> round(fat.freq_offset_hz)
    -434
    >>> fat.bandwidth_hz < 0.7 * abs(fat.freq_offset_hz)      # narrow enough to miss water
    True
    """

    def __init__(
        self,
        system: System,
        *,
        voxel_mm: float,
        flip_deg: float = 90.0,
        duration_us: float = 8000.0,
        time_bw_product: float = 1.6,
        shift_ppm: float = FAT_SHIFT_PPM,
        twists: float = 4.0,
        spoil_axes: tuple[str, ...] = ('x', 'y', 'z'),
        max_bandwidth_fraction: float = 0.7,
        regime: str = 'default',
    ) -> None:
        super().__init__(system, regime=regime)
        self.voxel_mm = float(voxel_mm)
        self.flip_deg = float(flip_deg)
        self.duration_us = float(duration_us)
        self.time_bw_product = float(time_bw_product)
        self.shift_ppm = float(shift_ppm)
        self.twists = float(twists)
        self.max_bandwidth_fraction = float(max_bandwidth_fraction)
        require_positive(self, 'voxel_mm', 'duration_us', 'twists', 'time_bw_product')

        offset = self.freq_offset_hz
        bandwidth = self.bandwidth_hz
        allowed = self.max_bandwidth_fraction * abs(offset)
        if bandwidth > allowed:
            needed_us = self.duration_us * bandwidth / allowed
            msg = format_error(
                f'a {self.duration_us:g} us saturation pulse with time_bw_product='
                f'{self.time_bw_product:g} has a {bandwidth:.0f} Hz bandwidth, above the '
                f'{allowed:.0f} Hz allowed against the {abs(offset):.0f} Hz fat-water separation '
                f'at {self.system.b0_T:.2f} T -- it would saturate water as well as fat.',
                {
                    'bandwidth_hz': round(bandwidth),
                    'separation_hz': round(abs(offset)),
                    'allowed_hz': round(allowed),
                    'b0_T': self.system.b0_T,
                },
                [
                    f'lengthen duration_us to at least {needed_us:.0f} us',
                    f'or reduce time_bw_product below '
                    f'{self.time_bw_product * allowed / bandwidth:.2f}',
                ],
            )
            raise ConfigurationError(msg)

        self.pulse = GaussSaturation(
            system,
            flip_deg=self.flip_deg,
            duration_us=self.duration_us,
            time_bw_product=self.time_bw_product,
            freq_offset_hz=offset,
            regime=regime,
        )
        self.spoiler = Spoiler(
            system, twists=self.twists, voxel_mm=self.voxel_mm, axes=spoil_axes, regime=regime
        )

    @property
    def freq_offset_hz(self) -> float:
        """Carrier frequency offset, ``shift_ppm * 1e-6 * gamma * B0``, in hertz."""
        return self.shift_ppm * 1e-6 * self.system.gamma * self.system.b0_T

    @property
    def bandwidth_hz(self) -> float:
        """Spectral width of the saturation pulse, ``time_bw_product / duration``, in hertz."""
        return self.time_bw_product / (self.duration_us / 1e6)

    @property
    def duration(self) -> float:
        """Seconds occupied by the pulse and its spoiler."""
        return self.pulse.duration + self.spoiler.duration

    def build(self) -> LogicBlock:
        """
        Return the saturation pulse followed by the spoiler.

        Examples
        --------
        >>> import seqcraft as sc
        >>> fat = FatSat(sc.System.preset('generic_3t'), voxel_mm=5)
        >>> len(fat.build())               # two nested blocks
        2
        """
        pulse = self.pulse.build()
        return LogicBlock('fatsat').add(0.0, pulse).add(pulse.duration, self.spoiler.build())


@register()
class InversionRecovery(Module):
    """
    Adiabatic inversion followed by a spoiler, for a T1-nulling preparation.

    Parameters
    ----------
    system
        The scanner.
    voxel_mm
        Voxel size the spoiler dephases across.
    duration_us
        Inversion pulse duration.
    twists
        Spoiler strength.
    slice_thickness_mm
        ``None`` for a non-selective inversion, which is what inversion recovery normally wants.

    Properties
    ----------
    duration
        Seconds occupied by the pulse and spoiler.  **The inversion time TI is measured from the
        pulse centre**, so place the excitation at ``t_inversion + inv.isodelay + ti``.

    Notes
    -----
    TI is conventionally measured from the *centre* of the inversion pulse to the centre of the
    excitation, which is why :attr:`isodelay` is exposed separately from :attr:`duration`.  A 10 ms
    adiabatic pulse means the two differ by 5 ms, which at short TI is a large fraction of the
    interval.

    To null a tissue of relaxation time ``T1``, use ``TI = T1 * ln(2)`` --
    about 0.693 T1 for a perfect inversion from full recovery.

    Examples
    --------
    >>> import math
    >>> import seqcraft as sc
    >>> ir = InversionRecovery(sc.System.preset('generic_3t'), voxel_mm=5)
    >>> round(1000 * math.log(2))          # TI to null a 1000 ms T1, in ms
    693
    >>> ir.isodelay < ir.duration
    True
    """

    def __init__(
        self,
        system: System,
        *,
        voxel_mm: float,
        duration_us: float = 10_000.0,
        twists: float = 4.0,
        slice_thickness_mm: float | None = None,
        regime: str = 'default',
    ) -> None:
        super().__init__(system, regime=regime)
        self.voxel_mm = float(voxel_mm)
        self.duration_us = float(duration_us)
        self.twists = float(twists)
        require_positive(self, 'voxel_mm', 'duration_us', 'twists')
        self.pulse = AdiabaticInversion(
            system,
            duration_us=self.duration_us,
            slice_thickness_mm=slice_thickness_mm,
            regime=regime,
        )
        self.spoiler = Spoiler(
            system, twists=self.twists, voxel_mm=self.voxel_mm, axes=('x', 'y', 'z'), regime=regime
        )

    @property
    def isodelay(self) -> float:
        """Seconds from the start of the built block to the inversion pulse's centre."""
        return self.pulse.isodelay

    @property
    def duration(self) -> float:
        """Seconds occupied by the pulse and its spoiler."""
        return ceil_to(
            self.pulse.duration + self.spoiler.duration, self.system.block_raster_s
        )

    def build(self) -> LogicBlock:
        """Return the inversion pulse followed by the spoiler."""
        pulse = self.pulse.build()
        return LogicBlock('inversion').add(0.0, pulse).add(pulse.duration, self.spoiler.build())
