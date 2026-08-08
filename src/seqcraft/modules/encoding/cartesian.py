"""
Cartesian encoding gradients: phase encode, partition encode, spoiler, crusher.

Each of these is a single trapezoid whose *area* carries the meaning, so each is parameterised
in the physical quantity it encodes rather than in an amplitude:

* a phase encode moves k-space by ``line * (1 / FOV)``, in 1/m;
* a spoiler dephases by a chosen number of cycles across a voxel, so its area is
  ``twists / voxel``.

Amplitude then falls out of the area and the available time, which is why the same module works
unchanged on a 40 mT/m system and a 200 mT/m one.

Examples
--------
>>> import seqcraft as sc
>>> system = sc.System.preset('generic_3t')
>>> pe = PhaseEncode(system, fov_pe_mm=250, matrix_pe=64)
>>> round(pe.dk_per_m, 3)                       # 1 / 250 mm
4.0
>>> round(pe.area_for(line=16), 3)              # 16 lines from the centre
64.0
>>> pe.build(line=16)
LogicBlock(pe, 1 node, 0.30 ms)
>>> pe.build(line=0)                            # k=0 still occupies the slot
LogicBlock(pe, 1 node, 0.30 ms)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pypulseq as pp

from ...core import events as ev
from ...core.errors import ConfigurationError, format_error
from ...core.logic import LogicBlock
from ...core.units import convert
from ...core.validate import Range, require_in, require_int_in, require_positive
from ..base import Module

if TYPE_CHECKING:
    from types import SimpleNamespace

    from ...core.system import System

__all__ = ['Crusher', 'PartitionEncode', 'PhaseEncode', 'Prephaser', 'Spoiler']

_FOV_RANGE = Range(0.5, 2000.0, 'mm', ((1e3, 'm'), (10.0, 'cm')))


class _AreaTrapezoid(Module):
    """
    Shared behaviour for a module that is one trapezoid on one or more axes.

    Designs at the maximum area it will ever need, so every variant has the **same duration**.
    That matters more than it looks: if a phase-encode blip got shorter for lines near the
    centre of k-space, the echo would move, TE would vary line by line, and the image would
    carry a phase ramp nothing in the reconstruction accounts for.

    Parameters
    ----------
    max_area_per_m
        The largest area this module will be asked for, 1/m.  Sets the duration.
    axes
        Which logical axes carry the gradient.
    min_duration_us
        Optional floor on the duration, to derate the slew for a quieter or PNS-safer waveform.
    """

    def __init__(
        self,
        system: System,
        *,
        max_area_per_m: float,
        axes: tuple[str, ...],
        min_duration_us: float | None = None,
        regime: str = 'default',
    ) -> None:
        super().__init__(system, regime=regime)
        self.max_area_per_m = float(max_area_per_m)
        self.axes = tuple(axes)
        self.min_duration_us = None if min_duration_us is None else float(min_duration_us)
        for axis in self.axes:
            if axis not in ('x', 'y', 'z'):
                msg = format_error(
                    f'{type(self).__name__}: axis must be one of x, y, z.',
                    {'got': axis, 'axes': ', '.join(self.axes)},
                )
                raise ConfigurationError(msg)

        area = abs(self.max_area_per_m) or 1.0
        kwargs = {'system': self.opts, 'area': area}
        if self.min_duration_us is not None:
            kwargs['duration'] = self.system.grad_raster.ceil(
                convert(self.min_duration_us, 'us', 's')
            )
        self._reference = {
            axis: pp.make_trapezoid(channel=axis, **kwargs) for axis in self.axes
        }

    @property
    def duration(self) -> float:
        """Seconds occupied, the same for every area this module can produce."""
        return self.system.block_raster.ceil(
            float(pp.calc_duration(*self._reference.values()))
        )

    def _scaled(self, area_per_m: float) -> tuple[SimpleNamespace, ...]:
        """Return the reference trapezoids scaled to `area_per_m`, keeping their duration."""
        reference = next(iter(self._reference.values()))
        factor = area_per_m / float(reference.area)
        return tuple(
            ev.derive(g, amplitude=float(g.amplitude) * factor, area=float(g.area) * factor)
            for g in self._reference.values()
        )

    def _block(self, tag: str, area_per_m: float) -> LogicBlock:
        """
        Return a block holding the scaled trapezoids, or a delay when the area is zero.

        A zero-area encode still has to occupy its slot, or the k=0 line would have a different
        TE from every other line.  A delay event does that with no gradient at all, which is
        both honest and cheaper than a zero-amplitude trapezoid.
        """
        out = LogicBlock(tag)
        if area_per_m == 0.0:
            return out.add(0.0, pp.make_delay(self.duration))
        return out.add(0.0, *self._scaled(area_per_m))


class PhaseEncode(_AreaTrapezoid):
    """
    In-plane phase-encode blip.

    Parameters
    ----------
    system
        The scanner.
    fov_pe_mm
        Field of view along the phase-encode direction, millimetres.
    matrix_pe
        Number of phase-encode lines.  Sets the largest blip, hence the duration.
    axis
        Which logical axis carries the blip.
    min_duration_us
        Optional duration floor, to derate the slew.

    Properties
    ----------
    dk_per_m
        k-space step per line, ``1 / FOV`` in 1/m.
    duration
        Seconds occupied, identical for every line.

    Build arguments
    ---------------
    line : int
        Signed line index, measured from the centre of k-space.  ``0`` is the k=0 line.
    scale : float, default 1.0
        Extra multiplier, for rewinding (``-1``) or fractional steps.

    Notes
    -----
    ``line`` is signed and centred, not a raw counter from zero.  That is what makes partial
    Fourier and parallel imaging expressible without arithmetic at the call site: the caller
    passes the line it means, and
    :attr:`~seqcraft.core.geometry.Geometry.pe_lines` produces the right set.

    Examples
    --------
    >>> import seqcraft as sc
    >>> pe = PhaseEncode(sc.System.preset('generic_3t'), fov_pe_mm=250, matrix_pe=64)
    >>> round(pe.area_for(line=-32), 1)             # the most negative line
    -128.0
    >>> round(pe.duration * 1e6)
    300
    >>> pe.build(line=-32, scale=-1.0)              # the rewinder for that line
    LogicBlock(pe, 1 node, 0.30 ms)
    """

    def __init__(
        self,
        system: System,
        *,
        fov_pe_mm: float,
        matrix_pe: int,
        axis: str = 'y',
        min_duration_us: float | None = None,
        regime: str = 'default',
    ) -> None:
        self.fov_pe_mm = float(fov_pe_mm)
        self.matrix_pe = int(matrix_pe)
        require_positive(self, 'fov_pe_mm')
        require_int_in(self, 'matrix_pe', lo=1, hi=8192)
        dk = 1e3 / self.fov_pe_mm
        super().__init__(
            system,
            max_area_per_m=dk * (self.matrix_pe / 2.0),
            axes=(axis,),
            min_duration_us=min_duration_us,
            regime=regime,
        )
        self.axis = axis
        require_in(self, 'axis', ('x', 'y', 'z'))

    @property
    def dk_per_m(self) -> float:
        """k-space step per line, 1/m."""
        return 1e3 / self.fov_pe_mm

    def area_for(self, line: float) -> float:
        """Return the gradient area for a signed, centred line index, in 1/m."""
        return self.dk_per_m * float(line)

    def build(self, *, line: float = 0.0, scale: float = 1.0) -> LogicBlock:
        """Return the blip for `line`."""
        return self._block('pe', self.area_for(line) * float(scale))


class PartitionEncode(PhaseEncode):
    """
    Through-slab partition encode for 3D imaging.

    Identical arithmetic to :class:`PhaseEncode` with the slab thickness in place of the FOV and
    the slice axis in place of the phase axis.  The separate name keeps 3D recipes readable and
    keeps the parameter names honest -- a partition step is ``1 / slab_thickness``, not
    ``1 / FOV``.

    Parameters
    ----------
    slab_thickness_mm
        Total excited slab thickness, millimetres.
    matrix_sl
        Number of partitions.

    Build arguments
    ---------------
    partition : int
        Signed partition index from the centre of k-space.
    scale : float, default 1.0

    Examples
    --------
    >>> import seqcraft as sc
    >>> par = PartitionEncode(sc.System.preset('generic_3t'), slab_thickness_mm=80, matrix_sl=16)
    >>> round(par.dk_per_m, 2)                  # 1 / 80 mm
    12.5
    >>> par.build(partition=4)
    LogicBlock(par, 1 node, 0.26 ms)
    """

    def __init__(
        self,
        system: System,
        *,
        slab_thickness_mm: float,
        matrix_sl: int,
        axis: str = 'z',
        min_duration_us: float | None = None,
        regime: str = 'default',
    ) -> None:
        super().__init__(
            system,
            fov_pe_mm=slab_thickness_mm,
            matrix_pe=matrix_sl,
            axis=axis,
            min_duration_us=min_duration_us,
            regime=regime,
        )
        self.slab_thickness_mm = float(slab_thickness_mm)
        self.matrix_sl = int(matrix_sl)

    def build(self, *, partition: float = 0.0, scale: float = 1.0) -> LogicBlock:
        """Return the partition encode for `partition`."""
        return self._block('par', self.area_for(partition) * float(scale))


class Prephaser(_AreaTrapezoid):
    """
    A gradient of an explicitly given area: readout prephaser, rewinder, or moment nuller.

    Parameters
    ----------
    area_per_m
        The area to produce, 1/m.  Negative is normal -- a readout prephaser is the negative
        half of the readout's own area.
    axis
        Which logical axis.
    min_duration_us
        Optional duration floor.

    Build arguments
    ---------------
    scale : float, default 1.0
        Multiplier on the designed area.

    Notes
    -----
    This is what replaces the old ``Prewinder``/``Rewinder`` pair.  A winder is not a distinct
    physical object -- it is a gradient of a chosen area placed before or after something else --
    and the tree already says where things go, so one module covers both.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> gx = pp.make_trapezoid('x', flat_area=256.0, flat_time=3.2e-3, system=system.default)
    >>> pre = Prephaser(system, area_per_m=-gx.area / 2, axis='x')
    >>> round(pre.area_per_m, 1)
    -128.8
    >>> pre.build()
    LogicBlock(prephaser, 1 node, 0.30 ms)
    """

    def __init__(
        self,
        system: System,
        *,
        area_per_m: float,
        axis: str = 'x',
        min_duration_us: float | None = None,
        regime: str = 'default',
    ) -> None:
        self.area_per_m = float(area_per_m)
        super().__init__(
            system,
            max_area_per_m=self.area_per_m,
            axes=(axis,),
            min_duration_us=min_duration_us,
            regime=regime,
        )
        self.axis = axis

    def build(self, *, scale: float = 1.0) -> LogicBlock:
        """Return the gradient, scaled by `scale`."""
        return self._block('prephaser', self.area_per_m * float(scale))


class Spoiler(_AreaTrapezoid):
    """
    Gradient spoiler: dephases the residual transverse magnetisation across a voxel.

    Parameters
    ----------
    system
        The scanner.
    twists
        Dephasing cycles across one voxel.  Four or more is usual for a spoiled gradient echo.
    voxel_mm
        The voxel size the twists are counted across.  For a 2D sequence that is the slice
        thickness, which is the dimension the residual signal is coherent over.
    axes
        Which axes carry the spoiler.  A single axis is normal; more spoils faster but adds
        eddy currents on more axes.
    min_duration_us
        Optional duration floor.

    Properties
    ----------
    area_per_m
        ``twists / voxel``, in 1/m.
    duration

    Notes
    -----
    Parameterised in twists per voxel rather than as a bare area because that is the physical
    intent: dephase by N cycles across a voxel, so residual signal cancels within it.  Quoting a
    raw area instead makes the same spoiler wrong the moment the resolution changes.

    Examples
    --------
    >>> import seqcraft as sc
    >>> spoil = Spoiler(sc.System.preset('generic_3t'), twists=4, voxel_mm=5)
    >>> round(spoil.area_per_m)                 # 4 cycles across 5 mm
    800
    >>> spoil.build()
    LogicBlock(spoil, 1 node, 0.74 ms)
    """

    def __init__(
        self,
        system: System,
        *,
        twists: float = 4.0,
        voxel_mm: float,
        axes: tuple[str, ...] = ('z',),
        min_duration_us: float | None = None,
        regime: str = 'default',
    ) -> None:
        self.twists = float(twists)
        self.voxel_mm = float(voxel_mm)
        require_positive(self, 'twists', 'voxel_mm')
        super().__init__(
            system,
            max_area_per_m=self.twists / convert(self.voxel_mm, 'mm', 'm'),
            axes=axes,
            min_duration_us=min_duration_us,
            regime=regime,
        )

    @property
    def area_per_m(self) -> float:
        """Spoiler area, ``twists / voxel``, in 1/m."""
        return self.twists / convert(self.voxel_mm, 'mm', 'm')

    def build(self, *, scale: float = 1.0) -> LogicBlock:
        """Return the spoiler."""
        return self._block('spoil', self.area_per_m * float(scale))


class Crusher(Spoiler):
    """
    Crusher pair lobe, for placing either side of a refocusing pulse.

    Physically a spoiler; the name records the intent, and
    :class:`~seqcraft.modules.rf.pulses.SincRefocusing` can own its own pair instead.  Use this
    when the crushers need to be placed independently -- for instance overlapping a diffusion lobe.

    Parameters
    ----------
    See :class:`Spoiler`: `twists`, `voxel_mm`, `axes` and `min_duration_us` all mean the same
    thing here.

    Notes
    -----
    The two lobes of a crusher pair are **equal**, not opposite: the refocusing pulse inverts the
    phase between them, so an equal pair leaves the wanted echo untouched while dephasing the FID
    an imperfect pulse creates.

    Examples
    --------
    >>> import seqcraft as sc
    >>> crush = Crusher(sc.System.preset('generic_3t'), twists=4, voxel_mm=5)
    >>> crush.build()
    LogicBlock(crush, 1 node, 0.74 ms)
    """

    def build(self, *, scale: float = 1.0) -> LogicBlock:
        """Return one crusher lobe."""
        return self._block('crush', self.area_per_m * float(scale))
