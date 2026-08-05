"""
Geometry: the single source of truth for field of view, matrix and sampling pattern.

One frozen dataclass owns FOV, matrix, slice/partition layout, partial Fourier and
acceleration -- and *derives* both the phase-encode table and the ``.seq``
``[DEFINITIONS]`` from the same arithmetic.  That shared derivation is the point.

In the reference implementation, ``write_seq(fov=220, Ny=80, Nx=80, ...)`` accepted the
geometry a second time at write time, separately from the values used to build the
gradients, so the metadata in the file could disagree with what was played.  Concretely it
wrote ``kSpaceCenterLine = Ny/2 = 73.0`` while its own navigator computed the centre line
as ``floor(146/2) - (1 - 0.75) * 146 = 36.5``.  Both cannot be right, and 36.5 is not even
an integer.  Here :attr:`Geometry.kspace_center_line` and :attr:`Geometry.pe_lines` come
from one expression in one index space, and ``Sequence.write()`` takes no geometry
arguments at all.

Index conventions
-----------------
Phase-encode lines are numbered on the **reconstruction grid**: ``0 .. matrix_pe - 1``,
with the k-space centre at ``matrix_pe // 2`` (the Siemens convention).  That is the index
space used for the ``LIN`` label, for ``kSpaceCenterLine``, and for the ``line=`` emit
parameter of :class:`~seqcraft.modules.encoding.phase_encode.PhaseEncode`.  Partial Fourier
truncates the *low* end; acceleration then skips forward, and the start index is nudged up
so the centre line is always sampled.

Examples
--------
>>> g = Geometry(fov_mm=(220, 220, 5), matrix=(128, 128, 1), slice_thickness_mm=5)
>>> g.fov_m
(0.22, 0.22, 0.005)
>>> g.dk_per_m[1] == 1 / 0.22
True
>>> g.kspace_center_line
64
>>> len(g.pe_lines)
128

Partial Fourier 0.75 with twofold acceleration on a 146 matrix reproduces the reference
implementation's own numbers exactly:

>>> g2 = Geometry(fov_mm=(220, 220, 15), matrix=(146, 146, 1), slice_thickness_mm=1.5,
...               n_slices=10, slice_gap_mm=4.5, partial_fourier_pe=0.75, accel_pe=2)
>>> g2.pe_first_index, g2.n_pe_acquired, g2.kspace_center_line
(37, 55, 73)
>>> g2.kspace_center_line in g2.pe_lines          # the centre line is always sampled
True
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal

from .errors import ConfigurationError
from .validate import Range, check_fields, require_in, require_in_range

if TYPE_CHECKING:
    from collections.abc import Sequence as Seq

__all__ = ['Geometry', 'round_half_up']

_FOV_RANGE = Range(0.5, 2000.0, 'mm', ((1e3, 'm'), (10.0, 'cm')))
_THICKNESS_RANGE = Range(0.05, 200.0, 'mm', ((1e3, 'm'),))


def round_half_up(x: float) -> int:
    """
    Round half away from zero, unlike Python's banker's rounding.

    ``round(109.5)`` is 110 in Python only by luck of the float representation;
    ``round(0.5)`` is 0.  Partial-Fourier line counts must not depend on that.

    Examples
    --------
    >>> round_half_up(0.5), round_half_up(1.5), round_half_up(109.5)
    (1, 2, 110)
    >>> round_half_up(-0.5)
    -1
    """
    return int(math.floor(x + 0.5)) if x >= 0 else -int(math.floor(-x + 0.5))


@dataclass(frozen=True)
class Geometry:
    """
    Field of view, matrix, slice layout and sampling pattern.

    Parameters
    ----------
    fov_mm
        Field of view in millimetres along the logical axes ``(RO, PE, SL)``.  For 2D
        multislice, the SL entry is the total slab covered by the slice stack; it is
        derived automatically if left as ``0``.
    matrix
        Matrix size ``(RO, PE, SL)``.  The SL entry is the number of partitions; it must
        be 1 for ``mode='2d'``.
    slice_thickness_mm
        Slice (2D) or partition (3D) thickness in millimetres.
    n_slices
        Number of 2D slices.  Must be 1 for ``mode='3d'``, where ``matrix[2]`` gives the
        partition count instead.
    slice_gap_mm
        Centre-to-centre spacing minus thickness.  ``0`` means contiguous slices.
    partial_fourier_pe, partial_fourier_sl
        Fraction of k-space sampled along PE / partition, in ``(0.5, 1.0]``.
    accel_pe, accel_sl
        Undersampling (skip) factors along PE / partition.
    n_shots
        Number of interleaved shots (segments) the PE table is split across.
    ro_oversampling
        Readout oversampling factor; increases ADC samples without changing k-space
        extent.
    mode
        ``'2d'`` for multislice, ``'3d'`` for slab-selective partition encoding.

    Notes
    -----
    Nothing here knows about gradients or timing.  Modules that genuinely need FOV and
    matrix to size a gradient (readouts, phase encoders) take those as their **own**
    fields, and :meth:`check_module` cross-validates them, so duplication is checked
    rather than trusted.
    """

    fov_mm: tuple[float, float, float]
    matrix: tuple[int, int, int]
    slice_thickness_mm: float = field(metadata={'range': _THICKNESS_RANGE})
    n_slices: int = 1
    slice_gap_mm: float = 0.0
    partial_fourier_pe: float = 1.0
    partial_fourier_sl: float = 1.0
    accel_pe: int = 1
    accel_sl: int = 1
    n_shots: int = 1
    ro_oversampling: float = 1.0
    mode: Literal['2d', '3d'] = '2d'

    # ------------------------------------------------------------------- validation
    def __post_init__(self) -> None:
        """Validate ranges, mode consistency, and the acceleration/matrix relationship."""
        if len(self.fov_mm) != 3 or len(self.matrix) != 3:
            msg = 'Geometry.fov_mm and Geometry.matrix must both have three entries (RO, PE, SL).'
            raise ConfigurationError(msg)
        check_fields(self)
        for axis, value in zip('RO PE SL'.split(), self.fov_mm):
            if value > 0:
                from .validate import check_value  # noqa: PLC0415

                check_value('Geometry', f'fov_mm[{axis}]', float(value), _FOV_RANGE)
        require_in(self, 'mode', ('2d', '3d'))
        require_in_range(self, 'partial_fourier_pe', 0.5, 1.0)
        require_in_range(self, 'partial_fourier_sl', 0.5, 1.0)
        require_in_range(self, 'ro_oversampling', 1.0, 16.0)
        for name in ('accel_pe', 'accel_sl', 'n_shots', 'n_slices'):
            if getattr(self, name) < 1:
                msg = f'Geometry.{name} must be >= 1, got {getattr(self, name)}.'
                raise ConfigurationError(msg)
        for i, name in enumerate(('matrix_ro', 'matrix_pe', 'matrix_sl')):
            if self.matrix[i] < 1:
                msg = f'Geometry.{name} must be >= 1, got {self.matrix[i]}.'
                raise ConfigurationError(msg)
        if self.mode == '2d' and self.matrix[2] != 1:
            msg = (
                f"Geometry.mode='2d' requires matrix[2] == 1, got {self.matrix[2]}. "
                f"Use mode='3d' for partition encoding."
            )
            raise ConfigurationError(msg)
        if self.mode == '3d' and self.n_slices != 1:
            msg = (
                f"Geometry.mode='3d' requires n_slices == 1 (partitions live in matrix[2]), "
                f'got {self.n_slices}.'
            )
            raise ConfigurationError(msg)
        if not self.pe_lines:
            msg = (
                f'the phase-encode table is empty for matrix_pe={self.matrix[1]}, '
                f'partial_fourier_pe={self.partial_fourier_pe}, accel_pe={self.accel_pe}.'
            )
            raise ConfigurationError(msg)

    # ------------------------------------------------------------ derived: geometry
    @property
    def fov_m(self) -> tuple[float, float, float]:
        """FOV in metres along (RO, PE, SL). The SL entry is derived when given as 0."""
        ro, pe, sl = self.fov_mm
        if sl <= 0:
            sl = self.slab_thickness_mm
        return (ro / 1e3, pe / 1e3, sl / 1e3)

    @property
    def dk_per_m(self) -> tuple[float, float, float]:
        """k-space sample spacing 1/FOV along (RO, PE, SL), in 1/m."""
        return tuple(1.0 / f for f in self.fov_m)  # type: ignore[return-value]

    @property
    def res_mm(self) -> tuple[float, float, float]:
        """Nominal voxel size in millimetres. The SL entry is the slice thickness."""
        ro, pe, _ = self.fov_mm
        return (ro / self.matrix[0], pe / self.matrix[1], self.slice_thickness_mm)

    @property
    def n_partitions(self) -> int:
        """Number of 3D partitions (1 for 2D multislice)."""
        return self.matrix[2]

    @property
    def slab_thickness_mm(self) -> float:
        """Total extent covered by the slice stack or the 3D slab."""
        if self.mode == '3d':
            return self.slice_thickness_mm * self.n_partitions
        pitch = self.slice_thickness_mm + self.slice_gap_mm
        return self.slice_thickness_mm + pitch * (self.n_slices - 1)

    @property
    def slice_positions_m(self) -> tuple[float, ...]:
        """
        Slice centre offsets along SL in metres, symmetric about zero.

        Examples
        --------
        >>> g = Geometry(fov_mm=(220, 220, 0), matrix=(64, 64, 1),
        ...              slice_thickness_mm=5, n_slices=3)
        >>> [round(z * 1e3, 3) for z in g.slice_positions_m]
        [-5.0, 0.0, 5.0]
        """
        if self.mode == '3d':
            return (0.0,)
        pitch_mm = self.slice_thickness_mm + self.slice_gap_mm
        first = -(self.n_slices - 1) / 2.0
        # Scale to metres last: (i * pitch_mm) / 1e3 is exact for the common integer
        # millimetre pitches, whereas i * (pitch_mm / 1e3) accumulates float noise that
        # then shows up in the SlicePositions definition.
        return tuple((first + i) * pitch_mm / 1e3 for i in range(self.n_slices))

    def slice_position_m(self, index: int) -> float:
        """Offset of slice `index` along SL, in metres."""
        return self.slice_positions_m[index]

    # ------------------------------------------------- derived: phase-encode table
    @property
    def kspace_center_line(self) -> int:
        """
        Recon-grid index of the k-space centre line: ``matrix_pe // 2``.

        This same value is written as the ``kSpaceCenterLine`` definition and used for the
        ``LIN`` label, so the two cannot disagree.
        """
        return self.matrix[1] // 2

    @property
    def kspace_center_partition(self) -> int:
        """Recon-grid index of the k-space centre partition: ``matrix_sl // 2``."""
        return self.matrix[2] // 2

    @property
    def pe_first_index(self) -> int:
        """
        Recon-grid index of the first acquired phase-encode line.

        Partial Fourier truncates the low end; the start is then nudged **up** by the
        residue ``(centre - skip) mod accel`` so the centre line is always sampled --
        without that, an even skip and an odd centre would miss k=0 entirely.
        """
        n, r = self.matrix[1], self.accel_pe
        c = self.kspace_center_line
        skip = n - round_half_up(self.partial_fourier_pe * n)
        skip += (c - skip) % r
        return skip

    @property
    def pe_lines(self) -> tuple[int, ...]:
        """Every acquired phase-encode line, in ascending recon-grid order."""
        return tuple(range(self.pe_first_index, self.matrix[1], self.accel_pe))

    @property
    def n_pe_acquired(self) -> int:
        """Number of acquired phase-encode lines across all shots."""
        return len(self.pe_lines)

    def pe_lines_for_shot(self, shot: int) -> tuple[int, ...]:
        """
        The subset of :attr:`pe_lines` belonging to interleaved `shot`.

        Shots interleave, so shot 0 takes lines 0, n_shots, 2*n_shots, ... of the table.

        Examples
        --------
        >>> g = Geometry(fov_mm=(220, 220, 5), matrix=(8, 8, 1),
        ...              slice_thickness_mm=5, n_shots=2)
        >>> g.pe_lines_for_shot(0), g.pe_lines_for_shot(1)
        ((0, 2, 4, 6), (1, 3, 5, 7))
        """
        if not 0 <= shot < self.n_shots:
            msg = f'shot must be in [0, {self.n_shots - 1}], got {shot}.'
            raise ConfigurationError(msg)
        return self.pe_lines[shot :: self.n_shots]

    @property
    def par_first_index(self) -> int:
        """Recon-grid index of the first acquired partition (3D)."""
        n, r = self.matrix[2], self.accel_sl
        c = self.kspace_center_partition
        skip = n - round_half_up(self.partial_fourier_sl * n)
        skip += (c - skip) % r
        return skip

    @property
    def par_lines(self) -> tuple[int, ...]:
        """Every acquired partition index, ascending (``(0,)`` for 2D)."""
        if self.mode == '2d':
            return (0,)
        return tuple(range(self.par_first_index, self.matrix[2], self.accel_sl))

    # -------------------------------------------------------------------- reporting
    def definitions(self) -> dict[str, Any]:
        """
        The ``.seq`` ``[DEFINITIONS]`` this geometry is responsible for.

        Merged by ``Sequence.write()`` with each readout's own definitions through a
        collision-checking merge, so a disagreement raises instead of being overwritten.
        """
        defs: dict[str, Any] = {
            'FOV': list(self.fov_m),
            'SliceThickness': self.slice_thickness_mm / 1e3,
            'SliceGap': self.slice_gap_mm / 1e3,
            'SlicePositions': list(self.slice_positions_m),
            'BaseResolution': self.matrix[0],
            'PhaseResolution': self.matrix[1] / self.matrix[0],
            'kSpaceCenterLine': self.kspace_center_line,
            'ReadoutOversamplingFactor': self.ro_oversampling,
        }
        if self.mode == '3d':
            defs['kSpaceCenterPartition'] = self.kspace_center_partition
        return defs

    def check_module(self, module: object, *, fov_field: str, matrix_field: str, axis: int) -> None:
        """
        Cross-validate a module's own FOV/matrix fields against this geometry.

        Readouts and phase encoders need FOV and matrix to size their gradients, so those
        values necessarily appear twice.  Duplication that is *checked* is safe;
        duplication that is unchecked is how the reference implementation's metadata drifted
        from what it played.

        Raises
        ------
        ConfigurationError
            On any mismatch, naming both values.
        """
        want_fov = self.fov_mm[axis]
        got_fov = getattr(module, fov_field)
        want_n = self.matrix[axis]
        got_n = getattr(module, matrix_field)
        label = 'RO PE SL'.split()[axis]
        if abs(float(got_fov) - float(want_fov)) > 1e-9:
            msg = (
                f'{type(module).__name__}.{fov_field} = {got_fov:g} mm disagrees with '
                f'Geometry.fov_mm[{label}] = {want_fov:g} mm.'
            )
            raise ConfigurationError(msg)
        if int(got_n) != int(want_n):
            msg = (
                f'{type(module).__name__}.{matrix_field} = {got_n} disagrees with '
                f'Geometry.matrix[{label}] = {want_n}.'
            )
            raise ConfigurationError(msg)

    def params(self) -> dict[str, Any]:
        """Return a JSON-safe description, for the provenance sidecar."""
        return {
            'fov_mm': list(self.fov_mm),
            'matrix': list(self.matrix),
            'slice_thickness_mm': self.slice_thickness_mm,
            'n_slices': self.n_slices,
            'slice_gap_mm': self.slice_gap_mm,
            'partial_fourier_pe': self.partial_fourier_pe,
            'partial_fourier_sl': self.partial_fourier_sl,
            'accel_pe': self.accel_pe,
            'accel_sl': self.accel_sl,
            'n_shots': self.n_shots,
            'ro_oversampling': self.ro_oversampling,
            'mode': self.mode,
            'derived': {
                'res_mm': list(self.res_mm),
                'slab_thickness_mm': self.slab_thickness_mm,
                'pe_first_index': self.pe_first_index,
                'n_pe_acquired': self.n_pe_acquired,
                'kspace_center_line': self.kspace_center_line,
            },
        }

    def describe(self) -> str:
        """Return a human-readable multi-line summary."""
        ro, pe, sl = self.res_mm
        lines = [
            f'Geometry {self.mode}  {self.matrix[0]}x{self.matrix[1]}'
            + (f'x{self.matrix[2]}' if self.mode == '3d' else f'  {self.n_slices} slices'),
            f'  FOV        {self.fov_mm[0]:g} x {self.fov_mm[1]:g} mm'
            f'   slab {self.slab_thickness_mm:g} mm',
            f'  resolution {ro:.3g} x {pe:.3g} x {sl:.3g} mm',
        ]
        if self.partial_fourier_pe < 1 or self.accel_pe > 1 or self.n_shots > 1:
            lines.append(
                f'  sampling   PF {self.partial_fourier_pe:g}  R {self.accel_pe}'
                f'  shots {self.n_shots}  ->  {self.n_pe_acquired} lines'
                f'  (first {self.pe_first_index}, centre {self.kspace_center_line})'
            )
        return '\n'.join(lines)
