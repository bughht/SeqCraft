"""
``Geometry``.  **Not part of the seqcraft package** -- see ``salvage/README.md``.

FOV, matrix, slice layout, partial Fourier and acceleration, in one frozen dataclass that derives
the ``.seq`` ``[DEFINITIONS]`` from them -- plus the unit-plausibility bands that catch a
millimetre passed as a metre.

Why it is not in the package
----------------------------
It never was on the compile path.  ``compile_sequence`` took ``geometry=`` only to call
``definitions()`` on it and merge the resulting eight keys, so ~450 lines of dataclass and range
framework sat inside the package to produce one dict.  The compiler now takes those keys directly::

    sc.compile(tree, opts, definitions=geometry.definitions())

which is the same information with none of the coupling, and it means a geometry can be *anything*
that produces a mapping -- a dataclass of your own, a dict read from JSON, a protocol object.

The deeper reason is that geometry is an **application** concern rather than a compiler one.  FOV,
matrix and slice order are decisions about the scan you are running; the compiler's job is to turn
a tree into legal pulseq blocks and it is indifferent to why the tree looks the way it does.
Keeping ``Geometry`` in the package made it look like a required input, which it never was.

What is worth keeping
---------------------
**One index space, stated once.**  In the reference implementation ``write_seq(fov=220, Ny=80,
...)`` accepted the geometry a second time at write time, separately from the values the gradients
were sized from, so the metadata in the file could disagree with what was played.  Concretely it
wrote ``kSpaceCenterLine = Ny/2 = 73.0`` while its own navigator computed the centre line as
``floor(146/2) - (1 - 0.75) * 146 = 36.5``.  Both cannot be right, and 36.5 is not an integer.
Deriving the definitions from the same fields the gradients are sized from is what prevents that.

**The plausibility bands.**  The failure they exist for is a *plausible number in the wrong unit* --
``fov=220`` meaning millimetres reaching code that reads metres.  That is not a type error, so no
type checker or coercion layer catches it; a range check does, and because every public float field
carries its unit as a suffix, :data:`DEFAULT_RANGES` covers most fields with no per-field code at
all.  It also says what the value probably meant::

    Geometry.fov_mm = 0.22 is outside the plausible range 5 .. 600 mm.
      got :  0.22
      hint:  0.22 looks like m. Did you mean fov_mm=220?

One test went with them and is worth re-establishing wherever this lands: every unit *name* in
:data:`DEFAULT_RANGES`, and every alias name, must be one :func:`seqcraft.convert` knows.  Without
it the two vocabularies drift, and an error message ends up quoting a unit the reader cannot pass
back in.  It was a four-line set comprehension over ``DEFAULT_RANGES`` against ``known_units()``.

Index conventions
-----------------
Phase-encode lines are numbered on the **reconstruction grid**: ``0 .. matrix_pe - 1``, with the
k-space centre at ``matrix_pe // 2`` (the Siemens convention).  That is the index space of the
``LIN`` label and of ``kSpaceCenterLine``.  Which lines an accelerated or partial-Fourier
acquisition actually takes is in ``geometry_pe.py`` beside this file.

Standalone on purpose: this module imports nothing from seqcraft, so it can be copied into a module
library or into user code as it stands.  Its two error types are plain ``ValueError`` subclasses
rather than ``seqcraft.errors`` ones for the same reason.

Where it is going
-----------------
A geometry of roughly this shape is wanted again when the module library gets its basic
infrastructure: a readout and a phase encoder both need FOV and matrix to size a gradient, and the
whole argument for deriving the definitions from the same fields is that the two cannot then
disagree.  This is the design to start from -- including the range bands, which are what turn a
metre passed as a millimetre into a sentence rather than a wrong image.  It is parked here rather
than kept in the package because that infrastructure does not exist yet, and a required-looking
parameter with no consumers is what the package has just finished removing.

Examples
--------
>>> g = Geometry(fov_mm=(220, 220, 5), matrix=(128, 128, 1), slice_thickness_mm=5)
>>> g.fov_m
(0.22, 0.22, 0.005)
>>> g.kspace_center_line
64
>>> g.definitions()['BaseResolution']
128

A 2D stack derives its slab from the slice pitch, so ``fov_mm[SL]`` may be left at zero:

>>> g2 = Geometry(fov_mm=(220, 220, 0), matrix=(146, 146, 1), slice_thickness_mm=1.5,
...               n_slices=10, slice_gap_mm=4.5, partial_fourier_pe=0.75, accel_pe=2)
>>> round(g2.slab_thickness_mm, 3), g2.kspace_center_line
(55.5, 73)

And the whole point of it -- what to hand the compiler:

>>> len(g.definitions()), g.definitions()['kSpaceCenterLine']
(8, 64)
"""

from __future__ import annotations

import dataclasses
import difflib
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal, NamedTuple

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping, Sequence

__all__ = ['ConfigurationError', 'Geometry', 'UnitSanityError']

#: The proton gyromagnetic ratio in Hz/T, for the gradient alias hints.  Spelled out rather than
#: imported, so this file stands alone.
GAMMA_1H = 42_576_384.74


class ConfigurationError(ValueError):
    """A parameter is missing, out of range, or inconsistent with another parameter."""


class UnitSanityError(ConfigurationError):
    """A value is outside the plausible range for its named unit (mm vs m, ms vs s, ...)."""


def format_error(
    headline: str,
    fields: Mapping[str, object] | None = None,
    fixes: Sequence[str] = (),
) -> str:
    """
    Build a message body in seqcraft's house shape: headline, aligned fields, then fixes.

    A four-line copy of :func:`seqcraft.errors.format_error`, so this module needs no import.
    """
    lines = [headline]
    if fields:
        width = max(len(str(k)) for k in fields)
        lines.extend(f'  {k!s:<{width}}:  {v}' for k, v in fields.items())
    if fixes:
        lines.append('  fix')
        lines.extend(f'    {f}' for f in fixes)
    return '\n'.join(lines)


class Range(NamedTuple):
    """
    A plausible interval for one named quantity, plus the units it is confused with.

    Parameters
    ----------
    lo, hi
        Inclusive bounds, expressed in `unit`.
    unit
        The unit the field is named for, e.g. ``'mm'``.
    aliases
        ``(factor, name)`` pairs.  A value `v` is suspected of being in unit `name` when
        ``v * factor`` lands inside ``[lo, hi]``.  Used only to build the hint.
    """

    lo: float
    hi: float
    unit: str
    aliases: tuple[tuple[float, str], ...] = ()


# 1 mT/m in Hz/m, for the gradient alias hints.  Approximate on purpose: the alias only
# has to land the suspected value inside the band, not reproduce it exactly.
_MT_PER_M_IN_HZ_PER_M = 1e-3 * GAMMA_1H

#: Ranges inferred from a field's unit suffix.  **Longest suffix wins**, which is what keeps
#: the reciprocal units separate from the plain ones: ``area_per_m`` is a k-space area in
#: 1/m and must not be matched by ``_m`` (metres), and ``b_value_s_per_mm2`` must not be
#: matched by ``_s`` (seconds).  Override per field with ``metadata={'range': Range(...)}``.
DEFAULT_RANGES: dict[str, Range] = {
    '_s_per_mm2': Range(0.0, 30_000.0, 's/mm^2', ((1e-6, 's/m^2'),)),
    '_Hz_per_m_per_s': Range(1e6, 1e11, 'Hz/m/s', ((GAMMA_1H, 'T/m/s'),)),
    '_Hz_per_m': Range(1e4, 4e7, 'Hz/m', ((_MT_PER_M_IN_HZ_PER_M, 'mT/m'),)),
    # Gradient moments and k-space positions.  Signed: a prewinder area is negative.
    '_per_m_per_s': Range(-1e9, 1e9, '1/m/s', ()),
    '_per_m2': Range(-1e12, 1e12, '1/m^2', ()),
    '_per_m': Range(-1e7, 1e7, '1/m', ()),
    '_mT_m': Range(0.5, 500.0, 'mT/m', ((1.0 / _MT_PER_M_IN_HZ_PER_M, 'Hz/m'),)),
    '_T_m_s': Range(1.0, 1000.0, 'T/m/s', ((1e-3, 'mT/m/ms'),)),
    '_mm': Range(0.01, 1000.0, 'mm', ((1e3, 'm'), (10.0, 'cm'))),
    '_um': Range(1.0, 1e6, 'um', ((1e6, 'm'),)),
    '_cm': Range(0.01, 100.0, 'cm', ((100.0, 'm'),)),
    '_ms': Range(1e-3, 60_000.0, 'ms', ((1e3, 's'), (1e-3, 'us'))),
    '_us': Range(0.1, 6e7, 'us', ((1e6, 's'), (1e3, 'ms'))),
    '_deg': Range(-720.0, 720.0, 'deg', ((57.29578, 'rad'),)),
    '_rad': Range(-13.0, 13.0, 'rad', ((0.0174533, 'deg'),)),
    '_uT': Range(0.1, 1000.0, 'uT', ((1e6, 'T'),)),
    '_m': Range(1e-5, 3.0, 'm', ((1e-3, 'mm'), (1e-2, 'cm'))),
    '_s': Range(1e-9, 6000.0, 's', ((1e-3, 'ms'), (1e-6, 'us'))),
    '_Hz': Range(-1e9, 1e9, 'Hz', ()),
}

#: Field names that are legitimately dimensionless.  Kept here rather than in the test
#: so that adding an exemption is a visible diff in library code, not in test code.
DIMENSIONLESS: frozenset[str] = frozenset({
    'accel',
    'accel_pe',
    'accel_sl',
    'apodization',
    'center_pos',
    'matrix',
    'matrix_pe',
    'matrix_ro',
    'matrix_sl',
    'moment_order',
    'n_adc',
    'n_averages',
    'n_bands',
    'n_echoes',
    'n_interleaves',
    'n_lines',
    'n_partitions',
    'n_samples',
    'n_shots',
    'n_slices',
    'n_twists',
    'oversampling',
    'partial_fourier',
    'ro_oversampling',
    'samples',
    'time_bw_product',
    'twists',
})

_SUFFIXES = tuple(sorted(DEFAULT_RANGES, key=len, reverse=True))


def _range_for(name: str, explicit: Range | None = None) -> Range | None:
    """
    Return the applicable :class:`Range` for a field name, or ``None`` if unconstrained.

    Examples
    --------
    >>> _range_for('te_ms').unit
    'ms'
    >>> _range_for('b_value_s_per_mm2').unit
    's/mm^2'
    >>> _range_for('n_slices') is None
    True
    """
    if explicit is not None:
        return explicit
    if name in DIMENSIONLESS:
        return None
    for suffix in _SUFFIXES:
        if name.endswith(suffix):
            return DEFAULT_RANGES[suffix]
    return None


def _hint(value: float, rng: Range, name: str) -> str | None:
    """Build a 'looks like <unit>' hint by testing each alias factor."""
    for factor, alias in rng.aliases:
        converted = value * factor
        if rng.lo <= converted <= rng.hi:
            shown = f'{converted:g}'
            return f'{value:g} looks like {alias}. Did you mean {name}={shown}?'
    return None


def _check_value(owner: str, name: str, value: float, rng: Range) -> None:
    """
    Raise :class:`~geometry.UnitSanityError` if `value` is outside `rng`.

    Examples
    --------
    >>> _check_value('EpiReadout', 'fov_mm', 220, Range(5.0, 600.0, 'mm'))
    >>> _check_value('EpiReadout', 'fov_mm', 0.22, Range(5.0, 600.0, 'mm', ((1e3, 'm'),)))
    Traceback (most recent call last):
        ...
    geometry.UnitSanityError: EpiReadout.fov_mm = 0.22 is outside ...
    """
    if rng.lo <= value <= rng.hi:
        return
    fields: dict[str, object] = {'got': f'{value:g}'}
    hint = _hint(value, rng, name)
    if hint:
        fields['hint'] = hint
    fields['note'] = 'seqcraft never auto-converts - a wrong unit is a wrong sequence.'
    msg = format_error(
        f'{owner}.{name} = {value:g} is outside the plausible range '
        f'{rng.lo:g} .. {rng.hi:g} {rng.unit}.',
        fields,
    )
    raise UnitSanityError(msg)


def _check_fields(obj: object) -> None:
    """
    Validate every numeric field of a dataclass instance against its range.

    Fields that are ``None``, boolean,
    or non-numeric are skipped; numpy arrays are skipped (they are waveforms, validated
    elsewhere).

    Parameters
    ----------
    obj
        Any dataclass instance.

    Raises
    ------
    UnitSanityError
        On the first field found outside its plausible range.
    """
    if not dataclasses.is_dataclass(obj):
        return
    owner = type(obj).__name__
    for f in dataclasses.fields(obj):
        value = getattr(obj, f.name, None)
        if value is None or isinstance(value, bool) or not isinstance(value, (int, float)):
            continue
        # Exactly zero is never a unit confusion: a wrong unit changes the magnitude, and
        # zero has no magnitude to get wrong.  0 mm is 0 m.  Skipping it lets gaps,
        # offsets and "derive this" sentinels keep their natural default without needing a
        # bespoke range each.  A zero that is genuinely invalid is a range check of its own.
        if value == 0:
            continue
        rng = _range_for(f.name, f.metadata.get('range'))
        if rng is not None:
            _check_value(owner, f.name, float(value), rng)


# ------------------------------------------------------------------ explicit helpers
def _require_in_range(obj: object, name: str, lo: float, hi: float, *, unit: str = '') -> None:
    """Require ``lo <= obj.name <= hi``."""
    owner = type(obj).__name__
    value = getattr(obj, name)
    if value is None:
        return
    if not lo <= value <= hi:
        suffix = f' {unit}' if unit else ''
        msg = f'{owner}.{name} must be in [{lo:g}, {hi:g}]{suffix}, got {value:g}{suffix}.'
        raise ConfigurationError(msg)


def _require_in(obj: object, name: str, allowed: Iterable[Any]) -> None:
    """
    Require the named field to be one of `allowed`, with a close-match suggestion.

    Examples
    --------
    >>> from dataclasses import dataclass
    >>> @dataclass
    ... class D:
    ...     axis: str
    >>> _require_in(D(axis='q'), 'axis', ('x', 'y', 'z'))
    Traceback (most recent call last):
        ...
    geometry.ConfigurationError: D.axis must be one of ('x', 'y', 'z'), got 'q'.
    """
    owner = type(obj).__name__
    value = getattr(obj, name)
    allowed = tuple(allowed)
    if value in allowed:
        return
    msg = f'{owner}.{name} must be one of {allowed!r}, got {value!r}.'
    if isinstance(value, str):
        close = difflib.get_close_matches(value, [str(a) for a in allowed], n=1)
        if close:
            msg += f' Did you mean {close[0]!r}?'
    raise ConfigurationError(msg)


_FOV_RANGE = Range(0.5, 2000.0, 'mm', ((1e3, 'm'), (10.0, 'cm')))
_THICKNESS_RANGE = Range(0.05, 200.0, 'mm', ((1e3, 'm'),))


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
    fields; comparing the two is the module's business, and a check that reached into a
    module by attribute name had no caller and made assumptions about how one is spelled.
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
        _check_fields(self)
        for axis, value in zip('RO PE SL'.split(), self.fov_mm):
            if value > 0:
                _check_value('Geometry', f'fov_mm[{axis}]', float(value), _FOV_RANGE)
        _require_in(self, 'mode', ('2d', '3d'))
        _require_in_range(self, 'partial_fourier_pe', 0.5, 1.0)
        _require_in_range(self, 'partial_fourier_sl', 0.5, 1.0)
        _require_in_range(self, 'ro_oversampling', 1.0, 16.0)
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

    # ------------------------------------------------------------ derived: geometry
    @property
    def fov_m(self) -> tuple[float, float, float]:
        """FOV in metres along (RO, PE, SL). The SL entry is derived when given as 0."""
        ro, pe, sl = self.fov_mm
        if sl <= 0:
            sl = self.slab_thickness_mm
        return (ro / 1e3, pe / 1e3, sl / 1e3)

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

    # ------------------------------------------------------- derived: k-space centre
    @property
    def kspace_center_line(self) -> int:
        """
        Recon-grid index of the k-space centre line: ``matrix_pe // 2``.

        Written as the ``kSpaceCenterLine`` definition.  A phase-encode table has to agree with
        it, which is why ``salvage/geometry_pe.py`` computes its centre with the same expression
        rather than a second one.
        """
        return self.matrix[1] // 2

    @property
    def kspace_center_partition(self) -> int:
        """Recon-grid index of the k-space centre partition: ``matrix_sl // 2``."""
        return self.matrix[2] // 2

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
