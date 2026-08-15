"""
Parameter validation: plausibility ranges, relational checks, and the "did you mean"
hint generator.

The failure this module exists to prevent is a *plausible number in the wrong unit*.  In
the reference implementation, ``write_seq(fov=220, ...)`` defaulted to 220 while callers
passed ``0.22``, and ``pSeq_GrOpt`` kept ``T_90``/``T_180``/``T_readout`` in ms beside
``dt`` in s and ``gmax`` in mT/m inside one dict.  Neither is a *type* error, so no type
checker or coercion layer would catch either.  A range check does, and it can say what
the value probably meant.

Two mechanisms:

1. **Automatic ranges from the field name.**  Because every public float field carries its
   unit as a suffix (see :mod:`seqcraft.core.units`), :data:`DEFAULT_RANGES` covers most
   fields with no per-field code at all.
2. **Explicit ranges via field metadata**, for fields where the default is wrong -- both
   ``fov_mm`` and ``slice_thickness_mm`` end in ``_mm`` but differ by two orders of
   magnitude::

       fov_mm: float = field(metadata={'range': Range(5.0, 600.0, 'mm')})

:func:`check_fields` takes an ordinary dataclass and asks nothing about its type, so a component
of your own gets the check by calling it.

This file also holds :func:`merge_definitions`, which is on the compile path proper: it is how the
compiler merges ``[DEFINITIONS]`` from several sources with a collision check.

Examples
--------
>>> from dataclasses import dataclass, field
>>> BAND = Range(5.0, 600.0, 'mm', ((1e3, 'm'),))     # aliases generate the hint
>>> @dataclass
... class Demo:
...     fov_mm: float = field(metadata={'range': BAND})
>>> check_fields(Demo(fov_mm=220))
>>> check_fields(Demo(fov_mm=0.22))
Traceback (most recent call last):
    ...
seqcraft.core.errors.UnitSanityError: Demo.fov_mm = 0.22 is outside the plausible range 5 .. 600 mm.
  got :  0.22
  hint:  0.22 looks like m. Did you mean fov_mm=220?
  note:  seqcraft never auto-converts - a wrong unit is a wrong sequence.
"""

from __future__ import annotations

import dataclasses
import difflib
from typing import TYPE_CHECKING, Any, NamedTuple

from .errors import ConfigurationError, UnitSanityError, format_error
from .units import GAMMA_1H

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

__all__ = [
    'DEFAULT_RANGES',
    'DIMENSIONLESS',
    'Range',
    'check_fields',
    'check_value',
    'merge_definitions',
    'range_for',
    'require_in',
    'require_in_range',
]


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


def range_for(name: str, explicit: Range | None = None) -> Range | None:
    """
    Return the applicable :class:`Range` for a field name, or ``None`` if unconstrained.

    Examples
    --------
    >>> range_for('te_ms').unit
    'ms'
    >>> range_for('b_value_s_per_mm2').unit
    's/mm^2'
    >>> range_for('n_slices') is None
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


def check_value(owner: str, name: str, value: float, rng: Range) -> None:
    """
    Raise :class:`~seqcraft.core.errors.UnitSanityError` if `value` is outside `rng`.

    Examples
    --------
    >>> check_value('EpiReadout', 'fov_mm', 220, Range(5.0, 600.0, 'mm'))
    >>> check_value('EpiReadout', 'fov_mm', 0.22, Range(5.0, 600.0, 'mm', ((1e3, 'm'),)))
    Traceback (most recent call last):
        ...
    seqcraft.core.errors.UnitSanityError: EpiReadout.fov_mm = 0.22 is outside ...
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


def check_fields(obj: object) -> None:
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
        rng = range_for(f.name, f.metadata.get('range'))
        if rng is not None:
            check_value(owner, f.name, float(value), rng)


# ------------------------------------------------------------------ explicit helpers
def require_in_range(obj: object, name: str, lo: float, hi: float, *, unit: str = '') -> None:
    """Require ``lo <= obj.name <= hi``."""
    owner = type(obj).__name__
    value = getattr(obj, name)
    if value is None:
        return
    if not lo <= value <= hi:
        suffix = f' {unit}' if unit else ''
        msg = f'{owner}.{name} must be in [{lo:g}, {hi:g}]{suffix}, got {value:g}{suffix}.'
        raise ConfigurationError(msg)


def require_in(obj: object, name: str, allowed: Iterable[Any]) -> None:
    """
    Require the named field to be one of `allowed`, with a close-match suggestion.

    Examples
    --------
    >>> from dataclasses import dataclass
    >>> @dataclass
    ... class D:
    ...     axis: str
    >>> require_in(D(axis='q'), 'axis', ('x', 'y', 'z'))
    Traceback (most recent call last):
        ...
    seqcraft.core.errors.ConfigurationError: D.axis must be one of ('x', 'y', 'z'), got 'q'.
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


def merge_definitions(sources: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    """
    Merge ``.seq`` definition dicts, raising on a genuine conflict.

    Parameters
    ----------
    sources
        ``source name -> {definition key: value}``.

    Returns
    -------
    dict
        The merged mapping.

    Raises
    ------
    DefinitionConflict
        If two sources claim the same key with different values.  Last-writer-wins is
        how ``.seq`` metadata came to disagree with what was actually played.
    """
    from .errors import DefinitionConflict  # noqa: PLC0415  (avoids a circular import)

    out: dict[str, Any] = {}
    owner: dict[str, str] = {}
    for source, mapping in sources.items():
        for key, value in mapping.items():
            if key in out and not _definitions_equal(out[key], value):
                msg = format_error(
                    f'definition {key!r} claimed twice with different values.',
                    {owner[key]: out[key], source: value},
                    [f'make {owner[key]} and {source} agree, or stop one of them emitting {key!r}'],
                )
                raise DefinitionConflict(msg)
            out[key] = value
            owner.setdefault(key, source)
    return out


def _definitions_equal(a: Any, b: Any) -> bool:
    """
    Compare two definition values, tolerating float noise and sequence-type mixing.

    A geometry writing ``FOV`` as a tuple and a module writing it as a numpy array are the same
    definition, and so are two floats that differ in the last bit after being derived by
    different arithmetic.  Neither should stop a compile.
    """
    import numpy as np  # noqa: PLC0415  (keeps `import seqcraft` light)

    seq_types = (list, tuple, np.ndarray)
    if isinstance(a, seq_types) or isinstance(b, seq_types):
        if not (isinstance(a, seq_types) and isinstance(b, seq_types)):
            return False
        return len(a) == len(b) and all(_definitions_equal(x, y) for x, y in zip(a, b))
    if isinstance(a, (int, float, np.number)) and isinstance(b, (int, float, np.number)):
        return abs(float(a) - float(b)) <= 1e-12 * max(1.0, abs(float(a)))
    return bool(a == b)
