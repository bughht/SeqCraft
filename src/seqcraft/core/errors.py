"""
Exception hierarchy and the shared error-message formatter.

seqcraft signals failure by raising, never by printing.  The reference implementation
this package replaces (``pSeq_Base``) contained no ``assert`` or ``raise`` at all: its
``get_report()`` printed the result of ``check_timing()`` and returned ``None``, so a
failing sequence could not be detected programmatically.

Every message follows one shape so it is scannable and so a test can assert on it::

    <ErrorClass>: <one line: what is wrong, with the number and unit>
      <key>:  <value>
      <key>:  <value>
      fix
        <option 1, with the numbers already filled in>
        <option 2>

Use :func:`format_error` to build the body.  Soft findings are *not* exceptions -- they
belong in a :class:`seqcraft.core.report.Report` as an ``Issue`` so they survive into the
provenance sidecar.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping, Sequence

__all__ = [
    'CompileError',
    'ConfigurationError',
    'DefinitionConflict',
    'HardwareLimitError',
    'MissingExtraError',
    'PurityError',
    'RasterError',
    'SeqCraftError',
    'UnitSanityError',
    'UnknownFieldError',
    'format_error',
]


class SeqCraftError(Exception):
    """Base class for every error raised by seqcraft."""


class ConfigurationError(SeqCraftError):
    """A parameter is missing, out of range, or inconsistent with another parameter."""


class UnitSanityError(ConfigurationError):
    """A value is outside the plausible range for its named unit (mm vs m, ms vs s, ...)."""


class RasterError(ConfigurationError):
    """A duration is not an exact multiple of the raster it must land on."""


class UnknownFieldError(ConfigurationError):
    """``with_()`` or ``override()`` was given a field the target does not declare."""


class HardwareLimitError(SeqCraftError):
    """A gradient exceeds an amplitude or slew limit, per axis or in vector norm."""


class CompileError(SeqCraftError):
    """
    A logic-block tree cannot be expressed as legal pulseq blocks.

    Raised for the conditions the compiler cannot resolve by scheduling: two RF or ADC events
    overlapping in time, a negative absolute start, or a block boundary falling inside a
    gradient that an ADC is sampling.  Amplitude and slew violations are reported through a
    :class:`~seqcraft.core.report.Report` instead, because there are usually several and
    seeing them all at once is more useful than stopping at the first.
    """


class DefinitionConflict(SeqCraftError):
    """Two sources claimed the same ``.seq`` definition key with different values."""


class PurityError(AssertionError):
    """A module mutated itself or a shared event. Raised by the test helpers."""


class MissingExtraError(SeqCraftError):
    """An optional dependency is needed. The message names the extra to install."""


def format_error(
    headline: str,
    fields: Mapping[str, object] | None = None,
    fixes: Iterable[str] = (),
    *,
    sections: Sequence[tuple[str, Mapping[str, object]]] = (),
) -> str:
    """
    Build a message body in the house shape.

    Parameters
    ----------
    headline
        One line stating what is wrong, including the offending number and its unit.
    fields
        ``key -> value`` pairs printed under the headline, aligned.
    fixes
        Concrete remedies, ideally with the resulting numbers already computed.
    sections
        Optional extra named blocks, printed after `fields`; used by the timing
        solvers to show one ledger per constraint chain.

    Returns
    -------
    str
        The formatted message, suitable as the sole argument to an exception.

    Examples
    --------
    >>> print(format_error('TE = 4.00 ms is too short.',
    ...                    {'minimum': '6.34 ms'},
    ...                    ["te_s='min'"]))
    TE = 4.00 ms is too short.
      minimum:  6.34 ms
      fix
        te_s='min'
    """
    lines = [headline]
    blocks: list[tuple[str, Mapping[str, object]]] = []
    if fields:
        blocks.append(('', fields))
    blocks.extend(sections)

    for title, mapping in blocks:
        if title:
            lines.append(f'  {title}')
        if not mapping:
            continue
        width = max(len(str(k)) for k in mapping)
        indent = '    ' if title else '  '
        lines.extend(f'{indent}{k!s:<{width}}:  {v}' for k, v in mapping.items())

    fix_list = list(fixes)
    if fix_list:
        lines.append('  fix')
        lines.extend(f'    {f}' for f in fix_list)
    return '\n'.join(lines)
