"""
The exception root, and the shared error-message formatter.

seqcraft signals failure by raising, never by printing.  The reference implementation
this package replaces (``pSeq_Base``) contained no ``assert`` or ``raise`` at all: its
``get_report()`` printed the result of ``check_timing()`` and returned ``None``, so a
failing sequence could not be detected programmatically.

**Only what more than one package raises lives here.**  Everything else stays with the code that
raises it -- :class:`~seqcraft.compiler.errors.CompileError` and its two siblings in
``compiler/``, :class:`~seqcraft.design.timing.RasterError` in ``design/timing.py``,
:class:`~seqcraft.scanner.opts.UnknownFieldError` in ``scanner/opts.py``.  Every one of them is
re-exported as ``sc.<Name>``, so the spelling a caller writes never depends on the layout::

    SeqCraftError                  # the base; catch this to catch everything
    +-- ConfigurationError         # a call is wrong: bad type, unknown unit, unusable Opts
    |   +-- RasterError            #   design/timing.py
    |   +-- UnknownFieldError      #   scanner/opts.py
    +-- CompileError               # compiler/errors.py
    +-- HardwareLimitError         # compiler/errors.py
    +-- DefinitionConflict         # compiler/errors.py
    +-- MissingExtraError          # an optional dependency is needed; names the extra

Every message follows one shape so it is scannable and so a test can assert on it::

    <ErrorClass>: <one line: what is wrong, with the number and unit>
      <key>:  <value>
      <key>:  <value>
      fix
        <option 1, with the numbers already filled in>
        <option 2>

Use :func:`format_error` to build the body.  What the compiler *did* rather than refused is not an
exception at all: it is a :class:`SeqCraftWarning`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping, Sequence

__all__ = [
    'ConfigurationError',
    'MissingExtraError',
    'SeqCraftError',
    'SeqCraftWarning',
    'format_error',
]


class SeqCraftError(Exception):
    """Base class for every error raised by seqcraft."""


class SeqCraftWarning(UserWarning):
    """
    The compiler changed a waveform or a time to make the sequence legal.

    Not a failure: it names something the compile *did* -- summed two gradients that shared an
    axis, resampled one onto the raster to join a boundary, snapped a reservation to the block
    raster -- rather than something it refused.  A ``UserWarning`` subclass, so the standard
    ``warnings`` machinery applies::

        with warnings.catch_warnings():
            warnings.simplefilter('error', SeqCraftWarning)   # treat any as fatal
            seq = sc.compile(tree, opts)
    """


class ConfigurationError(SeqCraftError):
    """
    A parameter is missing, out of range, or inconsistent with another parameter.

    Two packages subclass it, each keeping its subclass beside the code that raises it:
    :class:`seqcraft.design.timing.RasterError` and
    :class:`seqcraft.scanner.opts.UnknownFieldError`.
    """


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
