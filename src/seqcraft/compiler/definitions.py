"""
Merging ``[DEFINITIONS]`` from several sources, with a collision check.

One function, on the compile path.  A ``.seq`` file's definitions come from three places at once
-- the sequence name, the geometry, and whatever the caller passed -- and the reference
implementation let the last writer win.  That is how a file came to say
``kSpaceCenterLine = 73.0`` while its own navigator used 36.5: neither source was wrong about
itself, and nothing compared them.

Here two sources claiming one key with different values is a :class:`~seqcraft.errors.
DefinitionConflict`, and the message names *who* claimed it twice rather than reporting an
anonymous "already set".  Values that merely differ in representation -- a tuple against a numpy
array, two floats derived by different arithmetic -- are the same definition and merge silently;
that judgement is :func:`_definitions_equal`, and getting it wrong in the strict direction would
stop compiles for no reason.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ..errors import format_error

if TYPE_CHECKING:
    from collections.abc import Mapping

__all__ = ['merge_definitions']


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
    from ..errors import DefinitionConflict  # noqa: PLC0415  (avoids a circular import)

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
