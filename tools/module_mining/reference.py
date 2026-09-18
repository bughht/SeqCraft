"""
One record per runnable reference implementation, so later code can ask the same questions of
a SeqCraft tree, a PyPulseq script and a MATLAB-written ``.seq``.

Deliberately thin, and deliberately not frozen.  An adapter's only job is to produce a
``pypulseq.Sequence`` plus the provenance and the claims that go with it; it must not repair,
normalise or reinterpret what the source does.  Where a source is odd, the oddity is evidence.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pypulseq as pp


@dataclass
class ReferenceSequence:
    """A reference implementation, executed, with what it claims about itself."""

    #: Short identifier used in reports, e.g. ``'seqcraft_fse'``.
    name: str

    #: The executed sequence.  Every measurement is taken from this, not from the source.
    sequence: pp.Sequence

    #: The protocol the reference was run at, as the adapter passed it.
    parameters: dict[str, Any] = field(default_factory=dict)

    #: ``seq.definitions`` if the reference sets any.
    definitions: dict[str, Any] = field(default_factory=dict)

    #: Repository, commit, path, license, and any caveat.  Frozen by the adapter, not inferred.
    provenance: dict[str, Any] = field(default_factory=dict)

    #: What the implementation *claims*: echo spacing, effective TE, which sample is the echo,
    #: the phase-encode table.  Compared against what is measured, rather than trusted.
    semantic: dict[str, Any] = field(default_factory=dict)

    #: The SeqCraft ``LogicBlock`` when there is one.  Keeps the exact event digest (L1) reachable
    #: for SeqCraft-to-SeqCraft comparisons; ``None`` for external references.
    tree: Any = None
