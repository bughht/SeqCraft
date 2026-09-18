"""
R3 -- the Pulseq MATLAB ``writeTSE.m`` demo, through the ``.seq`` it writes.

No MATLAB parser.  The source is read for intent; the ``.seq`` is the executable artifact, and
``pypulseq.Sequence.read`` turns it into the same object every other adapter returns.

One property of that path was checked before relying on it, and it is the reason this adapter is
viable at all: ``pulseq`` writes file format **1.5.0**, whose ``[RF]`` section carries ``use``, and
pypulseq 1.5.1 reads it rather than guessing.  Without ``use``, ``calculate_kspacePP`` would not
conjugate k at the refocusing pulses and every k-space measurement here would be quietly wrong.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pypulseq as pp

from ....reference import ReferenceSequence

PROVENANCE = {
    'repo': 'pulseq/pulseq',
    'path': 'matlab/demoSeq/writeTSE.m',
    'license': 'MIT',
    'role': 'authoritative-external',
}


def read(path: str | Path, *, parameters: dict[str, Any] | None = None,
         name: str = 'pulseq_matlab_tse') -> ReferenceSequence:
    """Read a ``.seq`` written by the MATLAB reference."""
    path = Path(path)
    sequence = pp.Sequence()
    sequence.read(str(path))
    version = getattr(sequence, 'version_combined', None)
    return ReferenceSequence(
        name=name,
        sequence=sequence,
        parameters=dict(parameters or {}),
        definitions=dict(getattr(sequence, 'definitions', {}) or {}),
        provenance=PROVENANCE | {'seq_file': str(path), 'file_version_combined': version},
        # `TEeff` in the source is a *request*; what the file achieves has to be measured.
        semantic={'claims': 'none -- a .seq states definitions, not intent'},
    )
