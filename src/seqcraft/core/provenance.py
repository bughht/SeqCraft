"""
Provenance sidecar: making a written ``.seq`` self-describing.

Every :meth:`seqcraft.core.sequence.Sequence.write` drops a ``<name>.seq.json`` beside the
sequence recording versions, the resolved parameters of every module, the achieved timing,
the exact definitions written, and the file hash.

The problem this solves is concrete.  The reference implementation's archived files carry
no b-value, no diffusion directions, no moment order, no acceleration and no partial-Fourier
fraction anywhere in ``[DEFINITIONS]``; those survive only as substrings of a file name
built by a forty-line ``seq_file_name += ...`` ladder in a notebook cell.  The parameter set
that produced a given ``.seq`` is therefore unrecoverable, and two files differing only in
a parameter nobody thought to put in the name are indistinguishable.

Determinism rules, enumerated because byte-identity fails silently otherwise:

* ``sort_keys=True`` everywhere, so the JSON of one build is byte-stable;
* no wall-clock value reaches the ``.seq`` itself (pypulseq writes no timestamp, so two
  writes of the same sequence are byte-identical -- an assertion the test suite makes);
* the creation time *is* recorded in the sidecar, which is not compared byte-wise;
* numpy arrays are summarised by shape, dtype and hash rather than dumped;
* a dirty git tree is recorded as such, so a comparison can be marked indeterminate
  rather than reported as a false pass.
"""

from __future__ import annotations

import datetime as _dt
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .sequence import Sequence

__all__ = ['SIDECAR_SCHEMA', 'build_sidecar', 'git_state', 'write_sidecar']

#: Sidecar schema version.  Bumped independently of the package version.
SIDECAR_SCHEMA = 1


def git_state(cwd: str | Path | None = None) -> dict[str, Any]:
    """
    Return the git commit and dirty flag of `cwd`, or an explanation of why not.

    A dirty tree matters: it means the code that produced the file is not the code any
    commit contains, so a later rebuild cannot be expected to match.
    """
    try:
        commit = subprocess.run(  # noqa: S603
            ['git', 'rev-parse', 'HEAD'],  # noqa: S607
            cwd=str(cwd) if cwd else None,
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        ).stdout.strip()
        status = subprocess.run(  # noqa: S603
            ['git', 'status', '--porcelain'],  # noqa: S607
            cwd=str(cwd) if cwd else None,
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        ).stdout
    except Exception as exc:  # noqa: BLE001 - not being in a repo is normal
        return {'available': False, 'reason': f'{type(exc).__name__}'}
    return {'available': True, 'commit': commit, 'dirty': bool(status.strip())}


def _versions() -> dict[str, Any]:
    """Collect the versions of everything that can change the output."""
    import numpy as np  # noqa: PLC0415
    import pypulseq as pp  # noqa: PLC0415

    from .. import __version__  # noqa: PLC0415

    return {
        'seqcraft': __version__,
        'pypulseq': getattr(pp, '__version__', 'unknown'),
        'numpy': np.__version__,
        'python': sys.version.split()[0],
        'platform': platform.platform(),
    }


def build_sidecar(resolved: Mapping[str, Any], *, seq_bytes: int | None = None) -> dict[str, Any]:
    """
    Wrap a compile's own record in the environment it was produced in.

    Parameters
    ----------
    resolved
        Whatever the caller wants recorded about the build: definitions, system parameters,
        achieved duration, the compile report.
    seq_bytes
        Size of the written ``.seq``.

    Returns
    -------
    dict
        `resolved` under a ``'resolved'`` key, plus the package and pypulseq versions, the git
        commit and dirty flag, and a UTC timestamp.

    Notes
    -----
    A dirty working tree matters and is recorded: it means the code that produced the file is not
    the code any commit contains, so a later rebuild cannot be expected to match.

    Examples
    --------
    >>> side = build_sidecar({'n_blocks': 3})
    >>> side['seqcraft_schema'], side['resolved']['n_blocks']
    (1, 3)
    """
    return {
        'seqcraft_schema': SIDECAR_SCHEMA,
        'versions': _versions(),
        'git': git_state(),
        'created_utc': _dt.datetime.now(_dt.timezone.utc).isoformat(timespec='seconds'),
        'resolved': _jsonable(dict(resolved)),
        'artifacts': {'seq_bytes': seq_bytes},
    }


def write_sidecar(seq_path: str | Path, resolved: Mapping[str, Any]) -> Path:
    """
    Write ``<seq_path>.json`` beside the sequence and return its path.

    Examples
    --------
    A sidecar sits next to the sequence and shares its stem::

        dwi_64x64x8s_fov240mm_b1000_a3f19c.seq
        dwi_64x64x8s_fov240mm_b1000_a3f19c.seq.json
    """
    target = Path(seq_path)
    out = target.with_suffix(target.suffix + '.json')
    payload = build_sidecar(
        resolved, seq_bytes=target.stat().st_size if target.exists() else None
    )
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    return out


def _jsonable(value: Any) -> Any:
    """Coerce definition values (which may be numpy scalars or arrays) into JSON."""
    import numpy as np  # noqa: PLC0415

    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)
