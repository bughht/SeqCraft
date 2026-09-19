"""
The official PyPulseq radial references, loaded from disk and called.

Two scripts rather than one, because together they answer the question the fine scan exists to
settle:

``write_radial_gre.py``
    a full spoke through the centre of k-space.
``write_ute.py``
    the same family with a ``readout_asymmetry`` argument, which at ``1.0`` is centre-out.

Neither is edited.  Set ``MODULE_MINING_PYPULSEQ`` to choose a checkout; the default is the
official ``imr-framework/pypulseq``.
"""

from __future__ import annotations

import importlib.util
import os
import warnings
from pathlib import Path
from typing import Any

from ....reference import ReferenceSequence

CHECKOUT = Path(os.environ.get('MODULE_MINING_PYPULSEQ', '/Users/yiyund/Code_mgh/pypulseq'))

PROVENANCE = {
    'repo': 'imr-framework/pypulseq',
    'commit': 'f2c582b',
    'license': 'MIT',
    'role': 'authoritative-external',
}


def _load(script: str) -> Any:
    path = CHECKOUT / 'examples' / 'scripts' / script
    spec = importlib.util.spec_from_file_location(f'module_mining_{path.stem}', path)
    if spec is None or spec.loader is None:                           # pragma: no cover
        msg = f'cannot load the reference script at {path}'
        raise FileNotFoundError(msg)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build(script: str = 'write_radial_gre.py', **protocol: Any) -> ReferenceSequence:
    """Run one of the radial references at `protocol`, returning what it built."""
    module = _load(script)
    arguments = {'plot': False, 'write_seq': False, 'test_report': False, **protocol}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sequence = module.main(**arguments)

    return ReferenceSequence(
        name=f'pypulseq_{Path(script).stem}',
        sequence=sequence,
        parameters=arguments,
        definitions=dict(getattr(sequence, 'definitions', {}) or {}),
        provenance=PROVENANCE | {'path': f'examples/scripts/{script}'},
        # Neither script states an angle, a k-space extent or where its centre sample is.  Every
        # one of those has to be measured, which is the point.
        semantic={'claims': 'none -- angles are computed in the loop, not reported'},
    )
