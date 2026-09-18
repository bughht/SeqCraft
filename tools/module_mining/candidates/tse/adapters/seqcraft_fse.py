"""
R1 -- the SeqCraft notebook-local ``FSE2D``.

The notebook is the implementation under evaluation, so it is **executed, not copied**: its code
cells are run up to the class definition and the class is used from that namespace, exactly as
``tests/modules/test_se_notebooks.py`` already does.  Nothing here reimplements the notebook, and
a change to the notebook reaches this adapter without anyone editing it.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

import seqcraft as sc

from ....reference import ReferenceSequence

if TYPE_CHECKING:
    from collections.abc import Sequence

NOTEBOOK = Path(__file__).resolve().parents[5] / 'examples' / 'fse_2d' / '01_build.ipynb'

PROVENANCE = {
    'repo': 'bughht/SeqCraft',
    'path': 'examples/fse_2d/01_build.ipynb',
    'license': 'MIT',
    'role': 'primary',
}


def notebook_namespace(*, stop_at: str = 'class FSE2D(sc.Module):') -> dict[str, Any]:
    """Run the notebook's code cells up to `stop_at` and return what they defined."""
    import nbformat  # noqa: PLC0415

    cells = [c.source for c in nbformat.read(NOTEBOOK, as_version=4).cells
             if c.cell_type == 'code']
    end = next(i for i, source in enumerate(cells) if stop_at in source)
    namespace: dict[str, Any] = {'__name__': '__module_mining__'}
    here = os.getcwd()
    os.chdir(os.environ.get('TMPDIR', '/tmp'))                       # noqa: S108
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for index, source in enumerate(cells[: end + 1]):
                if 'plt.' in source:
                    continue
                exec(compile(source, f'{NOTEBOOK.name}:{index}', 'exec'), namespace)  # noqa: S102
    finally:
        os.chdir(here)
    return namespace


def interleaved(echoes: int, n_lines: int) -> list[list[int]]:
    """The notebook's own interleaved table: shot `s` takes every `shots`-th line."""
    shots = n_lines // echoes
    return [[s + n * shots for n in range(echoes)] for s in range(shots)]


def build(*, echoes: int = 16, segments: Sequence[Sequence[int]] | None = None,
          dummy_shots: int = 0, namespace: dict[str, Any] | None = None,
          **overrides: Any) -> ReferenceSequence:
    """
    Build and compile one ``FSE2D`` at the notebook's protocol, with `overrides` applied.

    `overrides` reaches ``FSE2D.__init__``; pass ``opts=`` to change the scanner, which is what
    the ringdown experiment does.
    """
    ns = namespace if namespace is not None else notebook_namespace()
    opts = overrides.pop('opts', ns['opts'])
    protocol = {
        'opts': opts, 'fov_mm': ns['FOV_MM'], 'matrix': ns['MATRIX'],
        'thickness_mm': ns['THICKNESS_MM'], 'tr_s': ns['TR_S'],
        'bandwidth_hz_px': ns['BANDWIDTH_HZ_PX'],
        'excitation_duration_s': ns['EXCITE_DURATION_S'],
        'refocus_duration_s': ns['REFOCUS_DURATION_S'],
        'refocus_thickness_factor': ns['REFOCUS_FACTOR'],
        'crush_cycles_slice': ns['CRUSH_CYCLES_SLICE'],
        'crush_cycles_readout': ns['CRUSH_CYCLES_READOUT'],
        'spoil_cycles_per_voxel': ns['SPOIL_CYCLES'], 'spoil_axis': ns['SPOIL_AXIS'],
        'echoes': echoes, **overrides,
    }
    module = ns['FSE2D'](**protocol)
    table = [list(s) for s in (segments if segments is not None
                               else interleaved(echoes, ns['NY']))]

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tree = module(segments=table, dummy_shots=dummy_shots)
        sequence = sc.compile(tree, opts, name='module_mining_fse')

    lines = [line for shot in table for line in shot]
    dky = 1e3 / float(ns['FOV_MM'][1])
    return ReferenceSequence(
        name='seqcraft_fse',
        sequence=sequence,
        tree=tree,
        parameters={k: v for k, v in protocol.items() if k != 'opts'} | {
            'rf_dead_time': float(opts.rf_dead_time),
            'rf_ringdown_time': float(opts.rf_ringdown_time),
            'dummy_shots': dummy_shots,
        },
        definitions=dict(getattr(sequence, 'definitions', {}) or {}),
        provenance=PROVENANCE,
        semantic={
            'echo_spacing_s': float(module.echo_spacing_s),
            'te_first_s': float(module.te_s()),
            'te_eff_s': module.te_eff_s(table),
            'bound': module.bound,
            'crush_s': float(module.crush_s),
            'shift_s': float(module.shift_s),
            'sample_of_echo': int(module.ro.echo_sample(0)),
            'segments': table,
            'expected_ky_per_m': [float((line - module.center_line) * dky) for line in lines],
            'dky_per_m': dky,
            'n_lines': int(ns['NY']),
            # The dummy shots play first and acquire nothing, so the acquired readouts are the
            # table's, in order.
            'dummy_shots': dummy_shots,
        },
    )


def expected_ky(ref: ReferenceSequence) -> list[float]:
    """The ky the table asks for, in acquisition order -- the input to invariant I2."""
    return list(np.asarray(ref.semantic['expected_ky_per_m'], dtype=float))
