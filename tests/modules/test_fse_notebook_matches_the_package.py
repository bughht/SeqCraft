"""
The build notebook's ``FSE2D`` and the package's pair are the same sequence, asserted.

``examples/fse_2d/01_build.ipynb`` writes its own ``FSE2D`` rather than importing one, for the
reason ``test_notebook_matches_the_package.py`` gives about the GRE notebook: a pass that imported
the answer would teach nothing.  What the package then ships is that composition **split in two**
-- :class:`~seqcraft.modules.TSEShot` for the physics of one shot, :class:`~seqcraft.modules.FSE2D`
for the table -- and the split is only trustworthy if it moved no event.

So this compares by event digest at four turbo factors, including the single-echo case that is a
conventional spin echo and the 72-echo single-shot case that is HASTE.  Absolute times and content
hashes of every leaf event, tags and nesting ignored: the trees are *not* identically shaped,
because the package nests a kernel where the notebook nests a method, and shape is not what has to
match.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

import pytest

import seqcraft as sc

nbformat = pytest.importorskip('nbformat', reason='needs seqcraft[dev]')

NOTEBOOK = Path(__file__).resolve().parents[2] / 'examples' / 'fse_2d' / '01_build.ipynb'


@pytest.fixture(scope='module')
def notebook(tmp_path_factory) -> dict:
    """The notebook's cells up to its ``FSE2D``, run in a scratch directory."""
    if not NOTEBOOK.exists():                                       # pragma: no cover
        pytest.skip(f'{NOTEBOOK} is not present')
    cells = [cell.source for cell in nbformat.read(NOTEBOOK, as_version=4).cells
             if cell.cell_type == 'code']
    end = next(i for i, s in enumerate(cells) if 'class FSE2D(sc.Module):' in s)
    sources = [s for s in cells[: end + 1] if 'plt.' not in s and 'sc.plot_block' not in s]

    namespace: dict = {'__name__': '__notebook__'}
    here = os.getcwd()
    os.chdir(tmp_path_factory.mktemp('fse_2d'))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', sc.SeqCraftWarning)
            for index, source in enumerate(sources):
                exec(compile(source, f'{NOTEBOOK.name}:{index}', 'exec'), namespace)  # noqa: S102
    finally:
        os.chdir(here)
    return namespace


def _digest(tree) -> list[tuple[float, str]]:
    """Every leaf event as (absolute time, content hash), tags and nesting ignored."""
    return sorted(
        (round(t, 12), sc.events.content_hash(event)) for t, event, _ in sc.flatten(tree)
    )


def _protocol(notebook: dict, **overrides) -> dict:
    """The notebook's own constants, so neither side restates the protocol."""
    return {
        'opts': notebook['opts'], 'fov_mm': notebook['FOV_MM'], 'matrix': notebook['MATRIX'],
        'thickness_mm': notebook['THICKNESS_MM'], 'tr_s': notebook['TR_S'],
        'bandwidth_hz_px': notebook['BANDWIDTH_HZ_PX'],
        'excitation_duration_s': notebook['EXCITE_DURATION_S'],
        'refocus_duration_s': notebook['REFOCUS_DURATION_S'],
        'refocus_thickness_factor': notebook['REFOCUS_FACTOR'],
        'crush_cycles_slice': notebook['CRUSH_CYCLES_SLICE'],
        'crush_cycles_readout': notebook['CRUSH_CYCLES_READOUT'],
        'spoil_cycles_per_voxel': notebook['SPOIL_CYCLES'],
        'spoil_axis': notebook['SPOIL_AXIS'], **overrides,
    }


def _interleaved(echoes: int, n_lines: int) -> list[list[int]]:
    shots = n_lines // echoes
    return [[s + n * shots for n in range(echoes)] for s in range(shots)]


CASES = [
    pytest.param(1, 0, {}, id='spin echo'),
    pytest.param(8, 1, {}, id='turbo 8'),
    pytest.param(16, 1, {}, id='turbo 16'),
    pytest.param(72, 0, {'partial_fourier': 0.625}, id='haste'),
]


@pytest.mark.parametrize(('echoes', 'dummy_shots', 'extra'), CASES)
def test_the_split_moved_no_event(notebook, echoes, dummy_shots, extra) -> None:
    """Event for event, acquired and dummy, at every turbo factor the notebook builds."""
    settings = _protocol(notebook, echoes=echoes, **extra)
    n_lines = notebook['NY']
    segments = (_interleaved(echoes, n_lines) if echoes != 72
                else [[n_lines // 2 - 8 + n for n in range(72)]])

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        from_notebook = notebook['FSE2D'](**settings)(
            segments=segments, dummy_shots=dummy_shots)
        from_package = sc.modules.FSE2D(**settings)(
            segments=segments, dummy_shots=dummy_shots)

    assert _digest(from_notebook) == _digest(from_package)
    assert from_notebook.duration == pytest.approx(from_package.duration, abs=1e-12)


@pytest.mark.parametrize(('echoes', 'dummy_shots', 'extra'), CASES)
def test_the_two_agree_on_the_numbers_they_report(notebook, echoes, dummy_shots, extra) -> None:
    """Echo spacing, the bound that set it, TR and the effective TE: what a caller acts on."""
    settings = _protocol(notebook, echoes=echoes, **extra)
    n_lines = notebook['NY']
    segments = (_interleaved(echoes, n_lines) if echoes != 72
                else [[n_lines // 2 - 8 + n for n in range(72)]])

    from_notebook = notebook['FSE2D'](**settings)
    from_package = sc.modules.FSE2D(**settings)

    assert from_notebook.echo_spacing_s == pytest.approx(from_package.shot.echo_spacing_s)
    assert from_notebook.bound == from_package.shot.bound
    assert from_notebook.tr_s == pytest.approx(from_package.shot.tr_s)
    assert from_notebook.center_line == from_package.center_line
    assert from_notebook.te_eff_s(segments) == pytest.approx(from_package.te_eff_s(segments))
    assert from_notebook.echo_bands(segments) == from_package.echo_bands(segments)


def test_the_leaf_now_owns_the_readout_geometry_the_notebook_derives(notebook) -> None:
    """
    The one piece of arithmetic the split moved *down* rather than sideways.

    The notebook writes ``(gx.area - 2*area_to_echo)/2`` inline, and so does ``se_2d/01``.  That
    quantity is a fact about the readout's own geometry, so ``CartesianLine`` states it and the
    kernel only decides what to do about it.  Equal to the last bit, because the extraction is
    only provable if nothing moved.
    """
    from_package = sc.modules.FSE2D(**_protocol(notebook, echoes=8))
    readout = from_package.shot.ro

    assert readout.area_to_echo_per_m + readout.area_after_echo_per_m == pytest.approx(
        float(readout.gx.area), abs=1e-12)
    assert readout.echo_moment_imbalance_per_m == (
        float(readout.gx.area) - 2 * readout.area_to_echo_per_m)
    assert from_package.shot.skew_per_m == readout.echo_moment_imbalance_per_m / 2.0
