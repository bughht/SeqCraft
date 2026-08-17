"""
The build notebook's ``GRE2DTR`` and the package's are the same class, asserted.

``examples/gre_2d/01_build.ipynb`` deliberately does **not** import ``GRE2DTR``.  It writes one,
because the third pass of that notebook is where a reader sees a working composition become
something reusable -- and a pass that imported the answer would teach nothing.

The cost of writing it twice is that the two can drift, and the drift is silent: a notebook still
runs, still prints plausible numbers, and still writes a ``.seq``.  So the assertion lives here
rather than in a cell, which keeps the notebook a tutorial and still makes CI notice.

Comparison is by **event digest** rather than by compiled bytes: absolute times and content hashes
of every leaf event, tags ignored.  The two trees are *not* identically shaped -- the notebook's
class nests one level differently in places -- and shape is not what has to match.  What has to
match is every event, at every instant.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import seqcraft as sc

nbformat = pytest.importorskip('nbformat', reason='needs seqcraft[dev]')

NOTEBOOK = Path(__file__).resolve().parents[2] / 'examples' / 'gre_2d' / '01_build.ipynb'


def _cells_up_to_the_class() -> list[str]:
    """
    Return the notebook's code cells up to and including the one defining ``GRE2DTR``.

    Plotting cells are dropped: they draw what the cells around them computed, so skipping them
    changes no value and saves the test a figure it would never look at.
    """
    notebook = nbformat.read(NOTEBOOK, as_version=4)
    sources = [cell.source for cell in notebook.cells if cell.cell_type == 'code']
    end = next(i for i, s in enumerate(sources) if 'class GRE2DTR(sc.Module):' in s)
    return [s for s in sources[: end + 1] if 'plt.' not in s and 'sc.plot_block' not in s]


@pytest.fixture(scope='module')
def notebook(tmp_path_factory) -> dict:
    """
    The notebook's build path, run up to its ``GRE2DTR``, in a scratch directory.

    Running the cells rather than re-deriving their contents is the point: the constants, the raw
    pass, the inline pass and the loop all come from the notebook itself, so a test that passes
    is a statement about the notebook a reader will actually run.

    A scratch cwd because the setup cell creates ``seq/`` for the files the notebook writes, and
    a test run should not leave one behind.
    """
    if not NOTEBOOK.exists():                                       # pragma: no cover
        pytest.skip(f'{NOTEBOOK} is not present')
    import os
    import warnings

    namespace: dict = {'__name__': '__notebook__'}
    here = os.getcwd()
    os.chdir(tmp_path_factory.mktemp('notebook'))
    try:
        with warnings.catch_warnings():
            # The vector-norm warning fires on every compile here and is expected: a blip on y
            # beside a prephaser on x reaches sqrt(2) of the per-axis slew in vector magnitude
            # while every amplifier axis stays inside its own limit.
            warnings.simplefilter('ignore', sc.SeqCraftWarning)
            for index, source in enumerate(_cells_up_to_the_class()):
                exec(compile(source, f'01_build.ipynb:{index}', 'exec'), namespace)  # noqa: S102
    finally:
        os.chdir(here)
    return namespace


def _digest(tree: sc.LogicBlock) -> list[tuple[float, str]]:
    """Every leaf event as (absolute time, content hash), tags and nesting ignored."""
    return sorted(
        (round(t, 12), sc.events.content_hash(event)) for t, event, _ in sc.flatten(tree)
    )


@pytest.fixture(scope='module')
def pair(notebook):
    """The notebook's repetition and the package's, built from the notebook's own constants."""
    from_notebook = notebook['gre_tr']
    from_package = sc.modules.GRE2DTR(
        opts=notebook['opts'],
        fov_mm=notebook['FOV_MM'],
        matrix=notebook['MATRIX'],
        thickness_mm=notebook['THICKNESS_MM'],
        flip_deg=notebook['FLIP_DEG'],
        bandwidth_hz_px=notebook['BANDWIDTH_HZ_PX'],
        tr_s=notebook['TR_S'],
        spoil_cycles_per_voxel=notebook['SPOIL_CYCLES'],
    )
    return from_notebook, from_package


@pytest.mark.parametrize('line', [0, 17, 32, 63])
@pytest.mark.parametrize('acquire', [True, False])
def test_the_two_repetitions_are_the_same_events(pair, line, acquire) -> None:
    """Every line and both build modes, because a drift can be in either branch."""
    from_notebook, from_package = pair

    assert _digest(from_notebook(line=line, phase_deg=117.0, acquire=acquire)) == _digest(
        from_package(line=line, phase_deg=117.0, acquire=acquire)
    )


def test_the_raw_pypulseq_pass_and_the_package_agree(notebook, pair) -> None:
    """
    The end of the chain, stated in one line.

    The notebook's own cells assert that the raw pass, the inline composition and the class it
    writes produce identical events -- and running the fixture runs those asserts.  This closes
    it: what the package builds is what raw ``pp.make_*`` calls build, which is the claim the
    whole module library rests on.
    """
    _, from_package = pair

    assert _digest(notebook['repetition_raw'](17, 117.0)) == _digest(
        from_package(line=17, phase_deg=117.0)
    )


def test_the_two_agree_on_the_numbers_they_report(pair) -> None:
    """TE, TR and the centre line: three answers a caller acts on."""
    from_notebook, from_package = pair

    assert from_notebook.min_te_s == pytest.approx(from_package.min_te_s)
    assert from_notebook.tr_s == pytest.approx(from_package.tr_s)
    assert from_notebook.center_line == from_package.center_line
    assert from_notebook.winder_s == pytest.approx(from_package.winder_s)


def test_a_whole_scan_compiles_to_the_same_thing(notebook, pair) -> None:
    """
    The notebook's ``scan`` loop over the package's repetition, against ``GRE2D``.

    This is the other half of the claim: the notebook writes the loop out so a reader sees the
    RF-spoiling counter running across the dummies, and ``GRE2D`` is that loop and nothing else.
    """
    _, from_package = pair
    lines = tuple(range(0, notebook['MATRIX'][1], 3)) + (from_package.center_line,)
    lines = tuple(sorted(set(lines)))

    by_hand = notebook['scan'](
        lambda line, phase_deg=0.0, acquire=True: from_package(
            line=line, phase_deg=phase_deg, acquire=acquire),
        lines,
    )
    from_module = sc.modules.GRE2D(
        opts=notebook['opts'],
        fov_mm=notebook['FOV_MM'],
        matrix=notebook['MATRIX'],
        thickness_mm=notebook['THICKNESS_MM'],
        flip_deg=notebook['FLIP_DEG'],
        bandwidth_hz_px=notebook['BANDWIDTH_HZ_PX'],
        tr_s=notebook['TR_S'],
        spoil_cycles_per_voxel=notebook['SPOIL_CYCLES'],
        dummies=notebook['DUMMIES'],
        rf_spoil_deg=notebook['RF_SPOIL_DEG'],
    )(lines=lines)

    assert _digest(by_hand) == _digest(from_module)


def test_the_notebook_does_not_import_the_classes_it_teaches(notebook) -> None:
    """
    The rule that makes the third pass worth reading, checked rather than trusted.

    ``sc.modules.Excitation`` and its three siblings are imported freely -- those are the leaves
    the notebook is composing.  ``GRE2DTR`` and ``GRE2D`` are what it writes.
    """
    source = NOTEBOOK.read_text(encoding='utf-8')

    for name in ('sc.modules.GRE2DTR', 'sc.modules.GRE2D('):
        assert name not in source, f'01_build.ipynb reaches for {name}; it is meant to write it'
