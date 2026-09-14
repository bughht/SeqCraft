"""
The README's opening figure, regenerated rather than trusted.

[`docs/figures/gre-tr.svg`](../docs/figures/gre-tr.svg) is drawn from the events
:class:`~seqcraft.modules.GRE2DTR` builds and the blocks ``sc.compile`` returns -- not sketched --
so the only way it stays true is if the drawing is *reproduced* here and compared.  That is what
this does: it runs ``tools/draw_sequence_figure.py`` in memory and asserts the committed file is
byte-identical.  Change the winder coupling, the spoiler tail or the boundary rule and the figure
becomes stale in the same commit, with a failure that says which command fixes it.

It is the argument ``tools/check_api_reference.py`` makes for ``docs/api_reference.md``, and the one
``tests/modules/test_notebook_matches_the_package.py`` makes for the tutorials: a picture nobody
executes is wrong within two commits, and the front page is the worst place for one.

The rest of the file asserts the two claims the figure makes *visually*, so that if the drawing is
ever rewritten the claims survive the rewrite: modules overlap on three axes, and one module's
events land in more than one pulseq block.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
FIGURE = ROOT / 'docs' / 'figures' / 'gre-tr.svg'
README = ROOT / 'README.md'
GENERATOR = ROOT / 'tools' / 'draw_sequence_figure.py'


def _generator():
    """Import the drawing script by path; ``tools/`` is not a package."""
    spec = importlib.util.spec_from_file_location('draw_sequence_figure', GENERATOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope='module')
def drawing():
    return _generator()


@pytest.fixture(scope='module')
def repetition(drawing):
    """The repetition the figure is drawn from, and the block edges beneath it."""
    return drawing.build()


def test_the_readme_still_shows_the_figure() -> None:
    """A figure nothing links to rots faster than one nothing checks."""
    assert FIGURE.exists()
    assert 'docs/figures/gre-tr.svg' in README.read_text(encoding='utf-8')


def test_the_committed_figure_is_the_one_the_code_draws(drawing) -> None:
    """Regenerate it and compare.  This is the whole guard; everything below documents a claim."""
    assert drawing.render() == FIGURE.read_text(encoding='utf-8'), (
        'docs/figures/gre-tr.svg is stale -- run: python tools/draw_sequence_figure.py'
    )


def test_modules_overlap_on_three_axes(drawing, repetition) -> None:
    """
    Three boxes across one window, which is the figure's central claim.

    The slice rephaser, the phase-encode blip and the readout prephaser are three modules' events
    on three axes at one instant.  Remove the winder coupling from ``GRE2DTR`` and they serialise,
    the figure loses its overlap, and this fails.
    """
    block, _ = repetition
    starts = {}
    for tag, events in drawing.instances(block):
        for t0, event in events:
            lane = drawing.lane_of(event)
            if lane in ('Gx', 'Gy', 'Gz'):
                starts.setdefault(round(drawing.shape_of(event, t0)[2][0], 9), set()).add((tag, lane))

    together = max(starts.values(), key=len)
    assert len({lane for _, lane in together}) == 3, 'three axes'
    assert len({tag for tag, _ in together}) == 3, 'three modules, none knowing the others'


def test_one_module_lands_in_more_than_one_block(drawing, repetition) -> None:
    """
    ``CartesianLine`` is written as one module and executed across two blocks.

    Nothing in the module says where a boundary falls -- that is the compiler's, and the figure
    shows it as a dashed line straight through the box.
    """
    block, edges = repetition
    spans = {
        tag: [drawing.shape_of(event, t0)[2] for t0, event in events]
        for tag, events in drawing.instances(block)
    }
    start = min(span[0] for span in spans['CartesianLine'])
    end = max(span[1] for span in spans['CartesianLine'])
    crossed = [edge for edge in edges if start < edge - 1e-12 and edge < end - 1e-12]
    assert crossed, 'the readout module should straddle at least one block boundary'


def test_every_module_in_the_tree_has_a_colour(drawing, repetition) -> None:
    """A module drawn in the fallback grey would silently stop being identifiable."""
    block, _ = repetition
    for tag, _ in drawing.instances(block):
        assert tag in drawing.COLOURS, f'{tag} has no colour of its own'
