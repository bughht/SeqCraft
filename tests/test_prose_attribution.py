"""
The attribution checker, checked.

A prose linter is only worth having if its failures are ones people agree with.  One that fires on
ordinary sentences teaches the next author to work around it -- rename the variable, split the
line, drop the word -- and then it is enforcing evasion rather than attribution.  So the
must-**pass** cases below matter at least as much as the must-fail ones: they are sentences this
repository should be able to write without argument.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

_TOOL = Path(__file__).resolve().parents[1] / 'tools' / 'check_prose_attribution.py'
_SPEC = importlib.util.spec_from_file_location('check_prose_attribution', _TOOL)
assert _SPEC and _SPEC.loader
prose = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(prose)


#: Formulations that must fail: an external artifact carries the defect, or a local observation
#: has been widened into a claim about everything.
MUST_FAIL = [
    'The public PNS reference model was wrong.',
    'the two worth writing down are the two the references get wrong',
    'pypulseq defaults all three to zero, which is wrong on every real scanner',
    "pypulseq's implementation is validated against vendor behaviour",
    "IRPrep keeps 'hypsec' as its default because that is right on a scanner",
    'The pulseq and pypulseq TSE demos get it right by accident.',
    'the same crusher trapezoid twice, which is what every reference implementation writes',
    'both pulseq reference implementations avoid it only by setting their dead time equal',
    'no vendor database can supply them',
    'a sequence built on those is refused or silently mangled at the console',
    'the kind of artefact that gets blamed on the scanner',
    'MRzero reports one sign and a scanner may use either',
    'our copy had drifted from it',
]

#: Sentences the repository should be able to write.  Several deliberately contain the review
#: words -- that is the point: those words are not what makes a sentence wrong.
MUST_PASS = [
    # The same facts, attributed to the artifact the evidence is about.
    "SeqCraft's synthetic_hardware() did not match safe_example_hw().",
    'pypulseq defaults all three to zero, and the PulseqSystems lookup does not carry them.',
    'seqcraft delegates PNS prediction to pypulseq SAFE implementation.',
    'writeTSE.m, as examined, sets dead time and ringdown to the same value; under that '
    'configuration the residual is zero.',
    # Ordinary use of the review words.
    'The compiler is correct here because every block is validated before emission.',
    'This is a bug in seqcraft: the spoiler area was computed from the wrong voxel size.',
    'A conservative derating policy is applied to the readout.',
    'The only axis that matters here is the readout axis.',
    'This test proves that k returns to zero at every echo.',
    'See the API reference for the full signature.',
    'Never mutate a component after it has been built.',
    'Upstream changes propagate automatically.',
    'This establishes the first-moment condition at every acquired echo.',
    # Near-misses that must not fire.
    'The reference frame is the RF effective centre.',
    'Every echo is a k = 0 crossing.',
    'It is wrong to assume the echo sits at the midpoint.',
]


def _scan_text(tmp_path: Path, text: str) -> list[str]:
    target = tmp_path / 'sample.md'
    target.write_text(text, encoding='utf-8')
    blocking, _ = prose.scan([str(target)])
    return blocking


@pytest.mark.parametrize('sentence', MUST_FAIL)
def test_unsupported_attribution_is_blocked(tmp_path, sentence: str) -> None:
    """Each of these is a formulation this repository actually contained, or a near variant."""
    assert _scan_text(tmp_path, sentence), f'not caught: {sentence!r}'


@pytest.mark.parametrize('sentence', MUST_PASS)
def test_ordinary_prose_is_not_blocked(tmp_path, sentence: str) -> None:
    """**False positives are the failure mode that matters**, so these are the real test."""
    found = _scan_text(tmp_path, sentence)
    assert not found, f'false positive on {sentence!r}:\n' + '\n'.join(found)


def test_a_phrase_split_across_a_line_break_is_still_caught(tmp_path) -> None:
    """
    Prose here is hard-wrapped at about a hundred characters.

    A checker that matched single lines would miss most of what it is looking for, and would miss
    it *unpredictably* -- the same sentence caught or not depending on where it happened to wrap.
    """
    wrapped = 'MRzero reports the opposite sign from that convention, and a scanner may\nuse either.'
    assert _scan_text(tmp_path, wrapped)


def test_a_hit_is_reported_once(tmp_path) -> None:
    """The two-line window must not report the same phrase from both windows that contain it."""
    assert len(_scan_text(tmp_path, 'one\nthe references get wrong\nthree')) == 1


def test_notebook_prose_is_read_from_the_cells(tmp_path) -> None:
    """A notebook is JSON; scanning it as text would both miss escaped newlines and match blobs."""
    notebook = {'cells': [{'cell_type': 'markdown', 'metadata': {},
                           'source': ['intro\n', 'which is wrong on every real scanner\n']}],
                'metadata': {}, 'nbformat': 4, 'nbformat_minor': 5}
    target = tmp_path / 'sample.ipynb'
    target.write_text(json.dumps(notebook), encoding='utf-8')

    blocking, _ = prose.scan([str(target)])

    assert len(blocking) == 1
    assert 'cell 1' in blocking[0]


def test_review_triggers_are_reported_but_never_block(tmp_path) -> None:
    """
    The two levels are the whole design.

    A broad word is a reason to look, not a defect; if these blocked, the honest response would be
    to stop using ordinary English rather than to fix an attribution.
    """
    target = tmp_path / 'sample.md'
    target.write_text('Every block is validated, and only the compiler proves it.', encoding='utf-8')

    blocking, review = prose.scan([str(target)])

    assert not blocking
    assert review


def test_the_exemption_list_is_exactly_these_four_files() -> None:
    """
    **The one way this checker can be defeated is by adding files here**, so the list is pinned.

    Each of the four quotes the blocked formulations as its subject: the checker defines them,
    its tests exercise them, and the two authoring guides teach the rule by showing what it
    rejects.  Anything else added here is a decision to stop checking a file, and should be an
    argument in a pull request rather than a quiet edit.
    """
    assert prose.SELF_REFERENTIAL == {
        'check_prose_attribution.py',
        'test_prose_attribution.py',
        'writing_a_module.md',
        'writing_examples.md',
    }


def test_the_repository_has_no_blocking_attribution_issues() -> None:
    """The gate itself, so a new one cannot be added without this failing."""
    blocking, _ = prose.scan([])

    assert not blocking, '\n\n'.join(blocking)
