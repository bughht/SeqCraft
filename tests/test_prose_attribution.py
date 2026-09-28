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
    # --- found by review after the first sweep: the same semantic families, different grammar.
    # Each of these was live on the branch while the checker reported the repository clean.
    'The four site constants are required because no vendor database has them.',
    'It depends on the coil loading, so no vendor database can know it.',
    'These are on no spec sheet and no vendor database has them.',
    'the configuration that makes writeTSE.m correct by accident',
    'writeTSE.m sets both to the same 100 us, so it gets this right by accident',
    'a sequence that compiles and validates cleanly before the console refuses it',
    'it writes a .seq the console refuses an hour later',
    'an unbalanced crusher simulates as perfectly fine and is wrong on a scanner',
    # --- second review pass: the same families, written with an indefinite article ------------
    'required keyword arguments, because a vendor database cannot supply them',
    'A vendor database supplies amplitudes.',
    'It cannot supply an installation dead times, because any vendor database lacks them.',
    'These are on no spec sheet, and the lookup does not return them.',
    'all three are ways a legal-looking tree produces a sequence a scanner refuses',
    'a scanner will reject a block like this',
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
    'A test asserts the gap, so the claim cannot become true by accident.',
    'Two parts cannot be built against different limits by accident.',
    'Every echo is a k = 0 crossing.',
    'It is wrong to assume the echo sits at the midpoint.',
    # An observed event with an exact message is evidence, not a population claim, and must
    # survive: this is the shape `compiler/verification.py` uses for the 67 388-sample readout.
    'Nothing checked them until a 67 388-sample spiral readout reached a scanner, which refused '
    'the block with `fRTEBFinish() failed for block type: ArbX ArbY ADC`.',
    'In the observed run the block was refused with an exact error naming the block type.',
    # The module-mining records define R1/R2/R3 as their reference set, so this is named.
    'The two worth writing down are where the references disagree about the crusher.',
    # --- the noun phrase is not the problem; the coverage claim is ---------------------------
    'This adapter reads a vendor database selected by the caller.',
    'The vendor database URL is configurable.',
    'A scanner database supplied by the site was used for this comparison.',
    'No spec sheet was provided for this installation, so the values were measured.',
    'The preset was taken from a vendor database the site maintains.',
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


def test_the_exemption_list_is_exactly_these_four_paths() -> None:
    """
    **The one way this checker can be defeated is by adding entries here**, so it is pinned.

    Repository-relative paths rather than basenames, or any future file sharing a name would
    inherit the exemption.  And the two guides are exempt only inside fenced blocks, so ordinary
    prose in the documents that define the rule is still held to it.
    """
    assert prose.TEACHING_FILES == {
        'tools/check_prose_attribution.py': 'whole file',
        'tests/test_prose_attribution.py': 'whole file',
        'docs/writing_a_module.md': 'fenced blocks only',
        'docs/writing_examples.md': 'fenced blocks only',
    }


def test_a_guide_is_exempt_inside_a_fence_and_not_outside_one() -> None:
    """
    **Structural rather than a list of sentences**, so re-wrapping a paragraph cannot silently
    widen or break the exemption.

    An earlier attempt matched the surrounding text, which meant each teaching example had to be
    quoted twice -- once in the guide and once in the checker -- and drifted the first time a line
    moved.
    """
    guide = prose.ROOT / 'docs' / 'writing_a_module.md'

    assert prose.exempt(guide, in_fence=True) is True
    assert prose.exempt(guide, in_fence=False) is False
    # The checker's own files quote every pattern by construction, fenced or not.
    assert prose.exempt(prose.ROOT / 'tests' / 'test_prose_attribution.py', in_fence=False) is True
    # Anything else is checked wherever the text sits.
    assert prose.exempt(prose.ROOT / 'README.md', in_fence=True) is False


def test_an_unsupported_claim_in_a_guide_outside_a_fence_still_fails(tmp_path) -> None:
    """The documents that define the rule are not the two documents exempt from it."""
    guide = prose.ROOT / 'docs' / 'writing_a_module.md'
    original = guide.read_text(encoding='utf-8')
    try:
        guide.write_text(original + '\n\nIt is right on a scanner, which settles it.\n',
                         encoding='utf-8')
        blocking, _ = prose.scan([str(guide)])
        assert blocking, 'prose outside a fence in the guide was not checked'
    finally:
        guide.write_text(original, encoding='utf-8')


def test_a_file_that_merely_shares_a_name_is_not_exempt(tmp_path) -> None:
    """A basename list would have let this through; a path list does not."""
    impostor = tmp_path / 'docs' / 'writing_examples.md'
    impostor.parent.mkdir(parents=True)
    impostor.write_text('that is right on a scanner\n', encoding='utf-8')

    blocking, _ = prose.scan([str(impostor)])

    assert blocking, 'a file outside the repository inherited a teaching exemption'


def test_the_repository_has_no_blocking_attribution_issues() -> None:
    """The gate itself, so a new one cannot be added without this failing."""
    blocking, _ = prose.scan([])

    assert not blocking, '\n\n'.join(blocking)


# --------------------------------------------------------------------- reading files safely
def test_a_notebook_that_will_not_parse_is_reported_rather_than_skipped(tmp_path) -> None:
    """
    **Fail closed.**

    Silently skipping an unparseable notebook makes it indistinguishable from a clean one, which
    is precisely the reading a lint must never offer.
    """
    broken = tmp_path / 'broken.ipynb'
    broken.write_text('{"cells": [ this is not json', encoding='utf-8')

    blocking, _ = prose.scan([str(broken)])

    assert len(blocking) == 1
    assert 'broken.ipynb' in blocking[0]
    assert 'unreadable' in blocking[0]


def test_a_cell_whose_source_is_one_string_is_read_as_lines(tmp_path) -> None:
    """
    nbformat permits ``source`` as a string as well as a list of them.

    Iterating a string walks it character by character, so every line would be one character long
    and nothing would ever match -- a false negative that looks exactly like a clean file.
    """
    notebook = {'cells': [{'cell_type': 'markdown', 'metadata': {},
                           'source': 'intro\nthat is right on a scanner\n'}],
                'metadata': {}, 'nbformat': 4, 'nbformat_minor': 5}
    target = tmp_path / 'one_string.ipynb'
    target.write_text(json.dumps(notebook), encoding='utf-8')

    blocking, _ = prose.scan([str(target)])

    assert len(blocking) == 1, 'a string source was not read as lines'


def test_review_can_be_scoped_to_the_lines_a_branch_changed() -> None:
    """
    Review over the whole tree is thousands of lines, and nobody reads thousands of lines.

    Scoping is what makes the second level usable at all, so the helper that does it is asserted
    rather than assumed: it must return a mapping, and a notebook must map to ``None`` meaning
    "all of its prose", because a notebook's diff is over JSON lines that have nothing to do with
    the cell positions its prose is reported at.
    """
    touched = prose.changed_lines('HEAD')

    assert touched is None or isinstance(touched, dict)
    if touched:
        for path, lines in touched.items():
            assert lines is None or isinstance(lines, set)
            if path.endswith('.ipynb'):
                assert lines is None
