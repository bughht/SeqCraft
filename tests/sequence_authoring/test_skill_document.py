"""
The skill document's load-bearing instructions, asserted so an edit cannot quietly drop one.

**These are document checks, not behaviour checks.** No test here establishes that an agent honours
an instruction; only a fresh-session end-to-end run with a human reviewer can show that. What they
do establish is that the instructions a reviewer agreed to are still present and still paired with
their obligations — which is the failure mode a prose edit actually has.

The pairs matter more than the sentences. "You may skip validation when asked" without "and label
it unvalidated" is worse than neither, so they are checked together.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

SKILL = (Path(__file__).resolve().parents[2]
         / '.claude' / 'skills' / 'sequence-authoring' / 'SKILL.md')


@pytest.fixture(scope='module')
def text() -> str:
    return SKILL.read_text(encoding='utf-8')


@pytest.fixture(scope='module')
def body(text: str) -> str:
    """The instructions, without the YAML frontmatter that only describes when to load them."""
    return re.sub(r'^---\n.*?\n---\n', '', text, count=1, flags=re.S)


@pytest.fixture(scope='module')
def flat(body: str) -> str:
    """
    The same prose with its hard wrapping collapsed.

    The document is wrapped at 110 characters, so a phrase this file cares about is as likely to
    straddle a newline as not -- matching the raw text makes a passing assertion depend on where
    the paragraph happened to break.
    """
    return ' '.join(body.split())


def sentences(flat_text: str) -> list[str]:
    return [s.strip() for s in re.split(r'(?<=[.!?])\s+', flat_text) if s.strip()]


def test_the_frontmatter_names_the_skill_and_describes_when_to_use_it(text: str) -> None:
    front = re.match(r'^---\n(.*?)\n---\n', text, flags=re.S)
    assert front, 'a skill needs YAML frontmatter to be discoverable'
    assert re.search(r'^name:\s*sequence-authoring\s*$', front.group(1), flags=re.M)
    assert re.search(r'^description:\s*\S', front.group(1), flags=re.M)


# --------------------------------------------------- validation is on by default, and waivable
def test_validation_is_stated_as_on_by_default(flat: str) -> None:
    assert 'on by default' in flat.lower()


def test_the_opt_out_comes_with_its_labelling_obligation(flat: str) -> None:
    """
    The pair, not the permission.

    An instruction allowing an unvalidated draft, without the instruction to say so, would turn a
    helpful shortcut into a misleading deliverable.
    """
    assert 'unvalidated' in flat
    assert re.search(r'mark it \*\*unvalidated\*\*', flat), 'the labelling instruction is missing'
    assert re.search(r'list the checks you did not run|checks that would normally run', flat), (
        'an unvalidated hand-off must still say what was not done'
    )


def test_unvalidated_means_independent_validation_was_not_performed(flat: str) -> None:
    """
    An opt-out declines independent validation; it does not switch off the built-in checks.

    Building, compiling and writing a file run checks of their own, and an unvalidated draft still
    compiles. "None were run" would then be untrue. Four things have to stay together: what the
    word means, that built-in checks are never bypassed, that they are reported truthfully, and
    that passing them is neither validation nor scanner readiness.
    """
    assert 'none were run' not in flat, 'an absolute claim that no checks ran'
    assert 'Unvalidated means independent validation was not performed' in flat
    assert re.search(r'Never bypass those', flat), 'built-in checks must not be bypassed'
    assert re.search(r'report truthfully that they ran', flat)
    assert re.search(r'Passing them is not independent validation', flat)


def test_an_unperformed_check_may_never_be_implied_to_have_passed(flat: str) -> None:
    assert 'Never imply an unperformed check passed' in flat


def test_compiling_is_denied_as_validation(flat: str) -> None:
    assert '"compiles" is never validation' in flat


def test_scanner_readiness_is_denied_rather_than_qualified(flat: str) -> None:
    """
    Every mention must sit in a denial; none may read as a condition under which it is claimed.

    Checked per sentence rather than per line: the denial and the phrase are routinely on opposite
    sides of a line break, and a per-line check fails on formatting instead of on meaning.
    """
    mentions = [s for s in sentences(flat) if 'scanner-ready' in s]

    assert mentions, 'the scanner-readiness denial has disappeared entirely'
    for sentence in mentions:
        assert re.search(r'\bnot\b|\bnever\b|\bnothing\b|must not', sentence, flags=re.I), (
            f'a scanner-readiness mention that is not a denial: {sentence!r}'
        )


# ----------------------------------------------- scope follows the request, and stays complete
def test_the_default_scope_is_the_smallest_complete_acquisition(flat: str) -> None:
    assert 'smallest complete' in flat
    assert 'not an acquisition' in flat, 'the "readout alone is not an acquisition" caveat is gone'


def test_optional_features_are_required_to_have_a_reason(flat: str) -> None:
    assert re.search(r'[Oo]ptional features need a reason', flat)
    assert 'fat saturation' in flat.lower(), 'the concrete cautionary example is gone'


def test_stopping_guidance_does_not_license_skipping_physics(body: str, flat: str) -> None:
    """Effort guidance that could be read as "cut checks to finish" would undo the rigour."""
    assert 'When to stop adding evidence' in body
    assert re.search(r'licenses skipping\s+a physical check', flat)
    assert re.search(r'rigou?r intact', flat)


# ------------------------------------------- saved simulation output, reused only while it applies
SIMULATION_INPUTS = (
    'the file actually simulated', 'surrogate or extraction', 'the phantom', 'tissue parameters',
    'initial magnetisation', 'the model assumptions', 'the simulation code and its version',
    'the numerical controls', 'discretisation', 'random seed',
)


def test_saved_simulation_output_is_reused_only_for_the_same_configuration(flat: str) -> None:
    """
    Reuse is the saving; the manifest of what produced the output is the safeguard.

    A `.seq` hash says nothing about the phantom, the model or the numerical controls, so reuse is
    keyed to every input of the simulation, and the record stays a manifest rather than a framework.
    """
    assert 'Reuse saved results for an unchanged simulation configuration' in flat
    assert 'When only the analysis, reconstruction or plotting changes' in flat
    assert re.search(r'A hash of the `\.seq` alone does not identify a simulation', flat)
    for item in SIMULATION_INPUTS:
        assert item in flat, f'the manifest no longer names {item!r}'
    assert 'not a caching framework' in flat


def test_a_changed_simulation_input_reopens_reuse_even_when_the_seq_is_unchanged(flat: str) -> None:
    for change in ('a pruning threshold', 'a state budget', 'an RF approximation', 'a surrogate'):
        assert change in flat, f'{change!r} is no longer named as a change that reopens reuse'
    assert re.search(r'even when the `\.seq` is unchanged', flat)
    assert 're-simulate unless you can show it does not depend on what changed' in flat


def test_reusing_a_result_replaces_no_convergence_check_and_validates_nothing(flat: str) -> None:
    """A saved result inherits exactly the evidence of the run that made it, and no more."""
    assert 'Caching does not replace convergence checks, independent controls or necessary repeats' \
        in flat
    assert 'no more validated than the run that produced it' in flat
    assert 'A number that moves with a numerical control is not a result' in flat, (
        'the convergence rule that reuse must not replace is gone'
    )


# ---------------------------------------------------------------- the boundaries, still stated
def test_the_opening_notice_is_marked_non_blocking(flat: str) -> None:
    assert 'Do not wait for an answer' in flat
    assert 'questionnaire' in flat


def test_no_mode_api_was_introduced(body: str) -> None:
    """
    v0 deliberately has no mode taxonomy; natural language controls both dimensions.

    A flag, a named tier or a preset registry would be the thing the maintainer declined.
    """
    for forbidden in ('--mode', 'mode=', 'ultra-light', 'MODE:'):
        assert forbidden not in body, f'a mode API crept in: {forbidden!r}'


def test_the_package_is_still_treated_as_immutable(flat: str) -> None:
    assert re.search(r'[Dd]o not edit, patch, or monkey-patch', flat)
    assert 'module-mining' in flat, 'the maintainer-workflow boundary is gone'


def test_module_reuse_is_still_the_first_preference(flat: str) -> None:
    order = [flat.index(s) for s in ('a shipped SeqCraft Module',
                                     'a sequence-local Module or function',
                                     'direct LogicBlock composition')]
    assert order == sorted(order), 'the implementation preference order has been disturbed'
