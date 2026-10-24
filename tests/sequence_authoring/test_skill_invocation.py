"""
The skill-invocation checker, which exists because a transcript cannot answer its question.

Searching a reply for the skill's vocabulary proves nothing either way: an agent can load the skill
and paraphrase none of it, and an agent can write about emitted ``.seq`` files because that is how
people talk about sequences. Only the session record shows an invocation.

An earlier version of the helper also reported ``offered`` whenever the skill's *name* appeared on
any log line, and ``absent`` when it did not. Neither holds. The name appears when a user types it,
when a path contains it, or when the agent merely discusses it -- in the session that developed this
skill it appears on hundreds of lines without a single invocation -- and absence from the searched
logs is absence of evidence, since a session started elsewhere logs elsewhere. So what is tested
here is that the strong claim rests on the ``tool_use`` record and the weak one is reported as a
hint and counted, not concluded from.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

_TOOL = (Path(__file__).resolve().parents[2]
         / 'tools' / 'sequence_authoring' / 'end_to_end' / 'skill_invocation.py')
_SPEC = importlib.util.spec_from_file_location('skill_invocation', _TOOL)
assert _SPEC is not None and _SPEC.loader is not None
checker = importlib.util.module_from_spec(_SPEC)
sys.modules['skill_invocation'] = checker
_SPEC.loader.exec_module(checker)

SKILL = 'sequence-authoring'


def _line(**message) -> str:
    return json.dumps({'type': 'assistant', 'message': message}) + '\n'


def _invocation(skill: str = SKILL, **extra) -> str:
    return _line(content=[{'type': 'tool_use', 'name': 'Skill', 'input': {'skill': skill, **extra}}])


# -------------------------------------------------------------------- where the logs are looked for
@pytest.mark.skipif(sys.platform == 'win32',
                    reason='the encoding was read off POSIX log directories; see the docstring')
def test_the_log_directory_encoding_converts_underscores_as_well_as_slashes() -> None:
    """
    A path with an underscore lands somewhere a slash-only substitution never looks.

    This is not guessable from the documented behaviour; it was read off the directories that
    exist on this machine. Getting it wrong makes every answer "no log directory", which looks
    like a missing session rather than a bug in the lookup.

    **Skipped on Windows, and that is a real limitation rather than a test convenience.** There
    `Path('/Users/...').resolve()` is drive-anchored with backslashes, so this substitution
    produces something else entirely — and what Claude Code actually names its log directory on
    Windows has not been checked. The helper is verified on POSIX only.
    """
    got = checker.project_log_dir(Path('/Users/someone/Code_mgh/sources/seqcraft_workspace'))

    assert got.name == '-Users-someone-Code-mgh-sources-seqcraft-workspace'
    assert got.parent == Path.home() / '.claude' / 'projects'


# ------------------------------------------------------------------------ the strong claim: invoked
def test_a_tool_use_record_is_counted_as_an_invocation(tmp_path) -> None:
    log = tmp_path / 's.jsonl'
    log.write_text(_invocation(args='some request'), encoding='utf-8')

    invocations, _ = checker.scan(log, SKILL)

    assert len(invocations) == 1
    assert invocations[0]['args'] == 'some request'


def test_a_namespaced_skill_name_still_counts(tmp_path) -> None:
    """Skills may be invoked as ``pack:name``; the record is still this skill being used."""
    log = tmp_path / 's.jsonl'
    log.write_text(_invocation(f'some-pack:{SKILL}'), encoding='utf-8')

    assert len(checker.scan(log, SKILL)[0]) == 1


def test_another_skill_being_invoked_is_not_counted(tmp_path) -> None:
    log = tmp_path / 's.jsonl'
    log.write_text(_invocation('dataviz'), encoding='utf-8')

    assert checker.scan(log, SKILL)[0] == []


# --------------------------------------------------------- the weak hint: a mention is not an offer
def test_mentioning_the_name_is_counted_but_is_not_an_invocation(tmp_path) -> None:
    """
    The correction this file exists for.

    Every line here contains the skill's name for a reason that has nothing to do with it being
    available: the user typed it, a path contains it, the agent discussed it. The old helper
    reported all three as `offered`.
    """
    log = tmp_path / 's.jsonl'
    log.write_text(
        _line(content=[{'type': 'text', 'text': f'let us talk about {SKILL}'}])
        + _line(content=[{'type': 'text', 'text': f'see .claude/skills/{SKILL}/SKILL.md'}])
        + _line(role='user', content=f'use {SKILL}'),
        encoding='utf-8',
    )

    invocations, mentions = checker.scan(log, SKILL)

    assert invocations == [], 'a mention was miscounted as an invocation'
    assert mentions == 3, 'mentions are counted so a human can follow them up'


def test_invocations_and_mentions_are_reported_separately(tmp_path) -> None:
    log = tmp_path / 's.jsonl'
    log.write_text(
        _line(content=[{'type': 'text', 'text': f'about to use {SKILL}'}]) + _invocation(),
        encoding='utf-8',
    )

    invocations, mentions = checker.scan(log, SKILL)

    assert len(invocations) == 1
    assert mentions == 2, 'the invocation line mentions the name too; the counts are independent'


# ------------------------------------------------------------------------------------- robustness
def test_a_malformed_line_does_not_stop_the_scan(tmp_path) -> None:
    """Logs are appended to live, so a truncated final line is ordinary rather than exceptional."""
    log = tmp_path / 's.jsonl'
    log.write_text(f'{{"broken": "{SKILL}"' + '\n' + _invocation(), encoding='utf-8')

    invocations, mentions = checker.scan(log, SKILL)

    assert len(invocations) == 1
    assert mentions == 2


def test_a_tool_use_without_an_input_does_not_raise(tmp_path) -> None:
    log = tmp_path / 's.jsonl'
    log.write_text(_line(content=[{'type': 'tool_use', 'name': 'Skill'}]), encoding='utf-8')

    assert checker.scan(log, SKILL)[0] == []


def test_a_missing_log_directory_is_reported_rather_than_crashing(tmp_path, caplog, monkeypatch) -> None:
    """
    The operator gets an actionable message, not a traceback, when a session logged elsewhere.

    The helper reports through `logging`, matching `tools/check_api_reference.py`, so the message
    is captured by `caplog` rather than by stdout or stderr.
    """
    monkeypatch.setattr(sys, 'argv', ['skill_invocation.py', str(tmp_path / 'never-ran')])

    assert checker.main() == 2
    assert 'no session log directory' in caplog.text
    assert 'expected' in caplog.text, 'the message must say where it looked'


@pytest.mark.parametrize('name', ['sequence-authoring', 'dataviz', 'module-mining'])
def test_the_scan_is_specific_to_the_skill_asked_about(tmp_path, name: str) -> None:
    log = tmp_path / 's.jsonl'
    log.write_text(''.join(_invocation(s) for s in
                           ('sequence-authoring', 'dataviz', 'module-mining')), encoding='utf-8')

    assert len(checker.scan(log, name)[0]) == 1
