#!/usr/bin/env python3
"""
Did the session actually invoke the skill?

Grepping a transcript for the skill's vocabulary cannot answer this: an agent can write about
emitted `.seq` files and validation layers because that is how people talk about sequences, and an
agent can load the skill and then paraphrase none of it. The question needs the session record.

Claude Code writes one JSONL file per session under ``~/.claude/projects/<encoded-cwd>/``, where
the encoding replaces both ``/`` and ``_`` in the working directory with ``-`` -- read off the
directories that exist, since the underscore rule is not guessable.

**That was read on POSIX, and the Windows form has not been checked.** A Windows path is
drive-anchored with backslashes, so this substitution would produce something else; if the lookup
reports no log directory there, suspect the encoding before the session. A real invocation appears
there as an assistant ``tool_use`` block whose ``name`` is ``Skill`` and whose ``input.skill`` is
the skill's name. That shape was confirmed against existing session logs before this was written.

Usage::

    python3 skill_invocation.py ~/skill-tests/epi-01          # the directory the session ran in
    python3 skill_invocation.py ~/skill-tests/epi-01 --skill sequence-authoring

What this can and cannot establish:

``invoked``      positive evidence.  An assistant ``tool_use`` block named ``Skill`` carrying this
                 skill's name was recorded.  Guided-authoring findings are meaningful.

``no invocation  the searched logs contain no such block.  That is the absence of evidence, not
found``          evidence of absence: the logs may not cover the whole session, a session may have
                 run from a different working directory, and a skill's content can reach a model by
                 routes this tool cannot see.  Report it as "no invocation found", and look at
                 where ``SKILL.md`` was placed before concluding anything about the skill itself.

The tool also reports whether the skill's **name appears anywhere** in the logs. That is a weak
hint for a human to follow, **not** proof the skill was offered: the name also appears when a user
types it, when a path contains it, or when the agent merely discusses it. It is printed as a
mention count and nothing is inferred from it.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path


def project_log_dir(workdir: Path) -> Path:
    """
    Claude Code's log directory for a session that ran in `workdir`.

    Both ``/`` and ``_`` become ``-``; a path with an underscore lands somewhere a naive
    slash-only substitution never looks.
    """
    encoded = str(workdir.resolve()).replace('/', '-').replace('_', '-')
    return Path.home() / '.claude' / 'projects' / encoded


def scan(log: Path, skill: str) -> tuple[list[dict], int]:
    """
    ``(invocations, name_mentions)`` for one session file.

    The first is positive evidence. The second is only a count of lines containing the skill's
    name, which is a hint worth following by hand and not a claim about availability.
    """
    invocations: list[dict] = []
    mentions = 0
    with log.open(encoding='utf-8', errors='replace') as handle:
        for line in handle:
            if skill not in line and '"Skill"' not in line:
                continue
            if skill in line:
                mentions += 1
            try:
                record = json.loads(line)
            except ValueError:
                continue
            content = (record.get('message') or {}).get('content')
            for block in content if isinstance(content, list) else []:
                if not isinstance(block, dict):
                    continue
                if block.get('type') == 'tool_use' and block.get('name') == 'Skill':
                    used = (block.get('input') or {}).get('skill')
                    if used == skill or skill in str(used):
                        invocations.append(block.get('input') or {})
    return invocations, mentions


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('workdir', type=Path, help='the directory the tested session ran in')
    parser.add_argument('--skill', default='sequence-authoring')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(message)s', stream=sys.stdout)

    logs = project_log_dir(args.workdir)
    if not logs.is_dir():
        logging.error('no session log directory for %s\n  expected %s', args.workdir, logs)
        return 2

    files = sorted(logs.glob('*.jsonl'))
    if not files:
        logging.error('no session files in %s', logs)
        return 2

    total, mentions_total = 0, 0
    for log in files:
        invocations, mentions = scan(log, args.skill)
        mentions_total += mentions
        total += len(invocations)
        state = (f'{len(invocations)} invocation(s)' if invocations
                 else f'no invocation found ({mentions} line(s) mention the name)')
        logging.info('%s  %s', log.name[:8], state)
        for inp in invocations:
            logging.info('            input: %s', json.dumps(inp)[:140])

    logging.info('')
    if total:
        logging.info('invoked  -- %d call(s) recorded; guided-authoring findings are meaningful',
                     total)
        return 0

    logging.info('no invocation found in %d searched log file(s).', len(files))
    logging.info('  This is the absence of evidence, not evidence of absence: the logs may not')
    logging.info('  cover the whole session, and a session started elsewhere logs elsewhere.')
    if mentions_total:
        logging.info('  The name appears on %d line(s) -- a hint to follow by hand, NOT proof the',
                     mentions_total)
        logging.info('  skill was offered; the name also appears when a user types it or a path')
        logging.info('  contains it.')
    logging.info('  Check where SKILL.md was placed before concluding anything about the skill.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
