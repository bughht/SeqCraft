#!/usr/bin/env python3
"""
Did the session actually invoke the skill?

Grepping a transcript for the skill's vocabulary cannot answer this: an agent can write about
emitted `.seq` files and validation layers because that is how people talk about sequences, and an
agent can load the skill and then paraphrase none of it. The question needs the session record.

Claude Code writes one JSONL file per session under ``~/.claude/projects/<encoded-cwd>/``, where
the encoding replaces both ``/`` and ``_`` in the working directory with ``-`` -- verified against
the directories that exist, since the underscore rule is not guessable. A real invocation appears
there as an assistant ``tool_use`` block whose ``name`` is ``Skill`` and whose ``input.skill`` is
the skill's name. That shape was confirmed against existing session logs before this was written.

Usage::

    python3 skill_invocation.py ~/skill-tests/epi-01          # the directory the session ran in
    python3 skill_invocation.py ~/skill-tests/epi-01 --skill sequence-authoring

What the answers mean:

``invoked``      the agent called the skill.  Guided-authoring findings are meaningful.
``offered``      the skill was listed as available and never called.  That is a discovery
                 finding, and says nothing about whether the skill improves authoring.
``absent``       the skill was never even listed.  Installation or discovery is broken; look at
                 where ``SKILL.md`` was placed before concluding anything about the skill itself.
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


def scan(log: Path, skill: str) -> tuple[list[dict], bool]:
    """``(invocations, was_offered)`` for one session file."""
    invocations: list[dict] = []
    offered = False
    with log.open(encoding='utf-8', errors='replace') as handle:
        for line in handle:
            if skill not in line and '"Skill"' not in line:
                continue
            if skill in line:
                offered = True                      # listed in the available-skills reminder
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
    return invocations, offered


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

    total, offered_anywhere = 0, False
    for log in files:
        invocations, offered = scan(log, args.skill)
        offered_anywhere |= offered
        total += len(invocations)
        state = f'{len(invocations)} invocation(s)' if invocations else (
            'offered, never invoked' if offered else 'not offered')
        logging.info('%s  %s', log.name[:8], state)
        for inp in invocations:
            logging.info('            input: %s', json.dumps(inp)[:140])

    logging.info('')
    if total:
        logging.info('invoked  -- %d call(s); guided-authoring findings are meaningful', total)
    elif offered_anywhere:
        logging.info('offered  -- the skill was available and never called; this is a DISCOVERY')
        logging.info('            finding, and says nothing about whether the skill improves')
        logging.info('            authoring')
    else:
        logging.info('absent   -- the skill was never listed; check where SKILL.md was placed')
        logging.info('            before concluding anything about the skill itself')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
