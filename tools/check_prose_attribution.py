"""
Attribution and provenance lint for prose: docstrings, docs, notebooks, plans.

**This is not a banned-word list.**  It enforces one rule --

    the grammatical subject of a defect claim must be the artifact for which the evidence exists

-- and it does that at two levels, because most of the words involved have perfectly good uses
inside SeqCraft and a checker that failed on all of them would only teach people to evade it.

``BLOCKING``
    Specific formulations that are unsupported or that attribute a defect to an external project
    without naming what was observed.  These fail the check.  Each one is a phrasing seen in this
    repository, not a guess at what might be bad.

``REVIEW``
    Broad words that are *often* fine and occasionally hide an unsupported claim.  These never
    fail.  They are printed on request so a human or an agent can ask the three questions in
    ``docs/writing_a_module.md``: who is the subject, what evidence exists, and is the sentence
    stronger than that evidence.

Usage::

    python tools/check_prose_attribution.py            # blocking only; exits 1 on a hit
    python tools/check_prose_attribution.py --review   # also list the review triggers
    python tools/check_prose_attribution.py PATH ...   # limit to some paths
"""

from __future__ import annotations

import json
import logging
import re
import sys
from collections.abc import Iterator
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: Directories whose prose is not maintained here.  `salvage/` is parked whole and frozen.
SKIP_DIRS = {'.git', '.venv', '__pycache__', 'node_modules', 'salvage', 'seq', '.mypy_cache',
             '.pytest_cache', '.ruff_cache', 'htmlcov', 'build', 'dist'}

SUFFIXES = {'.py', '.md', '.ipynb', '.rst', '.txt', '.toml', '.yaml', '.yml'}

#: Files whose *subject matter* is these patterns, so they quote them deliberately: the checker,
#: its tests, and the two authoring guides that teach the rule by showing what it rejects.
#:
#: This is the only exemption.  It is a list of **filenames**, not directories, so it cannot
#: quietly grow into somewhere to put prose that would otherwise fail -- adding a file here is a
#: visible edit to the checker, and `test_prose_attribution.py` pins the list.  The cost is real:
#: prose in the two guides is not checked by the tool that the guides describe, so it is checked
#: by whoever reviews them.
SELF_REFERENTIAL = {
    'check_prose_attribution.py',
    'test_prose_attribution.py',
    'writing_a_module.md',
    'writing_examples.md',
}

#: Formulations that fail.  Each is paired with what to write instead -- a checker that only says
#: "no" teaches nothing, and the replacement is the part that is actually hard.
BLOCKING: tuple[tuple[str, str, str], ...] = (
    (r'reference model was wrong',
     'names an external model as the defective artifact',
     'name the artifact the evidence is about, e.g. "our X did not match Y"'),
    (r'\breferences get (it )?wrong\b',
     'attributes a defect to unnamed "references"',
     'name the artifact and version, and state the observed behaviour'),
    (r'wrong on every real scanner',
     'a universal claim about all scanners',
     'say what the integration you use actually supplies'),
    (r'validated against vendor behaviou?r',
     'a validation claim the cited upstream source does not make',
     'say that SeqCraft delegates to that implementation'),
    (r'right on a scanner\b',
     'generalises a local observation to scanners in general',
     'say what was demonstrated, and on what'),
    (r'\b(get|gets|got|getting|got it|gets it)\s+(it\s+)?right by accident\b',
     'infers intent in an external project',
     'state the configuration observed and what it does or does not exercise'),
    (r'every reference implementation\b',
     'a universal claim about implementations that were not all examined',
     'name the ones actually tested, with versions'),
    (r'avoids? (it|this) only by\b',
     'infers intent, and implies the external artifact is otherwise wrong',
     'state the configuration and the consequence, without the motive'),
    (r'no (vendor|scanner) database can\b',
     'a universal claim about databases that were not all examined',
     'say what the lookup you use does not provide'),
    (r'silently mangle',
     'asserts console behaviour without naming platform, version or evidence',
     'name the platform and version, or drop the claim'),
    (r'gets? blamed on the scanner',
     'rhetoric in place of the measured consequence',
     'describe the artefact and its consequence directly'),
    (r'a scanner may use either',
     'generalises a convention to scanners in general',
     'say that the convention must be established for the pipeline in hand'),
    (r'\bdrifted from (it|upstream|the one it named)\b',
     'a historical claim; "drifted" says the two were once equal',
     'if history does not show that, write "did not match" or "differed"'),
)

#: Broad words worth a second look.  Never blocking.
REVIEW = (
    'wrong', 'bug', 'validated', 'correct', 'safe', 'robust', 'conservative',
    'every', 'never', 'only', 'proves', 'establishes', 'reference', 'upstream',
)
REVIEW_RE = re.compile(r'\b(' + '|'.join(REVIEW) + r')\b', re.IGNORECASE)


def prose_of(path: Path) -> Iterator[tuple[int, str]]:
    """Yield ``(line number, text)`` of the prose in one file.

    A notebook is JSON, so its prose is inside cell sources; scanning the raw file would both
    miss escaped newlines and match on base64 image payloads.
    """
    try:
        text = path.read_text(encoding='utf-8')
    except (UnicodeDecodeError, OSError):
        return
    if path.suffix != '.ipynb':
        yield from enumerate(text.splitlines(), start=1)
        return
    try:
        notebook = json.loads(text)
    except json.JSONDecodeError:
        return
    for index, cell in enumerate(notebook.get('cells', [])):
        for offset, line in enumerate(cell.get('source', [])):
            # Cell and line, which is what a reader can act on -- notebook files have no useful
            # line numbers of their own.
            yield (index + 1) * 100000 + offset + 1, line.rstrip('\n')


def files(paths: list[str]) -> Iterator[Path]:
    roots = [Path(p) for p in paths] if paths else [ROOT]
    for root in roots:
        if root.is_file():
            yield root
            continue
        for path in sorted(root.rglob('*')):
            if not path.is_file() or path.suffix not in SUFFIXES:
                continue
            if SKIP_DIRS & set(path.relative_to(ROOT).parts if path.is_relative_to(ROOT)
                               else path.parts):
                continue
            yield path


def where(path: Path, line: int) -> str:
    try:
        name = path.relative_to(ROOT)
    except ValueError:
        name = path
    if path.suffix == '.ipynb' and line > 100000:
        return f'{name}: cell {line // 100000}, line {line % 100000}'
    return f'{name}:{line}'


def scan(paths: list[str]) -> tuple[list[str], list[str]]:
    """Return ``(blocking, review)`` reports.

    Matching runs over a **two-line window**, because prose here is hard-wrapped at about a
    hundred characters and a phrase that falls across a line break is the same phrase.  A hit is
    reported once, at the line it starts on.
    """
    blocking, review = [], []
    checks = [(re.compile(p, re.IGNORECASE), why, fix) for p, why, fix in BLOCKING]
    for path in files(paths):
        if path.name in SELF_REFERENTIAL:
            continue                        # these quote the patterns as their subject matter
        lines = list(prose_of(path))
        for index, (line, text) in enumerate(lines):
            nxt = lines[index + 1][1] if index + 1 < len(lines) else ''
            window = f'{text} {nxt}'
            for pattern, why, fix in checks:
                found = pattern.search(window)
                # Only when it *starts* on this line, so a wrapped phrase is reported once.
                if found and found.start() <= len(text):
                    blocking.append(
                        f'{where(path, line)}\n'
                        f'    {window.strip()[:130]}\n'
                        f'    ^ {found.group(0)!r} -- {why}\n'
                        f'      instead: {fix}')
            hit = REVIEW_RE.search(text)
            if hit:
                review.append(f'{where(path, line)}  [{hit.group(0).lower()}]  {text.strip()[:96]}')
    return blocking, review


def main(argv: list[str]) -> int:
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    show_review = '--review' in argv
    paths = [a for a in argv if not a.startswith('--')]
    blocking, review = scan(paths)

    if blocking:
        logging.info('%d blocking attribution issue(s):\n', len(blocking))
        for item in blocking:
            logging.info('%s\n', item)
    else:
        logging.info('no blocking attribution issues')

    if show_review:
        logging.info('\n%d review trigger(s) -- these do NOT fail. Read them and ask who the '
                     'subject is, what evidence exists, and whether the sentence is stronger '
                     'than it:\n', len(review))
        for item in review:
            logging.info('  %s', item)
    else:
        logging.info('%d review trigger(s); pass --review to list them', len(review))

    return 1 if blocking else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
