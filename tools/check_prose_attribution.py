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

    python tools/check_prose_attribution.py                  # blocking, whole tree; exits 1
    python tools/check_prose_attribution.py --review         # + triggers in prose you changed
    python tools/check_prose_attribution.py --review=all     # + every trigger (thousands)
    python tools/check_prose_attribution.py --since=main     # changed prose against a ref
    python tools/check_prose_attribution.py PATH ...         # limit to some paths

**Blocking runs over the whole tree; review runs over what you changed.**  Repository-wide the
review level is some thousands of lines, which nobody reads, and a report nobody reads is the
thing this checker exists to argue against.  Scoped to the lines a branch touched it is a few
hundred, which is a review.

What this does **not** prove
---------------------------
A green run is evidence about the patterns below and nothing else.  In particular it does not
show that the repository contains no unsupported attribution:

- the blocking list is the set of formulations **already observed here**, generalised to the
  families they belong to.  A new way of saying the same thing passes until someone adds it.
- the review level finds words, not claims.  It cannot tell a supported "SeqCraft had a bug here"
  from an unsupported one; only a reader can.
- ``salvage/`` is excluded from the ongoing run.  It was audited once, on 2026-09-28, and is
  frozen; nothing re-checks it.

So this is a **maintained-tree floor**, not an attribution oracle.  The rule it partially enforces
lives in ``docs/writing_a_module.md`` and is the thing to apply; this file only catches the
mistakes that have already been made once.
"""

from __future__ import annotations

import json
import logging
import re
import subprocess
import sys
from collections.abc import Iterator
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: Directories whose prose is not maintained here.  `salvage/` is parked whole and frozen; it was
#: audited once on 2026-09-28 and is excluded from the ongoing run rather than never examined.
SKIP_DIRS = {'.git', '.venv', '__pycache__', 'node_modules', 'salvage', 'seq', '.mypy_cache',
             '.pytest_cache', '.ruff_cache', 'htmlcov', 'build', 'dist'}

SUFFIXES = {'.py', '.md', '.ipynb', '.rst', '.txt', '.toml', '.yaml', '.yml'}

#: The two authoring guides teach the rule by quoting what it rejects, so those quotations have to
#: be allowed somewhere.  They are allowed **inside fenced code blocks, in these two files, and
#: nowhere else** -- a structural rule rather than a list of sentences, so it does not rot the next
#: time a paragraph is re-wrapped.  Ordinary prose in both guides is still checked, which is the
#: property a per-file exemption would have thrown away.
#:
#: The rule is positional, so it exempts *anything* inside a fence in those two files, not only the
#: intended examples.  That is the price of not rotting: the sentence-matching version it replaced
#: had to repeat every example in this file and broke three times in one editing session.  Fenced
#: content in the two guides is reviewed by people rather than by this.
#:
#: The checker and its own tests are exempt wholesale: they quote every pattern by construction.
TEACHING_FILES: dict[str, str] = {
    'tools/check_prose_attribution.py': 'whole file',
    'tests/test_prose_attribution.py': 'whole file',
    'docs/writing_a_module.md': 'fenced blocks only',
    'docs/writing_examples.md': 'fenced blocks only',
}


def exempt(path: Path, in_fence: bool) -> bool:
    """Whether a hit here is a teaching example rather than a claim."""
    try:
        name = path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return False
    rule = TEACHING_FILES.get(name)
    if rule is None:
        return False
    return rule == 'whole file' or in_fence


BLOCKING: tuple[tuple[str, str, str], ...] = (
    # --- an external artifact made the subject of a defect claim -------------------------------
    (r'(reference|public|upstream|example) model was wrong',
     'names an external model as the defective artifact',
     'name the artifact the evidence is about, e.g. "our X did not match Y"'),
    (r'\bthe references\b[^.]{0,40}\b(get|got|have)\b[^.]{0,20}\bwrong\b',
     'attributes a defect to unnamed "references"',
     'name the artifact and version, and state the observed behaviour'),
    (r'every reference implementation\b',
     'a universal claim about implementations that were not all examined',
     'name the ones actually tested, with versions'),

    # --- inferred intent, rather than observed behaviour ---------------------------------------
    # Not every "by accident" -- "a test asserts the gap so the claim cannot become true by
    # accident" is about this repository and is fine.  What is not fine is saying that some other
    # implementation arrives at a *correct* result without meaning to.
    (r'\b(right|correct|wrong)\s+by accident\b',
     'infers intent: it says an implementation is right without meaning to be',
     'state the configuration observed and what it does or does not exercise'),
    (r'avoids? (it|this|them) only by\b',
     'infers intent, and implies the external artifact is otherwise wrong',
     'state the configuration and the consequence, without the motive'),

    # --- one observation widened into a population ---------------------------------------------
    (r'(wrong|right) on (a|every real|every) scanner\b',
     'generalises a local observation to scanners in general',
     'say what was demonstrated, and on what'),
    # The noun phrase is fine -- "this adapter reads a vendor database selected by the caller"
    # says nothing about the category.  What is not fine is a claim about what such a database
    # *can* or *does* hold, which is a claim about a population nobody examined.
    (r'\b(a|an|any|no|every|the) (vendor|scanner) database\s+'
     r'(can|cannot|can\'t|could|has|have|lacks?|knows?|supplies|supply|provides?|carries|'
     r'contains?|returns?)\b',
     'a coverage claim about vendor databases as a category, which was not examined',
     'name the lookup you actually use and say what it returns'),
    # Same shape: "on no spec sheet" is a coverage claim; "no spec sheet was provided" is not.
    (r'\b(on|in) no (spec ?sheet|preset|catalogue|catalog)\b',
     'a coverage claim about sources that were not all examined',
     'name the source you actually consulted'),
    (r'\bno (spec ?sheet|preset|catalogue|catalog)\s+'
     r'(can|could|has|have|lists?|carries|contains?|supplies|supply|provides?)\b',
     'a coverage claim about sources that were not all examined',
     'name the source you actually consulted'),

    # --- asserted behaviour of something not observed here -------------------------------------
    # Present or future tense with an article -- "a scanner refuses this", "the console will
    # reject it" -- claims what scanners do.  A past-tense report of one observed run ("was run,
    # and in that run the block was refused with <exact error>") is evidence, and passes.
    (r'\b(a|an|the) (console|scanner|interpreter) (refuses|rejects|will refuse|will reject|'
     r'silently|mangles?)\b',
     'asserts scanner or console behaviour as a general rule',
     'report one observed run with its exact message, or say where the failure surfaces'),
    (r'silently mangle',
     'asserts console behaviour without naming platform, version or evidence',
     'name the platform and version, or drop the claim'),
    (r'validated against vendor behaviou?r',
     'a validation claim the cited upstream source does not make',
     'say that SeqCraft delegates to that implementation'),

    # --- rhetoric and unscoped conventions ------------------------------------------------------
    (r'gets? blamed on the scanner',
     'rhetoric in place of the measured consequence',
     'describe the artefact and its consequence directly'),
    (r'a scanner may use either',
     'generalises a convention to scanners in general',
     'say that the convention must be established for the pipeline in hand'),

    # --- a historical claim that history has to support ------------------------------------------
    (r'\bdrifted from (it|upstream|the one it named|the model)\b',
     'a historical claim; "drifted" says the two were once equal',
     'if history does not show that, write "did not match" or "differed"'),
)

#: Broad words worth a second look.  Never blocking.
REVIEW = (
    'wrong', 'bug', 'validated', 'correct', 'safe', 'robust', 'conservative',
    'every', 'never', 'only', 'proves', 'establishes', 'reference', 'upstream',
)
REVIEW_RE = re.compile(r'\b(' + '|'.join(REVIEW) + r')\b', re.IGNORECASE)


class Unreadable(Exception):
    """A file this checker was asked to read and could not.

    Raised rather than skipped: for a lint, a notebook that fails to parse must not quietly
    become "contains no attribution issues".
    """

    def __init__(self, path: Path, why: str) -> None:
        super().__init__(f'{path}: {why}')


def prose_of(path: Path) -> Iterator[tuple[int, str, bool]]:
    """Yield ``(line number, text, inside a fenced block)`` for the prose in one file.

    A notebook is JSON, so its prose is inside cell sources; scanning the raw file would both
    miss escaped newlines and match on base64 image payloads.
    """
    try:
        text = path.read_text(encoding='utf-8')
    except (UnicodeDecodeError, OSError) as why:
        raise Unreadable(path, str(why)) from why
    if path.suffix != '.ipynb':
        fenced = False
        for number, line in enumerate(text.splitlines(), start=1):
            if line.lstrip().startswith('```'):
                fenced = not fenced
                continue
            yield number, line, fenced
        return
    try:
        notebook = json.loads(text)
    except json.JSONDecodeError as why:
        raise Unreadable(path, f'not valid notebook JSON -- {why}') from why
    for index, cell in enumerate(notebook.get('cells', [])):
        source = cell.get('source', [])
        # nbformat permits a single string as well as a list of them.  Iterating a string would
        # walk it character by character and silently find nothing.
        if isinstance(source, str):
            source = source.splitlines(keepends=True)
        for offset, line in enumerate(source):
            # Cell and line, which is what a reader can act on -- notebook files have no useful
            # line numbers of their own.
            yield (index + 1) * 100000 + offset + 1, line.rstrip('\n'), False


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
        try:
            lines = list(prose_of(path))
        except Unreadable as why:
            blocking.append(f'{why}\n    ^ unreadable, so it cannot be checked -- failing closed')
            continue
        for index, (line, text, fenced) in enumerate(lines):
            nxt = lines[index + 1][1] if index + 1 < len(lines) else ''
            window = f'{text} {nxt}'
            for pattern, why, fix in checks:
                found = pattern.search(window)
                # Only when it *starts* on this line, so a wrapped phrase is reported once.
                if found and found.start() <= len(text) and not exempt(path, fenced):
                    blocking.append(
                        f'{where(path, line)}\n'
                        f'    {window.strip()[:130]}\n'
                        f'    ^ {found.group(0)!r} -- {why}\n'
                        f'      instead: {fix}')
            hit = REVIEW_RE.search(text)
            if hit:
                review.append(f'{where(path, line)}  [{hit.group(0).lower()}]  {text.strip()[:96]}')
    return blocking, review


#: `@@ -old,n +new,n @@`, whose `+new,n` is the range this branch added or rewrote.
_HUNK = re.compile(r'^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@')


def changed_lines(since: str) -> dict[str, set[int] | None] | None:
    """
    Which lines this branch added or rewrote, per file, or ``None`` if git cannot say.

    Lines rather than files, because a one-word fix in ``CHANGELOG.md`` would otherwise put
    three thousand untouched lines in front of a reviewer -- and a review nobody reads is the
    failure this checker argues against everywhere else.

    A notebook maps to ``None``, meaning "all of it": its diff is over JSON lines, which have no
    relation to the cell-and-line positions reported for its prose.
    """
    try:
        base = subprocess.run(['git', 'merge-base', 'HEAD', since],
                              capture_output=True, text=True, check=True, cwd=ROOT).stdout.strip()
        diff = subprocess.run(['git', 'diff', '--unified=0', base],
                              capture_output=True, text=True, check=True, cwd=ROOT).stdout
    except (subprocess.CalledProcessError, OSError):
        return None

    touched: dict[str, set[int] | None] = {}
    current: str | None = None
    for line in diff.splitlines():
        if line.startswith('+++ b/'):
            name = line[len('+++ b/'):]
            path = ROOT / name
            current = None
            if path.is_file() and path.suffix in SUFFIXES:
                current = str(path)
                touched.setdefault(current, None if path.suffix == '.ipynb' else set())
            continue
        found = _HUNK.match(line) if current else None
        if found and touched.get(current) is not None:
            start, count = int(found.group(1)), int(found.group(2) or 1)
            touched[current].update(range(start, start + count))    # type: ignore[union-attr]
    return touched


def _within(item: str, touched: dict[str, set[int] | None]) -> bool:
    """Whether a review report points at a line this branch actually changed."""
    head = item.split('  ', 1)[0]
    name, _, position = head.rpartition(':')
    path = str(ROOT / name) if name else ''
    if path not in touched:
        return False
    lines = touched[path]
    if lines is None:
        return True                                     # a notebook: report all of its prose
    try:
        return int(position) in lines
    except ValueError:
        return True


def main(argv: list[str]) -> int:
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    show_review = any(a.startswith('--review') for a in argv)
    review_all = '--review=all' in argv
    since = next((a.split('=', 1)[1] for a in argv if a.startswith('--since=')), 'origin/main')
    paths = [a for a in argv if not a.startswith('--')]
    blocking, review = scan(paths)

    scope = 'the whole tree'
    if show_review and not review_all and not paths:
        touched = changed_lines(since)
        if touched is None:
            scope = f'the whole tree ({since} not available)'
        else:
            _, every = scan(list(touched))
            review = [item for item in every if _within(item, touched)]
            scope = f'lines changed against {since}'


    if blocking:
        logging.info('%d blocking attribution issue(s):\n', len(blocking))
        for item in blocking:
            logging.info('%s\n', item)
    else:
        logging.info('no blocking attribution issues')

    if show_review:
        logging.info('\n%d review trigger(s) in %s -- these do NOT fail. Read them and ask who '
                     'the subject is, what evidence exists, and whether the sentence is stronger '
                     'than it:\n', len(review), scope)
        for item in review:
            logging.info('  %s', item)
    else:
        logging.info('%d review trigger(s) tree-wide; --review lists the ones you changed',
                     len(review))

    return 1 if blocking else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
