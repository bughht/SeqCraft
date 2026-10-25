"""
The end-to-end setup script's refusals, which are the only part of it worth automating.

`tools/sequence_authoring/end_to_end/setup.sh` builds a throwaway environment, and an earlier
version began with ``rm -rf "$ROOT"``.  That is a destructive default in a script whose argument is
a path a person types, and the path check that was supposed to keep it away from the planning
repository failed open: it asked git about a path with a missing component, git errored, the
``if`` was false, and the script installed 224 MB into the repository it was meant to protect.

So what is tested here is not that the script builds anything -- that takes minutes and a network
-- but that **it refuses before it touches the disk**.  Every case below asserts both the refusal
and that the files that were already there survived it.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

SETUP = Path(__file__).resolve().parents[2] / 'tools' / 'sequence_authoring' / 'end_to_end' / 'setup.sh'

pytestmark = pytest.mark.skipif(
    sys.platform == 'win32' or shutil.which('bash') is None,
    reason='the setup script is bash; the behaviour it guards is the same everywhere',
)


def run(target: Path, mode: str = 'project') -> subprocess.CompletedProcess:
    """Invoke the script with a short timeout: a refusal is immediate, a build is not."""
    return subprocess.run(
        ['bash', str(SETUP), mode, str(target)],
        capture_output=True, text=True, timeout=60,
    )


def test_the_script_exists_and_is_executable() -> None:
    assert SETUP.is_file()
    assert SETUP.stat().st_mode & 0o111, 'setup.sh is not executable'


# ------------------------------------------------------------------- it refuses, and does no harm
def test_an_existing_directory_is_refused_and_its_contents_survive(tmp_path) -> None:
    """
    The case the destructive default would have destroyed.

    A directory with a user's file in it is exactly what someone types by accident, and the whole
    point of the fix is that this costs them nothing.
    """
    target = tmp_path / 'already_here'
    target.mkdir()
    precious = target / 'notes.txt'
    precious.write_text('work that predates the test run', encoding='utf-8')
    (target / 'subdir').mkdir()

    result = run(target)

    assert result.returncode != 0
    assert 'already exists' in result.stderr
    assert precious.read_text(encoding='utf-8') == 'work that predates the test run'
    assert (target / 'subdir').is_dir()


def test_an_existing_empty_directory_is_refused_too(tmp_path) -> None:
    """Each run gets a new directory, so 'empty' is not a licence to reuse one."""
    target = tmp_path / 'empty'
    target.mkdir()

    assert run(target).returncode != 0
    assert target.is_dir()


def test_a_symlink_is_refused_rather_than_followed(tmp_path) -> None:
    """
    A link is something that already exists, whatever it points at.

    Following one would let the real target be somewhere the path checks never looked -- and a
    dangling link is invisible to `-e`, which is why the script tests `-L` as well.
    """
    real = tmp_path / 'real'
    real.mkdir()
    (real / 'notes.txt').write_text('not to be touched', encoding='utf-8')
    link = tmp_path / 'link'
    link.symlink_to(real, target_is_directory=True)
    dangling = tmp_path / 'dangling'
    dangling.symlink_to(tmp_path / 'nowhere', target_is_directory=True)

    assert run(link).returncode != 0
    assert (real / 'notes.txt').read_text(encoding='utf-8') == 'not to be touched'
    assert run(dangling).returncode != 0
    assert dangling.is_symlink(), 'the dangling link was removed'


def test_a_missing_parent_is_refused_rather_than_created(tmp_path) -> None:
    """`mkdir -p` would invent the whole chain; a typo in a long path should stop here instead."""
    result = run(tmp_path / 'no' / 'such' / 'parent')

    assert result.returncode != 0
    assert 'parent directory does not exist' in result.stderr
    assert not (tmp_path / 'no').exists()


# ------------------------------------------------------- it stays out of the planning repository
def test_a_target_under_the_planning_repository_is_refused(tmp_path) -> None:
    """
    The marker is the criteria directory itself, not the enclosing git repository.

    Asking git which repository encloses a path is fooled by a nested one, which reports its own
    top level and hides the outer tree -- so the check looks for what it actually cares about.
    """
    planning = tmp_path / 'planning'
    (planning / 'docs' / 'plans' / 'sequence-authoring').mkdir(parents=True)
    (planning / 'docs' / 'plans' / 'sequence-authoring' / 'criteria.md').write_text('secret')
    target = planning / 'runs' / 'first'
    target.parent.mkdir(parents=True)

    result = run(target)

    assert result.returncode != 0
    assert 'review criteria' in result.stderr
    assert not target.exists()


def test_a_nested_git_repository_does_not_hide_the_planning_repository(tmp_path) -> None:
    """A repository checked out inside the planning tree must not launder the path check."""
    planning = tmp_path / 'planning'
    (planning / 'docs' / 'plans' / 'sequence-authoring').mkdir(parents=True)
    nested = planning / 'SeqCraft'
    (nested / '.git').mkdir(parents=True)
    target = nested / 'runs'

    result = run(target)

    assert result.returncode != 0
    assert 'review criteria' in result.stderr
    assert not target.exists()


def test_a_symlinked_parent_does_not_escape_the_check(tmp_path) -> None:
    """
    The parent is resolved with `cd -P`, so a link into the planning tree is seen through.

    Without that, `outside/link/runs` would pass a textual check while landing inside the tree the
    script exists to stay out of.
    """
    planning = tmp_path / 'planning'
    (planning / 'docs' / 'plans' / 'sequence-authoring').mkdir(parents=True)
    (planning / 'inner').mkdir()
    outside = tmp_path / 'outside'
    outside.mkdir()
    (outside / 'link').symlink_to(planning / 'inner', target_is_directory=True)

    result = run(outside / 'link' / 'runs')

    assert result.returncode != 0
    assert 'review criteria' in result.stderr
