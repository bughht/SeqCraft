"""
The notebook smoke runner's scheduling, with no kernels involved.

The runner is a CI gate, so the properties worth pinning are the ones a reader of a red build
depends on: every allowlisted notebook ran, the report says which ones failed, and it says so in
the same order every time.  Executing real notebooks is the `examples` job's work and is not
repeated here -- these tests stub execution and assert the scheduling around it.
"""

from __future__ import annotations

import importlib.util
import threading
from pathlib import Path

import pytest

_TOOL = Path(__file__).resolve().parents[1] / 'tools' / 'run_notebook_smoke.py'
_SPEC = importlib.util.spec_from_file_location('run_notebook_smoke', _TOOL)
assert _SPEC is not None and _SPEC.loader is not None
smoke = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(smoke)


class _CellFailed(RuntimeError):
    """What a failing notebook cell looks like to the runner: one exception carrying its name."""

    def __init__(self, notebook: Path) -> None:
        super().__init__(f'{notebook}: cell 3 raised')


@pytest.fixture
def notebooks() -> list[Path]:
    return [Path(f'dir_{i}/01_build.ipynb') for i in range(9)]


# ------------------------------------------------------------------ every notebook, exactly once
def test_every_notebook_is_executed_exactly_once(notebooks) -> None:
    seen: list[Path] = []
    lock = threading.Lock()

    def execute(relative: Path) -> None:
        with lock:
            seen.append(relative)

    assert smoke.run(notebooks, execute, max_workers=4) == []
    assert sorted(seen) == sorted(notebooks)
    assert len(seen) == len(notebooks), 'a notebook was submitted more than once'


def test_the_shipped_allowlist_is_what_runs() -> None:
    """
    ``_NOTEBOOKS`` is the source of truth, including the ``03`` notebooks in it.

    The count is deliberately not pinned: adding an example is ordinary work, and a test that
    fails on every such addition teaches people to edit the number rather than read the test.
    What does not change is that each entry is distinct and each one exists -- a duplicate would
    execute a notebook twice, and a stale path would silently stop covering one.

    Deliberately **not** asserted: that every file named ``01_build.ipynb`` appears here.  That
    would make the filename a second source of truth beside the allowlist, and it would also
    decide, by accident, that any future build notebook must be independent enough to run
    unordered -- a policy this repository has not adopted.  ``docs/testing.md`` leaves room for a
    dependent notebook to need ordering or another Layer-2 path instead.
    """
    assert len({*smoke._NOTEBOOKS}) == len(smoke._NOTEBOOKS), 'the allowlist repeats a notebook'
    examples = Path(__file__).resolve().parents[1] / 'examples'
    missing = [n for n in smoke._NOTEBOOKS if not (examples / n).exists()]
    assert missing == [], f'allowlisted notebooks that do not exist: {missing}'


# ------------------------------------------------------------------------------ the worker bound
@pytest.mark.parametrize(('cpus', 'expected'), [(1, 1), (2, 2), (4, 4), (8, 4), (64, 4)])
def test_the_worker_count_is_bounded(notebooks, cpus: int, expected: int) -> None:
    """Four is a ceiling, not a target: more CPUs must not buy more simultaneous kernels."""
    assert smoke._worker_count(notebooks, cpus) == expected
    assert smoke._worker_count(notebooks, cpus) <= smoke._MAX_WORKERS == 4


def test_a_short_list_does_not_start_idle_workers() -> None:
    assert smoke._worker_count([Path('only/01_build.ipynb')], 64) == 1


def test_no_more_than_the_bound_run_at_once(notebooks) -> None:
    """The bound is on *simultaneous* execution, which is what a runner's memory cares about."""
    lock = threading.Lock()
    live = 0
    high_water = 0

    def execute(_relative: Path) -> None:
        nonlocal live, high_water
        with lock:
            live += 1
            high_water = max(high_water, live)
        try:
            threading.Event().wait(0.02)
        finally:
            with lock:
                live -= 1

    smoke.run(notebooks, execute, max_workers=3)
    assert high_water <= 3


# ----------------------------------------------------------------------------------- failures
def test_one_failure_is_reported(notebooks) -> None:
    broken = notebooks[4]

    def execute(relative: Path) -> None:
        if relative == broken:
            raise _CellFailed(relative)

    failures = smoke.run(notebooks, execute, max_workers=4)

    assert [n for n, _ in failures] == [broken]
    assert 'cell 3 raised' in str(failures[0][1])


def test_a_failure_does_not_stop_the_others(notebooks) -> None:
    """
    The serial runner stopped at the first error.  With several kernels already running, where
    it stopped would depend on timing, so every scheduled notebook is allowed to finish.
    """
    ran: list[Path] = []
    lock = threading.Lock()

    def execute(relative: Path) -> None:
        with lock:
            ran.append(relative)
        if relative == notebooks[0]:
            raise _CellFailed(relative)

    failures = smoke.run(notebooks, execute, max_workers=4)

    assert [n for n, _ in failures] == [notebooks[0]]
    assert sorted(ran) == sorted(notebooks), 'a failure cancelled work that was already scheduled'


def test_failures_are_reported_in_allowlist_order_not_completion_order(notebooks) -> None:
    """
    The property the report exists for.

    The notebooks are made to fail in reverse order of their position, with the later ones
    returning first, so a report built from completion order would come out backwards.
    """
    doomed = [notebooks[1], notebooks[3], notebooks[7]]
    delay = {notebooks[1]: 0.06, notebooks[3]: 0.03, notebooks[7]: 0.0}

    def execute(relative: Path) -> None:
        if relative in doomed:
            threading.Event().wait(delay[relative])
            raise _CellFailed(relative)

    for _ in range(3):
        failures = smoke.run(notebooks, execute, max_workers=4)
        assert [n for n, _ in failures] == doomed


def test_a_clean_run_reports_nothing(notebooks) -> None:
    assert smoke.run(notebooks, lambda _relative: None, max_workers=4) == []
