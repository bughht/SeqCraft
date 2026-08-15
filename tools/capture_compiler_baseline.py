"""Capture reproducible structural and performance baselines for the compiler refactor."""

from __future__ import annotations

import argparse
import collections
import hashlib
import importlib
import importlib.metadata
import json
import logging
import platform
import runpy
import statistics
import subprocess
import tempfile
import time
import tracemalloc
import warnings
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any

import mypy.version
import numpy as np
import pypulseq
import pytest

import seqcraft as sc
from seqcraft.design.timing import EPS, to_ticks

if TYPE_CHECKING:
    from collections.abc import Callable

    from seqcraft.result import CompiledSequence

_ROOT = Path(__file__).resolve().parents[1]
_DEFAULT_OUTPUT = _ROOT / 'tests' / 'baselines' / 'compiler_phase0.json'
_CI_RUN = 'https://github.com/bughht/SeqCraft/actions/runs/31848625433'
_REGISTRATION_FIELDS = frozenset({'id', 'shape_IDs', '_pypulseq_sequence_event_cache'})
_EVENT_FIELDS = ('rf', 'adc', 'gx', 'gy', 'gz')


def _update_hash(digest: Any, value: Any) -> None:
    """Add a deterministic representation of a nested PyPulseq value to `digest`."""
    if isinstance(value, np.ndarray):
        array = np.ascontiguousarray(value, dtype=np.complex128 if np.iscomplexobj(value) else np.float64)
        digest.update(b'array')
        digest.update(repr(array.shape).encode())
        digest.update(array.tobytes())
    elif isinstance(value, np.generic):
        _update_hash(digest, value.item())
    elif isinstance(value, SimpleNamespace):
        digest.update(b'event')
        for key in sorted(vars(value)):
            if key not in _REGISTRATION_FIELDS:
                digest.update(key.encode())
                _update_hash(digest, getattr(value, key))
    elif isinstance(value, dict):
        digest.update(b'dict')
        for key in sorted(value, key=str):
            _update_hash(digest, key)
            _update_hash(digest, value[key])
    elif isinstance(value, (list, tuple)):
        digest.update(type(value).__name__.encode())
        for item in value:
            _update_hash(digest, item)
    else:
        digest.update(type(value).__name__.encode())
        digest.update(repr(value).encode())


def _sha256(value: Any) -> str:
    """Return the canonical digest used by Phase 0 structural snapshots."""
    digest = hashlib.sha256()
    _update_hash(digest, value)
    return digest.hexdigest()


def stable_summary(compiled: CompiledSequence) -> dict[str, Any]:
    """Return compiler output fields used to build the Phase 0 baseline."""
    counts: collections.Counter[str] = collections.Counter()
    content = hashlib.sha256()
    durations: list[int] = []
    for index in sorted(compiled.seq.block_events):
        block = compiled.seq.get_block(index)
        duration = to_ticks(float(compiled.seq.block_durations[index]))
        durations.append(duration)
        _update_hash(content, duration)
        for field in _EVENT_FIELDS:
            event = getattr(block, field, None)
            if event is None:
                continue
            counts[field] += 1
            _update_hash(content, field)
            _update_hash(content, event)
        for label in getattr(block, 'label', None) or ():
            counts[str(label.type)] += 1
            _update_hash(content, label)

    compile_issues = collections.Counter(f'{issue.severity}:{issue.kind}' for issue in compiled.report.issues)
    checked_issues = collections.Counter(
        f'{issue.severity}:{issue.kind}' for issue in compiled.check().issues
    )
    return {
        'n_blocks': compiled.n_blocks,
        'duration_ticks': to_ticks(compiled.duration_s),
        'block_duration_sha256': _sha256(durations),
        'emitted_content_sha256': content.hexdigest(),
        'origins_sha256': _sha256(compiled.origins),
        'event_counts': dict(sorted(counts.items())),
        'compile_issue_counts': dict(sorted(compile_issues.items())),
        'checked_issue_counts': dict(sorted(checked_issues.items())),
        'moments': {
            f'm{order}': {
                axis: float(f'{value:.12g}') for axis, value in sorted(compiled.moments(order).items())
            }
            for order in range(3)
        },
    }


def _builders() -> dict[str, Callable[[], CompiledSequence]]:
    """
    Load the integration recipes without copying them into this tool.

    Two recipes, both raw pypulseq.  The DTI and EPI-DWI entries went with the module library they
    were built on; they return when a module set is written to rebuild them against.
    """
    namespace = runpy.run_path(str(_ROOT / 'tests' / 'integration' / 'conftest.py'))

    return {
        'gre_2d': namespace['build_gre'],
        'se_2d': namespace['build_se'],
    }


def _tracked_build(builder: Callable[[], CompiledSequence]) -> tuple[CompiledSequence, dict[str, int]]:
    """Build one recipe while observing placement and split counts at existing seam points."""
    compiler = importlib.import_module('seqcraft.compiler')
    original_place = compiler._place
    original_axis_gradient = compiler._axis_gradient
    placed: list[Any] = []
    segments: collections.Counter[int] = collections.Counter()

    def track_place(*args: Any, **kwargs: Any) -> list[Any]:
        result = original_place(*args, **kwargs)
        placed.extend(result)
        return result

    def track_axis_gradient(
        axis: str,
        pieces: list[Any],
        start: float,
        end: float,
        *args: Any,
        **kwargs: Any,
    ) -> Any:
        for piece in pieces:
            if piece.start < end - EPS and piece.end > start + EPS:
                segments[id(piece)] += 1
        return original_axis_gradient(axis, pieces, start, end, *args, **kwargs)

    compiler._place = track_place
    compiler._axis_gradient = track_axis_gradient
    try:
        compiled = builder()
    finally:
        compiler._place = original_place
        compiler._axis_gradient = original_axis_gradient

    gradients = [piece for piece in placed if piece.kind in ('grad', 'trap')]
    return compiled, {
        'placed_events': len(placed),
        'input_gradients': len(gradients),
        'gradient_splits': sum(max(0, segments[id(piece)] - 1) for piece in gradients),
    }


def _file_size(compiled: CompiledSequence) -> int:
    """Return the written sequence size without leaving an artifact in the worktree."""
    with tempfile.TemporaryDirectory(prefix='seqcraft-phase0-') as scratch:
        path = Path(scratch) / 'baseline.seq'
        compiled.write(path, sidecar=False)
        return path.stat().st_size


def _environment() -> dict[str, Any]:
    """Describe the interpreter and pinned direct tools used for this capture."""
    commit = subprocess.run(
        ['git', 'rev-parse', 'HEAD'],
        cwd=_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    return {
        'commit': commit,
        'ci_run': _CI_RUN,
        'platform': platform.platform(),
        'machine': platform.machine(),
        'python': platform.python_version(),
        'seqcraft': sc.__version__,
        'pypulseq': pypulseq.__version__,
        'numpy': np.__version__,
        'pytest': pytest.__version__,
        'ruff': importlib.metadata.version('ruff'),
        'mypy': mypy.version.__version__,
    }


def capture(iterations: int) -> dict[str, Any]:
    """Capture deterministic outputs plus local timing and Python-allocation observations."""
    recipes: dict[str, Any] = {}
    for name, builder in _builders().items():
        logging.info('Capturing %s', name)
        samples: list[float] = []
        warning_counts: collections.Counter[str] = collections.Counter()
        stable: dict[str, Any] | None = None
        observed: dict[str, int] | None = None
        compiled: CompiledSequence | None = None
        for _ in range(iterations):
            started = time.perf_counter()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                current, current_observed = _tracked_build(builder)
            samples.append(time.perf_counter() - started)
            warning_counts.update(str(item.message) for item in caught)
            summary = stable_summary(current)
            if stable is not None and summary != stable:
                msg = f'{name} produced different structural output across identical builds'
                raise RuntimeError(msg)
            stable, observed, compiled = summary, current_observed, current

        tracemalloc.start()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            _tracked_build(builder)
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        assert stable is not None and observed is not None and compiled is not None
        emitted_content_sha256 = stable.pop('emitted_content_sha256')
        recipes[name] = {
            'stable': stable,
            'observed': {
                **observed,
                'emitted_content_sha256': emitted_content_sha256,
                'seq_size_bytes': _file_size(compiled),
                'python_warning_counts': dict(sorted(warning_counts.items())),
                'performance': {
                    'iterations': iterations,
                    'wall_s': [round(value, 6) for value in samples],
                    'median_wall_s': round(statistics.median(samples), 6),
                    'min_wall_s': round(min(samples), 6),
                    'max_wall_s': round(max(samples), 6),
                    'peak_python_bytes': peak,
                },
            },
        }
    return {'schema_version': 2, 'environment': _environment(), 'recipes': recipes}


def main() -> None:
    """Write a Phase 0 baseline artifact for review and later differential checks."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--iterations', type=int, default=3)
    parser.add_argument('--output', type=Path, default=_DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.iterations < 1:
        parser.error('--iterations must be positive')
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    payload = capture(args.iterations)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
    logging.info('Wrote %s', args.output)


if __name__ == '__main__':
    main()
