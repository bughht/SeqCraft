"""The process adapter has one job: LBTX in, Pulseq .seq out."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import seqcraft.cli as cli

FIXTURE = Path(__file__).parents[1] / 'fixtures' / 'lbtx' / 'minimal.lbtx.json'


def _run(*args: str):
    env = {
        **os.environ,
        'MPLCONFIGDIR': '/tmp/seqcraft-mpl-cache',
        'NUMBA_CACHE_DIR': '/tmp/seqcraft-numba-cache',
    }
    return subprocess.run(
        [sys.executable, '-m', 'seqcraft', *args], capture_output=True, text=True, env=env,
        check=False,
    )


def test_compile_writes_a_sequence_without_result_protocol(tmp_path) -> None:
    output = tmp_path / 'minimal.seq'
    completed = _run('compile-tree', str(FIXTURE), '--output', str(output))

    assert completed.returncode == 0
    assert completed.stdout == ''
    assert completed.stderr == ''
    assert output.read_text(encoding='utf-8').startswith('# Pulseq sequence file')


def test_invalid_input_has_exit_two_and_plain_error(tmp_path) -> None:
    source = tmp_path / 'bad.lbtx.json'
    source.write_text('{"format": "wrong"}\n', encoding='utf-8')

    completed = _run('compile-tree', str(source), '--output', str(tmp_path / 'never.seq'))

    assert completed.returncode == cli.EXIT_INPUT
    assert completed.stdout == ''
    assert '$:' in completed.stderr


def test_compiler_rejection_has_exit_three(tmp_path) -> None:
    source = tmp_path / 'off-raster.lbtx.json'
    document = json.loads(FIXTURE.read_text(encoding='utf-8'))
    document['root']['nodes'][1]['block']['nodes'][0]['start_s'] = 5e-6
    source.write_text(json.dumps(document), encoding='utf-8')

    completed = _run('compile-tree', str(source), '--output', str(tmp_path / 'never.seq'))

    assert completed.returncode == cli.EXIT_COMPILE
    assert completed.stdout == ''
    assert 'minimal.readout' in completed.stderr
