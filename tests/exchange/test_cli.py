"""The CLI is a stable machine boundary for MATLAB and other frontends."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import seqcraft.cli as cli

FIXTURE = Path(__file__).parents[1] / 'fixtures' / 'lbtx' / 'minimal.lb.json'


def _run(*args: str):
    env = {**os.environ, 'MPLCONFIGDIR': '/tmp/seqcraft-mpl-cache',
           'NUMBA_CACHE_DIR': '/tmp/seqcraft-numba-cache'}
    return subprocess.run(
        [sys.executable, '-m', 'seqcraft', *args], capture_output=True, text=True, env=env,
        check=False,
    )


def test_validate_json_contract() -> None:
    completed = _run('validate-tree', str(FIXTURE), '--diagnostics', 'json')
    result = json.loads(completed.stdout)

    assert completed.returncode == cli.EXIT_OK
    assert completed.stderr == ''
    assert result == {
        'schema': 'seqcraft.diagnostic', 'version': '0.1', 'command': 'validate-tree',
        'ok': True, 'diagnostics': [], 'input': str(FIXTURE),
    }


def test_compile_json_contract_writes_a_sequence(tmp_path) -> None:
    output = tmp_path / 'minimal.seq'
    completed = _run(
        'compile-tree', str(FIXTURE), '--output', str(output), '--diagnostics', 'json',
    )
    result = json.loads(completed.stdout)

    assert completed.returncode == cli.EXIT_OK
    assert completed.stderr == ''
    assert result['ok'] is True
    assert result['output'] == str(output)
    assert result['duration_s'] == 0.002
    assert result['block_count'] == 2
    assert output.read_text(encoding='utf-8').startswith('# Pulseq sequence file')


def test_invalid_input_has_exit_two_and_json_path(tmp_path) -> None:
    source = tmp_path / 'bad.lb.json'
    source.write_text('{"schema": "wrong"}\n', encoding='utf-8')

    completed = _run('validate-tree', str(source), '--diagnostics', 'json')
    result = json.loads(completed.stdout)

    assert completed.returncode == cli.EXIT_INPUT
    assert result['ok'] is False
    assert result['diagnostics'][0]['category'] == 'exchange'
    assert result['diagnostics'][0]['source_path'] == '$'


def test_compiler_rejection_has_exit_three_and_provenance_path(tmp_path) -> None:
    source = tmp_path / 'off-raster.lb.json'
    document = json.loads(FIXTURE.read_text(encoding='utf-8'))
    document['root']['nodes'][1]['item']['block']['nodes'][0]['start_s'] = '0.000005'
    source.write_text(json.dumps(document), encoding='utf-8')

    completed = _run(
        'compile-tree', str(source), '--output', str(tmp_path / 'never.seq'),
        '--diagnostics', 'json',
    )
    result = json.loads(completed.stdout)

    assert completed.returncode == cli.EXIT_COMPILE
    assert result['ok'] is False
    assert result['diagnostics'][0]['category'] == 'compiler'
    assert result['diagnostics'][0]['source_path'] == 'minimal.readout'


def test_compiler_warning_is_structured_with_category_and_source(tmp_path) -> None:
    source = tmp_path / 'merge.lb.json'
    document = json.loads(FIXTURE.read_text(encoding='utf-8'))
    node = document['root']['nodes'][1]['item']['block']['nodes'][0]
    document['root']['nodes'][1]['item']['block']['nodes'].append(node.copy())
    source.write_text(json.dumps(document), encoding='utf-8')

    completed = _run(
        'compile-tree', str(source), '--output', str(tmp_path / 'merge.seq'),
        '--diagnostics', 'json',
    )
    result = json.loads(completed.stdout)

    assert completed.returncode == cli.EXIT_OK
    assert result['ok'] is True
    assert result['diagnostics'][0]['severity'] == 'warning'
    assert result['diagnostics'][0]['category'] == 'merge'
    assert result['diagnostics'][0]['source_path'] == 'minimal.readout'
