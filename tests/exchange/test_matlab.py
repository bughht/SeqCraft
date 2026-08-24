"""Cross-language conformance through the repository MATLAB frontend."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc

ROOT = Path(__file__).parents[2]
MATLAB = Path('/Applications/MATLAB_R2024b.app/bin/matlab')
FRONTEND = ROOT / 'matlab'


def _quote(path: str | Path) -> str:
    return str(path).replace("'", "''")


def _matlab(code: str):
    env = {
        **os.environ,
        'PYTHONPATH': str(ROOT / 'src'),
        'MPLCONFIGDIR': '/tmp/seqcraft-mpl-cache',
        'NUMBA_CACHE_DIR': '/tmp/seqcraft-numba-cache',
    }
    return subprocess.run(
        [str(MATLAB), '-batch', f"addpath('{_quote(FRONTEND)}'); {code}"],
        cwd=ROOT, capture_output=True, text=True, env=env, check=False, timeout=120,
    )


@pytest.mark.crossval
@pytest.mark.skipif(not MATLAB.exists(), reason='MATLAB R2024b is not installed')
def test_python_writer_matlab_reader_writer_python_reader(tmp_path) -> None:
    opts = pp.Opts(
        max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s', B0=3.0,
        rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6, max_b1=1000,
    )
    rf = pp.make_block_pulse(
        flip_angle=0.5, duration=100e-6, delay=opts.rf_dead_time, use='excitation', system=opts,
    )
    gradient = pp.make_arbitrary_grad(
        'x', np.array([0.0, 100.0, 0.0]), first=0.0, last=0.0, system=opts,
    )
    source = tmp_path / 'python-all-events.lb.json'
    output = tmp_path / 'matlab-roundtrip.lb.json'
    root = sc.LogicBlock('cross-language').add([
        [0.0, rf], [0.3e-3, gradient], [0.4e-3, pp.make_delay(0.1e-3)],
    ])
    sc.write_logicblock(source, root, opts, definitions={'FOV': [0.22, 0.22, 0.005]})
    completed = _matlab(
        f"[root, opts, definitions] = seqcraft.readTree('{_quote(source)}'); "
        f"seqcraft.writeTree(root, opts, definitions, '{_quote(output)}');"
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    before = json.loads(source.read_text(encoding='utf-8'))
    after = json.loads(output.read_text(encoding='utf-8'))
    sc.validate_document(after)
    assert sc.semantic_hash(after) == sc.semantic_hash(before)


@pytest.mark.crossval
@pytest.mark.skipif(not MATLAB.exists(), reason='MATLAB R2024b is not installed')
def test_matlab_builder_validates_and_compiles_through_python(tmp_path) -> None:
    tree = tmp_path / 'matlab-built.lb.json'
    sequence = tmp_path / 'matlab-built.seq'
    python = Path(sys.executable)
    code = (
        "opts = seqcraft.scannerOpts(1703040, 6386400000, "
        "RFDeadTimeSeconds=100e-6, RFRingdownTimeSeconds=30e-6, "
        "ADCDeadTimeSeconds=10e-6, B0T=3); "
        "readout = seqcraft.LogicBlock('readout'); "
        "gx = seqcraft.trapezoid('x', 100000, 100e-6, 800e-6, 100e-6); "
        "adc = seqcraft.adc(64, 10e-6, DelaySeconds=100e-6, DeadTimeSeconds=10e-6); "
        "readout = readout.add(0, gx, adc); "
        "root = seqcraft.LogicBlock('matlab-minimal'); root = root.add(0, readout); "
        f"seqcraft.writeTree(root, opts, struct('FOV', [0.22 0.22 0.005]), '{_quote(tree)}'); "
        f"validation = seqcraft.validateTree('{_quote(tree)}', "
        f"PythonExecutable='{_quote(python)}'); assert(validation.ok); "
        f"compiled = seqcraft.compileTree('{_quote(tree)}', '{_quote(sequence)}', "
        f"PythonExecutable='{_quote(python)}'); assert(compiled.ok);"
    )
    completed = _matlab(code)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    root, opts, definitions = sc.read_logicblock(tree)
    assert root.tag == 'matlab-minimal'
    assert definitions['FOV'] == [0.22, 0.22, 0.005]
    assert sc.compile(root, opts, definitions=definitions).duration()[0] == 0.001
    assert sequence.exists()


@pytest.mark.crossval
@pytest.mark.skipif(not MATLAB.exists(), reason='MATLAB R2024b is not installed')
def test_matlab_unit_suite() -> None:
    completed = _matlab("results = runtests('matlab/tests/testLBTX.m'); assert(all([results.Passed]));")

    assert completed.returncode == 0, completed.stdout + completed.stderr
