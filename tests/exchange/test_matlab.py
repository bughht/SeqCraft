"""Cross-language checks for the small MATLAB-to-Python path."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from pypulseq import Sequence

from seqcraft.exchange import read_logicblock

ROOT = Path(__file__).parents[2]
MATLAB = Path('/Applications/MATLAB_R2024b.app/bin/matlab')
FRONTEND = ROOT / 'matlab'
MATLAB_PULSEQ = Path(os.environ.get('PULSEQ_MATLAB_PATH', '/Users/yiyund/Code_mgh/pulseq/matlab'))
MATLAB_EXAMPLES = ROOT / 'examples' / 'matlab'
HAS_MATLAB = MATLAB.exists() and MATLAB_PULSEQ.exists()


def _quote(path: str | Path) -> str:
    return str(path).replace("'", "''")


def _matlab(code: str):
    env = {
        **os.environ,
        'PYTHONPATH': str(ROOT / 'src'),
        'MPLCONFIGDIR': '/tmp/seqcraft-mpl-cache',
        'NUMBA_CACHE_DIR': '/tmp/seqcraft-numba-cache',
    }
    setup = f"addpath('{_quote(MATLAB_PULSEQ)}'); addpath('{_quote(FRONTEND)}');"
    return subprocess.run(
        [str(MATLAB), '-batch', f'{setup} {code}'],
        cwd=ROOT, capture_output=True, text=True, env=env, check=False, timeout=120,
    )


@pytest.mark.crossval
@pytest.mark.skipif(not HAS_MATLAB, reason='MATLAB R2024b or MATLAB Pulseq is not installed')
def test_matlab_native_events_are_read_generically(tmp_path) -> None:
    tree = tmp_path / 'native.lbtx.json'
    code = (
        "opts=mr.opts('MaxGrad',40,'GradUnit','mT/m','MaxSlew',150,'SlewUnit','T/m/s',"
        "'rfDeadTime',100e-6,'rfRingdownTime',30e-6,'adcDeadTime',10e-6,'B0',3); "
        "gx=mr.makeArbitraryGrad('x',[0 1000 0],opts,'first',0,'last',0); "
        "rf=mr.makeArbitraryRf([1+1i 1-1i],pi/1200,opts,'use','excitation'); "
        "root=seqcraft.LogicBlock('native'); root.add(0,gx); root.add(0.2e-3,rf); "
        f"seqcraft.writeLBTX(root,opts,'{_quote(tree)}',Definitions=struct('FOV',[.22 .22 .005]));"
    )
    completed = _matlab(code)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    root, opts, definitions = read_logicblock(tree)
    gradient = root.nodes[0].item
    rf = root.nodes[1].item
    assert gradient.type == 'grad'
    assert np.array_equal(gradient.waveform, np.array([0.0, 1000.0, 0.0]))
    assert np.iscomplexobj(rf.signal)
    assert rf.signal.size > 0
    assert opts.B0 == 3
    assert definitions['FOV'] == [0.22, 0.22, 0.005]


@pytest.mark.crossval
@pytest.mark.skipif(not HAS_MATLAB, reason='MATLAB R2024b or MATLAB Pulseq is not installed')
def test_matlab_compile_returns_official_sequence(tmp_path) -> None:
    output_dir = tmp_path / 'path with spaces'
    output_dir.mkdir()
    sequence = output_dir / 'matlab-built.seq'
    python = Path(sys.executable)
    code = (
        "opts=mr.opts('MaxGrad',40,'GradUnit','mT/m','MaxSlew',150,'SlewUnit','T/m/s',"
        "'adcDeadTime',10e-6,'B0',3); "
        "gx=mr.makeTrapezoid('x',opts,'Area',80,'Duration',1e-3); "
        "adc=mr.makeAdc(64,opts,'Dwell',10e-6,'Delay',gx.riseTime); "
        "root=seqcraft.LogicBlock('matlab-minimal'); root.add(0,gx,adc); "
        f"seq=seqcraft.compile(root,opts,'{_quote(sequence)}',"
        f"PythonExecutable='{_quote(python)}',Definitions=struct('FOV',[.22 .22 .005])); "
        "assert(isa(seq,'mr.Sequence')); [ok,report]=seq.checkTiming; "
        "assert(ok,strjoin(report,newline));"
    )
    completed = _matlab(code)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert sequence.exists()


@pytest.mark.crossval
@pytest.mark.skipif(not HAS_MATLAB, reason='MATLAB R2024b or MATLAB Pulseq is not installed')
def test_matlab_unit_suite() -> None:
    completed = _matlab("results=runtests('matlab/tests'); assert(all([results.Passed]));")

    assert completed.returncode == 0, completed.stdout + completed.stderr


@pytest.mark.crossval
@pytest.mark.skipif(not HAS_MATLAB, reason='MATLAB R2024b or MATLAB Pulseq is not installed')
def test_matlab_gre_examples_have_matching_sequence_semantics(tmp_path) -> None:
    logic_path = tmp_path / 'gre-logicblock.seq'
    module_path = tmp_path / 'gre-module.seq'
    python = Path(sys.executable)
    code = (
        f"addpath('{_quote(MATLAB_EXAMPLES)}'); "
        f"pythonExecutable='{_quote(python)}'; "
        "plotSequence=false; "
        f"outputPath='{_quote(logic_path)}'; "
        "run('examples/matlab/gre_2d_logicblock.m'); logicRoot=root; "
        f"outputPath='{_quote(module_path)}'; "
        "run('examples/matlab/gre_2d_module.m'); "
        "assert(abs(logicRoot.duration-root.duration)<1e-12); "
        "assert(logicRoot.duration==0.16); "
        "assert(numel(logicRoot.nodes)==8 && numel(root.nodes)==8);"
    )
    completed = _matlab(code)

    assert completed.returncode == 0, completed.stdout + completed.stderr

    sequences = []
    for path in (logic_path, module_path):
        sequence = Sequence()
        sequence.read(str(path))
        sequences.append(sequence)

    def semantics(sequence: Sequence):
        blocks = [sequence.get_block(index) for index in range(1, len(sequence.block_events) + 1)]
        return (
            len(blocks),
            sequence.duration()[0],
            sum(block.adc is not None for block in blocks),
            sum(bool(block.label) for block in blocks),
        )

    assert semantics(sequences[0]) == pytest.approx((40, 0.16, 8, 8))
    assert semantics(sequences[1]) == pytest.approx(semantics(sequences[0]))
    for sequence in sequences:
        assert np.allclose(sequence.definitions['FOV'], [0.22, 0.22, 0.005])
        assert np.array_equal(sequence.definitions['Matrix'], [16, 8])
