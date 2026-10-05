"""
``examples/bssfp_2d/01_build.ipynb``, asserted from outside it.

The notebook's own prints would catch a sequence that failed to compile.  What they would not
catch is a file that compiles, passes every k-space extent check and reconstructs into a plausible
image while being **unbalanced** -- because an image is not where a missing gradient lobe shows up
first, and because a train of identical repetitions can hide a z lobe carrying the wrong half.

So this file re-measures, on the written ``.seq`` files:

- **Zero net gradient area on every axis over every RF-centre-to-RF-centre interval.**  The
  physical definition, integrated from the compiled waveforms between consecutive RF effective
  centres -- not between block edges, which is a container boundary that merely resembles it.
- **The declared TR is the interval the magnetisation sees**, which holds here because every
  repetition in these files reports the same ``time_to_rf_center()``.
- **The segmented file encodes the same lines as the continuous one.**  Segmentation is a loop in
  the notebook; if it ever became something that changed the encoding, the two tables would part.
- **The ``TE`` and ``TR`` in each file's ``[DEFINITIONS]`` match the module's**, because ``02``
  reads those numbers back and measures against them.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc

nbformat = pytest.importorskip('nbformat', reason='needs seqcraft[dev]')

EXAMPLES = Path(__file__).resolve().parents[2] / 'examples'
NOTEBOOK = EXAMPLES / 'bssfp_2d' / '01_build.ipynb'
FILES = ('bssfp_2d', 'bssfp_2d_segmented')

#: Balance is measured here against **one k-space step**, ``1/FOV``, rather than against zero.
#:
#: The tree balances to about 1e-12 1/m -- ``tests/modules/test_bssfp_2d_tr.py`` measures that,
#: because it measures the design.  This file measures the written ``.seq``, and a Pulseq file
#: stores gradient amplitudes as rounded decimal text, so a round trip through it leaves a few
#: parts in 1e6 of the selection gradient's area behind.  Demanding 1e-12 here would be testing
#: the file format's decimal precision; demanding a fraction of the step that separates two
#: k-space lines is testing something a reconstruction would notice.
STEP_FRACTION = 1e-3


def _run(notebook: Path, tmp_path_factory) -> Path:
    """
    Execute the notebook's code cells into a scratch directory, and return its ``seq/``.

    The files are generated rather than read from the working tree, so this runs in CI on a
    clean checkout and cannot pass by reading an artefact somebody built months ago.  Figures
    are skipped; everything else runs, because the last cell is the one that writes the files.
    """
    if not notebook.exists():                                       # pragma: no cover
        pytest.skip(f'{notebook} is not present')
    sources = [cell.source for cell in nbformat.read(notebook, as_version=4).cells
               if cell.cell_type == 'code' and 'plt.subplots(' not in cell.source]

    namespace: dict = {'__name__': '__notebook__'}
    scratch = tmp_path_factory.mktemp('bssfp_2d')
    here = os.getcwd()
    os.chdir(scratch)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', sc.SeqCraftWarning)
            for index, source in enumerate(sources):
                exec(compile(source, f'{notebook.name}:{index}', 'exec'), namespace)  # noqa: S102
    finally:
        os.chdir(here)
    return scratch / 'seq'


@pytest.fixture(scope='module')
def seq_dir(tmp_path_factory) -> Path:
    return _run(NOTEBOOK, tmp_path_factory)


@pytest.fixture(scope='module')
def nominal(seq_dir):
    return np.load(seq_dir / 'bssfp_2d_nominal.npz')


@pytest.fixture(scope='module')
def written(seq_dir):
    out = {}
    for name in FILES:
        seq = pp.Sequence()
        seq.read(str(seq_dir / f'{name}.seq'))
        out[name] = seq
    return out


def test_the_notebook_writes_both_files(seq_dir) -> None:
    missing = [name for name in FILES if not (seq_dir / f'{name}.seq').exists()]

    assert not missing, missing


def rf_centres(seq) -> np.ndarray:
    """Absolute times of every RF effective centre, seconds."""
    out, t = [], 0.0
    for i in range(1, len(seq.block_events) + 1):
        block = seq.get_block(i)
        if getattr(block, 'rf', None) is not None:
            out.append(t + float(block.rf.delay) + float(pp.calc_rf_center(block.rf)[0]))
        t += float(seq.block_durations[i])
    return np.array(out)


def rf_to_rf_m0(seq) -> np.ndarray:
    """``M0`` over every RF-centre-to-RF-centre interval: shape ``(N, 3)``, 1/m."""
    waveforms = seq.waveforms_and_times()[0]
    centres = rf_centres(seq)
    rows = []
    for n in range(len(centres) - 1):
        a, b = centres[n], centres[n + 1]
        row = []
        for axis in range(3):
            t, g = waveforms[axis][0], waveforms[axis][1]
            if len(t) == 0:
                row.append(0.0)
                continue
            grid = np.unique(np.concatenate([t[(t > a) & (t < b)], [a, b]]))
            row.append(float(np.trapezoid(np.interp(grid, t, g, left=0.0, right=0.0), grid)))
        rows.append(row)
    return np.array(rows)


@pytest.mark.parametrize('name', FILES)
def test_every_rf_to_rf_interval_is_balanced(written, nominal, name) -> None:
    """The defining condition, on both written files, measured where it is defined."""
    m0 = rf_to_rf_m0(written[name])
    step_per_m = 1.0 / (float(nominal['fov_mm']) * 1e-3)

    assert len(m0) > 50          # the files are whole acquisitions, not a probe
    assert np.abs(m0).max() < STEP_FRACTION * step_per_m


def test_the_continuous_file_runs_at_the_declared_tr(written, nominal) -> None:
    """
    Every repetition has the same RF-centre offset, so stacking by block duration delivers the
    declared TR.  A notebook that mixed in a repetition of different pulse geometry would need
    the offsets instead, and this is where that would surface.
    """
    intervals = np.diff(rf_centres(written['bssfp_2d']))

    assert intervals == pytest.approx(float(nominal['tr_s']))


def test_the_segmented_file_is_the_interrupted_policy(written, nominal) -> None:
    """
    Interrupted, not maintained: the train stops between segments and the magnetisation recovers.

    Two policies go by "segmented" and they are different sequences.  This asserts which one the
    file is -- every interval is either TR or TR plus the recovery delay, and there are exactly
    as many long ones as there are segment boundaries.  A maintained acquisition would have no
    long intervals at all and would gate the ADC instead.
    """
    tr_s, recovery_s = float(nominal['tr_s']), float(nominal['recovery_s'])
    segments = int(nominal['segments'])
    intervals = np.diff(rf_centres(written['bssfp_2d_segmented']))
    long = np.abs(intervals - (tr_s + recovery_s)) < 1e-9
    short = np.abs(intervals - tr_s) < 1e-9

    assert np.all(long | short)
    assert long.sum() == segments - 1


@pytest.mark.parametrize('name', FILES)
def test_the_definitions_match_the_module(written, nominal, name) -> None:
    """``02`` reads TE and TR back out of the file and measures against them."""
    definitions = written[name].definitions

    assert float(definitions['TE']) == pytest.approx(float(nominal['te_s']), abs=1e-9)
    assert float(definitions['TR']) == pytest.approx(float(nominal['tr_s']), abs=1e-9)


def test_segmentation_changes_the_order_and_not_the_encoding(written, nominal) -> None:
    """
    The whole point of segmentation being a loop: it adds start-up repetitions and reorders
    nothing about what each repetition encodes.
    """
    tables = {}
    for name in FILES:
        labels = written[name].evaluate_labels(evolution='adc')
        tables[name] = np.sort(np.atleast_1d(np.asarray(labels['LIN'])))

    assert np.array_equal(tables['bssfp_2d'], tables['bssfp_2d_segmented'])
    assert len(tables['bssfp_2d']) == len(nominal['lines'])


def test_the_notebook_protocol_still_builds_the_same_repetition(nominal) -> None:
    """
    The module and the notebook agree on the protocol.

    Built here from the parameters the notebook saved, so a change to either side that moved TE or
    TR without the other noticing fails before the simulation notebook silently fits the wrong
    echo times.
    """
    opts = pp.Opts(
        max_grad=24, grad_unit='mT/m', max_slew=120, slew_unit='T/m/s', B0=3.0,
        rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
    )
    tr = sc.modules.bSSFP2DTR(
        opts=opts, fov_mm=float(nominal['fov_mm']),
        matrix=tuple(int(v) for v in nominal['matrix']),
        thickness_mm=float(nominal['thickness_mm']),
        flip_deg=float(nominal['flip_deg']),
        bandwidth_hz_px=float(nominal['bandwidth_hz_px']),
    )

    assert tr.te_s == pytest.approx(float(nominal['te_s']), abs=1e-12)
    assert tr.tr_s == pytest.approx(float(nominal['tr_s']), abs=1e-12)
    assert tr.time_to_rf_center() == pytest.approx(
        float(nominal['time_to_rf_center_s']), abs=1e-12)
    assert tr.symmetry_residual_s == pytest.approx(
        float(nominal['symmetry_residual_s']), abs=1e-12)
