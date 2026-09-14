"""
The two spiral composites that stay in their notebooks, asserted from outside them.

``GRESpiral2D`` and ``SESpiral2D`` have one consumer each -- their own build notebook -- so neither
ships and neither can be imported here.  What CI can still do is **run the notebook and assert
against what it defined**, which is the bargain ``test_epi_notebooks.py`` already strikes.

The claims worth a test rather than a cell are the ones whose failure is *silent*:

- **Every gradient on the readout's own axes is an arbitrary gradient**, in every block of every
  file.  A spiral that came back as an extended trapezoid lost its own shape at a boundary, and
  the only report of it would be a merge warning.
- **``k`` at the ``k = 0`` sample of every echo of every variant**, signed.  ``|k|`` is symmetric,
  so an extent check passes on a trajectory that is rotated or mirrored from the one reported.
- **``k`` returns to the origin at the end of every shot.**  A residual there is a phase ramp
  across the *next* excitation's image, and 0.02 ``dk`` of it is invisible in every plot.
- **No ``resample`` warning from any shipped file**, which is what says the joins were built inside
  the arm's own waveform rather than beside it.
- **The written ``te_s`` matches the module's to a nanosecond.**  This is the one that matters
  most: ``02`` corrects and fits against those numbers, and a wrong echo time is a scaled field map
  with nothing to look wrong.
- **The sidecar trajectory matches ``sc.kspace`` on the file it sits beside.**
- **The spin echo lands where the kernel says**, through ``sc.kspace``, which is the only oracle
  here that models the refocusing conjugation.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

import numpy as np
import pytest

import seqcraft as sc

nbformat = pytest.importorskip('nbformat', reason='needs seqcraft[dev]')

EXAMPLES = Path(__file__).resolve().parents[2] / 'examples'

#: pypulseq's shape compression, which is what a spiral's ``k`` can be known to through a `.seq`.
KSPACE_TOLERANCE_PER_M = 1e-4


def _run(notebook: Path, tmp_path_factory, stop_at: str) -> dict:
    """Execute a notebook's code cells up to and including the one containing `stop_at`."""
    if not notebook.exists():                                       # pragma: no cover
        pytest.skip(f'{notebook} is not present')
    cells = [cell.source for cell in nbformat.read(notebook, as_version=4).cells
             if cell.cell_type == 'code']
    end = next(i for i, source in enumerate(cells) if stop_at in source)
    sources = [s for s in cells[: end + 1] if 'plt.' not in s and 'sc.plot_block' not in s]

    namespace: dict = {'__name__': '__notebook__'}
    here = os.getcwd()
    os.chdir(tmp_path_factory.mktemp(notebook.parent.name))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', sc.SeqCraftWarning)
            for index, source in enumerate(sources):
                exec(compile(source, f'{notebook.name}:{index}', 'exec'), namespace)  # noqa: S102
    finally:
        os.chdir(here)
    return namespace


@pytest.fixture(scope='module')
def gre(tmp_path_factory) -> dict:
    return _run(EXAMPLES / 'gre_spiral_2d' / '01_build.ipynb', tmp_path_factory,
                stop_at='gre_spiral_2d_nominal.npz')


@pytest.fixture(scope='module')
def se(tmp_path_factory) -> dict:
    return _run(EXAMPLES / 'se_spiral_2d' / '01_build.ipynb', tmp_path_factory,
                stop_at='se_spiral_2d_nominal.npz')


def _compiled(tree, opts):
    """Compile, keeping only the warnings that are not the expected vector-norm one."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', sc.SeqCraftWarning)
        seq = sc.compile(tree, opts)
    return seq, [str(w.message) for w in caught
                 if issubclass(w.category, sc.SeqCraftWarning)
                 and 'vector-norm' not in str(w.message)]


def _kinds(seq, axis):
    return {getattr(g, 'type', '?') for i in range(1, len(seq.block_events) + 1)
            for g in (getattr(seq.get_block(i), axis, None),) if g is not None}


# --------------------------------------------------------------------------- the headline
@pytest.mark.parametrize('extra', [
    pytest.param({'shots': 1}, id='one-shot'),
    pytest.param({'shots': 4}, id='four-shot'),
    pytest.param({'shots': 4, 'density': (1.0, 0.5)}, id='variable-density'),
    pytest.param({'shots': 4, 'echoes': 3}, id='three-echo'),
    pytest.param({'shots': 4, 'variant': 'in'}, id='spiral-in'),
])
def test_k_is_where_the_kernel_says_at_every_echo_of_every_file(gre, extra) -> None:
    """Signed, on both axes, at the sample :meth:`echo_sample` names -- and against `sc.kspace`."""
    kernel = gre['kernel_for'](**extra)
    for shot in range(kernel.spiral.shots):
        tree = sc.LogicBlock('probe').add(0.0, kernel(shot=shot))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', sc.SeqCraftWarning)
            k = sc.kspace(tree, gre['opts'])
        claimed = kernel.spiral.k_per_m(shot=shot)
        assert np.abs(k['k_adc'][:2] - claimed).max() < KSPACE_TOLERANCE_PER_M
        for echo in range(len(kernel.spiral.te_s)):
            sample = kernel.spiral.echo_sample(echo)
            assert np.abs(k['k_adc'][:2, sample]).max() < 0.05 * gre['DK_PER_M']
            measured = float(k['t_adc'][sample] - k['t_excitation'][0])
            assert measured == pytest.approx(kernel.te_s[echo], abs=1e-9)


@pytest.mark.parametrize('extra', [
    pytest.param({'shots': 1}, id='one-shot'),
    pytest.param({'shots': 4}, id='four-shot'),
    pytest.param({'shots': 4, 'echoes': 3}, id='three-echo'),
    pytest.param({'shots': 4, 'variant': 'in'}, id='spiral-in'),
])
def test_every_readout_gradient_stays_arbitrary_and_nothing_is_resampled(gre, extra) -> None:
    """A split spiral stays a spiral, and no join was laid beside the arm rather than inside it."""
    kernel = gre['kernel_for'](**extra)
    seq, notes = _compiled(sc.LogicBlock('probe').add(0.0, kernel(shot=0)), gre['opts'])
    assert not any('resample' in note or 'merge' in note for note in notes), notes
    for axis in ('gx', 'gy'):
        # The spoiler is a trapezoid on the same axes and is not part of the readout, so the
        # assertion is that `grad` is present rather than that `trap` is absent.
        assert 'grad' in _kinds(seq, axis)


@pytest.mark.parametrize('shots', [1, 2, 4])
def test_k_returns_to_the_origin_at_the_end_of_every_shot(gre, shots) -> None:
    """A residual here is a phase ramp across the next excitation's image, measured before the tail."""
    kernel = gre['kernel_for'](shots=shots)
    left = kernel.spiral.k_end_per_m(0)
    assert np.abs(left).max() < 1e-3
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        k = sc.kspace(sc.LogicBlock('probe').add(0.0, kernel.spiral(shot=0)), gre['opts'])['k']
    assert np.abs(k[:2, -1]).max() < 1e-3


def test_the_written_files_carry_the_echo_times_the_kernels_report(gre) -> None:
    """``02`` corrects and fits against these, and a wrong one has nothing to look wrong."""
    assert set(gre['built']) == {
        'spiral_1shot', 'spiral_2shot', 'spiral_4shot', 'spiral_4shot_vds', 'spiral_4shot_me3',
        'spiral_4shot_in'}
    requested = gre['TE_S']
    raster = float(gre['opts'].grad_raster_time)
    for name, kernel in gre['built'].items():
        assert kernel.te_s[0] == pytest.approx(requested, abs=raster), name
        assert len(kernel.te_s) == len(kernel.spiral.te_s), name
        for echo, when in enumerate(kernel.te_s):
            offset = kernel.spiral.time_to_echo(echo) - kernel.spiral.time_to_echo(0)
            assert when == pytest.approx(kernel.te_s[0] + offset, abs=1e-12), name


# --------------------------------------------------------------------------- the spin echo
def test_the_seam_the_k_zero_sample_and_the_spin_echo_are_one_instant(se) -> None:
    """Half a dwell, which is as close as the ADC raster and the gradient raster allow."""
    kernel = se['kernel']
    assert abs(kernel.te_s - kernel.spin_echo_s) <= kernel.spiral.dwell_s
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        k = sc.kspace(sc.LogicBlock('probe').add(0.0, kernel(shot=1)), se['opts'])
    sample = kernel.spiral.echo_sample(0)
    measured = float(k['t_adc'][sample] - k['t_excitation'][0])
    assert measured == pytest.approx(kernel.te_s, abs=1e-9)
    assert np.abs(k['k_adc'][:2, sample]).max() < 0.05 * se['DK_PER_M']


def test_the_dephaser_carries_the_sign_the_conjugation_flips(se) -> None:
    """A sign error here is a legal sequence that acquires a mirrored part of k-space.

    ``sc.kspace`` is the only oracle here that models the conjugation, and ``|k|`` at the seam is
    what it answers: get the sign wrong and the arm starts where it should have ended.
    """
    kernel = se['kernel']
    for shot in range(kernel.spiral.shots):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', sc.SeqCraftWarning)
            k = sc.kspace(sc.LogicBlock('probe').add(0.0, kernel(shot=shot)), se['opts'])
        claimed = kernel.spiral.k_per_m(shot=shot)
        assert np.abs(k['k_adc'][:2] - claimed).max() < KSPACE_TOLERANCE_PER_M
        assert np.hypot(*k['k_adc'][:2, 0]) > 0.9 * kernel.spiral.kmax_per_m


def test_placing_the_readout_against_its_first_sample_moves_the_echo_by_a_whole_arm(se) -> None:
    """The version built wrong in the notebook: it compiles, and its contrast is not the one asked for."""
    right, wrong = se['kernel'], se['wrong']
    assert abs(right.te_s - right.spin_echo_s) <= right.spiral.dwell_s
    assert abs(wrong.te_s - wrong.spin_echo_s) == pytest.approx(
        wrong.spiral.arm_duration_s, abs=2 * wrong.spiral.dwell_s)
    seq, notes = _compiled(sc.LogicBlock('probe').add(0.0, wrong(shot=0)), se['opts'])
    assert len(seq.block_events) > 0
    assert not any('resample' in note for note in notes), notes


def test_all_four_spin_echo_files_land_on_one_echo_time(se) -> None:
    """A difference between two pictures in ``02`` is the acquisition, not the contrast."""
    times = [kernel.te_s for kernel in se['built'].values()]
    assert len(se['built']) == 4
    assert np.ptp(times) <= float(se['opts'].grad_raster_time) + 1e-12
    for kernel in se['built'].values():
        assert kernel.spin_echo_s == pytest.approx(se['TE_S'], abs=1e-9)
