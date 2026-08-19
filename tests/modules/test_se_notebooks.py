"""
The two spin-echo composites that stay in their notebooks, asserted from outside them.

``SE2D`` and ``FSE2D`` have one consumer each -- their own build notebook -- so neither ships and
neither can be imported here.  What CI can still do is **run the notebook and assert against what
it defined**, which is the bargain ``test_mprage_notebooks.py`` already strikes: the tutorial stays
a tutorial, and it cannot drift without something noticing.

The claims worth a test rather than a cell are the ones whose failure is *silent*:

- **k at the echo sample, signed, on all three axes.**  A wrong sign or an unbalanced crusher is a
  legal block structure that meets every limit, and against a symmetric phantom it produces a
  plausible image.  ``|k|`` is symmetric, so a k-space *extent* check passes on the mirrored one.
- **Uniform pulse spacing, and the echo at the midpoint.**  An FSE signal is a mixture of primary
  and stimulated echoes, and they coincide only if both hold.  Nothing downstream complains if they
  do not; the image simply loses signal.
- **`z` quiet while any ADC is open.**  Through-slice dephasing is signal loss rather than an
  illegal block, so it shows up in none of the k-space checks above.
- **The two expected compiler warnings, pinned by content and count**, because a real warning
  hiding inside an expected one is the failure this costs nothing to prevent.
- **``FSE2D(echoes=1)`` reproduces ``SE2D``**, event for event.  Two notebooks are only safe if
  their divergence is something a test says out loud.
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
            # The two the notebooks themselves pin below: same-axis merges of gradients that do
            # not overlap in time, and three axes ramping together over the vector-norm limit.
            warnings.simplefilter('ignore', sc.SeqCraftWarning)
            for index, source in enumerate(sources):
                exec(compile(source, f'{notebook.name}:{index}', 'exec'), namespace)  # noqa: S102
    finally:
        os.chdir(here)
    return namespace


@pytest.fixture(scope='module')
def se(tmp_path_factory) -> dict:
    return _run(EXAMPLES / 'se_2d' / '01_build.ipynb', tmp_path_factory,
                stop_at='class SE2D(sc.Module):')


@pytest.fixture(scope='module')
def fse(tmp_path_factory) -> dict:
    return _run(EXAMPLES / 'fse_2d' / '01_build.ipynb', tmp_path_factory,
                stop_at='class FSE2D(sc.Module):')


def _digest(tree) -> list[tuple[float, str]]:
    """Every leaf event as (absolute time, content hash), tags and nesting ignored."""
    return sorted(
        (round(t, 12), sc.events.content_hash(event)) for t, event, _ in sc.flatten(tree)
    )


def _instance(fse: dict, *, echoes: int, **kwargs):
    """An ``FSE2D`` on the notebook's own protocol, at a turbo factor of `echoes`."""
    return fse['FSE2D'](
        opts=fse['opts'], fov_mm=fse['FOV_MM'], matrix=fse['MATRIX'],
        thickness_mm=fse['THICKNESS_MM'], echoes=echoes, tr_s=fse['TR_S'],
        bandwidth_hz_px=fse['BANDWIDTH_HZ_PX'],
        excitation_duration_s=fse['EXCITE_DURATION_S'],
        refocus_duration_s=fse['REFOCUS_DURATION_S'],
        refocus_thickness_factor=fse['REFOCUS_FACTOR'],
        crush_cycles_slice=fse['CRUSH_CYCLES_SLICE'],
        crush_cycles_readout=fse['CRUSH_CYCLES_READOUT'],
        spoil_cycles_per_voxel=fse['SPOIL_CYCLES'], spoil_axis=fse['SPOIL_AXIS'], **kwargs,
    )


def _interleaved(ny: int, echoes: int) -> list[list[int]]:
    shots = ny // echoes
    return [[s + n * shots for n in range(echoes)] for s in range(shots)]


# --------------------------------------------------------------------------- the headline
@pytest.mark.parametrize('echoes', [1, 8, 16, 72])
def test_k_is_exact_at_every_echo_of_every_train_length(fse, echoes) -> None:
    """
    The one assertion the whole sequence family exists to satisfy, at four turbo factors.

    ``1e-6`` 1/m is a real bound rather than a shrug: the measured worst case is 2.9e-7 over eight
    echoes and 2.6e-7 over seventy-two, so the tolerance is three times the float drift and a
    factor of ten million below one ``dk``.
    """
    ny = fse['NY']
    module = _instance(fse, echoes=echoes)
    segments = _interleaved(ny, echoes)
    lines = [line for segment in segments for line in segment]

    tree = module(segments=segments)
    k = sc.kspace(tree, fse['opts'])['k_adc'].reshape(3, len(lines), module.ro.num_samples)
    at_echo = k[:, :, module.ro.pre_echo_samples]
    wanted = (np.asarray(lines) - module.center_line) * (1e3 / fse['FOV_MM'][1])

    assert np.abs(at_echo[0]).max() < 1e-6, 'k_x: the dephaser sign and the two readout lobes'
    assert np.abs(at_echo[2]).max() < 1e-6, 'k_z: the crusher balance about the effective centre'
    assert np.abs(at_echo[1] - wanted).max() < 1e-6, 'k_y: the blip and the rewinder'


def test_the_refocusing_pulses_are_uniformly_spaced(fse) -> None:
    """
    Measured spread: 0 ps.  The stimulated-echo pathways depend on it, and a train whose pulses
    drift produces a smaller signal rather than a wrong one -- which is why arithmetic has to say
    so.  Per shot, because the TR gap between shots is not an echo spacing.
    """
    module = _instance(fse, echoes=16)
    segments = _interleaved(fse['NY'], 16)

    centres = sc.kspace(module(segments=segments), fse['opts'])['t_refocusing']
    per_shot = np.asarray(centres).reshape(len(segments), 16)

    assert np.ptp(np.diff(per_shot, axis=1)) < 1e-12
    assert np.diff(per_shot, axis=1).mean() == pytest.approx(module.echo_spacing_s, abs=1e-9)


def test_every_echo_is_within_one_gradient_raster_of_the_midpoint(fse) -> None:
    """
    Measured +1.95 us, and it is **not** zero on purpose: an ADC's delay lives on the 1 us RF
    raster and gradient starts live on the 10 us gradient raster, so the best available is a
    placement within one of them.  The test asserts the raster, not zero.
    """
    for echoes in (16, 72):
        module = _instance(fse, echoes=echoes)
        segments = _interleaved(fse['NY'], echoes)
        lines = [line for segment in segments for line in segment]

        k = sc.kspace(module(segments=segments), fse['opts'])
        centres = np.asarray(k['t_refocusing']).reshape(len(segments), echoes)
        echo_times = np.asarray(k['t_adc']).reshape(len(lines), module.ro.num_samples)[
            :, module.ro.pre_echo_samples].reshape(len(segments), echoes)

        midpoints = (centres[:, :-1] + centres[:, 1:]) / 2
        offset = echo_times[:, :-1] - midpoints

        assert np.abs(offset).max() <= float(fse['opts'].grad_raster_time), f'{echoes} echoes'


def test_z_is_quiet_while_every_adc_is_open(fse) -> None:
    """
    The one invariant here whose failure is **signal loss rather than a k-space error**, so
    nothing else in this file would catch it: a slice gradient on during the readout dephases
    through the slice, and every k check above still passes.
    """
    module = _instance(fse, echoes=16)
    shot = module(segments=_interleaved(fse['NY'], 16)[:1])

    grid, grads, marks = sc.sample(shot, fse['opts'])
    open_adc = np.zeros_like(grid, dtype=bool)
    for kind, start, end, _ in marks:
        if kind == 'adc':
            open_adc |= (grid >= start) & (grid < end)

    assert open_adc.any(), 'the probe has to contain an ADC for this to mean anything'
    assert np.abs(grads['z'][open_adc]).max() < 1e-9


# ------------------------------------------------------------------------------- the pin
def test_one_echo_of_the_train_is_the_spin_echo_the_other_notebook_writes(se, fse) -> None:
    """
    Event for event, acquired and dummy.  Two example directories are only safe if their
    divergence is something a test says out loud -- and this is what makes ``echoes=1`` a
    *generalisation* rather than a claim.
    """
    lines = list(range(se['NY']))
    single = _instance(fse, echoes=1)

    assert _digest(se['se'](lines=lines)) == _digest(
        single(segments=[[line] for line in lines]))
    assert _digest(se['se'](lines=lines[:4], dummy_shots=2)) == _digest(
        single(segments=[[line] for line in lines[:4]], dummy_shots=2))


def test_the_class_reproduces_the_loop_written_beside_it(se) -> None:
    """
    The notebook asserts this too; here it cannot be skipped by not running the cell.

    It is what keeps the class honest: if it grows a decision the cells above it do not make, this
    is what notices -- and those cells are what a reader actually learns from.
    """
    lines = list(range(se['NY']))

    assert _digest(se['scan_inline'](lines)) == _digest(se['se'](lines=lines))


# ------------------------------------------------------------------------- the sampling
def test_a_full_table_covers_k_space_exactly_once_and_the_labels_say_so(fse) -> None:
    """
    The round trip a reconstruction depends on: ``LIN`` is a permutation of ``range(ny)`` and each
    readout's measured ``k_y`` matches the label it carries.  Acquisition order is not k-space
    order in a segmented scan, and the file is the only record of which is which.
    """
    ny = fse['NY']
    module = _instance(fse, echoes=16)
    segments = _interleaved(ny, 16)
    tree = module(segments=segments)

    seq = sc.compile(tree, fse['opts'])
    labels = [int(v) for v in np.atleast_1d(np.asarray(
        seq.evaluate_labels(evolution='adc')['LIN']))]
    k = sc.kspace(tree, fse['opts'])['k_adc'].reshape(3, ny, module.ro.num_samples)

    assert sorted(labels) == list(range(ny))
    measured = k[1, :, module.ro.pre_echo_samples] / (1e3 / fse['FOV_MM'][1]) + module.center_line
    assert np.allclose(measured, labels, atol=1e-6)


def test_seg_splits_the_shots_and_no_echo_counter_is_emitted(fse) -> None:
    """
    ``SEG`` per shot; **``ECO`` deliberately absent**.  Every echo of an FSE shot goes into *one*
    image at a different ``ky``, so an echo counter would tell a reconstruction to split what
    belongs together -- and it is the right label for a multi-echo GRE, so the two look identical
    in a file.
    """
    module = _instance(fse, echoes=16)
    segments = _interleaved(fse['NY'], 16)

    seq = sc.compile(module(segments=segments), fse['opts'])
    labels = seq.evaluate_labels(evolution='adc')
    per_readout = np.asarray(labels['SEG']).reshape(len(segments), 16)

    assert (per_readout == np.arange(len(segments))[:, None]).all()
    assert 'ECO' not in labels or not np.any(np.asarray(labels['ECO']))


def test_te_eff_agrees_with_the_compiled_adc_for_two_orderings(fse) -> None:
    """
    A query that can disagree with its own block is worse than no query.  Two orderings, because
    they put ``center_line`` at opposite ends of the train: one that read the first echo's TE off
    ``self`` would be right for centric and wrong by a train for interleaved.
    """
    ny, echoes = fse['NY'], 16
    module = _instance(fse, echoes=echoes)
    shots = ny // echoes
    order = sorted(range(ny), key=lambda line: (abs(line - ny // 2), line))
    tables = {
        'interleaved': _interleaved(ny, echoes),
        'centric': [[order[n * shots + s] for n in range(echoes)] for s in range(shots)],
    }

    for name, segments in tables.items():
        lines = [line for segment in segments for line in segment]
        tree = module(segments=segments)
        k = sc.kspace(tree, fse['opts'])
        row = lines.index(module.center_line)
        t_echo = float(np.asarray(k['t_adc']).reshape(len(lines), module.ro.num_samples)[
            row, module.ro.pre_echo_samples])
        shot = row // echoes
        t_excitation = float(np.asarray(k['t_excitation'])[shot])

        assert t_echo - t_excitation == pytest.approx(module.te_eff_s(segments), abs=1e-5), name

    assert tables['interleaved'] != tables['centric']
    assert module.te_eff_s(tables['interleaved']) != module.te_eff_s(tables['centric'])


# ------------------------------------------------------------------------- the refusals
def test_the_table_refusals(fse) -> None:
    """Four things about a table the compiler cannot see, each with its own message."""
    ny, echoes = fse['NY'], 16
    module = _instance(fse, echoes=echoes)

    with pytest.raises(sc.ConfigurationError, match='non-empty'):
        module(segments=[])
    with pytest.raises(sc.ConfigurationError, match=f'do not have {echoes} lines'):
        module(segments=[list(range(echoes - 1))])
    with pytest.raises(sc.ConfigurationError, match='outside'):
        module(segments=[[*range(echoes - 1), ny]])
    with pytest.raises(sc.ConfigurationError, match='more than once'):
        module(segments=[list(range(echoes)), list(range(echoes))])


def test_a_scattered_echo_band_warns_and_still_builds(fse) -> None:
    """
    Warn, do not refuse: a table that scatters one echo index across k-space is a legitimate
    experiment as well as a common mistake, and the symptom -- ghosting along the phase-encode
    direction -- reads as motion rather than as a table problem.
    """
    ny, echoes = fse['NY'], 16
    module = _instance(fse, echoes=echoes)
    linear = [list(range(s * echoes, (s + 1) * echoes)) for s in range(ny // echoes)]

    with pytest.warns(sc.SeqCraftWarning, match='widely separated'):
        assert module(segments=linear)

    # A centric ordering puts its last echo at *both* edges of k-space, and must not warn: its
    # point-spread function is symmetric, so it blurs rather than ghosts, and a warning a reader
    # learns to ignore is worse than none.
    shots = ny // echoes
    order = sorted(range(ny), key=lambda line: (abs(line - ny // 2), line))
    centric = [[order[n * shots + s] for n in range(echoes)] for s in range(shots)]
    with warnings.catch_warnings():
        warnings.simplefilter('error', sc.SeqCraftWarning)
        assert module(segments=centric)


# ---------------------------------------------------------------------- the two warnings
def test_the_expected_warnings_are_pinned_by_content_and_count(fse) -> None:
    """
    Pinning is what stops a real warning hiding inside an expected one.

    The first is the readout's trailing lobe sharing a block with the tail spoiler on the same
    axis, and echo *n*'s rewinder with echo *n+1*'s blip -- exact superpositions of gradients that
    do not overlap in time.  The second is three axes ramping together at per-axis slew limits,
    which is legal on a real amplifier and over the vector norm.
    """
    module = _instance(fse, echoes=16)
    segments = _interleaved(fse['NY'], 16)[:2]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', sc.SeqCraftWarning)
        sc.compile(module(segments=segments), fse['opts'])

    messages = [str(w.message).splitlines()[0] for w in caught]

    assert len(messages) == 2, messages
    assert 'same-axis gradient merges' in messages[0]
    assert 'vector-norm limit' in messages[1]
    assert 'legal on real amplifiers' in messages[1]


def test_the_haste_shot_is_one_excitation_and_seventy_two_echoes(fse) -> None:
    """Partial Fourier is 0.625 rather than 0.6 because ``adc_samples_divisor`` wants a
    multiple of four, and 0.6 of 128 is 77."""
    module = _instance(fse, echoes=72, partial_fourier=0.625)
    segments = [[fse['CENTER_LINE'] - 8 + n for n in range(72)]]

    tree = module(segments=segments)
    excitations = sc.kspace(tree, fse['opts'])['t_excitation']

    assert module.ro.num_samples == 80
    assert module.ro.num_samples % 4 == 0
    assert len(excitations) == 1
    assert module.shot_s < 0.8, 'a breath hold, which is the whole point'
