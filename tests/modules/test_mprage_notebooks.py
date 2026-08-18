"""
The two composites that stay in their notebooks, asserted from outside them.

``MPRAGE2D`` and ``MP2RAGE2D`` have one consumer each -- their own build notebook -- so neither
ships, and neither can be imported here.  What CI can still do is **run the notebook and assert
against what it defined**, which is the same bargain
``test_notebook_matches_the_package.py`` strikes for ``GRE2DTR``: the tutorial stays a tutorial,
and it cannot drift without something noticing.

The claims worth a test rather than a cell are the ones whose failure is *silent*:

- TI is inversion centre to the acquisition of ``k = 0``, for **both** orderings.  A placement
  that ignored ``time_to_center_line`` passes for centric and fails for linear.
- Every shot's *first* train carries ``SET=0``.  Pulseq labels are stateful, so a single ``SET=0``
  at the start of the sequence mislabels every shot after the first -- and the file stays legal.
- The second flip angle of a pair reaches the scanner.  A pair collapsing to its first element is
  the failure mode of tuple arguments and is invisible in the tree.
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
            # Every compile here warns about the vector-norm limit, which is legal on real
            # amplifiers and is the same warning gre_2d/ produces.
            warnings.simplefilter('ignore', sc.SeqCraftWarning)
            for index, source in enumerate(sources):
                exec(compile(source, f'{notebook.name}:{index}', 'exec'), namespace)  # noqa: S102
    finally:
        os.chdir(here)
    return namespace


@pytest.fixture(scope='module')
def mprage(tmp_path_factory) -> dict:
    return _run(EXAMPLES / 'mprage_2d' / '01_build.ipynb', tmp_path_factory,
                stop_at='class MPRAGE2D(sc.Module):')


@pytest.fixture(scope='module')
def mp2rage(tmp_path_factory) -> dict:
    return _run(EXAMPLES / 'mp2rage_2d' / '01_build.ipynb', tmp_path_factory,
                stop_at='class MP2RAGE2D(sc.Module):')


def _digest(tree) -> list[tuple[float, str]]:
    """Every leaf event as (absolute time, content hash), tags and nesting ignored."""
    return sorted(
        (round(t, 12), sc.events.content_hash(event)) for t, event, _ in sc.flatten(tree)
    )


def _labels(tree, opts, key: str) -> list[int]:
    """One label per acquired readout, read off the compiled sequence."""
    seq = sc.compile(tree, opts)
    return [int(v) for v in np.atleast_1d(np.asarray(
        seq.evaluate_labels(evolution='adc')[key]))]


# ------------------------------------------------------------------------------ MPRAGE
@pytest.mark.parametrize('ordering', ['linear', 'centric'])
def test_the_centre_line_lands_exactly_ti_after_the_inversion(mprage, ordering) -> None:
    """
    Both orderings, because that is what distinguishes a correct placement from a lucky one.

    A sequence that referenced TI to the train's *start* rather than to the acquisition of
    ``k = 0`` gives the right answer for centric -- where they nearly coincide -- and is wrong by
    half a train for linear.
    """
    module, segments = mprage['mprage'], mprage['segmentations'][ordering]
    tree = module(segments=segments)

    acquired = [line for segment in segments for line in segment]
    k = sc.kspace(tree, mprage['opts'])
    nx, center = mprage['NX'], module.center_line
    t_echo = float(k['t_adc'].reshape(len(acquired), nx)[acquired.index(center), nx // 2])

    shot = acquired.index(center) // mprage['TURBO'] + module.dummy_shots
    interval = float(mprage['block_raster'].ceil(module.shot_interval_s))
    t_inversion = shot * interval + module.inversion.time_to_center()

    assert t_echo - t_inversion == pytest.approx(module.ti_s, abs=1e-5)


def test_the_class_and_the_loop_beside_it_are_the_same_events(mprage) -> None:
    """
    The notebook asserts this too; here it cannot be skipped by not running the cell.

    It is the assertion that keeps the class honest: if it grows a decision the loop above it does
    not make, this is what notices -- and the loop is what a reader actually learns from.
    """
    for segments in mprage['segmentations'].values():
        assert _digest(mprage['scan'](segments)) == _digest(
            mprage['mprage'](segments=segments)
        )


def test_a_dummy_shot_is_a_whole_shot_and_samples_nothing(mprage) -> None:
    module, segments = mprage['mprage'], mprage['segmentations']['centric']

    tree = module(segments=segments)
    adcs = [e for _, e, _ in sc.flatten(tree) if getattr(e, 'type', '') == 'adc']
    rfs = [e for _, e, _ in sc.flatten(tree) if getattr(e, 'type', '') == 'rf']

    shots = len(segments) + module.dummy_shots
    assert len(adcs) == sum(len(segment) for segment in segments)
    assert len(rfs) == shots * (1 + len(segments[0])), 'each dummy shot inverts and excites too'


def test_the_provenance_path_is_four_levels_and_reads(mprage) -> None:
    """The deepest nesting the package has produced, and it still names what it came from."""
    tree = mprage['mprage'](segments=mprage['segmentations']['centric'])

    paths = {path for _, _, path in sc.flatten(tree)}

    assert ('MPRAGE2D', 'GRE2D', 'GRE2DTR', 'CartesianLine') in paths
    assert ('MPRAGE2D', 'IRPrep', 'spoiler') in paths


def test_a_ti_below_the_floor_is_refused(mprage) -> None:
    """Not merely tight: below it the first excitation overlaps the inversion crusher."""
    segments = mprage['segmentations']['linear']
    floor = mprage['mprage'].min_ti_s(segments=segments)
    too_short = mprage['MPRAGE2D'](
        opts=mprage['opts'], fov_mm=mprage['FOV_MM'], matrix=mprage['MATRIX'],
        thickness_mm=mprage['THICKNESS_MM'], ti_s=floor / 2,
        shot_interval_s=mprage['SHOT_INTERVAL_S'], inversion_pulse=mprage['INVERSION_PULSE'],
        inversion_duration_s=mprage['INVERSION_DURATION_S'],
    )

    with pytest.raises(sc.ConfigurationError, match='below the minimum'):
        too_short(segments=segments)


# ----------------------------------------------------------------------------- MP2RAGE
def test_both_inversion_times_land_where_they_were_asked_to(mp2rage) -> None:
    module, segments = mp2rage['mp2rage'], mp2rage['segments']
    tree = module(segments=segments)

    k = sc.kspace(tree, mp2rage['opts'])
    nx, center = mp2rage['NX'], module.center_line
    lines = [line for segment in segments for line in segment]
    per_readout = k['t_adc'].reshape(2 * len(lines), nx)

    interval = float(mp2rage['block_raster'].ceil(module.shot_interval_s))
    shot = lines.index(center) // mp2rage['TURBO'] + module.dummy_shots
    t_inversion = shot * interval + module.inversion.time_to_center()
    # Readouts run train 0 then train 1 within each shot, both over the same segment.
    within = lines.index(center) % mp2rage['TURBO']
    base = (lines.index(center) // mp2rage['TURBO']) * 2 * mp2rage['TURBO']

    for contrast, ti in enumerate(module.ti_s):
        row = base + contrast * mp2rage['TURBO'] + within
        measured = float(per_readout[row, nx // 2]) - t_inversion
        assert measured == pytest.approx(ti, abs=1e-5), f'TI{contrast}'


def test_every_shot_labels_its_first_train_set_zero(mp2rage) -> None:
    """
    The stateful-label bug, and the reason it needs **more than one shot** to show up.

    ``type='SET'`` writes a value that persists.  Emitting ``SET=0`` once at the start of the
    sequence and ``SET=1`` before each second train leaves shot 2's first train wearing the
    ``SET=1`` shot 1's second train wrote -- and the file is still legal.
    """
    module, segments = mp2rage['mp2rage'], mp2rage['segments']
    assert len(segments) > 1, 'a single-shot check passes with the bug present'

    values = np.asarray(_labels(module(segments=segments), mp2rage['opts'], 'SET'))
    per_shot = values.reshape(len(segments), 2, len(segments[0]))

    assert (per_shot[:, 0, :] == 0).all(), 'a later shot inherited SET=1'
    assert (per_shot[:, 1, :] == 1).all()
    assert sorted(np.bincount(values).tolist()) == [len(values) // 2] * 2


def test_the_stateful_version_is_refused_by_the_compiler(mp2rage) -> None:
    """
    The label-address check catches it: 96 readouts sharing one ``(LIN, SET)`` address.

    That is the same check that catches a repeated phase-encode line, doing a job nobody wrote it
    for -- and it is why the notebook builds the broken version rather than describing it.
    """
    with pytest.raises(sc.CompileError, match='repeat a k-space address'):
        sc.compile(mp2rage['buggy_scan'](mp2rage['segments']), mp2rage['opts'])


def test_both_trains_acquire_the_same_lines(mp2rage) -> None:
    """Which is what makes the ratio in ``02`` defined pixel by pixel."""
    module, segments = mp2rage['mp2rage'], mp2rage['segments']

    lines = np.asarray(_labels(module(segments=segments), mp2rage['opts'], 'LIN'))
    per_shot = lines.reshape(len(segments), 2, len(segments[0]))

    assert (per_shot[:, 0, :] == per_shot[:, 1, :]).all()


def test_the_second_flip_angle_reaches_the_scanner(mp2rage) -> None:
    """A pair collapsing to its first element is the failure mode of tuple arguments."""
    module = mp2rage['mp2rage']
    tree = module(segments=mp2rage['segments'])

    amplitudes = sorted({round(float(np.abs(e.signal).max()), 6)
                         for _, e, _ in sc.flatten(tree) if getattr(e, 'type', '') == 'rf'})

    assert len(amplitudes) == 3, 'the inversion and two distinct excitations'
    assert module.flip_deg[0] != module.flip_deg[1]
    assert amplitudes[1] / amplitudes[0] == pytest.approx(
        module.flip_deg[1] / module.flip_deg[0], rel=1e-3,
    )


def test_ti_s_must_be_increasing_and_far_enough_apart(mp2rage) -> None:
    def variant(**kwargs):
        base = dict(opts=mp2rage['opts'], fov_mm=mp2rage['FOV_MM'], matrix=mp2rage['MATRIX'],
                    thickness_mm=mp2rage['THICKNESS_MM'], ti_s=mp2rage['TI_S'],
                    flip_deg=mp2rage['FLIP_DEG'],
                    shot_interval_s=mp2rage['SHOT_INTERVAL_S'],
                    inversion_pulse=mp2rage['INVERSION_PULSE'],
                    inversion_duration_s=mp2rage['INVERSION_DURATION_S'])
        return mp2rage['MP2RAGE2D'](**{**base, **kwargs})

    with pytest.raises(sc.ConfigurationError, match='increasing'):
        variant(ti_s=(mp2rage['TI_S'][1], mp2rage['TI_S'][0]))
    with pytest.raises(sc.ConfigurationError, match='overlap'):
        variant(ti_s=(0.70, 0.75))(segments=mp2rage['segments'])
