"""
``FSE2D`` -- the table, and everything that is a property of it rather than of a waveform.

The physics is ``TSEShot``'s and is tested there.  What is left here is acquisition policy: that a
table covers k-space exactly once and the labels say which line went where, that dummies sample
nothing while loading the same gradients, that the effective TE follows the echo carrying the
centre line, and that the two mistakes a table can make are caught -- one refused, one warned
about, and the difference between them is deliberate.
"""

from __future__ import annotations

import numpy as np
import pytest

import seqcraft as sc

MATRIX = (32, 16)
ECHOES = 4


@pytest.fixture(scope='module')
def fse(opts):
    """Four shots of four echoes over sixteen lines."""
    return sc.modules.FSE2D(opts=opts, fov_mm=200.0, matrix=MATRIX, thickness_mm=5.0,
                            echoes=ECHOES, tr_s=0.5)


def _interleaved(echoes: int = ECHOES, n_lines: int = MATRIX[1]) -> list[list[int]]:
    shots = n_lines // echoes
    return [[s + n * shots for n in range(echoes)] for s in range(shots)]


def _labels(tree, name: str) -> list[int]:
    return [event.value for _, event, _ in sc.flatten(tree)
            if getattr(event, 'type', '') == 'labelset' and event.label == name]


# ------------------------------------------------------------------------- the table
def test_a_full_table_covers_k_space_once_and_the_labels_say_so(opts, fse) -> None:
    """Signed ``ky`` at each echo against the line its own ``LIN`` label claims."""
    segments = _interleaved()
    tree = fse(segments=segments)
    lines = [line for segment in segments for line in segment]

    k = sc.kspace(tree, opts)['k_adc'].reshape(3, len(lines), fse.shot.ro.num_samples)
    at_echo = k[1, :, fse.shot.ro.echo_sample(0)]
    want = (np.asarray(lines) - fse.center_line) * (1e3 / 200.0)

    assert _labels(tree, 'LIN') == lines
    assert sorted(lines) == list(range(MATRIX[1]))
    assert np.abs(at_echo - want).max() < 1e-3


def test_one_seg_label_per_acquired_shot(opts, fse) -> None:
    """Counted from the first acquired shot, so a dummy does not consume an index."""
    tree = fse(segments=_interleaved(), dummy_shots=2)
    assert _labels(tree, 'SEG') == list(range(MATRIX[1] // ECHOES))


def test_dummies_sample_nothing_and_load_the_same_gradients(opts, fse) -> None:
    """
    Shot 0 starts from equilibrium and the rest do not, and with interleaved segmentation that
    one difference is a comb in k-space -- which is a replica of the object, not a blur.
    """
    segments = _interleaved()
    plain = fse(segments=segments)
    with_dummies = fse(segments=segments, dummy_shots=1)

    kinds = [getattr(event, 'type', '') for _, event, _ in sc.flatten(with_dummies)]
    assert kinds.count('adc') == len(segments) * ECHOES
    assert with_dummies.duration == pytest.approx(plain.duration + fse.shot.tr_s)


def test_shots_are_stacked_at_tr(opts, fse) -> None:
    starts = sorted(node.start for node in fse(segments=_interleaved()))
    assert np.allclose(np.diff(starts), fse.shot.tr_s)


# ---------------------------------------------------------------- the effective echo time
def test_the_effective_te_follows_the_echo_carrying_the_centre_line(opts, fse) -> None:
    """
    The whole point of an ordering: the same lines and the same waveforms, a different contrast.
    """
    interleaved = _interleaved()
    linear = [list(range(s * ECHOES, (s + 1) * ECHOES)) for s in range(MATRIX[1] // ECHOES)]

    assert fse.echo_of_center_line(interleaved) == (0, 2)
    assert fse.te_eff_s(interleaved) == pytest.approx(fse.shot.te_s(2))
    assert fse.te_eff_s(linear) == pytest.approx(fse.shot.te_s(fse.echo_of_center_line(linear)[1]))
    assert fse.te_eff_s(interleaved) != fse.te_eff_s(linear)


def test_the_effective_te_is_the_measured_time_to_the_centre_of_k_space(opts, fse) -> None:
    """The claim against the sequence, rather than against the expression that produced it."""
    segments = _interleaved()
    tree = fse(segments=segments)
    lines = [line for segment in segments for line in segment]
    k = sc.kspace(tree, opts)
    samples = fse.shot.ro.num_samples

    at_echo = k['k_adc'].reshape(3, len(lines), samples)[1, :, fse.shot.ro.echo_sample(0)]
    centre = int(np.argmin(np.abs(at_echo)))
    t_echo = k['t_adc'].reshape(len(lines), samples)[centre, fse.shot.ro.echo_sample(0)]
    shot_start = (centre // ECHOES) * fse.shot.tr_s

    assert t_echo - shot_start == pytest.approx(
        fse.te_eff_s(segments) + fse.shot.exc.time_to_center(), abs=1e-9)


def test_echo_bands_measure_distance_from_the_centre_not_signed_ky(opts, fse) -> None:
    """
    Centric deals lines outward from ``k = 0``, so its last echo sits at *both* edges: a span of
    the whole matrix in signed ``ky`` and a narrow band in distance.  Its point spread is
    symmetric and it blurs, which is not what the warning is for.
    """
    order = sorted(range(MATRIX[1]), key=lambda line: (abs(line - MATRIX[1] // 2), line))
    shots = MATRIX[1] // ECHOES
    centric = [[order[n * shots + s] for n in range(ECHOES)] for s in range(shots)]
    assert max(fse.echo_bands(centric).values()) <= max(1, MATRIX[1] // ECHOES)


# -------------------------------------------------------------------------- the refusals
def test_a_line_acquired_twice_raises(opts, fse) -> None:
    """Two readouts writing one k-space address; the compiler says less about it downstream."""
    with pytest.raises(sc.ConfigurationError, match='more than once'):
        fse(segments=[[0, 1, 2, 3], [0, 4, 5, 6]])


def test_an_empty_table_raises(opts, fse) -> None:
    with pytest.raises(sc.ConfigurationError, match='non-empty'):
        fse(segments=[])


def test_a_short_segment_raises_with_the_shots_own_message(opts, fse) -> None:
    """One contract, asked once: length and range are the shot's to check."""
    with pytest.raises(sc.ConfigurationError, match='refocusing pulses'):
        fse(segments=[[0, 1, 2]])


def test_negative_dummies_raise(opts, fse) -> None:
    with pytest.raises(sc.ConfigurationError, match='negative'):
        fse(segments=_interleaved(), dummy_shots=-1)


def test_a_table_without_the_centre_line_warns(opts, fse) -> None:
    """No DC term, and ``te_eff_s`` undefined -- but a legitimate thing to build."""
    segments = [s for s in _interleaved() if fse.center_line not in s]
    with pytest.warns(sc.SeqCraftWarning, match='centre line'):
        fse(segments=segments)
    assert fse.te_eff_s(segments) is None


def test_a_table_that_scatters_one_echo_across_k_space_warns(opts, fse) -> None:
    """
    Warned rather than refused, because it is the same mechanism the dummy shot addresses at a
    different period, and a deliberate experiment is a legitimate thing to build.
    """
    linear = [list(range(s * ECHOES, (s + 1) * ECHOES)) for s in range(MATRIX[1] // ECHOES)]
    with pytest.warns(sc.SeqCraftWarning, match='widely separated'):
        fse(segments=linear)


def test_the_component_contract(opts, fse, component_checks) -> None:
    component_checks.all(fse, segments=_interleaved())
