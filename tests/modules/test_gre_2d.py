"""
``GRE2D`` -- the loop, the schedule, and the validation of what is being acquired.

The three things a whole acquisition owns and a single repetition cannot: how many repetitions
there are, what phase each one gets, and which lines the table contains.  The phase counter is
the one with a quiet failure mode, so it is checked against the closed form *and* against the
compiled sequence.
"""

from __future__ import annotations

import numpy as np
import pytest

import seqcraft as sc

MATRIX = (64, 32)
NY = MATRIX[1]


@pytest.fixture(scope='module')
def gre(opts):
    """A small but complete scan, with no dummies, so a test can add its own."""
    return sc.modules.GRE2D(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0)


def _rf_phases(block) -> list[float]:
    """Carrier phase of every RF event in the tree, degrees, in time order."""
    return [
        float(np.rad2deg(e.phase_offset))
        for _, e, _ in sorted(sc.flatten(block), key=lambda item: item[0])
        if getattr(e, 'type', '') == 'rf'
    ]


# ---------------------------------------------------------------------- the loop
def test_one_repetition_per_line_stacked_at_tr(opts, gre) -> None:
    scan = gre(lines=range(NY))

    assert len(scan) == NY
    assert [n.start for n in scan] == pytest.approx(
        [n * gre.tr.tr_s for n in range(NY)], abs=1e-12,
    )
    assert scan.duration == pytest.approx(NY * gre.tr.tr_s, abs=1e-9)


def test_the_pattern_is_a_build_argument_so_one_instance_serves_all_of_them(opts, gre) -> None:
    """
    The ``__init__`` designs / ``build`` assembles split doing its job.

    The waveforms are created once; a three-way comparison is three cheap calls on one object,
    which is what makes an undersampling study a study rather than three sequences.
    """
    acs = set(range(gre.center_line - 4, gre.center_line + 4))
    patterns = {
        'full': tuple(range(NY)),
        'uniform': tuple(sorted(acs | set(range(0, NY, 3)))),
        'reversed': tuple(reversed(range(NY))),
    }

    built = {name: gre(lines=lines) for name, lines in patterns.items()}

    assert len(built['full']) == NY
    assert len(built['uniform']) == len(patterns['uniform'])
    assert built['reversed'].duration == built['full'].duration


def test_acquisition_order_is_the_order_of_the_list(opts, gre) -> None:
    """``lines=reversed(...)`` is line ordering, with no new argument."""
    scan = gre(lines=reversed(range(NY)))

    seq = sc.compile(scan, opts)
    seen = [int(v) for v in np.atleast_1d(np.asarray(seq.evaluate_labels(evolution='adc')['LIN']))]

    assert seen == list(reversed(range(NY)))


def test_lines_has_no_default(opts, gre) -> None:
    """Fully sampled is ``range(ny)`` written out, so undersampling reads as a peer choice."""
    with pytest.raises(TypeError, match='lines'):
        gre()


# ------------------------------------------------------------------ RF spoiling
def test_the_phase_follows_the_quadratic_schedule(opts, gre) -> None:
    assert [gre.phase_deg(n) for n in range(5)] == [0.0, 117.0, 351.0, 702.0, 1170.0]


def test_the_phase_counter_runs_across_the_dummies(opts) -> None:
    """
    Repetition *n* is counted from the first dummy, not from the first acquired line.

    Restarting it at the acquisition is a real and quiet bug: the steady state the dummies
    established is a steady state *of a particular phase sequence*, and resetting the counter
    discards it at the moment it starts to matter.
    """
    with_dummies = sc.modules.GRE2D(opts=opts, fov_mm=250.0, matrix=MATRIX,
                                    thickness_mm=5.0, dummies=3)

    phases = _rf_phases(with_dummies(lines=range(NY)))

    assert len(phases) == NY + 3
    expected = [(0.5 * 117.0 * n * (n + 1)) % 360.0 for n in range(NY + 3)]
    assert [p % 360.0 for p in phases] == pytest.approx(expected, abs=1e-9)
    assert phases[3] != pytest.approx(0.0), (
        'the first acquired repetition continues the schedule rather than restarting it'
    )


def test_rf_spoiling_off_gives_every_repetition_the_same_phase(opts) -> None:
    """Gradient spoiling is a separate mechanism, and is unaffected."""
    plain = sc.modules.GRE2D(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                             rf_spoil=False)

    assert set(_rf_phases(plain(lines=range(NY)))) == {0.0}
    assert sc.moments(plain(lines=[16]))['z'] > 0.0, 'the gradient spoiler is still there'


def test_a_spoil_phase_that_cannot_take_effect_raises(opts) -> None:
    """The same shape as ``axis`` alongside a non-selective ``Excitation``."""
    with pytest.raises(sc.ConfigurationError, match='rf_spoil_deg'):
        sc.modules.GRE2D(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                         rf_spoil=False, rf_spoil_deg=117.0)


def test_a_custom_spoil_increment_is_used(opts) -> None:
    gre = sc.modules.GRE2D(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                           rf_spoil_deg=50.0)

    assert gre.phase_deg(2) == pytest.approx(150.0)


# ---------------------------------------------------------------------- dummies
def test_dummies_come_first_and_sample_nothing(opts) -> None:
    dummied = sc.modules.GRE2D(opts=opts, fov_mm=250.0, matrix=MATRIX,
                               thickness_mm=5.0, dummies=4)

    scan = dummied(lines=range(NY))
    seq = sc.compile(scan, opts)

    assert len(scan) == NY + 4
    adcs = [e for _, e, _ in sc.flatten(scan) if getattr(e, 'type', '') == 'adc']
    assert len(adcs) == NY, 'four repetitions played and none of them sampled'
    labels = [int(v) for v in np.atleast_1d(
        np.asarray(seq.evaluate_labels(evolution='adc')['LIN']))]
    assert labels == list(range(NY)), 'and none of them wrote a line index'


def test_a_dummy_is_gradient_identical_to_the_repetition_it_precedes(opts) -> None:
    """What it is establishing has to be the steady state that then gets acquired."""
    dummied = sc.modules.GRE2D(opts=opts, fov_mm=250.0, matrix=MATRIX,
                               thickness_mm=5.0, dummies=1, rf_spoil=False)

    scan = dummied(lines=[7])
    first, second = (n.item for n in scan)

    assert sc.moments(first) == pytest.approx(sc.moments(second))
    assert first.duration == pytest.approx(second.duration)


# -------------------------------------------------------------------- validation
def test_an_out_of_range_line_raises(opts, gre) -> None:
    with pytest.raises(sc.ConfigurationError, match=f'0 ... {NY - 1}'):
        gre(lines=[0, NY])


def test_a_repeated_line_raises(opts, gre) -> None:
    """Two readouts writing one k-space address; the compiler would reject it with less to say."""
    with pytest.raises(sc.ConfigurationError, match='more than once'):
        gre(lines=[0, 16, 16])


def test_an_empty_table_raises(opts, gre) -> None:
    with pytest.raises(sc.ConfigurationError, match='nothing to acquire'):
        gre(lines=[])


def test_a_pattern_missing_the_centre_of_k_space_warns(opts, gre) -> None:
    """
    Not illegal, and warned rather than refused.

    An image reconstructed with no DC term reads as a windowing problem rather than as a
    sampling bug, and with no generators shipped this warning is the only thing between a caller
    and that.
    """
    with pytest.warns(sc.SeqCraftWarning, match='centre of k-space'):
        gre(lines=[n for n in range(NY) if n != gre.center_line])


def test_a_pattern_carrying_the_centre_is_silent(opts, gre) -> None:
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter('error', sc.SeqCraftWarning)
        gre(lines=[gre.center_line])


def test_negative_dummies_raise(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='dummies'):
        sc.modules.GRE2D(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0, dummies=-1)


# ---------------------------------------------------------------- the whole thing
def test_the_scan_matches_a_hand_written_loop(opts, gre) -> None:
    """
    ``GRE2D`` is the loop and nothing else, so writing the loop out must give the same tree.

    This is the assertion that keeps the class honest: if it ever grows a decision the notebook
    version does not make, this is what notices.
    """
    lines = tuple(range(0, NY, 2)) + (gre.center_line,)
    lines = tuple(sorted(set(lines)))

    by_hand = sc.LogicBlock('GRE2D')
    for n, line in enumerate(lines):
        by_hand.add(n * gre.tr.tr_s, gre.tr(line=line, phase_deg=0.5 * 117.0 * n * (n + 1)))

    from_module = gre(lines=lines)
    left = sorted(sc.flatten(by_hand), key=lambda item: (item[0], str(item[2])))
    right = sorted(sc.flatten(from_module), key=lambda item: (item[0], str(item[2])))

    assert len(left) == len(right)
    for (t_a, a, path_a), (t_b, b, path_b) in zip(left, right, strict=True):
        assert t_a == pytest.approx(t_b, abs=1e-12)
        assert path_a == path_b
        assert sc.events.content_hash(a) == sc.events.content_hash(b)


def test_the_sampled_trajectory_covers_the_requested_lines(opts, gre) -> None:
    """
    ``sc.kspace`` over the whole scan, reduced to the ky each readout sat on.

    The arithmetic check the build notebook makes, in the suite so it cannot rot: the lines that
    were asked for are the lines that ended up in the file.
    """
    lines = tuple(sorted({*range(0, NY, 3), gre.center_line}))
    k = sc.kspace(gre(lines=lines), opts)

    per_readout = k['k_adc'][1].reshape(len(lines), -1)
    encoded = np.round(per_readout[:, 0] * (250.0 / 1e3)).astype(int) + gre.center_line

    assert list(encoded) == list(lines)


def test_it_is_pure_and_compiles_alone(opts, gre, component_checks) -> None:
    component_checks.all(gre, lines=range(0, NY, 4))
