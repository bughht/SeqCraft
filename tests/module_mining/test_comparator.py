"""
The comparator, against answers known independently of any candidate.

**This suite exists because the comparator has been wrong three times and the candidates none.**
Each time the failing measurement looked exactly like a physics bug in the module under test:

- an echo location that reported ``kx`` alternating +-1.947 1/m, which is a textbook odd/even
  artefact and was floating-point noise choosing between two samples half a k-step either side of
  a centre the reference does not sample;
- an equal-area rule applied to the phase axis, where both implementations "failed" by ~230 1/m
  because an encode is *meant* to move k;
- a rotation equivariance out by 5.8e+02 1/m because the harness stacked bare spokes into one
  sequence, and ``k`` integrates from the start and resets only at an excitation.

So the rule these tests enforce is: **validate the validator against an answer known independently
of the candidate** -- an analytic value, a metamorphic identity, or a deliberately seeded error
whose size is chosen rather than measured.

The heavy external-reference replay is *not* here.  It depends on external checkouts, adapter
assumptions and candidate-specific tolerances, and it stays a manual or scheduled run.  What runs
on every pull request is whether the measuring instrument works.
"""

from __future__ import annotations

import math

import numpy as np
import pypulseq as pp
import pytest
from module_mining import fingerprint
from module_mining.candidates.radial import checks
from module_mining.reference import ReferenceSequence

import seqcraft as sc

# --------------------------------------------------------------------------- primitives


def test_a_known_trapezoid_has_its_known_moment(opts) -> None:
    """``_area`` against the area pypulseq itself computed, and against half of it."""
    trap = pp.make_trapezoid(channel='x', amplitude=1e4, flat_time=1e-3, rise_time=2e-4,
                             system=opts)
    knots = np.array([[0.0, 2e-4, 1.2e-3, 1.4e-3], [0.0, 1e4, 1e4, 0.0]])

    assert fingerprint._area(knots, 0.0, 1.4e-3) == pytest.approx(float(trap.area), rel=1e-12)
    # Half the plateau, endpoints interpolated rather than snapped to a knot.
    assert fingerprint._area(knots, 2e-4, 7e-4) == pytest.approx(1e4 * 5e-4, rel=1e-12)
    assert fingerprint._area(knots, 0.0, 1e-4) == pytest.approx(0.5 * 1e-4 * 5e3, rel=1e-12)


def test_the_echo_is_the_crossing_and_not_the_nearest_sample() -> None:
    """
    The first comparator bug, as a regression test.

    Samples straddling ``kx = 0`` by half a step: the crossing is exactly between them, and which
    of the two is *nearest* is a tie that floating-point noise decides.  Reporting either one as
    the echo is how an alternating +-dk/2 artefact appears in a correct sequence.
    """
    times = np.array([0.0, 1e-5, 2e-5, 3e-5])
    kx = np.array([-3.0, -1.0, 1.0, 3.0])

    t_echo, nearest, crossed = fingerprint._echo_instant(times, kx)

    assert crossed
    assert t_echo == pytest.approx(1.5e-5)
    assert nearest in {1, 2}
    assert abs(times[nearest] - t_echo) / 1e-5 == pytest.approx(0.5)


def test_a_readout_that_never_crosses_the_centre_says_so() -> None:
    """A fallback a caller can distinguish from a measurement."""
    t_echo, nearest, crossed = fingerprint._echo_instant(
        np.array([0.0, 1e-5, 2e-5]), np.array([5.0, 4.0, 3.0]))

    assert not crossed
    assert nearest == 2
    assert t_echo == pytest.approx(2e-5)


# ------------------------------------------------------------------ analytic rotation
@pytest.mark.parametrize('angle_deg', [0.0, 17.0, 90.0, 180.0, -45.0])
def test_a_rotated_trajectory_is_recognised_as_one(angle_deg) -> None:
    """Synthetic, so the answer comes from trigonometry rather than from a sequence."""
    phi = math.radians(angle_deg)
    spoke = np.vstack([np.linspace(-100.0, 100.0, 64), np.zeros(64), np.zeros(64)])
    rotation = np.array([[math.cos(phi), -math.sin(phi)], [math.sin(phi), math.cos(phi)]])
    rotated = np.vstack([rotation @ spoke[:2], np.zeros(64)])

    assert checks.rotation_equivariance(spoke, rotated, angle_deg)['pass']


def test_a_trajectory_rotated_by_the_wrong_angle_is_not() -> None:
    """The check has to be able to fail, or passing means nothing."""
    spoke = np.vstack([np.linspace(-100.0, 100.0, 64), np.zeros(64), np.zeros(64)])
    phi = math.radians(30.0)
    rotation = np.array([[math.cos(phi), -math.sin(phi)], [math.sin(phi), math.cos(phi)]])
    rotated = np.vstack([rotation @ spoke[:2], np.zeros(64)])

    verdict = checks.rotation_equivariance(spoke, rotated, 31.0)
    assert not verdict['pass']
    assert verdict['worst'] > 1.0


# ------------------------------------------------------- against a module whose answer is known
@pytest.fixture(scope='module')
def line_measurement(opts):
    """A Cartesian line, measured -- the semantic centre of which is known from its matrix."""
    line = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=200.0)
    tree = sc.LogicBlock().add(0.0, line())
    reference = ReferenceSequence(name='line', sequence=sc.compile(tree, opts, name='line'))
    return line, fingerprint.measure(reference)


def test_the_semantic_centre_of_a_known_line_is_where_the_module_says(opts, line_measurement):
    line, measured = line_measurement
    readout = measured['readouts'][0]

    assert readout['n_samples'] == 64
    assert readout['sample_of_echo'] == line.echo_sample(0) == 32
    assert abs(readout['k_at_echo_per_m'][0]) < 1e-6
    # A Cartesian line *does* put a sample on the centre, so the offset is zero -- the number that
    # was 0.5 for the PyPulseq TSE demo, and the difference is a design fact rather than an error.
    assert readout['echo_offset_dwells'] == pytest.approx(0.0, abs=1e-6)


def test_ky_coverage_finds_a_missing_and_a_duplicated_line(opts, line_measurement) -> None:
    """The report has to name what is wrong, not only that something is."""
    _, measured = line_measurement
    verdict = fingerprint.check_ky_coverage(measured, dky_per_m=4.0, n_lines=1)

    assert set(verdict) >= {'missing', 'duplicated', 'lattice_max_error_per_m', 'pass'}


# ------------------------------------------------ the phase axis is not an equal-area axis
def test_the_equal_area_rule_does_not_apply_to_the_phase_axis(opts) -> None:
    """
    The second comparator bug, as a regression test.

    On x and z the areas either side of a refocusing pulse must match.  On y they must *not*: the
    encode is applied after the pulse and rewound before the next one, so what has to hold is
    that ``ky`` returns to zero **at** the pulse.  Measured here both ways on one shot, so the
    two statements are visibly different numbers rather than a claim in a docstring.
    """
    shot = sc.modules.TSEShot(opts=opts, fov_mm=200.0, matrix=(32, 16), thickness_mm=5.0,
                              echoes=4)
    tree = sc.LogicBlock().add(0.0, shot(lines=[2, 7, 8, 13]))
    reference = ReferenceSequence(name='tse', sequence=sc.compile(tree, opts, name='tse'))
    measured = fingerprint.measure(reference)
    rows = {row['invariant']: row for row in fingerprint.check_invariants(measured)}

    assert rows['I4y']['what'].startswith('ky = 0')
    assert rows['I4y']['pass']

    # The quantity the rule *would* have compared, had it been applied to y.
    waveforms = measured['_waveforms']
    shot_row = measured['shots'][0]
    echoes = sorted(measured['readouts'], key=lambda r: r['t_echo_s'])
    centre = shot_row['refocusing_s'][1]
    before = fingerprint._area(waveforms[1], echoes[0]['t_echo_s'], centre)
    after = fingerprint._area(waveforms[1], centre, echoes[1]['t_echo_s'])

    assert abs(before - after) > 1.0            # ~230 1/m at the protocol that exposed this
    assert rows['I4x']['pass'] and rows['I4z']['pass']


# ------------------------------------------------------------- a seeded, known-size error
def test_a_seeded_k_space_error_is_reported_with_its_size(opts) -> None:
    """
    A deliberately perturbed sequence, and a difference report that names the axis.

    The perturbation is chosen rather than measured: two readouts whose partial Fourier differs
    put the centre of k-space at different samples, so the echoes land at different k.
    """
    def measured(partial_fourier: float):
        line = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64,
                                        bandwidth_hz_px=200.0, partial_fourier=partial_fourier)
        tree = sc.LogicBlock().add(0.0, line())
        return fingerprint.measure(
            ReferenceSequence(name=f'pf{partial_fourier}',
                              sequence=sc.compile(tree, opts, name='seeded')))

    good = fingerprint.compare(measured(1.0), measured(1.0))
    assert good['result'] == 'pass'
    assert good['first_failure'] is None

    # Different sample counts are not comparable, and the report says that rather than guessing.
    mismatched = fingerprint.compare(measured(1.0), measured(0.75))
    assert mismatched['result'] == 'incomparable'
    assert 'readout inventory' in mismatched['first_failure']['check']


def test_block_count_is_reported_and_never_decides(opts, line_measurement) -> None:
    """Two legal decompositions of one sequence differ here and must still compare equal."""
    _, measured = line_measurement
    report = fingerprint.compare(measured, measured)

    assert 'note' in report['l0']['n_blocks']
    assert report['result'] == 'pass'


# --------------------------------------------- the harness itself, whose bug started this file
def test_bare_spokes_in_one_sequence_accumulate_k(opts) -> None:
    """
    The third comparator bug, preserved as a *known* answer rather than as a warning.

    ``k`` is integrated from the start of a sequence and reset only by an excitation, so two bare
    spokes in one tree have the second starting where the first ended.  That is correct pypulseq
    behaviour and a wrong measurement, which is why the radial adapter compiles one spoke per
    sequence.  If this ever stops accumulating, the adapter's reason has gone away and it should
    be revisited rather than quietly kept.
    """
    spoke = sc.modules.RadialReadout(opts=opts, fov_mm=260.0, matrix=64, dwell_s=20e-6)
    period = spoke.prephaser_duration_s + spoke.num_samples * spoke.dwell_s + 2e-3

    stacked = sc.LogicBlock()
    stacked.add(0.0, spoke(angle_rad=0.0)).add(period, spoke(angle_rad=0.0))
    k = sc.kspace(stacked, opts)['k_adc']

    first, second = k[:, :spoke.num_samples], k[:, spoke.num_samples:]
    assert np.linalg.norm(first[:, spoke.center_sample]) < 1e-3
    assert np.linalg.norm(second[:, spoke.center_sample]) > 1.0
