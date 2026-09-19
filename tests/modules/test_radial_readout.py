"""
``RadialReadout`` -- the trajectory, measured rather than asserted about.

A spoke's failures are quiet ones.  A wrong angle is a rotated image; a wrong *increment* is
streaks; a centre sample half a step out puts every spoke's centre somewhere slightly different
and the reconstruction blurs rather than complains.  None of it shows up in a k-space extent
check, so everything here is measured against pypulseq's own integration of the compiled
waveform.

The rotation test is the one worth having twice: it was written against the official PyPulseq
reference in ``tools/module_mining/`` and passed there **before this module existed**, so it
cannot be encoding this module's own behaviour.  This is the same statement, in the package.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import seqcraft as sc

FOV_MM = 260.0
MATRIX = 64
DWELL_S = 20e-6


@pytest.fixture(scope='module')
def spoke(opts):
    """A full spoke on the official reference's geometry."""
    return sc.modules.RadialReadout(opts=opts, fov_mm=FOV_MM, matrix=MATRIX, dwell_s=DWELL_S)


def _k(readout, opts, angle_rad: float = 0.0) -> np.ndarray:
    """The compiled trajectory at the ADC samples, for one spoke alone."""
    return sc.kspace(sc.LogicBlock().add(0.0, readout(angle_rad=angle_rad)), opts)['k_adc']


# ------------------------------------------------------------------- the trajectory
def test_a_sample_lands_on_the_centre_of_k_space(opts, spoke) -> None:
    """
    Every spoke shares this one point, and it is the most densely sampled in the acquisition.

    A design whose samples straddle it puts each spoke's centre somewhere slightly different --
    which no extent check and no single-spoke inspection would show.
    """
    k = _k(spoke, opts)
    assert np.linalg.norm(k[:, spoke.center_sample]) < 1e-3
    assert spoke.center_sample == 32                       # not 31, and not the ADC midpoint


def test_the_module_says_where_its_own_trajectory_goes(opts, spoke) -> None:
    """The claim against the sequence, rather than against the expression that produced it."""
    k = _k(spoke, opts)
    assert k[0, 0] == pytest.approx(spoke.k_first_per_m, abs=1e-6)
    assert k[0, -1] == pytest.approx(spoke.k_last_per_m, abs=1e-6)
    assert np.mean(np.diff(k[0])) == pytest.approx(spoke.dk_per_m, abs=1e-9)
    assert spoke.dk_per_m == pytest.approx(1e3 / FOV_MM)
    assert spoke.k_max_per_m == pytest.approx(abs(spoke.k_first_per_m))


def test_time_to_center_matches_the_sequence(opts, spoke) -> None:
    """Not measurable from the tree: a block knows when events play, not which one is k = 0."""
    tree = sc.LogicBlock().add(0.0, spoke(angle_rad=0.4))
    k = sc.kspace(tree, opts)
    assert k['t_adc'][spoke.center_sample] == pytest.approx(spoke.time_to_center(), abs=1e-9)


@pytest.mark.parametrize(('partial_fourier', 'samples', 'centre'),
                         [(1.0, 64, 32), (0.75, 48, 16), (0.5, 32, 0)])
def test_partial_fourier_walks_the_centre_to_the_first_sample(
    opts, partial_fourier, samples, centre,
) -> None:
    """
    Full spoke to centre-out is one parameter, and ``dk`` does not move while it happens.

    The official UTE example spans the same continuum with an argument running the other way;
    this is SeqCraft's existing vocabulary saying the same thing.
    """
    readout = sc.modules.RadialReadout(opts=opts, fov_mm=FOV_MM, matrix=MATRIX, dwell_s=DWELL_S,
                                       partial_fourier=partial_fourier)
    k = _k(readout, opts)

    assert (readout.num_samples, readout.center_sample) == (samples, centre)
    assert np.linalg.norm(k[:, readout.center_sample]) < 1e-3
    assert readout.dk_per_m == pytest.approx(1e3 / FOV_MM)
    assert k[0, -1] == pytest.approx(readout.k_last_per_m, abs=1e-6)


def test_the_full_spoke_is_asymmetric_by_one_sample(opts, spoke) -> None:
    """
    ``-32*dk ... +31*dk`` at an even matrix, which is ``PhaseEncode``'s ``matrix // 2`` convention
    and the official reference's.  Deliberately not tidied into something symmetric.
    """
    assert spoke.k_first_per_m == pytest.approx(-32 * spoke.dk_per_m)
    assert spoke.k_last_per_m == pytest.approx(+31 * spoke.dk_per_m)


# ----------------------------------------------------------------------- the rotation
@pytest.mark.parametrize('angle_deg', [0.0, 22.5, 90.0, 137.5, -40.0])
def test_the_spoke_points_where_it_was_asked_to(opts, spoke, angle_deg) -> None:
    k = _k(spoke, opts, math.radians(angle_deg))
    measured = math.degrees(math.atan2(k[1, -1] - k[1, 0], k[0, -1] - k[0, 0]))
    assert (measured - angle_deg + 90) % 180 - 90 == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize('angle_deg', [22.5, 45.0, 90.0, 137.5, -40.0])
def test_a_spoke_is_the_spoke_at_zero_rotated(opts, spoke, angle_deg) -> None:
    """
    The metamorphic test, and the defining property of a radial acquisition.

    It is not a restatement of the construction: nothing here knows *how* the module orients a
    spoke, only that the trajectory it produces is a rotation of the one it produces at zero.
    """
    phi = math.radians(angle_deg)
    rotation = np.array([[math.cos(phi), -math.sin(phi)], [math.sin(phi), math.cos(phi)]])

    assert np.abs(rotation @ _k(spoke, opts)[:2] - _k(spoke, opts, phi)[:2]).max() < 1e-3


def test_the_samples_lie_on_a_straight_line_and_stay_in_plane(opts, spoke) -> None:
    """A curved or out-of-plane spoke still grids into a plausible image."""
    k = _k(spoke, opts, 0.7)
    dx, dy = math.cos(0.7), math.sin(0.7)
    # The 2D cross product spelled out: np.cross on 2-vectors is deprecated in NumPy 2.
    assert np.abs(dx * k[1] - dy * k[0]).max() < 1e-3
    assert np.abs(k[2]).max() < 1e-3
    assert np.ptp(np.diff(np.array([dx, dy]) @ k[:2])) < 1e-3


def test_an_axis_aligned_spoke_is_one_gradient_not_two(opts, spoke) -> None:
    """A direction cosine of zero emits nothing, rather than a 1e-11 Hz/m ghost beside the real one."""
    axes = [event.channel for _, event, _ in sc.flatten(spoke(angle_rad=0.0))
            if getattr(event, 'type', '') in {'trap', 'grad'}]
    oblique = [event.channel for _, event, _ in sc.flatten(spoke(angle_rad=math.pi / 4))
               if getattr(event, 'type', '') in {'trap', 'grad'}]

    assert set(axes) == {'x'}
    assert set(oblique) == {'x', 'y'}


def test_a_rotated_gradient_carries_a_rotated_area(opts, spoke) -> None:
    """
    The stored ``area`` is scaled with the amplitude.

    Both radial references that extend a readout's flat time warn in a comment that their event's
    area is then wrong.  Anything asking this module what it plays would read the same defect.
    """
    for _, event, _ in sc.flatten(spoke(angle_rad=math.pi / 3)):
        if getattr(event, 'type', '') == 'trap':
            expected = float(event.amplitude) * (
                event.rise_time / 2 + event.flat_time + event.fall_time / 2)
            assert float(event.area) == pytest.approx(expected, rel=1e-9)


# ---------------------------------------------------------------------- the contract
def test_acquire_false_drops_the_adc_and_nothing_else(opts, spoke) -> None:
    kinds = [getattr(e, 'type', '') for _, e, _ in sc.flatten(spoke(angle_rad=0.3))]
    dummy = [getattr(e, 'type', '') for _, e, _ in sc.flatten(spoke(angle_rad=0.3, acquire=False))]

    assert kinds.count('adc') == 1
    assert dummy.count('adc') == 0
    assert [k for k in kinds if k != 'adc'] == dummy


def test_the_component_contract(opts, spoke, component_checks) -> None:
    """Pure, deterministic, and legal on its own -- including at an oblique angle."""
    component_checks.all(spoke, angle_rad=0.9)


def test_it_exposes_the_prephaser_duration_a_caller_has_to_overlap(opts, spoke) -> None:
    assert spoke.prephaser_duration_s > 0.0
    assert spoke.prephaser_duration_s == pytest.approx(spoke.line.prephaser_duration_s)


@pytest.mark.parametrize('value', [0.0, 0.25, 0.49, 1.01, 2.0])
def test_a_partial_fourier_outside_the_evidence_raises(opts, value) -> None:
    """Closed at both ends: 0.5 is centre-out, which is a family rather than an edge case."""
    with pytest.raises(sc.ConfigurationError, match=r'\[0.5, 1.0\]'):
        sc.modules.RadialReadout(opts=opts, fov_mm=FOV_MM, matrix=MATRIX, dwell_s=DWELL_S,
                                 partial_fourier=value)


def test_the_receiver_phase_reaches_the_adc(opts, spoke) -> None:
    """An RF-spoiling schedule has to move the receiver with the transmitter."""
    adc = next(e for _, e, _ in sc.flatten(spoke(angle_rad=0.2, phase_deg=117.0))
               if getattr(e, 'type', '') == 'adc')
    assert float(adc.phase_offset) == pytest.approx(math.radians(117.0))


def test_an_illegal_sample_count_is_refused_by_the_line_it_is_built_from(opts) -> None:
    """
    Inherited, not duplicated.

    ``CartesianLine`` computes the sample count, so it is the layer that knows whether the
    receiver can digitise it, and a second copy of the rule here would be a second place to keep
    it right.  What this test pins is that composing the line really does carry the refusal --
    the error names the count and the divisor, and arrives at construction rather than from the
    compiler.
    """
    with pytest.raises(sc.ConfigurationError, match='divisible by 4'):
        sc.modules.RadialReadout(opts=opts, fov_mm=FOV_MM, matrix=65, dwell_s=DWELL_S)

    # ...and the same geometry is fine where the count works out.
    assert sc.modules.RadialReadout(opts=opts, fov_mm=FOV_MM, matrix=68,
                                    dwell_s=DWELL_S).num_samples == 68
