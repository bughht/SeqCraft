"""
``TSEShot`` -- the three things a train has that a single refocused echo does not.

Everything here is a property of the *coupling* rather than of any leaf, which is the reason the
module exists.  Each is measured against pypulseq's own k-space integration or against the
compiled sequence, because those are the two independent oracles available -- and each has the
same failure mode, which is that nothing downstream complains.  A wrong crusher, a drifting pulse
spacing and an echo off the midpoint are all legal block structures that meet every limit; what
they cost is signal, and against a symmetric phantom the image still looks plausible.
"""

from __future__ import annotations

import numpy as np
import pytest

import seqcraft as sc

MATRIX = (32, 16)
ECHOES = 4


@pytest.fixture(scope='module')
def shot(opts):
    """A small but complete turbo-spin-echo shot."""
    return sc.modules.TSEShot(opts=opts, fov_mm=200.0, matrix=MATRIX, thickness_mm=5.0,
                              echoes=ECHOES)


def _lines() -> list[int]:
    return [2, 7, 8, 13]


def _kspace(shot, lines):
    return sc.kspace(shot(lines=lines), shot.opts)


# ------------------------------------------------------------------ k at every echo
def test_k_is_exact_at_every_echo(opts, shot) -> None:
    """
    Signed, on all three axes.  ``|k|`` is symmetric, so an extent check passes on a mirrored
    echo and a sign error survives it; this does not.
    """
    lines = _lines()
    k = _kspace(shot, lines)
    at_echo = k['k_adc'].reshape(3, ECHOES, shot.ro.num_samples)[:, :, shot.ro.echo_sample(0)]
    want_y = (np.asarray(lines) - shot.center_line) * (1e3 / 200.0)

    assert np.abs(at_echo[0]).max() < 1e-3
    assert np.abs(at_echo[2]).max() < 1e-3
    assert np.abs(at_echo[1] - want_y).max() < 1e-3


def test_ky_returns_to_zero_at_every_refocusing_pulse(opts, shot) -> None:
    """
    The blip goes after the 180 and the rewinder before the next one.

    The alternative -- blip first, which a single-echo spin echo may legally do -- leaves ``ky``
    non-zero at the pulse, and a refocusing pulse conjugates whatever it finds there, so every
    later echo carries the wrong line by a growing amount.
    """
    k = _kspace(shot, _lines())
    at_pulse = np.interp(k['t_refocusing'], k['t_k'], k['k'][1])
    assert np.abs(at_pulse).max() < 1e-3


def test_the_refocusing_pulses_are_uniformly_spaced(opts, shot) -> None:
    """An FSE signal is a mixture of primary and stimulated echoes; drift dephases the mixture."""
    spacing = np.diff(_kspace(shot, _lines())['t_refocusing'])
    assert np.ptp(spacing) < 1e-9
    assert spacing[0] == pytest.approx(shot.echo_spacing_s, abs=1e-9)


def test_every_echo_sits_at_the_midpoint_between_its_two_pulses(opts, shot) -> None:
    """
    Within one gradient raster, and not zero: the echo instant lives on the ADC raster while
    gradient starts live on the gradient raster, so the two cannot coincide in general.
    """
    k = _kspace(shot, _lines())
    centres = np.asarray(k['t_refocusing'])
    echoes = np.asarray(k['t_adc']).reshape(ECHOES, shot.ro.num_samples)[
        :, shot.ro.echo_sample(0)]
    midpoints = (centres[:-1] + centres[1:]) / 2
    assert np.abs(echoes[:-1] - midpoints).max() <= float(opts.grad_raster_time)


def test_the_first_interval_is_half_the_echo_spacing(opts, shot) -> None:
    """Measured from the excitation's *effective* centre, which is not its midpoint."""
    k = _kspace(shot, _lines())
    first = float(k['t_refocusing'][0]) - float(k['t_excitation'][0])
    assert first == pytest.approx(shot.echo_spacing_s / 2, abs=1e-9)


# --------------------------------------------------------------- what the kernel couples
def test_the_two_readout_lobes_differ_by_the_readouts_own_imbalance(opts, shot) -> None:
    """
    The one number this module takes from ``CartesianLine`` rather than deriving.

    The leaf states how asymmetric it is about its echo; the kernel decides that the difference
    is split evenly between the lobe before the pulse and the lobe after it.  Deriving the
    imbalance here too would be the duplication the property exists to remove.
    """
    difference = float(shot.lobe_in.area) - float(shot.lobe_out.area)
    assert difference == pytest.approx(shot.ro.echo_moment_imbalance_per_m, abs=1e-12)
    assert shot.skew_per_m == pytest.approx(shot.ro.echo_moment_imbalance_per_m / 2, abs=1e-12)


def test_the_crusher_window_is_the_longest_of_three(opts, shot) -> None:
    """Slice crusher on z, compensating lobes on x, phase blip on y: one window, three axes."""
    assert set(shot.crush_min_s) == {'x', 'y', 'z'}
    assert shot.crush_s >= max(shot.crush_min_s.values()) - 1e-9
    assert shot.pe(line=0).duration == pytest.approx(shot.crush_s)


def test_a_requested_echo_spacing_is_delivered(opts) -> None:
    """Widening the window by ``d`` lengthens the spacing by ``2d``, so the request inverts."""
    shortest = sc.modules.TSEShot(opts=opts, fov_mm=200.0, matrix=MATRIX, thickness_mm=5.0,
                                  echoes=ECHOES)
    wanted = shortest.echo_spacing_s + 2e-3
    longer = sc.modules.TSEShot(opts=opts, fov_mm=200.0, matrix=MATRIX, thickness_mm=5.0,
                                echoes=ECHOES, echo_spacing_s=wanted)
    assert longer.echo_spacing_s >= wanted - 1e-9
    assert longer.echo_spacing_s <= wanted + float(opts.grad_raster_time) * 2
    # A request below the minimum is not an error and does not shorten anything.
    assert sc.modules.TSEShot(opts=opts, fov_mm=200.0, matrix=MATRIX, thickness_mm=5.0,
                              echoes=ECHOES, echo_spacing_s=1e-4).echo_spacing_s == pytest.approx(
        shortest.echo_spacing_s)


def test_time_to_echo_matches_the_sequence_it_builds(opts, shot) -> None:
    """The claim against the measurement -- the failure no k-space check can see."""
    k = _kspace(shot, _lines())
    measured = np.asarray(k['t_adc']).reshape(ECHOES, shot.ro.num_samples)[
        :, shot.ro.echo_sample(0)]
    claimed = [shot.time_to_echo(n) for n in range(ECHOES)]
    assert np.abs(measured - claimed).max() < 1e-9
    assert shot.te_s(0) == pytest.approx(shot.time_to_echo(0) - shot.exc.time_to_center())


def test_the_block_is_exactly_tr_long(opts, shot) -> None:
    """A caller stacks shots at ``tr_s``, so the block has to fill it."""
    assert shot(lines=_lines()).duration == pytest.approx(shot.tr_s, abs=1e-12)
    assert shot.tr_s >= shot.shot_s - 1e-12


def test_a_dummy_shot_plays_the_same_gradients(opts, shot) -> None:
    """``acquire=False`` drops the ADC and the labels and nothing else."""
    kinds = [getattr(event, 'type', '') for _, event, _ in sc.flatten(shot(lines=_lines()))]
    dummy = [getattr(event, 'type', '')
             for _, event, _ in sc.flatten(shot(lines=_lines(), acquire=False))]
    assert kinds.count('adc') == ECHOES
    assert dummy.count('adc') == 0
    assert [k for k in kinds if k not in {'adc', 'labelset'}] == [
        k for k in dummy if k not in {'adc', 'labelset'}]


def test_the_shot_places_the_segment_label_the_caller_values(opts, shot) -> None:
    """The value is the acquisition's; where it goes is the shot's, because it knows."""
    labels = [(t, e.label, e.value) for t, e, _ in sc.flatten(shot(lines=_lines(), segment_index=3))
              if getattr(e, 'type', '') == 'labelset']
    assert [v for _, name, v in labels if name == 'SEG'] == [3]
    assert [v for _, name, v in labels if name == 'LIN'] == _lines()
    assert not [name for _, name, _ in labels if name == 'ECO']
    assert not [name for _, _, name in
                [(t, e, getattr(e, 'label', '')) for t, e, _ in
                 sc.flatten(shot(lines=_lines()))] if name == 'SEG']


def test_the_component_contract(opts, shot, component_checks) -> None:
    """Pure, deterministic, and legal on its own."""
    component_checks.all(shot, lines=_lines())


# ------------------------------------------------------------------------ the refusals
def test_a_tr_shorter_than_the_train_raises(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='shorter than'):
        sc.modules.TSEShot(opts=opts, fov_mm=200.0, matrix=MATRIX, thickness_mm=5.0,
                           echoes=ECHOES, tr_s=1e-3)


def test_no_echoes_raises(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='not a train'):
        sc.modules.TSEShot(opts=opts, fov_mm=200.0, matrix=MATRIX, thickness_mm=5.0, echoes=0)


def test_the_wrong_number_of_lines_raises(opts, shot) -> None:
    """Every shot plays ``echoes`` pulses whether or not they sample."""
    with pytest.raises(sc.ConfigurationError, match='refocusing pulses'):
        shot(lines=[0, 1])


def test_a_line_outside_the_matrix_raises(opts, shot) -> None:
    with pytest.raises(sc.ConfigurationError, match='outside'):
        shot(lines=[0, 1, 2, MATRIX[1]])


def test_a_phase_encode_offset_raises(opts, shot) -> None:
    with pytest.raises(sc.ConfigurationError, match='per-line ADC'):
        shot(lines=_lines(), center_mm=(0.0, 5.0, 0.0))


def test_an_echo_outside_the_train_raises(opts, shot) -> None:
    with pytest.raises(sc.ConfigurationError, match='outside'):
        shot.time_to_echo(ECHOES)
