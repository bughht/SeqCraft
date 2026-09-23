"""
The peak-B1 contract: SeqCraft does not emit RF over ``opts.max_b1``.

Two layers, and the second is what makes it a contract rather than four habits:

* every RF-producing module checks its own pulse at **design time**, where it still knows what to
  suggest changing;
* :func:`seqcraft.compiler.verification.check_rf_amplitude` measures the **emitted** waveform, so a
  raw pypulseq event, or a future module that forgets its own check, is caught anyway.

The limit is live by default.  ``pp.Opts()`` sets ``max_b1`` to 851.52 Hz -- 20 uT -- so these
tests do not have to opt in to it; ``tests/conftest.py``'s ``unbounded_b1`` is how a test opts
*out*, and why.
"""

from __future__ import annotations

import copy

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.modules._support import duration_for_peak_b1, peak_b1_hz

#: Pulses are quoted against this throughout, because it is what a bare ``pp.Opts()`` carries.
DEFAULT_MAX_B1_HZ = 851.52


def _peak(module) -> float:
    return peak_b1_hz(module.rf)


# --------------------------------------------------------------- the limit, and what disables it
def test_the_pypulseq_default_is_a_real_limit(opts) -> None:
    """
    The fact that makes this contract worth enforcing rather than opt-in.

    ``pp.Opts()`` does not leave ``max_b1`` unset: it defaults to 20 uT, which an ordinary 1 ms
    90 degree sinc exceeds.  A user who never thinks about B1 still has a limit, so a module that
    never checks one is silently emitting pulses that limit forbids.
    """
    assert opts.max_b1 == pytest.approx(DEFAULT_MAX_B1_HZ)
    assert sc.convert(opts.max_b1, 'Hz', 'uT') == pytest.approx(20.0, rel=1e-3)


def test_zero_disables_the_check_everywhere(unbounded_b1) -> None:
    """``max_b1 = 0`` is pypulseq's "no limit", the same convention ``adc_samples_limit`` uses."""
    assert unbounded_b1.max_b1 == 0.0
    loud = sc.modules.Excitation(opts=unbounded_b1, flip_deg=90.0, thickness_mm=5.0,
                                 duration_s=1e-3)

    assert _peak(loud) > DEFAULT_MAX_B1_HZ, 'the pulse really is over the default limit'
    sc.compile(sc.LogicBlock('loud').add(0.0, loud()), unbounded_b1)


# ------------------------------------------------------------------------- accepted, and refused
@pytest.mark.parametrize('duration_s', [3e-3, 5e-3])
def test_a_pulse_below_the_limit_is_accepted(opts, duration_s: float) -> None:
    module = sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0,
                                   duration_s=duration_s)

    assert _peak(module) < opts.max_b1


def test_a_pulse_just_below_the_limit_is_accepted(opts) -> None:
    """
    The boundary, approached from the side that must work.

    ``duration_for_peak_b1`` rounds up onto a tenth of a millisecond, so the duration it returns
    lands just under the limit rather than exactly on it -- and *that* is the number a refusal
    quotes, so it has to be accepted when a caller uses it.
    """
    over = 1.0e-3
    peak = peak_b1_hz(sc.modules.Excitation(opts=copy.copy(_no_limit(opts)), flip_deg=90.0,
                                            thickness_mm=5.0, duration_s=over).rf)
    floor_s = duration_for_peak_b1(over, peak, opts.max_b1)

    module = sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0,
                                   duration_s=floor_s)

    assert _peak(module) <= opts.max_b1
    assert _peak(module) > 0.9 * opts.max_b1, 'the quoted floor is not wastefully conservative'


def _no_limit(opts):
    relaxed = copy.copy(opts)
    relaxed.max_b1 = 0.0
    return relaxed


def test_excitation_over_the_limit_is_refused(opts) -> None:
    """The case that motivated this contract: an ordinary 1 ms 90, previously emitted silently."""
    with pytest.raises(sc.ConfigurationError, match='max_b1') as caught:
        sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0, duration_s=1e-3)

    assert 'duration_s' in str(caught.value), 'the remedy names what to change'


def test_refocusing_over_the_limit_is_refused(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='max_b1') as caught:
        sc.modules.Refocusing(opts=opts, thickness_mm=5.0, duration_s=2e-3)

    assert 'duration_s' in str(caught.value)


def test_the_two_rf_leaves_refuse_the_same_excess_consistently(opts) -> None:
    """
    A 180 needs exactly twice a 90's peak at a fixed shape, so the same excess is reached at twice
    the duration.  Both leaves must refuse it, and before this contract only one did.
    """
    for make in (
        lambda: sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0,
                                      duration_s=1e-3),
        lambda: sc.modules.Refocusing(opts=opts, thickness_mm=5.0, duration_s=2e-3),
    ):
        with pytest.raises(sc.ConfigurationError, match='max_b1'):
            make()


# ------------------------------------------------------------------------ the preparation layer
def test_ir_prep_shaped_over_the_limit_is_refused(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='max_b1') as caught:
        sc.modules.IRPrep(opts=opts, thickness_mm=None, pulse='sinc', duration_s=2e-3,
                          spoil_voxel_mm=5.0)

    assert 'duration_s' in str(caught.value), 'a shaped pulse scales as 1 / duration'


def test_ir_prep_adiabatic_is_refused_with_the_remedy_it_actually_has(opts) -> None:
    """
    The reason the remedy is the module's and not the shared helper's.

    An adiabatic pulse's peak is set by its frequency sweep, not by its length: a hyperbolic
    secant's peak does not move with duration at all, and WURST's falls only as its square root.
    Quoting "make it longer" here would be advice that does not work.
    """
    with pytest.raises(sc.ConfigurationError, match='max_b1') as caught:
        sc.modules.IRPrep(opts=opts, thickness_mm=None, pulse='wurst', duration_s=4e-3,
                          spoil_voxel_mm=5.0)

    message = str(caught.value)
    assert 'bandwidth' in message and 'adiabaticity' in message
    assert 'does not help much' in message


def test_the_hyperbolic_secant_default_fits_and_does_not_move_with_duration(opts) -> None:
    """The default inversion is playable, and the claim about adiabatic pulses is measured."""
    peaks = {
        d: _peak(sc.modules.IRPrep(opts=opts, thickness_mm=None, pulse='hypsec', duration_s=d,
                                   spoil_voxel_mm=5.0))
        for d in (5e-3, 10e-3, 20e-3)
    }

    assert max(peaks.values()) < opts.max_b1
    assert len(set(round(v, 6) for v in peaks.values())) == 1, 'unchanged across a factor of four'


def test_saturation_prep_over_the_limit_is_refused(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='max_b1'):
        sc.modules.SaturationPrep(opts=opts, shift_ppm=-3.4, flip_deg=110.0, duration_s=0.15e-3,
                                  bandwidth_hz=200.0, spoil_voxel_mm=5.0)


def test_a_realistic_fat_saturation_is_far_below_the_limit(opts) -> None:
    """The protocol `examples/fat_sat/` uses, so the check is not in a working caller's way."""
    module = sc.modules.SaturationPrep(opts=opts, shift_ppm=-3.4, flip_deg=110.0, duration_s=8e-3,
                                       bandwidth_hz=200.0, spoil_voxel_mm=5.0)

    assert _peak(module) < 0.2 * opts.max_b1


# ------------------------------------------------------------------------------- the backstop
def test_a_raw_rf_event_is_caught_by_the_compiler(opts) -> None:
    """
    The backstop, on the path no module guards: a pypulseq event added straight to a tree.

    This is what makes the contract hold for a future RF module that forgets its own check.
    """
    rf = pp.make_sinc_pulse(flip_angle=np.pi, duration=2e-3, time_bw_product=4, system=opts,
                            slice_thickness=5e-3, use='refocusing', return_gz=True)[0]

    with pytest.raises(sc.HardwareLimitError, match='max_b1') as caught:
        sc.compile(sc.LogicBlock('raw').add(0.0, rf), opts)

    assert 'block 1' in str(caught.value), 'the failure names where it came from'


def test_the_backstop_reports_the_worst_event_not_the_first(opts) -> None:
    quiet = pp.make_sinc_pulse(flip_angle=np.pi / 2, duration=4e-3, time_bw_product=4, system=opts,
                               slice_thickness=5e-3, use='excitation', return_gz=True)[0]
    loud = pp.make_sinc_pulse(flip_angle=np.pi, duration=2e-3, time_bw_product=4, system=opts,
                              slice_thickness=5e-3, use='refocusing', return_gz=True)[0]
    tree = sc.LogicBlock('two').add(0.0, quiet).add(10e-3, loud)

    with pytest.raises(sc.HardwareLimitError) as caught:
        sc.compile(tree, opts)

    assert f'{peak_b1_hz(loud):.2f}' in str(caught.value)


def test_a_raw_rf_event_under_the_limit_compiles(opts) -> None:
    rf = pp.make_sinc_pulse(flip_angle=np.pi / 2, duration=4e-3, time_bw_product=4, system=opts,
                            slice_thickness=5e-3, use='excitation', return_gz=True)[0]

    sc.compile(sc.LogicBlock('quiet').add(0.0, rf), opts)
