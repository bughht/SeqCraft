"""
The peak-B1 contract: SeqCraft does not emit RF over ``opts.max_b1``.

Two layers, and the second is what makes it a contract rather than four habits:

* every RF-producing module checks its own pulse at **design time**, where it still knows what to
  suggest changing;
* :func:`seqcraft.compiler.verification.check_rf_amplitude` measures the **emitted** waveform, so a
  raw pypulseq event, or a future module that forgets its own check, is caught anyway.

The limit is live by default.  ``pp.Opts()`` sets ``max_b1`` to 851.52 Hz -- 20 uT -- so these
tests do not have to opt in to it; ``tests/conftest.py``'s ``unbounded_b1`` is how a test opts
*out*, and it uses ``inf`` rather than zero for the reason
:func:`test_zero_is_not_a_disabling_value` records.
"""

from __future__ import annotations

import copy
import re

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.modules._support import peak_b1_hz

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


def test_infinity_disables_enforcement_cleanly(unbounded_b1) -> None:
    """
    ``inf`` is how a caller designs without a transmit limit, and it needs no sentinel branch.

    The invariant stays literal -- ``peak <= opts.max_b1`` -- and an infinite limit satisfies it.
    Both layers are exercised: the module builds, and the emitted sequence compiles.
    """
    assert unbounded_b1.max_b1 == np.inf
    loud = sc.modules.Excitation(opts=unbounded_b1, flip_deg=90.0, thickness_mm=5.0,
                                 duration_s=1e-3)

    assert _peak(loud) > DEFAULT_MAX_B1_HZ, 'the pulse really is over the default limit'
    sc.compile(sc.LogicBlock('loud').add(0.0, loud()), unbounded_b1)


@pytest.mark.parametrize('limit', [0.0, float('nan'), -np.inf])
def test_only_positive_infinity_means_unlimited(opts, limit: float) -> None:
    """
    Zero, NaN and negative infinity are **not** spellings of "unlimited", at either layer.

    NaN is the dangerous one: every comparison against it is False, so a check written as
    ``worst > limit`` would pass it silently -- which is exactly the failure this contract exists
    to prevent.  Both layers therefore test ``limit > 0.0``, which is False for all three.
    """
    broken = copy.copy(opts)
    broken.max_b1 = limit

    with pytest.raises(sc.ConfigurationError, match='not a transmit limit'):
        sc.modules.Excitation(opts=broken, flip_deg=90.0, thickness_mm=5.0, duration_s=3e-3)

    quiet = pp.make_sinc_pulse(flip_angle=np.pi / 2, duration=4e-3, time_bw_product=4,
                               system=opts, slice_thickness=5e-3, use='excitation',
                               return_gz=True)[0]
    with pytest.raises(sc.HardwareLimitError, match='not a transmit limit'):
        sc.compile(sc.LogicBlock('raw').add(0.0, quiet), broken)


def test_the_backstop_allows_positive_infinity(opts) -> None:
    """The one non-finite value that does disable enforcement, on the compiler path."""
    unlimited = copy.copy(opts)
    unlimited.max_b1 = np.inf
    loud = pp.make_sinc_pulse(flip_angle=np.pi, duration=2e-3, time_bw_product=4, system=opts,
                              slice_thickness=5e-3, use='refocusing', return_gz=True)[0]

    sc.compile(sc.LogicBlock('raw').add(0.0, loud), unlimited)


def test_zero_is_not_a_disabling_value(opts) -> None:
    """
    Why ``unbounded_b1`` uses ``inf``, recorded as a test rather than as a comment.

    ``adc_samples_limit = 0`` is pypulseq's documented "no limit"; ``max_b1`` is **not** the same.
    ``make_sinc_pulse`` compares ``rf_amplitude > system.max_b1`` unconditionally, so at zero it
    warns at ``inf %``.  SeqCraft therefore does not treat zero as disabling either -- it refuses,
    and says the limit itself is the thing to set.
    """
    zeroed = copy.copy(opts)
    zeroed.max_b1 = 0.0

    with pytest.raises(sc.ConfigurationError) as caught:
        sc.modules.Excitation(opts=zeroed, flip_deg=90.0, thickness_mm=5.0, duration_s=3e-3)

    message = str(caught.value)
    assert 'max_b1 is 0' in message
    assert 'max_b1=inf' in message, 'it points at the value that does disable it'


# ------------------------------------------------------------------------- accepted, and refused
@pytest.mark.parametrize('duration_s', [3e-3, 5e-3])
def test_a_pulse_below_the_limit_is_accepted(opts, duration_s: float) -> None:
    module = sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0,
                                   duration_s=duration_s)

    assert _peak(module) < opts.max_b1


def _no_limit(opts):
    relaxed = copy.copy(opts)
    relaxed.max_b1 = np.inf
    return relaxed


def _quoted_duration_ms(message: str) -> float:
    """The duration a refusal quotes, in milliseconds, parsed out of the message it printed."""
    found = re.search(r'duration_s >= ([\d.]+) ms', message)
    assert found, f'no duration quoted in:\n{message}'
    return float(found.group(1))


@pytest.mark.parametrize('pulse', ['sinc', 'gauss'])
def test_a_quoted_floor_really_does_fit(opts, pulse: str) -> None:
    """
    **The rule for quoting a number: rebuild at it and check.**

    A sinc or gauss at a fixed time--bandwidth product is merely stretched by a longer duration,
    so the inverted ``1 / duration`` relation is a floor and the message says so.  This rebuilds
    at exactly the duration the refusal printed, and requires that it passes.
    """
    over = 0.6e-3 if pulse == 'gauss' else 1.0e-3
    with pytest.raises(sc.ConfigurationError) as caught:
        sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0, duration_s=over,
                              pulse=pulse)

    message = str(caught.value)
    assert 'which is where this shape fits' in message, 'a floor is claimed'
    at_the_floor = sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0, pulse=pulse,
                                         duration_s=_quoted_duration_ms(message) / 1e3)

    assert _peak(at_the_floor) <= opts.max_b1
    assert _peak(at_the_floor) > 0.8 * opts.max_b1, 'and it is not wastefully conservative'


def test_gauss_with_an_explicit_bandwidth_does_not_claim_a_floor(opts) -> None:
    """
    The pulse **name** is not enough to know whether lengthening helps.

    ``make_gauss_pulse`` takes either a time--bandwidth product or an explicit ``bandwidth``.  With
    a bandwidth supplied the design is pinned to it, and a longer duration buys almost nothing:
    measured at 4 kHz and 90 degrees, the peak is 1000.0 Hz at 1 ms **and** 1000.0 Hz at 2 ms.

    So this configuration must not be told a floor -- rebuilding at the quoted duration would
    still be over the limit, which is advice that does not work.
    """
    fixed_bandwidth = {'bandwidth': 4000.0}
    with pytest.raises(sc.ConfigurationError) as caught:
        sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0, duration_s=1e-3,
                              pulse='gauss', pulse_opts=fixed_bandwidth)

    message = str(caught.value)
    assert 'which is where this shape fits' not in message, 'no floor may be claimed here'
    assert 'starting point rather than a floor' in message

    # And the reason: rebuilding at the number it printed really would still be over.
    still_over = sc.modules.Excitation(opts=_no_limit(opts), flip_deg=90.0, thickness_mm=5.0,
                                       duration_s=_quoted_duration_ms(message) / 1e3,
                                       pulse='gauss', pulse_opts=fixed_bandwidth)

    assert _peak(still_over) > opts.max_b1


def test_gauss_without_a_bandwidth_does_claim_a_floor(opts) -> None:
    """The same shape, pinned by its time--bandwidth product instead, stretches exactly."""
    with pytest.raises(sc.ConfigurationError) as caught:
        sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0, duration_s=0.6e-3,
                              pulse='gauss')

    assert 'which is where this shape fits' in str(caught.value)


def test_slr_does_not_claim_a_floor(opts) -> None:
    """
    An SLR filter is **recomputed** when the duration changes, so the same arithmetic is a guide.

    The number is still worth printing -- it is right for the default filter -- but the wording
    must not promise what is not guaranteed for every ``filter_type``.
    """
    with pytest.raises(sc.ConfigurationError) as caught:
        sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0, duration_s=1e-3,
                              pulse='slr')

    message = str(caught.value)
    assert 'starting point rather than a floor' in message
    assert 'which is where this shape fits' not in message


def test_saturation_prep_does_not_claim_a_floor(opts) -> None:
    """
    `bandwidth_hz` is part of this module's public contract and is held fixed, so changing the
    duration changes the time--bandwidth product and therefore the design.  No floor is promised.
    """
    with pytest.raises(sc.ConfigurationError) as caught:
        sc.modules.SaturationPrep(opts=opts, shift_ppm=-3.4, flip_deg=110.0, duration_s=0.3e-3,
                                  bandwidth_hz=200.0, spoil_voxel_mm=5.0)

    assert 'starting point rather than a floor' in str(caught.value)


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


def test_wurst_remedy_names_bandwidth_because_bandwidth_works(opts) -> None:
    """
    WURST's peak goes as ``sqrt(bandwidth)``: 2523 Hz at 40 kHz and 564 at 2 kHz.

    A factor of 4.47 for a factor of 20, which is ``sqrt(20)`` -- so the remedy is real but the
    message must not promise a linear one.
    """
    with pytest.raises(sc.ConfigurationError, match='max_b1') as caught:
        sc.modules.IRPrep(opts=opts, thickness_mm=None, pulse='wurst', duration_s=4e-3,
                          spoil_voxel_mm=5.0)

    message = str(caught.value)
    assert 'bandwidth' in message and 'adiabaticity' in message
    assert 'square root' in message, 'the scaling is not claimed to be linear'
    narrowed = sc.modules.IRPrep(opts=opts, thickness_mm=None, pulse='wurst', duration_s=4e-3,
                                 spoil_voxel_mm=5.0, pulse_opts={'bandwidth': 2000})

    assert _peak(narrowed) <= opts.max_b1, 'the remedy it names actually works'


def test_hypsec_remedy_does_not_name_bandwidth_because_bandwidth_does_nothing(opts) -> None:
    """
    The reason the two adiabatic shapes do not share a message.

    ``make_adiabatic_pulse`` ignores ``bandwidth`` for a hyperbolic secant -- 563.7 Hz at 40, 10
    and 2 kHz alike -- so recommending it would be advice that cannot work.  ``beta`` and ``mu``
    are what set the peak here.
    """
    loud = {'beta': 3000.0}
    with pytest.raises(sc.ConfigurationError, match='max_b1') as caught:
        sc.modules.IRPrep(opts=opts, thickness_mm=None, pulse='hypsec', duration_s=10e-3,
                          spoil_voxel_mm=5.0, pulse_opts=loud)

    message = str(caught.value)
    assert 'beta' in message and 'mu' in message
    assert 'bandwidth is not a remedy' in message
    assert 'adiabaticity' in message and 'robustness' in message, 'and it is not free'


def test_bandwidth_really_does_nothing_to_a_hyperbolic_secant(unbounded_b1) -> None:
    """The measurement the previous test's message rests on."""
    peaks = {
        bw: _peak(sc.modules.IRPrep(opts=unbounded_b1, thickness_mm=None, pulse='hypsec',
                                    duration_s=10e-3, spoil_voxel_mm=5.0,
                                    pulse_opts={'bandwidth': bw}))
        for bw in (40000, 10000, 2000)
    }

    assert len({round(v, 6) for v in peaks.values()}) == 1


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


# ------------------------------------------------------------------------------ other nuclei
def test_the_reported_microtesla_follow_opts_gamma(opts) -> None:
    """
    SeqCraft stores and compares ``max_b1`` in hertz; the microtesla it prints follow ``gamma``.

    The hertz value is what the caller configured, so it is unchanged by the nucleus.  The **same
    field strength** in microtesla is a different number of hertz for a different nucleus, which
    is exactly why the conversion cannot use a hard-coded proton gamma: the µT printed beside a
    carbon limit would otherwise be wrong by the ratio of the two.
    """
    carbon = copy.copy(opts)
    carbon.gamma = 10.7084e6                          # 13C
    assert carbon.max_b1 == opts.max_b1, 'the configured limit is in hertz and is not converted'

    with pytest.raises(sc.ConfigurationError) as proton:
        sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0, duration_s=1e-3)
    with pytest.raises(sc.ConfigurationError) as other:
        sc.modules.Excitation(opts=carbon, flip_deg=90.0, thickness_mm=5.0, duration_s=1e-3)

    # The configured limit and the measured peak are the same hertz on both sides ...
    assert 'peak_b1_hz:  1107.63' in str(proton.value)
    assert 'peak_b1_hz:  1107.63' in str(other.value)
    # ... and the microtesla reported for them follow the nucleus.
    assert 'max_b1_uT :  20.0' in str(proton.value)
    assert 'max_b1_uT :  79.52' in str(other.value)
