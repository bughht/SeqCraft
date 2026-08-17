"""
``CartesianLine`` -- where k = 0 actually lands, measured rather than asserted about.

The module carries the hardest arithmetic in the library, so it is tested against pypulseq's own
k-space calculation rather than against the formula that produced it: ``calculate_kspacePP``
integrates the compiled waveform at the true ADC sample times, which is an independent oracle.

The ramp-heavy case is the one that matters.  ``-gx.flat_area / 2`` and ``-gx.area / 2`` differ
by the ramp area alone, so on a fast-slewing system the wrong formula hides inside tolerance;
against a deliberately slow slew rate it is off by whole ``dk`` steps.
"""

from __future__ import annotations

import numpy as np
import pytest
from pypulseq.opts import Opts

import seqcraft as sc
from seqcraft.modules._support import area_until


@pytest.fixture(scope='module')
def slow_opts() -> Opts:
    """
    A deliberately ramp-heavy scanner: the same amplitude, a sixth of the slew.

    The readout's ramps are then a large fraction of its total area, which is exactly where the
    two candidate prephaser formulas separate.
    """
    return Opts(
        max_grad=40, grad_unit='mT/m', max_slew=25, slew_unit='T/m/s', B0=3.0,
        rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
    )


def _k_at_samples(readout: sc.modules.CartesianLine, opts: Opts) -> np.ndarray:
    """k along the readout axis at every ADC sample, from pypulseq's own integration."""
    axis = 'xyz'.index(readout.axis)
    return sc.kspace(sc.LogicBlock('probe').add(0.0, readout()), opts)['k_adc'][axis]


# ------------------------------------------------------------------- the prephaser
#: A scanner and a bandwidth: the second pair makes the ramps a third of the readout, which is
#: where the two candidate prephaser formulas separate by 31 ``dk`` rather than by half of one.
SCANNERS = [('opts', 250.0), ('slow_opts', 2000.0)]


@pytest.mark.parametrize(('scanner', 'bandwidth'), SCANNERS)
def test_k_is_zero_at_the_echo_sample(request, scanner, bandwidth) -> None:
    """
    The single assertion the module exists to satisfy, on a fast and a ramp-heavy system.

    ``-gx.flat_area / 2`` passes neither, and misses by the ramp-up area.
    """
    opts = request.getfixturevalue(scanner)
    ro = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=bandwidth)

    k = _k_at_samples(ro, opts)

    assert k[ro.pre_echo_samples] == pytest.approx(0.0, abs=1e-4 * ro.dk_per_m)


@pytest.mark.parametrize(('scanner', 'bandwidth'), SCANNERS)
def test_the_samples_span_the_intended_k_range(request, scanner, bandwidth) -> None:
    """Adjacent samples are one ``dk`` apart, and the extent is ``matrix * dk``."""
    opts = request.getfixturevalue(scanner)
    ro = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=bandwidth)

    k = _k_at_samples(ro, opts)

    assert np.allclose(np.diff(k), ro.dk_per_m, rtol=1e-6)
    assert k[0] == pytest.approx(-32 * ro.dk_per_m, rel=1e-6)


def test_the_prephaser_is_not_half_the_flat_area(slow_opts) -> None:
    """
    The classic error, named so a future reader sees what the number is *not*.

    Here the ramps are a third of the readout, so the two differ by 31 ``dk`` -- the echo would
    land a third of the way through the ADC instead of at its centre, and the image would show
    ringing and signal loss that read as a hardware fault.  With a fast ramp the same error is
    half a ``dk`` and hides inside any tolerance loose enough to pass.
    """
    ro = sc.modules.CartesianLine(opts=slow_opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=2000.0)

    wrong = -float(ro.gx.flat_area) / 2

    assert abs(ro.prephaser_area_per_m - wrong) > 30 * ro.dk_per_m


def test_the_closed_form_is_right_only_up_to_the_half_dwell_offset(opts) -> None:
    """
    ``-gx.area / 2`` is the *nearly* right answer, and the residue is worth naming.

    k = 0 belongs at sample index ``matrix // 2``, whose centre time is half a dwell after the
    midpoint of the ADC window -- that is the convention a centred ``ifft`` assumes, and it is
    what makes the sampled positions land on ``-N/2 ... N/2-1`` times ``dk`` rather than on a
    half-integer grid.  So the exact prephaser is half the total area *plus half a* ``dk``, and
    a closed form derived from the trapezoid's geometry alone cannot know that.

    Integrating the gradient's own knots up to the echo produces it for free, which is the whole
    argument for computing the prephaser that way rather than from a formula per case.
    """
    ro = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0)

    assert float(ro.gx.flat_time) == pytest.approx(ro.num_samples * ro.dwell_s), (
        'this scanner takes the whole ADC window on the flat top with nothing left over, so '
        'the half dwell is the only thing between the two numbers'
    )
    assert ro.prephaser_area_per_m == pytest.approx(
        -float(ro.gx.area) / 2 - 0.5 * ro.dk_per_m, abs=1e-9,
    )


def test_the_prephaser_cancels_the_measured_moment(opts) -> None:
    """Stated as the identity rather than as a formula, which is why the property is public."""
    ro = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=128, bandwidth_hz_px=200.0)

    accumulated = area_until(ro.gx, ro.time_to_echo() - ro.prephaser_duration_s)

    assert ro.prephaser_area_per_m == pytest.approx(-accumulated, abs=1e-9)


# ------------------------------------------------------------------- partial echo
def test_partial_fourier_one_is_the_symmetric_case_exactly(opts) -> None:
    """A factor that reduces to the existing constant adds no branch, and this proves it."""
    full = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0)
    explicit = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64,
                                        bandwidth_hz_px=250.0, partial_fourier=1.0)

    assert explicit.num_samples == full.num_samples
    assert explicit.prephaser_area_per_m == full.prephaser_area_per_m
    assert explicit.time_to_echo() == full.time_to_echo()


def test_partial_fourier_three_quarters_puts_the_echo_a_third_in(opts) -> None:
    """``(2*pf - 1) / (2*pf)`` is 1/3 at 0.75 -- and it is the *pre-echo* side that is dropped."""
    pf = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0,
                                  partial_fourier=0.75)

    assert pf.num_samples == 48
    assert pf.pre_echo_samples == 16
    assert pf.pre_echo_samples / pf.num_samples == pytest.approx(1 / 3)

    k = _k_at_samples(pf, opts)
    assert k[pf.pre_echo_samples] == pytest.approx(0.0, abs=1e-4 * pf.dk_per_m)
    assert k[-1] == pytest.approx(31 * pf.dk_per_m, rel=1e-6), 'the post-echo side stays full'


def test_time_to_echo_tracks_partial_fourier(opts) -> None:
    """A shorter pre-echo side is a shorter prephaser and an earlier echo, and TE follows."""
    full = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0)
    pf = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0,
                                  partial_fourier=0.75)

    assert pf.time_to_echo() < full.time_to_echo()
    assert abs(pf.prephaser_area_per_m) < abs(full.prephaser_area_per_m)


@pytest.mark.parametrize('value', [0.0, -0.5, 1.01, 2.0])
def test_partial_fourier_out_of_range_raises(opts, value) -> None:
    with pytest.raises(sc.ConfigurationError, match='partial_fourier'):
        sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64,
                                 bandwidth_hz_px=250.0, partial_fourier=value)


# ------------------------------------------------------------------ the sampling rate
def test_bandwidth_and_dwell_are_one_number_given_either_way(opts) -> None:
    by_bandwidth = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64,
                                            bandwidth_hz_px=250.0)
    by_dwell = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64,
                                        dwell_s=by_bandwidth.dwell_s)

    assert by_dwell.dwell_s == by_bandwidth.dwell_s
    assert by_dwell.bandwidth_hz_px == pytest.approx(by_bandwidth.bandwidth_hz_px)


@pytest.mark.parametrize('kwargs', [
    {},
    {'bandwidth_hz_px': 250.0, 'dwell_s': 6e-6},
])
def test_both_or_neither_raises(opts, kwargs) -> None:
    """Two ways to say one thing is fine; saying it twice is a question about which one wins."""
    with pytest.raises(sc.ConfigurationError, match='exactly one'):
        sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, **kwargs)


def test_a_bandwidth_the_amplifier_cannot_reach_raises_with_the_one_that_works(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='Hz/m'):
        sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=20000.0)


# --------------------------------------------------------------------- the two builds
def test_acquire_false_drops_the_adc_and_nothing_else(opts) -> None:
    """
    A dummy has to load the gradients exactly as a real repetition does.

    If it does not, the steady state it establishes is not the one that gets acquired -- and
    nothing about the file looks wrong.
    """
    ro = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0)

    real, dummy = ro(), ro(acquire=False)

    assert dummy.duration == real.duration
    assert [n.item.type for n in dummy] == ['trap', 'trap']
    assert sorted(n.item.type for n in real) == ['adc', 'trap', 'trap']


def test_the_receiver_phase_follows_the_transmitter(opts) -> None:
    """
    The fix the simulation forced, and the reason ``phase_deg`` is on this signature at all.

    The receiver is phase-locked to the transmitter.  With an RF-spoiling schedule on the pulse
    and a receiver left at zero, the schedule's quadratic phase lands in ``ky``: a single
    off-centre voxel came back scattered across the whole phase-encode direction with its peak
    thirteen pixels out, while the readout direction stayed perfectly correct -- which is what
    makes it look like a phase-encode bug rather than a demodulation one.
    """
    ro = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0)

    adc = next(n.item for n in ro(phase_deg=117.0) if getattr(n.item, 'type', '') == 'adc')

    assert adc.phase_offset == pytest.approx(np.deg2rad(117.0))
    assert ro.adc.phase_offset == 0.0, 'the design was not mutated'


def test_the_receiver_phase_and_the_fov_offset_compose(opts) -> None:
    """Two independent terms on one field, and neither may quietly replace the other."""
    ro = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0)

    both = next(n.item for n in ro(phase_deg=117.0, offset_mm=30.0)
                if getattr(n.item, 'type', '') == 'adc')
    shift_only = next(n.item for n in ro(offset_mm=30.0)
                      if getattr(n.item, 'type', '') == 'adc')

    assert both.freq_offset == pytest.approx(shift_only.freq_offset)
    assert both.phase_offset == pytest.approx(shift_only.phase_offset + np.deg2rad(117.0))


def test_an_offset_fov_demodulates_at_a_shifted_frequency(opts) -> None:
    """
    ``adc.freq_offset = gx.amplitude * offset_m``, with the phase referenced to the echo.

    The phase term is the one that gets omitted, and omitting it is a linear phase across the
    image that reads as a gradient-delay problem rather than as a missing reference.
    """
    ro = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0)

    adc = next(n.item for n in ro(offset_mm=30.0) if getattr(n.item, 'type', '') == 'adc')

    expected_freq = float(ro.gx.amplitude) * 0.030
    assert adc.freq_offset == pytest.approx(expected_freq)
    t_ref = ro.time_to_echo() - ro.prephaser_duration_s - float(adc.delay)
    assert adc.phase_offset == pytest.approx(-2 * np.pi * expected_freq * t_ref)
    assert ro().nodes[-1].item.phase_offset == 0.0, 'the design itself was not touched'


def test_the_offset_is_symmetric_about_isocentre(opts) -> None:
    """A sign error mirrors the FOV, which against a symmetric phantom looks entirely correct."""
    ro = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0)

    plus = next(n.item for n in ro(offset_mm=20.0) if getattr(n.item, 'type', '') == 'adc')
    minus = next(n.item for n in ro(offset_mm=-20.0) if getattr(n.item, 'type', '') == 'adc')

    assert plus.freq_offset == pytest.approx(-minus.freq_offset)
    assert plus.phase_offset == pytest.approx(-minus.phase_offset)


def test_the_adc_delay_carries_the_offset_not_the_node(opts) -> None:
    """
    ``pp.make_adc`` raises a delay below ``adc_dead_time`` up to it, and seqcraft preserves an
    event's own delay -- so an offset written into the node would be added twice.
    """
    ro = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0)

    block = ro()
    # Two events share the x channel here -- the prephaser and the readout -- so the readout is
    # the one that starts where the prephaser ends, not simply the first x gradient found.
    gx_node = next(n for n in block if n.item is ro.gx)
    adc_node = next(n for n in block if getattr(n.item, 'type', '') == 'adc')

    assert gx_node.start == pytest.approx(ro.prephaser_duration_s)
    assert adc_node.start == gx_node.start
    assert float(adc_node.item.delay) == pytest.approx(float(ro.gx.rise_time))


# --------------------------------------------------------------------- the winder
def test_a_stretched_prephaser_keeps_the_same_area(opts) -> None:
    """Stretching the shorter of the two winders keeps TE at its minimum, which a delay would not."""
    short = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0)
    stretched = sc.modules.CartesianLine(
        opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0,
        prephaser_duration_s=short.prephaser_duration_s + 200e-6,
    )

    assert stretched.prephaser_duration_s == pytest.approx(short.prephaser_duration_s + 200e-6)
    assert stretched.prephaser_area_per_m == pytest.approx(short.prephaser_area_per_m, rel=1e-9)
    assert stretched.time_to_echo() == pytest.approx(short.time_to_echo() + 200e-6)


def test_a_prephaser_shorter_than_its_minimum_raises(opts) -> None:
    short = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0)

    with pytest.raises(sc.ConfigurationError, match='min_prephaser_duration_s'):
        sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0,
                                 prephaser_duration_s=short.prephaser_duration_s / 2)


def test_it_is_pure_and_compiles_alone(opts, component_checks) -> None:
    component_checks.all(
        sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0),
    )


def test_the_line_index_is_not_an_argument(opts) -> None:
    """
    A Cartesian line is identical every TR; which line it encodes belongs to ``PhaseEncode``.

    Stated as a test because the tempting API is the wrong one: putting ``line=`` here would
    make the readout depend on the phase-encode convention, and would put the ``LIN`` label in
    a module that does not know which line it is.
    """
    ro = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0)

    with pytest.raises(TypeError, match='line'):
        ro(line=3)
