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

import warnings

import numpy as np
import pypulseq as pp
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


# ------------------------------------------------------------------ prephase=False
#: The whole argument for the boolean is that **nothing else moves**.  Every designed number is
#: listed here rather than spot-checked, because "nothing else moves" is the claim.
DESIGNED = ('dwell_s', 'bandwidth_hz_px', 'num_samples', 'pre_echo_samples', 'dk_per_m')


@pytest.mark.parametrize('partial_fourier', [1.0, 0.75, 0.625])
def test_prephase_false_changes_no_designed_number(opts, partial_fourier) -> None:
    """A second readout module would have been a copy of all of this, kept in step by hand."""
    kwargs = dict(opts=opts, fov_mm=250.0, matrix=128, bandwidth_hz_px=200.0,
                  partial_fourier=partial_fourier)
    winding = sc.modules.CartesianLine(**kwargs)
    spin_echo = sc.modules.CartesianLine(**kwargs, prephase=False)

    for name in DESIGNED:
        assert getattr(spin_echo, name) == getattr(winding, name), name
    assert float(spin_echo.gx.amplitude) == float(winding.gx.amplitude)
    assert float(spin_echo.gx.flat_time) == float(winding.gx.flat_time)
    assert sc.events.content_hash(spin_echo.gx) == sc.events.content_hash(winding.gx)
    assert sc.events.content_hash(spin_echo.adc) == sc.events.content_hash(winding.adc)


@pytest.mark.parametrize('partial_fourier', [1.0, 0.75, 0.625])
def test_area_to_echo_is_minus_the_prephaser(opts, partial_fourier) -> None:
    """
    The physics number and the event's number, stated as the identity between them.

    ``area_to_echo_per_m`` is what a spin echo's crusher pair is balanced around, and it is
    available with **no** prephaser -- which is what lets one module serve both sequences.
    """
    kwargs = dict(opts=opts, fov_mm=250.0, matrix=128, bandwidth_hz_px=200.0,
                  partial_fourier=partial_fourier)
    winding = sc.modules.CartesianLine(**kwargs)
    spin_echo = sc.modules.CartesianLine(**kwargs, prephase=False)

    assert winding.area_to_echo_per_m == pytest.approx(-winding.prephaser_area_per_m, abs=1e-9)
    assert spin_echo.area_to_echo_per_m == pytest.approx(winding.area_to_echo_per_m, abs=1e-12)
    assert spin_echo.area_to_echo_per_m == pytest.approx(
        area_until(spin_echo.gx, spin_echo.time_to_echo()), abs=1e-12)
@pytest.mark.parametrize('partial_fourier', [1.0, 0.75, 0.625])
def test_the_two_halves_of_the_lobe_sum_to_it_and_their_difference_is_the_imbalance(
    opts, partial_fourier,
) -> None:
    """
    The geometry a spin-echo caller needs, stated by the readout rather than re-derived by each.

    A refocusing pulse conjugates k, so a caller straddling this readout with a lobe either side
    of the pulse needs the area before the echo to equal the area after it -- and the two lobes
    therefore differ by exactly this imbalance.  What the readout cannot decide is how to split
    that difference between them, because that depends on the crusher window and the pulse it
    surrounds; ``TSEShot`` does.

    Under ``partial_fourier`` the echo is nowhere near the middle of the lobe, which is the case
    that separates this from "half the area".
    """
    line = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=128, bandwidth_hz_px=200.0,
                                    partial_fourier=partial_fourier, prephase=False)

    assert line.area_to_echo_per_m + line.area_after_echo_per_m == pytest.approx(
        float(line.gx.area), abs=1e-12)
    assert line.echo_moment_imbalance_per_m == pytest.approx(
        line.area_after_echo_per_m - line.area_to_echo_per_m, abs=1e-12)
    # Never zero: the flat top is rounded up onto the gradient raster with the slack after the
    # last sample, and the echo sample sits half a dwell off the window's centre.
    assert line.echo_moment_imbalance_per_m != 0.0
    if partial_fourier < 1.0:
        # A partial-Fourier readout starts closer to k = 0, so far more area follows the echo.
        assert line.area_after_echo_per_m > line.area_to_echo_per_m


def test_prephase_false_is_the_gradient_and_its_adc_and_nothing_else(opts) -> None:
    """The block starts at zero, so the caller's own lobe is what precedes it."""
    ro = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=128, bandwidth_hz_px=200.0,
                                  prephase=False)

    block = ro()

    assert sorted(n.item.type for n in block) == ['adc', 'trap']
    assert block.duration == pytest.approx(float(pp.calc_duration(ro.gx)))
    assert 0.0 < ro.time_to_echo() < block.duration
    assert ro(acquire=False).duration == block.duration
    assert [n.item.type for n in ro(acquire=False)] == ['trap']


def test_time_to_echo_measures_from_the_readout_gradient_when_there_is_no_prephaser(opts) -> None:
    """Both are 'from the start of the block build returns' -- and the block starts elsewhere."""
    winding = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=128, bandwidth_hz_px=200.0)
    spin_echo = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=128,
                                         bandwidth_hz_px=200.0, prephase=False)

    assert spin_echo.time_to_echo() == pytest.approx(
        winding.time_to_echo() - winding.prephaser_duration_s, abs=1e-12)


@pytest.mark.parametrize('name', ['prephaser_duration_s', 'prephaser_area_per_m'])
def test_asking_for_a_prephaser_that_does_not_exist_raises(opts, name) -> None:
    """The shape the library already uses for an argument that cannot take effect."""
    ro = sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=128, bandwidth_hz_px=200.0,
                                  prephase=False)

    with pytest.raises(sc.ConfigurationError, match='prephase=False'):
        getattr(ro, name)


def test_a_prephaser_duration_with_prephase_false_names_both_arguments(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='prephase') as caught:
        sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=128, bandwidth_hz_px=200.0,
                                 prephase=False, prephaser_duration_s=600e-6)

    assert 'prephaser_duration_s' in str(caught.value)


def test_the_spin_echo_readout_is_pure_and_compiles_alone(opts, component_checks) -> None:
    """Starting at zero is what makes that true; a prephaser is not what made it legal."""
    component_checks.all(
        sc.modules.CartesianLine(opts=opts, fov_mm=250.0, matrix=64, bandwidth_hz_px=250.0,
                                 prephase=False),
    )


# =============================================================================================
# The multi-echo train
# =============================================================================================
#
# Four numbers here are wrong *quietly*, and every test below isolates one of them:
#
# 1. **The reverse lobe's grid.**  Build a bipolar train out of the flat-top-rounded lobe and the
#    two polarities sample grids 0.2051 dk apart.  Every k-space extent check passes, the image is
#    an image, and the error surfaces only as a wrong value in the map the sequence exists to
#    produce.  ``test_sc_kspace_sees_a_reversed_lobe_correctly`` and
#    ``test_the_reverse_lobe_lands_on_the_forward_grid``.
# 2. **The fly-back area.**  Minus the lobe's *total* area.  The plausible wrong answer,
#    ``-area_to_echo_per_m``, is a number the module already exposes and it is only a factor of
#    two out.  ``test_the_flyback_cancels_the_whole_lobe``.
# 3. **The dwell.**  The guard is an ADC ``delay``, and pypulseq divides that by
#    ``rf_raster_time``, not ``adc_raster_time``.  ``test_an_adc_delay_off_the_rf_raster_is_refused``
#    is the oracle and ``test_the_bipolar_dwell_lands_on_the_quantum`` the consequence.
# 4. **The echo times.**  A bipolar train's dTE alternates by two dwells, and its partial-Fourier
#    version by 1953 us.  Nothing checks it; a reconstruction that assumes a uniform dTE simply
#    returns a scaled field map.  ``test_the_echo_times_alternate_by_the_closed_form``.

REFERENCE = dict(fov_mm=220.0, matrix=128, axis='x')


@pytest.fixture(scope='module')
def fast_opts() -> Opts:
    """The reference scanner of ``docs/multi_echo_api.md``: 40 mT/m, 180 T/m/s, 10 us dead time."""
    return Opts(
        max_grad=40, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s', B0=3.0,
        rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
    )


def _line(opts, **kwargs) -> sc.modules.CartesianLine:
    return sc.modules.CartesianLine(opts=opts, **REFERENCE, **kwargs)


def _k_rows(line: sc.modules.CartesianLine, opts: Opts) -> np.ndarray:
    """``k`` along the readout axis, reshaped ``(echo, sample)``, from pypulseq's integration."""
    axis = 'xyz'.index(line.axis)
    k = sc.kspace(sc.LogicBlock('probe').add(0.0, line()), opts)['k_adc'][axis]
    return k.reshape(line.echoes, line.num_samples)


def _gradient_kinds(seq) -> set[str]:
    return {getattr(g, 'type', '?') for i in range(1, len(seq.block_events) + 1)
            for g in (getattr(seq.get_block(i), 'gx', None),) if g is not None}


def _events(tree, kind: str) -> list:
    return [e for _, e, _ in sc.flatten(tree) if getattr(e, 'type', '') == kind]


# ------------------------------------------------------------- Phase 0: the two oracles
#
# Both are cheap, both invalidate everything after them if they fail, and both are written as
# tests because the oracle for every physics assertion below is ``sc.kspace`` and pypulseq's own
# timing check -- and an oracle nothing checks is an assumption.

def _raw_two_echo(opts: Opts, *, centred: bool) -> tuple[sc.LogicBlock, int, float]:
    """A prephaser, a forward lobe, a barrier and a reverse lobe, out of raw pypulseq."""
    from seqcraft.design.timing import from_ticks, to_ticks
    from seqcraft.modules._support import ceil_raster, halve_onto

    n, dk = 128, 1e3 / 220.0
    dwell = 31.5e-6 if centred else 31.2e-6
    amplitude = dk / dwell
    sample_s = from_ticks(n * to_ticks(dwell))
    if centred:
        rise_s = ceil_raster(amplitude / float(opts.max_slew), opts.grad_raster_time)
        total_s = ceil_raster(sample_s + 2 * max(float(opts.adc_dead_time), rise_s),
                              opts.grad_raster_time)
        guard_s = halve_onto(total_s, sample_s, opts.rf_raster_time, opts=opts)
        gx = pp.make_trapezoid(channel='x', amplitude=amplitude, rise_time=rise_s,
                               flat_time=total_s - 2 * rise_s, system=opts)
    else:
        # Today's lobe: the flat top rounded UP onto the gradient raster, the slack left after the
        # last sample.  Correct for one line -- the slack changes no k position -- and wrong the
        # moment a reverse lobe has to sample the same grid.
        gx = pp.make_trapezoid(channel='x', amplitude=amplitude, system=opts,
                               flat_time=ceil_raster(sample_s, opts.grad_raster_time))
        guard_s = float(gx.rise_time)
    adc = pp.make_adc(num_samples=n, dwell=dwell, delay=guard_s, system=opts)
    lobe_s = float(pp.calc_duration(gx))
    pre = pp.make_trapezoid(channel='x', system=opts,
                            area=-area_until(gx, guard_s + (n // 2 + 0.5) * dwell))
    out = sc.LogicBlock('raw').add(0.0, pre)
    t0 = float(pp.calc_duration(pre))
    out.add(t0, gx).add(t0, adc)
    out.add(t0 + lobe_s, sc.barrier('seam'))
    out.add(t0 + lobe_s, sc.events.derive(gx, amplitude=-amplitude)).add(t0 + lobe_s, adc)
    return out, n, dk


def _signed_mismatch(block, opts: Opts, n: int) -> float:
    """Forward samples against the reverse lobe's reversed, in 1/m.  Zero is one grid."""
    kx = np.asarray(sc.kspace(block, opts)['k_adc'])[0]
    return float(np.max(np.abs(kx[:n] - kx[n:2 * n][::-1])))


@pytest.mark.parametrize(('centred', 'worst_dk'), [(True, 1e-9), (False, 0.2051)])
def test_sc_kspace_sees_a_reversed_lobe_correctly(fast_opts, centred, worst_dk) -> None:
    """
    **Oracle 1**, and *both halves are the test*: the second is what says the first is not vacuous.

    A centred window puts the two polarities on one grid to 1.7e-13 1/m.  The flat-top-rounded
    lobe -- which is what :class:`CartesianLine` correctly builds for a *line*, and what a
    monopolar train reuses unchanged -- puts them **0.2051 dk** apart, and the whole design of
    ``polarity='bipolar'`` exists to stop that.
    """
    block, n, dk = _raw_two_echo(fast_opts, centred=centred)

    mismatch = _signed_mismatch(block, fast_opts, n) / dk
    if centred:
        assert mismatch < worst_dk, f'{mismatch:.4f} dk apart'
    else:
        assert mismatch == pytest.approx(worst_dk, abs=1e-4), f'{mismatch:.4f} dk apart'


def test_a_k_space_extent_check_would_pass_on_both(fast_opts) -> None:
    """
    **This is the check a k-space extent test would pass**, written down so that nobody adds one
    and believes it.  ``max|k|`` is the same to a part in a thousand either way; only the *signed*
    difference between the two polarities separates them, and it separates them by 1e12.
    """
    good, n, _ = _raw_two_echo(fast_opts, centred=True)
    bad, _, _ = _raw_two_echo(fast_opts, centred=False)
    extent = lambda b: float(np.abs(np.asarray(  # noqa: E731
        sc.kspace(b, fast_opts)['k_adc'])[0]).max())

    assert extent(good) == pytest.approx(extent(bad), rel=2e-3), 'the extent check cannot tell'
    assert _signed_mismatch(bad, fast_opts, n) > 1e10 * _signed_mismatch(good, fast_opts, n), (
        'and the signed difference can, which is why every k assertion here is signed'
    )


def test_an_adc_delay_off_the_rf_raster_is_refused(fast_opts) -> None:
    """
    **Oracle 2**, and the fact :func:`~seqcraft.modules._support.dwell_quantum` exists for.

    ``pypulseq.Sequence.check_timing._check_timing_block`` selects ``rf_raster_time`` for ``adc``,
    ``rf`` and ``output`` events and ``grad_raster_time`` for everything else.  That is **not in
    pulseq's specification**: it is two lines of pypulseq source, and it is the whole reason a
    dwell cannot simply be snapped onto the ADC raster when the guard has to be symmetric.

    23.2 us is a whole number of 0.1 us ADC rasters and *not* of 1 us RF rasters.  23.0 us is both.
    """
    def probe(delay_s: float) -> sc.LogicBlock:
        gx = pp.make_trapezoid(channel='x', amplitude=1e5, flat_time=1e-3, rise_time=2e-4,
                               system=fast_opts)
        adc = pp.make_adc(num_samples=64, dwell=1e-5, delay=delay_s, system=fast_opts)
        return sc.LogicBlock('probe').add(0.0, gx).add(0.0, adc)

    with pytest.raises(sc.SeqCraftError) as caught:
        sc.compile(probe(23.2e-6), fast_opts)
    assert 'adc' in str(caught.value).lower(), str(caught.value)

    sc.compile(probe(23.0e-6), fast_opts)                # the same block, one RF raster over


# ------------------------------------------------------------------ the echoes=1 baseline
#
# The test that makes the whole revision safe.  It goes first because everything after it is
# allowed to move and this is not.

#: ``(bandwidth_hz_px, partial_fourier) -> [(time, content_hash), ...]``, captured from the
#: single-echo module **as it stood before the train existed** and pinned here.  Regenerate only
#: for a deliberate change to the *line*, never for one to the train: a train that moved these has
#: moved every GRE, MPRAGE and SE example in the repository.
_BASELINE: dict[tuple[float, float], list[tuple[float, str]]] = {
    (250.0, 1.0): [
        (0.0, 'a195270dc2fda1ce028a793de261e379978572162915eefdbf9dd8cb2a451ea0'),
        (0.0004, '7db8250c09220a682064152eea0a0174c19070da7143bdd282b2e5ab55d34388'),
        (0.0004, '7e1c0d6d87d3d41459249c73b6022b15984e7f346b9d71aa4229c711f9909d0b'),
    ],
    (250.0, 0.75): [
        (0.0, '8b410d955fc24f7568c9f1ce03e7d5ba964a94df0db74e3e3c80d2b2e632cf71'),
        (0.00028, '27fc857cbe83ba98a7b0f3aa81e1c00391cd30f805e3d23d62f52082599ebd91'),
        (0.00028, '3e2302f145046bba05ea5d9de4e95db5f12fa4342648662be91be31cd45690c0'),
    ],
    (250.0, 0.625): [
        (0.0, 'd648e7af31c5e5fc7188d0225b8ae4f63edf6815e47f720b96a6f48cc8870a64'),
        (0.0002, '1a59527e0f868d00da6ecdc5fae2fa8e9656eb44ee5dfb3237c131987f63f847'),
        (0.0002, '3b9356d8096ff3a98735e5a1bd66c44945bfd721070068cfd496e6529380c41d'),
    ],
    (500.0, 1.0): [
        (0.0, '8c16425e7fe3f608f2a6e90902f777b6ecfb7e57da326cfc865e57d631344b87'),
        (0.0004, '523e4903ab6c53f98b15f497aa5b790b514393c47c27052a88b1c0c0d47a4c72'),
        (0.0004, '8ce3490e529f363fdc8bf058560d634f2b7a7bbb68256d72264553e8ccbee1fc'),
    ],
    (500.0, 0.75): [
        (0.0, 'f8db23990fedef601c17109dab5a4ef3b4c37bc0ba47748982c56fa762d513fe'),
        (0.0003, '4339ba3227d3de86cc6ed1704977b3428485c977b2345c1b86e5c742fbf333e9'),
        (0.0003, '6f698797b09af632b6b59b81394352d426e691286bd524d1b0f15a7f72f8af50'),
    ],
    (500.0, 0.625): [
        (0.0, '2b75ed259cb2fdda8ac2c32136e81e0cdd86a62c7d4425c61d8196cadaf34677'),
        (0.00022, '7c21a98e70dd70eb8205a20f3aa587def18bb2c040e4221190b768740ffca79e'),
        (0.00022, '92f8a3114d965c82405cccbfc501e4dfd09d1acfb3efe77a71423d765ae35835'),
    ],
    (700.0, 1.0): [
        (0.0, 'c66500e7128c628283a8caa94bc8cbc8a39db9be670e0e530a80d8a813353e9f'),
        (0.0004, 'aa490d5ef8624cfc82689bdc382af1bf0c20b1d43677e392e8bf09423c607d72'),
        (0.0004, 'f62c7dd313729729f90b3c2a373438220cb949854837115e9e7061ed72980f75'),
    ],
    (700.0, 0.75): [
        (0.0, 'c24f9c5af196ce54419914af4dd8aed3cf86927640eea45c34a9fa0ef080340b'),
        (0.0003, '9e40344e640684ce43bd2bc4a21fa09edcebed2e208fe6ce82da8c458d7a50b8'),
        (0.0003, 'fa6a475834e0a39d9cb8164d8a72aa3c2014946d0154c1ee934e0e3ffcdd3727'),
    ],
    (700.0, 0.625): [
        (0.0, '524884d59b37be94726cab906fe8eaecc050925f3401c1682d1212269ec3aa91'),
        (0.00022, 'd7019069baca3ff5afb774bf353e7da69e972a57db86d1020fcea50273fa4d8f'),
        (0.00022, 'dda9a5fb8cfbd339856f9c3ae990c09124e6613523d734361a643be3113201a6'),
    ],
}


@pytest.mark.parametrize('bandwidth', [250.0, 500.0, 700.0])
@pytest.mark.parametrize('partial_fourier', [1.0, 0.75, 0.625])
def test_one_echo_is_the_module_that_shipped(fast_opts, bandwidth, partial_fourier) -> None:
    """
    ``echoes=1`` is byte-identical: same dwell, same lobe, same prephaser, same block, same
    duration, same ``time_to_echo()``.

    Pinned against digests captured from the module **before** the train existed, so a failure
    here is a claim about every existing example rather than about this feature.
    """
    line = _line(fast_opts, bandwidth_hz_px=bandwidth, partial_fourier=partial_fourier)

    digest = sorted((round(t, 12), sc.events.content_hash(e))
                    for t, e, _ in sc.flatten(line()))
    assert digest == _BASELINE[bandwidth, partial_fourier]
    assert (line.echoes, line.polarity) == (1, None)


# ----------------------------------------------------------------------- the refusals
@pytest.mark.parametrize('polarity', ['monopolar', 'bipolar'])
def test_polarity_at_one_echo_raises_naming_both(fast_opts, polarity) -> None:
    """An argument that cannot take effect is reported rather than ignored."""
    with pytest.raises(sc.ConfigurationError, match='echoes=1') as caught:
        _line(fast_opts, bandwidth_hz_px=250.0, polarity=polarity)
    assert 'polarity' in str(caught.value)


def test_a_train_without_a_polarity_raises_and_says_what_each_costs(fast_opts) -> None:
    """
    No default, on purpose.  The two modes produce files that differ in the dwell, the echo times,
    the k ordering of every second echo and the block count -- and look identical in a protocol
    printout, which is why the message has to say what each one buys.
    """
    with pytest.raises(sc.ConfigurationError, match='polarity is required') as caught:
        _line(fast_opts, bandwidth_hz_px=250.0, echoes=8)

    message = str(caught.value)
    assert 'monopolar' in message and 'bipolar' in message
    assert 'fly-back' in message


def test_an_unknown_polarity_lists_the_two(fast_opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='monopolar') as caught:
        _line(fast_opts, bandwidth_hz_px=250.0, echoes=4, polarity='alternating')
    assert 'bipolar' in str(caught.value)


def test_echo_spacing_at_one_echo_raises(fast_opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='echoes=1') as caught:
        _line(fast_opts, bandwidth_hz_px=250.0, echo_spacing_s=5e-3)
    assert 'echo_spacing_s' in str(caught.value)


@pytest.mark.parametrize('name', ['echo_spacing_s', 'min_echo_spacing_s', 'flyback_area_per_m',
                                  'flyback_duration_s'])
def test_a_train_question_at_one_echo_raises(fast_opts, name) -> None:
    """The shape ``prephaser_area_per_m`` already uses with ``prephase=False``."""
    line = _line(fast_opts, bandwidth_hz_px=250.0)
    with pytest.raises(sc.ConfigurationError):
        getattr(line, name)


@pytest.mark.parametrize('echo', [-1, 8, 99])
def test_an_echo_this_train_does_not_have_raises(fast_opts, echo) -> None:
    train = _line(fast_opts, bandwidth_hz_px=250.0, echoes=8, polarity='monopolar')
    with pytest.raises(sc.ConfigurationError, match='outside this train'):
        train.time_to_echo(echo)


@pytest.mark.parametrize('value', [0, -1])
def test_a_non_positive_echo_count_raises(fast_opts, value) -> None:
    with pytest.raises(sc.ConfigurationError, match='echoes'):
        _line(fast_opts, bandwidth_hz_px=250.0, echoes=value)


def test_bipolar_with_prephase_false_raises(fast_opts) -> None:
    """
    Legal in principle, unmeasured in fact: a multi-echo *spin* echo with alternating polarity has
    to balance the reversed lobes across the refocusing pulse as well, which is ``se_api.md``'s
    "which pair is a pair" question with a sign flip in it.  It raises until somebody measures it.
    """
    with pytest.raises(sc.ConfigurationError, match='prephase=False') as caught:
        _line(fast_opts, bandwidth_hz_px=250.0, echoes=4, polarity='bipolar', prephase=False)
    assert 'monopolar' in str(caught.value)


def test_a_monopolar_spin_echo_train_builds(fast_opts) -> None:
    """``prephase=False`` works with a train, and ``area_to_echo_per_m`` still means echo one."""
    train = _line(fast_opts, bandwidth_hz_px=250.0, echoes=4, polarity='monopolar',
                  prephase=False)
    one = _line(fast_opts, bandwidth_hz_px=250.0, prephase=False)

    assert train.area_to_echo_per_m == pytest.approx(one.area_to_echo_per_m, abs=1e-12)
    assert train.te_s[0] == pytest.approx(one.time_to_echo(), abs=1e-15)
    assert train().duration == pytest.approx(one().duration + 3 * train.echo_spacing_s, abs=1e-12)


# -------------------------------------------------------------------------- monopolar
def test_the_monopolar_lobe_is_the_shipped_lobe(fast_opts) -> None:
    """
    Byte-identical, asserted rather than intended.  The whole argument for one module rather than
    two is that the shared half does not move, and this is what makes that a claim.
    """
    one = _line(fast_opts, bandwidth_hz_px=250.0)
    train = _line(fast_opts, bandwidth_hz_px=250.0, echoes=8, polarity='monopolar')

    assert sc.events.content_hash(train.gx) == sc.events.content_hash(one.gx)
    assert sc.events.content_hash(train.adc) == sc.events.content_hash(one.adc)
    assert (train.dwell_s, train.bandwidth_hz_px, train.num_samples, train.pre_echo_samples) == (
        one.dwell_s, one.bandwidth_hz_px, one.num_samples, one.pre_echo_samples)
    assert train.area_to_echo_per_m == one.area_to_echo_per_m


def test_the_flyback_cancels_the_whole_lobe(fast_opts) -> None:
    """
    **The dangerous wrong answer is excluded by name.**  ``-area_to_echo_per_m`` is the pre-echo
    part only -- a number this module already computes and already exposes -- and it is
    -294.638695 against -585.664336 1/m.  Each echo would start half a k-space further along than
    the last, and the reconstruction would sum eight partial acquisitions believing they were
    eight echoes of one line.
    """
    train = _line(fast_opts, bandwidth_hz_px=250.0, echoes=6, polarity='monopolar')
    lobe = area_until(train.gx, float(pp.calc_duration(train.gx)))

    assert train.flyback_area_per_m == pytest.approx(-lobe, abs=1e-9)
    assert round(train.flyback_area_per_m, 6) == -585.664336
    assert train.flyback_area_per_m != pytest.approx(-train.area_to_echo_per_m, rel=0.1)
    assert train.flyback_area_per_m != pytest.approx(-float(train.gx.flat_area), rel=1e-3)
    assert train.flyback_area_per_m != pytest.approx(-float(train.gx.area) / 2, rel=0.1)


def test_the_monopolar_echo_times_are_uniform_to_the_nanosecond(fast_opts) -> None:
    """And every gap is exactly ``echo_spacing_s``, which is what makes the period the dTE here."""
    train = _line(fast_opts, bandwidth_hz_px=500.0, echoes=8, polarity='monopolar')
    gaps = np.diff(train.te_s)

    assert float(np.ptp(gaps)) * 1e9 == pytest.approx(0.0, abs=1e-3)
    assert np.allclose(gaps, train.echo_spacing_s, atol=1e-15)
    assert len(train.te_s) == train.echoes == 8
    assert train.time_to_echo() == train.te_s[0]
    assert train.time_to_echo(7) == train.te_s[7]


def test_monopolar_polarity_and_echo_sample_never_reverse(fast_opts) -> None:
    train = _line(fast_opts, bandwidth_hz_px=500.0, echoes=6, polarity='monopolar',
                  partial_fourier=0.75)

    assert [train.polarity_of(n) for n in range(6)] == [1] * 6
    assert [train.echo_sample(n) for n in range(6)] == [train.pre_echo_samples] * 6


# --------------------------------------------------------------------------- bipolar
def test_the_bipolar_dwell_lands_on_the_quantum(fast_opts) -> None:
    """
    The same request produces a different bandwidth in the two modes, by up to 2.5 %.  That is a
    fact about the file, it belongs in the protocol table, and ``bandwidth_hz_px`` reads it back
    -- which is why that attribute has always been the achieved value rather than the requested.
    """
    from seqcraft.modules._support import dwell_quantum

    mono = _line(fast_opts, bandwidth_hz_px=500.0, echoes=8, polarity='monopolar')
    bi = _line(fast_opts, bandwidth_hz_px=500.0, echoes=8, polarity='bipolar')
    quantum = dwell_quantum(128, opts=fast_opts)

    assert round(quantum * 1e9) == 500
    assert round(bi.dwell_s * 1e12) % round(quantum * 1e12) == 0
    assert bi.dwell_s != mono.dwell_s
    assert bi.dwell_s > 1.0 / (500.0 * 128), 'rounded up, never down'
    assert (round(mono.bandwidth_hz_px, 1), round(bi.bandwidth_hz_px, 1)) == (500.8, 488.3)


@pytest.mark.parametrize('bandwidth', [200.0, 350.0, 500.0, 700.0, 1000.0])
def test_rule_one_holds_with_equality(fast_opts, bandwidth) -> None:
    """``T - 2*guard - N*dwell == 0``, which is the invariant that names Rule 1.  Measured 0 ns."""
    bi = _line(fast_opts, bandwidth_hz_px=bandwidth, echoes=4, polarity='bipolar')
    guard = float(bi.adc.delay)

    residue = bi.echo_spacing_s - 2 * guard - bi.num_samples * bi.dwell_s
    assert residue * 1e9 == pytest.approx(0.0, abs=1e-6), f'{residue * 1e9:.3f} ns'
    assert guard >= float(fast_opts.adc_dead_time)
    assert guard >= float(bi.gx.rise_time), 'the window must be entirely on the flat top'


def test_the_reverse_lobe_lands_on_the_forward_grid(fast_opts) -> None:
    """
    Rule 1 as a statement about the **compiled waveform**, not about the arithmetic that produced
    it.  Measured 1.7e-13 1/m against a dk of 4.5455; the negative half of this claim is
    ``test_sc_kspace_sees_a_reversed_lobe_correctly``, which measures 0.2051 dk on the lobe this
    module correctly builds for a single line.
    """
    bi = _line(fast_opts, bandwidth_hz_px=250.0, echoes=6, polarity='bipolar')
    rows = _k_rows(bi, fast_opts)

    grid = rows[0]
    for n in range(bi.echoes):
        got = rows[n] if bi.polarity_of(n) > 0 else rows[n][::-1]
        assert np.abs(got - grid).max() < 1e-9, f'echo {n} is on a different grid'


@pytest.mark.parametrize('polarity', ['monopolar', 'bipolar'])
@pytest.mark.parametrize('echoes', [1, 2, 6, 8])
@pytest.mark.parametrize('partial_fourier', [1.0, 0.75])
def test_k_is_zero_at_the_echo_sample_of_every_echo(
    fast_opts, polarity, echoes, partial_fourier,
) -> None:
    """
    Signed, at the sample :meth:`echo_sample` names -- which is the assertion that says the mirror
    is right.  The wrong branch is not an error; it is a mirrored echo, and against a symmetric
    phantom that looks entirely correct.
    """
    extra = {} if echoes == 1 else {'polarity': polarity}
    line = _line(fast_opts, bandwidth_hz_px=250.0, echoes=echoes,
                 partial_fourier=partial_fourier, **extra)
    rows = _k_rows(line, fast_opts)

    at_echo = np.array([rows[n][line.echo_sample(n)] for n in range(echoes)])
    assert np.abs(at_echo).max() < 1e-9, f'{np.abs(at_echo).max():.3e} 1/m'


@pytest.mark.parametrize(('partial_fourier', 'echoes'), [(1.0, 8), (0.75, 4)])
def test_the_echo_times_alternate_by_the_closed_form(fast_opts, partial_fourier, echoes) -> None:
    """
    ``2 * |N - 1 - 2*pre_echo| * dwell``, asserted against the closed form rather than against the
    measured constant -- with the constant printed in the failure message, because that is the
    number a protocol table quotes.  Measured **32.00 us** at pf 1.0 and **1953.00 us** at 0.75,
    which is comparable to the spacing itself.

    The data is exact either way.  What is not merely imprecise but *meaningless* is the uniform
    dTE a caller might assume, which is why ``te_s`` is the only sanctioned source of echo times.
    """
    bi = _line(fast_opts, bandwidth_hz_px=500.0, echoes=echoes, polarity='bipolar',
               partial_fourier=partial_fourier)
    predicted = 2 * abs(bi.num_samples - 1 - 2 * bi.pre_echo_samples) * bi.dwell_s
    spread = float(np.ptp(np.diff(bi.te_s)))

    assert spread == pytest.approx(predicted, abs=1e-12), (
        f'{spread * 1e6:.2f} us measured against {predicted * 1e6:.2f} us predicted'
    )
    gaps = sorted({round(g, 12) for g in np.diff(bi.te_s)})
    assert len(gaps) == 2, f'the gaps should take exactly two values, got {gaps}'
    assert float(np.mean(gaps)) == pytest.approx(bi.echo_spacing_s, abs=1e-12), (
        'the period is the mean of the two, and it is not either of them'
    )


def test_echo_sample_alternates_under_partial_echo(fast_opts) -> None:
    """32, 63, 32, 63 -- measured, and the reason ``echo_sample`` is a method and not a constant."""
    bi = _line(fast_opts, bandwidth_hz_px=500.0, echoes=4, polarity='bipolar',
               partial_fourier=0.75)

    assert (bi.num_samples, bi.pre_echo_samples) == (96, 32)
    assert [bi.echo_sample(n) for n in range(4)] == [32, 63, 32, 63]
    assert [bi.polarity_of(n) for n in range(4)] == [1, -1, 1, -1]


def test_bipolar_forces_an_even_sample_count(fast_opts) -> None:
    """An odd count leaves the guard half a raster short, which is a grid the reversal misses."""
    odd = sc.modules.CartesianLine(opts=fast_opts, fov_mm=220.0, matrix=127, axis='x',
                                   bandwidth_hz_px=500.0, echoes=4, polarity='bipolar')
    assert odd.num_samples % 2 == 0


def test_the_bipolar_period_is_shorter_than_the_monopolar_one(fast_opts) -> None:
    """11.5 % at 250 Hz/px, which is what the fly-back costs and the whole trade in one number."""
    mono = _line(fast_opts, bandwidth_hz_px=250.0, echoes=6, polarity='monopolar')
    bi = _line(fast_opts, bandwidth_hz_px=250.0, echoes=6, polarity='bipolar')

    assert bi.echo_spacing_s < mono.echo_spacing_s
    assert round(mono.echo_spacing_s * 1e6, 1) == 4610.0
    assert round(bi.echo_spacing_s * 1e6, 1) == 4080.0
    assert round(mono.flyback_duration_s * 1e6, 1) == 570.0


@pytest.mark.parametrize('name', ['flyback_area_per_m', 'flyback_duration_s'])
def test_the_flyback_is_refused_under_bipolar(fast_opts, name) -> None:
    """Naming ``echo_spacing_s`` as what a caller wanting to lengthen the period should use."""
    bi = _line(fast_opts, bandwidth_hz_px=250.0, echoes=4, polarity='bipolar')
    with pytest.raises(sc.ConfigurationError, match='fly-back') as caught:
        getattr(bi, name)
    assert 'echo_spacing_s' in str(caught.value)


# -------------------------------------------------------------------------- the barriers
@pytest.mark.parametrize(('polarity', 'blocks'), [('monopolar', 16), ('bipolar', 9)])
def test_the_train_compiles_to_trapezoids_and_no_merge_warning(
    fast_opts, polarity, blocks,
) -> None:
    """
    **The barriers' test.**  Left to choose, the compiler cuts the midpoint of the legal gap
    between consecutive ADCs, which lands inside a gradient, and splits it into two arbitrary
    waveforms.  Nothing else fails when that happens: it compiles, it passes every k-space check,
    and the only report of it is a ``merge`` warning naming the readout axis.

    The block counts are asserted rather than the warning's absence alone, because the counts are
    what say the *second* monopolar barrier is there.
    """
    train = _line(fast_opts, bandwidth_hz_px=250.0, echoes=8, polarity=polarity)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', sc.SeqCraftWarning)
        seq = sc.compile(train(), fast_opts)

    assert len(seq.block_events) == blocks
    assert _gradient_kinds(seq) == {'trap'}
    assert not [w for w in caught if 'merge' in str(w.message)], [str(w.message) for w in caught]


@pytest.mark.parametrize('polarity', ['monopolar', 'bipolar'])
@pytest.mark.parametrize('bandwidth', [250.0, 500.0, 1000.0])
def test_the_barriers_change_nothing_against_this_compiler(fast_opts, polarity, bandwidth) -> None:
    """
    **The negative half, and it says the opposite of what the design document predicted.**

    Two barriers, one, and none produce the *identical* compiled sequence: the same block count,
    every gradient a ``trap``, and no ``merge`` warning, at three bandwidths and both polarities.
    :func:`~seqcraft.compiler.boundaries.find_boundaries` already offers every gradient edge as a
    boundary candidate and accepts one that falls strictly inside no gradient -- and every seam
    here is exactly such an edge, because the lobe ends where the fly-back starts and the fly-back
    ends where the next lobe starts, with nothing overlapping either join.

    That is the difference from :class:`~seqcraft.modules.EPI2D`, whose blip is *centred on* the
    seam and therefore covers it.  The barriers stay because the opportunistic rule is documented
    as a preference rather than a guarantee; this test is what stops anyone claiming they are what
    keeps the lobes whole.
    """
    train = _line(fast_opts, bandwidth_hz_px=bandwidth, echoes=8, polarity=polarity)
    period, lobe = train.echo_spacing_s, float(pp.calc_duration(train.gx))

    def compiled(barriers: int):
        out = sc.LogicBlock(f'{barriers}-barriers').add(0.0, train._prephaser)
        start = train.prephaser_duration_s
        for n in range(8):
            t0 = start + n * period
            out.add(t0, train.gx if train.polarity_of(n) > 0 else train._gx_reverse)
            out.add(t0, train.adc)
            if n + 1 < 8:
                if barriers >= 1:
                    out.add(t0 + lobe, sc.barrier('seam'))
                if train._flyback is not None:
                    out.add(t0 + lobe, train._flyback)
                    if barriers >= 2:
                        out.add(start + (n + 1) * period, sc.barrier('lobe'))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always', sc.SeqCraftWarning)
            seq = sc.compile(out, fast_opts)
        merges = [w for w in caught if 'merge' in str(w.message)]
        return len(seq.block_events), _gradient_kinds(seq), len(merges)

    two = compiled(2)
    assert two == compiled(1) == compiled(0), f'the barriers are load-bearing after all: {two}'
    assert two[1] == {'trap'} and two[2] == 0


# --------------------------------------------------------------------- a requested period
def test_a_requested_spacing_lengthens_the_flyback_and_not_the_lobe(fast_opts) -> None:
    """
    A longer fly-back at the same area is a lower amplitude and a lower slew rate, which is the
    eddy-current-cheaper of the two ways to spend the time.  A gap inserted after it would be
    neither, and would leave the lobe's own rounding where it is.
    """
    short = _line(fast_opts, bandwidth_hz_px=250.0, echoes=4, polarity='monopolar')
    stretched = _line(fast_opts, bandwidth_hz_px=250.0, echoes=4, polarity='monopolar',
                      echo_spacing_s=6.0e-3)

    assert sc.events.content_hash(stretched.gx) == sc.events.content_hash(short.gx)
    assert stretched.echo_spacing_s == pytest.approx(6.0e-3, abs=1e-12)
    assert stretched.flyback_duration_s > short.flyback_duration_s
    assert stretched.flyback_area_per_m == pytest.approx(short.flyback_area_per_m, abs=1e-9)
    assert abs(float(stretched._flyback.amplitude)) < abs(float(short._flyback.amplitude))
    assert stretched().duration == pytest.approx(
        short().duration + 3 * (6.0e-3 - short.echo_spacing_s), abs=1e-12)


def test_a_requested_spacing_keeps_rule_one_under_bipolar(fast_opts) -> None:
    """Lengthening ``T`` grows both guards symmetrically, which keeps Rule 1 rather than breaking it."""
    bi = _line(fast_opts, bandwidth_hz_px=250.0, echoes=4, polarity='bipolar',
               echo_spacing_s=5.0e-3)
    guard = float(bi.adc.delay)

    assert bi.echo_spacing_s == pytest.approx(5.0e-3, abs=1e-12)
    assert (bi.echo_spacing_s - 2 * guard - bi.num_samples * bi.dwell_s) * 1e9 == pytest.approx(
        0.0, abs=1e-6)
    rows = _k_rows(bi, fast_opts)
    assert np.abs(rows[0] - rows[1][::-1]).max() < 1e-9


@pytest.mark.parametrize('polarity', ['monopolar', 'bipolar'])
def test_a_spacing_below_the_minimum_raises_naming_it(fast_opts, polarity) -> None:
    train = _line(fast_opts, bandwidth_hz_px=250.0, echoes=4, polarity=polarity)

    with pytest.raises(sc.ConfigurationError, match='echo_spacing_s') as caught:
        _line(fast_opts, bandwidth_hz_px=250.0, echoes=4, polarity=polarity,
              echo_spacing_s=train.min_echo_spacing_s - 100e-6)
    assert f'{train.min_echo_spacing_s:.6g}' in str(caught.value)


# ---------------------------------------------------------------------------- the labels
@pytest.mark.parametrize(('polarity', 'reversed_'), [('monopolar', [0] * 8),
                                                     ('bipolar', [0, 1] * 4)])
def test_eco_and_rev_are_emitted_on_every_acquired_echo(fast_opts, polarity, reversed_) -> None:
    """
    Both labels, unconditionally, including echo 0 and including monopolar -- because pulseq
    labels are **stateful**.  A skipped ``ECO`` on echo 0 leaves it wearing the previous
    repetition's last value; an omitted ``REV`` under monopolar inherits whatever an EPI train, a
    bipolar train or a navigator left standing, and a monopolar train that inherits ``REV = 1`` is
    mirrored by the reconstruction -- which a symmetric phantom does not reveal.
    """
    train = _line(fast_opts, bandwidth_hz_px=250.0, echoes=8, polarity=polarity)
    labels = [(e.label, e.value) for e in _events(train(), 'labelset')]

    assert [v for name, v in labels if name == 'ECO'] == list(range(8))
    assert [v for name, v in labels if name == 'REV'] == reversed_


def test_the_labels_are_read_back_off_the_compiled_sequence(fast_opts) -> None:
    """Off the file rather than off the tree, which is what a reconstruction actually sees."""
    train = _line(fast_opts, bandwidth_hz_px=250.0, echoes=4, polarity='bipolar')
    seq = sc.compile(train(), fast_opts)

    seen = []
    for index in range(1, len(seq.block_events) + 1):
        label = getattr(seq.get_block(index), 'label', None)
        if label is None:
            continue
        events = label if isinstance(label, (list, tuple)) else [label]
        seen.append({e.label: e.value for e in events})

    assert [d.get('ECO') for d in seen] == [0, 1, 2, 3]
    assert [d.get('REV') for d in seen] == [0, 1, 0, 1]


def test_a_single_echo_emits_no_labels(fast_opts) -> None:
    """Which is what keeps ``echoes=1`` byte-identical: two extra events are not identical."""
    assert not _events(_line(fast_opts, bandwidth_hz_px=250.0)(), 'labelset')


@pytest.mark.parametrize('polarity', ['monopolar', 'bipolar'])
def test_acquire_false_drops_the_adc_and_both_labels(fast_opts, polarity) -> None:
    """
    Same duration, same gradients, no ADC and **no echo index** -- for the reason ``GRE2DTR``
    suppresses ``LIN``: a dummy that emitted one would put an echo index in the file for a readout
    that was never sampled.
    """
    train = _line(fast_opts, bandwidth_hz_px=250.0, echoes=6, polarity=polarity)
    live, dummy = train(), train(acquire=False)
    hashes = lambda tree, kind: [sc.events.content_hash(e) for e in _events(tree, kind)]  # noqa: E731

    assert dummy.duration == pytest.approx(live.duration, abs=1e-15)
    assert hashes(dummy, 'trap') == hashes(live, 'trap')
    assert hashes(dummy, 'adc') == []
    assert hashes(dummy, 'labelset') == []
    assert len(hashes(live, 'adc')) == 6


# --------------------------------------------------------------------------- the contract
@pytest.mark.parametrize('polarity', ['monopolar', 'bipolar'])
def test_the_train_is_pure_and_compiles_alone(fast_opts, polarity, component_checks) -> None:
    """Two calls, one train: the readout gradient is reused ``echoes`` times and never mutated."""
    component_checks.all(_line(fast_opts, bandwidth_hz_px=250.0, echoes=4, polarity=polarity))
