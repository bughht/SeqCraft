"""
``WaveReadout`` -- the corkscrew's contract, measured on what is emitted and sampled.

Every k below is integrated from emitted or compiled gradients and evaluated at the actual ADC
sample times, ``adc.delay + (n + 1/2) * dwell``.  The analytic Wave is an oracle for the shape,
never for the contract.

The reference protocol is the sagittal Wave-MPRAGE/Wave-GRE one: readout z, sine y, cosine x,
250 pixels at 4x oversampling, a 5 us dwell and a 5 ms window, at the reference amplitudes for 10,
20 and 25 cycles.  The 25-cycle case is the odd one, where a centring computed for a different
realisation leaves the cosine axis off target at the echo.
"""

from __future__ import annotations

import warnings

import numpy as np
import pypulseq as pp
import pytest
from pypulseq.opts import Opts

import seqcraft as sc
from seqcraft.design.events import check_limits, knots_of, pwl_moment
from seqcraft.modules.readout.wave_readout import _integral_at

#: (cycles, amplitude mT/m): the coupled reference cases, odd and even.
REFERENCE = [(10, 12.732), (20, 6.3662), (25, 5.093)]


@pytest.fixture(scope='module')
def scanner() -> Opts:
    """The reference protocols' 80 mT/m, 200 T/m/s scanner, with stated dead times."""
    return Opts(max_grad=80, grad_unit='mT/m', max_slew=200, slew_unit='T/m/s', B0=3.0,
                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=20e-6)


def wave_readout(opts: Opts, **changes) -> sc.modules.WaveReadout:
    """The sagittal reference readout, with `changes` applied."""
    kwargs = dict(opts=opts, fov_mm=250.0, matrix=250, dwell_s=5e-6, cycles=25,
                  amplitude_mT_m=5.093, axis='z', sine_axis='y', cosine_axis='x')
    return sc.modules.WaveReadout(**{**kwargs, **changes})


def compiled_k(wave: sc.modules.WaveReadout, opts: Opts) -> tuple[np.ndarray, np.ndarray]:
    """
    ``(t_adc, k_adc)`` of the block compiled alone, with the stated moment added when the block
    leaves it to its caller -- which is exactly what a caller is required to supply.
    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = sc.kspace(sc.LogicBlock('probe').add(0.0, wave()), opts)
    k = out['k_adc'].copy()
    if not wave.prephase:
        for row, axis in enumerate('xyz'):
            k[row] += wave.prephaser_area_per_m(axis)
    return out['t_adc'], k


def wave_event(wave: sc.modules.WaveReadout, axis: str) -> tuple[np.ndarray, np.ndarray]:
    """The emitted knots of the Wave on `axis`, on the block clock."""
    found = [knots_of(event, start) for start, event, _ in sc.flatten(wave())
             if getattr(event, 'channel', None) == axis and event.type == 'grad']
    assert len(found) == 1
    return found[0]


def adc_start(wave: sc.modules.WaveReadout) -> float:
    return float(wave.sample_times_s()[0] - wave.line.dwell_s / 2)


# ------------------------------------------------------------------ the echo contract
@pytest.mark.parametrize('prephase', [True, False])
@pytest.mark.parametrize(('cycles', 'amplitude'), [*REFERENCE, (7, 'max'), (8, 'max')])
def test_sinzero_puts_the_echo_sample_on_target_on_both_wave_axes(
        scanner, cycles, amplitude, prephase) -> None:
    """
    For odd and even cycles alike, measured on the compiled block at the true sample times.
    At 25 cycles this is the regression case: the cosine axis closes to zero, not to ``G0 * dt``.
    """
    wave = wave_readout(scanner, cycles=cycles, amplitude_mT_m=amplitude, prephase=prephase)
    _, k = compiled_k(wave, scanner)
    echo = wave.line.echo_sample(0)
    dk = wave.line.dk_per_m

    assert np.max(np.abs(k[:, echo])) < 1e-6 * dk
    assert np.max(np.abs(wave.k_per_m()[:, echo])) < 1e-9 * dk


@pytest.mark.parametrize('prephase', [True, False])
@pytest.mark.parametrize(('cycles', 'amplitude'), [*REFERENCE, (7, 'max')])
def test_sinctr_centres_the_sine_corkscrew_and_keeps_the_cosine_echo_on_target(
        scanner, cycles, amplitude, prephase) -> None:
    wave = wave_readout(scanner, cycles=cycles, amplitude_mT_m=amplitude, centering='sinctr',
                        prephase=prephase)
    _, k = compiled_k(wave, scanner)
    echo = wave.line.echo_sample(0)
    dk = wave.line.dk_per_m

    assert abs(np.mean(k[1])) < 1e-6 * dk                        # sine: centred on the line
    assert abs(k[0, echo]) < 1e-6 * dk                           # cosine: echo on target
    assert k[1, echo] == pytest.approx(wave.echo_offset_per_m('y'), abs=1e-6 * dk)
    assert wave.echo_offset_per_m('x') == wave.echo_offset_per_m('z') == 0.0


@pytest.mark.parametrize(('cycles', 'amplitude'), REFERENCE)
def test_the_sinctr_offset_is_the_closed_form_for_the_emitted_wave(scanner, cycles,
                                                                   amplitude) -> None:
    """
    The analytic oracle.  For the continuous Wave, k at the echo minus the corkscrew's centre is
    ``(-1)**(cycles + 1) * G0/omega * cos(omega * dwell / 2)``.  The emitted Wave interpolates the
    sinusoid linearly between gradient-raster edges, which scales its oscillation by
    ``1 - (omega * dt)**2 / 12`` to leading order, and the module reports that.
    """
    wave = wave_readout(scanner, cycles=cycles, amplitude_mT_m=amplitude, centering='sinctr')
    g0 = amplitude * 1e-3 * float(scanner.gamma)
    omega = wave.angular_frequency_rad_s
    continuous = (-1) ** (cycles + 1) * g0 / omega * np.cos(omega * wave.line.dwell_s / 2)
    interpolated = 1 - (omega * float(scanner.grad_raster_time)) ** 2 / 12
    assert wave.echo_offset_per_m('y') == pytest.approx(continuous * interpolated, rel=5e-4)


@pytest.mark.parametrize('prephase', [True, False])
def test_centering_never_changes_an_active_sample(scanner, prephase) -> None:
    zero = wave_readout(scanner, prephase=prephase)
    centred = wave_readout(scanner, prephase=prephase, centering='sinctr')
    start, end = adc_start(zero), adc_start(zero) + zero.window_s
    for axis in ('y', 'x'):
        (t0, g0), (t1, g1) = wave_event(zero, axis), wave_event(centred, axis)
        inside0, inside1 = (t0 >= start) & (t0 <= end), (t1 >= start) & (t1 <= end)
        np.testing.assert_array_equal(t0[inside0], t1[inside1])
        np.testing.assert_array_equal(g0[inside0], g1[inside1])


# --------------------------------------------------------------------- the timing anchor
@pytest.mark.parametrize(('matrix', 'oversampling', 'dwell_s', 'cycles'), [
    (250, 4, 5e-6, 25),            # the reference: a 5 ms window, a knot on every period edge
    (96, 4, 6.5e-6, 7),            # a 2496 us window: it ends between two gradient rasters
    (100, 3, 3.3e-6, 7),           # a 141.43 us period: not a whole number of rasters
])
def test_the_wave_is_anchored_to_the_adc_dwell_window(scanner, matrix, oversampling, dwell_s,
                                                      cycles) -> None:
    """
    Every emitted knot inside the window is the continuous prescription sampled there: phase
    anchored at the ADC start, and exactly `cycles` periods over ``num_samples * dwell_s``.
    """
    wave = wave_readout(scanner, matrix=matrix, oversampling=oversampling, dwell_s=dwell_s,
                        cycles=cycles, amplitude_mT_m='max')
    window = wave.line.num_samples * wave.line.dwell_s
    assert wave.window_s == pytest.approx(window, abs=1e-15)
    omega = 2 * np.pi * cycles / window
    assert wave.angular_frequency_rad_s == pytest.approx(omega, rel=1e-12)

    g0 = wave.amplitude_mT_m * 1e-3 * float(scanner.gamma)
    start = adc_start(wave)
    for axis, shape in (('y', np.sin), ('x', np.cos)):
        t, g = wave_event(wave, axis)
        inside = (t >= start - 1e-12) & (t < start + window - 1e-12)
        assert inside.sum() > 0
        np.testing.assert_allclose(g[inside], g0 * shape(omega * (t[inside] - start)),
                                   atol=1e-9 * g0)


def test_sample_times_are_the_compiled_adc_times(scanner) -> None:
    wave = wave_readout(scanner, prephase=False)
    t_adc, _ = compiled_k(wave, scanner)
    np.testing.assert_allclose(t_adc, wave.sample_times_s(), atol=1e-12)
    assert wave.time_to_echo() == pytest.approx(t_adc[wave.line.echo_sample(0)], abs=1e-12)
    assert wave.time_to_echo() == pytest.approx(
        float(wave.sample_times_s()[0]) + wave.line.echo_sample(0) * wave.line.dwell_s,
        abs=1e-12)


# ------------------------------------------------- independent ADC and gradient rasters
@pytest.fixture(scope='module')
def late_adc() -> Opts:
    """An ADC dead time of 25 us: the ADC opens 25 us into a 20 us-rise lobe, between rasters."""
    return Opts(max_grad=80, grad_unit='mT/m', max_slew=200, slew_unit='T/m/s', B0=3.0,
                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=25e-6)


#: Independently legal rasters whose sample points do not coincide: dwells of 6.5, 3.3 and 5 us
#: against a 10 us gradient raster; windows that end between rasters; periods that are not whole
#: rasters; and, with `late_adc`, a window that starts between rasters too.
NON_COMMON = [
    ('scanner', dict(matrix=96, oversampling=4, dwell_s=6.5e-6, cycles=7)),
    ('scanner', dict(matrix=100, oversampling=3, dwell_s=3.3e-6, cycles=9)),
    ('scanner', dict(matrix=77, oversampling=4, dwell_s=4.1e-6, cycles=4)),   # echo off N/2
    ('late_adc', dict(matrix=64, oversampling=4, dwell_s=5e-6, cycles=5, fov_mm=300.0)),
]


@pytest.mark.parametrize('centering', ['sinzero', 'sinctr'])
@pytest.mark.parametrize('prephase', [True, False])
@pytest.mark.parametrize(('which', 'changes'), NON_COMMON)
def test_non_common_rasters_still_meet_the_contract(request, which, changes, prephase,
                                                    centering) -> None:
    opts = request.getfixturevalue(which)
    wave = wave_readout(opts, amplitude_mT_m='max', prephase=prephase, centering=centering,
                        **changes)
    t_adc, k = compiled_k(wave, opts)
    echo = wave.line.echo_sample(0)
    dk = wave.line.dk_per_m
    raster = float(opts.grad_raster_time)

    # The samples are not on the gradient raster ...
    phase = (t_adc / raster) % 1.0
    assert np.min(np.minimum(phase, 1.0 - phase)) > 1e-6
    # ... and the emitted Wave still covers the window and meets the contract on its samples.
    start = adc_start(wave)
    for axis in ('y', 'x'):
        t, _ = wave_event(wave, axis)
        assert t[0] <= start + 1e-12 and t[-1] >= start + wave.window_s - 1e-12
    assert abs(k[0, echo]) < 1e-6 * dk
    assert abs(np.mean(k[1]) if centering == 'sinctr' else k[1, echo]) < 1e-6 * dk
    np.testing.assert_allclose(k, wave.k_per_m(), atol=1e-5 * dk)


@pytest.mark.parametrize(('which', 'changes'), NON_COMMON)
def test_the_sampled_corkscrew_is_the_prescribed_one_to_within_the_raster(request, which,
                                                                          changes) -> None:
    """
    The rasterised Wave is an approximation of the prescription, and the size of the difference
    is predictable: linear interpolation between raster edges scales a sinusoid by
    ``1 - (omega * dt)**2 / 12``.  Aligned at the echo, the analytic trajectory moves up to
    ``2 * G0 / omega`` from its echo value, so the sampled one stays within 1.2 times that
    fraction of ``2 * G0 / omega`` -- including where the window starts or ends between rasters.
    """
    opts = request.getfixturevalue(which)
    wave = wave_readout(opts, amplitude_mT_m='max', **changes)
    omega = wave.angular_frequency_rad_s
    g0 = wave.amplitude_mT_m * 1e-3 * float(opts.gamma)
    since = wave.sample_times_s() - adc_start(wave)
    analytic = {'y': g0 / omega * (1 - np.cos(omega * since)),
                'x': g0 / omega * np.sin(omega * since)}
    k = wave.k_per_m()
    echo = wave.line.echo_sample(0)
    attenuation = (omega * float(opts.grad_raster_time)) ** 2 / 12
    for row, axis in ((1, 'y'), (0, 'x')):
        oracle = analytic[axis] - analytic[axis][echo]
        assert np.max(np.abs(k[row] - oracle)) < 1.2 * attenuation * 2 * g0 / omega


@pytest.mark.parametrize('prephase', [True, False])
def test_an_off_raster_adc_start_anchors_the_prescription_not_the_emitted_value(
        late_adc, tmp_path, prephase) -> None:
    """
    With the ADC opening between gradient-raster edges, the analytic phase anchor is still the ADC
    start: every emitted knot inside the window is ``G0 * sin|cos(omega * (t - t_a))``.  The emitted
    waveform is that prescription's interpolant between raster edges, so **at** ``t_a`` it takes the
    chord between the two surrounding edges.  Where the waveform before the window is the module's
    own fixed one (``prephase=False``), that chord is measurably neither the sine's ``0`` nor the
    cosine's ``G0``.  The contract is met by the trajectory at the actual ADC samples, which is what
    the compiled and the written-and-read-back sequence play.
    """
    wave = wave_readout(late_adc, matrix=64, oversampling=4, dwell_s=5e-6, cycles=5, fov_mm=300.0,
                        amplitude_mT_m='max', prephase=prephase)
    raster = float(late_adc.grad_raster_time)
    t_a = adc_start(wave)
    assert 1e-3 < (t_a / raster) % 1.0 < 1 - 1e-3                    # between two edges
    edges = np.array([np.floor(t_a / raster), np.ceil(t_a / raster)]) * raster
    g0 = wave.amplitude_mT_m * 1e-3 * float(late_adc.gamma)
    omega = wave.angular_frequency_rad_s

    for axis, shape, endpoint in (('y', np.sin, 0.0), ('x', np.cos, g0)):
        t, g = wave_event(wave, axis)
        inside = (t > t_a) & (t < t_a + wave.window_s)
        np.testing.assert_allclose(g[inside], g0 * shape(omega * (t[inside] - t_a)),
                                   atol=1e-9 * g0)                      # anchored at t_a
        emitted = np.interp(t_a, t, g)
        chord = np.interp(t_a, edges, np.interp(edges, t, g))
        assert emitted == pytest.approx(chord, abs=1e-9 * g0)          # rasterised: the chord
        if not prephase:
            assert abs(emitted - endpoint) > 1e-2 * g0                 # not the analytic value

    t_adc, k = compiled_k(wave, late_adc)
    dk = wave.line.dk_per_m
    assert np.max(np.abs(k[:, wave.line.echo_sample(0)])) < 1e-6 * dk
    np.testing.assert_allclose(k, wave.k_per_m(), atol=1e-5 * dk)
    np.testing.assert_allclose(t_adc, wave.sample_times_s(), atol=1e-12)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        seq = sc.compile(sc.LogicBlock('probe').add(0.0, wave()), late_adc)
        seq.write(str(tmp_path / 'late.seq'))
        back = pp.Sequence(system=late_adc)
        back.read(str(tmp_path / 'late.seq'))
        k_back = sc.sequence_kspace(back)['k_adc'].copy()
    if not prephase:
        k_back += np.array([wave.prephaser_area_per_m(a) for a in 'xyz'])[:, None]
    module = wave.k_per_m()
    np.testing.assert_allclose(k_back, module, atol=1e-5 * np.max(np.abs(module)))


# -------------------------------------------------------------- transitions and prephase
@pytest.mark.parametrize('prephase', [True, False])
def test_the_cosine_ramps_are_emitted_whatever_prephase_is(scanner, prephase) -> None:
    wave = wave_readout(scanner, prephase=prephase)
    g0 = wave.amplitude_mT_m * 1e-3 * float(scanner.gamma)
    start, end = adc_start(wave), adc_start(wave) + wave.window_s
    t, g = wave_event(wave, 'x')

    assert g[0] == 0.0 and g[-1] == 0.0                          # from rest, to rest
    assert np.interp(start, t, g) == pytest.approx(g0, rel=1e-12)    # at G0 when sampling starts
    assert np.interp(end, t, g) == pytest.approx(g0, rel=1e-12)      # and when it stops
    assert np.all(np.abs(g[t > end]) < g0)                       # then it ramps straight down


def test_the_exit_ramp_starts_when_the_window_closes(scanner) -> None:
    """No hold through the ADC dead time: the shortest legal ramp, from the window's end."""
    wave = wave_readout(scanner, prephase=False)
    g0 = wave.amplitude_mT_m * 1e-3 * float(scanner.gamma)
    end = adc_start(wave) + wave.window_s
    t, g = wave_event(wave, 'x')
    ramp = np.ceil(g0 / float(scanner.max_slew) / float(scanner.grad_raster_time) - 1e-9) \
        * float(scanner.grad_raster_time)
    assert t[-1] == pytest.approx(end + ramp, abs=1e-12)


def test_prephase_false_states_exactly_what_prephase_true_plays(scanner) -> None:
    """
    The two give one trajectory on the samples: what one block realises, the other states, and
    the caller then supplies.
    """
    plays = wave_readout(scanner, prephase=True, centering='sinctr')
    states = wave_readout(scanner, prephase=False, centering='sinctr')
    np.testing.assert_allclose(plays.k_per_m(), states.k_per_m(), atol=1e-9 * plays.line.dk_per_m)
    for axis in 'xyz':
        assert plays.prephaser_area_per_m(axis) == pytest.approx(
            states.prephaser_area_per_m(axis), abs=1e-9)


def test_prephase_true_merges_the_cosine_moment_and_entrance_into_one_event(scanner) -> None:
    """One event per Wave axis; the cosine's leaves rest at the block's start and is at G0 by t_a."""
    wave = wave_readout(scanner, prephase=True)
    t, g = wave_event(wave, 'x')
    assert t[0] == 0.0
    assert np.interp(adc_start(wave), t, g) == pytest.approx(
        wave.amplitude_mT_m * 1e-3 * float(scanner.gamma), rel=1e-12)
    events = [e for _, e, _ in sc.flatten(wave()) if getattr(e, 'channel', None) == 'x']
    assert len(events) == 1


def test_the_block_leaves_k_where_area_after_echo_says(scanner) -> None:
    wave = wave_readout(scanner, prephase=True, centering='sinctr')
    for axis in 'xyz':
        total = sum(pwl_moment(*knots_of(e, s))
                    for s, e, _ in sc.flatten(wave()) if getattr(e, 'channel', None) == axis)
        assert total == pytest.approx(
            wave.echo_offset_per_m(axis) + wave.area_after_echo_per_m(axis), abs=1e-6)


@pytest.mark.parametrize('flow', [False, True])
@pytest.mark.parametrize('centering', ['sinzero', 'sinctr'])
def test_a_scope_composition_delivers_the_stated_moment(scanner, centering, flow) -> None:
    """
    ``prephase=False`` under :class:`~seqcraft.PhysicalDesignScope`: slab rephasing, the LIN/PAR
    encoding, the Wave's own requirement and optionally flow compensation, merged into one lobe
    per axis, measured on the compiled sequence over a centre and a corner state.
    """
    wave = wave_readout(scanner, fov_mm=200.0, matrix=64, dwell_s=7.8e-6, cycles=7,
                        amplitude_mT_m='max', prephase=False, centering=centering)
    exc = sc.modules.Excitation(opts=scanner, flip_deg=15.0, thickness_mm=150.0, axis='x',
                                duration_s=1e-3)
    dk = {'y': 1e3 / 220.0, 'x': 1e3 / 150.0}
    states = ((16, 12), (0, 0), (31, 23))

    def encoded(state, axis):
        return {'y': (state[0] - 16) * dk['y'], 'x': (state[1] - 12) * dk['x'], 'z': 0.0}[axis]

    scope = sc.PhysicalDesignScope(
        origin_s=exc.time_to_center(), before=exc(rephase=False), after=wave(),
        echo_in_after_s=wave.time_to_echo(), axes=('x', 'y', 'z'), states=states,
        k_at_echo=lambda state, axis: encoded(state, axis) + wave.echo_offset_per_m(axis))
    design = scope.design(opts=scanner,
                          flow_comp=sc.FlowCompensation(axis=('x', 'y', 'z')) if flow else None)
    echo = wave.line.echo_sample(0)
    for state in states:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            k = sc.kspace(design.build(state), scanner)['k_adc']
        assert abs(k[2, echo]) < 1e-6 * wave.line.dk_per_m
        assert abs(k[0, echo] - encoded(state, 'x')) < 1e-6 * dk['x']
        sine = np.mean(k[1]) if centering == 'sinctr' else k[1, echo]
        assert abs(sine - encoded(state, 'y')) < 1e-6 * dk['y']


# ------------------------------------------------------------------------- amplitude
def test_max_is_the_slew_bound_when_the_slew_binds(scanner) -> None:
    wave = wave_readout(scanner, cycles=25, amplitude_mT_m='max')
    g0 = wave.amplitude_mT_m * 1e-3 * float(scanner.gamma)
    assert g0 * wave.angular_frequency_rad_s == pytest.approx(float(scanner.max_slew), rel=1e-12)
    assert check_limits([e for _, e, _ in sc.flatten(wave())], scanner) == []


def test_max_is_the_gradient_bound_when_the_gradient_binds() -> None:
    weak = Opts(max_grad=6, grad_unit='mT/m', max_slew=200, slew_unit='T/m/s', B0=3.0,
                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=20e-6)
    wave = wave_readout(weak, cycles=2, amplitude_mT_m='max')
    assert wave.amplitude_mT_m == pytest.approx(6.0, rel=1e-9)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sc.compile(sc.LogicBlock('probe').add(0.0, wave()), weak)


def test_a_numeric_amplitude_is_played_exactly(scanner) -> None:
    wave = wave_readout(scanner, amplitude_mT_m=5.093, prephase=False)
    assert wave.amplitude_mT_m == pytest.approx(5.093, rel=1e-12)
    t, g = wave_event(wave, 'x')
    assert np.interp(adc_start(wave), t, g) == pytest.approx(
        5.093e-3 * float(scanner.gamma), rel=1e-12)


@pytest.mark.parametrize(('cycles', 'amplitude', 'bound'), [(25, 6.5, 'max_slew'),
                                                            (1, 90.0, 'max_grad')])
def test_an_infeasible_amplitude_is_refused_not_reduced(scanner, cycles, amplitude,
                                                        bound) -> None:
    with pytest.raises(sc.ConfigurationError) as raised:
        wave_readout(scanner, cycles=cycles, amplitude_mT_m=amplitude)
    message = str(raised.value)
    assert bound in message and 'the largest that fits is' in message


# ----------------------------------------------------------------------- oversampling
@pytest.mark.parametrize('oversampling', [1, 2, 3, 4, 5])
def test_any_oversampling_samples_the_same_contract(scanner, oversampling) -> None:
    wave = wave_readout(scanner, matrix=100, dwell_s=20e-6 / oversampling,
                        oversampling=oversampling, cycles=9, amplitude_mT_m='max')
    k = wave.k_per_m()
    echo = wave.line.echo_sample(0)
    assert wave.line.num_samples == 100 * oversampling
    assert echo == wave.line.num_samples // 2
    np.testing.assert_allclose(np.diff(k[2]), wave.line.dk_per_m / oversampling, rtol=1e-9)
    assert np.max(np.abs(k[:2, echo])) < 1e-9 * wave.line.dk_per_m


def test_the_default_oversampling_is_four(scanner) -> None:
    wave = sc.modules.WaveReadout(opts=scanner, fov_mm=250.0, matrix=250, dwell_s=5e-6,
                                  cycles=10)
    assert wave.line.oversampling == 4 and wave.line.num_samples == 1000


# --------------------------------------------------------------------- the written file
@pytest.mark.parametrize('prephase', [True, False])
@pytest.mark.parametrize(('cycles', 'changes'), [(25, {}), (7, NON_COMMON[0][1])])
def test_the_written_file_plays_the_reported_trajectory(scanner, tmp_path, cycles, changes,
                                                        prephase) -> None:
    """
    Write, read back, and integrate the **read-back** gradients exactly at the **read-back** ADC
    times: the module's trajectory is the file's, sample for sample.
    """
    wave = wave_readout(scanner, amplitude_mT_m='max', prephase=prephase,
                        **{**changes, 'cycles': cycles})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        seq = sc.compile(sc.LogicBlock('probe').add(0.0, wave()), scanner)
        path = tmp_path / 'wave.seq'
        seq.write(str(path))
        back = pp.Sequence(system=scanner)
        back.read(str(path))
        trajectory = sc.sequence_kspace(back)
        k_file, t_adc = trajectory['k_adc'], trajectory['t_adc']

    starts, clock = {}, 0.0
    for block in back.block_events:
        starts[block] = clock
        clock += back.block_durations[block]
    k = np.zeros((3, len(t_adc)))
    for block in back.block_events:
        for row, axis in enumerate('xyz'):
            event = getattr(back.get_block(block), f'g{axis}', None)
            if event is not None:
                k[row] += _integral_at(np.asarray(t_adc), *knots_of(event, starts[block]))
    k_file = k_file.copy()
    if not prephase:
        stated = np.array([wave.prephaser_area_per_m(axis) for axis in 'xyz'])[:, None]
        k += stated
        k_file += stated

    # Writing rounds every amplitude to six significant digits, up to 5e-6 of each event, so the
    # file's trajectory is held to twice that fraction of the largest excursion.
    module = wave.k_per_m()
    rounding = 1e-5 * np.max(np.abs(module))
    np.testing.assert_allclose(np.asarray(t_adc), wave.sample_times_s(), atol=1e-9)
    np.testing.assert_allclose(k, module, atol=rounding)
    np.testing.assert_allclose(k_file, module, atol=rounding)


def test_it_compiles_late_in_a_sequence(scanner) -> None:
    """The per-block moment checks, with the block a second into the sequence."""
    wave = wave_readout(scanner)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sc.compile(sc.LogicBlock('late').add(1.0, wave()), scanner)


# ----------------------------------------------------------------------------- shapes
@pytest.mark.parametrize('missing', ['sine_axis', 'cosine_axis'])
def test_a_single_axis_wave_meets_its_contract(scanner, missing) -> None:
    wave = wave_readout(scanner, **{missing: None})
    k = wave.k_per_m()
    idle = 'xyz'.index('y' if missing == 'sine_axis' else 'x')
    assert not np.any(k[idle])
    assert np.max(np.abs(k[:, wave.line.echo_sample(0)])) < 1e-9 * wave.line.dk_per_m


def test_acquire_false_plays_the_same_gradients(scanner) -> None:
    wave = wave_readout(scanner)
    real = [e for _, e, _ in sc.flatten(wave()) if getattr(e, 'type', None) in ('trap', 'grad')]
    dummy = [e for _, e, _ in sc.flatten(wave(acquire=False))]
    assert not [e for e in dummy if e.type == 'adc']
    assert [sc.events.content_hash(e) for e in real] == [
        sc.events.content_hash(e) for e in dummy if e.type in ('trap', 'grad')]


@pytest.mark.parametrize('prephase', [True, False])
def test_it_is_pure_and_compiles_alone(scanner, component_checks, prephase) -> None:
    component_checks.all(wave_readout(scanner, prephase=prephase))


# ---------------------------------------------------------------------------- refusals
@pytest.mark.parametrize('changes', [
    dict(sine_axis='z'),                                   # on the readout axis
    dict(sine_axis='x'),                                   # on the cosine axis
    dict(sine_axis=None, cosine_axis=None),                # no Wave at all
    dict(centering='sinctr', sine_axis=None),              # nothing to centre
    dict(centering='centred'),
    dict(cycles=0),
    dict(amplitude_mT_m='min'),
    dict(cycles=300),                                      # a 16.7 us period on a 10 us raster
    dict(prephase=False, prephaser_duration_s=1e-3),
])
def test_what_it_refuses(scanner, changes) -> None:
    with pytest.raises(sc.ConfigurationError):
        wave_readout(scanner, **changes)


def test_the_wave_lobes_can_set_the_preparation_window(scanner) -> None:
    """
    Two cycles at 40.7 mT/m need a 600 us window, where the readout's own prephaser needs 260:
    the readout prephaser is lengthened to match, and a shorter requested window is refused.
    """
    changes = dict(axis='x', sine_axis='y', cosine_axis='z', matrix=64, dwell_s=10e-6,
                   cycles=2, amplitude_mT_m='max')
    alone = sc.modules.CartesianLine(opts=scanner, fov_mm=250.0, matrix=64, dwell_s=10e-6,
                                     oversampling=4)
    wave = wave_readout(scanner, **changes)
    assert wave.prephaser_duration_s > alone.prephaser_duration_s
    assert wave.line.prephaser_area_per_m == pytest.approx(alone.prephaser_area_per_m, rel=1e-9)

    with pytest.raises(sc.ConfigurationError) as raised:
        wave_readout(scanner, **changes, prephaser_duration_s=4e-4)
    assert f'prephaser_duration_s={wave.prephaser_duration_s:.6g} is the shortest' in str(
        raised.value)
    longer = wave_readout(scanner, **changes,
                          prephaser_duration_s=wave.prephaser_duration_s + 1e-4)
    assert longer.prephaser_duration_s == pytest.approx(wave.prephaser_duration_s + 1e-4,
                                                        abs=1e-12)
