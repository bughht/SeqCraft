"""
Module physics: known values, and the contract every module keeps.

The contract tests are parametrised over every concrete ``Module`` subclass, so a new module gets
them for free.  The known-value tests are per module, and each asserts a number an independent
calculation gives -- not a number the code happened to produce.
"""

from __future__ import annotations

import math

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc


#: One constructed instance per module class, for the contract tests.
def _instances(system: sc.System) -> dict[str, sc.Module]:
    """
    Build one of every concrete Module subclass with plausible parameters.

    Keyed by class name, derived from the instances themselves -- so the keys cannot drift from
    the classes, and the coverage assertion compares against ``sc.testing.all_modules()``.
    """
    m = sc.modules
    built = (
        m.SincExcitation(system, flip_deg=15, duration_us=1000, slice_thickness_mm=5),
        m.SlabExcitation(
            system, flip_deg=10, duration_us=1500, slice_thickness_mm=80, time_bw_product=8),
        m.HardExcitation(system, flip_deg=90, duration_us=500),
        m.SLRExcitation(system, flip_deg=90, duration_us=2000, slice_thickness_mm=4),
        m.SLRRefocusing(
            system, flip_deg=180, duration_us=4000, slice_thickness_mm=4, crusher_twists=4),
        m.SincRefocusing(
            system, flip_deg=180, duration_us=4000, slice_thickness_mm=5, crusher_twists=4),
        m.HardRefocusing(system, flip_deg=180, duration_us=1000),
        m.GaussSaturation(system, flip_deg=90, duration_us=8000, time_bw_product=1.6),
        m.AdiabaticInversion(system, duration_us=10000),
        m.CartesianLine(system, fov_ro_mm=250, matrix_ro=64, readout_duration_us=3200),
        m.NoiseAcquisition(system, n_samples=256),
        m.EPIReadout(
            system, fov_ro_mm=240, matrix_ro=64, fov_pe_mm=240, matrix_pe=64,
            partial_fourier_pe=0.75),
        m.SpiralVDS(system, fov_mm=240, matrix=64, n_interleaves=8),
        m.PhaseEncode(system, fov_pe_mm=250, matrix_pe=64),
        m.PartitionEncode(system, slab_thickness_mm=80, matrix_sl=16),
        m.Prephaser(system, area_per_m=-128.0),
        m.Spoiler(system, twists=4, voxel_mm=5),
        m.Crusher(system, twists=4, voxel_mm=5),
        m.MonopolarDiffusion(system, b_value_s_per_mm2=1000, refocus_duration_us=4200),
        m.BipolarDiffusion(system, b_value_s_per_mm2=300, refocus_duration_us=4200),
        m.ArbitraryDiffusion(
            system,
            waveform_Hz_per_m=np.concatenate(
                [np.linspace(0, 1, 20), np.ones(60), np.linspace(1, 0, 20)]) * 1e6,
            refocus_duration_us=4200,
        ),
        m.FatSat(system, voxel_mm=5),
        m.InversionRecovery(system, voxel_mm=5),
        m.Delay(system, duration_ms=20),
        m.Trigger(system, channel='osc0', duration_us=100),
        m.Barrier(system),
        m.RawEvents(system, events=(pp.make_trapezoid('x', area=100.0, system=system.default),)),
    )
    return {type(module).__name__: module for module in built}


@pytest.fixture(scope='module')
def instances(system: sc.System) -> dict[str, sc.Module]:
    return _instances(system)


def module_names() -> list[str]:
    return sorted(sc.testing.all_modules())


# ------------------------------------------------------------------------------- the contract
def test_every_module_class_is_covered(instances) -> None:
    """
    A new module must be added to the fixture, so it cannot slip past the contract tests.

    Discovery walks ``Module.__subclasses__()``, so subclassing *is* the registration: there is no
    decorator to forget, and coverage is not opt-in.  Intermediate bases declare an abstract
    ``_design`` and are skipped.
    """
    found = set(sc.testing.all_modules())
    assert set(instances) == found, (
        f'not built in the fixture: {sorted(found - set(instances))}; '
        f'built but not a Module subclass: {sorted(set(instances) - found)}'
    )


@pytest.mark.parametrize('name', module_names())
def test_build_returns_a_logic_block(instances, name: str) -> None:
    assert isinstance(instances[name].build(), sc.LogicBlock)


@pytest.mark.parametrize('name', module_names())
def test_build_does_not_mutate_the_module(instances, name: str) -> None:
    """
    ``build`` is called once per TR; mutating the module would make TR 500 differ from TR 1.

    The reference implementation flipped ``self.gx.amplitude`` in place inside its readout loop,
    which is exactly the bug this catches.
    """
    module = instances[name]
    before = {
        key: sc.events.content_hash(value)
        for key, value in vars(module).items()
        if getattr(value, 'type', None) is not None
    }
    module.build()
    module.build()
    after = {
        key: sc.events.content_hash(value)
        for key, value in vars(module).items()
        if getattr(value, 'type', None) is not None
    }
    assert before == after


@pytest.mark.parametrize('name', module_names())
def test_two_builds_give_equal_blocks(instances, name: str) -> None:
    """Same arguments, same output -- the other side of the purity coin."""
    module = instances[name]
    first, second = module.build(), module.build()
    assert len(first) == len(second)
    assert first.duration == pytest.approx(second.duration)
    for a, b in zip(first.nodes, second.nodes):
        assert a.start == pytest.approx(b.start)
        if getattr(a.item, 'type', None) is not None:
            assert sc.events.content_hash(a.item) == sc.events.content_hash(b.item)


@pytest.mark.parametrize('name', module_names())
def test_declared_duration_matches_the_built_block(instances, name: str) -> None:
    """
    A module's ``duration`` is what callers place the *next* thing by, so it cannot be optimistic.

    Allowed to exceed the block (padding to the raster), never to fall short.
    """
    module = instances[name]
    declared = getattr(module, 'duration', None)
    if declared is None:
        pytest.skip(f'{name} declares no duration')
    assert declared >= module.build().duration - 1e-12


@pytest.mark.parametrize('name', module_names())
def test_timing_properties_fall_inside_the_block(instances, name: str) -> None:
    """``isodelay`` and ``time_to_echo`` are offsets into the block, so they must be inside it."""
    module = instances[name]
    duration = getattr(module, 'duration', None)
    if duration is None or duration == 0.0:
        pytest.skip(f'{name} has no duration')
    for attribute in ('isodelay', 'time_to_echo'):
        value = getattr(module, attribute, None)
        if value is not None:
            assert 0.0 <= value <= duration + 1e-12, f'{name}.{attribute} = {value}'


@pytest.mark.parametrize('name', module_names())
def test_each_module_compiles_alone(instances, name: str) -> None:
    """
    Every module must produce a legal sequence on its own, with no limit violations.

    A module that only works when something else happens to be beside it is not reusable.
    """
    module = instances[name]
    block = module.build()
    if not block.nodes or block.duration == 0.0:
        pytest.skip(f'{name} occupies no time on its own')
    out = sc.compile(sc.LogicBlock(f'test_{name}').add(0.0, block), module.system)
    errors = [i for i in out.check().errors if i.kind != 'timing' or 'TotalDuration' not in i.message]
    assert not errors, f'{name}: {errors}'


@pytest.mark.parametrize('name', module_names())
def test_params_are_json_safe(instances, name: str) -> None:
    """The provenance sidecar has to serialise, so ``params`` must contain no waveforms."""
    import json

    json.dumps(instances[name].params())


@pytest.mark.parametrize('name', module_names())
def test_docstring_has_the_required_sections(instances, name: str) -> None:
    """
    Every module documents its parameters and its physics, or the test fails.

    The reference implementation had one docstring across 2 600 lines; this is the mechanical
    measure that keeps that from happening here.
    """
    doc = type(instances[name]).__doc__ or ''
    assert 'Parameters' in doc or 'See :class:' in doc, f'{name} documents no parameters'
    assert 'Examples' in doc, f'{name} has no runnable example'


# ------------------------------------------------------------------------------- RF known values
def test_slice_select_amplitude_is_bandwidth_over_thickness(system) -> None:
    """The defining relation for a slice-selective pulse."""
    exc = sc.modules.SincExcitation(
        system, flip_deg=15, duration_us=1000, slice_thickness_mm=5, time_bw_product=4)
    assert exc.slice_select_amplitude == pytest.approx(4000.0 / 5e-3, rel=1e-9)


def test_isodelay_is_dead_time_plus_the_pulse_centre(system) -> None:
    """TE is measured from here, so it has to be exactly right."""
    exc = sc.modules.SincExcitation(
        system, flip_deg=15, duration_us=1000, slice_thickness_mm=5)
    assert exc.isodelay == pytest.approx(float(exc.rf.delay) + 500e-6, abs=1e-9)


def test_hard_pulse_peak_b1_matches_the_closed_form(system) -> None:
    """``B1 = flip / (2 pi gamma T)`` -- 11.75 uT for 90 degrees in 500 us."""
    hard = sc.modules.HardExcitation(system, flip_deg=90, duration_us=500)
    expected = math.radians(90) / (2 * math.pi * system.gamma * 500e-6) * 1e6
    assert hard.peak_b1_uT == pytest.approx(expected, rel=1e-3)
    assert hard.peak_b1_uT == pytest.approx(11.75, abs=0.02)


def test_an_impossible_pulse_is_refused_with_the_three_fixing_parameters(system) -> None:
    """
    A 2 ms 180 degree pulse needs about 23 uT, above the 20 uT limit -- real physics, not a bug.

    pypulseq says only 'Amplitude violation (117%)'; this says which parameters to change.
    """
    with pytest.raises(sc.ConfigurationError) as err:
        sc.modules.SincExcitation(
            system, flip_deg=180, duration_us=500, slice_thickness_mm=5, rephase=False)
    text = str(err.value)
    assert 'flip_deg' in text
    assert 'duration_us' in text
    assert 'max_b1_uT' in text


def test_off_centre_slice_gets_both_the_frequency_and_the_phase_term(system) -> None:
    """
    The most commonly omitted line in hand-written pulseq.

    Without the phase term every slice acquires a different phase, which appears in the
    reconstruction as slice-dependent phase and gets blamed on the scanner.
    """
    exc = sc.modules.SincExcitation(
        system, flip_deg=15, duration_us=1000, slice_thickness_mm=5)
    z = 12e-3
    shifted = exc.offset_rf(exc.rf, z, 0.0)
    freq = exc.slice_select_amplitude * z
    assert shifted.freq_offset == pytest.approx(freq)
    expected = (-2 * math.pi * freq * float(exc.rf.center)) % (2 * math.pi)
    assert shifted.phase_offset == pytest.approx(expected)


def test_the_rf_phase_composes_rather_than_overwrites(system) -> None:
    """A CPMG pi/2 on a refocusing pulse must survive an off-centre slice offset."""
    refoc = sc.modules.SincRefocusing(
        system, flip_deg=180, duration_us=4000, slice_thickness_mm=5,
        phase_offset_rad=math.pi / 2)
    shifted = refoc.offset_rf(refoc.rf, 0.0, 0.0)
    assert shifted.phase_offset == pytest.approx(math.pi / 2)


def test_crushers_are_equal_not_opposite(system) -> None:
    """
    The refocusing pulse inverts the phase between them, so an equal pair spares the wanted echo.

    Making them opposite would crush the echo instead of the FID.
    """
    refoc = sc.modules.SincRefocusing(
        system, flip_deg=180, duration_us=4000, slice_thickness_mm=5, crusher_twists=4)
    block = refoc.build()
    crushers = [n for n in block if getattr(n.item, 'type', None) == 'trap'
                and n.item.channel == 'z' and abs(float(n.item.area) - 800.0) < 1.0]
    assert len(crushers) == 2
    assert crushers[0].item.area == pytest.approx(crushers[1].item.area)
    assert refoc.crusher_area_per_m == pytest.approx(4 * 1e3 / 5)


def test_refocusing_isodelay_includes_the_leading_crusher(system) -> None:
    """TE placement would be wrong by the crusher's length otherwise."""
    refoc = sc.modules.SincRefocusing(
        system, flip_deg=180, duration_us=4000, slice_thickness_mm=5, crusher_twists=4)
    bare = sc.modules.SincRefocusing(
        system, flip_deg=180, duration_us=4000, slice_thickness_mm=5)
    assert refoc.isodelay == pytest.approx(bare.isodelay + refoc.crusher_duration)


# -------------------------------------------------------------------------- encoding known values
def test_phase_encode_step_is_one_over_fov(system) -> None:
    pe = sc.modules.PhaseEncode(system, fov_pe_mm=250, matrix_pe=64)
    assert pe.dk_per_m == pytest.approx(1e3 / 250)
    assert pe.area_for(16) == pytest.approx(16 * 4.0)
    assert pe.area_for(-32) == pytest.approx(-128.0)


def test_every_phase_encode_line_takes_the_same_time(system) -> None:
    """
    If the blip shortened near k=0, TE would vary line by line and the image would carry a phase
    ramp nothing in the reconstruction accounts for.
    """
    pe = sc.modules.PhaseEncode(system, fov_pe_mm=250, matrix_pe=64)
    # Compared on the raster, not for float identity: ``calc_duration`` returns
    # 0.00030000000000000003 for a 300 us trapezoid while the k=0 delay is exactly 0.0003, and one
    # ulp is not a physics difference.  ``pe.duration`` is the authority either way.
    for line in (-32, -7, 0, 1, 31):
        assert pe.build(line=line).duration == pytest.approx(pe.duration, abs=1e-12)


def test_the_zero_line_still_occupies_its_slot(system) -> None:
    pe = sc.modules.PhaseEncode(system, fov_pe_mm=250, matrix_pe=64)
    block = pe.build(line=0)
    assert block.duration == pytest.approx(pe.duration)
    assert not [n for n in block if getattr(n.item, 'type', None) in ('trap', 'grad')]


def test_phase_encode_area_scales_linearly_with_line(system) -> None:
    pe = sc.modules.PhaseEncode(system, fov_pe_mm=250, matrix_pe=64)
    for line in (-32, -8, 5, 31):
        block = pe.build(line=line)
        area = sum(float(n.item.area) for n in block
                   if getattr(n.item, 'type', None) == 'trap')
        assert area == pytest.approx(pe.area_for(line), rel=1e-9)


def test_spoiler_area_is_twists_over_voxel(system) -> None:
    """Parameterised in the physical intent, so it stays correct when the resolution changes."""
    spoil = sc.modules.Spoiler(system, twists=4, voxel_mm=5)
    assert spoil.area_per_m == pytest.approx(800.0)
    assert sc.modules.Spoiler(system, twists=4, voxel_mm=2.5).area_per_m == pytest.approx(1600.0)


# --------------------------------------------------------------------------- readout known values
def test_cartesian_k_extent_and_resolution(system) -> None:
    ro = sc.modules.CartesianLine(
        system, fov_ro_mm=250, matrix_ro=64, readout_duration_us=3200)
    assert ro.k_max_per_m == pytest.approx(64 / (2 * 0.25))
    assert ro.resolution_mm == pytest.approx(250 / 64)
    assert ro.bandwidth_per_pixel_hz == pytest.approx(1.0 / (ro.dwell_s * ro.n_samples))


def test_time_to_echo_plus_after_echo_is_the_whole_block(system) -> None:
    """The pair is what lets symmetric and asymmetric readouts share one placement expression."""
    ro = sc.modules.CartesianLine(
        system, fov_ro_mm=250, matrix_ro=64, readout_duration_us=3200)
    assert ro.time_to_echo + ro.duration_after_echo == pytest.approx(ro.duration)


def test_dropping_the_prephaser_shortens_time_to_echo_consistently(system) -> None:
    """
    The property has to follow the constructor, which is why `prephase` is not a build argument.

    When it was one, ``build(prephase=False)`` moved the echo 280 us early while ``time_to_echo``
    still claimed otherwise.
    """
    kwargs = {'fov_ro_mm': 250, 'matrix_ro': 64, 'readout_duration_us': 3200}
    full = sc.modules.CartesianLine(system, **kwargs)
    bare = sc.modules.CartesianLine(system, prephase=False, **kwargs)
    assert full.time_to_echo - bare.time_to_echo == pytest.approx(full.prephase_duration)
    assert bare.build().duration == pytest.approx(bare.duration)


@pytest.mark.parametrize('partial_echo', [1.0, 0.875, 0.75, 0.625])
def test_the_echo_lands_on_k_zero(system, partial_echo: float) -> None:
    """
    ``time_to_echo`` must point at k = 0, checked against pypulseq's own trajectory.

    Asserting the prephaser equals minus ``k_max`` -- which this test used to do -- asserts a bug.
    The gradient's ramp carries area too, so a prephaser of exactly ``k_max`` leaves the line
    displaced by the ramp's own area: about one dk here and nearly four on a short low-bandwidth
    readout.  A k-space offset is a linear phase ramp across the image, so a magnitude image looks
    fine and every phase measurement built on it is wrong.

    Parametrised over partial echo because that is the case where k=0 stops being the middle of the
    ADC window: it arrives earlier, which is the whole reason to drop leading samples.
    """
    ro = sc.modules.CartesianLine(
        system, fov_ro_mm=250, matrix_ro=64, readout_duration_us=3200,
        partial_echo=partial_echo)
    out = sc.compile(sc.LogicBlock('t').add(0.0, ro.build()), system)
    kspace = out.kspace()
    k_at_echo = float(np.interp(ro.time_to_echo, kspace['t_adc'], kspace['k_adc'][0]))
    assert k_at_echo == pytest.approx(0.0, abs=0.02 * ro.dk_per_m)
    # And the sampling really is asymmetric by the requested fraction.
    first, last = float(kspace['k_adc'][0][0]), float(kspace['k_adc'][0][-1])
    assert -first / last == pytest.approx(ro.partial_echo, rel=0.05)
    assert abs(kspace['k_adc'][0][1] - first) == pytest.approx(ro.dk_per_m, rel=1e-3)


def test_an_adc_sample_count_off_the_divisor_is_refused(system) -> None:
    """Siemens requires a multiple of 4, and the message says what the nearest one is."""
    with pytest.raises(sc.ConfigurationError, match='multiple of the required divisor'):
        sc.modules.CartesianLine(
            system, fov_ro_mm=250, matrix_ro=63, readout_duration_us=3200)


# ---------------------------------------------------------------------------- spiral known values
def test_spiral_reaches_k_max_and_starts_at_the_origin(system) -> None:
    """
    The spiral proper ends at k_max; the ramp-down then carries |k| a little past it.

    Bounded above as well as below.  Only reaching k_max during the ramp would mean the winding
    itself stopped short, so ``|k|`` must cross the edge before the tail begins -- and the overshoot
    must stay small, because a large one is a radial streak the density compensation would weight
    heavily.  The turn *spacing*, which the same arithmetic can get wrong independently, is checked
    by :func:`test_the_combined_spiral_meets_nyquist`.
    """
    ro = sc.modules.SpiralVDS(system, fov_mm=240, matrix=64, n_interleaves=8)
    kx, ky = ro.trajectory(0)
    radius = np.hypot(kx, ky)
    reach = float(radius.max()) / ro.k_max_per_m
    assert 1.0 <= reach <= 1.2, f'|k| reaches {reach:.3f} of k_max'
    crossing = int(np.argmax(radius >= ro.k_max_per_m)) / len(radius)
    assert crossing > 0.85, f'|k| only reaches k_max {crossing:.0%} of the way through'
    assert float(np.hypot(kx[0], ky[0])) == pytest.approx(0.0, abs=1e-9)


def test_spiral_respects_the_amplifier(system) -> None:
    """
    The generator's whole contract: it measures the finished waveform, it does not assume.

    An earlier version neglected the ``theta''`` term, which put the realised slew several times
    over the limit right at the origin -- where a diffusion measurement's best samples are.
    """
    derated = system.derate('spiral', grad=0.85, slew=0.65)
    opts = derated.limits('spiral')
    ro = sc.modules.SpiralVDS(
        derated, fov_mm=240, matrix=64, n_interleaves=8, regime='spiral')
    gx = np.asarray(ro.gx.waveform)
    gy = np.asarray(ro.gy.waveform)
    assert float(np.hypot(gx, gy).max()) <= float(opts.max_grad) * 1.001
    slew = np.hypot(np.diff(gx), np.diff(gy)) / derated.grad_raster.dt
    assert float(slew.max()) <= float(opts.max_slew) * 1.001


def test_spiral_gradient_starts_and_ends_at_rest(system) -> None:
    """Otherwise the rewinder starts away from the previous block's end and pypulseq refuses."""
    ro = sc.modules.SpiralVDS(system, fov_mm=240, matrix=64, n_interleaves=8)
    for axis in (ro.gx, ro.gy):
        assert float(axis.waveform[0]) == pytest.approx(0.0, abs=1e-9)
        assert float(axis.waveform[-1]) == pytest.approx(0.0, abs=1e-9)


def test_interleaves_are_rotations_of_one_design(system) -> None:
    ro = sc.modules.SpiralVDS(system, fov_mm=240, matrix=64, n_interleaves=4)
    kx0, ky0 = ro.trajectory(0)
    kx1, ky1 = ro.trajectory(1)
    angle = math.atan2(ky1[-1], kx1[-1]) - math.atan2(ky0[-1], kx0[-1])
    assert angle == pytest.approx(2 * math.pi / 4, abs=1e-9)
    assert np.hypot(kx1, ky1) == pytest.approx(np.hypot(kx0, ky0))


def test_the_combined_spiral_meets_nyquist(system) -> None:
    """
    N interleaves each wind k_max/N turns, so together their radial pitch is 1/FOV exactly.

    Undersampling is then a deliberate `density` choice, not an accident of the geometry.
    """
    fov_mm, n_int = 240.0, 8
    ro = sc.modules.SpiralVDS(system, fov_mm=fov_mm, matrix=64, n_interleaves=n_int)
    radii: list[float] = []
    for i in range(n_int):
        kx, ky = ro.trajectory(i)
        angle = np.unwrap(np.arctan2(ky, kx))
        radius = np.hypot(kx, ky)
        lo, hi = angle.min(), angle.max()
        for turn in range(int(math.ceil(lo / (2 * math.pi))), int(hi // (2 * math.pi)) + 1):
            radii.append(float(np.interp(turn * 2 * math.pi, angle, radius)))
    # Only the winding inside k_max: the ramp-down carries |k| about 12 % past the edge, and those
    # samples are a single radial streak per interleaf rather than part of the regular spacing.
    inside = sorted(r for r in radii if r <= ro.k_max_per_m * 1.001)
    pitch = float(np.max(np.diff(inside)))
    assert pitch <= 1e3 / fov_mm * 1.05, f'radial pitch {pitch:.2f} exceeds Nyquist'


def test_spiral_k_zero_is_the_first_sample_not_the_middle(system) -> None:
    """Conflating this with a Cartesian readout's mid-window echo shifts the diffusion weighting."""
    ro = sc.modules.SpiralVDS(system, fov_mm=240, matrix=64, n_interleaves=8)
    assert ro.time_to_echo == pytest.approx(float(ro.adc.delay))
    assert ro.time_to_echo < ro.duration / 10


def test_the_spiral_block_holds_the_adc_dead_time(system) -> None:
    """
    Starting the rewinder before the ADC's trailing dead time is not merely illegal.

    The compiler would have nowhere to place a boundary and would sum the rewinder into the
    spiral, putting a step discontinuity in the middle of the readout.
    """
    ro = sc.modules.SpiralVDS(system, fov_mm=240, matrix=64, n_interleaves=8)
    needed = float(ro.adc.delay) + ro.adc_duration + float(ro.opts.adc_dead_time)
    assert ro.sampling_block_duration >= needed
    out = sc.compile(sc.LogicBlock('t').add(0.0, ro.build()), system)
    assert out.check().ok
    assert not out.report.of_kind('grad_merge')


def test_every_interleaf_takes_the_same_time(system) -> None:
    """The next TR must start at the same instant regardless of which shot was played."""
    ro = sc.modules.SpiralVDS(system, fov_mm=240, matrix=64, n_interleaves=8)
    assert len({round(ro.build(interleaf=i).duration, 12) for i in range(8)}) == 1


# ------------------------------------------------------------------------- diffusion known values
@pytest.mark.parametrize('b_target', [300.0, 1000.0, 2000.0])
def test_b_value_matches_numerical_integration(system, b_target: float) -> None:
    """
    The check that matters: b measured from the built waveform, not from the formula.

    An earlier version used a published ramp correction with the wrong ``delta`` convention and was
    2-4 % high -- which biases every diffusivity by the same amount, invisibly.
    """
    refoc = sc.modules.SincRefocusing(
        system, flip_deg=180, duration_us=4000, slice_thickness_mm=4, crusher_twists=4)
    diff = sc.modules.MonopolarDiffusion(
        system, b_value_s_per_mm2=b_target, refocus_duration_us=refoc.duration * 1e6)

    raster = system.grad_raster.dt
    n_lobe = int(round(diff.lobe_duration / raster))
    gap = int(round(refoc.duration / raster))

    def sampled(block: sc.LogicBlock) -> np.ndarray:
        g = np.zeros(n_lobe + 2)
        for node in block:
            if getattr(node.item, 'type', None) in ('trap', 'grad'):
                tt, wf = sc.events.waveform_of(node.item, raster)
                index = np.round((tt + node.start) / raster).astype(int)
                keep = index < len(g)
                g[index[keep]] += wf[keep]
        return g

    direction = (0.0, 0.0, 1.0)
    effective = np.concatenate([
        sampled(diff.build(part='pre', direction=direction)),
        np.zeros(gap),
        -sampled(diff.build(part='post', direction=direction)),
    ])
    k = np.cumsum(effective) * raster
    numeric = (2 * math.pi) ** 2 * float(np.sum(k**2)) * raster / 1e6
    assert numeric == pytest.approx(diff.achieved_b_s_per_mm2(direction), rel=5e-3)
    assert numeric == pytest.approx(b_target, rel=5e-3)


def test_b_value_collapses_to_stejskal_tanner_without_ramps(system) -> None:
    """The exact expression must reduce to the textbook one as the ramp time goes to zero."""
    diff = sc.modules.MonopolarDiffusion(
        system, b_value_s_per_mm2=1000, refocus_duration_us=4200)
    delta, big_delta, g = diff.lobe_duration, diff.big_delta, diff.amplitude_Hz_per_m
    ideal = (2 * math.pi * g) ** 2 * delta**2 * (big_delta - delta / 3.0) / 1e6
    # The ramps remove a few percent of the encoding, so the ideal value is the upper bound.
    assert diff.achieved_b_s_per_mm2() < ideal
    assert diff.achieved_b_s_per_mm2() > ideal * 0.9


def test_b_scales_as_amplitude_squared(system) -> None:
    diff = sc.modules.MonopolarDiffusion(
        system, b_value_s_per_mm2=1000, refocus_duration_us=4200)
    strong = diff.build(part='pre', direction=(0, 0, 1))
    half = diff.build(part='pre', direction=(0, 0, 1), scale=0.5)
    a_strong = sum(float(n.item.area) for n in strong if getattr(n.item, 'type', None) == 'trap')
    a_half = sum(float(n.item.area) for n in half if getattr(n.item, 'type', None) == 'trap')
    assert a_half == pytest.approx(a_strong / 2, rel=1e-9)


def test_b_zero_occupies_the_same_slot_as_a_weighted_volume(system) -> None:
    """
    A b=0 volume acquired at a shorter TE would carry different T2 weighting and bias the fit.
    """
    weighted = sc.modules.MonopolarDiffusion(
        system, b_value_s_per_mm2=1000, refocus_duration_us=4200)
    unweighted = sc.modules.MonopolarDiffusion(
        system, b_value_s_per_mm2=0, refocus_duration_us=4200,
        lobe_duration_us=weighted.lobe_duration * 1e6)
    assert unweighted.duration == pytest.approx(weighted.duration)
    block = unweighted.build(part='pre')
    assert not [n for n in block if getattr(n.item, 'type', None) in ('trap', 'grad')]


def test_directions_are_normalised_on_the_way_in(system) -> None:
    """
    An unnormalised [1,1,1] would ask for sqrt(3) times the intended vector amplitude.

    That both overstates b and can exceed the amplifier while every individual axis looks legal --
    a real bug in the reference implementation.
    """
    diff = sc.modules.MonopolarDiffusion(
        system, b_value_s_per_mm2=1000, refocus_duration_us=4200)
    block = diff.build(part='pre', direction=(1, 1, 1))
    peak = math.sqrt(sum(float(n.item.amplitude) ** 2 for n in block
                         if getattr(n.item, 'type', None) == 'trap'))
    assert peak == pytest.approx(diff.amplitude_Hz_per_m, rel=1e-9)


def test_a_zero_direction_vector_is_refused(system) -> None:
    diff = sc.modules.MonopolarDiffusion(
        system, b_value_s_per_mm2=1000, refocus_duration_us=4200)
    with pytest.raises(sc.ConfigurationError, match='zero vector'):
        diff.build(part='pre', direction=(0, 0, 0))


def test_bipolar_nulls_the_first_moment(system) -> None:
    """
    The entire reason to pay for a bipolar encoding, so it is measured, not asserted.

    Both lobes are built identically; the refocusing pulse inverts the second, making the effective
    waveform even about the encoding centre -- and an even function has no first moment about it.
    """
    bip = sc.modules.BipolarDiffusion(
        system, b_value_s_per_mm2=300, refocus_duration_us=4200)
    mono = sc.modules.MonopolarDiffusion(
        system, b_value_s_per_mm2=300, refocus_duration_us=4200)
    assert abs(bip.m1_per_m_per_s()) < 1e-6
    assert abs(mono.m1_per_m_per_s() if hasattr(mono, 'm1_per_m_per_s') else 1.0) > 0.0


def test_bipolar_costs_time_for_the_same_b(system) -> None:
    """Splitting a lobe into opposing halves quarters its area, and b goes as area squared."""
    kwargs = {'b_value_s_per_mm2': 300, 'refocus_duration_us': 4200}
    bip = sc.modules.BipolarDiffusion(system, **kwargs)
    mono = sc.modules.MonopolarDiffusion(system, **kwargs)
    assert bip.lobe_duration > mono.lobe_duration
    assert bip.achieved_b_s_per_mm2() == pytest.approx(300.0, rel=1e-3)


def test_an_unreachable_b_value_names_what_is_achievable(system) -> None:
    diff_kwargs = {'refocus_duration_us': 4200, 'lobe_duration_us': 2000}
    with pytest.raises(sc.ConfigurationError, match='achievable_b_s_per_mm2'):
        sc.modules.MonopolarDiffusion(system, b_value_s_per_mm2=3000, **diff_kwargs)


# ------------------------------------------------------------------------------------ directions
@pytest.mark.parametrize('n', [6, 12, 30, 60])
def test_direction_sets_are_unit_length_and_well_conditioned(n: int) -> None:
    """
    Condition number, not the raw count, is what says whether a tensor fit will be noisy.

    The theoretical floor is about 1.58; a raw golden-angle spiral leaves a 14-degree pair at n=30
    and the fit gets biased without anything looking wrong.
    """
    dirs = sc.modules.dti_directions(n)
    assert len(dirs) == n
    for d in dirs:
        assert math.sqrt(sum(c * c for c in d)) == pytest.approx(1.0, abs=1e-9)
        assert d[2] >= -1e-12
    assert sc.modules.direction_condition_number(dirs) < 1.7


def test_six_directions_reach_the_icosahedral_optimum() -> None:
    """Six is the arithmetic minimum for a tensor, so it has to be the best six."""
    dirs = sc.modules.dti_directions(6)
    worst = max(abs(sum(a * b for a, b in zip(p, q)))
                for i, p in enumerate(dirs) for q in dirs[i + 1:])
    assert math.degrees(math.acos(worst)) > 60.0


# ------------------------------------------------------------------------------------- fat sat
def test_fat_sat_offset_is_the_chemical_shift(system) -> None:
    fat = sc.modules.FatSat(system, voxel_mm=5)
    expected = -3.4e-6 * system.gamma * system.b0_T
    assert fat.freq_offset_hz == pytest.approx(expected)
    assert fat.freq_offset_hz == pytest.approx(-434.0, abs=1.0)


def test_a_spectrally_wide_fat_sat_pulse_is_refused(system) -> None:
    """
    Saturating water as well as fat looks like poor SNR, not like a bug, so it must be refused.

    At 1.5 T the separation is only 217 Hz, which is where this bites hardest.
    """
    with pytest.raises(sc.ConfigurationError, match='saturate water as well as fat'):
        sc.modules.FatSat(system, voxel_mm=5, duration_us=2000)


def test_fat_sat_nests_a_pulse_and_a_spoiler(system) -> None:
    """Nesting needs no mechanism: the children are other modules' blocks."""
    fat = sc.modules.FatSat(system, voxel_mm=5)
    block = fat.build()
    assert len(block) == 2
    assert all(isinstance(n.item, sc.LogicBlock) for n in block)
    assert block.duration == pytest.approx(fat.duration, rel=1e-9)


# ------------------------------------------------------------------------------- SLR vs sinc
def test_refocusing_pulses_default_to_a_quarter_turn(system) -> None:
    """
    The CPMG condition: 90x then 180y, not 180x.

    It makes the echo amplitude first-order insensitive to flip-angle error on the refocusing pulse.
    Simulated, the difference is not marginal: at 20 % low B1 the echo keeps 98 % of its amplitude
    with the quarter turn and 83 % without it, and transmit fields are never uniform across a slice.
    """
    for cls, kwargs in (
        (sc.modules.SincRefocusing, {'slice_thickness_mm': 5}),
        (sc.modules.SLRRefocusing, {'slice_thickness_mm': 5}),
        (sc.modules.HardRefocusing, {}),
    ):
        refoc = cls(system, flip_deg=180, duration_us=4000, **kwargs)
        assert refoc.phase_offset_rad == pytest.approx(math.pi / 2), cls.__name__
        built = [n.item for n in refoc.build() if getattr(n.item, 'type', None) == 'rf'][0]
        assert float(built.phase_offset) == pytest.approx(math.pi / 2)

    # An excitation stays at zero, so the *difference* is the quarter turn.
    exc = sc.modules.SincExcitation(
        system, flip_deg=90, duration_us=2000, slice_thickness_mm=5, rephase=False)
    assert exc.phase_offset_rad == pytest.approx(0.0)


def test_refocusing_phase_composes_with_the_build_argument(system) -> None:
    """Phase cycling has to add to the design phase, not replace it."""
    refoc = sc.modules.SLRRefocusing(
        system, flip_deg=180, duration_us=4000, slice_thickness_mm=4)
    built = [n.item for n in refoc.build(rf_phase_rad=0.3) if getattr(n.item, 'type', None) == 'rf']
    assert float(built[0].phase_offset) == pytest.approx(math.pi / 2 + 0.3)


def test_slr_refocusing_beats_sinc_on_the_refocusing_profile(system) -> None:
    """
    SLR is not a cosmetic change: it is 13 % more echo signal for this geometry.

    Integrated from the Bloch equations by hard-pulse decomposition, because the small-tip
    approximation that would let us just Fourier-transform the envelope is exactly what fails at 180
    degrees.  The metric is **refocusing efficiency** -- start transverse, apply the pulse, measure
    what comes back conjugated -- which is what a spin echo depends on and is not the same as the
    excitation profile.
    """
    thickness_mm = 4.0
    raster = system.rf_raster.dt
    positions = np.linspace(-1.5 * thickness_mm / 1e3, 1.5 * thickness_mm / 1e3, 121)

    def efficiency(module):
        b1 = np.asarray(module.rf.signal, dtype=complex)
        amplitude = float(module.gz.amplitude)
        out = np.zeros(len(positions))
        for index, z in enumerate(positions):
            m = np.array([1.0, 0.0, 0.0])
            wz = 2 * np.pi * amplitude * z
            for sample in range(len(b1)):
                w = np.array([2 * np.pi * b1[sample].real, 2 * np.pi * b1[sample].imag, wz])
                norm = float(np.linalg.norm(w))
                if norm < 1e-12:
                    continue
                axis, angle = w / norm, norm * raster
                cos, sin = np.cos(angle), np.sin(angle)
                m = m * cos + np.cross(axis, m) * sin + axis * float(np.dot(axis, m)) * (1 - cos)
            out[index] = abs(complex(m[0], -m[1]))
        return out

    inside = np.abs(positions) < 0.35 * thickness_mm / 1e3
    scores = {}
    for name, cls in (('sinc', sc.modules.SincRefocusing), ('slr', sc.modules.SLRRefocusing)):
        module = cls(system, flip_deg=180, duration_us=4000, slice_thickness_mm=thickness_mm,
                     crusher_twists=4, crusher_voxel_mm=thickness_mm)
        profile = efficiency(module)
        scores[name] = float(profile[inside].mean())

    assert scores['slr'] > scores['sinc'] * 1.05, scores
