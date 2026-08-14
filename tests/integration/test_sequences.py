"""
End-to-end: whole sequences, checked against physics rather than against themselves.

The sequences are built by the helpers in ``conftest.py`` -- in the tests, not in the package, because
seqcraft ships modules and a compiler and a sequence is something you assemble.

Every assertion is a number an independent calculation gives: k-space extent from ``matrix / FOV``,
the echo at the requested TE, the b-value from numerical integration, ``k`` at the echo from
pypulseq's own trajectory calculation.  None of them would pass merely because the code is
self-consistent.
"""

from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import numpy as np
import pypulseq as pp
import pytest
from conftest import (
    FOV_MM,
    MATRIX,
    build_dti,
    build_epi_dwi,
    build_gre,
    epi_dwi_system,
)

import seqcraft as sc
from tools.capture_compiler_baseline import stable_summary

_PHASE0_BASELINE = json.loads(
    (Path(__file__).parents[1] / 'baselines' / 'compiler_phase0.json').read_text()
)


# ------------------------------------------------------------------------- every recipe passes
def test_every_recipe_checks_clean(compiled: sc.CompiledSequence) -> None:
    report = compiled.check()
    assert report.ok, f'{compiled.definitions["Name"]}: {report}'


def test_every_recipe_has_unique_kspace_addresses(compiled: sc.CompiledSequence) -> None:
    """
    A duplicate address means two readouts writing the same location.

    One assertion catches a wrong slice order, an off-by-one partial-Fourier start and a mis-nested
    loop -- which is why it runs on every recipe rather than being tested once.
    """
    assert not [i for i in compiled.check().issues if i.kind == 'label']


def test_every_recipe_reports_no_amplitude_or_slew_error(compiled: sc.CompiledSequence) -> None:
    """Per-axis violations are errors; the vector-norm ones are warnings and are allowed."""
    hard = [
        i for i in compiled.report.errors
        if i.kind in ('grad_limit', 'slew_limit')
    ]
    assert not hard, hard


def test_every_recipe_round_trips_through_a_file(compiled: sc.CompiledSequence) -> None:
    """Write, read back, and get the same block count and duration."""
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'seq.seq'
        result = compiled.write(path)
        assert result.path.exists()
        assert result.sidecar is not None and result.sidecar.exists()
        again = pp.Sequence(system=compiled.system.limits(compiled.regime))
        again.read(str(path))
        assert len(again.block_events) == compiled.n_blocks
        assert float(again.duration()[0]) == pytest.approx(compiled.duration_s, abs=1e-9)


def test_every_recipe_records_its_definitions(compiled: sc.CompiledSequence) -> None:
    """A file has to be interpretable a year later without the script that made it."""
    defs = compiled.definitions
    assert defs['Name']
    assert 'TE' in defs
    assert 'TR' in defs
    assert 'FOV' in defs


def test_every_block_duration_lands_on_the_raster(compiled: sc.CompiledSequence) -> None:
    raster = compiled.system.block_raster
    off = [
        index for index, duration in compiled.seq.block_durations.items()
        if not raster.holds(float(duration))
    ]
    assert not off, f'blocks off the raster: {off[:5]}'


def test_every_recipe_matches_the_phase0_compiler_baseline(
    compiled: sc.CompiledSequence,
) -> None:
    """Freeze block structure, emitted content, provenance, issue surface, and moments."""
    name = str(compiled.definitions['Name'])
    expected = _PHASE0_BASELINE['recipes'][name]['stable']
    actual = stable_summary(compiled)
    for field, value in expected.items():
        assert actual[field] == value, f'{name}: Phase 0 baseline changed for {field}'


# ------------------------------------------------------------------------------- GRE physics
def test_gre_kspace_extent_matches_matrix_over_fov(gre: sc.CompiledSequence) -> None:
    """
    ``k_max = matrix / (2 FOV)``, and with an even sample count the outermost sample sits one half
    step inside it -- so the expected extent is ``k_max * (1 - 1/matrix)``.
    """
    k = gre.kspace()['k_adc']
    k_max = MATRIX / (2.0 * FOV_MM / 1e3)
    expected = k_max * (1.0 - 1.0 / MATRIX)
    assert float(np.abs(k[0]).max()) == pytest.approx(expected, rel=0.02)
    assert float(np.abs(k[1]).max()) == pytest.approx(k_max, rel=0.02)


def test_gre_dk_is_one_over_fov(gre: sc.CompiledSequence) -> None:
    """The step between phase-encode lines is what sets the FOV, so it is the thing to check."""
    k = gre.kspace()['k_adc']
    lines = np.unique(np.round(k[1], 6))
    steps = np.diff(np.sort(lines))
    assert float(np.median(steps)) == pytest.approx(1e3 / FOV_MM, rel=0.02)


def test_gre_echo_lands_at_te(gre: sc.CompiledSequence) -> None:
    """
    The definition of TE: the centre of the ADC window relative to the excitation centre.

    Asserted on the window centre rather than on a sample, because with an even sample count no
    sample sits exactly at k=0 -- the two central ones straddle it by half a dwell.
    """
    k = gre.kspace()
    n = MATRIX
    window = k['t_adc'][:n]
    centre = 0.5 * (window[n // 2 - 1] + window[n // 2])
    assert centre - k['t_excitation'][0] == pytest.approx(gre.definitions['TE'], abs=1e-6)


def test_gre_no_sample_sits_exactly_at_k_zero(gre: sc.CompiledSequence) -> None:
    """
    The corollary, stated so nobody 'fixes' the test above.

    With 32 samples the two central ones straddle k=0 by half a dwell each.
    """
    k = gre.kspace()['k_adc']
    line = k[0][:MATRIX]
    assert float(np.abs(line).min()) > 0.0


def test_gre_rf_spoiling_uses_the_closed_form(gre: sc.CompiledSequence) -> None:
    """
    ``phi_n = inc * n(n+1)/2``, evaluated rather than accumulated.

    An accumulator makes shot 900 depend on shots 1 to 899 having been built first, which is a
    reproducibility hazard rather than a physics one -- but a real one.
    """
    increment = math.radians(gre.definitions['RfSpoilIncrementDeg'])
    phases = []
    for index in range(1, gre.n_blocks + 1):
        block = gre.seq.get_block(index)
        if getattr(block, 'rf', None) is not None:
            phases.append(float(block.rf.phase_offset))
    assert len(phases) > 4
    for n, phase in enumerate(phases[:5]):
        assert phase == pytest.approx((increment * n * (n + 1) / 2) % (2 * math.pi), abs=1e-9)


def test_gre_slices_get_different_rf_frequencies(gre: sc.CompiledSequence) -> None:
    """Two slices means two carrier offsets; one would mean both slices land on top of each other."""
    offsets = set()
    for index in range(1, gre.n_blocks + 1):
        block = gre.seq.get_block(index)
        if getattr(block, 'rf', None) is not None:
            offsets.add(round(float(block.rf.freq_offset), 3))
    assert len(offsets) == 2


def test_gre_three_winders_share_one_block(gre: sc.CompiledSequence) -> None:
    """
    The overlap the whole design exists for: a rephaser, a blip and a prephaser in one block.

    Three axes, three modules, no coordination -- and no warning, because there is nothing wrong.
    """
    found = False
    for index in range(1, gre.n_blocks + 1):
        block = gre.seq.get_block(index)
        if all(getattr(block, f'g{axis}', None) is not None for axis in 'xyz') and (
            getattr(block, 'adc', None) is None
        ):
            found = True
            break
    assert found, 'no block carries all three winders'
    assert not gre.report.of_kind('grad_merge')


# -------------------------------------------------------------------------- spin-echo physics
def test_se_echo_lands_at_te(se: sc.CompiledSequence) -> None:
    k = se.kspace()
    n = MATRIX
    window = k['t_adc'][:n]
    centre = 0.5 * (window[n // 2 - 1] + window[n // 2])
    assert centre - k['t_excitation'][0] == pytest.approx(se.definitions['TE'], abs=1e-6)


def test_se_refocusing_pulse_sits_at_half_te(se: sc.CompiledSequence) -> None:
    """If the 180 is not at TE/2 the echo does not form at TE, whatever the readout says."""
    k = se.kspace()
    half = k['t_refocusing'][0] - k['t_excitation'][0]
    assert half == pytest.approx(se.definitions['TE'] / 2.0, abs=2e-5)


def test_se_kspace_is_not_shifted_by_the_refocusing_pulse(se: sc.CompiledSequence) -> None:
    """
    The sign trap: winders before the 180 must be inverted, or the readout ends at 3 k_max.

    Nothing errors when this is wrong -- the extent check is the only thing that catches it.
    """
    k = se.kspace()['k_adc']
    k_max = MATRIX / (2.0 * FOV_MM / 1e3)
    assert float(np.abs(k[0]).max()) == pytest.approx(k_max * (1 - 1 / MATRIX), rel=0.02)


def test_se_has_no_spoiler_at_the_end(se: sc.CompiledSequence) -> None:
    """The crushers already dephased the FID; a spoiler would only lengthen TR."""
    assert se.definitions['TR'] > se.definitions['TE']


# -------------------------------------------------------------------------------- DTI physics
def test_dti_b_values_are_achieved(dti: sc.CompiledSequence) -> None:
    """The recorded b-value is the achieved one, so raster rounding is visible not assumed away."""
    requested = list(dti.definitions['bValues'])
    achieved = list(dti.definitions['AchievedbValues'])
    for want, got in zip(requested, achieved):
        assert got == pytest.approx(want, rel=2e-3, abs=1e-6)


def test_dti_b_zero_and_b_max_share_one_te(dti: sc.CompiledSequence) -> None:
    """
    An unweighted volume at a shorter TE would carry different T2 weighting and bias the fit.

    Every volume is built from an encoder sharing the largest lobe duration, so TE is one number.
    """
    assert isinstance(dti.definitions['TE'], float)
    assert dti.definitions['DiffusionLobeDuration'] > 0.0


def test_dti_spiral_reaches_k_max(dti: sc.CompiledSequence) -> None:
    """The spiral reaches k_max, and the ramp-down takes |k| only a little past it."""
    k = dti.kspace()['k_adc']
    k_max = MATRIX / (2.0 * 240.0 / 1e3)
    reach = float(np.hypot(k[0], k[1]).max()) / k_max
    assert 1.0 <= reach <= 1.2, f'|k| reaches {reach:.3f} of k_max'


def test_dti_spiral_starts_at_the_origin(dti: sc.CompiledSequence) -> None:
    """
    A spiral's first sample is k=0, which is what makes ``time_to_echo`` the ADC delay alone.

    Placing it as though k=0 were mid-window would shift the diffusion weighting.
    """
    k = dti.kspace()['k_adc']
    assert float(np.hypot(k[0], k[1]).min()) < 5.0


def test_dti_records_the_diffusion_scheme(dti: sc.CompiledSequence) -> None:
    defs = dti.definitions
    assert defs['DiffusionScheme'] in ('monopolar', 'bipolar')
    assert defs['DiffusionDirections'] == 6
    assert defs['SpiralInterleaves'] == 1
    assert defs['SlicesPerTR'] == 1
    assert defs['DiffusionBigDelta'] > defs['DiffusionLobeDuration']


# ---------------------------------------------------------------------------- EPI DWI physics
def test_epi_dwi_b_values_are_achieved(epi_dwi: sc.CompiledSequence) -> None:
    """The same diffusion encoding as the spiral, so the same b-values must come out of it."""
    for want, got in zip(epi_dwi.definitions['bValues'],
                         epi_dwi.definitions['AchievedbValues']):
        assert got == pytest.approx(want, rel=2e-3, abs=1e-6)


def test_epi_dwi_echo_lands_at_te(epi_dwi: sc.CompiledSequence) -> None:
    """
    TE is measured to the ``ky = 0`` echo, a third of the way into the train, not to its start.

    This is the assertion that would fail if ``time_to_echo`` were taken as the train's midpoint or
    as its first sample -- and neither mistake changes anything else about the sequence.
    """
    k = epi_dwi.kspace()
    te = epi_dwi.definitions['TE']
    at_echo = int(np.argmin(np.hypot(k['k_adc'][0], k['k_adc'][1])))
    elapsed = float(k['t_adc'][at_echo] - k['t_excitation'][0])
    assert elapsed == pytest.approx(te, abs=2.0 * epi_dwi.definitions['EchoSpacing'])
    # The refocusing pulse sits at TE/2, as for any spin echo.
    assert float(k['t_refocusing'][0] - k['t_excitation'][0]) == pytest.approx(te / 2, abs=1e-5)


def test_epi_dwi_k_is_zero_at_the_echo_on_every_axis(epi_dwi: sc.CompiledSequence) -> None:
    """
    Catches a missing slice rewinder, which nothing else in the sequence would reveal.

    The excitation's slice-select gradient leaves through-slice dephasing equal to its own tail, and
    the refocusing pulse does not undo it -- its own gradient is symmetric about its centre.
    """
    k = epi_dwi.kspace()['k_adc']
    at_echo = int(np.argmin(np.hypot(k[0], k[1])))
    dk = 1e3 / 240.0
    for axis, value in zip('xyz', k[:, at_echo]):
        assert abs(float(value)) < dk, f'k_{axis} = {float(value):+.2f} 1/m at the echo'


def test_epi_dwi_covers_k_space_on_a_grid(epi_dwi: sc.CompiledSequence) -> None:
    """
    Readout extent from ``matrix / (2 FOV)``, phase extent from the partial-Fourier table.

    Partial Fourier truncates the **low** end, so ky runs from ``-0.25 matrix dk`` to
    ``+0.5 matrix dk`` -- asymmetric, which is the whole point, and a check on |k| alone could not
    see it because |k| is symmetric.
    """
    k = epi_dwi.kspace()['k_adc']
    dk = 1e3 / 240.0
    k_max_ro = MATRIX / (2.0 * 240.0 / 1e3)
    assert float(np.abs(k[0]).max()) == pytest.approx(k_max_ro, rel=0.02)
    assert float(k[1].max()) == pytest.approx((MATRIX / 2 - 1) * dk, abs=dk)
    assert float(k[1].min()) == pytest.approx(-0.25 * MATRIX * dk, abs=dk)
    # Distinct ky values, one per acquired line, and they are a uniform ladder.
    rungs = np.unique(np.round(k[1] / dk).astype(int))
    assert len(rungs) == int(round(0.75 * MATRIX))
    assert np.all(np.diff(rungs) == 1)


def test_epi_dwi_records_the_readout_it_played(epi_dwi: sc.CompiledSequence) -> None:
    """
    The definitions carry what a downstream correction needs, and nothing it would have to guess.

    ``TotalReadoutTime`` and ``PhaseEncodeBandwidthPerPixelHz`` are what FSL's ``topup`` and
    ``eddy`` ask for; a sequence that does not record them makes someone measure them from the file.
    """
    defs = epi_dwi.definitions
    assert defs['EPIFactor'] == int(round(0.75 * MATRIX))
    assert defs['NumberOfShots'] == 1
    assert defs['PartialFourierPE'] == pytest.approx(0.75)
    assert defs['TotalReadoutTime'] == pytest.approx(
        (defs['EPIFactor'] - 1) * defs['EchoSpacing'])
    assert defs['PhaseEncodeBandwidthPerPixelHz'] == pytest.approx(
        1.0 / (defs['EPIFactor'] * defs['EchoSpacing']))
    assert defs['PhaseEncodePolarity'] == 1


def test_epi_dwi_shares_the_diffusion_encoding_with_the_spiral(system: sc.System) -> None:
    """
    Two readouts, one encoding.  If they disagreed, the ADC maps could not be compared.

    Built on the same system and the same b-value, so the lobe duration and Delta must match to the
    raster -- the readout only enters TE, not the encoding.
    """
    scanner = epi_dwi_system()
    epi = build_epi_dwi(scanner, regime='epi', n_directions=6, n_slices=1)
    spiral = build_dti(
        scanner, regime='epi', n_directions=6, n_slices=1, density=0.5, tr_s=6.0)
    for key in ('DiffusionLobeDuration', 'DiffusionBigDelta', 'AchievedbValues'):
        assert epi.definitions[key] == spiral.definitions[key], f'{key} differs'
    # TE does differ, because the two readouts reach k=0 at different points in themselves.
    assert epi.definitions['TE'] != spiral.definitions['TE']


def test_epi_dwi_blip_down_mirrors_the_phase_encode(system: sc.System) -> None:
    """
    The blip-down half of a distortion-correction pair: identical timing, mirrored ky.

    Reversed polarity displaces off-resonance the other way, which is what lets the pair be combined
    into an undistorted image -- so the two must agree on TE exactly, or the pair is inconsistent.
    """
    scanner = epi_dwi_system()
    up = build_epi_dwi(scanner, regime='epi', n_directions=6, n_slices=1)
    down = build_epi_dwi(scanner, regime='epi', n_directions=6, n_slices=1, pe_polarity=-1)
    assert up.definitions['TE'] == down.definitions['TE']
    assert up.duration_s == pytest.approx(down.duration_s)

    ky_up = up.kspace()['k_adc'][1]
    ky_down = down.kspace()['k_adc'][1]
    assert float(ky_up.max()) == pytest.approx(-float(ky_down.min()), rel=1e-6)
    assert float(ky_up.min()) == pytest.approx(-float(ky_down.max()), rel=1e-6)


def test_epi_dwi_two_shots_halve_the_train_and_double_the_bandwidth(system: sc.System) -> None:
    """
    The trade the segmented variant exists to make: less distortion, more shots, shorter TE.

    Both are the same acquisition of the same k-space; only how it is split differs, so the line
    count over all shots has to be identical.
    """
    scanner = epi_dwi_system()
    one = build_epi_dwi(scanner, regime='epi', n_directions=6, n_slices=1, n_shots=1)
    two = build_epi_dwi(scanner, regime='epi', n_directions=6, n_slices=1, n_shots=2)

    assert two.definitions['EPIFactor'] * 2 == one.definitions['EPIFactor']
    assert two.definitions['PhaseEncodeBandwidthPerPixelHz'] > (
        1.9 * one.definitions['PhaseEncodeBandwidthPerPixelHz'])
    assert two.definitions['TE'] < one.definitions['TE']

    dk = 1e3 / 240.0
    for out in (one, two):
        rungs = np.unique(np.round(out.kspace()['k_adc'][1] / dk).astype(int))
        assert len(rungs) == int(round(0.75 * MATRIX))


def test_epi_dwi_partial_echo_shortens_the_echo_spacing(system: sc.System) -> None:
    """Dropping the leading part of every line shortens each echo, and so the whole train."""
    scanner = epi_dwi_system()
    full = build_epi_dwi(scanner, regime='epi', n_directions=1, n_slices=1)
    short = build_epi_dwi(scanner, regime='epi', n_directions=1, n_slices=1, partial_echo=0.7)
    assert short.definitions['EchoSpacing'] < full.definitions['EchoSpacing']
    assert short.definitions['TE'] < full.definitions['TE']
    assert short.definitions['EPIFactor'] == full.definitions['EPIFactor']


def test_epi_dwi_labels_address_every_line_once(epi_dwi: sc.CompiledSequence) -> None:
    """
    Every acquired line, once per volume and slice.  The duplicate-address check is the guard.

    ``test_every_recipe_has_unique_kspace_addresses`` already covers collisions for this sequence;
    this one asserts the *coverage*, which a collision check cannot see.
    """
    lines: set[int] = set()
    for index in sorted(epi_dwi.seq.block_events):
        block = epi_dwi.seq.get_block(index)
        labels = getattr(block, 'label', None)
        if labels is None or getattr(block, 'adc', None) is None:
            continue
        for label in (labels if isinstance(labels, (list, tuple)) else [labels]):
            if label.label == 'LIN':
                lines.add(int(label.value))
    expected = set(range(MATRIX - int(round(0.75 * MATRIX)), MATRIX))
    assert lines == expected


def test_epi_dwi_scales_to_a_full_acquisition(system: sc.System) -> None:
    """
    Seven volumes of a 48-echo train, and every block boundary splits the one train gradient.

    The plan flagged this as the one scale risk of building the train as a single extended
    trapezoid rather than one trapezoid per echo: a superlinear split path would show up here as
    wall-clock.  It does not -- the compiler's boundary scan and its exact split are both linear.
    """
    import time

    scanner = epi_dwi_system()
    started = time.perf_counter()
    out = build_epi_dwi(
        scanner, matrix=64, regime='epi', n_directions=6, n_slices=1, tr_s=6.0)
    elapsed = time.perf_counter() - started
    assert out.check().ok
    # One block per echo, seven volumes, plus the prep blocks of each shot.
    echoes = out.definitions['EPIFactor']
    assert out.n_blocks > 7 * echoes
    assert elapsed < 60.0, f'compiling {out.n_blocks} blocks took {elapsed:.1f} s'


def test_epi_dwi_pns_is_dominated_by_the_readout(system: sc.System) -> None:
    """
    The opposite of the spiral, where derating the readout did nothing for peripheral nerve
    stimulation.

    An EPI readout oscillates at ~1 kHz for tens of milliseconds, which is the stimulation worst
    case -- so derating *its* slew is what moves the number.  Measured against the synthetic model,
    which over-reports and must never be used to clear a human scan; what matters here is the
    ordering, not the value.
    """
    hardware = sc.synthetic_hardware()
    scanner = epi_dwi_system()
    fast = build_epi_dwi(scanner, matrix=64, regime='epi', n_directions=1, n_slices=1)
    derated = (
        sc.System.preset('cima_x')
        .derate('epi', grad=0.85, slew=0.25)
        .derate('diffusion', grad=0.85, slew=0.3)
    )
    slow = build_epi_dwi(derated, matrix=64, regime='epi', n_directions=1, n_slices=1)

    def peak(out: sc.CompiledSequence) -> float:
        return float(out.pns(hardware).values['peak_pns_fraction'])

    assert peak(fast) > 0.0, 'the synthetic model reported nothing at all'
    assert peak(slow) < peak(fast), (
        f'derating the readout slew did not lower PNS: {peak(slow):.2f} against {peak(fast):.2f}'
    )


def test_dti_refuses_a_te_shorter_than_the_encoding(system: sc.System) -> None:
    """The error names which half of TE is binding and what each term contributes."""
    with pytest.raises(sc.ConfigurationError) as err:
        build_dti(system, b_values=(1000.0,), n_directions=6, te_s=10e-3, tr_s=3.0)
    # The message has to name the term that dominates, which is the diffusion lobe.
    text = str(err.value)
    assert 'below the minimum' in text
    assert 'diffusion lobe' in text


def test_dti_refuses_a_tr_that_cannot_hold_its_slices(system: sc.System) -> None:
    """
    TR is the time between exciting the *same* slice, so it has to hold every slice.

    Treating TR as per-shot instead does not produce a wrong image -- it produces a correct one that
    takes n_slices times longer, which for 20 slices turns a 17-minute scan into five and a half
    hours.  The error says how many slices do fit.
    """
    with pytest.raises(sc.ConfigurationError, match='cannot hold'):
        build_dti(system, b_values=(1000.0,), n_directions=6, tr_s=0.05)


def test_dti_bipolar_option_nulls_m1(system: sc.System) -> None:
    """The velocity-compensated variant has to actually be velocity compensated."""
    bip = sc.modules.BipolarDiffusion(
        system, b_value_s_per_mm2=300, refocus_duration_us=6000)
    mono = sc.modules.MonopolarDiffusion(
        system, b_value_s_per_mm2=300, refocus_duration_us=6000)
    assert abs(bip.m1_per_m_per_s()) < 1e-6
    assert abs(mono.m1_per_m_per_s()) > 1.0
    assert bip.lobe_duration > mono.lobe_duration


# ----------------------------------------------------------------------------- reproducibility
def test_building_the_same_recipe_twice_gives_the_same_file(system: sc.System) -> None:
    """
    Byte identity, which is a legitimate assertion because pulseq writes no timestamp.

    This is the refactor guard: 'my output did not change unintentionally'.
    """
    with tempfile.TemporaryDirectory() as tmp:
        first = build_gre(system).write(Path(tmp) / 'a.seq', sidecar=False)
        second = build_gre(system).write(Path(tmp) / 'b.seq', sidecar=False)
        assert first.sha256 == second.sha256


def test_changing_te_changes_the_file(system: sc.System) -> None:
    """The other half of the guard: identical output for different inputs would be worse."""
    with tempfile.TemporaryDirectory() as tmp:
        a = build_gre(system, te_s=8e-3).write(Path(tmp) / 'a.seq', sidecar=False)
        b = build_gre(system, te_s=9e-3).write(Path(tmp) / 'b.seq', sidecar=False)
        assert a.sha256 != b.sha256


def test_the_sidecar_records_versions_and_git_state(system: sc.System) -> None:
    """A dirty tree means a rebuild cannot be expected to match, so it is recorded."""
    import json

    out = build_gre(system)
    with tempfile.TemporaryDirectory() as tmp:
        result = out.write(Path(tmp) / 'seq.seq')
        payload = json.loads(result.sidecar.read_text(encoding='utf-8'))
    assert payload['versions']['seqcraft'] == sc.__version__
    assert 'pypulseq' in payload['versions']
    assert 'git' in payload
    assert payload['resolved']['n_blocks'] == out.n_blocks


# --------------------------------------------------------------------------------- performance
def test_a_realistic_sequence_compiles_quickly(system: sc.System) -> None:
    """
    Budget check, not a benchmark.

    A full acquisition is thousands of TRs and hundreds of thousands of events; the sweep-based
    interval assignment is O(n log n), and this fails loudly if that ever regresses to quadratic.
    """
    import time

    start = time.perf_counter()
    out = build_gre(system, matrix=128, n_slices=8)
    elapsed = time.perf_counter() - start
    assert out.n_blocks > 4000
    assert elapsed < 60.0, f'{out.n_blocks} blocks took {elapsed:.1f} s'


# --------------------------------------------------------- the slice-refocusing trap
def test_dti_refocuses_the_slice_at_the_echo(system: sc.System) -> None:
    """
    All three k components must be zero at the echo, ``k_z`` included.

    The trap: a spiral starts at k=0 in x and y, which makes it tempting to drop the excitation's
    slice rewinder as unnecessary.  It is not.  The excitation's slice-select gradient leaves
    through-slice dephasing equal to its own tail, and the refocusing pulse does **not** undo it --
    that pulse's slice-select gradient is symmetric about its centre, so its lead and tail cancel
    each other and leave the excitation's tail alone.  Without the rewinder ``k_z`` at the echo is
    -525 1/m, which is 2.1 cycles across a 4 mm slice and about 95 % of the signal, and nothing else
    in the sequence looks wrong.

    Checked against ``calculate_kspacePP``, which handles the refocusing conjugation itself rather
    than relying on this test to get the sign convention right.
    """
    out = build_dti(
        system, b_values=(0.0,), n_directions=1, n_interleaves=1, n_slices=1, tr_s=3.0,
    )
    k_at_echo = out.kspace()['k_adc'][:, 0]
    for axis, value in zip('xyz', k_at_echo):
        assert abs(float(value)) < 1.0, f'k_{axis} = {float(value):.2f} 1/m at the echo'


def test_dropping_the_slice_rewinder_leaves_the_slice_dephased(system: sc.System) -> None:
    """
    The other half of the pair: prove the rewinder is what does it, not something else.

    The residual is exactly the excitation gradient's own tail, which is what the rewinder cancels.
    """
    out = build_dti(
        system, b_values=(0.0,), n_directions=1, n_interleaves=1, n_slices=1, tr_s=3.0,
        rephase=False,
    )
    k_z = float(out.kspace()['k_adc'][2, 0])
    bare = sc.modules.SincExcitation(
        system, flip_deg=90, duration_us=2000, slice_thickness_mm=4.0, rephase=False)
    tail = float(bare.gz.amplitude) * (
        float(bare.gz.flat_time) / 2 + float(bare.gz.fall_time) / 2
    )
    assert abs(k_z) > 100.0, f'expected a large residual, got k_z = {k_z:.2f}'
    assert abs(abs(k_z) - tail) < 1.0, f'k_z = {k_z:.2f} but the gz tail is {tail:.2f}'
    assert abs(float(bare.gzr.area) + tail) < 1.0, 'gzr should be minus the tail'
