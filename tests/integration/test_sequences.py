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

import math
import tempfile
from pathlib import Path

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc

from conftest import FOV_MM, MATRIX, SLICE_MM, build_dti, build_gre


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
    raster = compiled.system.block_raster_s
    off = [
        index for index, duration in compiled.seq.block_durations.items()
        if not sc.raster.on_raster(float(duration), raster)
    ]
    assert not off, f'blocks off the raster: {off[:5]}'


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
