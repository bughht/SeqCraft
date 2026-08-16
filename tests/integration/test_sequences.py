"""
End-to-end: whole sequences, checked against physics rather than against themselves.

The sequences are built by the helpers in ``conftest.py``, out of raw pypulseq events.  That is
deliberate and is the point of this tier: it demonstrates that the compile path stands alone, with
no module layer involved at all.

Every assertion is a number an independent calculation gives -- k-space extent from
``matrix / FOV``, the echo at the requested TE, ``k`` at the echo from pypulseq's own trajectory
calculation.  None of them would pass merely because the code is self-consistent.

The DTI and EPI-DWI tiers that used to live here went with the module library they were built on.
They are the acceptance test for whatever module set is written next, so they come back when it
does -- rebuilt against it rather than ported.
"""

from __future__ import annotations

import hashlib
import json
import math
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np
import pypulseq as pp
import pytest
from conftest import FOV_MM, MATRIX, OPTS, build_gre

import seqcraft as sc
from tools.capture_compiler_baseline import stable_summary

_BASELINE = json.loads(
    (Path(__file__).parents[1] / 'baselines' / 'compiler_phase0.json').read_text()
)


# ------------------------------------------------------------------------- every recipe passes
def test_every_recipe_compiles(compiled: pp.Sequence) -> None:
    """
    The whole of what ``check()`` used to assert, in three lines.

    Every legality failure raises now -- limits on the summed waveform, duplicate k-space
    addresses, per-event sample counts, pypulseq's own timing audit, and the four
    against-the-tree invariants -- so a recipe that reaches this line has passed all of them.
    What is left to say is that the thing handed back is a bare pypulseq sequence with no report
    bolted to it.
    """
    assert type(compiled).__module__.startswith('pypulseq')
    assert not hasattr(compiled, 'report')
    assert len(compiled.block_events) > 0


def test_every_recipe_round_trips_through_a_file(compiled: pp.Sequence) -> None:
    """Write, read back, and get the same block count and duration."""
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'seq.seq'
        compiled.write(str(path))
        assert path.exists()
        again = pp.Sequence(system=OPTS)
        again.read(str(path))
        assert len(again.block_events) == len(compiled.block_events)
        assert float(again.duration()[0]) == pytest.approx(compiled.duration()[0], abs=1e-9)


def test_every_recipe_records_its_definitions(compiled: pp.Sequence) -> None:
    """
    A file has to be interpretable a year later without the script that made it.

    Set by the compile rather than by ``write()``, which is what makes the returned sequence
    self-sufficient: no wrapper has to survive until write time to carry them.
    """
    defs = compiled.definitions
    assert defs['Name']
    assert 'TE' in defs
    assert 'TR' in defs
    assert 'FOV' in defs
    assert defs['TotalDuration'] == pytest.approx(compiled.duration()[0])


def test_every_block_duration_lands_on_the_raster(compiled: pp.Sequence) -> None:
    raster = sc.Raster(OPTS.block_duration_raster, 'block')
    off = [
        index for index, duration in compiled.block_durations.items()
        if not raster.holds(float(duration))
    ]
    assert not off, f'blocks off the raster: {off[:5]}'


def test_every_recipe_matches_the_frozen_compiler_baseline(compiled: pp.Sequence) -> None:
    """
    Freeze the compiler fields that are exact across supported CI platforms.

    This is the revision's invariant I1, asserted on every run: what the compiler *emits* --
    block durations, event content, origins, moments -- does not move while what it *returns*
    does.
    """
    name = str(compiled.definitions['Name'])
    expected = _BASELINE['recipes'][name]['stable']
    actual = stable_summary(compiled)
    for field, value in expected.items():
        assert actual[field] == value, f'{name}: compiler baseline changed for {field}'


# ------------------------------------------------------------------------------- GRE physics
def test_gre_kspace_extent_matches_matrix_over_fov(gre_k) -> None:
    """
    ``k_max = matrix / (2 FOV)``, and with an even sample count the outermost sample sits one half
    step inside it -- so the expected extent is ``k_max * (1 - 1/matrix)``.
    """
    k = gre_k['k_adc']
    k_max = MATRIX / (2.0 * FOV_MM / 1e3)
    expected = k_max * (1.0 - 1.0 / MATRIX)
    assert float(np.abs(k[0]).max()) == pytest.approx(expected, rel=0.02)
    assert float(np.abs(k[1]).max()) == pytest.approx(k_max, rel=0.02)


def test_gre_dk_is_one_over_fov(gre_k) -> None:
    """The step between phase-encode lines is what sets the FOV, so it is the thing to check."""
    k = gre_k['k_adc']
    lines = np.unique(np.round(k[1], 6))
    steps = np.diff(np.sort(lines))
    assert float(np.median(steps)) == pytest.approx(1e3 / FOV_MM, rel=0.02)


def test_gre_echo_lands_at_te(gre: pp.Sequence, gre_k) -> None:
    """
    The definition of TE: the centre of the ADC window relative to the excitation centre.

    Asserted on the window centre rather than on a sample, because with an even sample count no
    sample sits exactly at k=0 -- the two central ones straddle it by half a dwell.
    """
    n = MATRIX
    window = gre_k['t_adc'][:n]
    centre = 0.5 * (window[n // 2 - 1] + window[n // 2])
    assert centre - gre_k['t_excitation'][0] == pytest.approx(gre.definitions['TE'], abs=1e-5)


def test_gre_no_sample_sits_exactly_at_k_zero(gre_k) -> None:
    """
    The corollary, stated so nobody 'fixes' the test above.

    With 32 samples the two central ones straddle k=0 by half a dwell each.
    """
    line = gre_k['k_adc'][0][:MATRIX]
    assert float(np.abs(line).min()) > 0.0


def test_gre_rf_spoiling_uses_the_closed_form(gre: pp.Sequence) -> None:
    """
    ``phi_n = inc * n(n+1)/2``, evaluated rather than accumulated.

    An accumulator makes shot 900 depend on shots 1 to 899 having been built first, which is a
    reproducibility hazard rather than a physics one -- but a real one.
    """
    increment = math.radians(gre.definitions['RfSpoilIncrementDeg'])
    phases = [
        float(gre.get_block(i).rf.phase_offset)
        for i in sorted(gre.block_events)
        if getattr(gre.get_block(i), 'rf', None) is not None
    ]
    assert len(phases) > 4
    # The slice offset contributes a constant of its own, so compare against shot 0 -- whose own
    # quadratic term is zero, leaving exactly the closed form.
    for n in range(1, 5):
        expected = increment * n * (n + 1) / 2
        got = phases[n] - phases[0]
        error = abs(got - expected) % (2 * math.pi)
        assert min(error, 2 * math.pi - error) < 1e-6


def test_gre_slices_get_different_rf_frequencies(gre: pp.Sequence) -> None:
    """Two slices means two carrier offsets; one would put both slices on top of each other."""
    offsets = {
        round(float(gre.get_block(i).rf.freq_offset), 3)
        for i in sorted(gre.block_events)
        if getattr(gre.get_block(i), 'rf', None) is not None
    }
    assert len(offsets) == 2


def test_gre_three_winders_share_one_block(gre: pp.Sequence) -> None:
    """
    The overlap the whole design exists for: a rephaser, a blip and a prephaser in one block.

    Three axes, three unrelated gradients, no coordination between them -- and no warning, because
    there is nothing wrong.
    """
    found = any(
        all(getattr(gre.get_block(i), f'g{axis}', None) is not None for axis in 'xyz')
        and getattr(gre.get_block(i), 'adc', None) is None
        for i in sorted(gre.block_events)
    )
    assert found, 'no block carries all three winders'
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        build_gre()
    assert not [w for w in caught if 'merge' in str(w.message)], (
        'three unrelated gradients on three axes are not a merge'
    )


# -------------------------------------------------------------------------- spin-echo physics
def test_se_echo_lands_at_te(se: pp.Sequence, se_k) -> None:
    n = MATRIX
    window = se_k['t_adc'][:n]
    centre = 0.5 * (window[n // 2 - 1] + window[n // 2])
    assert centre - se_k['t_excitation'][0] == pytest.approx(se.definitions['TE'], abs=1e-5)


def test_se_refocusing_pulse_sits_at_half_te(se: pp.Sequence, se_k) -> None:
    """If the 180 is not at TE/2 the echo does not form at TE, whatever the readout says."""
    half = se_k['t_refocusing'][0] - se_k['t_excitation'][0]
    assert half == pytest.approx(se.definitions['TE'] / 2.0, abs=5e-5)


def test_se_kspace_is_not_shifted_by_the_refocusing_pulse(se_k) -> None:
    """
    The sign trap: winders before the 180 must be inverted, or the readout ends at 3 k_max.

    Nothing errors when this is wrong -- the extent check is the only thing that catches it.
    """
    k = se_k['k_adc']
    k_max = MATRIX / (2.0 * FOV_MM / 1e3)
    assert float(np.abs(k[0]).max()) == pytest.approx(k_max * (1 - 1 / MATRIX), rel=0.02)


# ----------------------------------------------------------------------------- reproducibility
def _written(seq, path: Path) -> str:
    """Write `seq` and return the file's sha256."""
    seq.write(str(path))
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_building_the_same_recipe_twice_gives_the_same_file() -> None:
    """
    Byte identity, which is a legitimate assertion because pulseq writes no timestamp.

    This is the refactor guard: 'my output did not change unintentionally'.  There is no sidecar
    to suppress any more -- the JSON that recorded versions and git state went with the result
    wrapper, and nothing records them for the moment.  See docs/serialization.md.
    """
    with tempfile.TemporaryDirectory() as tmp:
        first = _written(build_gre(), Path(tmp) / 'a.seq')
        second = _written(build_gre(), Path(tmp) / 'b.seq')
        assert first == second


def test_changing_te_changes_the_file() -> None:
    """The other half of the guard: identical output for different inputs would be worse."""
    with tempfile.TemporaryDirectory() as tmp:
        a = _written(build_gre(te_s=8e-3), Path(tmp) / 'a.seq')
        b = _written(build_gre(te_s=9e-3), Path(tmp) / 'b.seq')
        assert a != b


# ------------------------------------------------------------------------------- the scanner
def test_a_derated_opts_does_not_lose_the_site_constants() -> None:
    """
    The multi-regime ``System`` guaranteed this with a cross-regime consistency check.

    With one ``Opts`` per design there is nothing to check, because there is nothing that can
    disagree -- provided a derated scanner is *derived* rather than rebuilt.  A sequence compiled
    against one that lost its dead times passes every check here and is refused at the console.
    """
    epi = sc.opts.derate(OPTS, grad=0.85, slew=0.6)

    assert epi.rf_dead_time == OPTS.rf_dead_time
    assert epi.adc_dead_time == OPTS.adc_dead_time
    assert epi.grad_raster_time == OPTS.grad_raster_time
    assert epi.gamma == OPTS.gamma


# --------------------------------------------------------------------------------- performance
@pytest.mark.slow
def test_a_realistic_sequence_compiles_quickly() -> None:
    """
    Budget check, not a benchmark.

    A full acquisition is thousands of TRs and hundreds of thousands of events; the sweep-based
    interval assignment is O(n log n), and this fails loudly if that ever regresses to quadratic.
    """
    start = time.perf_counter()
    seq = build_gre(matrix=128, n_slices=8)
    elapsed = time.perf_counter() - start

    n_blocks = len(seq.block_events)
    assert n_blocks > 4000
    assert elapsed < 60.0, f'{n_blocks} blocks took {elapsed:.1f} s'
