"""
The authoring contract, checked on the emitted file rather than on the authoring script.

`sc.compile` already refuses an illegal sequence -- limits on the summed waveform, rasters,
duplicate k-space addresses, pypulseq's timing audit -- and `tests/compiler/` covers that.  None of
it answers the question this directory exists for:

    is the file the acquisition the author said they were building?

So every expectation below comes from the *specification* -- the matrix and field of view the
protocol asked for -- and from physics: ``k_y = (line - n/2) / FOV_y``, TE from the RF centre to
the echo, TR between successive excitations.  All of it is measured on the sequence **after it has
been written and read back**.

The distinction is not cosmetic.  An earlier draft of this file compared the emitted table against
the view list the authoring code had just built, which agreed with itself: deleting a line from the
protocol left all of these tests passing.  A check that consults the thing under test for its own
expected answer is not a check.

These are deterministic artifact checks.  Whether an agent *chose* a sensible acquisition is a
different question, reviewed by a human against the calibration prompts, and is deliberately not
asserted here.
"""

from __future__ import annotations

import numpy as np
import pytest
import reference_acquisition as ref

import seqcraft as sc

#: The specification -- what was asked for -- as opposed to anything the build reports about
#: itself.  A fully sampled Cartesian acquisition of this matrix collects every line of it exactly
#: once, and that expectation needs no help from the authoring code's own view list.
NX, NY = ref.MATRIX
FOV_Y_M = ref.FOV_MM[1] * 1e-3
EXPECTED_KY = sorted((line - NY // 2) / FOV_Y_M for line in range(NY))


def _readout_count(seq) -> int:
    return sum(1 for i in range(1, len(seq.block_events) + 1)
               if getattr(seq.get_block(i), 'adc', None) is not None)


# ------------------------------------------------------------------ what the file actually walks
def test_the_emitted_file_covers_the_whole_phase_encode_grid(emitted, trajectory) -> None:
    """
    The acquisition-level property the compiler cannot check: the *whole table* was collected.

    A train that silently drops, repeats or truncates a view fails here even though each individual
    repetition is still legal and the file still compiles.
    """
    k_y = np.asarray(trajectory['k_adc'])[1]
    readouts = _readout_count(emitted)
    acquired = sorted(float(np.mean(k_y[i * NX:(i + 1) * NX])) for i in range(readouts))

    assert readouts == NY, f'{readouts} readouts for a {NY}-line phase-encode grid'
    assert acquired == pytest.approx(EXPECTED_KY, abs=1e-3), 'the grid is not covered exactly once'


def test_the_dummy_repetitions_excite_without_acquiring(emitted, trajectory) -> None:
    """
    A preparation policy that quietly acquires is a different acquisition.

    One excitation per repetition including the dummies; one readout per acquired view only.
    """
    n_excitations = np.asarray(trajectory['t_excitation']).size

    assert _readout_count(emitted) == NY
    assert n_excitations == ref.DUMMIES + NY


# --------------------------------------------------------------------------- timing, as emitted
def test_te_measured_from_the_emitted_echo(built, trajectory) -> None:
    """
    TE is the RF effective centre to the echo, measured on the file -- not the module's `te_s`.

    The echo is located as the ADC sample whose `k_x` is smallest, so the comparison can be no
    sharper than half a dwell; that is a property of sampling, not slack in the check.
    """
    kernel = built['kernel']
    k_x = np.asarray(trajectory['k_adc'])[0]
    t_adc = np.asarray(trajectory['t_adc'])
    n_acquired = t_adc.size // NX
    t_exc = np.asarray(trajectory['t_excitation'])[-n_acquired:]

    assert kernel.te_s >= kernel.min_te_s, 'the declared TE is below this kernel own minimum'

    dwell = float(np.median(np.diff(t_adc[:NX])))
    for i, t0 in enumerate(t_exc):
        window = slice(i * NX, (i + 1) * NX)
        echo_at = t_adc[window][int(np.argmin(np.abs(k_x[window])))]
        assert echo_at - t0 == pytest.approx(kernel.te_s, abs=dwell / 2)


def test_tr_is_the_interval_between_emitted_excitations(built, trajectory) -> None:
    """Constant, and the one the acquisition declared -- measured RF centre to RF centre."""
    kernel = built['kernel']
    intervals = np.diff(np.asarray(trajectory['t_excitation']))

    assert kernel.tr_s >= kernel.min_tr_s, 'the declared TR is below this kernel own minimum'
    assert intervals == pytest.approx(kernel.tr_s, abs=1e-9), 'TR is not constant across the train'


def test_the_echo_sits_where_the_readout_says_it_does(trajectory) -> None:
    """
    `k_x` must sweep through the origin inside every readout, at the sample the protocol puts it on.

    A prephaser wrong by a lobe still produces a legal file; it produces the wrong image.
    """
    k_x = np.asarray(trajectory['k_adc'])[0]

    for i in range(NY):
        window = k_x[i * NX:(i + 1) * NX]
        assert int(np.argmin(np.abs(window))) == NX // 2
        assert window[0] < 0 < window[-1], 'the readout does not sweep through the k-space origin'


# ------------------------------------------------------------ the contract a consumer reads back
def test_the_written_definitions_describe_this_acquisition(built, emitted) -> None:
    """
    Downstream work reads these, so they are part of the artifact rather than metadata.

    Geometry is checked against the specification; TE and TR against what the kernel computed, which
    is the one thing here that tests the emission path rather than the protocol.
    """
    kernel = built['kernel']
    fov = np.asarray(emitted.definitions['FOV'])

    assert list(np.asarray(emitted.definitions['Matrix'], dtype=int)) == [NX, NY]
    assert fov[:2] * 1e3 == pytest.approx(ref.FOV_MM, abs=1e-6)
    assert float(emitted.definitions['TE']) == pytest.approx(kernel.te_s, abs=1e-12)
    assert float(emitted.definitions['TR']) == pytest.approx(kernel.tr_s, abs=1e-12)


def test_the_build_raises_nothing_but_the_known_informational_warning(built) -> None:
    """
    Not warning-free, and the distinction matters.

    Several gradient axes ramping together exceed the scalar per-axis norm while each amplifier
    axis stays inside its own limit; the compiler says so and it is informational.  Any *other*
    `SeqCraftWarning` from this build would be a physics or legality signal and should fail here
    rather than scroll past.
    """
    kinds = {str(w.message).split(':')[0] for w in built['warnings']
             if issubclass(w.category, sc.SeqCraftWarning)}

    unexpected = {k for k in kinds if 'vector-norm limit' not in k}
    assert unexpected == set(), f'unexpected warnings from the authoring build: {unexpected}'
