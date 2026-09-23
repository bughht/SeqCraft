"""
``T2Prep`` -- measured on the emitted sequence, never on the module's own arithmetic.

The organising rule: every claim below is read back out of the events the module emits, or out of
the compiled sequence, rather than from an attribute the module computed and could be wrong about
in the same way twice.  `prep_time_s` is the sharp case -- it is *measured* off the placement, so
a test that compared it with itself would assert nothing.  These compare it with the RF centres
pypulseq reports for the compiled file.
"""

from __future__ import annotations

import copy

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.design.logic import flatten

#: The voxel the crusher counts cycles across.  Required: the preparation is non-selective.
VOXEL_MM = 4.0

#: MLEV-4 over the four composites, each 90x-180y-90x.  Degrees, modulo 360.
EXPECTED_PHASES = [0.0,
                   0.0, 90.0, 0.0,
                   0.0, 90.0, 0.0,
                   180.0, 270.0, 180.0,
                   180.0, 270.0, 180.0,
                   0.0, 180.0]

#: Flip angles in the same order: the tip-down, four composites, and the two-pulse tip-up.
EXPECTED_FLIPS = [90.0, *(90.0, 180.0, 90.0) * 4, 270.0, 360.0]


def prep(opts, **overrides):
    kwargs = {'opts': opts, 'prep_time_s': 50e-3, 'spoil_voxel_mm': VOXEL_MM}
    kwargs.update(overrides)
    return sc.modules.T2Prep(**kwargs)


def _rf_events(tree):
    """Every RF event in the emitted tree, with its start, in time order."""
    found = [(start, event) for start, event, _ in flatten(tree)
             if getattr(event, 'type', None) == 'rf']
    return sorted(found, key=lambda pair: pair[0])


def _flip_deg(rf) -> float:
    """The flip angle the emitted waveform actually delivers, in degrees."""
    return float(np.abs(np.trapezoid(rf.signal, rf.t)) * 360.0)


def _compiled_rf_centres(tree, opts) -> list[float]:
    """
    Effective RF centres read back out of the **compiled** sequence, in seconds.

    Deliberately not ``sc.kspace``: it reports ``t_excitation`` and ``t_refocusing``, and a
    preparation pulse is neither, so it sees nothing here.  Walking the emitted blocks and
    accumulating their durations is a different path from the arithmetic that placed them, which
    is the point -- it is what makes the timing assertions evidence rather than a restatement.
    """
    seq = sc.compile(tree, opts)
    centres, elapsed = [], 0.0
    for index in sorted(seq.block_events):
        block = seq.get_block(index)
        rf = getattr(block, 'rf', None)
        if rf is not None:
            centres.append(elapsed + float(rf.delay) + float(pp.calc_rf_center(rf)[0]))
        elapsed += float(seq.block_durations[index])
    return centres


# ------------------------------------------------------------------------ the mode structure
def test_the_block_is_fifteen_hard_pulses_in_the_mlev4_pattern(opts) -> None:
    """
    One tip-down, four composites of three, and a two-pulse tip-up.

    The count is the mode: a train of a different length is `mlev2-hard` or something else, and
    this module ships one mode.
    """
    module = prep(opts)

    events = _rf_events(module())

    assert len(events) == 15
    assert [round(_flip_deg(rf), 1) for _, rf in events] == EXPECTED_FLIPS


def test_every_pulse_declares_itself_a_preparation(opts) -> None:
    """``rf.use`` is what a downstream tool reads to know this is not an imaging pulse."""
    events = _rf_events(prep(opts)())

    assert {rf.use for _, rf in events} == {'preparation'}


def test_the_phases_are_mlev4_over_the_composites(opts) -> None:
    """
    ``+ + - -``: the third and fourth composites are the first two turned through 180 degrees.

    Read off the emitted events, because a phase that is right in the constructor and lost on the
    way out is the failure this is for.
    """
    events = _rf_events(prep(opts)())

    measured = [round(float(np.rad2deg(rf.phase_offset)) % 360.0, 1) for _, rf in events]
    assert measured == EXPECTED_PHASES


def test_the_tip_up_reverses_the_tip_down(opts) -> None:
    """
    270 about +x then 360 about -x is a net -90 about +x, which is the tip-down undone.

    Get the sign wrong and the second pulse continues the rotation instead of reversing it: the
    prepared magnetisation is destroyed rather than stored, and the sequence still compiles.
    """
    events = _rf_events(prep(opts)())
    tip_down, (up_a, up_b) = events[0][1], (events[13][1], events[14][1])

    net = _flip_deg(up_a) * np.cos(up_a.phase_offset) + _flip_deg(up_b) * np.cos(up_b.phase_offset)
    assert net == pytest.approx(-_flip_deg(tip_down), abs=0.5)


def test_one_b1_amplitude_for_the_whole_block(opts) -> None:
    """A flip angle is a duration here, which is what makes the composites composites."""
    events = _rf_events(prep(opts)())

    peaks = {round(float(np.abs(rf.signal).max()), 6) for _, rf in events}
    assert len(peaks) == 1
    assert peaks.pop() == pytest.approx(prep(opts).b1_hz)


# ----------------------------------------------------------------------------- the timing
@pytest.mark.parametrize('prep_time_s', [26e-3, 30e-3, 40e-3, 50e-3, 80e-3])
def test_the_preparation_interval_is_realised_on_the_compiled_sequence(opts,
                                                                       prep_time_s: float) -> None:
    """
    **The claim the module exists to make good on.**

    `prep_time_s` is the tip-down to tip-up interval, and it is checked against the RF centres
    ``pp.calculate_kspacePP`` reports for the compiled file -- a different code path from the one
    that placed them.  The tip-up is a composite, so its effective instant is the
    amplitude-weighted centre of the group; the first and last centres below bracket the whole
    preparation.
    """
    module = prep(opts, prep_time_s=prep_time_s)

    centres = _compiled_rf_centres(module(), opts)
    weights = [_flip_deg(rf) for _, rf in _rf_events(module())]
    assert len(centres) == 15, 'every pulse reached the file'
    tip_down = centres[0]
    tip_up = float(np.average(centres[13:], weights=weights[13:]))

    assert tip_up - tip_down == pytest.approx(module.prep_time_s, abs=1e-9)
    assert module.prep_time_s == pytest.approx(prep_time_s, abs=float(opts.grad_raster_time) / 2)


@pytest.mark.parametrize('prep_time_s', [26e-3, 45e-3, 50e-3])
def test_the_gaps_either_side_of_every_refocusing_pulse_are_equal(opts,
                                                                  prep_time_s: float) -> None:
    """
    The refocusing condition, and the reason the fractions are 1/8, 3/8, 5/8 and 7/8.

    Static dephasing accumulated before a refocusing pulse is only undone if the same time passes
    after it, so the four inner intervals must be equal and the two outer ones half of each --
    to within the raster the placement is quantised onto.
    """
    module = prep(opts, prep_time_s=prep_time_s)

    first, *middle, last = module.interval_gaps_s
    raster = float(opts.grad_raster_time)

    assert max(middle) - min(middle) <= raster
    assert first == pytest.approx(middle[0] / 2.0, abs=raster)
    assert last == pytest.approx(middle[-1] / 2.0, abs=raster)


def test_a_composite_refocusing_pulse_is_centred_on_its_180(opts) -> None:
    """
    The amplitude-weighted centre reduces to the obvious answer for a symmetric composite.

    90x-180y-90x is symmetric, so its weighted centre is the middle pulse's centre.  Stating the
    weighted rule rather than "the middle one" is what lets the asymmetric tip-up use it too.
    """
    module = prep(opts)
    events = _rf_events(module())

    for index in range(4):
        _, middle = events[1 + 3 * index + 1]
        start, _ = events[1 + 3 * index + 1]
        centre = start + float(middle.delay) + float(pp.calc_rf_center(middle)[0])
        assert module.time_to_refocus(index) == pytest.approx(centre, abs=1e-12)


# ------------------------------------------------------------------------------ the gradients
def test_no_gradient_plays_during_the_preparation(opts) -> None:
    """
    A gradient inside the preparation would dephase what the module is storing.

    This mode has no selection gradient and no encoding, so the only gradient in the block is the
    crusher -- and it must lie wholly outside the transverse period.
    """
    module = prep(opts)
    tree = module()

    gradients = [(start, event) for start, event, _ in flatten(tree)
                 if getattr(event, 'type', None) in ('trap', 'grad')]
    assert gradients, 'the spoiler is a gradient and should be here'
    for start, _ in gradients:
        assert start >= module.time_to_tip_up()


def test_the_spoiler_is_after_the_tip_up_and_carries_the_requested_moment(opts) -> None:
    """The crusher destroys what the tip-up left transverse, and is not part of the weighting."""
    module = prep(opts, spoil_cycles_per_voxel=6.0)

    assert module.time_to_spoiler() >= module.time_to_tip_up()
    # 6 cycles across a 4 mm voxel is 1500 1/m, and it is the only gradient in the block.
    assert sc.moments(module(), 0)['z'] == pytest.approx(6.0 / (VOXEL_MM / 1e3), rel=1e-9)


def test_the_spoiler_axis_is_the_callers(opts) -> None:
    module = prep(opts, spoil_axis='x')

    assert set(sc.moments(module(), 0)) == {'x'}


# ---------------------------------------------------------------------------- metamorphic
def test_doubling_the_preparation_doubles_the_interval_and_changes_nothing_else(opts) -> None:
    """
    The preparation time is a length, not a redesign: the pulses are identical either side of it.

    What must move is the placement; what must not is the pulse count, the flip angles, the
    phases or the crusher.

    Doubling holds to within the raster the placement is quantised onto.  The two realised
    intervals carry the same sub-raster residue, so doubling one doubles that residue too.
    """
    short, long = prep(opts, prep_time_s=30e-3), prep(opts, prep_time_s=60e-3)

    assert long.prep_time_s == pytest.approx(2.0 * short.prep_time_s,
                                             abs=float(opts.grad_raster_time))
    for a, b in zip(_rf_events(short()), _rf_events(long())):
        assert _flip_deg(a[1]) == pytest.approx(_flip_deg(b[1]))
        assert a[1].phase_offset == pytest.approx(b[1].phase_offset)
    assert sc.moments(short(), 0) == pytest.approx(sc.moments(long(), 0))


def test_a_shorter_refocus_duration_raises_b1_and_shortens_every_pulse(opts) -> None:
    """One B1 for the block means the 180's duration sets every other pulse's."""
    slow, fast = prep(opts, refocus_duration_s=2e-3), prep(opts, refocus_duration_s=1e-3)

    assert fast.b1_hz == pytest.approx(2.0 * slow.b1_hz)
    assert fast().duration < slow().duration


# ------------------------------------------------------------------------------ the refusals
def test_a_preparation_shorter_than_the_pulses_is_refused(opts) -> None:
    """
    Naming the floor, because "too short" without a number is a search.

    The binding constraint is the last composite against the tip-up, and it is solved rather than
    found by trying.
    """
    module = prep(opts)
    with pytest.raises(sc.ConfigurationError, match='prep_time_s') as caught:
        prep(opts, prep_time_s=module.min_prep_time_s / 2.0)

    assert f'{module.min_prep_time_s * 1e3:g}' in str(caught.value)


def test_the_quoted_floor_is_accepted(opts) -> None:
    """The rule for quoting a number: build at it and check."""
    floor = prep(opts).min_prep_time_s

    assert prep(opts, prep_time_s=floor).prep_time_s >= floor - float(opts.grad_raster_time)


def test_a_pulse_over_max_b1_is_refused(opts) -> None:
    """
    `T2Prep` is an RF-producing path, so it holds the same peak-B1 contract as every other.

    One B1 for the whole block means the remedy is the 180's duration, and that is what the
    message names.
    """
    with pytest.raises(sc.ConfigurationError, match='max_b1') as caught:
        prep(opts, refocus_duration_s=0.2e-3)

    assert 'refocus_duration_s' in str(caught.value)


def test_an_unusable_max_b1_is_refused_like_every_other_rf_module(opts) -> None:
    unset = copy.copy(opts)
    unset.max_b1 = None

    with pytest.raises(sc.ConfigurationError, match='max_b1 is not set'):
        prep(unset)


def test_the_compiler_accepts_the_whole_block(opts) -> None:
    """The contract is the emitted file, so the block has to survive compilation on its own."""
    seq = sc.compile(prep(opts)(), opts)

    assert len(seq.block_events) > 15


# ----------------------------------------------------------------------------- composition
def test_the_carrier_phase_rotates_the_whole_preparation(opts) -> None:
    """
    Phase-cycling the preparation across shots must not disturb the pattern inside it.

    Every *difference* -- the tip pair's inverse relationship, the MLEV progression -- has to
    survive, so only the common offset moves.
    """
    plain = _rf_events(prep(opts)())
    turned = _rf_events(prep(opts)(phase_deg=117.0))

    for (_, a), (_, b) in zip(plain, turned):
        difference = float(np.rad2deg(b.phase_offset - a.phase_offset)) % 360.0
        assert difference == pytest.approx(117.0, abs=1e-9)
