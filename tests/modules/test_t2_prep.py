"""
``T2Prep`` -- measured on the emitted sequence, never on the module's own arithmetic.

The organising rule: every claim below is read back out of the events the module emits, or out of
the compiled sequence, rather than from an attribute the module computed and could be wrong about
in the same way twice.  `prep_time_s` is the sharp case -- it is *measured* off the placement, so
a test that compared it with itself would assert nothing.  These compare it with the RF centres
pypulseq reports for the compiled file.

The second organising rule, and the reason for :data:`REALISATIONS`: **the contract is one thing
and the tip-up is a realisation of it.**  Everything that is true of a T2 preparation -- the
structure, the MLEV-4 pattern, where `prep_time_s` is measured, the single B1, the spoiler outside
the interval -- is asserted for *both* realisations, by parametrising rather than by testing the
default and trusting the other.  What legitimately differs between them is the tip-up's pulses and
its cost in time, and only that.
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

#: Everything up to the tip-up, which is the part the two realisations share.  MLEV-4 over the
#: four composites, each 90x-180y-90x: ``+ + - -``.  Degrees, modulo 360.
SHARED_PHASES = [0.0,
                 0.0, 90.0, 0.0,
                 0.0, 90.0, 0.0,
                 180.0, 270.0, 180.0,
                 180.0, 270.0, 180.0]

#: Flip angles in the same order: the tip-down, then four composites.
SHARED_FLIPS = [90.0, *(90.0, 180.0, 90.0) * 4]

#: The two tip-up realisations, as (flip, phase) pairs the emitted waveform must deliver.
#:
#: ``'simple'`` is a single hard -90, written as a 90 about -x.  ``'composite_270_360'`` is 270
#: about +x then 360 about -x, also a net -90.  Both end `prep_time_s` at the centre of their
#: **first** pulse -- the -90 in one case, the 270 in the other -- so the endpoint is
#: ``centres[13]`` either way, and no test below needs to know which realisation it is looking at
#: to find it.
REALISATIONS = {
    'simple': ((90.0, 180.0),),
    'composite_270_360': ((270.0, 0.0), (360.0, 180.0)),
}

#: The index of the pulse whose centre ends the preparation: the first pulse of the tip-up group,
#: after the tip-down and the four composites of three.
TIP_UP_INDEX = len(SHARED_FLIPS)


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


def _phase_deg(rf) -> float:
    """The carrier phase the emitted event actually carries, in degrees modulo 360."""
    return float(np.rad2deg(rf.phase_offset)) % 360.0


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


# --------------------------------------------------------------- what is contract, what is not
def test_the_default_realisation_is_the_simple_tip_up(opts) -> None:
    """
    The canonical ``-90`` is what a caller who says nothing gets.

    This is the correction the module's design rests on: the composite tip-up is an established
    way to *build* a tip-up, not part of what a T2 preparation is, so it cannot be the default
    merely because one implementation reached for it.
    """
    assert prep(opts).tip_up == 'simple'
    assert len(_rf_events(prep(opts)())) == TIP_UP_INDEX + 1


def test_the_two_realisations_differ_only_in_the_tip_up(opts) -> None:
    """
    Choosing a tip-up must not quietly redesign the rest of the preparation.

    The tip-down, the four composites and the MLEV-4 pattern are the contract; they are emitted
    identically either way, and the placement of the refocusing pulses tracks the realised
    interval rather than the realisation.
    """
    simple = _rf_events(prep(opts)())
    composite = _rf_events(prep(opts, tip_up='composite_270_360')())

    for (_, a), (_, b) in zip(simple[:TIP_UP_INDEX], composite[:TIP_UP_INDEX]):
        assert _flip_deg(a) == pytest.approx(_flip_deg(b))
        assert _phase_deg(a) == pytest.approx(_phase_deg(b))
    assert len(simple) == TIP_UP_INDEX + 1
    assert len(composite) == TIP_UP_INDEX + 2


def test_the_simple_tip_up_is_the_cheaper_one(opts) -> None:
    """
    Two pulses instead of one cost time, and time in a preparation is paid in relaxation.

    Stated as a measurement rather than as a claim about which realisation is *better*: the module
    makes no quantitative robustness claim for the composite form, so the only difference it is
    entitled to report is this one.
    """
    simple, composite = prep(opts), prep(opts, tip_up='composite_270_360')

    assert simple().duration < composite().duration
    assert simple.prep_time_s == pytest.approx(composite.prep_time_s, abs=1e-9)


def test_an_unknown_tip_up_is_refused_by_name(opts) -> None:
    """A realisation this module does not have is a configuration error, not a silent default."""
    with pytest.raises(sc.ConfigurationError, match='tip_up') as caught:
        prep(opts, tip_up='adiabatic')

    assert "'simple'" in str(caught.value)
    assert "'composite_270_360'" in str(caught.value)


# ------------------------------------------------------------------------ the mode structure
@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_the_block_is_the_mlev4_train_and_the_chosen_tip_up(opts, tip_up: str) -> None:
    """
    One tip-down, four composites of three, and the tip-up the caller named.

    The train's length is the mode: a train of a different length is `mlev2-hard` or something
    else, and this module ships one refocusing structure.
    """
    events = _rf_events(prep(opts, tip_up=tip_up)())

    expected = SHARED_FLIPS + [flip for flip, _ in REALISATIONS[tip_up]]
    assert [round(_flip_deg(rf), 1) for _, rf in events] == expected


@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_the_phases_are_mlev4_over_the_composites(opts, tip_up: str) -> None:
    """
    ``+ + - -``: the third and fourth composites are the first two turned through 180 degrees.

    Read off the emitted events, because a phase that is right in the constructor and lost on the
    way out is the failure this is for.  MLEV-4 governs the refocusing train and says nothing
    about the tip-up, which is why the pattern is identical under both realisations.
    """
    events = _rf_events(prep(opts, tip_up=tip_up)())

    expected = SHARED_PHASES + [phase for _, phase in REALISATIONS[tip_up]]
    assert [round(_phase_deg(rf), 1) for _, rf in events] == expected


@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_the_tip_up_reverses_the_tip_down(opts, tip_up: str) -> None:
    """
    Both realisations are a net ``-90`` about +x, which is the tip-down undone.

    Get a sign wrong -- in the single pulse, or between the 270 and the 360 -- and the tip-up
    continues the rotation instead of reversing it: the prepared magnetisation is destroyed rather
    than stored, and the sequence still compiles.  Summed off the emitted waveform, so it is the
    rotation the scanner would play rather than the one the table declares.
    """
    events = _rf_events(prep(opts, tip_up=tip_up)())
    tip_down = events[0][1]

    net = sum(_flip_deg(rf) * np.cos(rf.phase_offset) for _, rf in events[TIP_UP_INDEX:])
    assert net == pytest.approx(-_flip_deg(tip_down), abs=0.5)


@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_every_pulse_declares_itself_a_preparation(opts, tip_up: str) -> None:
    """``rf.use`` is what a downstream tool reads to know this is not an imaging pulse."""
    events = _rf_events(prep(opts, tip_up=tip_up)())

    assert {rf.use for _, rf in events} == {'preparation'}


# ------------------------------------------------------------------------------- the single B1
@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_one_b1_amplitude_for_the_whole_block(opts, tip_up: str) -> None:
    """A flip angle is a duration here, which is what makes the composites composites."""
    module = prep(opts, tip_up=tip_up)

    peaks = {round(float(np.abs(rf.signal).max()), 6) for _, rf in _rf_events(module())}
    assert len(peaks) == 1
    assert peaks.pop() == pytest.approx(module.b1_hz)


@pytest.mark.parametrize('requested_s', [1e-3, 1.001e-3, 0.777e-3, 1.5e-3, 2.345e-3])
@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_one_b1_survives_a_refocus_duration_that_does_not_divide_evenly(opts, tip_up: str,
                                                                       requested_s: float) -> None:
    """
    **The one-B1 claim is public, so it may not hold only for convenient numbers.**

    At one amplitude a hard pulse's flip angle is its duration, so 90 : 180 : 270 : 360 must be
    exactly 1 : 2 : 3 : 4 -- and each of those four durations has to land on the RF raster.
    Rounding each pulse separately would break the ratio for any `refocus_duration_s` that is not
    twice a raster multiple; the module rounds the *base* 90 instead and reports the 180 it
    realised.  These requests are deliberately awkward: 1.001 ms and 0.777 ms are not even
    multiples of the 1 us raster when halved.
    """
    module = prep(opts, tip_up=tip_up, refocus_duration_s=requested_s)
    raster = float(opts.rf_raster_time)

    peaks = {round(float(np.abs(rf.signal).max()), 9) for _, rf in _rf_events(module())}
    assert len(peaks) == 1, 'every pulse in the block plays at one amplitude'

    durations = sorted({round(float(rf.shape_dur), 12) for _, rf in _rf_events(module())})
    base = durations[0]
    assert durations == pytest.approx([base * n for n in (1, 2, 3, 4)][:len(durations)])
    for duration in durations:
        assert round(duration / raster, 6) == pytest.approx(round(duration / raster), abs=1e-9)
    assert module.refocus_duration_s == pytest.approx(2.0 * base)
    assert module.requested_refocus_duration_s == requested_s
    assert module.refocus_duration_s >= requested_s


# ----------------------------------------------------------------------------- the timing
@pytest.mark.parametrize('prep_time_s', [26e-3, 30e-3, 40e-3, 50e-3, 80e-3])
@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_the_preparation_interval_is_realised_on_the_compiled_sequence(opts, tip_up: str,
                                                                       prep_time_s: float) -> None:
    """
    **The claim the module exists to make good on.**

    `prep_time_s` runs between two plain RF centres -- the tip-down's, and the tip-up group's
    first pulse -- and it is checked here against the centres read back out of the *compiled*
    sequence, which is a different code path from the one that placed them.

    Both ends are ordinary pulse centres under both realisations, which is what makes the interval
    something a compiled file can be measured for.  For the composite form the endpoint is the
    270's centre: that is the convention the implementation that uses this realisation defines its
    own preparation time by, and the module adopts it rather than inventing a derived instant.
    """
    module = prep(opts, tip_up=tip_up, prep_time_s=prep_time_s)

    centres = _compiled_rf_centres(module(), opts)
    assert len(centres) == len(module.pulses), 'every pulse reached the file'

    measured = centres[TIP_UP_INDEX] - centres[0]
    assert measured == pytest.approx(module.prep_time_s, abs=1e-9)
    assert module.prep_time_s == pytest.approx(prep_time_s, abs=float(opts.grad_raster_time) / 2)


@pytest.mark.parametrize('prep_time_s', [26e-3, 45e-3, 50e-3])
@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_the_gaps_either_side_of_every_refocusing_pulse_are_equal(opts, tip_up: str,
                                                                  prep_time_s: float) -> None:
    """
    The refocusing condition, and the reason the fractions are 1/8, 3/8, 5/8 and 7/8.

    Static dephasing accumulated before a refocusing pulse is only undone if the same time passes
    after it, so the four inner intervals must be equal and the two outer ones half of each --
    to within the raster the placement is quantised onto.
    """
    module = prep(opts, tip_up=tip_up, prep_time_s=prep_time_s)

    first, *middle, last = module.interval_gaps_s
    raster = float(opts.grad_raster_time)

    assert max(middle) - min(middle) <= raster
    assert first == pytest.approx(middle[0] / 2.0, abs=raster)
    assert last == pytest.approx(middle[-1] / 2.0, abs=raster)


def test_a_composite_refocusing_pulse_is_centred_on_its_180(opts) -> None:
    """
    The refocusing instant is the 180 the composite is built around, not the group's midpoint.

    90x-180y-90x is symmetric, so for the refocusing train the two coincide -- but naming the 180
    is what makes the rule checkable rather than a coincidence of this particular composite.
    """
    module = prep(opts)
    events = _rf_events(module())

    for index in range(4):
        start, middle = events[1 + 3 * index + 1]
        centre = start + float(middle.delay) + float(pp.calc_rf_center(middle)[0])
        assert module.time_to_refocus(index) == pytest.approx(centre, abs=1e-12)


# ------------------------------------------------------------------------------ the gradients
@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_no_gradient_plays_during_the_preparation(opts, tip_up: str) -> None:
    """
    A gradient inside the preparation would dephase what the module is storing.

    This module has no selection gradient and no encoding, so the only gradient in the block is
    the crusher -- and it must lie wholly outside the transverse period.
    """
    module = prep(opts, tip_up=tip_up)

    gradients = [(start, event) for start, event, _ in flatten(module())
                 if getattr(event, 'type', None) in ('trap', 'grad')]
    assert gradients, 'the spoiler is a gradient and should be here'
    for start, _ in gradients:
        assert start >= module.time_to_tip_up()


@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_the_spoiler_is_after_the_tip_up_and_carries_the_requested_moment(opts,
                                                                         tip_up: str) -> None:
    """The crusher destroys what the tip-up left transverse, and is not part of the weighting."""
    module = prep(opts, tip_up=tip_up, spoil_cycles_per_voxel=6.0)

    assert module.time_to_spoiler() >= module.time_to_tip_up()
    # 6 cycles across a 4 mm voxel is 1500 1/m, and it is the only gradient in the block.
    assert sc.moments(module(), 0)['z'] == pytest.approx(6.0 / (VOXEL_MM / 1e3), rel=1e-9)


def test_the_spoiler_axis_is_the_callers(opts) -> None:
    module = prep(opts, spoil_axis='x')

    assert set(sc.moments(module(), 0)) == {'x'}


# ---------------------------------------------------------------------------- metamorphic
@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_doubling_the_preparation_doubles_the_interval_and_changes_nothing_else(opts,
                                                                                tip_up: str
                                                                                ) -> None:
    """
    The preparation time is a length, not a redesign: the pulses are identical either side of it.

    What must move is the placement; what must not is the pulse count, the flip angles, the
    phases or the crusher.

    Doubling holds to within the raster the placement is quantised onto.
    """
    short = prep(opts, tip_up=tip_up, prep_time_s=30e-3)
    long = prep(opts, tip_up=tip_up, prep_time_s=60e-3)

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
@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_a_preparation_shorter_than_the_pulses_is_refused(opts, tip_up: str) -> None:
    """
    Naming the floor, because "too short" without a number is a search.

    The floor depends on the realisation -- a two-pulse tip-up leaves less room for the last
    composite -- so the number the message quotes is the one this configuration computed.
    """
    module = prep(opts, tip_up=tip_up)
    with pytest.raises(sc.ConfigurationError, match='prep_time_s') as caught:
        prep(opts, tip_up=tip_up, prep_time_s=module.min_prep_time_s / 2.0)

    assert f'{module.min_prep_time_s * 1e3:g}' in str(caught.value)


@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_the_quoted_floor_is_accepted(opts, tip_up: str) -> None:
    """The rule for quoting a number: build at it and check."""
    floor = prep(opts, tip_up=tip_up).min_prep_time_s
    built = prep(opts, tip_up=tip_up, prep_time_s=floor)

    assert built.prep_time_s >= floor - float(opts.grad_raster_time)


@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_the_quoted_floor_is_the_true_one_under_unequal_dead_and_ringdown(opts,
                                                                         tip_up: str) -> None:
    """
    **`min_prep_time_s` is a public number, so it must be the real floor.**

    The floor is set by how close the outer composites can come to the tip pulses, and that is
    measured from the composite's refocusing instant -- the centre of its 180 -- to each end of
    the group.  Those two distances are not equal: a group is a run of raster-ceiled slots and
    every pulse carries a dead time in front and a ringdown behind, so the refocusing instant is
    the midpoint only by coincidence.  At 300 us dead and 20 us ringdown the two differ by 280 us,
    and solving from half the group quotes a floor about 1.1 ms too long.

    The check is that the binding seam has almost no slack left at the quoted floor -- which is
    what "true minimum" means, and what half a group would not give.
    """
    lopsided = copy.copy(opts)
    lopsided.rf_dead_time = 300e-6
    lopsided.rf_ringdown_time = 20e-6
    raster = sc.Raster(float(lopsided.grad_raster_time))

    floor = prep(lopsided, tip_up=tip_up).min_prep_time_s
    events = _rf_events(prep(lopsided, tip_up=tip_up, prep_time_s=floor)())
    ends = [start + float(raster.ceil(float(pp.calc_duration(rf)))) for start, rf in events]

    seams = [events[index + 1][0] - ends[index] for index in range(len(events) - 1)]
    assert min(seams) > -1e-12, 'nothing overlaps at the floor'
    outer = (seams[0], seams[TIP_UP_INDEX - 1])
    assert min(outer) < 2.0 * float(lopsided.grad_raster_time), 'the floor is tight, not padded'


@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_a_pulse_over_max_b1_is_refused(opts, tip_up: str) -> None:
    """
    `T2Prep` is an RF-producing path, so it holds the same peak-B1 contract as every other.

    One B1 for the whole block means the remedy is the 180's duration, and that is what the
    message names -- under either realisation, since the tip-up shares that amplitude too.
    """
    with pytest.raises(sc.ConfigurationError, match='max_b1') as caught:
        prep(opts, tip_up=tip_up, refocus_duration_s=0.2e-3)

    assert 'refocus_duration_s' in str(caught.value)


def test_an_unusable_max_b1_is_refused_like_every_other_rf_module(opts) -> None:
    unset = copy.copy(opts)
    unset.max_b1 = None

    with pytest.raises(sc.ConfigurationError, match='max_b1 is not set'):
        prep(unset)


@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_the_compiler_accepts_the_whole_block(opts, tip_up: str) -> None:
    """The contract is the emitted file, so the block has to survive compilation on its own."""
    seq = sc.compile(prep(opts, tip_up=tip_up)(), opts)

    assert len(seq.block_events) > len(prep(opts, tip_up=tip_up).pulses)


# ----------------------------------------------------------------------------- composition
@pytest.mark.parametrize('tip_up', list(REALISATIONS))
def test_the_carrier_phase_rotates_the_whole_preparation(opts, tip_up: str) -> None:
    """
    Phase-cycling the preparation across shots must not disturb the pattern inside it.

    Every *difference* -- the tip pair's inverse relationship, the MLEV progression, and the 180
    between the 270 and the 360 -- has to survive, so only the common offset moves.
    """
    plain = _rf_events(prep(opts, tip_up=tip_up)())
    turned = _rf_events(prep(opts, tip_up=tip_up)(phase_deg=117.0))

    for (_, a), (_, b) in zip(plain, turned):
        difference = float(np.rad2deg(b.phase_offset - a.phase_offset)) % 360.0
        assert difference == pytest.approx(117.0, abs=1e-9)
