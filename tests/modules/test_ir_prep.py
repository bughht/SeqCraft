"""
``IRPrep`` -- ``time_to_center`` first, because every inversion time downstream is measured from it.

The module is two events in fixed relationship, so almost nothing here is about the waveforms.  It
is about the *number*: TI is the inversion's effective centre to the acquisition of k = 0, and a
10 ms hyperbolic secant puts that centre 5 ms after the block starts.  A ``time_to_center`` that
quietly used a different origin -- the block's midpoint, say, which the spoiler makes wrong -- would
be a five-millisecond error in the one quantity an inversion-prepared sequence exists to control,
and every timeline built from it would be wrong by the same amount with nothing looking odd.

So the first two tests measure the number against the block rather than against the formula.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc

VOXEL_MM = 5.0


@pytest.fixture(scope='module')
def inv(opts):
    """The default: a 10 ms non-selective hyperbolic secant, as MPRAGE uses it."""
    return sc.modules.IRPrep(opts=opts, thickness_mm=None, spoil_voxel_mm=VOXEL_MM)


def _leaves(block):
    """Every leaf event in a built block as ``(absolute time, event)``, in time order."""
    return sorted(((t, e) for t, e, _ in sc.flatten(block)), key=lambda item: item[0])


def _rf_of(block):
    """The one RF event in a built block, reached through the nesting."""
    rf, = [e for _, e in _leaves(block) if getattr(e, 'type', '') == 'rf']
    return rf


def _gradients_of(block):
    """Every gradient event in a built block, in time order.  The spoiler nests, so flatten."""
    return [e for _, e in _leaves(block) if getattr(e, 'channel', None) is not None]


# --------------------------------------------------------------------- the number
def test_the_centre_is_where_the_pulse_peak_actually_is(opts, inv) -> None:
    """
    ``time_to_center()`` against the built block, not against the expression that computes it.

    The convention -- every ``time_to_*`` is measured from the start of the block that module's
    ``build`` returns -- is only true if something checks it, and this is the check.  The route
    here is the tree's own absolute times and the waveform's peak sample; the module's route is
    ``rf.delay + calc_rf_center(rf)``.  A hyperbolic secant's envelope peaks at its centre, so
    the two must land on the same instant.
    """
    block = inv()
    t_node, rf = next((t, e) for t, e in _leaves(block) if getattr(e, 'type', '') == 'rf')

    peak_s = t_node + float(rf.delay) + float(rf.t[int(np.argmax(np.abs(rf.signal)))])

    assert inv.time_to_center() == pytest.approx(peak_s, abs=2e-5)


def test_the_centre_falls_inside_the_pulse_the_sampler_sees(opts, inv) -> None:
    """
    The same number again, against ``sc.sample``'s RF span rather than against the waveform.

    Cheap, and it catches the one failure the test above cannot: a ``time_to_center`` that is
    self-consistent with the event but describes a block the module does not actually build.
    """
    _grid, _grads, marks = sc.sample(inv(), opts)

    spans = [(start, end) for kind, start, end, *_ in marks if kind == 'rf']

    assert len(spans) == 1
    start, end = spans[0]
    assert start < inv.time_to_center() < end


def test_the_centre_is_not_half_the_block(opts, inv) -> None:
    """
    The spoiler is *after* the pulse, so anything assuming symmetry is wrong by milliseconds.

    Stated as its own test because "halfway through the block" is the plausible wrong answer: it
    is right for a bare pulse, and this module is deliberately not one.
    """
    block = inv()

    assert inv.time_to_center() < 0.5 * block.duration - 0.5e-3
    assert inv.time_to_center() == pytest.approx(
        float(inv.rf.delay) + float(pp.calc_rf_center(inv.rf)[0]),
    )


def test_it_matches_excitations_answer_to_the_same_question(opts) -> None:
    """
    Two modules, one origin.  They compose by addition only because of that.

    An MPRAGE timeline reads ``t_shot + inv.time_to_center() + ti_s - gre.time_to_center_line()``,
    and a term measured from somewhere else would not merely be wrong -- it would be wrong in a
    way no other number in the sequence could contradict.
    """
    exc = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0)
    prep = sc.modules.IRPrep(opts=opts, thickness_mm=5.0, pulse='sinc', duration_s=3e-3)

    assert prep.time_to_center() == pytest.approx(exc.time_to_center())


# -------------------------------------------------------------------- the two events
def test_the_block_is_the_pulse_and_then_the_spoiler(opts, inv) -> None:
    block = inv()
    leaves = _leaves(block)

    assert [getattr(e, 'type', '') for _, e in leaves] == ['rf', 'trap']
    assert leaves[1][0] >= float(pp.calc_duration(inv.rf)) - 1e-12, 'the crusher is after it'
    assert {p for _, _, p in sc.flatten(block)} == {('IRPrep',), ('IRPrep', 'spoiler')}


def test_the_spoiler_winds_the_cycles_it_was_asked_for(opts, inv) -> None:
    """``area = cycles_per_voxel / voxel_m``, which is what ``spoiler`` is for."""
    crusher, = _gradients_of(inv())

    assert float(crusher.area) == pytest.approx(8.0 / (VOXEL_MM / 1e3))
    assert crusher.channel == 'z'


def test_the_spoiler_defaults_higher_than_a_readout_spoiler(opts) -> None:
    """
    Eight cycles rather than four, and the reason is the whole inversion time.

    Anything the pulse leaves behind has a TI of recovery to rephase in before it reaches the
    first readout, where it appears as a stripe that changes with TI -- which reads as a contrast
    problem rather than as a spoiling one.
    """
    assert sc.modules.IRPrep(
        opts=opts, thickness_mm=None, spoil_voxel_mm=VOXEL_MM,
    ).spoil_cycles_per_voxel == 8.0


def test_the_spoiler_axis_and_length_are_separate_from_the_selection(opts) -> None:
    prep = sc.modules.IRPrep(opts=opts, thickness_mm=None, spoil_axis='x',
                             spoil_voxel_mm=2.0, spoil_cycles_per_voxel=6.0)

    crusher, = _gradients_of(prep())

    assert crusher.channel == 'x'
    assert float(crusher.area) == pytest.approx(6.0 / 2e-3)


def test_a_non_selective_inversion_has_no_length_to_spoil_against(opts) -> None:
    """
    Refused rather than guessed.  A wrong voxel under-spoils by the ratio of the two lengths,
    and the symptom is faint residual banding that changes with TI -- exactly what the crusher
    exists to prevent, so a default would be a guess disguised as a decision.
    """
    with pytest.raises(sc.ConfigurationError, match='spoil_voxel_mm'):
        sc.modules.IRPrep(opts=opts, thickness_mm=None)


def test_a_selective_inversion_supplies_it_from_the_slab(opts) -> None:
    prep = sc.modules.IRPrep(opts=opts, thickness_mm=8.0)

    assert prep.spoil_voxel_mm == 8.0
    assert float(_gradients_of(prep())[-1].area) == pytest.approx(8.0 / 8e-3)


# --------------------------------------------------------------------- selectivity
def test_the_folder_rule_is_the_pulses_use(opts, inv) -> None:
    """``use='inversion'`` is what puts this file in ``preparation/`` rather than in ``rf/``."""
    assert inv.rf.use == 'inversion'
    assert _rf_of(inv(phase_deg=90.0)).use == 'inversion'


def test_a_non_selective_inversion_is_the_pulse_and_the_crusher(opts, inv) -> None:
    assert inv.gz is None
    assert len(inv()) == 2


def test_a_selective_inversion_gets_a_gradient_and_no_rephaser(opts) -> None:
    """
    An inversion is **not** rephased, unlike an excitation.

    ``make_adiabatic_pulse(return_gz=True)`` hands one back and it is dropped on purpose: no
    transverse magnetisation is meant to survive this block, so rephasing would restore exactly
    what the crusher in the next event is there to destroy.
    """
    prep = sc.modules.IRPrep(opts=opts, thickness_mm=8.0)

    gradients = _gradients_of(prep())

    assert prep.gz is not None
    assert len(gradients) == 2, 'the selection gradient and the spoiler, and nothing between them'
    assert float(gradients[0].amplitude) > 0.0
    assert float(gradients[1].area) == pytest.approx(
        float(prep.spoiler.nodes[0].item.area),
    ), 'the second one is the crusher, not a rephaser'


def test_the_selection_axis_defaults_to_z_and_can_move(opts) -> None:
    on_x = sc.modules.IRPrep(opts=opts, thickness_mm=8.0, axis='x', spoil_axis='y')

    assert on_x.gz.channel == 'x'
    assert [g.channel for g in _gradients_of(on_x())] == ['x', 'y']


def test_an_axis_that_cannot_take_effect_raises(opts) -> None:
    """The same shape as ``Excitation``: an argument that cannot act is reported, not ignored."""
    with pytest.raises(sc.ConfigurationError, match='axis'):
        sc.modules.IRPrep(opts=opts, thickness_mm=None, axis='z', spoil_voxel_mm=VOXEL_MM)


def test_a_selective_block_pulse_raises(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='slab-selective'):
        sc.modules.IRPrep(opts=opts, thickness_mm=8.0, pulse='block', duration_s=1e-3)


# --------------------------------------------- the identity value, and the sentinel
def test_a_zero_offset_is_valid_on_a_non_selective_pulse(opts, inv) -> None:
    """
    ``position_mm=0.0`` is what a non-selective pulse *delivers*, so it is not a refusal.

    The rule, which ``Excitation`` set: ``None`` for an argument with no identity value, and the
    identity value itself when there is one.  There is no axis meaning "no axis", so `axis` needs
    a sentinel; there is a position meaning "no offset", so making `position_mm` ``None`` would
    create two spellings of one thing and force every caller to branch on selectivity.
    """
    assert len(inv(position_mm=0.0)) == 2


def test_a_non_zero_offset_on_a_non_selective_pulse_raises(opts, inv) -> None:
    with pytest.raises(sc.ConfigurationError, match='non-selective'):
        inv(position_mm=5.0)


def test_a_slab_offset_carries_the_phase_reference(opts) -> None:
    """
    The second line of ``shift_slice``, and the one that gets omitted.

    For a 10 ms adiabatic pulse the reference is taken 5 ms into the waveform, so the term it
    contributes is large -- this is the case where dropping it is least survivable.
    """
    prep = sc.modules.IRPrep(opts=opts, thickness_mm=8.0)

    rf = _rf_of(prep(position_mm=20.0))

    expected_freq = float(prep.gz.amplitude) * 0.020
    assert rf.freq_offset == pytest.approx(expected_freq)
    assert rf.phase_offset == pytest.approx(
        -2 * np.pi * expected_freq * float(pp.calc_rf_center(prep.rf)[0]),
    )


def test_the_carrier_phase_is_a_build_argument_in_degrees(opts, inv) -> None:
    rf = _rf_of(inv(phase_deg=117.0))

    assert rf.phase_offset == pytest.approx(np.deg2rad(117.0))
    assert inv.rf.phase_offset == 0.0, 'the design was not mutated'


# ------------------------------------------------------------------- the vocabulary
@pytest.mark.parametrize(('pulse', 'duration_s'), [
    ('hypsec', 10e-3),
    ('wurst', 4e-3),
    ('sinc', 3e-3),
    ('block', 1e-3),
])
def test_every_shape_builds_an_inversion(opts, pulse, duration_s) -> None:
    prep = sc.modules.IRPrep(opts=opts, thickness_mm=None, pulse=pulse,
                             duration_s=duration_s, spoil_voxel_mm=VOXEL_MM)

    block = prep()

    assert _rf_of(block).use == 'inversion'
    assert prep.time_to_center() < block.duration


def test_the_adiabatic_pair_is_the_default_and_the_other_three_are_not(opts) -> None:
    """
    B1 insensitivity is the entire reason to use an adiabatic pulse.

    A nominal 180 that delivers 150 inverts nothing like as well, and the error is spatially
    varying -- so it appears as shading that survives every uniformity correction downstream.
    The signature says so by defaulting rather than by documenting.
    """
    assert sc.modules.IRPrep(opts=opts, thickness_mm=None,
                             spoil_voxel_mm=VOXEL_MM).pulse == 'hypsec'


def test_a_shaped_inversion_delivers_180_degrees(opts) -> None:
    """The adiabatic pair have no flip angle at all; the other three are pinned to pi."""
    prep = sc.modules.IRPrep(opts=opts, thickness_mm=None, pulse='sinc', duration_s=3e-3,
                             spoil_voxel_mm=VOXEL_MM)

    flip_rad = float(2 * np.pi * np.trapezoid(np.real(prep.rf.signal), prep.rf.t))

    assert flip_rad == pytest.approx(np.pi, rel=1e-3)


def test_an_unknown_pulse_lists_the_five(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match="'hypsec', 'wurst', 'sinc', 'slr', 'block'"):
        sc.modules.IRPrep(opts=opts, thickness_mm=None, pulse='hermite',
                          spoil_voxel_mm=VOXEL_MM)


def test_pulse_opts_reach_the_factory_and_cannot_reach_past_it(opts) -> None:
    plain = sc.modules.IRPrep(opts=opts, thickness_mm=None, spoil_voxel_mm=VOXEL_MM)
    swept = sc.modules.IRPrep(opts=opts, thickness_mm=None, spoil_voxel_mm=VOXEL_MM,
                              pulse_opts={'beta': 1200.0})

    assert not np.allclose(plain.rf.signal, swept.rf.signal)
    with pytest.raises(sc.ConfigurationError, match='accepted'):
        sc.modules.IRPrep(opts=opts, thickness_mm=None, spoil_voxel_mm=VOXEL_MM,
                          pulse_opts={'use': 'excitation'})


# -------------------------------------------------------------------- the contract
def test_it_is_pure_and_compiles_alone(opts, inv, component_checks) -> None:
    component_checks.all(inv, phase_deg=0.0)


def test_a_selective_one_is_pure_and_compiles_alone(opts, component_checks) -> None:
    component_checks.all(
        sc.modules.IRPrep(opts=opts, thickness_mm=8.0), position_mm=10.0,
    )
