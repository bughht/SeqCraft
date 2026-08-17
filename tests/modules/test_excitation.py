"""
``Excitation`` -- the rephaser invariant first, because it decides what the module is.

The plan for this module opened with a question rather than a specification: does pypulseq's own
``gzr`` rewind the selection gradient's area *after the pulse's effective centre*, or half its
total?  The two agree for a symmetric pulse and disagree for every other kind, so the answer
decided whether ``Excitation`` computes its own rephaser or coordinates pypulseq's.

It coordinates.  ``test_the_rephaser_matches_the_area_after_the_effective_centre`` is the
measurement that settled it, and the minimum-phase case is why a symmetric sinc alone would not
have: its effective centre is at the *end* of the pulse, so the correct rephaser is 6.8 1/m and
"half the area" is 406.7 -- a factor of sixty, and an obvious slab of signal loss.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.modules._support import area_until


def _flip_rad(rf) -> float:
    """The flip angle a pulse actually delivers, integrated from its own waveform."""
    return float(2 * np.pi * np.trapezoid(np.real(rf.signal), rf.t))


def _area_after(gz, t_from: float) -> float:
    """Exact area of `gz` from `t_from` to its end, both measured from the block start."""
    return float(gz.area) - area_until(gz, t_from)


# ------------------------------------------------------------- the deciding measurement
@pytest.mark.parametrize(('pulse', 'pulse_opts', 'flip_deg'), [
    ('sinc', None, 15.0),
    ('gauss', None, 15.0),
    ('slr', None, 90.0),                              # ptype 'ex': a large-tip SLR design
    ('slr', {'filter_type': 'min'}, 90.0),            # minimum phase: centre at the very end
])
def test_the_rephaser_matches_the_area_after_the_effective_centre(
    opts, pulse, pulse_opts, flip_deg,
) -> None:
    """
    ``gzr.area == -(area of gz after calc_rf_center(rf))``.  Not half the total.

    The minimum-phase row is the one that distinguishes the two: for a symmetric sinc they are
    the same number, so a suite that tested only that would pass on the wrong formula.
    """
    exc = sc.modules.Excitation(opts=opts, flip_deg=flip_deg, thickness_mm=5.0,
                                pulse=pulse, pulse_opts=pulse_opts)

    after = _area_after(exc.gz, exc.time_to_center())

    assert float(exc.gzr.area) == pytest.approx(-after, rel=1e-6)


def test_the_minimum_phase_case_is_nowhere_near_half_the_area(opts) -> None:
    """The measurement that makes the previous test worth having, stated as a number."""
    asymmetric = sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0,
                                       pulse='slr', pulse_opts={'filter_type': 'min'})
    symmetric = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0)

    assert asymmetric.time_to_center() > symmetric.time_to_center() + 1e-3
    assert abs(float(asymmetric.gzr.area)) < 0.1 * abs(float(symmetric.gzr.area))


def test_the_rephaser_reports_where_it_starts_and_how_long_it_is(opts) -> None:
    """
    The two numbers a composite needs to overlap its own gradients with the rephaser.

    The rephaser is on ``z`` and occupies the tail of this block; an encode on ``y`` and a
    prephaser on ``x`` play beside it perfectly happily.  A caller that waits for the whole block
    instead waits through the rephaser for nothing, and pays for it in TE.
    """
    exc = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0)

    assert exc.time_to_rephaser() == pytest.approx(float(pp.calc_duration(exc.gz)))
    assert exc.rephaser_duration_s == pytest.approx(float(pp.calc_duration(exc.gzr)))
    assert exc.time_to_rephaser() + exc.rephaser_duration_s == pytest.approx(exc().duration)
    assert exc.time_to_rephaser() > exc.time_to_center(), 'the tail is after the pulse centre'


def test_a_non_selective_pulse_has_no_rephaser_to_overlap(opts) -> None:
    """Zero length, and the block ends where the pulse does -- so a caller needs no branch."""
    hard = sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=None,
                                 pulse='block', duration_s=500e-6)

    assert hard.rephaser_duration_s == 0.0
    assert hard.time_to_rephaser() == pytest.approx(hard().duration)


def test_the_slice_is_refocused_by_the_end_of_the_block(opts) -> None:
    """
    The whole point of the rephaser, checked through the compiler rather than the arithmetic.

    ``calculate_kspacePP`` starts k at the excitation's effective centre, which is the physics:
    spins are not transverse before it, so the selection gradient's *leading* half is not
    dephasing to undo.  A raw moment over the whole block is therefore **not** zero, and
    expecting it to be is how a correct rephaser gets "fixed".

    Without the rephaser ``k_z`` is left at the selection gradient's tail -- for a 5 mm slice a
    couple of cycles across it and most of the signal, with nothing else looking wrong.
    """
    exc = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0)
    block = exc()
    probe = sc.LogicBlock('probe').add(0.0, block).add(
        block.duration, pp.make_adc(num_samples=4, dwell=10e-6, system=opts),
    )

    k_z = sc.kspace(probe, opts)['k_adc'][2]

    assert k_z[0] == pytest.approx(0.0, abs=1e-6)


# ------------------------------------------------------------------------ selectivity
def test_a_non_selective_pulse_is_the_rf_event_alone(opts) -> None:
    """``thickness_mm=None`` drops the gradient and the rephaser, and that is all it does."""
    hard = sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=None,
                                 pulse='block', duration_s=500e-6)

    block = hard()

    assert len(block) == 1
    assert block.nodes[0].item.type == 'rf'
    assert _flip_rad(block.nodes[0].item) == pytest.approx(np.pi / 2, rel=1e-3)
    assert hard.gz is None


@pytest.mark.parametrize('pulse', ['sinc', 'gauss', 'slr'])
def test_a_shaped_pulse_without_a_gradient_is_spectrally_selective(opts, pulse) -> None:
    """A sinc or a gauss played with no gradient is how water excitation works."""
    exc = sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=None, pulse=pulse)

    assert len(exc()) == 1
    assert _flip_rad(exc().nodes[0].item) == pytest.approx(np.pi / 2, rel=1e-2)


def test_a_selective_block_pulse_raises_rather_than_surfacing_a_keyword_error(opts) -> None:
    """
    ``pp.make_block_pulse`` has neither ``slice_thickness`` nor ``return_gz``.

    A rectangular envelope's slice profile is a sinc with sidelobes nobody wants, so the refusal
    says that rather than letting pypulseq report an unexpected keyword.
    """
    with pytest.raises(sc.ConfigurationError, match='slice-selective'):
        sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0, pulse='block')


def test_the_selection_axis_defaults_to_z_and_can_move(opts) -> None:
    on_z = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0)
    on_x = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0, axis='x')

    assert {n.item.channel for n in on_z() if hasattr(n.item, 'channel')} == {'z'}
    assert {n.item.channel for n in on_x() if hasattr(n.item, 'channel')} == {'x'}
    assert float(on_x.gz.area) == pytest.approx(float(on_z.gz.area))


def test_an_axis_that_cannot_take_effect_raises(opts) -> None:
    """The same shape as ``rf_spoil_deg`` alongside ``rf_spoil=False``: reported, not ignored."""
    with pytest.raises(sc.ConfigurationError, match='axis'):
        sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=None, axis='z')


def test_a_slice_offset_needs_a_slice(opts) -> None:
    exc = sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=None)

    with pytest.raises(sc.ConfigurationError, match='non-selective'):
        exc(position_mm=20.0)


# ------------------------------------------------------------------------ the vocabulary
def test_an_unknown_pulse_lists_the_four(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match="'sinc', 'slr', 'gauss', 'block'"):
        sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0, pulse='hermite')


def test_an_unrecognised_pulse_opts_key_names_what_is_accepted(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='accepted'):
        sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0,
                              pulse='sinc', pulse_opts={'filter_type': 'pm'})


def test_pulse_opts_cannot_reach_past_the_module(opts) -> None:
    """The escape hatch is for pulse *design*; the module owns the rest of the signature."""
    with pytest.raises(sc.ConfigurationError, match='accepted'):
        sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0,
                              pulse_opts={'flip_angle': 1.0})


def test_pulse_opts_that_are_accepted_reach_the_factory(opts) -> None:
    plain = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0)
    apodized = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0,
                                     pulse_opts={'apodization': 0.5})

    assert not np.allclose(plain.rf.signal, apodized.rf.signal)


def test_time_bw_product_with_a_block_pulse_raises(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='time_bw_product'):
        sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=None,
                              pulse='block', time_bw_product=4.0)


def test_time_bw_product_sets_the_slice_bandwidth(opts) -> None:
    """A higher time-bandwidth product is a sharper profile and a stronger gradient."""
    narrow = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0,
                                   time_bw_product=2.0)
    wide = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0,
                                 time_bw_product=8.0)

    assert float(wide.gz.amplitude) == pytest.approx(4 * float(narrow.gz.amplitude), rel=1e-6)


# ------------------------------------------------------------------------ per-call state
def test_the_carrier_phase_is_a_build_argument_in_degrees(opts) -> None:
    """117 degrees is how RF spoiling is quoted, so degrees is what the argument takes."""
    exc = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0)

    rf = exc(phase_deg=117.0).nodes[0].item

    assert rf.phase_offset == pytest.approx(np.deg2rad(117.0))
    assert exc.rf.phase_offset == 0.0, 'the design was not mutated'


def test_a_slice_offset_carries_the_phase_reference(opts) -> None:
    """
    ``freq_offset = gz.amplitude * z``, and ``phase_offset -= 2*pi * freq * calc_rf_center(rf)``.

    The second line is the one that gets omitted, and multi-slice is silently wrong without it.
    """
    exc = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0)

    rf = exc(position_mm=20.0).nodes[0].item

    expected_freq = float(exc.gz.amplitude) * 0.020
    assert rf.freq_offset == pytest.approx(expected_freq)
    assert rf.phase_offset == pytest.approx(
        -2 * np.pi * expected_freq * float(pp.calc_rf_center(exc.rf)[0]),
    )


def test_an_asymmetric_pulse_gets_a_different_phase_reference(opts) -> None:
    """
    The assertion that proves the correction term is present rather than cancelling.

    At the midpoint of a symmetric pulse the term is one specific number; for a minimum-phase
    pulse at the same position it is a different one, because its effective centre is elsewhere.
    A missing term would give both the same phase (zero).
    """
    sinc = sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0)
    min_phase = sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0,
                                      pulse='slr', pulse_opts={'filter_type': 'min'})

    a = sinc(position_mm=20.0).nodes[0].item
    b = min_phase(position_mm=20.0).nodes[0].item

    assert a.freq_offset == pytest.approx(b.freq_offset), 'same gradient, same frequency'
    assert abs(a.phase_offset - b.phase_offset) > 1e-3, 'different centre, different reference'


def test_two_slices_at_plus_and_minus_z_are_symmetric(opts) -> None:
    exc = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0)

    up = exc(position_mm=20.0).nodes[0].item
    down = exc(position_mm=-20.0).nodes[0].item

    assert up.freq_offset == pytest.approx(-down.freq_offset)
    assert up.phase_offset == pytest.approx(-down.phase_offset)


def test_it_is_pure_and_compiles_alone(opts, component_checks) -> None:
    component_checks.all(
        sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0), phase_deg=117.0,
    )
