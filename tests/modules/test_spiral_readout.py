"""
``SpiralReadout`` -- the family relationships first, because that is what makes it one module.

The four variants are one arm traversed differently, and the tests say so directly: ``'in'`` is
tested *as the reversal of* ``'out'`` rather than by restating every number independently.  A
suite that checked each variant against its own expected values would pass just as happily on four
unrelated implementations sharing a class, which is the thing the mode contract exists to rule out.

Evidence tiers are kept apart, because only one of them has an external witness:

    out                     comparable against pulseq/pulseq writeSpiral.m -- path, turn spacing,
                            extent, hardware-limited traversal.  NOT the tail: the endpoint policy
                            differs deliberately.
    in, in-out, out-in      analytic and metamorphic only.  No independent implementation exists.

Several tests are named for the historical failure they guard, from
``tools/module_mining/plans/spiral/pr23_archaeology.md``.  Those failures all compiled and looked
right, which is why they get tests rather than comments.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc

VARIANTS = ('out', 'in', 'in-out', 'out-in')


@pytest.fixture
def opts() -> pp.Opts:
    """A scanner with a real ADC sample limit, which a spiral is the first readout here to meet."""
    return pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
                   adc_samples_limit=8192)


def spiral(opts: pp.Opts, **overrides: object) -> sc.modules.SpiralReadout:
    """The reference-ish protocol, with `overrides` applied."""
    kwargs: dict = {'opts': opts, 'fov_mm': 220.0, 'matrix': 64, 'shots': 4, 'dwell_s': 4e-6}
    kwargs.update(overrides)
    return sc.modules.SpiralReadout(**kwargs)


# ------------------------------------------------------- the family, as relationships
@pytest.mark.parametrize('variant', VARIANTS)
def test_every_arm_begins_and_ends_at_rest(opts: pp.Opts, variant: str) -> None:
    """
    The endpoint policy, which is the whole reason there is a family.

    An arm ending at full gradient time-reverses into one starting at full gradient, which no
    block can begin with -- and that is why the independent reference has no spiral-in.
    """
    module = spiral(opts, variant=variant)
    # The policy is about the gradient at the block's edges, which the emitted event carries as
    # its `first` and `last` knots.  An arbitrary gradient's *samples* sit at raster centres, so
    # the first sample is legitimately non-zero -- it is one raster into the ramp away from rest.
    for times, amplitudes in module._knots:
        assert amplitudes[0] == 0.0
        assert amplitudes[-1] == 0.0
        assert times[0] == 0.0
    waveform = np.diff(module._waveform, axis=1) / opts.grad_raster_time
    reachable = float(opts.max_slew) * float(opts.grad_raster_time)
    assert np.hypot(*waveform[:, 0]) <= reachable
    assert np.hypot(*waveform[:, -1]) <= reachable


def test_in_is_out_reversed(opts: pp.Opts) -> None:
    """
    The metamorphic check the whole family rests on, and the one with no external witness.

    Stated as a relationship rather than as numbers: if this holds, ``'in'`` needs no independent
    validation of its geometry, because it *is* ``'out'``.
    """
    out, inward = spiral(opts, variant='out'), spiral(opts, variant='in')
    assert np.allclose(inward._waveform, out._waveform[:, ::-1], atol=1e-9)


@pytest.mark.parametrize(('variant', 'arms'), [('in-out', (-1, +1)), ('out-in', (+1, -1))])
def test_two_arm_variants_are_the_same_arm_twice(opts: pp.Opts, variant: str,
                                                 arms: tuple[int, int]) -> None:
    """Each half is the one arm, forwards or backwards -- not a second design."""
    module, out = spiral(opts, variant=variant), spiral(opts, variant='out')
    half = out._waveform.shape[1]
    first, second = module._waveform[:, :half], module._waveform[:, half - 1:]
    for piece, direction in zip((first, second), arms, strict=True):
        expected = out._waveform if direction > 0 else out._waveform[:, ::-1]
        assert np.allclose(piece, expected, atol=1e-9)


@pytest.mark.parametrize('variant', ('in-out', 'out-in'))
def test_the_join_is_continuous_and_at_rest(opts: pp.Opts, variant: str) -> None:
    """
    Two arms meet at zero gradient, so they need no connector.

    That is what the endpoint policy buys, and it is why multi-echo costs a two-arm variant
    nothing.
    """
    module = spiral(opts, variant=variant)
    half = spiral(opts, variant='out')._waveform.shape[1]
    waveform = np.diff(module._waveform, axis=1) / opts.grad_raster_time
    # Within one raster of rest on both sides of the seam -- which is what "at rest" can mean on
    # a sampled waveform, and is why the two arms need no connector between them.
    reachable = float(opts.max_slew) * float(opts.grad_raster_time)
    assert np.hypot(*waveform[:, half - 2]) <= reachable
    assert np.hypot(*waveform[:, half - 1]) <= reachable


# --------------------------------------------------- origin crossings, plural from the start
#: Matrix, shots and FOV chosen to move the *tail* margin, which is what the crossing count
#: turned out to depend on.  `fov_mm=400` halves `dk` and so halves the half-step the last
#: sample has to land inside.
PROTOCOLS = [
    {'matrix': 16, 'shots': 1}, {'matrix': 16, 'shots': 4},
    {'matrix': 64, 'shots': 4}, {'matrix': 128, 'shots': 4},
    {'matrix': 64, 'shots': 1, 'fov_mm': 400.0},
    {'matrix': 128, 'shots': 4, 'fov_mm': 400.0},
]


@pytest.mark.parametrize(('variant', 'count'),
                         [('out', 1), ('in', 1), ('in-out', 1), ('out-in', 2)])
@pytest.mark.parametrize('protocol', PROTOCOLS, ids=lambda p: '-'.join(map(str, p.values())))
def test_origin_crossing_count(opts: pp.Opts, variant: str, count: int,
                               protocol: dict) -> None:
    """
    ``'out-in'`` crosses twice; the durable model is a sequence because of this one case.

    **Parametrized over protocols after a near miss.**  This assertion held at the default
    protocol and failed elsewhere: ``'in'`` reaches the origin at the *end* of the arm, and
    `_plan_adc` was sizing the acquisition by shrinking a conservative estimate and never growing
    it back, so it discarded the last ~20 us -- which for a braking arm is the last dk.  At
    `matrix=128` the final sample landed at 2.281 1/m against a half-step of 2.273, and the
    module truthfully reported no crossing for a variant whose contract promises one.

    One protocol could not have caught that, because the margin it turns on is not one the
    default protocol stresses.
    """
    module = spiral(opts, variant=variant, **protocol)
    assert len(module.origin_crossing_samples) == count


def test_crossings_are_a_sequence_even_when_there_is_one(opts: pp.Opts) -> None:
    """
    The shape that must not change meaning later.

    A scalar here would work for three variants and break on the fourth -- and it would break by
    changing what an existing attribute means, which is worse than breaking loudly.
    """
    module = spiral(opts, variant='out')
    assert isinstance(module.origin_crossing_samples, tuple)
    assert isinstance(module.origin_crossing_times, tuple)
    assert module.echo_sample(0) == module.origin_crossing_samples[0]
    with pytest.raises(sc.ConfigurationError, match='does not exist'):
        module.echo_sample(1)


@pytest.mark.parametrize('variant', VARIANTS)
def test_every_crossing_is_within_half_a_step_of_the_origin(opts: pp.Opts, variant: str) -> None:
    """Measured on the reported trajectory, which is measured on the emitted events."""
    module = spiral(opts, variant=variant)
    k = module.k_per_m()
    for index in range(len(module.origin_crossing_samples)):
        radius = float(np.hypot(*k[:, module.echo_sample(index)]))
        assert radius < 0.5 * module.dk_per_m


@pytest.mark.parametrize(('variant', 'where'),
                         [('out', 'first'), ('in', 'last'), ('in-out', 'middle')])
def test_the_crossing_falls_where_the_mode_contract_says(opts: pp.Opts, variant: str,
                                                         where: str) -> None:
    """Each variant moves k = 0, and that is the mode-specific consequence the table records."""
    module = spiral(opts, variant=variant)
    sample = module.echo_sample(0)
    last = module.num_samples - 1
    if where == 'first':
        assert sample == 0
    elif where == 'last':
        assert sample == last
    else:
        assert 0.4 * last < sample < 0.6 * last


def test_out_in_reports_two_crossings_and_a_derived_spacing(opts: pp.Opts) -> None:
    """The plural machinery is exercised at ``echoes=1``, not merely declared."""
    module = spiral(opts, variant='out-in')
    first, second = module.origin_crossing_samples
    assert first == 0
    assert second == module.num_samples - 1
    assert module.echo_spacing_s == pytest.approx(
        module.time_to_echo(1) - module.time_to_echo(0),
    )
    assert spiral(opts, variant='out').echo_spacing_s is None


def test_no_crossing_lands_in_an_adc_guard(opts: pp.Opts) -> None:
    """
    PR23 failure D: split per arm rather than per acquisition region and ``k = 0`` falls into the
    dead time between two ADC events.

    Here the arms join at rest, so they are one contiguous gradient and one acquisition region --
    and every crossing is at a real sample, which is what "inside a window" means.
    """
    module = spiral(opts, variant='in-out', matrix=96, shots=2)
    assert len(module._segments) >= 1
    for index in range(len(module.origin_crossing_samples)):
        assert 0 <= module.echo_sample(index) < module.num_samples


# ------------------------------------------------- the trajectory that plays (failure C)
@pytest.mark.parametrize('variant', VARIANTS)
def test_reported_k_matches_an_independent_measurement(opts: pp.Opts, variant: str) -> None:
    """
    PR23 failure C: report the design trajectory and a NUFFT is told the wrong sample positions.

    ``sc.kspace`` shares no code with the module and walks the compiled tree, so agreement is a
    check rather than a tautology.
    """
    module = spiral(opts, variant=variant)
    block = module()
    measured = sc.kspace(block, opts)['k_adc']
    reported = module.k_per_m()
    count = min(measured.shape[1], reported.shape[1])
    error = float(np.abs(measured[:2, :count] - reported[:, :count]).max())
    assert error < 1e-9 * module.dk_per_m, f'{error:.3g} 1/m'


def test_reported_k_follows_the_rotation(opts: pp.Opts) -> None:
    """The block is already rotated, so the reported trajectory has to be too."""
    module = spiral(opts)
    angle = 0.7
    measured = sc.kspace(module(angle_rad=angle), opts)['k_adc']
    reported = module.k_per_m(angle_rad=angle)
    count = min(measured.shape[1], reported.shape[1])
    assert np.abs(measured[:2, :count] - reported[:, :count]).max() < 1e-9 * module.dk_per_m


def test_rotation_is_a_rotation(opts: pp.Opts) -> None:
    """Metamorphic, and the same property ``RadialReadout`` asserts."""
    module = spiral(opts)
    base, turned = module.k_per_m(), module.k_per_m(angle_rad=np.pi / 3)
    cos, sin = np.cos(np.pi / 3), np.sin(np.pi / 3)
    expected = np.vstack([cos * base[0] - sin * base[1], sin * base[0] + cos * base[1]])
    assert np.abs(turned - expected).max() < 1e-9 * module.dk_per_m


# --------------------------------------------------------- hardware, measured not inferred
@pytest.mark.parametrize('variant', VARIANTS)
@pytest.mark.parametrize(('max_grad', 'max_slew'), [(40, 150), (20, 80), (80, 200)])
def test_the_emitted_waveform_respects_both_limits(variant: str, max_grad: float,
                                                   max_slew: float) -> None:
    """
    Both limits, measured off what is emitted.

    Not inferred from the traversal's constraints: the solver designs against a margin and the
    resampling onto the gradient raster can still round a segment over, which is why the module
    measures and tightens rather than trusting itself.
    """
    scanner = pp.Opts(max_grad=max_grad, grad_unit='mT/m', max_slew=max_slew, slew_unit='T/m/s',
                      rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
                      adc_samples_limit=8192)
    module = spiral(scanner, variant=variant)
    gradient, slew = module.limits()
    assert gradient <= float(scanner.max_grad)
    assert slew <= float(scanner.max_slew)


def test_the_traversal_actually_uses_the_hardware(opts: pp.Opts) -> None:
    """A design that respected the limits by staying far below them would pass the test above."""
    gradient, slew = spiral(opts).limits()
    assert slew > 0.5 * float(opts.max_slew)


# ------------------------------------------------------------- the path, and the reference
def test_k_extent_matches_the_requested_resolution(opts: pp.Opts) -> None:
    """Reference-supported: ``k_max = matrix / (2 fov)``, whatever the traversal does."""
    module = spiral(opts)
    reached = float(np.hypot(*module.k_per_m()).max())
    assert reached == pytest.approx(module.k_max_per_m, rel=2e-3)


def test_turn_spacing_follows_the_density(opts: pp.Opts) -> None:
    """
    Reference-supported, and the check PR23's failure B is about.

    ``density`` is sampled field-of-view multipliers: halving it at the edge halves the field of
    view there, which *doubles* the turn spacing.  Under the polynomial-coefficient convention
    the same tuple means the opposite, and nothing raises.
    """
    uniform = spiral(opts, density=(1.0,))
    tapered = spiral(opts, density=(1.0, 0.5))
    assert tapered.duration_s < uniform.duration_s


def test_shots_enter_the_path_not_the_schedule(opts: pp.Opts) -> None:
    """One shot of a four-shot spiral is a different curve, not the same curve rotated."""
    one, four = spiral(opts, shots=1), spiral(opts, shots=4)
    assert four.duration_s < one.duration_s
    assert one.k_max_per_m == pytest.approx(four.k_max_per_m)


# ---------------------------------------------------------------- prephase and rewind
@pytest.mark.parametrize(('variant', 'starts_out', 'ends_out'),
                         [('out', False, True), ('in', True, False),
                          ('in-out', True, True), ('out-in', False, False)])
def test_which_end_is_away_from_the_origin(opts: pp.Opts, variant: str, starts_out: bool,
                                           ends_out: bool) -> None:
    """
    One rule, evaluated per mode: whichever endpoint is away from the origin needs a moment.

    Asserted as magnitudes rather than as booleans, because "at the origin" is not exactly zero
    on an emitted waveform -- a variant that finishes at the centre finishes 0.015 `dk` away, and
    the rewinder that removes *that* is failure E's correction rather than a design requirement.
    """
    module = spiral(opts, variant=variant)
    far, near = 0.5 * module.k_max_per_m, 0.05 * module.dk_per_m
    assert (module.start_k_per_m > far) if starts_out else (module.start_k_per_m < near)
    assert (module.end_k_per_m > far) if ends_out else (module.end_k_per_m < near)
    assert module.needs_prephase is starts_out


def test_the_residual_at_a_centre_ending_variant_is_still_corrected(opts: pp.Opts) -> None:
    """
    ``'in'`` ends at the origin by design and 0.015 `dk` away in fact, and the fact wins.

    PR23 measured the same residual at 0.024 `dk` and called it invisible in every plot and a
    phase ramp across the next excitation's image.  So the correction is applied even where the
    mode table says no rewinder is required -- the table describes the design, and the emitted
    waveform is what the next excitation sees.
    """
    module = spiral(opts, variant='in')
    assert 0.0 < module.end_k_per_m < 0.05 * module.dk_per_m
    moments = sc.moments(module(), order=0)
    assert max(abs(moments.get(axis, 0.0)) for axis in 'xy') < 1e-3 * module.dk_per_m


def test_the_requirement_is_stated_even_when_not_realised(opts: pp.Opts) -> None:
    """
    ``prephase=False`` hands the realisation to a caller; it does not mean none is needed.

    The same split ``CartesianLine`` and ``Excitation`` already use.
    """
    module = spiral(opts, variant='in')
    assert module.needs_prephase
    assert module(prephase=False).duration < module(prephase=True).duration


@pytest.mark.parametrize('variant', VARIANTS)
def test_the_rewinder_returns_k_to_the_origin(opts: pp.Opts, variant: str) -> None:
    """
    PR23 failure E: a join designed in isolation leaves residual area, invisible in every plot.

    Checked on the *assembled* block's total moment, which is where the residual would be.
    """
    module = spiral(opts, variant=variant)
    moments = sc.moments(module(), order=0)
    for axis in ('x', 'y'):
        assert abs(moments.get(axis, 0.0)) < 0.05 * module.dk_per_m


# ------------------------------------------------------------------ ADC segmentation
@pytest.mark.parametrize('matrix', [48, 64, 96, 112, 128])
def test_segmentation_respects_the_sample_limit(opts: pp.Opts, matrix: int) -> None:
    """Every ADC event fits, and the segments account for every sample."""
    module = spiral(opts, matrix=matrix, shots=4, variant='in-out')
    assert sum(module._segments) == module.num_samples
    for count in module._segments:
        assert 0 < count <= opts.adc_samples_limit


def test_crossing_stays_correct_across_a_segmentation_boundary(opts: pp.Opts) -> None:
    """
    The sweep the acceptance plan asks for: walk the protocol across a segment boundary.

    Segment count is driven by the sample limit alone, so a boundary may fall anywhere --
    including on a crossing.  The contract requires only that the crossing be a real sample.
    """
    seen = set()
    for matrix in range(56, 136, 8):
        module = spiral(opts, matrix=matrix, shots=2, variant='in-out')
        seen.add(len(module._segments))
        block = module()
        sc.compile(block, opts, name=f'spiral_{matrix}')
        sample = module.echo_sample(0)
        assert 0 <= sample < module.num_samples
        assert float(np.hypot(*module.k_per_m()[:, sample])) < 0.5 * module.dk_per_m
    assert len(seen) > 1, f'no segmentation boundary was crossed; saw {seen}'


def test_a_refusal_names_the_way_out(opts: pp.Opts) -> None:
    """``adc_segments=1`` refuses rather than silently splitting, and says what would fit."""
    with pytest.raises(sc.ConfigurationError, match='will not fit'):
        spiral(opts, matrix=128, shots=1, variant='in-out', adc_segments=1)


@pytest.mark.parametrize('variant', VARIANTS)
@pytest.mark.parametrize(('matrix', 'shots'), [(96, 1), (112, 1), (128, 2), (128, 4)])
def test_the_acquisition_finishes_inside_its_gradient(opts: pp.Opts, variant: str, matrix: int,
                                                      shots: int) -> None:
    """
    The regression for the defect that took the longest to find.

    Each ADC event spends a lead delay and a trailing dead time **and** is rounded up to the
    gradient raster, so the overhead per segment is up to ``lead + trail + raster``.  Sizing the
    sample budget by subtracting ``lead + trail`` under-estimates it by one raster, the
    acquisition outlasts the gradient, the compiler holds the last block open, and the waveform
    is padded with area nobody designed -- 0.045 1/m on ``matrix=112``, refused by the m0
    contract check.

    ``matrix=112, shots=1`` was the minimal failing case: its neighbours at 111 and 113 both
    passed, and it was the only one in the range whose acquisition ran past its gradient.
    """
    module = spiral(opts, matrix=matrix, shots=shots, variant=variant)
    assert module.acquisition_end_s <= module.duration_s
    sc.compile(module(), opts, name='inside')


# ----------------------------------------------------------------------- it compiles
@pytest.mark.parametrize('variant', VARIANTS)
def test_it_compiles_to_legal_pulseq(opts: pp.Opts, variant: str) -> None:
    """Layer 2: the compiler emits it, which is a different question from the physics above."""
    module = spiral(opts, variant=variant)
    seq = sc.compile(module(), opts, name=f'spiral_{variant}')
    assert len(seq.block_events) >= 1
    assert seq.duration()[0] > 0.0


@pytest.mark.parametrize('variant', VARIANTS)
def test_a_dummy_shot_plays_the_gradients_without_acquiring(opts: pp.Opts,
                                                            variant: str) -> None:
    """Steady state needs the gradients; it does not need the samples."""
    module = spiral(opts, variant=variant)
    quiet = sc.compile(module(acquire=False), opts, name='quiet')
    assert not any(getattr(quiet.get_block(n), 'adc', None) is not None
                   for n in quiet.block_events)


# ------------------------------------------------------------------------ the refusals
def test_multi_echo_refuses_and_names_the_contract(opts: pp.Opts) -> None:
    """Not implemented, and the refusal says why rather than pretending the option exists."""
    with pytest.raises(sc.ConfigurationError, match='not implemented yet'):
        spiral(opts, echoes=2)


def test_an_unknown_variant_lists_the_four(opts: pp.Opts) -> None:
    """The four are the contract; anything else is a typo."""
    with pytest.raises(sc.ConfigurationError, match='variant must be one of'):
        spiral(opts, variant='spiral-out')


def test_density_must_be_positive_multipliers(opts: pp.Opts) -> None:
    """A zero or negative field of view is not a density, it is a typo with a plausible shape."""
    with pytest.raises(sc.ConfigurationError, match='field-of-view multipliers'):
        spiral(opts, density=(1.0, 0.0))
