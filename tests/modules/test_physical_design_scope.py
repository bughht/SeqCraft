"""
The physical-design scope, from both sides: a packaged kernel and a user's own composition.

The substrate exists because those two are the same problem.  A packaged kernel has a class and a
public API; a user composing `Excitation` + a readout + a spoiler has neither, and should not have
to promote that composition into a kernel to get physical design.  Every test here asserts that
the middle does not know which side it was handed.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.design import joint, scope

FOV_MM, MATRIX, SHOTS, THICKNESS_MM = 220.0, 64, 8, 5.0


# ----------------------------------------------------------------- the user composition
def spiral_scope(opts, *, flow_comp_axes=('x', 'y')):
    """
    Adapt a hand-written spiral repetition to the shared designer.  **A function, not a class.**

    The composition is the one `examples/gre_spiral_2d` writes -- excite, play one interleaf,
    spoil -- with `variant='in'` so the arm plays before the echo, and `prephase=False`, which is
    `SpiralReadout`'s existing way of saying "I state the moment I need at my start, you realise
    it".  That handed-over prephaser **is** the adjustable window.

    The state is an interleaf angle.  Nothing in `scope` or `joint` learns that.
    """
    exc = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=THICKNESS_MM,
                                duration_s=1e-3)
    arm = sc.modules.SpiralReadout(opts=opts, fov_mm=FOV_MM, matrix=MATRIX, shots=SHOTS,
                                   dwell_s=4e-6, variant='in')
    raster = float(opts.grad_raster_time)
    states = tuple(2.0 * np.pi * i / SHOTS for i in range(SHOTS))

    geometry = scope.ScopeGeometry(
        origin_s=exc.time_to_center(),
        window_start_s=float(sc.Raster(raster).ceil(exc().duration)),
        # `time_to_echo` is reported from a block that includes the module's OWN prephaser, and
        # this scope has taken that prephaser over -- so its duration comes back off.  See the
        # invariant in `scope.ScopeGeometry`.
        tail_s=arm.time_to_echo(0) - arm.prephaser_duration_s,
    )

    def requirement(axis: str) -> scope.AxisRequirement:
        def fixed(state, schedule):
            played = sc.LogicBlock()
            played.add(0.0, exc())
            played.add(schedule.window_start_s + schedule.window_s,
                       arm(angle_rad=float(state), prephase=False, rewind=False, acquire=False))
            return tuple(
                joint.measure_moment(played, order, axis, origin_s=schedule.origin_s,
                                     start_s=schedule.origin_s, end_s=schedule.endpoint_s)
                for order in joint.ORDERS
            )

        # The arm is one trajectory rotated by the state, so on `x` the requirement is a
        # state-independent magnitude times cos(angle) and on `y` times sin(angle).  Each is
        # bounded by that magnitude and each bound is attained -- on `x` at angle 0, on `y` at
        # pi/2 -- so the axis-aligned states are the ones worth designing against.  A proposal,
        # which `design_scope` then verifies across the whole family.
        wanted = 0.0 if axis == 'x' else np.pi / 2.0
        nearest = min(states, key=lambda a: abs(abs(np.cos(a - wanted)) - 1.0))
        return scope.AxisRequirement(
            axis=axis, states=states, design_states=(nearest,),
            target=lambda _state: (0.0, 0.0 if axis in flow_comp_axes else None),
            fixed=fixed,
        )

    def materialise(design, state, *, acquire=True):
        """One repetition, built ONCE from the finished design.  Nothing is patched afterwards."""
        out = sc.LogicBlock('spiral_tr')
        out.add(0.0, exc())
        start = design.schedule.window_start_s
        for axis in ('x', 'y'):
            # `events_for` reports instants in the schedule's frame, and this block shares it.
            for at, event in design.events_for(axis, state):
                out.add(at, event)
        out.add(start + design.window_s,
                arm(angle_rad=float(state), prephase=False, rewind=True, acquire=acquire))
        return out

    return geometry, [requirement(a) for a in ('x', 'y')], materialise, arm, states


def residual(shot, geometry, design, order: int, axis: str) -> float:
    """One moment of the whole emitted repetition, between the two semantic instants."""
    end = geometry.schedule_at(design.window_s, design.fill_s).endpoint_s
    return joint.measure_moment(shot, order, axis, origin_s=geometry.origin_s,
                                start_s=geometry.origin_s, end_s=end)


def test_a_user_composition_reaches_the_designer_without_a_kernel(opts) -> None:
    """
    The seam the substrate exists for: no `GRESpiralTR`, no Module, no public API.

    This stresses something no Cartesian kernel does.  The spiral's fixed contribution **rotates
    with the state on both axes at once**, so `x` and `y` are neither independent of the state
    nor independent of each other -- and the designer still never learns what a state is.
    """
    geometry, requirements, materialise, arm, states = spiral_scope(opts)
    design = scope.design_scope(geometry, requirements, opts,
                                min_window_s=arm.prephaser_duration_s)

    for state in states:
        shot = materialise(design, state)
        for axis in ('x', 'y'):
            assert abs(residual(shot, geometry, design, 0, axis)) < 1e-9, f'k != 0 on {axis}'
            assert abs(residual(shot, geometry, design, 1, axis)) < 1e-11, f'm1 != 0 on {axis}'


def test_the_spiral_scope_compiles_end_to_end(opts) -> None:
    """Designed, materialised once, compiled.  No finished block is rewritten on the way."""
    geometry, requirements, materialise, arm, states = spiral_scope(opts)
    design = scope.design_scope(geometry, requirements, opts,
                                min_window_s=arm.prephaser_duration_s)

    tr_s = float(sc.Raster(float(opts.block_duration_raster)).ceil(
        max(materialise(design, s).duration for s in states)))
    scan = sc.LogicBlock('gre_spiral_2d_fc')
    for index, state in enumerate(states):
        scan.add(index * tr_s, materialise(design, state))

    seq = sc.compile(scan, opts, name='gre_spiral_2d_fc')
    assert len(seq.block_events) > 0


def test_the_owned_region_must_be_taken_out_of_the_borrowed_timing(opts) -> None:
    """
    The adapter invariant, as the number it costs when it is broken.

    Semantic instants, fixed contributions and the adjustable region must describe **one**
    physical decomposition.  `SpiralReadout.time_to_echo` is measured from a block containing its
    own prephaser; a scope that takes that prephaser over and reuses the number unchanged puts
    the endpoint past the end of the arm and inside the rewinder, where `fixed` and the emitted
    block disagree about what is playing.

    Not an exception -- a residual, which is why it is worth a test rather than a comment.
    """
    geometry, requirements, materialise, arm, states = spiral_scope(opts)
    honest = scope.design_scope(geometry, requirements, opts,
                                min_window_s=arm.prephaser_duration_s)
    assert abs(residual(materialise(honest, states[0]), geometry, honest, 0, 'x')) < 1e-9

    broken = scope.ScopeGeometry(origin_s=geometry.origin_s,
                                 window_start_s=geometry.window_start_s,
                                 tail_s=geometry.tail_s + arm.prephaser_duration_s)
    wrong = scope.design_scope(broken, requirements, opts,
                               min_window_s=arm.prephaser_duration_s)
    off = abs(residual(materialise(wrong, states[0]), broken, wrong, 0, 'x'))
    assert off > 1e-3, 'the broken decomposition should be visible, not subtle'


# ------------------------------------------------------- design_states, and verification
def test_a_poor_design_state_proposal_costs_time_and_not_correctness(opts) -> None:
    """
    `design_states` is a proposal the designer checks, not metadata it trusts.

    Proposing the *easiest* state on each axis -- where the rotated requirement nearly vanishes --
    would size a window far too short.  Verification across the family widens it back to the same
    schedule the honest proposal reaches.
    """
    geometry, requirements, _materialise, arm, states = spiral_scope(opts)
    honest = scope.design_scope(geometry, requirements, opts,
                               min_window_s=arm.prephaser_duration_s)

    poor = [
        scope.AxisRequirement(axis=r.axis, states=r.states,
                              design_states=(states[2] if r.axis == 'x' else states[0],),
                              target=r.target, fixed=r.fixed)
        for r in requirements
    ]
    corrected = scope.design_scope(geometry, poor, opts, min_window_s=arm.prephaser_duration_s)

    assert corrected.window_s == honest.window_s, 'verification reaches the same schedule'


def test_the_family_is_re_realisable_at_a_longer_schedule(opts) -> None:
    """
    A protocol holding two families takes the longer echo time and asks both to design again.

    Going back through realisation rather than stretching what exists is the point: every moment
    is measured to the echo, and the echo moved.
    """
    geometry, requirements, materialise, arm, states = spiral_scope(opts)
    design = scope.design_scope(geometry, requirements, opts,
                                min_window_s=arm.prephaser_duration_s)

    raster = float(opts.grad_raster_time)
    for extra in (50, 200):
        longer = scope.realise_scope_at(geometry, requirements, opts,
                                        window_s=design.window_s + extra * raster)
        assert longer is not None
        assert longer.window_s > design.window_s
        for state in states:
            shot = materialise(longer, state)
            for axis in ('x', 'y'):
                assert abs(residual(shot, geometry, longer, 1, axis)) < 1e-11


# --------------------------------------------------------------- the packaged-kernel side
@pytest.mark.parametrize('axes', ['x', 'y', 'z', ('x', 'y'), ('x', 'y', 'z')])
def test_gre2dtr_flow_compensates_every_axis_it_owns(opts, axes) -> None:
    """
    One public spelling, three different owners, measured on the emitted repetition.

    `x` is `CartesianLine`'s own analytic solve, `y` is the shared designer, and `z` is the slice
    rephasing that the scope takes over by asking `Excitation` not to emit it.  A caller writes
    the same thing for all three.
    """
    kernel = sc.modules.GRE2DTR(opts=opts, fov_mm=FOV_MM, matrix=(32, 32),
                                thickness_mm=THICKNESS_MM,
                                flow_comp=sc.FlowCompensation(axis=axes))
    shot = kernel(line=6)
    origin, echo = kernel.exc.time_to_center(), kernel.time_to_echo()
    for axis in ((axes,) if isinstance(axes, str) else axes):
        got = joint.measure_moment(shot, 1, axis, origin_s=origin, start_s=origin, end_s=echo)
        assert abs(got) < 1e-11, f'm1 on {axis} is {got}'


def test_taking_the_slice_rephasing_over_still_rephases_the_slice(opts) -> None:
    """
    `z` compensation must not buy `m1 = 0` by leaving `k_z` somewhere other than zero.

    The two are one problem on that axis and the scope solves both, which is the reason it can
    own `z` at all: `Excitation.build(rephase=False)` hands the whole winder over rather than
    leaving a realised rephaser to be worked around.
    """
    shared = dict(opts=opts, fov_mm=FOV_MM, matrix=(32, 32), thickness_mm=THICKNESS_MM)
    plain = sc.modules.GRE2DTR(**shared)
    compensated = sc.modules.GRE2DTR(**shared, flow_comp=sc.FlowCompensation(axis='z'))

    def moments(kernel):
        origin, echo = kernel.exc.time_to_center(), kernel.time_to_echo()
        return tuple(joint.measure_moment(kernel(line=6), order, 'z', origin_s=origin,
                                          start_s=origin, end_s=echo) for order in joint.ORDERS)

    plain_m0, plain_m1 = moments(plain)
    done_m0, done_m1 = moments(compensated)

    assert abs(plain_m0) < 1e-9 and abs(done_m0) < 1e-9, 'rephased either way'
    assert abs(plain_m1) > 1e-3, 'and uncompensated it really does have a first moment'
    assert abs(done_m1) < 1e-11
    assert compensated.te_s > plain.te_s, 'which costs echo time'


def test_the_kernel_sizes_its_schedule_from_two_lines_not_every_line(opts) -> None:
    """
    Family-level timing: the search runs over the signed extremes, every line is realised at the
    schedule that comes out.

    Both ends rather than "the largest |k|", because a signed fixed contribution makes one end
    harder than the other.
    """
    kernel = sc.modules.GRE2DTR(opts=opts, fov_mm=FOV_MM, matrix=(64, 64),
                                thickness_mm=THICKNESS_MM,
                                flow_comp=sc.FlowCompensation(axis='y'))
    designed = kernel._joint['y']

    lines = {state[0] for state in designed.realisations}
    assert lines == set(range(64)), 'every line is realised'

    origin, echo = kernel.exc.time_to_center(), kernel.time_to_echo()
    for line in (0, 17, 32, 63):
        got = joint.measure_moment(kernel(line=line), 1, 'y', origin_s=origin,
                                   start_s=origin, end_s=echo)
        assert abs(got) < 1e-11, f'line {line} was realised at a schedule it does not meet'


# --------------------------------------------------- the public custom-composition surface
def spiral_pieces(opts):
    """The composition `examples/gre_spiral_2d/03_flow_comp` writes, as plain leaves."""
    exc = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=THICKNESS_MM,
                                duration_s=1e-3)
    arm = sc.modules.SpiralReadout(opts=opts, fov_mm=FOV_MM, matrix=MATRIX, shots=SHOTS,
                                   dwell_s=4e-6, variant='in')
    angles = tuple(2.0 * np.pi * i / SHOTS for i in range(SHOTS))
    return exc, arm, angles


def public_scope(opts, *, axes=('x', 'y', 'z')):
    exc, arm, angles = spiral_pieces(opts)
    return sc.PhysicalDesignScope(
        origin_s=exc.time_to_center(),
        before=exc(rephase=False),
        after=lambda angle: arm(angle_rad=angle, prephase=False),
        echo_in_after_s=arm.time_to_echo(0) - arm.prephaser_duration_s,
        axes=axes,
        states=angles,
        design_states=(angles[0], angles[SHOTS // 4]),
        min_window_s=arm.prephaser_duration_s,
    ), exc, arm, angles


def public_residual(design, scope_, shot, order, axis):
    end = scope_.origin_s + design.te_s
    return joint.measure_moment(shot, order, axis, origin_s=scope_.origin_s,
                                start_s=scope_.origin_s, end_s=end)


def test_the_public_surface_designs_a_composition_that_has_no_kernel(opts) -> None:
    """
    The whole point, in the vocabulary a user is expected to have.

    Nothing here names a schedule, a claim, a realisation family or a moment target.  The scope
    says where the spins were excited, what plays either side, which axes may be used and what
    states exist; the intent says what physics is wanted.
    """
    scope_, _exc, _arm, angles = public_scope(opts)
    design = scope_.design(opts=opts, flow_comp=sc.FlowCompensation(axis=('x', 'y', 'z')))

    for angle in angles:
        shot = design.build(angle)
        for axis in ('x', 'y', 'z'):
            assert abs(public_residual(design, scope_, shot, 0, axis)) < 1e-9
            assert abs(public_residual(design, scope_, shot, 1, axis)) < 1e-11


def test_both_paths_end_at_the_same_machinery(opts) -> None:
    """
    A packaged kernel and a user composition are two doors into one designer, not two designers.

    Asserted structurally: both produce the internal `ScopeDesign`, from the same module, with the
    same shape -- one schedule shared across axes, and a realisation per state.
    """
    kernel = sc.modules.GRE2DTR(opts=opts, fov_mm=FOV_MM, matrix=(32, 32),
                                thickness_mm=THICKNESS_MM,
                                flow_comp=sc.FlowCompensation(axis='y'))
    scope_, _exc, _arm, angles = public_scope(opts, axes=('x', 'y'))
    custom = scope_.design(opts=opts, flow_comp=sc.FlowCompensation(axis=('x', 'y')))

    assert type(kernel._joint['y']).__module__ == 'seqcraft.design.joint'
    assert type(custom._designed).__module__ == 'seqcraft.design.scope'
    assert type(custom._designed.designs['y']) is type(kernel._joint['y'])

    schedules = {d.schedule for d in custom._designed.designs.values()}
    assert len(schedules) == 1, 'one schedule shared across axes, as for the kernel'


def test_the_public_scope_re_realises_at_a_longer_echo_time(opts) -> None:
    """A protocol takes the longer of two families' echo times and asks both for it."""
    scope_, _exc, _arm, angles = public_scope(opts)
    design = scope_.design(opts=opts, flow_comp=sc.FlowCompensation(axis=('x', 'y', 'z')))

    longer = design.at(te_s=design.te_s + 1.0e-3)
    assert longer is not None and longer.te_s > design.te_s
    for angle in angles:
        shot = longer.build(angle)
        for axis in ('x', 'y', 'z'):
            assert abs(public_residual(longer, scope_, shot, 1, axis)) < 1e-11

    assert design.at(te_s=-1.0) is None, 'an impossible request is None, not a wrong answer'


def test_the_public_scope_refuses_an_axis_it_was_not_given(opts) -> None:
    """The scope owns what `axes` lists.  An intent elsewhere is a refusal, not an omission."""
    scope_, _exc, _arm, _angles = public_scope(opts, axes=('x', 'y'))

    with pytest.raises(sc.errors.ConfigurationError) as raised:
        scope_.design(opts=opts, flow_comp=sc.FlowCompensation(axis='z'))

    said = str(raised.value)
    assert 'cannot apply flow compensation' in said and "'z'" in said
    assert '`axes`' in said


def test_the_public_scope_carries_velocity_encoding_states(opts) -> None:
    """
    The other intent, on a composition with a state family of its own.

    The two multiply: eight interleaves times two encoding states, all at one echo time, and the
    difference between the states is what was asked for.
    """
    scope_, _exc, _arm, angles = public_scope(opts, axes=('x', 'y'))
    venc = sc.VelocityEncoding(venc_m_s=1.5, axis='y')
    design = scope_.design(opts=opts, velocity_encode=venc)

    for angle in angles[:3]:
        got = {state: public_residual(design, scope_,
                                      design.build(angle, encoding_state=state), 1, 'y')
               for state in venc.states}
        assert got[+1] - got[-1] == pytest.approx(venc.delta_m1_s_per_m, rel=1e-9)

    with pytest.raises(sc.errors.ConfigurationError, match='which state'):
        design.build(angles[0])


def test_a_scope_declaration_is_a_value(opts) -> None:
    """Frozen, comparable, and reusable -- the re-realisation case needs it to outlive one call."""
    first, _exc, _arm, _angles = public_scope(opts)
    second, *_ = public_scope(opts)

    assert first.axes == ('x', 'y', 'z')
    assert first.states == second.states
    with pytest.raises(Exception):
        first.origin_s = 0.0                                          # noqa: B018 -- frozen


def test_a_two_lobe_split_lands_on_the_gradient_raster(opts) -> None:
    """
    An odd number of raster steps cannot be halved onto the raster, and the compiler knows it.

    The two-lobe family split the window at exactly 0.5, so a 1450 us window asked for two 725 us
    lobes and the second began 5 us off-raster.  Nothing in the designer noticed -- the moments
    were right -- and it surfaced as a `CompileError` from the emitted sequence, which is the
    compiler doing its job but rather late.
    """
    raster = float(opts.grad_raster_time)
    for steps in (144, 145, 146, 147):                      # even and odd alike
        schedule = joint.Schedule(origin_s=0.5e-3, endpoint_s=6.0e-3,
                                  window_start_s=1.0e-3, window_s=steps * raster)
        made = joint.realise_two_lobes('y', (145.0, 0.2), (0.0, 0.0), schedule, opts)
        if made is None:
            continue
        at = schedule.window_start_s
        for event in made.events:
            assert at / raster == pytest.approx(round(at / raster), abs=1e-9), (
                f'{steps} raster steps: a lobe starts at {at * 1e6:.3f} us, off the raster'
            )
            at += float(pp.calc_duration(event))


def test_the_published_spiral_scope_compiles_on_all_three_axes(opts) -> None:
    """x + y + z through the public surface, assembled and compiled the way the notebook does."""
    scope_, _exc, _arm, angles = public_scope(opts)
    design = scope_.design(opts=opts, flow_comp=sc.FlowCompensation(axis=('x', 'y', 'z')))

    scan = sc.LogicBlock('gre_spiral_2d_flow_comp')
    for index, angle in enumerate(angles):
        scan.add(index * 30e-3, design.build(angle))

    seq = sc.compile(scan, opts, name='gre_spiral_2d_flow_comp')
    assert len(seq.block_events) > 0


# ----------------------------------------------- the user's state is never inspected
def tuple_state_scope(opts, *, axes=('x', 'y')):
    """
    A family whose state is **itself a tuple**, which is the case tuple-shape sniffing breaks.

    `(ky, kz)`, `(partition, line)` and `(angle, echo)` are all natural states.  Any internal
    scheme that reads a tuple as "state plus encoding state" gets them wrong, silently, by
    handing `before` and `after` half of what the caller passed.
    """
    exc, arm, angles = spiral_pieces(opts)
    states = tuple((angle, index) for index, angle in enumerate(angles))
    seen: list = []

    def after(state):
        seen.append(state)
        angle, _index = state
        return arm(angle_rad=angle, prephase=False)

    scope_ = sc.PhysicalDesignScope(
        origin_s=exc.time_to_center(),
        before=exc(rephase=False),
        after=after,
        echo_in_after_s=arm.time_to_echo(0) - arm.prephaser_duration_s,
        axes=axes,
        states=states,
        min_window_s=arm.prephaser_duration_s,
    )
    return scope_, states, seen


def test_a_tuple_valued_state_is_handed_back_exactly(opts) -> None:
    """The callback sees the caller's object, not a piece of it."""
    scope_, states, seen = tuple_state_scope(opts)
    design = scope_.design(opts=opts, flow_comp=sc.FlowCompensation(axis=('x', 'y')))

    assert seen, 'after was never called'
    for state in seen:
        assert isinstance(state, tuple) and len(state) == 2
        assert state in states, f'after got {state!r}, which is not one of the states given'

    design.build(states[3])


def test_a_tuple_valued_state_survives_velocity_encoding(opts) -> None:
    """
    The case that breaks tuple sniffing hardest: a tuple state **and** two encoding states.

    Velocity encoding multiplies the family internally.  If that multiplication were spelled as a
    tuple, a caller whose state is already a tuple would have their state silently unpacked --
    `after` would be handed an angle where it expected `(angle, index)`, or the encoding state
    where it expected an index.
    """
    scope_, states, seen = tuple_state_scope(opts, axes=('x', 'y'))
    venc = sc.VelocityEncoding(venc_m_s=1.5, axis='y')
    design = scope_.design(opts=opts, velocity_encode=venc)

    for state in seen:
        assert state in states, f'after got {state!r} once velocity encoding was added'

    got = {sign: public_residual(design, scope_, design.build(states[2], encoding_state=sign),
                                 1, 'y')
           for sign in venc.states}
    assert got[+1] - got[-1] == pytest.approx(venc.delta_m1_s_per_m, rel=1e-9)


def test_an_unhashable_state_is_refused_at_declaration(opts) -> None:
    """A state is a key, so it has to be one -- said plainly rather than as a later `TypeError`."""
    exc, arm, _angles = spiral_pieces(opts)
    with pytest.raises(sc.errors.ConfigurationError, match='hashable'):
        sc.PhysicalDesignScope(
            origin_s=exc.time_to_center(), before=exc(), after=arm(prephase=False),
            echo_in_after_s=arm.time_to_echo(0) - arm.prephaser_duration_s,
            axes=('x',), states=([0.0, 1.0], [2.0, 3.0]),
        )


# ------------------------------------------------------- a non-zero, state-dependent k
def encoding_scope(opts):
    """A family that **encodes** with the region it owns: a Cartesian phase encode."""
    exc = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=THICKNESS_MM,
                                duration_s=1e-3)
    pe = sc.modules.PhaseEncode(opts=opts, fov_mm=FOV_MM, matrix=32, axis='y')
    ro = sc.modules.CartesianLine(opts=opts, fov_mm=FOV_MM, matrix=32, bandwidth_hz_px=500.0,
                                  prephase=False)
    lines = tuple(range(32))
    scope_ = sc.PhysicalDesignScope(
        origin_s=exc.time_to_center(),
        before=exc(),
        after=ro(),
        echo_in_after_s=ro.time_to_echo(),
        axes=('y',),
        states=lines,
        design_states=(lines[0], lines[-1]),
        k_at_echo=lambda line, _axis: pe.k_per_m(line),
    )
    return scope_, pe, lines


def test_a_custom_scope_can_encode_with_the_region_it_owns(opts) -> None:
    """
    `k_at_echo` is the zeroth-moment half of the requirement, and it may depend on the state.

    Without it a scope could only ask for `k = 0` -- fine for a prephaser or a rephaser, and
    useless for anything that encodes.  The intents stay first-moment-only, so the two kinds of
    requirement stay in separate places.
    """
    scope_, pe, lines = encoding_scope(opts)
    design = scope_.design(opts=opts, flow_comp=sc.FlowCompensation(axis='y'))

    for line in (0, 7, 16, 31):
        shot = design.build(line)
        assert public_residual(design, scope_, shot, 0, 'y') == pytest.approx(
            pe.k_per_m(line), abs=1e-6), f'line {line} did not land on its own k'
        assert abs(public_residual(design, scope_, shot, 1, 'y')) < 1e-11


def test_an_encoding_family_still_shares_one_schedule(opts) -> None:
    """
    Every line is realised at one echo time, sized from the two signed extremes.

    A per-line timing search would give 32 different echo times for one image, which is the thing
    family-level design exists to prevent.
    """
    scope_, _pe, lines = encoding_scope(opts)
    design = scope_.design(opts=opts, flow_comp=sc.FlowCompensation(axis='y'))

    realised = design._designed.designs['y'].realisations
    assert len({key.state for key in realised}) == len(lines), 'every line realised'
    assert len({d.schedule for d in design._designed.designs.values()}) == 1


# ------------------------------------------------------------- explicit longer echo times
def test_a_longer_echo_time_goes_through_the_design_search_not_a_stretch(opts) -> None:
    """
    `.at(te_s=...)` redesigns; it does not spend the whole increase on the region.

    The extra time may become fill in front of the window, a wider window, or both -- whichever
    the family can actually be realised at.  Forcing it all into the region would be a different
    and worse answer.
    """
    scope_, _exc, _arm, angles = public_scope(opts)
    design = scope_.design(opts=opts, flow_comp=sc.FlowCompensation(axis=('x', 'y', 'z')))

    wanted = design.te_s + 0.5e-3
    longer = design.at(te_s=wanted)
    assert longer is not None

    grew = longer.window_s - design.window_s
    assert 0.0 <= grew < 0.5e-3, 'not all of the extra time went into the region'
    assert longer._designed.fill_s > 0.0, 'some of it became fill, which is the point'

    for angle in angles:
        shot = longer.build(angle)
        for axis in ('x', 'y', 'z'):
            assert abs(public_residual(longer, scope_, shot, 1, axis)) < 1e-11


def test_an_achieved_echo_time_is_never_shorter_than_the_one_asked_for(opts) -> None:
    """
    Requests between two legal instants come back at the first one at or above.

    Rounding the fill to nearest would sometimes hand back an echo earlier than requested, which
    a caller harmonising two families would then quietly violate.
    """
    scope_, _exc, _arm, _angles = public_scope(opts)
    design = scope_.design(opts=opts, flow_comp=sc.FlowCompensation(axis=('x', 'y', 'z')))

    raster = float(opts.grad_raster_time)
    for offset in (0.1 * raster, 0.5 * raster, 0.9 * raster, 7.3 * raster):
        wanted = design.te_s + offset
        longer = design.at(te_s=wanted)
        assert longer is not None
        assert longer.te_s >= wanted - 1e-12, f'{longer.te_s} < {wanted}'
        assert longer.te_s - wanted < raster, 'and not overshooting by more than one step'

    assert design.at(te_s=design.te_s - 1e-3) is None, 'a shorter request is None, not a lie'


# ------------------------------------------------------------- structural validation
def test_a_before_that_varies_in_duration_is_refused(opts) -> None:
    """
    The family shares one window start, so `before` may vary in waveform but not in length.

    Otherwise one state overlaps the designed region and another leaves a gap, and neither shows
    up as anything but a wrong moment.
    """
    exc, arm, angles = spiral_pieces(opts)
    long_exc = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=THICKNESS_MM,
                                     duration_s=2e-3)

    scope_ = sc.PhysicalDesignScope(
        origin_s=exc.time_to_center(),
        before=lambda angle: exc() if angle < 1.0 else long_exc(),
        after=lambda angle: arm(angle_rad=angle, prephase=False),
        echo_in_after_s=arm.time_to_echo(0) - arm.prephaser_duration_s,
        axes=('x', 'y'), states=angles,
    )
    with pytest.raises(sc.errors.ConfigurationError) as raised:
        scope_.design(opts=opts, flow_comp=sc.FlowCompensation(axis=('x', 'y')))

    said = str(raised.value)
    assert 'differ in duration between states' in said
    assert 'vary in waveform between states, but not in duration' in said


def test_a_before_that_varies_only_in_waveform_is_fine(opts) -> None:
    """The permitted half of the same rule, so the refusal above is not over-broad."""
    exc, arm, angles = spiral_pieces(opts)
    scope_ = sc.PhysicalDesignScope(
        origin_s=exc.time_to_center(),
        before=lambda angle: exc(phase_deg=float(np.degrees(angle)), rephase=False),
        after=lambda angle: arm(angle_rad=angle, prephase=False),
        echo_in_after_s=arm.time_to_echo(0) - arm.prephaser_duration_s,
        axes=('x', 'y', 'z'), states=angles, min_window_s=arm.prephaser_duration_s,
    )
    design = scope_.design(opts=opts, flow_comp=sc.FlowCompensation(axis=('x', 'y', 'z')))
    assert abs(public_residual(design, scope_, design.build(angles[1]), 1, 'z')) < 1e-11


@pytest.mark.parametrize('bad, match', [
    (dict(origin_s=-1e-3), 'origin_s'),
    (dict(echo_in_after_s=-1e-3), 'echo_in_after_s'),
    (dict(min_window_s=-1e-3), 'min_window_s'),
    (dict(axes=()), 'at least one axis'),
    (dict(states=()), 'at least one state'),
    (dict(design_states=()), 'nothing would size the schedule'),
    (dict(design_states=(99.0,)), 'states this scope does not have'),
])
def test_an_invalid_declaration_is_refused_where_it_is_written(opts, bad, match) -> None:
    """Refused at declaration, not discovered later as a residual or a compiler error."""
    exc, arm, angles = spiral_pieces(opts)
    good = dict(origin_s=exc.time_to_center(), before=exc(),
                after=lambda a: arm(angle_rad=a, prephase=False),
                echo_in_after_s=arm.time_to_echo(0) - arm.prephaser_duration_s,
                axes=('x', 'y'), states=angles)

    with pytest.raises(sc.errors.ConfigurationError, match=match):
        sc.PhysicalDesignScope(**{**good, **bad})


# ------------------------------------------------------- the cascade, on any window length
@pytest.mark.parametrize('steps', [60, 61, 100, 101, 145, 146, 199])
def test_whatever_the_cascade_returns_lands_on_the_gradient_raster(opts, steps) -> None:
    """
    The contract that matters, at the level that protects it.

    A family may only emit boundaries the sequence can name.  Fixing one helper is not the same
    as holding the invariant: the cascade picks whichever family fits first, so the assertion
    belongs where the choice is made.  Odd window lengths are the interesting ones -- they cannot
    be halved onto the raster.
    """
    raster = float(opts.grad_raster_time)
    schedule = joint.Schedule(origin_s=0.5e-3, endpoint_s=12.0e-3,
                              window_start_s=1.0e-3, window_s=steps * raster)
    found = 0
    for target in ((40.0, 0.005), (145.0, 0.02), (20.0, -0.01), (80.0, 0.05)):
        problem = joint.JointProblem(axis='y', targets={'only': target},
                                     fixed=lambda _s, _sched: (0.0, 0.0))
        design = joint.attempt(problem, schedule, opts)
        if design is None:
            continue
        found += 1
        at = schedule.window_start_s
        for event in design.realisations['only'].events:
            assert at / raster == pytest.approx(round(at / raster), abs=1e-9), (
                f'{steps} steps, family {design.family}: a boundary at {at * 1e6:.3f} us'
            )
            at += float(pp.calc_duration(event))
        assert at / raster == pytest.approx(round(at / raster), abs=1e-9)
    assert found, f'{steps} raster steps realised nothing, so this proved nothing'
