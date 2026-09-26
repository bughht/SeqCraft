"""
The physical-design scope, from both sides: a packaged kernel and a user's own composition.

The substrate exists because those two are the same problem.  A packaged kernel has a class and a
public API; a user composing `Excitation` + a readout + a spoiler has neither, and should not have
to promote that composition into a kernel to get physical design.  Every test here asserts that
the middle does not know which side it was handed.
"""

from __future__ import annotations

import numpy as np
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
