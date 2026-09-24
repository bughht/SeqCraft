"""
The repetition-level joint design prototype -- **provisional, internal, and not a public API**.

What these tests are for: showing that one piece of machinery serves two repetition families and
two augmentations, with no `N x M` branches.  Neither kernel is modified; each supplies facts it
already computes, which is the point -- the boundary is tested before anything is rewired.

Every assertion integrates the **emitted** waveform.  Nothing asks the solver what it thinks it
built, because a solver that reports its own residual is checking its arithmetic rather than its
physics.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.modules._joint import (
    CommonModeClaim,
    DifferenceClaim,
    JointProblem,
    measure_moment,
    realise,
    resolve,
    solve,
)

#: The two encoding states of a phase-contrast acquisition.
POLARITIES = (+1, -1)


def facts(kernel, axis: str, state_keys, base_m0, *, venc_m_s=None, flow_comp=True):
    """
    Build a `JointProblem` from what a kernel already exposes.

    This is the per-family adapter the design document asks for, written out longhand so it is
    visible how little it is: a semantic origin, a semantic endpoint, a window, the moments
    already accumulated between those instants, and the targets.  **No gradient event crosses
    it.**
    """
    origin = kernel.exc.time_to_center()
    endpoint = kernel.time_to_echo()
    window_start = endpoint - kernel.ro.time_to_echo()

    claims = []
    if flow_comp:
        claims.append(CommonModeClaim(axis, 1))
    if venc_m_s is not None:
        claims.append(DifferenceClaim(
            axis, 1, sc.modules.VelocityEncode.delta_m1_for(venc_m_s), POLARITIES))

    # `resolve` works over the states a claim actually relates -- here the two encoding
    # polarities -- so it is called once per encode index.  A state key is `(index, polarity)`,
    # and the encode index is what the sequence varies while the augmentation does not see it.
    m1_targets = {}
    for index in dict.fromkeys(key[0] for key in state_keys):
        group = [key for key in state_keys if key[0] == index]
        resolved = resolve(claims, [key[1] for key in group], axis=axis, order=1,
                           base={key[1]: 0.0 for key in group})
        m1_targets.update({key: resolved[key[1]] for key in group})

    return JointProblem(
        axis=axis,
        origin_s=origin,
        endpoint_s=endpoint,
        window_start_s=window_start,
        # Nothing else plays on this axis in the probe composition, so the fixed part is zero.
        # A real kernel would integrate its own events here, exactly as section 6 of the spike
        # findings measures the excitation's `ms`.
        fixed=(0.0, 0.0),
        targets={key: (base_m0[key], m1_targets[key]) for key in state_keys},
    )


def emitted(problem, solution, state, opts):
    """The lobes for one state, placed in a block at the window start."""
    out = sc.LogicBlock('winder')
    at = problem.window_start_s
    for lobe in realise(problem, solution, state, opts):
        out.add(at, lobe)
        at += float(pp.calc_duration(lobe))
    return out


def check(problem, solution, state, opts):
    """The realised (m0, m1) for one state, integrated off the emitted events."""
    block = emitted(problem, solution, state, opts)
    return tuple(
        measure_moment(block, order, problem.axis, origin_s=problem.origin_s,
                       start_s=0.0, end_s=problem.endpoint_s)
        for order in problem.orders
    )


# ------------------------------------------------------- the same machinery, two kernel families
def test_gre2dtr_supplies_a_problem_the_shared_solver_realises(opts) -> None:
    """
    **Demonstration A, half one.** A 2D kernel's own facts, a velocity encoding, a flow
    compensation, and one shared window across every state.
    """
    kernel = sc.modules.GRE2DTR(opts=opts, fov_mm=220.0, matrix=(64, 64),
                                thickness_mm=5.0, flip_deg=15.0)
    lines = (0, 31, -32)
    states = [(line, polarity) for line in lines for polarity in POLARITIES]
    base = {state: state[0] * kernel.pe.dk_per_m for state in states}

    problem = facts(kernel, 'y', states, base, venc_m_s=1.5)
    solution = solve(problem, opts)

    for state in states:
        m0, m1 = check(problem, solution, state, opts)
        assert m0 == pytest.approx(problem.targets[state][0], abs=1e-9 * max(1.0, abs(m0)))
        assert m1 == pytest.approx(problem.targets[state][1], abs=1e-12)


def test_gre3dtr_supplies_the_same_shape_of_problem(opts) -> None:
    """
    **Demonstration A, half two.** A 3D kernel, whose extra partition term enters as part of its
    own base requirement rather than as anything the augmentation knows about.
    """
    kernel = sc.modules.GRE3DTR(opts=opts, fov_mm=(220.0, 220.0, 120.0),
                                matrix=(64, 64, 16), flip_deg=15.0)
    partitions = (0, 8, 15)
    states = [(part, polarity) for part in partitions for polarity in POLARITIES]
    base = {state: kernel.combined_z_area_per_m(state[0]) for state in states}

    problem = facts(kernel, 'z', states, base, venc_m_s=1.5)
    solution = solve(problem, opts)

    for state in states:
        m0, m1 = check(problem, solution, state, opts)
        assert m0 == pytest.approx(problem.targets[state][0], abs=1e-9 * max(1.0, abs(m0)))
        assert m1 == pytest.approx(problem.targets[state][1], abs=1e-12)


def test_neither_kernel_adapter_contains_augmentation_physics(opts) -> None:
    """
    The `N x M` test, stated as an assertion rather than as a claim in a document.

    The same adapter serves flow compensation alone, velocity encoding alone, and both together,
    with only the claim list changing.  If an augmentation needed kernel-specific physics, this
    is where a third branch would have to appear.
    """
    kernel = sc.modules.GRE2DTR(opts=opts, fov_mm=220.0, matrix=(64, 64),
                                thickness_mm=5.0, flip_deg=15.0)
    states = [(0, polarity) for polarity in POLARITIES]
    base = {state: 0.0 for state in states}

    combinations = {
        'flow comp only': dict(flow_comp=True, venc_m_s=None),
        'venc only': dict(flow_comp=False, venc_m_s=1.5),
        'both': dict(flow_comp=True, venc_m_s=1.5),
    }
    realised = {}
    for name, kwargs in combinations.items():
        problem = facts(kernel, 'y', states, base, **kwargs)
        solution = solve(problem, opts)
        realised[name] = {s: check(problem, solution, s, opts)[1] for s in states}

    delta = sc.modules.VelocityEncode.delta_m1_for(1.5)
    # Flow compensation alone: both states nulled.
    assert realised['flow comp only'][(0, +1)] == pytest.approx(0.0, abs=1e-12)
    assert realised['flow comp only'][(0, -1)] == pytest.approx(0.0, abs=1e-12)
    # Both: the difference is delta and the common mode is zero, which is the derived pair.
    both = realised['both']
    assert both[(0, +1)] - both[(0, -1)] == pytest.approx(delta, rel=1e-9)
    assert both[(0, +1)] + both[(0, -1)] == pytest.approx(0.0, abs=1e-12)


# --------------------------------------------------------------------- the requirement layer
def test_the_two_augmentations_compose_without_a_precedence_rule(opts) -> None:
    """
    Velocity encoding claims the **difference**, flow compensation the **common mode**.

    Different components of the same state-indexed vector, so composing them needs no ordering
    and the symmetric pair is derived rather than assumed.
    """
    delta = sc.modules.VelocityEncode.delta_m1_for(1.5)
    claims = [DifferenceClaim('y', 1, delta, POLARITIES), CommonModeClaim('y', 1)]

    targets = resolve(claims, POLARITIES, axis='y', order=1,
                      base=dict.fromkeys(POLARITIES, 0.0))

    assert targets[+1] == pytest.approx(+delta / 2.0)
    assert targets[-1] == pytest.approx(-delta / 2.0)
    # Order of the claims is not a precedence rule.
    assert resolve(list(reversed(claims)), POLARITIES, axis='y', order=1,
                   base=dict.fromkeys(POLARITIES, 0.0)) == targets


def test_velocity_encoding_alone_leaves_the_common_mode_free(opts) -> None:
    """
    Without flow compensation only the difference is claimed, so a non-zero base survives it.

    That is the distinction the C3 record insists on: the public quantity is the change between
    states, and the symmetric split is a realisation of it rather than the contract.
    """
    delta = sc.modules.VelocityEncode.delta_m1_for(1.5)
    base = dict.fromkeys(POLARITIES, 0.02)

    targets = resolve([DifferenceClaim('y', 1, delta, POLARITIES)], POLARITIES,
                      axis='y', order=1, base=base)

    assert targets[+1] - targets[-1] == pytest.approx(delta)
    assert np.mean(list(targets.values())) == pytest.approx(0.02), 'common mode untouched'


@pytest.mark.parametrize('duplicate', (
    [CommonModeClaim('y', 1), CommonModeClaim('y', 1)],
    [DifferenceClaim('y', 1, 0.1, POLARITIES), DifferenceClaim('y', 1, 0.2, POLARITIES)],
))
def test_two_claims_on_one_component_refuse(opts, duplicate: list) -> None:
    """A disagreement about physics is not something to resolve by precedence."""
    with pytest.raises(sc.ConfigurationError, match='claim'):
        resolve(duplicate, POLARITIES, axis='y', order=1, base=dict.fromkeys(POLARITIES, 0.0))


def test_a_claim_on_another_axis_or_order_is_not_read(opts) -> None:
    """Claims are per axis and per order, so an unrelated one changes nothing."""
    base = dict.fromkeys(POLARITIES, 0.0)
    targets = resolve([CommonModeClaim('x', 1), CommonModeClaim('y', 2)], POLARITIES,
                      axis='y', order=1, base=base)

    assert targets == base


# ---------------------------------------------------------------------------- the window
def test_one_window_serves_every_state_and_comes_from_the_limiting_one(opts) -> None:
    """
    A window that varied with the state would put TE on a k-space axis.

    The limiting state is found by enumeration, never assumed to be an index extreme -- which is
    `GRE3DTR`'s own finding, and it holds in a larger state space too.
    """
    kernel = sc.modules.GRE2DTR(opts=opts, fov_mm=220.0, matrix=(64, 64),
                                thickness_mm=5.0, flip_deg=15.0)
    lines = (-32, -16, 0, 16, 31)
    states = [(line, polarity) for line in lines for polarity in POLARITIES]
    base = {state: state[0] * kernel.pe.dk_per_m for state in states}
    problem = facts(kernel, 'y', states, base, venc_m_s=1.5)

    solution = solve(problem, opts)
    raster = float(opts.grad_raster_time)

    assert all(len(solution.areas[state]) == 2 for state in states)
    # Minimal: one raster shorter, some state no longer fits.
    from seqcraft.modules._joint import _fits
    assert _fits(problem, opts, solution.half_s)
    assert not _fits(problem, opts, solution.half_s - raster)


def test_the_solve_is_independent_of_where_the_target_came_from(opts) -> None:
    """
    The realisation layer sees a number, not an augmentation.

    Two problems with identical targets realise identically even though one was assembled from a
    velocity-encoding claim and the other written down by hand.  That is what keeps the solver
    free of `N x M` knowledge.
    """
    kernel = sc.modules.GRE2DTR(opts=opts, fov_mm=220.0, matrix=(64, 64),
                                thickness_mm=5.0, flip_deg=15.0)
    states = [(0, polarity) for polarity in POLARITIES]
    delta = sc.modules.VelocityEncode.delta_m1_for(1.5)

    from_claims = facts(kernel, 'y', states, dict.fromkeys(states, 0.0), venc_m_s=1.5)
    by_hand = JointProblem(
        axis=from_claims.axis, origin_s=from_claims.origin_s,
        endpoint_s=from_claims.endpoint_s, window_start_s=from_claims.window_start_s,
        fixed=from_claims.fixed,
        targets={state: (0.0, state[1] * delta / 2.0) for state in states},
    )

    assert solve(from_claims, opts).areas == pytest.approx(solve(by_hand, opts).areas)


def test_a_difference_claim_needs_both_of_its_states(opts) -> None:
    """
    A phase difference is between two acquisitions, so claiming one is a configuration error.

    Found by the spike rather than designed in: an acquisition that plays only one polarity
    used to raise `KeyError` from inside the resolver.
    """
    delta = sc.modules.VelocityEncode.delta_m1_for(1.5)

    with pytest.raises(sc.ConfigurationError, match='difference claim') as caught:
        resolve([DifferenceClaim('y', 1, delta, POLARITIES)], [+1],
                axis='y', order=1, base={+1: 0.0})

    assert 'both states' in str(caught.value)


# ------------------------------------------------------------------------- the timing edge
@pytest.mark.parametrize(('kernel_kind', 'axis'), (('2d', 'y'), ('3d', 'z')))
def test_the_window_growth_is_reported_and_costs_echo_time(opts, kernel_kind: str,
                                                           axis: str) -> None:
    """
    Retasking the window for a first moment costs TE, and the cost is a number, not a surprise.

    Measured rather than asserted loosely, because this is the quantity a protocol trades
    against: a compensated, velocity-encoded winder is longer than one that only winds `k`.
    """
    if kernel_kind == '2d':
        kernel = sc.modules.GRE2DTR(opts=opts, fov_mm=220.0, matrix=(64, 64),
                                    thickness_mm=5.0, flip_deg=15.0)
        indices, base_of = (0, 31, -32), lambda i: i * kernel.pe.dk_per_m
    else:
        kernel = sc.modules.GRE3DTR(opts=opts, fov_mm=(220.0, 220.0, 120.0),
                                    matrix=(64, 64, 16), flip_deg=15.0)
        indices, base_of = (0, 8, 15), kernel.combined_z_area_per_m

    states = [(index, polarity) for index in indices for polarity in POLARITIES]
    problem = facts(kernel, axis, states, {s: base_of(s[0]) for s in states}, venc_m_s=1.5)
    solution = solve(problem, opts)

    assert solution.duration_s > kernel.winder_s, 'a first moment needs a longer window'
    for state in states:
        m0, m1 = check(problem, solution, state, opts)
        assert m1 == pytest.approx(problem.targets[state][1], abs=1e-12)


@pytest.mark.parametrize('kind', ('te_s', 'tr_s'))
def test_the_kernels_already_refuse_a_time_they_cannot_reach(opts, kind: str) -> None:
    """
    The refusal the longer window would trip is the one the kernels already have.

    Retasking moves `min_te_s` and `min_tr_s` up; it does not need a new refusal path, and this
    pins the existing one so that a rewiring which silently lengthened TE instead would fail.
    """
    shared = dict(opts=opts, fov_mm=220.0, matrix=(64, 64), thickness_mm=5.0, flip_deg=15.0)
    kernel = sc.modules.GRE2DTR(**shared)
    floor = kernel.min_te_s if kind == 'te_s' else kernel.min_tr_s

    with pytest.raises(sc.ConfigurationError, match=kind) as caught:
        sc.modules.GRE2DTR(**shared, **{kind: floor / 2.0})

    assert 'shorter than' in str(caught.value)
    # And the floor itself is accepted, which is what makes it a floor rather than a guess.
    assert sc.modules.GRE2DTR(**shared, **{kind: floor}) is not None

