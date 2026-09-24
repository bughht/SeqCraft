"""
Repetition-level joint moment design -- **provisional and internal**, not a public API.

Three things are under test, and they are deliberately separate files' worth of concern in one
place only because the whole thing is provisional:

* the requirement layer, where augmentations claim components and never meet a waveform;
* the realisation families, which own their own maths and report utilisation rather than a bool;
* the two rewired kernels, measured on the **complete emitted repetition** from the semantic
  origin to the achieved echo -- not on the joint lobes alone, which is how a leaked contribution
  from somewhere else in the repetition would hide.
"""

from __future__ import annotations

import numpy as np
import pytest

import seqcraft as sc
from seqcraft.modules._joint import (
    CommonModeClaim,
    DifferenceClaim,
    Schedule,
    measure_moment,
    realise_base_plus_bipolar,
    realise_single_lobe,
    realise_split_search,
    realise_two_lobes,
    resolve_claims,
)

#: The two encoding states of a phase-contrast acquisition.
POLARITIES = (+1, -1)

#: One VENC, used throughout so the numbers can be compared across tests.
VENC_M_S = 1.5


def delta() -> float:
    return sc.modules.VelocityEncode.delta_m1_for(VENC_M_S)


def claims(axis: str, *, flow_comp: bool = True, venc: bool = True) -> tuple:
    out = []
    if flow_comp:
        out.append(CommonModeClaim(axis, 1))
    if venc:
        out.append(DifferenceClaim(axis, 1, delta(), POLARITIES))
    return tuple(out)


def residual(kernel, axis: str, shot, wanted: tuple[float, float]) -> tuple[float, float]:
    """How far the **whole repetition** is from its target, at the achieved echo."""
    origin, echo = kernel.exc.time_to_center(), kernel.time_to_echo()
    return tuple(
        abs(measure_moment(shot, order, axis, origin_s=origin, start_s=0.0, end_s=echo)
            - wanted[order])
        for order in (0, 1)
    )


# ------------------------------------------------------------------------ the requirement layer
def test_the_two_augmentations_compose_without_a_precedence_rule() -> None:
    """Velocity encoding claims the difference, flow compensation the common mode."""
    resolved = resolve_claims(claims('y'), POLARITIES, axis='y', order=1,
                              base=dict.fromkeys(POLARITIES, 0.0))

    assert resolved[+1] == pytest.approx(+delta() / 2.0)
    assert resolved[-1] == pytest.approx(-delta() / 2.0)
    assert resolve_claims(tuple(reversed(claims('y'))), POLARITIES, axis='y', order=1,
                          base=dict.fromkeys(POLARITIES, 0.0)) == resolved


def test_velocity_encoding_alone_leaves_the_common_mode_free() -> None:
    """The public quantity is the change between states; the symmetric split is a realisation."""
    base = dict.fromkeys(POLARITIES, 0.02)

    resolved = resolve_claims(claims('y', flow_comp=False), POLARITIES, axis='y', order=1,
                              base=base)

    assert resolved[+1] - resolved[-1] == pytest.approx(delta())
    assert np.mean(list(resolved.values())) == pytest.approx(0.02)


@pytest.mark.parametrize('duplicate', (
    (CommonModeClaim('y', 1), CommonModeClaim('y', 1)),
    (DifferenceClaim('y', 1, 0.1, POLARITIES), DifferenceClaim('y', 1, 0.2, POLARITIES)),
))
def test_two_claims_on_one_component_refuse(duplicate: tuple) -> None:
    with pytest.raises(sc.ConfigurationError, match='claim'):
        resolve_claims(duplicate, POLARITIES, axis='y', order=1,
                       base=dict.fromkeys(POLARITIES, 0.0))


def test_a_difference_claim_needs_both_of_its_states() -> None:
    """A phase difference is between two acquisitions."""
    with pytest.raises(sc.ConfigurationError, match='difference claim') as caught:
        resolve_claims(claims('y', flow_comp=False), [+1], axis='y', order=1, base={+1: 0.0})

    assert 'both states' in str(caught.value)


def test_claims_are_read_per_axis_and_per_order() -> None:
    base = dict.fromkeys(POLARITIES, 0.0)

    assert resolve_claims((CommonModeClaim('x', 1), CommonModeClaim('y', 2)), POLARITIES,
                          axis='y', order=1, base=base) == base


# --------------------------------------------------------------------- the realisation families
def a_schedule(window_s: float) -> Schedule:
    return Schedule(origin_s=0.5e-3, endpoint_s=3.0e-3, window_start_s=1.0e-3,
                    window_s=window_s)


def test_a_single_lobe_serves_an_unconstrained_first_moment_and_declines_a_wrong_one(
        opts) -> None:
    """
    One degree of freedom, so it serves an `m0` target and declines to guess at `m1`.

    **A target of zero is not the same as no target**, and this is where conflating them would
    show: the lobe an `m0` target implies has a first moment of ``area * centroid``, which is not
    zero, so a family that only asked "was an m1 requested?" would report success while emitting
    the wrong waveform.  Asking what it *would* produce is the difference.
    """
    schedule = a_schedule(600e-6)

    unconstrained = realise_single_lobe('y', (50.0, None), (0.0, 0.0), schedule, opts)
    assert unconstrained is not None and unconstrained.feasible

    assert realise_single_lobe('y', (50.0, 0.0), (0.0, 0.0), schedule, opts) is None

    # The one m1 it can serve is the one it makes anyway.
    centre = schedule.window_start_s + schedule.window_s / 2.0 - schedule.origin_s
    exact = realise_single_lobe('y', (50.0, 50.0 * centre), (0.0, 0.0), schedule, opts)
    assert exact is not None and exact.feasible


def test_the_two_bases_span_the_same_set(opts) -> None:
    """
    Adjacent lobes and base-plus-bipolar reach the same targets in the same window.

    They are two coordinate systems for one two-dimensional space, so the decoupled basis is
    better conditioned rather than more capable -- worth knowing before treating it as a fallback.
    """
    schedule = a_schedule(760e-6)
    target = (31 * 1e3 / 220.0, +delta() / 2.0)

    adjacent = realise_two_lobes('y', target, (0.0, 0.0), schedule, opts)
    decoupled = realise_base_plus_bipolar('y', target, (0.0, 0.0), schedule, opts)

    assert adjacent is not None and decoupled is not None
    assert adjacent.feasible == decoupled.feasible
    assert adjacent.peak_grad == pytest.approx(decoupled.peak_grad, rel=1e-9)


def test_the_shape_search_finds_more_headroom_than_a_fixed_split(opts) -> None:
    """
    The richer family is what keeps a narrow one from inventing a feasibility cliff.

    Equal durations are one point of a one-parameter family; letting the split move buys real
    headroom at the same window, which is the difference between a true physical limit and an
    artefact of the shape that was picked.
    """
    schedule = a_schedule(740e-6)
    target = (31 * 1e3 / 220.0, +delta() / 2.0)

    fixed = realise_two_lobes('y', target, (0.0, 0.0), schedule, opts)
    searched = realise_split_search('y', target, (0.0, 0.0), schedule, opts)

    assert searched is not None and searched.feasible
    assert max(searched.peak_grad, searched.peak_slew) <= max(fixed.peak_grad, fixed.peak_slew)
    assert searched.family != fixed.family, 'the searched split is reported, not hidden'


def test_a_family_reports_utilisation_and_what_limits_it(opts) -> None:
    """A bool cannot distinguish "failed by two per cent" from "failed everywhere"."""
    target = (31 * 1e3 / 220.0, +delta() / 2.0)

    tight = realise_two_lobes('y', target, (0.0, 0.0), a_schedule(300e-6), opts)

    assert tight is None or tight.limiting in ('gradient', 'slew')
    if tight is not None:
        assert tight.peak_grad > 0.0 and tight.peak_slew > 0.0


# -------------------------------------------------------------- the kernels, whole repetition
def test_gre2dtr_meets_its_targets_on_the_complete_repetition(opts) -> None:
    """
    **Demonstration A, 2D.** Measured from the excitation's instant to the achieved echo, over
    every gradient the repetition plays on the axis -- not over the joint lobes alone.
    """
    kernel = sc.modules.GRE2DTR(opts=opts, fov_mm=220.0, matrix=(32, 32), thickness_mm=5.0,
                                flip_deg=15.0, joint_claims=claims('y'),
                                encoding_states=POLARITIES)

    for line in (0, 16, 31):
        for state in POLARITIES:
            shot = kernel(line=line, encoding_state=state)
            wanted = (kernel.pe.k_per_m(line), state * delta() / 2.0)
            off_m0, off_m1 = residual(kernel, 'y', shot, wanted)
            assert off_m0 < 1e-6 * max(1.0, abs(wanted[0]))
            assert off_m1 < 1e-11


@pytest.mark.parametrize('slab_thickness_mm', (None, 120.0))
def test_gre3dtr_meets_its_targets_including_the_slab_contribution(opts,
                                                                   slab_thickness_mm) -> None:
    """
    **Demonstration A, 3D**, and the case where `fixed` earns its keep.

    A selective slab plays its own lobe on `z` and carries a first moment from the excitation's
    effective instant -- the fine scan's `ms`. It enters as a fixed contribution the kernel
    integrates from events it already places, so no leaf publishes a moment.
    """
    shared = dict(opts=opts, fov_mm=(220.0, 220.0, 120.0), matrix=(32, 32, 8), flip_deg=15.0)
    if slab_thickness_mm is not None:
        shared['slab_thickness_mm'] = slab_thickness_mm
    kernel = sc.modules.GRE3DTR(**shared, joint_claims=claims('z'),
                                encoding_states=POLARITIES)

    for partition in range(kernel.matrix[2]):
        for state in POLARITIES:
            shot = kernel(line=16, partition=partition, encoding_state=state)
            wanted = (kernel.pe_z.k_per_m(partition), state * delta() / 2.0)
            off_m0, off_m1 = residual(kernel, 'z', shot, wanted)
            assert off_m0 < 1e-6 * max(1.0, abs(wanted[0]))
            assert off_m1 < 1e-11


def test_the_kernel_adapter_does_not_branch_on_the_augmentation(opts) -> None:
    """
    The `N x M` test as an assertion: three claim sets, one code path.

    Flow compensation alone, velocity encoding alone, and both together differ only in what is
    handed to the constructor. A kernel that had to know which augmentation was in play would
    need a branch here to pass.
    """
    shared = dict(opts=opts, fov_mm=220.0, matrix=(32, 32), thickness_mm=5.0, flip_deg=15.0)
    realised = {}
    for name, kwargs in (('flow comp', dict(flow_comp=True, venc=False)),
                         ('venc', dict(flow_comp=False, venc=True)),
                         ('both', dict(flow_comp=True, venc=True))):
        kernel = sc.modules.GRE2DTR(**shared, joint_claims=claims('y', **kwargs),
                                    encoding_states=POLARITIES)
        origin, echo = kernel.exc.time_to_center(), kernel.time_to_echo()
        realised[name] = {
            state: measure_moment(kernel(line=8, encoding_state=state), 1, 'y',
                                  origin_s=origin, start_s=0.0, end_s=echo)
            for state in POLARITIES
        }

    assert realised['flow comp'][+1] == pytest.approx(0.0, abs=1e-11)
    assert realised['flow comp'][-1] == pytest.approx(0.0, abs=1e-11)
    assert realised['both'][+1] - realised['both'][-1] == pytest.approx(delta(), rel=1e-9)
    assert realised['both'][+1] + realised['both'][-1] == pytest.approx(0.0, abs=1e-11)


# ------------------------------------------------------------------------ product-like timing
def test_enabling_the_feature_moves_the_automatic_minima(opts) -> None:
    """
    AUTO timing recomputes; it does not quietly keep the old TE and hope.

    A longer winder is a later echo, and TR follows TE -- both are reported, and both are what a
    protocol trades against.
    """
    shared = dict(opts=opts, fov_mm=220.0, matrix=(32, 32), thickness_mm=5.0, flip_deg=15.0)
    plain = sc.modules.GRE2DTR(**shared)
    joint = sc.modules.GRE2DTR(**shared, joint_claims=claims('y'),
                               encoding_states=POLARITIES)

    assert joint.winder_s > plain.winder_s
    assert joint.min_te_s > plain.min_te_s
    assert joint.min_tr_s > plain.min_tr_s
    assert joint.te_s == pytest.approx(joint.min_te_s), 'te_s=None means the new minimum'


@pytest.mark.parametrize('kind', ('te_s', 'tr_s'))
def test_an_explicit_time_that_the_feature_made_infeasible_is_refused(opts, kind: str) -> None:
    """
    **The one thing that must not happen silently.**

    A TE that was legal without the feature and is not legal with it has to raise, naming the new
    minimum -- never be rounded up behind the caller's back.  The old minimum is exactly such a
    request, which is what makes this a real case rather than an arbitrary small number.
    """
    shared = dict(opts=opts, fov_mm=220.0, matrix=(32, 32), thickness_mm=5.0, flip_deg=15.0)
    plain = sc.modules.GRE2DTR(**shared)
    joint = dict(joint_claims=claims('y'), encoding_states=POLARITIES)
    was_legal = plain.min_te_s if kind == 'te_s' else plain.min_tr_s

    with pytest.raises(sc.ConfigurationError, match=kind) as caught:
        sc.modules.GRE2DTR(**shared, **joint, **{kind: was_legal})

    assert 'shorter than' in str(caught.value)


@pytest.mark.parametrize('kind', ('te_s', 'tr_s'))
def test_the_new_minimum_itself_is_accepted(opts, kind: str) -> None:
    """The rule for quoting a number: build at it and check."""
    shared = dict(opts=opts, fov_mm=220.0, matrix=(32, 32), thickness_mm=5.0, flip_deg=15.0)
    joint = dict(joint_claims=claims('y'), encoding_states=POLARITIES)
    floor = getattr(sc.modules.GRE2DTR(**shared, **joint), f'min_{kind[:2]}_s')

    built = sc.modules.GRE2DTR(**shared, **joint, **{kind: floor})

    assert getattr(built, kind) == pytest.approx(floor, abs=1e-9)


def test_nothing_changes_when_no_claim_is_made(opts) -> None:
    """
    The local path is untouched, which is the other half of the architecture.

    A kernel with no claims designs exactly what it designed before -- same winder, same TE, same
    events -- so the shared machinery costs nothing where it is not used.
    """
    shared = dict(opts=opts, fov_mm=220.0, matrix=(32, 32), thickness_mm=5.0, flip_deg=15.0)

    plain = sc.modules.GRE2DTR(**shared)
    also_plain = sc.modules.GRE2DTR(**shared, joint_claims=(), encoding_states=())

    assert also_plain.winder_s == plain.winder_s
    assert also_plain.te_s == plain.te_s
    assert also_plain(line=3).duration == plain(line=3).duration


def test_two_axes_are_designed_against_one_common_schedule(opts) -> None:
    """
    **Multi-axis.** Velocity encoding on `z`, flow compensation on `y`, one window, one echo.

    Each axis is its own physical problem -- one `JointProblem` each, never a single
    multidimensional solver -- but the *schedule* is shared, because they share a winder and an
    echo.  Both are then measured on the complete repetition at that one achieved echo, which is
    the check that the common schedule really was common.
    """
    shared = dict(opts=opts, fov_mm=(220.0, 220.0, 120.0), matrix=(16, 16, 8), flip_deg=15.0,
                  slab_thickness_mm=120.0)
    kernel = sc.modules.GRE3DTR(
        **shared, encoding_states=POLARITIES,
        joint_claims=(DifferenceClaim('z', 1, delta(), POLARITIES),
                      CommonModeClaim('z', 1), CommonModeClaim('y', 1)))

    assert sorted(kernel._joint) == ['y', 'z'], 'one design per claimed axis'
    schedules = {designed.schedule for designed in kernel._joint.values()}
    assert len(schedules) == 1, 'and one schedule shared between them'

    origin, echo = kernel.exc.time_to_center(), kernel.time_to_echo()
    for axis, count, wanted_m1 in (('y', kernel.matrix[1], lambda state: 0.0),
                                   ('z', kernel.matrix[2], lambda state: state * delta() / 2.0)):
        for index in range(count):
            for state in POLARITIES:
                shot = kernel(line=index if axis == 'y' else 8,
                              partition=index if axis == 'z' else 4, encoding_state=state)
                m1 = measure_moment(shot, 1, axis, origin_s=origin, start_s=0.0, end_s=echo)
                assert m1 == pytest.approx(wanted_m1(state), abs=1e-11)


def test_the_compiler_accepts_a_jointly_designed_repetition(opts) -> None:
    """The contract is the emitted file, so the repetition has to survive compilation."""
    kernel = sc.modules.GRE2DTR(opts=opts, fov_mm=220.0, matrix=(32, 32), thickness_mm=5.0,
                                flip_deg=15.0, joint_claims=claims('y'),
                                encoding_states=POLARITIES)

    seq = sc.compile(kernel(line=8, encoding_state=+1), opts)

    assert len(seq.block_events) > 3


@pytest.mark.parametrize('extra_ms', (0.0, 0.5, 2.0, 5.0))
def test_an_explicit_te_above_the_minimum_does_not_stale_the_first_moment(opts,
                                                                          extra_ms: float
                                                                          ) -> None:
    """
    **The schedule-shift bug, pinned.**

    An explicit TE inserts fill *in front of* the winder, so a waveform solved at the minimum
    schedule would be translated bodily -- and translating changes the first moment whenever
    ``m0`` is non-zero, by ``dt * m0``.  Before the common-schedule fix this was wrong by
    0.264 s/m at TE + 2 ms against a target of 0.167.

    A non-centre line is the case that shows it; at ``ky = 0`` the error vanishes and the bug
    hides.
    """
    shared = dict(opts=opts, fov_mm=220.0, matrix=(32, 32), thickness_mm=5.0, flip_deg=15.0)
    joint = dict(joint_claims=claims('y'), encoding_states=POLARITIES)
    floor = sc.modules.GRE2DTR(**shared, **joint).min_te_s

    kernel = sc.modules.GRE2DTR(**shared, **joint, te_s=floor + extra_ms * 1e-3)
    origin, echo = kernel.exc.time_to_center(), kernel.time_to_echo()

    assert kernel.te_s == pytest.approx(floor + extra_ms * 1e-3, abs=1e-9)
    for line in (0, 5, 31):
        for state in POLARITIES:
            shot = kernel(line=line, encoding_state=state)
            assert measure_moment(shot, 1, 'y', origin_s=origin, start_s=0.0,
                                  end_s=echo) == pytest.approx(state * delta() / 2.0, abs=1e-11)
            assert measure_moment(shot, 0, 'y', origin_s=origin, start_s=0.0,
                                  end_s=echo) == pytest.approx(kernel.pe.k_per_m(line),
                                                               abs=1e-6 * max(1.0, abs(
                                                                   kernel.pe.k_per_m(line))))


def test_the_design_schedule_is_the_schedule_the_kernel_emits(opts) -> None:
    """
    The invariant behind both schedule bugs, asserted directly rather than via a protocol.

    A design made against one schedule and emitted at another is stale whenever ``m0`` is
    non-zero.  I could not construct a GRE protocol where *axis-local-then-max* bit -- the joint
    axis wanted the longest window in every case tried -- so the property is pinned here instead
    of relying on a case that happens to expose it.
    """
    shared = dict(opts=opts, fov_mm=220.0, matrix=(32, 32), thickness_mm=5.0, flip_deg=15.0)
    for extra in (0.0, 1.5e-3):
        floor = sc.modules.GRE2DTR(**shared, joint_claims=claims('y'),
                                   encoding_states=POLARITIES).min_te_s
        kernel = sc.modules.GRE2DTR(**shared, joint_claims=claims('y'),
                                    encoding_states=POLARITIES, te_s=floor + extra)
        designed = kernel._joint['y'].schedule

        assert designed.window_s == pytest.approx(kernel.winder_s, abs=1e-12)
        assert designed.endpoint_s == pytest.approx(kernel.time_to_echo(), abs=1e-9)
        assert designed.origin_s == pytest.approx(kernel.exc.time_to_center(), abs=1e-12)

