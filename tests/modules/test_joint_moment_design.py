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
    require_owned_axes,
    resolve_claims,
    utilisation,
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


def test_the_diagnostic_agrees_with_what_the_emitter_can_actually_do(opts) -> None:
    """
    `u <= 1` has to mean the same thing as "a family fitted", or it is not a diagnostic.

    It reports the best the two-lobe family can do, and that family emits on the gradient
    raster -- so the split has to be searched on the raster too.  A coarser fixed grid is wrong
    in both directions: it proposes splits nothing can emit, and it steps over the narrow minima
    where `m1 / m0` lands on a lobe centre and the solve degenerates to one comfortable lobe.
    A stress sweep measured the old fixed grid overstating `u` by up to x1.82 -- a diagnostic
    reporting a cliff that the emitter walks straight past.

    The residual disagreement is the boundary itself: `_lobe_utilisation` scans ramp times
    analytically while `make_trapezoid` picks one shape, so the two differ by a fraction of a
    per cent at `u == 1`.  Anything further apart than that is a real divergence.
    """
    target = (31 * 1e3 / 220.0, +delta() / 2.0)
    disagreements = []
    for steps in range(30, 140, 3):
        schedule = a_schedule(steps * 10e-6)
        fitted = realise_split_search('y', target, (0.0, 0.0), schedule, opts) is not None
        u = utilisation('y', target, (0.0, 0.0), schedule, opts)['utilisation']
        if fitted != (u <= 1.0):
            disagreements.append((steps * 10e-6, u))

    assert all(abs(u - 1.0) < 0.05 for _, u in disagreements), (
        f'the diagnostic and the emitter disagree away from the boundary: {disagreements}'
    )


def test_the_coarse_grid_the_diagnostic_used_to_sample_gets_the_boundary_wrong(opts) -> None:
    """
    The concrete window where the old fixed grid disagreed with the emitter, kept as the reason.

    At a 40 1/m, 0.197 s/m target the split search fits a 930 us window.  Sampling 48 fixed
    fractions reports `u = 1.014` there -- infeasible -- while the raster lattice reports 0.994.
    One window is all it takes: a schedule search that trusts `u` walks past a window the
    emitter would have accepted, and the next one it takes is longer TE bought for nothing.
    Sweeping this schedule found 148 such windows.

    The gap widens when the endpoint moves with the window, as it does in a real repetition,
    because then the lobe centres sweep too and the narrow minima move between samples; there
    the old grid was measured overstating `u` by up to x1.82.
    """
    target = (40.0, 0.1971)
    schedule = a_schedule(930e-6)

    assert realise_split_search('y', target, (0.0, 0.0), schedule, opts) is not None
    assert utilisation('y', target, (0.0, 0.0), schedule, opts)['utilisation'] <= 1.0
    assert utilisation('y', target, (0.0, 0.0), schedule, opts,
                       splits=48)['utilisation'] > 1.0


def test_a_family_reports_utilisation_and_what_limits_it(opts) -> None:
    """A bool cannot distinguish "failed by two per cent" from "failed everywhere"."""
    target = (31 * 1e3 / 220.0, +delta() / 2.0)

    tight = realise_two_lobes('y', target, (0.0, 0.0), a_schedule(300e-6), opts)

    assert tight is None or tight.limiting in ('gradient', 'slew')
    if tight is not None:
        assert tight.peak_grad > 0.0 and tight.peak_slew > 0.0


# ----------------------------------------------------------------------- capability
def _emits_on(shot, axis: str) -> bool:
    """Does the emitted tree actually play a gradient on `axis`, from the joint block?"""
    return any(getattr(event, 'channel', None) == axis and 'joint' in path
               for _, event, path in sc.design.logic.flatten(shot))


@pytest.mark.parametrize('kernel_of, build_kwargs', [
    (lambda opts, ax: sc.modules.GRE2DTR(
        opts=opts, fov_mm=220.0, matrix=(32, 32), thickness_mm=5.0,
        joint_claims=(CommonModeClaim(ax, 1),), encoding_states=('only',)), {'line': 4}),
    (lambda opts, ax: sc.modules.GRE3DTR(
        opts=opts, fov_mm=(220.0, 220.0, 120.0), matrix=(32, 32, 8),
        joint_claims=(CommonModeClaim(ax, 1),), encoding_states=('only',)),
     {'line': 4, 'partition': 1}),
])
def test_every_advertised_axis_is_designed_emitted_and_met(opts, kernel_of, build_kwargs) -> None:
    """
    An advertised capability is three claims, and the middle one is the easy one to lose.

    A design that is produced but never materialised is worse than a refusal: the winder grows,
    TE goes out with it, and the events are dropped -- the caller pays for a compensation they do
    not get.  So this asserts the design exists, that something is emitted on the axis, and that
    the **whole** repetition meets the target, for every axis the kernel says it owns.

    Off-centre indices, deliberately: at the centre line the encode area is zero, the design has
    nothing to realise, and the emitted block would be the ordinary blip either way -- which
    would pass this test without exercising the joint path at all.
    """
    probe = kernel_of(opts, 'y')
    for axis in probe.joint_axes:
        kernel = kernel_of(opts, axis)

        assert axis in kernel._joint, f'{axis} is advertised but no design was produced'

        shot = kernel(encoding_state='only', **build_kwargs)
        assert _emits_on(shot, axis), f'{axis} is advertised but nothing is emitted on it'

        off_m1 = measure_moment(shot, 1, axis, origin_s=kernel.exc.time_to_center(),
                                start_s=0.0, end_s=kernel.time_to_echo())
        assert abs(off_m1) < 1e-11, f'{axis} is advertised but the repetition emits m1 = {off_m1}'


@pytest.mark.parametrize('kernel_of, unowned', [
    (lambda opts, ax: sc.modules.GRE2DTR(
        opts=opts, fov_mm=220.0, matrix=(32, 32), thickness_mm=5.0,
        joint_claims=(CommonModeClaim(ax, 1),), encoding_states=('only',)), ('x', 'z', 'q')),
    (lambda opts, ax: sc.modules.GRE3DTR(
        opts=opts, fov_mm=(220.0, 220.0, 120.0), matrix=(32, 32, 8),
        joint_claims=(CommonModeClaim(ax, 1),), encoding_states=('only',)), ('x', 'q')),
])
def test_an_axis_the_repetition_does_not_own_is_refused_before_realisation(
        opts, kernel_of, unowned) -> None:
    """
    `x` and `z` are legal axes these kernels play gradients on, and neither owns their winder.

    Before the check existed, `GRE2DTR` accepted a claim on `z`, designed it, grew the winder
    from 520 to 2400 us -- 1.9 ms of TE the caller paid for -- and then emitted nothing, leaving
    `m1 = -0.736 s/m` on the axis it reported as compensated.  So the refusal has to come from
    what the repetition **materialises**, not from whether the axis letter is legal.
    """
    for axis in unowned:
        with pytest.raises(sc.errors.ConfigurationError) as raised:
            kernel_of(opts, axis)
        assert axis in str(raised.value)
        assert 'adjustable window' in str(raised.value)


def test_running_out_of_candidate_windows_is_not_called_infeasible(opts, monkeypatch) -> None:
    """
    The window search stops at a ceiling, and the ceiling is an implementation number.

    Reporting "these moments cannot be realised" when the search simply stopped would assert
    something it never established -- the next window along was never tried.  The two failures
    need different words because they need different fixes: relax the physics, or raise the
    ceiling.

    The ceiling is lowered here for speed, not to make the refusal fire: the split search tries
    every raster-aligned split at every candidate window, so walking to the real limit of 4000 is
    quadratic and takes the better part of a minute.  That cost is the search's, not this test's.
    """
    monkeypatch.setattr(sc.modules._joint, 'SEARCH_LIMIT_WINDOWS', 200)
    absurd = DifferenceClaim('y', 1, 5000.0, POLARITIES)

    with pytest.raises(sc.errors.ConfigurationError) as raised:
        sc.modules.GRE2DTR(opts=opts, fov_mm=220.0, matrix=(8, 8), thickness_mm=5.0,
                           joint_claims=(absurd,), encoding_states=POLARITIES)

    said = str(raised.value)
    assert 'design search limit' in said or 'search reached its limit' in said
    assert 'not a proof of infeasibility' in said
    assert 'no candidate schedule realises' not in said


def test_capability_is_a_property_of_the_repetition_not_of_the_claim(opts) -> None:
    """
    The same claim is fine on one repetition and refused on another, and nothing inspects a type.

    `z` is owned by `GRE3DTR`, whose z winder carries the partition encode, and not by `GRE2DTR`,
    whose z gradient is the slice rephaser that `Excitation` owns.  A check keyed on the claim --
    its class, or the augmentation that made it -- could not tell those apart.
    """
    claim = CommonModeClaim('z', 1)

    require_owned_axes([claim], ('y', 'z'), component='GRE3DTR')
    with pytest.raises(sc.errors.ConfigurationError):
        require_owned_axes([claim], ('y',), component='GRE2DTR')


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

