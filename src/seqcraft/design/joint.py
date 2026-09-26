r"""
**Provisional, internal.** Repetition-level joint moment design: requirements, then realisation.

Nothing here is public.  It is the shared part of the architecture in
``plans/current/2026-09-24_repetition_physical_design_architecture.md`` -- the boundary and the
orchestration, **not** one solver that grows to cover every MRI design problem.

The stages, composed as functions rather than as a class hierarchy:

.. code-block:: text

    claims + states        resolve_claims   -> one absolute target per state
    targets + fixed        JointProblem     -> what must be true, and what is already there
    a candidate schedule   Schedule         -> semantic origin, endpoint, and the free window
    a realisation family   realise_*        -> events, or None, plus utilisation
    the cascade            design           -> the first feasible schedule and family

Two separations matter.

**The realisation layer never learns where a target came from.**  It sees a number.  That is what
keeps velocity encoding's state semantics -- a *difference* between encoding states -- out of the
waveform arithmetic.

**Timing is part of the physical design, not a fixed frame around it.**  A longer window moves the
echo, which moves every fixed contribution measured to it, so a candidate schedule is re-evaluated
in full rather than patched.  :class:`JointProblem` therefore takes `fixed` as a **callable** of
state and schedule.

Scope of the moment arithmetic here: **orders 0 and 1 only**.  For a symmetric lobe ``m0`` is its
area and ``m1`` is its area times its centroid, exactly; that shortcut does not extend to ``m2``,
where the lobe's own internal shape contributes.  A second-order family would have to compute its
own exact response and carry enough independent degrees of freedom.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ..errors import ConfigurationError, format_error
from .events import knots_of, pwl_moment
from .logic import flatten

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Mapping, Sequence

    from pypulseq.opts import Opts

    from .events import Event
    from .logic import LogicBlock

__all__ = [
    'ORDERS',
    'CommonModeClaim',
    'DifferenceClaim',
    'JointDesign',
    'JointProblem',
    'Realisation',
    'Schedule',
    'attempt',
    'claimed_axis',
    'design',
    'group_by_axis',
    'placed',
    'refuse_infeasible',
    'refuse_search_exhausted',
    'require_owned_axes',
    'utilisation',
    'events_for',
    'measure_moment',
    'resolve_claims',
]

#: Moments below this are zero for the purpose of choosing a lobe, in the axis's own units.
_NEGLIGIBLE = 1e-9

#: The orders this module's realisation families compute exactly.  See the module docstring.
ORDERS = (0, 1)

#: How many candidate windows a kernel's schedule search will try before giving up.  An
#: implementation guard, not a physical statement: 4000 gradient raster steps is 40 ms at the
#: usual 10 us raster, far past any repetition that would be useful, and reaching it is reported
#: as a search limit rather than as infeasibility.  See :func:`refuse_search_exhausted`.
SEARCH_LIMIT_WINDOWS = 4000


# ------------------------------------------------------------------ the independent validator
def measure_moment(tree: LogicBlock, order: int, axis: str, *,
                   origin_s: float, start_s: float, end_s: float) -> float:
    """
    The `order`-th moment of **every** gradient on `axis` over ``[start_s, end_s]``, about
    `origin_s`.

    Integrated from each event's own piecewise-linear knots and summed across events, which is the
    part no per-event check reaches.  For tests, not for the design: a validator that asked the
    designer for its own residual would be checking arithmetic rather than physics.
    """
    total = 0.0
    for at, event, _ in flatten(tree):
        if getattr(event, 'type', None) not in ('trap', 'grad'):
            continue
        if getattr(event, 'channel', None) != axis:
            continue
        times, amps = knots_of(event, at)
        if times.size < 2:
            continue
        low, high = max(float(times[0]), start_s), min(float(times[-1]), end_s)
        if high <= low:
            continue
        if low > times[0] or high < times[-1]:
            inside = (times > low) & (times < high)
            edges = (float(np.interp(low, times, amps)), float(np.interp(high, times, amps)))
            times = np.concatenate(([low], times[inside], [high]))
            amps = np.concatenate(([edges[0]], amps[inside], [edges[1]]))
        total += pwl_moment(times - origin_s, amps, order)
    return float(total)


# ------------------------------------------------------------------------ the requirement layer
@dataclass(frozen=True)
class DifferenceClaim:
    """``m_order(state_a) - m_order(state_b) = value`` -- what velocity encoding asks for."""

    axis: str
    order: int
    value: float
    states: tuple[object, object]


@dataclass(frozen=True)
class CommonModeClaim:
    """``mean over states of m_order = value`` -- what flow compensation asks for."""

    axis: str
    order: int
    value: float = 0.0


def group_by_axis(claims: Iterable[object]) -> dict[str, tuple[object, ...]]:
    """Claims grouped by the axis they constrain.  One physical problem per axis."""
    grouped: dict[str, list[object]] = {}
    for claim in claims:
        grouped.setdefault(str(getattr(claim, 'axis')), []).append(claim)
    return {axis: tuple(group) for axis, group in grouped.items()}


def placed(event: Event, at_s: float) -> LogicBlock:
    """One event in a block at `at_s` -- for measuring a fixed contribution where it will play."""
    from .logic import LogicBlock as _LogicBlock
    return _LogicBlock().add(at_s, event)


def refuse_infeasible(problem: JointProblem) -> None:
    """Public spelling of the refusal, for a kernel whose candidate schedules all failed."""
    _refuse_infeasible(problem)


def refuse_search_exhausted(problem: JointProblem, *, steps: int, raster_s: float) -> None:
    """
    The schedule search hit its own ceiling, which is **not** a proof that nothing would work.

    A window search has to stop somewhere, and the number it stops at is an implementation guard.
    Reporting that as "these moments cannot be realised" would be a claim the search never
    established -- the next window along was simply never tried.  The two need different words
    because they need different fixes: one is relax the physics, the other is raise the ceiling.
    """
    _refuse_search_exhausted(problem, steps, raster_s)


def require_owned_axes(claims: Iterable[object], owned: Sequence[str], *, component: str) -> None:
    """
    Refuse a claim on an axis this repetition cannot carry -- **before** any realisation.

    The question is physical and it belongs to the repetition, which is the only thing that knows
    which of its axes have an adjustable window between the semantic origin and the endpoint that
    it will actually **materialise**.  Not a class-name test, and not a check that the letter is
    one of ``x``, ``y``, ``z``: a kernel can own a legal axis it never emits on, and designing a
    waveform for that axis would lengthen the winder, push out TE, and then drop the events.

    `owned` is that list, and the caller is responsible for it naming the axes it emits rather
    than the axes it can imagine.
    """
    allowed = tuple(owned)
    for claim in claims:
        axis = str(getattr(claim, 'axis', ''))
        if axis not in allowed:
            _refuse_unowned_axis(axis, allowed, component)


def claimed_axis(claims: Sequence[object]) -> str:
    """The single axis a set of claims refers to, or a refusal."""
    axes = {getattr(claim, 'axis', None) for claim in claims}
    if len(axes) == 1:
        return str(next(iter(axes)))
    msg = format_error(
        f'these claims name {len(axes)} axes, and one coupled design serves one axis.',
        {'axes': sorted(str(a) for a in axes)},
        ['claim one axis per design -- each logical axis is encoded independently'],
    )
    raise ConfigurationError(msg)


def resolve_claims(claims: Iterable[object], states: Sequence[object], *,
                   axis: str, order: int,
                   base: Mapping[object, float]) -> dict[object, float]:
    """
    Turn claims into one **absolute** target per state, or refuse.

    `base` is what the sequence itself wants at this order -- the encode area at order 0, and
    nothing at order 1.  A claim on the difference and a claim on the common mode constrain
    different components of the same state-indexed vector, so they compose without a precedence
    rule; two claims on the *same* component do not, and raise.

    >>> resolve_claims([DifferenceClaim('y', 1, 0.5, (1, -1)), CommonModeClaim('y', 1)],
    ...                [1, -1], axis='y', order=1, base={1: 0.0, -1: 0.0})
    {1: 0.25, -1: -0.25}
    """
    mine = [c for c in claims
            if getattr(c, 'axis', None) == axis and getattr(c, 'order', None) == order]
    difference = [c for c in mine if isinstance(c, DifferenceClaim)]
    common = [c for c in mine if isinstance(c, CommonModeClaim)]
    for kind, found in (('difference', difference), ('common mode', common)):
        if len(found) > 1:
            _refuse_duplicate(kind, axis, order, len(found))

    out = {state: float(base[state]) for state in states}
    if common:
        shift = common[0].value - float(np.mean([out[s] for s in states]))
        out = {state: value + shift for state, value in out.items()}
    if difference:
        claim = difference[0]
        first, second = claim.states
        if first not in out or second not in out:
            _refuse_missing_state(claim, states)
        half = (claim.value - (out[first] - out[second])) / 2.0
        out[first] += half
        out[second] -= half
    return out


# ---------------------------------------------------------------------- the physical problem
@dataclass(frozen=True)
class Schedule:
    """
    One candidate repetition schedule for one axis.

    The window and the endpoint move together, which is the point: lengthening a winder pushes
    the echo out, and every fixed contribution is measured *to* the echo.
    """

    origin_s: float
    endpoint_s: float
    window_start_s: float
    window_s: float


@dataclass(frozen=True)
class JointProblem:
    """
    What must be true on one axis, and what is already there.

    `fixed` is a callable rather than a tuple because the contributions depend on both the state
    and the schedule -- a wave gradient's moment about the echo changes when the echo moves.  It
    returns ``(m0, m1)`` already accumulated over the schedule's interval.
    """

    axis: str
    targets: Mapping[object, tuple[float, float | None]]
    fixed: Callable[[object, Schedule], tuple[float, float]]


@dataclass(frozen=True)
class Realisation:
    """One state's waveform at one schedule, with why it did or did not fit."""

    events: tuple[Event, ...]
    family: str
    peak_grad: float
    peak_slew: float
    limiting: str = ''

    @property
    def feasible(self) -> bool:
        return not self.limiting


@dataclass(frozen=True)
class JointDesign:
    """The chosen schedule and family, and one realisation per state."""

    schedule: Schedule
    family: str
    realisations: Mapping[object, Realisation]

    @property
    def peak_grad(self) -> float:
        return max(r.peak_grad for r in self.realisations.values())

    @property
    def peak_slew(self) -> float:
        return max(r.peak_slew for r in self.realisations.values())


# ------------------------------------------------------------------- the realisation families
def _lobe(axis: str, area: float, duration_s: float, opts: Opts) -> Event | None:
    """One trapezoid, or ``None`` when the area is negligible or it will not fit."""
    if abs(area) <= _NEGLIGIBLE:
        return None
    try:
        return pp.make_trapezoid(channel=axis, area=area, duration=duration_s, system=opts)
    except (ValueError, ZeroDivisionError):
        return _UNFITTABLE


#: Sentinel: pypulseq refused to build this lobe at this duration.
_UNFITTABLE: Event = object()  # type: ignore[assignment]


def _assemble(axis: str, pieces: Sequence[tuple[float, float]], opts: Opts,
              family: str) -> Realisation | None:
    """
    Build ``(area, duration)`` pieces in order and report utilisation, or ``None`` if impossible.

    Utilisation is reported rather than only feasibility, because a family that fails everywhere
    over a wide interval and one that fails by two per cent are different problems -- section 20
    of the architecture document.
    """
    events: list[Event] = []
    peak_grad = peak_slew = 0.0
    limiting = ''
    for area, duration_s in pieces:
        lobe = _lobe(axis, area, duration_s, opts)
        if lobe is _UNFITTABLE:
            return None
        if lobe is None:
            continue
        amplitude = abs(float(lobe.amplitude))
        peak_grad = max(peak_grad, amplitude / float(opts.max_grad))
        peak_slew = max(peak_slew, amplitude / float(lobe.rise_time) / float(opts.max_slew))
        events.append(lobe)
    if peak_grad > 1.0 + 1e-9:
        limiting = 'gradient'
    elif peak_slew > 1.0 + 1e-9:
        limiting = 'slew'
    return Realisation(tuple(events), family, peak_grad, peak_slew, limiting)


def realise_single_lobe(axis: str, target: tuple[float, float | None],
                        fixed: tuple[float, float], schedule: Schedule,
                        opts: Opts) -> Realisation | None:
    """
    One lobe filling the window -- the shape every kernel already uses, written as a family.

    **It serves an m0 target only.**  One degree of freedom cannot satisfy two conditions, so it
    declines unless the first moment is either unconstrained (``None``) or happens to equal the
    ``area * centroid`` this lobe would produce anyway.  Checking what it *would* produce rather
    than only whether a target was supplied is the difference between declining and being quietly
    wrong: a target of zero is not the same as no target.
    """
    area = target[0] - fixed[0]
    if target[1] is not None:
        centre = schedule.window_start_s + schedule.window_s / 2.0 - schedule.origin_s
        if abs(fixed[1] + area * centre - target[1]) > _NEGLIGIBLE:
            return None
    return _assemble(axis, [(area, schedule.window_s)], opts, 'single-lobe')


def _two_lobe_areas(target: tuple[float, float], fixed: tuple[float, float],
                    centres: tuple[float, float]) -> tuple[float, float]:
    """Solve ``[[1, 1], [c0, c1]] a = target - fixed`` -- exact for symmetric lobes at m0/m1."""
    wanted = np.asarray(target, dtype=float) - np.asarray(fixed, dtype=float)
    matrix = np.array([[1.0, 1.0], list(centres)], dtype=float)
    return tuple(float(v) for v in np.linalg.solve(matrix, wanted))  # type: ignore[return-value]


def realise_two_lobes(axis: str, target: tuple[float, float], fixed: tuple[float, float],
                      schedule: Schedule, opts: Opts, *,
                      split: float = 0.5) -> Realisation | None:
    """
    Two adjacent lobes filling the window, the first taking `split` of it.

    ``split = 0.5`` is PR #39's equal-duration shape.  Letting the split move is one extra degree
    of freedom and costs nothing but a search over a raster-aligned fraction -- which is what
    :func:`realise_split_search` does with it.

    **The boundary between the two lobes is snapped onto the gradient raster**, because it is an
    instant the sequence has to be able to name.  Halving an odd number of raster steps does not
    land on one: a 1450 us window split evenly asks for two 725 us lobes, and the compiler
    rightly refuses the second for starting off-raster.  The snap moves the boundary by less than
    one step and the areas are solved against the durations that result, so the moments stay
    exact rather than being corrected afterwards.
    """
    if target[1] is None:
        return None       # under-determined: nothing constrains the second degree of freedom
    raster = float(opts.grad_raster_time)
    steps = max(int(round(schedule.window_s * split / raster)), 1)
    first_s = steps * raster
    second_s = schedule.window_s - first_s
    if min(first_s, second_s) <= 0.0:
        return None
    centres = (schedule.window_start_s + first_s / 2.0 - schedule.origin_s,
               schedule.window_start_s + first_s + second_s / 2.0 - schedule.origin_s)
    areas = _two_lobe_areas(target, fixed, centres)
    return _assemble(axis, [(areas[0], first_s), (areas[1], second_s)], opts,
                     'two-lobe' if split == 0.5 else f'two-lobe/{split:.3f}')


def realise_base_plus_bipolar(axis: str, target: tuple[float, float],
                              fixed: tuple[float, float], schedule: Schedule,
                              opts: Opts) -> Realisation | None:
    """
    A base lobe carrying the area, plus a **zero-area** bipolar carrying pure first moment.

    The decoupled basis `wave-gre-flow-comp` uses.  Same span as :func:`realise_two_lobes` and the
    same two degrees of freedom, but the columns are orthogonal in what they do, which changes
    which target values push a lobe over the limit first.
    """
    if target[1] is None:
        return None       # under-determined: nothing constrains the second degree of freedom
    half = schedule.window_s / 2.0
    if half <= 0.0:
        return None
    area = target[0] - fixed[0]
    base_centre = schedule.window_start_s + schedule.window_s / 2.0 - schedule.origin_s
    # A zero-area bipolar of lobes +-b over the window has m1 = -b * half.
    residual = target[1] - fixed[1] - area * base_centre
    bipolar = -residual / half
    return _assemble(axis, [(area / 2.0 + bipolar, half), (area / 2.0 - bipolar, half)],
                     opts, 'base+bipolar')


def realise_split_search(axis: str, target: tuple[float, float], fixed: tuple[float, float],
                         schedule: Schedule, opts: Opts) -> Realisation | None:
    """
    The two-lobe family with its split searched on the gradient raster.

    A small, low-dimensional shape family -- the same kind of freedom `wave-gre-flow-comp` searches
    inside a fixed duration.  It returns the feasible split with the most headroom, so that a
    schedule which only just fits is not chosen over one that fits comfortably.
    """
    raster = float(opts.grad_raster_time)
    steps = int(round(schedule.window_s / raster))
    best: Realisation | None = None
    for first in range(1, steps):
        found = realise_two_lobes(axis, target, fixed, schedule, opts, split=first / steps)
        if found is None or not found.feasible:
            continue
        if best is None or max(found.peak_grad, found.peak_slew) < max(best.peak_grad,
                                                                      best.peak_slew):
            best = found
    return best


#: The cascade: cheapest first, then richer.  Order is the policy; each entry owns its own maths.
FAMILIES: tuple[Callable[..., Realisation | None], ...] = (
    realise_single_lobe,
    realise_two_lobes,
    realise_base_plus_bipolar,
    realise_split_search,
)


# ------------------------------------------------------------------------- the orchestration
def attempt(problem: JointProblem, schedule: Schedule, opts: Opts, *,
            families: Sequence[Callable[..., Realisation | None]] = FAMILIES,
            ) -> JointDesign | None:
    """
    Try each family at one candidate schedule, and return the first that serves **every** state.

    Every state, because a window that served only some would put TE on a k-space axis.
    """
    for family in families:
        realisations = {}
        for state, target in problem.targets.items():
            found = family(problem.axis, target, problem.fixed(state, schedule), schedule, opts)
            if found is None or not found.feasible:
                realisations = {}
                break
            realisations[state] = found
        if realisations:
            name = next(iter(realisations.values())).family
            return JointDesign(schedule, name, realisations)
    return None


def design(problem: JointProblem, schedules: Iterable[Schedule], opts: Opts, *,
           families: Sequence[Callable[..., Realisation | None]] = FAMILIES,
           ) -> JointDesign:
    """
    Walk candidate schedules in order and return the first that any family can serve.

    **The schedule search is the caller's**, supplied as an iterable, because what may move -- the
    window, TE, the echo spacing -- and in what order is repetition timing policy, which belongs
    to the kernel.  This function owns only "try them in the order given".

    No duration-search algorithm is baked in here.  A family may bisect where it has a
    monotonicity proof, enumerate a shape, or scan; the generic layer does not assume which.
    """
    for schedule in schedules:
        found = attempt(problem, schedule, opts, families=families)
        if found is not None:
            return found
    _refuse_infeasible(problem)


def events_for(designed: JointDesign, state: object) -> tuple[tuple[float, Event], ...]:
    """The state's events with the instants they start at, ready to place in a block."""
    at = designed.schedule.window_start_s
    placed = []
    for event in designed.realisations[state].events:
        placed.append((at, event))
        at += float(pp.calc_duration(event))
    return tuple(placed)


# ----------------------------------------------------------------------- the diagnostics
def _lobe_utilisation(area: float, duration_s: float, opts: Opts, *,
                      steps: int = 64) -> tuple[float, float]:
    """
    The best ``(G / G_limit, S / S_limit)`` a symmetric trapezoid of this area and duration can do.

    Parameterised by the ramp time `r`, which is the shape's only freedom once the area and the
    duration are fixed: ``g = A / (T - r)`` and ``slew = g / r``.  A short ramp lowers the peak
    amplitude and raises the slew; a long one does the reverse.  Scanning `r` and keeping the best
    ``max`` of the two fractions is what turns "pypulseq refused" into "it needed 5.5x the slew" --
    which is the difference between a family that nearly fits and one that is nowhere near.

    Not a design path: nothing here is emitted.
    """
    if abs(area) <= _NEGLIGIBLE:
        return 0.0, 0.0
    best = (float('inf'), float('inf'))
    for step in range(1, steps + 1):
        rise_s = duration_s * 0.5 * step / steps
        flat_s = duration_s - rise_s
        if flat_s <= 0.0:
            continue
        amplitude = abs(area) / flat_s
        fractions = (amplitude / float(opts.max_grad),
                     amplitude / rise_s / float(opts.max_slew))
        if max(fractions) < max(best):
            best = fractions
    return best


def utilisation(axis: str, target: tuple[float, float | None], fixed: tuple[float, float],
                schedule: Schedule, opts: Opts, *, splits: int | None = None) -> dict[str, float]:
    """
    How close the best shape comes at this schedule, **even when it does not fit**.

    ``pp.make_trapezoid`` refuses an over-limit design, so a feasibility search reads as
    ``None, None, None, success`` and cannot say whether a candidate missed by two per cent or by
    a factor of ten.  This reports

    .. code-block:: text

        u(T) = min over allowed shapes  max(G_required / G_limit, S_required / S_limit)

    over the two-lobe family's split and each lobe's ramp.  Diagnostic only; the emitter never
    takes this path, so no illegal waveform can reach a sequence through it.

    **The split is searched on the gradient raster**, which is the lattice
    :func:`realise_split_search` can actually emit on.  It used to be a fixed 48 fractions of the
    window, and that is the wrong lattice in both directions: it proposes splits no raster can
    hold, and it misses the ones that fit.  ``u(T)`` is not smooth -- when ``m1 / m0`` happens to
    equal a lobe centre the two-lobe solve degenerates to a single comfortable lobe, so the curve
    has narrow minima -- and a coarse fixed grid steps over them.  Measured against the raster
    lattice on a 145 1/m, 0.055 s/m target at 20 mT/m and 45 T/m/s, the old grid overstated ``u``
    by up to **x1.82**, reporting 1.171 where the lattice reaches 0.653.  That is a diagnostic
    that manufactures a cliff, which is the one thing this function exists to rule out.
    """
    best = {'utilisation': float('inf'), 'grad': float('inf'), 'slew': float('inf'),
            'split': float('nan')}
    if target[1] is None:
        return best
    steps = splits if splits is not None else int(round(schedule.window_s /
                                                        float(opts.grad_raster_time)))
    for step in range(1, max(steps, 2)):
        split = step / max(steps, 2)
        first_s = schedule.window_s * split
        second_s = schedule.window_s - first_s
        if min(first_s, second_s) <= 0.0:
            continue
        centres = (schedule.window_start_s + first_s / 2.0 - schedule.origin_s,
                   schedule.window_start_s + first_s + second_s / 2.0 - schedule.origin_s)
        areas = _two_lobe_areas((target[0], target[1]), fixed, centres)
        grad = slew = 0.0
        for area, duration_s in zip(areas, (first_s, second_s)):
            lobe_grad, lobe_slew = _lobe_utilisation(area, duration_s, opts)
            grad, slew = max(grad, lobe_grad), max(slew, lobe_slew)
        if max(grad, slew) < best['utilisation']:
            best = {'utilisation': max(grad, slew), 'grad': grad, 'slew': slew, 'split': split}
    return best


# -------------------------------------------------------------------------------- the refusals
def _refuse_duplicate(component: str, axis: str, order: int, count: int) -> None:
    msg = format_error(
        f'{count} augmentations claim the {component} of moment {order} on {axis}.',
        {'axis': axis, 'order': order, 'component': component, 'claims': count},
        ['one augmentation owns one component of one moment on one axis',
         'a second claim is a disagreement about physics, not something to order by precedence'],
    )
    raise ConfigurationError(msg)


def _refuse_search_exhausted(problem: JointProblem, steps: int, raster_s: float) -> None:
    msg = format_error(
        f'the schedule search reached its limit of {steps} windows on {problem.axis} without '
        f'finding one that fits.  This is a design search limit, not a proof of infeasibility: '
        f'windows longer than {steps * raster_s * 1e3:.1f} ms were never tried.',
        {'axis': problem.axis, 'states': len(problem.targets),
         'search_limit_windows': steps, 'longest_window_tried_s': steps * raster_s},
        [f'the requirement may still be realisable beyond {steps * raster_s * 1e3:.1f} ms; '
         'if that is plausible here, raise the search ceiling',
         'more usually it is not -- a target needing a window this long is worth rechecking '
         'against the units it was stated in',
         'relaxing the target or widening the gradient limits shortens the window it needs'],
    )
    raise ConfigurationError(msg)


def _refuse_unowned_axis(axis: str, owned: Sequence[str], component: str) -> None:
    msg = format_error(
        f'{component} does not own an adjustable window on {axis!r}, so it cannot carry a '
        f'moment claim there.',
        {'axis': axis, 'owned_axes': tuple(owned), 'component': component},
        [f'claim one of {tuple(owned)}, which is what this repetition designs and emits',
         'an axis whose gradients this repetition does not own is realised by the leaf that '
         'does own them, and cannot be jointly designed here',
         'this is about what gets emitted, not about whether the axis name is legal'],
    )
    raise ConfigurationError(msg)


def _refuse_missing_state(claim: DifferenceClaim, states: Sequence[object]) -> None:
    msg = format_error(
        f'a difference claim on moment {claim.order} of {claim.axis} needs states '
        f'{claim.states}, and this acquisition has {tuple(states)}.',
        {'axis': claim.axis, 'order': claim.order,
         'claim_states': claim.states, 'acquired_states': tuple(states)},
        ['acquire both states -- a phase difference needs two acquisitions to subtract',
         'or drop the claim, which leaves this moment to whatever else constrains it'],
    )
    raise ConfigurationError(msg)


def _refuse_infeasible(problem: JointProblem) -> None:
    msg = format_error(
        f'no candidate schedule realises these moments on {problem.axis}.',
        {'axis': problem.axis, 'states': len(problem.targets)},
        ['allow a longer window, which is what a longer TE buys',
         'or relax the target -- a larger venc needs less first moment',
         'or widen the gradient limits if the scanner allows it'],
    )
    raise ConfigurationError(msg)
