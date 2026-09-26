r"""
**Internal.**  The region of a sequence whose physical degrees of freedom a designer may own.

A physical-design scope is **opt-in**.  Building events, modules or functions into a
``LogicBlock`` and compiling it stays valid with no scope anywhere, and not every useful
composition has to be promoted into a packaged kernel to get physical design.  A scope says one
thing:

    *for this part of my sequence, design the physics and the timing together before anything is
    materialised.*

It is not a Module, not a layer of the sequence hierarchy, not necessarily one literal TR, and
not something the compiler knows about.  It is also **not** a ``LogicBlock`` that already exists
and is then patched: the adapter materialises once, from a finished design.

.. code-block:: text

    packaged kernel ---adapter--.
                                 >--- ScopeGeometry + AxisRequirement --> design_scope --> events
    user composition ---adapter-'                                                            |
                                                                                             v
                                                                        LogicBlock --> compiler

The middle never learns which side it came from.

**What the adapter owns, and what this owns.**  An adapter knows its own physics: where the
semantic instants are, what already plays, which axes it may reshape, and which states are worth
designing against.  This module owns the search, the one-schedule-across-axes rule, and the
verification that turns an adapter's *guess* about states into a *correct* schedule.

**Two representational limits**, current rather than fundamental, and both refused rather than
worked around:

.. code-block:: text

    one adjustable window        a family needing two independently adjustable regions -- a
                                 prephaser and a rewinder designed together, say -- does not fit,
                                 because two scopes cannot share one schedule

    a constant tail              what plays between the window and the echo may not depend on the
                                 state.  A rotated spiral arm is fine; an arm whose length varies
                                 per state is not

Also current: ``m0`` and ``m1`` only, which is :mod:`seqcraft.design.joint`'s scope.

**What a scope does not see.**  It guarantees physics from the contributions it is given.  If a
preparation, a gradient history or an RF pathway *outside* the scope affects the target moment,
the two principled answers are to enlarge the scope or to pass that contribution in explicitly as
part of `fixed`.  Nothing here inspects the surrounding sequence, and there is no pathway model.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

from . import joint

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping, Sequence

    from pypulseq.opts import Opts

__all__ = [
    'AxisRequirement',
    'ScopeDesign',
    'ScopeGeometry',
    'design_scope',
    'realise_scope_at',
]


@dataclass(frozen=True)
class ScopeGeometry:
    """
    Where the adjustable region sits between two semantic instants.

    Three numbers, all measured from the start of whatever the adapter will emit:

    .. code-block:: text

        0                origin_s        window_start_s      +window_s      +tail_s
        |                   |                  |                 |             |
        [--- excitation ----X------ fill ------[==== window ====]---- lead ----X
                            ^                                                  ^
                     semantic origin                                   semantic endpoint
                     (moments measured                                 (the echo)
                      from here)

    `tail_s` is what plays between the end of the window and the echo without moving when the
    window grows -- a Cartesian readout's pre-echo lobe, or a spiral arm.  The endpoint *does*
    move with the window, which is why a candidate schedule is re-evaluated in full rather than
    patched: ``m1' = m1 - dt * m0``.

    **The invariant an adapter must hold.**  Semantic instants, fixed contributions and the
    adjustable region must all describe *the same physical decomposition*.  If a scope takes over
    a component's prephaser, winder or rephaser, then timing metadata that component reports for
    its **own** configuration -- one that still includes the region now owned elsewhere -- cannot
    be reused unchanged.  `SpiralReadout.time_to_echo` is measured from a block containing its own
    prephaser; an adapter that sets ``prephase=False`` and designs that prephaser itself must
    subtract `prephaser_duration_s` from the tail.  Not doing so moves the endpoint past the arm
    and into the rewinder, where `fixed` and the emitted block disagree about what is playing --
    which is a residual, not an exception.
    """

    origin_s: float
    window_start_s: float
    tail_s: float

    def schedule_at(self, window_s: float, fill_s: float = 0.0) -> joint.Schedule:
        """The schedule this geometry gives for one window length, plus any explicit-TE fill."""
        start = self.window_start_s + max(fill_s, 0.0)
        return joint.Schedule(
            origin_s=self.origin_s,
            endpoint_s=start + window_s + self.tail_s,
            window_start_s=start,
            window_s=window_s,
        )

    def echo_at(self, window_s: float, fill_s: float = 0.0) -> float:
        """The echo time this geometry reaches, measured from the semantic origin."""
        sched = self.schedule_at(window_s, fill_s)
        return sched.endpoint_s - sched.origin_s


@dataclass(frozen=True)
class AxisRequirement:
    """
    What one axis of one state family needs, and which states are worth designing against.

    `design_states` is a **proposal, not a proof**.  It is the subset the adapter thinks will size
    the schedule -- signed extremes of a phase encode, the axis-aligned angles of a rotated
    spiral, corners of a multidimensional space, a representative handful, or everything when
    nothing better is available.  Correctness does not depend on it being right:

    .. code-block:: text

        design_states  ->  propose a candidate schedule
                       ->  verify it across the supported family
                       ->  widen if necessary

    A good proposal saves search time.  A bad one costs search time.  Neither can ship an
    infeasible state, because :func:`design_scope` realises every state before accepting a
    schedule -- and it has to, since the feasible set has small raster holes, so "this state needs
    the most area" does not imply "the window that serves it serves everything".

    `states` is what will be realised, and nothing here or in :mod:`seqcraft.design.joint` ever
    looks inside one.  That is what lets a state be a phase-encode index, a signed encoding state,
    a spiral interleaf angle, or a tuple of them.
    """

    axis: str
    states: Sequence[object]
    design_states: Sequence[object]
    target: Callable[[object], tuple[float, float | None]]
    fixed: Callable[[object, joint.Schedule], tuple[float, float]]

    def problem_for(self, states: Sequence[object]) -> joint.JointProblem:
        """The `joint` problem restricted to some states -- the proposal, or all of them."""
        return joint.JointProblem(
            axis=self.axis,
            targets={state: self.target(state) for state in states},
            fixed=self.fixed,
        )


@dataclass(frozen=True)
class ScopeDesign:
    """One schedule, and one realisation per state per axis."""

    window_s: float
    fill_s: float
    schedule: joint.Schedule
    designs: Mapping[str, joint.JointDesign]

    def events_for(self, axis: str, state: object) -> tuple[tuple[float, object], ...]:
        """
        The placed events one state needs on one axis.

        Instants are in the **schedule's own frame**, so an adapter whose emitted block starts at
        the same zero can place them directly.  One whose block starts elsewhere subtracts
        ``schedule.window_start_s`` and adds its own offset.
        """
        return joint.events_for(self.designs[axis], state)


def design_scope(geometry: ScopeGeometry, requirements: Sequence[AxisRequirement], opts: Opts, *,
                 min_window_s: float = 0.0, te_request_s: float | None = None,
                 admissible: Callable[[ScopeDesign], bool] | None = None) -> ScopeDesign:
    """
    Size one schedule from the proposed states, then realise and verify the whole family at it.

    Two rules are kept from the kernel loop this generalises, because both of them were bugs
    once: a candidate is a window **and** the fill an explicit TE needs at that window, so an
    explicit TE cannot translate a solved waveform out from under its own moment; and every axis
    is solved against the same schedule, so one axis' minimum cannot widen the window after
    another axis has already solved at it.

    Parameters
    ----------
    admissible
        **The seam a local PNS criterion will attach to, and nothing uses it yet.**  It is given
        a fully realised candidate -- the schedule *and* the waveforms every state would emit --
        because admissibility depends on the realised gradient, not on the timing alone.  Two
        realisation families at the same schedule may differ in whether they pass, and rejecting
        one candidate leaves longer schedules free to be tried, which is why this filters
        candidates rather than constraining targets.

        It is deliberately not another physical intent.  `FlowCompensation` and `VelocityEncoding`
        say *what relation must hold*; this says *which realisations are acceptable*.  Whole-
        sequence PNS stays a later, global evaluator with temporal history across repetitions,
        whose feedback returns through design rather than stretching emitted events.
    """
    raster = float(opts.grad_raster_time)
    proposed = [r.problem_for(r.design_states) for r in requirements]

    # Every echo this geometry can reach sits an exact number of raster steps above the one a
    # zero-length window reaches, so an explicit request is quantised **up** onto that grid once,
    # here, rather than rounded inside the loop.  Two things follow: the fill is then an exact
    # multiple of the raster at every candidate, and the echo time handed back is never earlier
    # than the one asked for.
    wanted_s = te_request_s
    if wanted_s is not None:
        floor = geometry.echo_at(0.0)
        wanted_s = floor + max(math.ceil((wanted_s - floor) / raster - 1e-9), 0) * raster

    exhausted = True
    for steps in range(max(int(round(min_window_s / raster)), 1), joint.SEARCH_LIMIT_WINDOWS):
        window = steps * raster
        fill = 0.0
        if wanted_s is not None:
            fill = wanted_s - geometry.echo_at(window)
            if fill < -1e-12:
                # Already past the requested echo; a longer window only overshoots further.
                exhausted = False
                break
            fill = max(fill, 0.0)
        schedule = geometry.schedule_at(window, fill)
        if not all(joint.attempt(p, schedule, opts) is not None for p in proposed):
            continue
        designs = _realise_every_state(requirements, schedule, opts)
        if designs is None:
            continue
        candidate = ScopeDesign(window, max(fill, 0.0), schedule, designs)
        if admissible is None or admissible(candidate):
            return candidate
    first = proposed[0]
    if exhausted:
        return joint.refuse_search_exhausted(
            first, steps=joint.SEARCH_LIMIT_WINDOWS, raster_s=raster)
    return joint.refuse_infeasible(first)


def realise_scope_at(geometry: ScopeGeometry, requirements: Sequence[AxisRequirement], opts: Opts,
                     *, window_s: float, fill_s: float = 0.0) -> ScopeDesign | None:
    """
    Realise the whole family again at a schedule someone else chose, or ``None`` if it will not fit.

    This is what keeps a designed family **re-realisable**.  A protocol holding two families with
    different minima takes the longer echo time and asks both to design again at it, with no
    acquisition-wide timing manager anywhere.

    Going back through realisation is the point.  Stretching the waveform that already exists
    would be a different operation with a different answer: every moment is measured to the echo,
    and the echo moved.
    """
    schedule = geometry.schedule_at(window_s, fill_s)
    designs = _realise_every_state(requirements, schedule, opts)
    return None if designs is None else ScopeDesign(window_s, fill_s, schedule, designs)


def _realise_every_state(requirements: Sequence[AxisRequirement], schedule: joint.Schedule,
                         opts: Opts) -> dict[str, joint.JointDesign] | None:
    """Every state of every axis at one fixed schedule, or ``None`` when one of them will not fit."""
    out: dict[str, joint.JointDesign] = {}
    for requirement in requirements:
        found = joint.attempt(requirement.problem_for(requirement.states), schedule, opts)
        if found is None:
            return None
        out[requirement.axis] = found
    return out
