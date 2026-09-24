r"""
**Provisional, internal.** Repetition-level joint moment design -- requirements, then realisation.

Nothing here is public and nothing here is final.  It exists to test one boundary:

.. code-block:: text

    fixed physical contribution
    + adjustable waveform family in a window
    + a target at a semantic instant
            -> one linear solve
            -> areas

A requirement is a **moment at an instant, measured from a named origin**.  Both instants are
required, because under a shift of origin :math:`m_1' = m_1 - \Delta t\, m_0` -- so a first moment
is not a number until the origin is named.  ``order == 0`` does not read the origin; that is not a
reason to leave it out of the semantics.

Two layers, deliberately separated:

.. code-block:: text

    resolve()   turns per-augmentation claims into ONE absolute target per repetition state
    solve()     realises one state's absolute target, and never learns where it came from

That separation is what keeps velocity encoding's state semantics -- a *difference* between
encoding states -- out of the waveform solve, which only ever sees a number.

What this is not: a plugin bus, a constraint language, a compiler change, or an optimiser.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ..design.events import knots_of, pwl_moment
from ..design.logic import flatten
from ..errors import ConfigurationError, format_error
from ._support import ceil_raster

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping, Sequence

    from pypulseq.opts import Opts

    from ..design.events import Event
    from ..design.logic import LogicBlock

__all__ = [
    'CommonModeClaim',
    'DifferenceClaim',
    'JointProblem',
    'JointSolution',
    'measure_moment',
    'resolve',
    'solve',
]

#: Moments below this are zero for the purpose of choosing a lobe, in the axis's own units.
_NEGLIGIBLE = 1e-9


# ------------------------------------------------------------------ the independent validator
def measure_moment(tree: LogicBlock, order: int, axis: str, *,
                   origin_s: float, start_s: float, end_s: float) -> float:
    """
    The `order`-th moment of **every** gradient on `axis` over ``[start_s, end_s]``, about
    `origin_s`.

    Integrated from each event's own piecewise-linear knots and summed across events, which is the
    part no per-event check reaches.  Used by tests rather than by the solve: the validator must
    not ask the designer what it thinks it built.
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
            times = np.concatenate(([low], times[inside], [high]))
            amps = np.concatenate((
                [np.interp(low, *knots_of(event, at))], amps[inside],
                [np.interp(high, *knots_of(event, at))],
            ))
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


def resolve(claims: Iterable[object], states: Sequence[object], *,
            axis: str, order: int, base: Mapping[object, float]) -> dict[object, float]:
    """
    Turn claims into one **absolute** target per state, or refuse.

    `base` is what the sequence itself wants at this order -- the encode area at order 0, and
    nothing at order 1.  A claim on the difference and a claim on the common mode constrain
    different components of the same state-indexed vector, so they compose without a precedence
    rule; two claims on the *same* component do not, and raise.

    Composing the two for a two-state acquisition gives the symmetric pair, which is therefore
    derived rather than assumed:

    >>> resolve([DifferenceClaim('y', 1, 0.5, (1, -1)), CommonModeClaim('y', 1)],
    ...         [1, -1], axis='y', order=1, base={1: 0.0, -1: 0.0})
    {1: 0.25, -1: -0.25}
    """
    mine = [c for c in claims if getattr(c, 'axis', None) == axis and getattr(c, 'order', None) == order]
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
        a, b = claim.states
        if a not in out or b not in out:
            _refuse_missing_state(claim, states)
        half = (claim.value - (out[a] - out[b])) / 2.0
        out[a] += half
        out[b] -= half
    return out


def _refuse_missing_state(claim: DifferenceClaim, states: Sequence[object]) -> None:
    """A difference is between two states, so both have to be acquired."""
    msg = format_error(
        f'a difference claim on moment {claim.order} of {claim.axis} needs states '
        f'{claim.states}, and this acquisition has {tuple(states)}.',
        {'axis': claim.axis, 'order': claim.order,
         'claim_states': claim.states, 'acquired_states': tuple(states)},
        ['acquire both states -- a phase difference needs two acquisitions to subtract',
         'or drop the claim, which leaves this moment to whatever else constrains it'],
    )
    raise ConfigurationError(msg)


def _refuse_duplicate(component: str, axis: str, order: int, count: int) -> None:
    msg = format_error(
        f'{count} augmentations claim the {component} of moment {order} on {axis}.',
        {'axis': axis, 'order': order, 'component': component, 'claims': count},
        ['one augmentation owns one component of one moment on one axis',
         'a second claim is a disagreement about physics, not something to order by precedence'],
    )
    raise ConfigurationError(msg)


# ------------------------------------------------------------------------ the realisation layer
@dataclass(frozen=True)
class JointProblem:
    """
    One axis of one repetition: what is already there, what is wanted, and where it may be built.

    Every field is a number the repetition already knows.  No gradient event crosses this
    boundary in either direction.
    """

    axis: str
    origin_s: float
    endpoint_s: float
    window_start_s: float
    fixed: tuple[float, ...]
    targets: Mapping[object, tuple[float, ...]]

    @property
    def orders(self) -> tuple[int, ...]:
        return tuple(range(len(self.fixed)))


@dataclass(frozen=True)
class JointSolution:
    """The shared window, and one set of lobe areas per state."""

    half_s: float
    areas: Mapping[object, tuple[float, ...]]

    @property
    def duration_s(self) -> float:
        return 2.0 * self.half_s


def response(problem: JointProblem, half_s: float) -> np.ndarray:
    """
    The moment response of the basis: column `j` is lobe `j`'s moments about the origin.

    Two lobes of duration `half_s` played back to back from the window start.  A symmetric
    trapezoid's `order`-th moment about an external origin is its area times its centroid to that
    power, exactly, so the matrix needs no waveform -- which is what makes the solve independent
    of how the lobe is later realised.
    """
    centres = [problem.window_start_s + (0.5 + n) * half_s - problem.origin_s for n in range(2)]
    return np.array([[c ** order for c in centres] for order in problem.orders], dtype=float)


def solve(problem: JointProblem, opts: Opts, *, half_s: float | None = None) -> JointSolution:
    """
    Solve every state in one shared window, and return the areas.

    The window is shared because a repetition whose TE moved with the state would put a contrast
    gradient across the acquisition that no k-space check would show.  Its length is the shortest
    that serves the **limiting** state, found by enumeration -- never assumed to be an extreme of
    any index, because the limit is a result of the coupling.
    """
    raster = float(opts.grad_raster_time)
    shortest = half_s if half_s is not None else _shortest_half_s(problem, opts, raster)
    return JointSolution(
        half_s=shortest,
        areas={state: _areas_for(problem, state, shortest) for state in problem.targets},
    )


def realise(problem: JointProblem, solution: JointSolution, state: object,
            opts: Opts) -> list[Event]:
    """Build one state's lobes.  Zero-area lobes are dropped rather than handed to pypulseq."""
    return [
        pp.make_trapezoid(channel=problem.axis, area=area, duration=solution.half_s, system=opts)
        for area in solution.areas[state] if abs(area) > _NEGLIGIBLE
    ]


def _areas_for(problem: JointProblem, state: object, half_s: float) -> tuple[float, ...]:
    """The lobe areas for one state: ``M a = target - fixed``."""
    wanted = np.asarray(problem.targets[state], dtype=float) - np.asarray(problem.fixed, dtype=float)
    return tuple(float(a) for a in np.linalg.solve(response(problem, half_s), wanted))


def _fits(problem: JointProblem, opts: Opts, half_s: float) -> bool:
    """Whether every state's lobes fit a window of this length, inside the gradient limits."""
    if half_s <= 0.0:
        return False
    for state in problem.targets:
        for area in _areas_for(problem, state, half_s):
            if abs(area) <= _NEGLIGIBLE:
                continue
            try:
                lobe = pp.make_trapezoid(
                    channel=problem.axis, area=area, duration=half_s, system=opts)
            except (ValueError, ZeroDivisionError):
                return False
            peak = abs(float(lobe.amplitude))
            if peak > float(opts.max_grad) * (1.0 + 1e-9):
                return False
            if peak > float(opts.max_slew) * float(lobe.rise_time) * (1.0 + 1e-9):
                return False
    return True


def _shortest_half_s(problem: JointProblem, opts: Opts, raster: float) -> float:
    """Bracket by doubling, then bisect, then confirm the step below does not fit."""
    high = raster
    for _ in range(40):
        if _fits(problem, opts, high):
            break
        high *= 2.0
    else:  # pragma: no cover - unreachable on any real gradient system
        _refuse_infeasible(problem)
    low = high / 2.0
    while high - low > raster:
        middle = ceil_raster((low + high) / 2.0, raster)
        if middle >= high:
            break
        if _fits(problem, opts, middle):
            high = middle
        else:
            low = middle
    return float(high)


def _refuse_infeasible(problem: JointProblem) -> None:  # pragma: no cover - unreachable
    msg = format_error(
        f'no window on {problem.axis} carries these moments.',
        {'axis': problem.axis, 'targets': dict(problem.targets)},
        ['relax the target', 'or widen the gradient limits if the scanner allows it'],
    )
    raise ConfigurationError(msg)
