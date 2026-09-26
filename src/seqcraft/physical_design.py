r"""
Ask SeqCraft to design part of your own composition, instead of writing that part yourself.

A :class:`PhysicalDesignScope` names a region of a repetition you are willing to hand over:

.. code-block:: text

    0            origin_s                  window                     echo
    |               |                        |                          |
    [---- before ---X------------[==== designed here ====]---- after ---X

You say where the spins were excited, what plays on either side of the region, which axes the
designer may use inside it, and which states your repetition comes in.  You then attach the same
physical intents a packaged kernel takes -- :class:`~seqcraft.FlowCompensation`,
:class:`~seqcraft.VelocityEncoding` -- and get back a family of fully timed blocks sharing one
echo time.

**This is opt-in and it is not a layer.**  Building events, modules and functions into a
``LogicBlock`` and compiling it needs none of this, and a composition does not have to become a
packaged module to use it.  `GRE2DTR` reaches the same designer through its own adapter; this is
the door for a composition that has no class.

.. code-block:: python

    exc = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0)
    arm = sc.modules.SpiralReadout(opts=opts, fov_mm=240.0, matrix=64, shots=8, variant='in')

    scope = sc.PhysicalDesignScope(
        origin_s=exc.time_to_center(),
        before=exc(rephase=False),                       # hand the slice rephasing over
        after=lambda angle: arm(angle_rad=angle, prephase=False),
        echo_in_after_s=arm.time_to_echo(0) - arm.prephaser_duration_s,
        axes=('x', 'y', 'z'),
        states=angles,
    )
    design = sc.design_repetition(scope, opts=opts,
                                  flow_comp=sc.FlowCompensation(axis=('x', 'y', 'z')))

    for angle in angles:
        scan.add(at, design.repetition(angle))

**What the designed region is asked for.**  On every axis it owns, it brings `k` back to the
origin at the echo -- which is what a prephaser, a winder and a slice rephaser all do -- and it
meets whatever the intents add on top.  A composition needing a *non-zero* `k` at the echo on a
designed axis, as a Cartesian phase encode does, is not expressible here yet; that is what the
packaged kernels are for.

**Handing a region over is explicit and it happens before anything is built.**  `before=exc(
rephase=False)` is you saying the slice rephasing is the designer's now.  Nothing inspects a
finished block and rewrites it, and nothing outside the declared region is read or changed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from .design import scope as _scope
from .design.joint import ORDERS, Schedule, measure_moment
from .design.logic import LogicBlock
from .design.timing import Raster
from .design.validation import require_axis
from .errors import ConfigurationError, format_error
from .modules import _augment

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

    from pypulseq.opts import Opts

    from .augmentation import FlowCompensation, VelocityEncoding

__all__ = ['PhysicalDesignScope', 'RepetitionDesign', 'design_repetition']

#: The state of a repetition that comes in only one.
_ONLY = _augment.ONE_STATE


def _blocks(piece: LogicBlock | Callable[[Any], LogicBlock], state: Any) -> LogicBlock:
    """`before` and `after` may be one block, or a function of the state."""
    return piece(state) if callable(piece) else piece


@dataclass(frozen=True)
class PhysicalDesignScope:
    """
    The region of your composition a designer may own, and what surrounds it.

    Parameters
    ----------
    origin_s
        When the spins were excited, measured from the start of the blocks you give below.  For
        an `Excitation` this is ``exc.time_to_center()``.  Moments are measured **from** this
        instant, because a gradient playing before it dephases nothing -- there is no transverse
        magnetisation yet.
    before
        What plays before the designed region: a ``LogicBlock``, or a function of the state
        returning one.  This is where a component hands its own winder over —
        ``exc(rephase=False)``.
    after
        What plays after it, on the same terms.  The echo is inside this.
    echo_in_after_s
        Where the echo sits inside `after`, measured from `after`'s own start.

        **Measure it for the configuration you are emitting.**  A module reports timings for its
        own block, and if you asked it to leave out a region -- ``prephase=False`` -- that
        region's duration is still in the number it reports.  `SpiralReadout.time_to_echo` is
        measured from a block containing its own prephaser, so a scope designing that prephaser
        subtracts `prephaser_duration_s`.  Getting this wrong does not raise; it moves the echo
        to an instant your sequence does not have, and the moments come out non-zero.
    axes
        The logical gradient channels the designer may use inside the region.  It brings `k` back
        to the origin at the echo on each of them, plus whatever the intents ask.
    states
        The states your repetition comes in -- interleaf angles, shot indices, whatever they are.
        SeqCraft never looks inside one.  Default: a single state.
    design_states
        Which states to **propose** the schedule from.  A proposal, not a proof: every state is
        realised and checked before a schedule is accepted, so a poor choice costs search time
        and never correctness.  Default: all of them, which is always safe and sometimes slow.
    min_window_s
        A floor on the designed region, when you already know one -- a module's own
        ``prephaser_duration_s``, say.  Default: none.
    """

    origin_s: float
    before: LogicBlock | Callable[[Any], LogicBlock]
    after: LogicBlock | Callable[[Any], LogicBlock]
    echo_in_after_s: float
    axes: str | Sequence[str]
    states: Sequence[Any] = (_ONLY,)
    design_states: Sequence[Any] | None = None
    min_window_s: float = 0.0

    def __post_init__(self) -> None:
        names = (self.axes,) if isinstance(self.axes, str) else tuple(self.axes)
        if not names:
            raise ConfigurationError(format_error(
                'a physical-design scope needs at least one axis to design on.',
                {'axes': self.axes},
                ["pass axes='y', or a tuple such as axes=('x', 'y')"],
            ))
        object.__setattr__(self, 'axes', tuple(require_axis(n, 'axes') for n in names))
        if not tuple(self.states):
            raise ConfigurationError(format_error(
                'a physical-design scope needs at least one state.', {'states': self.states},
                ['omit `states` for a repetition that comes in only one'],
            ))
        object.__setattr__(self, 'states', tuple(self.states))
        if self.design_states is not None:
            unknown = [s for s in self.design_states if s not in self.states]
            if unknown:
                raise ConfigurationError(format_error(
                    'design_states names states this scope does not have.',
                    {'unknown': tuple(unknown), 'states': self.states},
                    ['design_states proposes which of `states` sizes the schedule, so each one '
                     'has to be in `states`',
                     'omit it to propose all of them'],
                ))
            object.__setattr__(self, 'design_states', tuple(self.design_states))


@dataclass(frozen=True)
class RepetitionDesign:
    """
    A designed family: one schedule, and a block for each state.

    Attributes
    ----------
    te_s
        The echo time this family reached, from the excitation instant.  One value, shared by
        every state -- that is the point of designing the family rather than each state.
    window_s
        How long the designed region came out.
    """

    te_s: float
    window_s: float
    _scope: PhysicalDesignScope = field(repr=False)
    _opts: Opts = field(repr=False)
    _designed: _scope.ScopeDesign = field(repr=False)
    _states: tuple[Any, ...] = field(repr=False)
    _claims: tuple[Any, ...] = field(repr=False, default=())
    _encoding: tuple[Any, ...] = field(repr=False, default=(_ONLY,))

    def repetition(self, state: Any = _ONLY, *, encoding_state: int | None = None) -> LogicBlock:
        """
        One state's repetition: what you gave, with the designed region in the middle.

        Built **once**, from the finished design.  Nothing here edits a block that already
        exists.
        """
        key = self._key(state, encoding_state)
        out = LogicBlock('scope')
        out.add(0.0, _blocks(self._scope.before, state))
        for axis in self._scope.axes:
            for at, event in self._designed.events_for(axis, key):
                out.add(at, event)
        out.add(self._designed.schedule.window_start_s + self.window_s,
                _blocks(self._scope.after, state))
        return out

    def at(self, *, te_s: float) -> RepetitionDesign | None:
        """
        The same family designed again at a longer echo time, or ``None`` if it will not fit.

        For a protocol holding two families with different minima: take the longer echo time and
        ask both for it.  The design goes back through realisation rather than stretching what
        exists, because every moment is measured to the echo and the echo moved.
        """
        raster = float(self._opts.grad_raster_time)
        geometry, requirements = _translate(self._scope, self._opts, self._states,
                                            self._claims, self._encoding)
        grown = te_s - geometry.echo_at(0.0)
        steps = int(round(grown / raster))
        if steps < 1:
            return None
        found = _scope.realise_scope_at(geometry, requirements, self._opts,
                                        window_s=steps * raster)
        if found is None:
            return None
        return RepetitionDesign(geometry.echo_at(found.window_s), found.window_s,
                                self._scope, self._opts, found, self._states,
                                self._claims, self._encoding)

    def _key(self, state: Any, encoding_state: int | None) -> Any:
        """The internal key for one state, refusing an encoding state in either direction."""
        if self._encoding == (_ONLY,):
            if encoding_state is not None:
                _augment.refuse_state_mismatch('this scope', self._encoding, encoding_state)
            return state
        if encoding_state not in self._encoding:
            _augment.refuse_state_mismatch('this scope', self._encoding, encoding_state)
        return (state, encoding_state)


def _translate(scope: PhysicalDesignScope, opts: Opts, states: Sequence[Any],
               claims: Sequence[Any], encoding: Sequence[Any],
               ) -> tuple[_scope.ScopeGeometry, list[_scope.AxisRequirement]]:
    """
    The whole translation: a public declaration in, the internal substrate out.

    This is the only place the two vocabularies meet.  Nothing below learns that a scope existed,
    and nothing above learns what a `JointProblem` is.
    """
    raster = float(opts.grad_raster_time)
    probe = states[0][0] if isinstance(states[0], tuple) else states[0]
    geometry = _scope.ScopeGeometry(
        origin_s=scope.origin_s,
        window_start_s=float(Raster(raster).ceil(_blocks(scope.before, probe).duration)),
        tail_s=scope.echo_in_after_s,
    )

    def requirement(axis: str, targets: dict[Any, tuple[float, float | None]],
                    proposal: Sequence[Any]) -> _scope.AxisRequirement:
        def fixed(key: Any, schedule: Schedule) -> tuple[float, float]:
            """What the user's own blocks already put on this axis, where the schedule puts them."""
            state = key[0] if isinstance(key, tuple) else key
            played = LogicBlock()
            played.add(0.0, _blocks(scope.before, state))
            played.add(schedule.window_start_s + schedule.window_s,
                       _blocks(scope.after, state))
            return tuple(  # type: ignore[return-value]
                measure_moment(played, order, axis, origin_s=schedule.origin_s,
                               start_s=schedule.origin_s, end_s=schedule.endpoint_s)
                for order in ORDERS
            )

        return _scope.AxisRequirement(axis=axis, states=tuple(targets), design_states=proposal,
                                      target=lambda key: targets[key], fixed=fixed)

    return geometry, [requirement(axis, *_targets(scope, axis, states, claims, encoding))
                      for axis in scope.axes]


def _targets(scope: PhysicalDesignScope, axis: str, states: Sequence[Any],
             claims: Sequence[Any], encoding: Sequence[Any],
             ) -> tuple[dict[Any, tuple[float, float | None]], tuple[Any, ...]]:
    """
    What each state wants on one axis, and which states propose the schedule.

    The zeroth moment is always zero: the designed region brings `k` back to the origin at the
    echo.  The first is whatever the intents resolved to, or unconstrained when nothing asked.
    """
    wants_m1 = any(getattr(c, 'axis', None) == axis and getattr(c, 'order', None) == 1
                   for c in claims)
    targets: dict[Any, tuple[float, float | None]] = {}
    resolved = _joint_resolve(claims, encoding, axis) if wants_m1 else None
    for key in states:
        enc = key[1] if isinstance(key, tuple) else _ONLY
        targets[key] = (0.0, None if resolved is None else resolved[enc])
    proposal = scope.design_states
    if proposal is None:
        chosen: tuple[Any, ...] = tuple(targets)
    else:
        chosen = tuple(k for k in targets
                       if (k[0] if isinstance(k, tuple) else k) in proposal)
    return targets, chosen


def _joint_resolve(claims: Sequence[Any], encoding: Sequence[Any], axis: str) -> dict[Any, float]:
    """One absolute first-moment target per encoding state, from the claims."""
    from .design.joint import resolve_claims
    return resolve_claims(claims, tuple(encoding), axis=axis, order=1,
                          base={key: 0.0 for key in encoding})


def design_repetition(scope: PhysicalDesignScope, *, opts: Opts,
                      flow_comp: FlowCompensation | None = None,
                      velocity_encode: VelocityEncoding | None = None) -> RepetitionDesign:
    """
    Design one repetition family, and return a block-maker for every state.

    The same physical intents a packaged kernel takes.  They say *what relation must hold*; where
    it is realised is this scope's business, and which waveform realises it is the designer's.

    Raises
    ------
    ConfigurationError
        If an intent names an axis the scope does not own, or if no schedule can serve every
        state.  The scope owns what `axes` lists and nothing else, which is why an intent on
        another axis is a refusal rather than a silent omission.
    """
    from .augmentation import FlowCompensation, VelocityEncoding

    _augment.require_intent_type(flow_comp, FlowCompensation, 'flow_comp')
    _augment.require_intent_type(velocity_encode, VelocityEncoding, 'velocity_encode')
    owners = dict.fromkeys(scope.axes, 'joint')
    for intent, axis in _augment.axes_claimed((flow_comp, velocity_encode)):
        if axis not in owners:
            _augment.refuse_unowned_axis(
                'this physical-design scope', intent, axis, owners,
                f'{axis!r} is not in this scope\'s `axes`, so it owns no adjustable region there')

    claims, encoding = _augment.claims_and_states(flow_comp, velocity_encode)
    states: tuple[Any, ...] = (
        tuple(scope.states) if encoding == (_ONLY,)
        else tuple((state, key) for state in scope.states for key in encoding))
    geometry, requirements = _translate(scope, opts, states, claims, encoding)
    designed = _scope.design_scope(geometry, requirements, opts,
                                   min_window_s=scope.min_window_s)
    return RepetitionDesign(geometry.echo_at(designed.window_s), designed.window_s,
                            scope, opts, designed, states, claims, encoding)
