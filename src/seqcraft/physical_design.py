r"""
Ask SeqCraft to design part of your own composition, instead of writing that part yourself.

A :class:`PhysicalDesignScope` names a region you are willing to hand over:

.. code-block:: text

    0            origin_s                  window                     echo
    |               |                        |                          |
    [---- before ---X------------[==== designed here ====]---- after ---X

You say where the pathway starts, what plays on either side of the region, which axes the designer
may use inside it, and which states your family comes in.  You then attach the same physical
intents a packaged kernel takes -- :class:`~seqcraft.FlowCompensation`,
:class:`~seqcraft.VelocityEncoding` -- and get back one design covering every state, at one echo
time.

**This is opt-in and it is not a layer.**  Building events, modules and functions into a
``LogicBlock`` and compiling it needs none of this, and a composition does not have to become a
packaged module to use it.  `GRE2DTR` reaches the same designer through its own adapter; this is
the door for a composition that has no class.

.. code-block:: python

    scope = sc.PhysicalDesignScope(
        origin_s=exc.time_to_center(),
        before=exc(rephase=False),                       # hand the slice rephasing over
        after=lambda angle: arm(angle_rad=angle, prephase=False),
        echo_in_after_s=arm.time_to_echo(0) - arm.prephaser_duration_s,
        axes=('x', 'y', 'z'),
        states=angles,
    )
    design = scope.design(opts=opts, flow_comp=sc.FlowCompensation(axis=('x', 'y', 'z')))

    for angle in angles:
        scan.add(at, design.build(angle))

**What the designed region is asked for.**  On each axis it owns it puts `k` where
:attr:`~PhysicalDesignScope.k_at_echo` says at the echo -- zero by default, which is what a
prephaser, a winder and a slice rephaser all do -- and meets whatever the intents add on top.

**Handing a region over is explicit and it happens before anything is built.**  ``before=exc(
rephase=False)`` is you saying the slice rephasing is the designer's now.  Nothing inspects a
finished block and rewrites it, and nothing outside the declared region is read or changed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from .design import scope as _scope
from .design.joint import ORDERS, Schedule, measure_moment, resolve_claims
from .design.logic import LogicBlock
from .design.timing import Raster
from .design.validation import require_axis
from .errors import ConfigurationError, format_error
from .modules import _augment

if TYPE_CHECKING:
    from collections.abc import Callable, Hashable, Sequence

    from pypulseq.opts import Opts

    from .augmentation import FlowCompensation, VelocityEncoding

__all__ = ['PhysicalDesign', 'PhysicalDesignScope']

#: The encoding state of a family that has only one.
_ONLY = _augment.ONE_STATE


@dataclass(frozen=True)
class _Key:
    """
    **Internal.**  One state of the designed family: the caller's state, and an encoding state.

    A wrapper rather than a tuple, and that is the whole point of it.  A caller's state may itself
    be a tuple -- ``(ky, kz)``, ``(partition, line)`` -- so tuple shape cannot be read as a type
    tag.  The caller's object is carried through untouched and handed back to `before` and
    `after` exactly as it was given.
    """

    state: Hashable
    encoding: Any = _ONLY


def _piece(piece: LogicBlock | Callable[[Any], LogicBlock], state: Any) -> LogicBlock:
    """`before` and `after` may be one block, or a function of the caller's own state."""
    return piece(state) if callable(piece) else piece


def _refuse(headline: str, fields: dict[str, Any], fixes: Sequence[str]) -> None:
    raise ConfigurationError(format_error(headline, fields, list(fixes)))


@dataclass(frozen=True)
class PhysicalDesignScope:
    """
    The region of your composition a designer may own, and what surrounds it.

    Parameters
    ----------
    origin_s
        The **semantic origin** of the pathway whose moments this scope designs, measured from the
        start of the blocks below.  Moments are integrated from it.

        For the excitation-to-echo pathway of a gradient echo -- which is what the shipped
        examples are -- that origin is the RF effective centre, ``exc.time_to_center()``, and
        gradient played before it does not act on the signal this scope is designing for, because
        the transverse magnetisation does not exist yet.  **That reasoning is specific to this
        pathway.**  This first public scope assumes the pathway begins at `origin_s`; RF history
        outside the scope, and the several coherence pathways a refocusing pulse creates, are not
        modelled.
    before
        What plays before the designed region: a ``LogicBlock``, or a function of your state
        returning one.  This is where a component hands its own winder over --
        ``exc(rephase=False)``.  Its waveform may vary with the state; its **duration may not**,
        because the family shares one window start.
    after
        What plays after it, on the same terms.  The echo is inside this.
    echo_in_after_s
        Where the echo sits inside `after`, measured from `after`'s own start.

        **Measure it for the configuration you are emitting.**  A module reports timings for its
        own block, and if you asked it to leave out a region -- ``prephase=False`` -- that
        region's duration is still in the number it reports. `SpiralReadout.time_to_echo` is
        measured from a block containing its own prephaser, so a scope designing that prephaser
        subtracts `prephaser_duration_s`.
    axes
        The logical gradient channels the designer may use inside the region.
    states
        The states your family comes in -- interleaf angles, ``(ky, kz)`` pairs, whatever they
        are.  SeqCraft never looks inside one; they need only be **hashable**, and they are handed
        back to `before` and `after` exactly as given.  Default: a single state.
    design_states
        Which states to **propose** the schedule from.  A proposal, not a proof: every state is
        realised and checked before a schedule is accepted, so a poor choice costs search time and
        never correctness.  Default: all of them, which is always safe and sometimes slow.
    k_at_echo
        Where the designed region should leave `k` at the echo, in ``1/m``, as
        ``k_at_echo(state, axis)``.  Default: zero on every axis, which is what a prephaser, a
        winder and a slice rephaser all do.

        A family that *encodes* with the region it owns says so here -- a phase encode is
        ``lambda state, axis: pe.k_per_m(state) if axis == 'y' else 0.0``.  This is the zeroth
        moment of the whole scope at the echo; the first moment is what the physical intents are
        for.
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
    k_at_echo: Callable[[Any, str], float] | None = None
    min_window_s: float = 0.0

    def __post_init__(self) -> None:
        names = (self.axes,) if isinstance(self.axes, str) else tuple(self.axes)
        if not names:
            _refuse('a physical-design scope needs at least one axis to design on.',
                    {'axes': self.axes}, ["pass axes='y', or a tuple such as axes=('x', 'y')"])
        object.__setattr__(self, 'axes', tuple(require_axis(n, 'axes') for n in names))

        states = tuple(self.states)
        if not states:
            _refuse('a physical-design scope needs at least one state.', {'states': self.states},
                    ['omit `states` for a family that comes in only one'])
        for state in states:
            try:
                hash(state)
            except TypeError:
                _refuse(
                    'a state has to be hashable, because the design is stored per state.',
                    {'state': repr(state), 'type': type(state).__name__},
                    ['use a tuple rather than a list, or a frozen dataclass rather than a dict',
                     'a state is an identity SeqCraft never looks inside -- it only has to be '
                     'usable as a key'],
                )
        object.__setattr__(self, 'states', states)

        for name in ('origin_s', 'echo_in_after_s', 'min_window_s'):
            value = float(getattr(self, name))
            if value < 0.0:
                _refuse(f'{name} = {value:g} s is negative.', {name: value},
                        [f'{name} is measured forward from the start of what you emit'])

        if self.design_states is not None:
            proposal = tuple(self.design_states)
            if not proposal:
                _refuse('design_states is empty, so nothing would size the schedule.',
                        {'states': states},
                        ['omit design_states to propose every state',
                         'or name the few you expect to be hardest'])
            unknown = [s for s in proposal if s not in states]
            if unknown:
                _refuse('design_states names states this scope does not have.',
                        {'unknown': tuple(unknown), 'states': states},
                        ['design_states proposes which of `states` sizes the schedule, so each '
                         'one has to be in `states`', 'omit it to propose all of them'])
            object.__setattr__(self, 'design_states', proposal)

    # ------------------------------------------------------------------ designing
    def design(self, *, opts: Opts, flow_comp: FlowCompensation | None = None,
               velocity_encode: VelocityEncoding | None = None) -> PhysicalDesign:
        """
        Design this scope, and return one design covering every state.

        The physical intents are the same objects a packaged kernel takes.  They say *what
        relation must hold*; where it applies is this scope's business, and which waveform
        realises it is the designer's.

        Raises
        ------
        ConfigurationError
            If an intent names an axis this scope does not own, if `before` varies in duration
            across the family, or if no schedule can serve every state.
        """
        return _design(self, opts, flow_comp, velocity_encode, te_request_s=None)


@dataclass(frozen=True)
class PhysicalDesign:
    """
    A designed family: one schedule, and a block for each state.

    Attributes
    ----------
    te_s
        The echo time this family reached, from the semantic origin.  One value shared by every
        state -- that is the point of designing the family rather than each state.
    window_s
        How long the designed region came out.
    """

    te_s: float
    window_s: float
    scope: PhysicalDesignScope
    _opts: Opts = field(repr=False)
    _designed: _scope.ScopeDesign = field(repr=False)
    _keys: tuple[_Key, ...] = field(repr=False)
    _intents: tuple[Any, Any] = field(repr=False)
    _encoding: tuple[Any, ...] = field(repr=False)

    def build(self, state: Any = _ONLY, *, encoding_state: int | None = None) -> LogicBlock:
        """
        One state's block: what you declared, with the designed region in the middle.

        Built **once**, from the finished design.  Nothing here edits a block that already
        exists, and what comes back is exactly ``before + designed + after`` -- a spoiler, a TR
        fill or anything else outside the scope is yours to add.
        """
        key = self._key(state, encoding_state)
        out = LogicBlock('scope')
        out.add(0.0, _piece(self.scope.before, key.state))
        for axis in self.scope.axes:
            for at, event in self._designed.events_for(axis, key):
                out.add(at, event)
        out.add(self._designed.schedule.window_start_s + self.window_s,
                _piece(self.scope.after, key.state))
        return out

    def __call__(self, state: Any = _ONLY, *, encoding_state: int | None = None) -> LogicBlock:
        """:meth:`build`, so a design reads like the modules beside it."""
        return self.build(state, encoding_state=encoding_state)

    def at(self, *, te_s: float) -> PhysicalDesign | None:
        """
        The same family designed again at a longer echo time, or ``None`` if it will not fit.

        This **redesigns**; it does not stretch.  The requested echo time goes back through the
        same search as an explicit TE, so the designer is free to keep the window it already had
        and let the extra time become fill in front of it, rather than being forced to spend all
        of it on a wider region.  Every moment is measured to the echo, and the echo moved.

        The achieved :attr:`te_s` is **never shorter than requested**: a request between two
        legal instants comes back at the first one at or above it.

        For a protocol holding two families with different minima: take the longer echo time and
        ask both for it.
        """
        if te_s < self.te_s:
            return None
        try:
            return _design(self.scope, self._opts, *self._intents, te_request_s=te_s)
        except ConfigurationError:
            return None

    # ------------------------------------------------------------------ internals
    def _key(self, state: Any, encoding_state: int | None) -> _Key:
        """The internal key for one state, refusing an encoding state in either direction."""
        if self._encoding == (_ONLY,):
            if encoding_state is not None:
                _augment.refuse_state_mismatch('this scope', self._encoding, encoding_state)
            return _Key(state)
        if encoding_state not in self._encoding:
            _augment.refuse_state_mismatch('this scope', self._encoding, encoding_state)
        return _Key(state, encoding_state)


# ------------------------------------------------------------------ the translation
def _design(scope: PhysicalDesignScope, opts: Opts, flow_comp: FlowCompensation | None,
            velocity_encode: VelocityEncoding | None, *,
            te_request_s: float | None) -> PhysicalDesign:
    """Public declaration in, internal substrate out, and a design back.  The only crossing."""
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
    keys = tuple(_Key(state, enc) for state in scope.states for enc in encoding)

    geometry = _geometry(scope, opts)
    requirements = [_requirement(scope, axis, keys, claims, encoding, geometry)
                    for axis in scope.axes]
    designed = _scope.design_scope(geometry, requirements, opts,
                                   min_window_s=scope.min_window_s,
                                   te_request_s=te_request_s)
    return PhysicalDesign(geometry.echo_at(designed.window_s, designed.fill_s),
                          designed.window_s, scope, opts, designed, keys,
                          (flow_comp, velocity_encode), encoding)


def _geometry(scope: PhysicalDesignScope, opts: Opts) -> _scope.ScopeGeometry:
    """
    Where the region sits, and the one structural invariant that has to hold to place it.

    `before` may differ in waveform between states, because the designed region is solved per
    state.  It may **not** differ in duration: the family shares one window start, so a state
    whose `before` is longer would either overlap the designed region or leave a gap that nothing
    accounts for.  Checked here rather than discovered as a residual.
    """
    raster = Raster(float(opts.grad_raster_time))
    durations = {state: _piece(scope.before, state).duration for state in scope.states}
    spread = max(durations.values()) - min(durations.values())
    if spread > float(opts.grad_raster_time) * 1e-6:
        longest = max(durations, key=lambda s: durations[s])
        shortest = min(durations, key=lambda s: durations[s])
        _refuse(
            'the blocks before the designed region differ in duration between states, and the '
            'family shares one window start.',
            {'longest state': repr(longest), 'longest': f'{durations[longest] * 1e6:.3f} us',
             'shortest state': repr(shortest), 'shortest': f'{durations[shortest] * 1e6:.3f} us'},
            ['`before` may vary in waveform between states, but not in duration',
             'pad the shorter ones to a common duration, or move what varies into `after`'],
        )
    return _scope.ScopeGeometry(
        origin_s=scope.origin_s,
        window_start_s=float(raster.ceil(max(durations.values()))),
        tail_s=scope.echo_in_after_s,
    )


def _requirement(scope: PhysicalDesignScope, axis: str, keys: Sequence[_Key],
                 claims: Sequence[Any], encoding: Sequence[Any],
                 geometry: _scope.ScopeGeometry) -> _scope.AxisRequirement:
    """One axis' requirement, with `target` and `fixed` derived rather than asked for."""
    wants_m1 = any(getattr(c, 'axis', None) == axis and getattr(c, 'order', None) == 1
                   for c in claims)
    k_at_echo = scope.k_at_echo or (lambda _state, _axis: 0.0)

    targets: dict[_Key, tuple[float, float | None]] = {}
    for key in keys:
        base = float(k_at_echo(key.state, axis))
        first = None
        if wants_m1:
            first = resolve_claims(claims, tuple(encoding), axis=axis, order=1,
                                   base=dict.fromkeys(encoding, 0.0))[key.encoding]
        targets[key] = (base, first)

    def fixed(key: _Key, schedule: Schedule) -> tuple[float, float]:
        """What the caller's own blocks already put on this axis, where the schedule puts them."""
        played = LogicBlock()
        played.add(0.0, _piece(scope.before, key.state))
        played.add(schedule.window_start_s + schedule.window_s, _piece(scope.after, key.state))
        return tuple(  # type: ignore[return-value]
            measure_moment(played, order, axis, origin_s=schedule.origin_s,
                           start_s=schedule.origin_s, end_s=schedule.endpoint_s)
            for order in ORDERS
        )

    if scope.design_states is None:
        proposal: tuple[_Key, ...] = tuple(targets)
    else:
        wanted = list(scope.design_states)
        proposal = tuple(k for k in targets if any(k.state is w or k.state == w for w in wanted))
    return _scope.AxisRequirement(axis=axis, states=tuple(targets), design_states=proposal,
                                  target=lambda key: targets[key], fixed=fixed)
