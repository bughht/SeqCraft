r"""
**Internal.**  Public intent in, internal claims out.  The only place the two vocabularies meet.

.. code-block:: text

    public intent      FlowCompensation / VelocityEncoding
           |   claims_and_states()     <- the one dispatch point, above the numerics
           v
    internal claims    CommonModeClaim / DifferenceClaim
           |   _joint.resolve_claims
           v
    absolute targets   one (m0, m1) per state
           |
           v
    realisation        sees numbers.  Never learns an augmentation existed

Keeping that single point is what lets a third augmentation be added without touching the schedule
search, the realisation families or either kernel's design method.  Nothing below this module
imports :mod:`seqcraft.augmentation`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ..augmentation import FlowCompensation, VelocityEncoding
from ..errors import ConfigurationError, format_error
from ._joint import CommonModeClaim, DifferenceClaim

if TYPE_CHECKING:
    from collections.abc import Sequence

__all__ = [
    'ONE_STATE',
    'claims_and_states',
    'flow_comp_for',
    'refuse_state_mismatch',
    'refuse_unowned_axis',
    'refuse_unsupported_route',
    'require_intent_type',
]

#: The state key of a repetition that acquires one state, which is every repetition that is not
#: velocity encoding.  A name rather than ``None`` so that the state-indexed maps downstream have
#: one shape, and so a refusal can print something a reader recognises.
ONE_STATE = 'only'


def claims_and_states(
    flow_comp: FlowCompensation | None,
    velocity_encode: VelocityEncoding | None,
) -> tuple[tuple[object, ...], tuple[object, ...]]:
    """
    Translate the public intents into claims, and say which states they are stated over.

    >>> claims, states = claims_and_states(FlowCompensation(axis='y'), None)
    >>> claims
    (CommonModeClaim(axis='y', order=1, value=0.0),)
    >>> states
    ('only',)

    Several axes are several independent requirements, one per waveform:

    >>> [c.axis for c in claims_and_states(FlowCompensation(axis=('y', 'z')), None)[0]]
    ['y', 'z']

    Velocity encoding is what makes a repetition two-state:

    >>> claims, states = claims_and_states(None, VelocityEncoding(venc_m_s=1.5, axis='z'))
    >>> states
    (1, -1)
    >>> round(claims[0].value, 6)
    0.333333
    """
    claims: list[object] = []
    states: tuple[object, ...] = (ONE_STATE,)

    if velocity_encode is not None:
        states = velocity_encode.states
        claims.append(DifferenceClaim(velocity_encode.axis, 1,
                                      velocity_encode.delta_m1_s_per_m, states))
    if flow_comp is not None:
        # CommonModeClaim, never a per-state zero.  That is the whole of FlowCompensation's
        # contract: with a difference claim beside it the two resolve to +-delta/2, and nothing
        # here computes that half.  There is deliberately no branch on how many states there are
        # -- such a branch is exactly how a per-state zero would creep back and break the
        # composed case.
        claims.extend(CommonModeClaim(axis, 1) for axis in flow_comp.axes)

    return tuple(claims), states


def require_intent_type(value: object, expected: type, keyword: str) -> None:
    """Refuse the wrong intent in a keyword, naming where it should have gone."""
    if value is None or isinstance(value, expected):
        return
    other = {'flow_comp': 'velocity_encode', 'velocity_encode': 'flow_comp'}.get(keyword, '')
    fixes = [f'pass a {expected.__name__}, or None']
    if other and isinstance(value, (FlowCompensation, VelocityEncoding)):
        fixes.insert(0, f'{type(value).__name__} goes in {other}=')
    raise ConfigurationError(format_error(
        f'{keyword} expects a {expected.__name__}, and got {type(value).__name__}.',
        {keyword: type(value).__name__, 'expected': expected.__name__},
        fixes,
    ))


def flow_comp_for(flow_comp: FlowCompensation | None,
                  axes: set[str]) -> FlowCompensation | None:
    """The part of a flow-compensation intent that `axes` covers, or ``None``."""
    if flow_comp is None:
        return None
    wanted = tuple(a for a in flow_comp.axes if a in axes)
    return FlowCompensation(axis=wanted) if wanted else None


def refuse_unowned_axis(component: str, axis: str, owners: dict[str, str], why: str) -> None:
    """
    Refuse an axis this repetition has no owner for, in physical terms.

    `why` is the sequence-specific sentence -- which gradient that axis carries here, and who
    realises it -- because "not supported" is not a reason a caller can act on.
    """
    raise ConfigurationError(format_error(
        f'{component} cannot apply flow compensation on {axis!r}.',
        {'axis': axis, 'axes it can use': tuple(owners)},
        [why,
         f'this repetition designs {tuple(owners)}; an axis it does not own is realised by the '
         f'leaf that does own it, and cannot be augmented from here'],
    ))


def refuse_unsupported_route(component: str, axis: str, joint_axes: Sequence[str]) -> None:
    """
    Refuse velocity encoding on an axis whose owner realises one repetition at a time.

    The reason is the current local realisation's scope, **not** an impossibility: the readout's
    solve implements first-moment nulling, one fixed target, and a difference between two
    acquisitions needs a state-resolved one.  Extending it would make this work with no change to
    the routing, so the message must not say it cannot be done.
    """
    raise ConfigurationError(format_error(
        f'{component} cannot velocity encode on {axis!r}.',
        {'axis': axis, 'axes it can velocity encode on': tuple(joint_axes)},
        ['the readout axis is designed by its readout, which solves one repetition at a time and '
         'nulls its first moment rather than aiming it at a value',
         f'velocity encoding is a difference between two acquisitions, so it needs an axis this '
         f'repetition designs jointly: {tuple(joint_axes)}'],
    ))


def refuse_state_mismatch(component: str, states: Sequence[object], given: object) -> None:
    """
    Refuse a state that means nothing here, in either direction.

    Both directions, because silence is worse than a refusal in both: a missing state realises an
    arbitrary one, and an ignored state hands the caller two identical repetitions that they find
    out about when the subtraction comes back zero.
    """
    if tuple(states) == (ONE_STATE,):
        raise ConfigurationError(format_error(
            f'{component} does not velocity encode, so encoding_state={given!r} has no meaning '
            f'here.',
            {'encoding_state': given},
            ['drop encoding_state -- this repetition has one state and build() already realises '
             'it',
             'or pass velocity_encode=... if this acquisition is meant to be phase contrast'],
        ))
    if given is None:
        raise ConfigurationError(format_error(
            f'{component} velocity encodes, so build() needs to be told which state it is '
            f'realising.',
            {'states': tuple(states)},
            [f'pass encoding_state={states[0]!r} or encoding_state={states[1]!r}',
             'acquire both, and subtract their phase to get velocity'],
        ))
    raise ConfigurationError(format_error(
        f'encoding_state={given!r} is not one of this repetition\'s states.',
        {'encoding_state': given, 'states': tuple(states)},
        ['iterate the states the intent defines:  for state in venc.states:'],
    ))
