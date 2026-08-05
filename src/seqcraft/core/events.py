"""
Safe derivation of pypulseq events.

pypulseq events are plain :class:`types.SimpleNamespace` objects and seqcraft does **not**
wrap them: wrapping would be a large API surface for no benefit and would break interop
with the rest of the pulseq ecosystem.  What seqcraft does add is one sanctioned way to
produce a modified copy, because two properties of pypulseq 1.5 make ad-hoc copying
dangerous.

1. ``Sequence.set_block`` writes a registration cache **onto the event object**,
   ``event._pypulseq_sequence_event_cache``, keyed by ``id(sequence)``.  A
   :func:`copy.deepcopy` carries that stale key along, and CPython reuses addresses -- so
   a copied event can be silently mis-attributed to a *different* ``Sequence`` that
   happens to land at the same address.  That is a wrong-sequence hazard, not a
   performance note.  **Never deepcopy an event.**
2. Events carry an ``id`` attribute once registered.  A copy that keeps it will be taken
   for the original.

:func:`derive` therefore shallow-copies, strips both, and applies the requested changes.
Waveform arrays are shared by reference on purpose -- they are never mutated in place
(the contract tests assert this), and sharing is what lets pypulseq's own per-event
registration cache hit for repeated events.

Examples
--------
>>> import pypulseq as pp
>>> from pypulseq.opts import Opts
>>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=170, slew_unit='T/m/s')
>>> g = pp.make_trapezoid(channel='x', area=100.0, system=o)
>>> g2 = derive(g, delay=1e-4)
>>> g2.delay, g.delay
(0.0001, 0.0)
"""

from __future__ import annotations

import copy
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any

import numpy as np
import pypulseq as pp

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from pypulseq.opts import Opts

__all__ = [
    'Event',
    'GRAD_TYPES',
    'channels_of',
    'content_hash',
    'derive',
    'duration_of',
    'kinds_of',
    'moment_of',
    'sanitise',
    'trapz',
    'waveform_of',
]

#: A pypulseq event.  Alias rather than a class: seqcraft never wraps these.
Event = SimpleNamespace

#: Event ``type`` values that carry a gradient on a channel.
GRAD_TYPES = frozenset({'trap', 'grad'})

_STRIP = ('id', 'shape_IDs', '_pypulseq_sequence_event_cache')

# numpy renamed trapz -> trapezoid in 2.0; support both so numpy>=1.24 stays valid.
trapz = getattr(np, 'trapezoid', None) or np.trapz


def derive(event: Event, **changes: Any) -> Event:
    """
    Return a modified shallow copy of `event`.

    Parameters
    ----------
    event
        Any pypulseq event.
    **changes
        Attributes to overwrite on the copy.

    Returns
    -------
    SimpleNamespace
        A new event.  `event` itself is never touched.

    Raises
    ------
    AttributeError
        If `changes` names an attribute the event does not already have.  Silently
        adding a field would hide a typo such as ``phase_offest=``.
    """
    out = copy.copy(event)
    for name in _STRIP:
        if hasattr(out, name):
            delattr(out, name)
    for key, value in changes.items():
        if not hasattr(out, key):
            msg = (
                f'{getattr(event, "type", "event")!r} event has no attribute {key!r}; '
                f'available: {sorted(vars(event))}'
            )
            raise AttributeError(msg)
        setattr(out, key, value)
    return out


def sanitise(event: Event) -> Event:
    """
    Strip registration state from an event obtained from ``Sequence.get_block``.

    Needed when importing an existing sequence: events read back out of a ``Sequence``
    carry that sequence's ids and cache, which must not travel into a new one.
    """
    return derive(event)


def duration_of(events: Sequence[Event], *, explicit: float | None = None) -> float:
    """
    Return the duration a block of `events` will occupy.

    Parameters
    ----------
    events
        The block's events.  May be empty.
    explicit
        An explicit block duration.  When given it is returned as-is; pulseq treats it
        as a floor and raises if the events are longer, which is the semantics wanted.

    Returns
    -------
    float
        Duration in seconds.  Zero for an empty block with no explicit duration.
    """
    if explicit is not None:
        return float(explicit)
    if not events:
        return 0.0
    return float(pp.calc_duration(*events))


def kinds_of(events: Iterable[Event]) -> frozenset[str]:
    """Return the set of ``type`` strings present in `events`."""
    return frozenset(getattr(e, 'type', '?') for e in events)


def channels_of(events: Iterable[Event]) -> dict[str, Event]:
    """
    Return ``channel -> event`` for the gradient events in `events`.

    Raises
    ------
    ValueError
        If two gradients share a channel, which pulseq cannot represent in one block.
    """
    out: dict[str, Event] = {}
    for e in events:
        if getattr(e, 'type', None) in GRAD_TYPES:
            ch = e.channel
            if ch in out:
                msg = f'two gradients on channel {ch!r} in one block'
                raise ValueError(msg)
            out[ch] = e
    return out


def waveform_of(event: Event, raster: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Sample a gradient event onto a uniform raster.

    Parameters
    ----------
    event
        A ``trap`` or ``grad`` event.
    raster
        Gradient raster in seconds.

    Returns
    -------
    t, g
        Times in seconds from the start of the block, and amplitudes in Hz/m.  Both
        include the event's own ``delay`` as leading zeros so that events on different
        channels can be summed directly.
    """
    kind = getattr(event, 'type', None)
    if kind == 'trap':
        rise, flat, fall = event.rise_time, event.flat_time, event.fall_time
        corners_t = np.array([0.0, rise, rise + flat, rise + flat + fall])
        corners_g = np.array([0.0, event.amplitude, event.amplitude, 0.0])
        n = int(round((rise + flat + fall) / raster))
        t = np.arange(n + 1) * raster
        g = np.interp(t, corners_t, corners_g)
    elif kind == 'grad':
        t = np.asarray(event.tt, dtype=float)
        g = np.asarray(event.waveform, dtype=float)
    else:
        msg = f'not a gradient event: type={kind!r}'
        raise ValueError(msg)

    delay = float(getattr(event, 'delay', 0.0) or 0.0)
    if delay > 0:
        n_pad = int(round(delay / raster))
        t = np.concatenate([np.arange(n_pad) * raster, t + delay])
        g = np.concatenate([np.zeros(n_pad), g])
    return t, g


def moment_of(event: Event, raster: float, order: int = 0) -> float:
    """
    Return the `order`-th gradient moment of a single event.

    Parameters
    ----------
    event
        A ``trap`` or ``grad`` event.
    raster
        Gradient raster in seconds, used for trapezoid sampling.
    order
        ``0`` for m0 (area, 1/m), ``1`` for m1 (s/m), ``2`` for m2 (s^2/m).  Moments are
        referenced to the start of the event's block.

    Returns
    -------
    float
        The moment, in pulseq units (amplitudes are Hz/m, so m0 is already 1/m, i.e.
        k-space units, and no gamma appears anywhere).

    Notes
    -----
    ``trap`` events are integrated analytically where possible via ``event.area``; higher
    orders and arbitrary waveforms use the trapezoidal rule on the raster, which is exact
    for piecewise-linear gradients when `order` is 0 and accurate to O(raster^2) above.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=170, slew_unit='T/m/s')
    >>> g = pp.make_trapezoid(channel='x', area=100.0, system=o)
    >>> round(moment_of(g, o.grad_raster_time, 0), 6)
    100.0
    """
    if order == 0 and getattr(event, 'type', None) == 'trap':
        return float(event.area)
    t, g = waveform_of(event, raster)
    return float(trapz(g * t**order, t)) if order else float(trapz(g, t))


def content_hash(event: Event) -> str:
    """
    Return a stable hash of an event's numeric content.

    Used by the purity tests to detect in-place mutation, and by the self-golden digests.
    Registration state (``id``, shape ids, the pypulseq cache) is excluded so that
    registering an event does not change its hash.
    """
    import hashlib  # noqa: PLC0415  (local: keeps import seqcraft light)

    h = hashlib.sha256()
    for key in sorted(vars(event)):
        if key in _STRIP:
            continue
        value = getattr(event, key)
        h.update(key.encode())
        if isinstance(value, np.ndarray):
            h.update(np.ascontiguousarray(value, dtype=np.float64).tobytes())
        else:
            h.update(repr(value).encode())
    return h.hexdigest()


def check_limits(
    events: Iterable[Event],
    opts: Opts,
    raster: float,
) -> list[tuple[str, str, float, float]]:
    """
    Check gradient amplitude and slew against `opts`, per axis and in vector norm.

    pypulseq checks neither: ``Sequence.check_timing`` validates raster, duration and
    continuity only, and ``rotate_3d`` performs no limit check at all.  So this is
    seqcraft's job.

    Two different invariants are reported, and the distinction matters:

    ``'grad'`` / ``'slew'`` on a channel
        The **per-axis** limit, which is what a gradient amplifier actually enforces.  These
        are hard violations.

    ``'grad_norm'`` / ``'slew_norm'``
        The **vector-norm** limit across simultaneous axes.  Two axes ramping together reach
        ``sqrt(2)`` times the per-axis slew in vector magnitude, and three reach ``sqrt(3)``
        -- so a prewinder on x concurrent with a phase-encode blip on y exceeds it routinely
        while being perfectly legal on the hardware.  It is therefore reported as a *warning*
        by default.

        It becomes a hard constraint as soon as the sequence itself applies a rotation,
        because a rotation can concentrate the whole vector onto one physical axis.  That is
        the general form of a real bug in the reference implementation: it never normalised
        its diffusion direction vectors, so ``[1, 1, 0]`` silently requested
        ``sqrt(2) * Gmax`` and ``[1, 1, 1]`` requested ``sqrt(3) * Gmax``.

    Returns
    -------
    list of (kind, where, achieved, limit)
        One entry per violation; empty when everything is within limits.  ``where`` is a
        channel name for per-axis entries and ``'norm'`` for the vector ones.
    """
    events = [e for e in events if getattr(e, 'type', None) in GRAD_TYPES]
    if not events:
        return []

    out: list[tuple[str, str, float, float]] = []
    sampled: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for e in events:
        t, g = waveform_of(e, raster)
        sampled[e.channel] = (t, g)
        peak = float(np.max(np.abs(g))) if g.size else 0.0
        if peak > opts.max_grad * (1 + 1e-6):
            out.append(('grad', e.channel, peak, float(opts.max_grad)))
        if g.size > 1:
            slew = float(np.max(np.abs(np.diff(g) / np.diff(t))))
            if slew > opts.max_slew * (1 + 1e-6):
                out.append(('slew', e.channel, slew, float(opts.max_slew)))

    if len(sampled) > 1:
        n = max(len(g) for _, g in sampled.values())
        stack = np.zeros((len(sampled), n))
        for i, (t, g) in enumerate(sampled.values()):
            stack[i, : len(g)] = g
        norm = np.linalg.norm(stack, axis=0)
        peak = float(np.max(norm))
        if peak > opts.max_grad * (1 + 1e-6):
            out.append(('grad_norm', 'norm', peak, float(opts.max_grad)))
        if n > 1:
            slew = float(np.max(np.abs(np.diff(norm)) / raster))
            if slew > opts.max_slew * (1 + 1e-6):
                out.append(('slew_norm', 'norm', slew, float(opts.max_slew)))
    return out
