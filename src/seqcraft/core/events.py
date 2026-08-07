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
import itertools
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any

import numpy as np
import pypulseq as pp

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from pypulseq.opts import Opts

#: Times closer than this are the same knot.  One nanosecond: 100x finer than the finest
#: pulseq raster, and ~1e6x coarser than float64 resolution at millisecond magnitudes.
_EDGE_EPS = 1e-9

__all__ = [
    'Event',
    'GRAD_TYPES',
    'channels_of',
    'content_hash',
    'derive',
    'duration_of',
    'kinds_of',
    'knots_of',
    'moment_of',
    'pwl_moment',
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
    Return a gradient event as a curve to plot or interpolate, with its own time axis.

    A trapezoid is sampled onto `raster`.  A ``grad`` event is returned **at its own sample
    times**, which are the raster *centres* for ``make_arbitrary_grad`` and the *knots* for
    ``make_extended_trapezoid`` -- neither of which lies on `raster`, and the latter of which is
    not uniformly spaced at all.

    So **never assume the returned times are uniform**: use the ``t`` that comes back.  Every
    caller here does (``np.interp(grid, tt, wf)`` in the MRzero bridge, ``trapz(wf * tt, tt)`` in
    the diffusion moment check).  :func:`check_limits` did not, and divided by `raster` where it
    should have divided by ``diff(t)``; on an EPI train's extended trapezoid that overstated the
    slew by 26x.  For arithmetic -- moments, splits, sums -- use :func:`knots_of` and
    :func:`pwl_moment`, which are exact.

    Parameters
    ----------
    event
        A ``trap`` or ``grad`` event.
    raster
        Gradient raster in seconds.  Used for ``trap`` only; a ``grad`` event carries its own
        sample times and they are returned unchanged.

    Returns
    -------
    t, g
        Times in seconds from the start of the block, and amplitudes in Hz/m.  Both
        include the event's own ``delay`` as leading zeros so that events on different
        channels can be summed directly.

    Examples
    --------
    >>> import numpy as np
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=170, slew_unit='T/m/s')
    >>> t, g = waveform_of(pp.make_trapezoid(channel='x', area=100.0, system=o), 1e-5)
    >>> bool(np.allclose(np.diff(t), 1e-5))                 # a trapezoid: uniform
    True
    >>> extended = pp.make_extended_trapezoid(
    ...     'x', amplitudes=np.array([0.0, 1e6, 0.0]),
    ...     times=np.array([0.0, 2.6e-4, 5.2e-4]), system=o)
    >>> t, g = waveform_of(extended, 1e-5)
    >>> [round(v * 1e6) for v in t]                         # its knots, not the 10 us raster
    [0, 260, 520]
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


def knots_of(event: Event, t0: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """
    Return the exact piecewise-linear knots of a gradient event, in absolute time.

    Every pulseq gradient *is* a PWL function, so this loses nothing -- and it is the only
    representation in which a moment, a split or a sum is exact.  :func:`waveform_of` samples
    onto a uniform raster instead, which is what plotting wants and what arithmetic does not:
    an arbitrary gradient's samples sit at raster *centres*, so they never land on the grid.

    The subtlety is ``grad``: its samples sit at raster centres when it came from
    ``make_arbitrary_grad``, and the amplitudes at the two edges live outside the sample array
    in ``first`` and ``last``.  Omitting them truncates the waveform by half a raster at each
    end.  For an extended trapezoid ``tt`` already reaches both edges, so they coincide with
    existing knots and are skipped.

    Parameters
    ----------
    event
        A ``trap`` or ``grad`` event.
    t0
        Absolute time of the event's node.  The event's own ``delay`` is added to it.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=170, slew_unit='T/m/s')
    >>> flat = pp.make_trapezoid(channel='x', flat_area=100.0, flat_time=1e-3, system=o)
    >>> t, a = knots_of(flat)
    >>> len(t), float(a[0]), float(a[-1])           # up, along, down
    (4, 0.0, 0.0)
    >>> t, a = knots_of(pp.make_trapezoid(channel='x', area=100.0, system=o))
    >>> len(t)                                      # no flat top: a triangle has three
    3
    """
    kind = getattr(event, 'type', None)
    start = t0 + float(getattr(event, 'delay', 0.0) or 0.0)
    if kind == 'trap':
        rise = float(event.rise_time)
        flat = float(event.flat_time)
        fall = float(event.fall_time)
        amp = float(event.amplitude)
        if flat > 0:
            return (
                np.array([start, start + rise, start + rise + flat, start + rise + flat + fall]),
                np.array([0.0, amp, amp, 0.0]),
            )
        return (
            np.array([start, start + rise, start + rise + fall]),
            np.array([0.0, amp, 0.0]),
        )
    if kind != 'grad':
        msg = f'not a gradient event: type={kind!r}'
        raise ValueError(msg)

    tt = np.asarray(event.tt, dtype=float)
    wf = np.asarray(event.waveform, dtype=float)
    shape_dur = float(getattr(event, 'shape_dur', tt[-1] if len(tt) else 0.0))
    times = [start + t for t in tt]
    amps = list(wf)
    first = getattr(event, 'first', None)
    last = getattr(event, 'last', None)
    if first is not None and len(tt) and tt[0] > _EDGE_EPS:
        times.insert(0, start)
        amps.insert(0, float(first))
    if last is not None and len(tt) and shape_dur - tt[-1] > _EDGE_EPS:
        times.append(start + shape_dur)
        amps.append(float(last))
    return np.asarray(times, dtype=float), np.asarray(amps, dtype=float)


def pwl_moment(times: np.ndarray, amps: np.ndarray, order: int = 0) -> float:
    """
    Return the exact ``integral g(t) * t**order dt`` of a piecewise-linear waveform.

    Closed form, not quadrature of a sampled curve.  On one segment ``g(t) * t**order`` is a
    polynomial of degree ``order + 1``, and `k`-node Gauss-Legendre is exact to degree
    ``2k - 1`` -- so a handful of nodes integrates it with no discretisation error at all.

    Why this matters rather than being a nicety: the trapezoidal rule on raster samples is
    exact for ``order == 0`` and *only* then.  Using it for m1 made the tree side and the
    compiled side disagree by ~0.1 % on nothing more than a split, because they were sampled
    at different points -- useless as an invariant.  And using it for m0 on both sides made the
    errors cancel, which is worse: the invariant then cannot see a resampling that moved the
    waveform, which is exactly how a 2.5 % amplitude loss stayed hidden.

    Examples
    --------
    >>> import numpy as np
    >>> t = np.array([0.0, 1.0, 2.0])
    >>> g = np.array([0.0, 1.0, 0.0])
    >>> round(pwl_moment(t, g, 0), 12)          # area of a unit triangle
    1.0
    >>> round(pwl_moment(t, g, 1), 12)          # centroid at t = 1
    1.0
    """
    if len(times) < 2:
        return 0.0
    t0, t1 = times[:-1], times[1:]
    g0, g1 = amps[:-1], amps[1:]
    h = t1 - t0
    nodes, weights = np.polynomial.legendre.leggauss(max(2, (order + 3) // 2))
    total = 0.0
    for x, w in zip(nodes, weights):
        s = 0.5 * (x + 1.0)  # map the Gauss nodes from [-1, 1] onto the segment's [0, 1]
        total += float(np.sum(w * 0.5 * h * (g0 + s * (g1 - g0)) * (t0 + s * h) ** order))
    return total


def moment_of(event: Event, raster: float = 0.0, order: int = 0) -> float:
    """
    Return the `order`-th gradient moment of a single event, exactly.

    Parameters
    ----------
    event
        A ``trap`` or ``grad`` event.
    raster
        Unused, and kept only so existing calls keep working.  The moment is computed from the
        event's exact knots, which needs no raster -- passing one never made it more accurate,
        only less.
    order
        ``0`` for m0 (area, 1/m), ``1`` for m1 (s/m), ``2`` for m2 (s^2/m).  Referenced to the
        start of the event's block, so an event's own ``delay`` counts toward orders above 0.

    Returns
    -------
    float
        The moment, in pulseq units (amplitudes are Hz/m, so m0 is already 1/m, i.e.
        k-space units, and no gamma appears anywhere).

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=170, slew_unit='T/m/s')
    >>> g = pp.make_trapezoid(channel='x', area=100.0, system=o)
    >>> round(moment_of(g, o.grad_raster_time, 0), 6)
    100.0
    """
    del raster
    return pwl_moment(*knots_of(event), order)


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


def _merge_times(times: np.ndarray) -> np.ndarray:
    """Return `times` sorted with knots closer than one nanosecond collapsed into one."""
    ordered = np.sort(np.asarray(times, dtype=float))
    if ordered.size < 2:
        return ordered
    keep = np.concatenate([[True], np.diff(ordered) > _EDGE_EPS])
    return ordered[keep]


def _on_grid(pieces: Sequence[tuple[np.ndarray, np.ndarray]]) -> tuple[np.ndarray, np.ndarray]:
    """
    Return `pieces` summed on the union of their knots.

    Exact, because a sum of piecewise-linear functions is piecewise linear with knots at the
    union of theirs -- the same fact :func:`seqcraft.core.compiler._superpose` rests on.
    """
    grid = _merge_times(np.concatenate([t for t, _ in pieces]))
    total = np.zeros(len(grid))
    for t, a in pieces:
        total += np.interp(grid, t, a, left=0.0, right=0.0)
    return grid, total


def check_limits(
    events: Iterable[Event],
    opts: Opts,
    raster: float = 0.0,
    *,
    starts: Sequence[float] | None = None,
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

    Parameters
    ----------
    events
        Gradient events.  Non-gradient events are ignored.
    opts
        The limits to check against.
    raster
        Unused, and kept so existing positional calls work.  Everything here is computed from
        exact knots; see the note below for why it was a liability rather than a parameter.
    starts
        Optional node time of each event, seconds.  Needed only when two events **on the same
        axis** are given: without it they are both taken to play at their own ``delay``, and
        two lobes of a bipolar pair then land on top of each other.  The compiler passes
        compiled-block events, at most one per axis, so it needs nothing here.

    Returns
    -------
    list of (kind, where, achieved, limit)
        One entry per violation; empty when everything is within limits.  ``where`` is a
        channel name for per-axis entries and ``'norm'`` for the vector ones.

    Notes
    -----
    **Everything is measured on the union of the events' knots**, where a piecewise-linear
    function's slope is piecewise constant, so ``diff(a) / diff(t)`` is the slew rather than an
    approximation of it.  Three bugs lived in the sampled version this replaces, all of them
    reachable only through an *extended trapezoid* -- whose knots are neither uniform nor on the
    raster, and which is what an EPI readout train has to be if it is to be one event per axis
    rather than one per echo:

    * the vector-norm slew divided ``diff(norm)`` by the **raster** rather than by ``diff(t)``.
      On one EPI echo junction, whose knots are 260 us apart, that reported **2441 % of the slew
      limit where the truth is 94 %** -- 4882.9 T/m/s being exactly the peak amplitude over one
      10 us raster.  The per-axis path divided by ``diff(t)`` and was right, which is why the
      error surfaced only as a warning.
    * the axes were stacked by sample **index**, so ``gx`` at 260 us was combined with ``gy`` at
      490 us, and the shorter axis was zero-padded at the end instead of held at its ``last``.
      That can hide a real violation as easily as invent one.
    * a second event on an axis **overwrote** the first in the norm's dictionary, so only the
      last was ever combined.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=170, slew_unit='T/m/s')
    >>> gx = pp.make_trapezoid(channel='x', amplitude=0.9 * o.max_grad, duration=1e-3, system=o)
    >>> gy = pp.make_trapezoid(channel='y', amplitude=0.9 * o.max_grad, duration=1e-3, system=o)
    >>> [kind for kind, *_ in check_limits([gx], o)]              # legal on its own
    []
    >>> [kind for kind, *_ in check_limits([gx, gy], o)]          # sqrt(2) in vector norm
    ['grad_norm', 'slew_norm']
    """
    offsets = itertools.repeat(0.0) if starts is None else starts
    chosen = [
        (event, float(t0))
        for event, t0 in zip(events, offsets)
        if getattr(event, 'type', None) in GRAD_TYPES
    ]
    if not chosen:
        return []

    by_channel: dict[str, list[tuple[np.ndarray, np.ndarray]]] = {}
    for event, t0 in chosen:
        by_channel.setdefault(event.channel, []).append(knots_of(event, t0))

    out: list[tuple[str, str, float, float]] = []
    curves: list[tuple[np.ndarray, np.ndarray]] = []
    for channel, pieces in by_channel.items():
        grid, amps = _on_grid(pieces)
        curves.append((grid, amps))
        peak = float(np.max(np.abs(amps))) if amps.size else 0.0
        if peak > opts.max_grad * (1 + 1e-6):
            out.append(('grad', channel, peak, float(opts.max_grad)))
        if grid.size > 1:
            slew = float(np.max(np.abs(np.diff(amps) / np.diff(grid))))
            if slew > opts.max_slew * (1 + 1e-6):
                out.append(('slew', channel, slew, float(opts.max_slew)))

    if len(curves) > 1:
        grid = _merge_times(np.concatenate([t for t, _ in curves]))
        stack = np.stack([np.interp(grid, t, a, left=0.0, right=0.0) for t, a in curves])
        norm = np.linalg.norm(stack, axis=0)
        peak = float(np.max(norm))
        if peak > opts.max_grad * (1 + 1e-6):
            out.append(('grad_norm', 'norm', peak, float(opts.max_grad)))
        if grid.size > 1:
            slew = float(np.max(np.abs(np.diff(norm) / np.diff(grid))))
            if slew > opts.max_slew * (1 + 1e-6):
                out.append(('slew_norm', 'norm', slew, float(opts.max_slew)))
    return out
