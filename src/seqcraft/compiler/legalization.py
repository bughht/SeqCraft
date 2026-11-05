"""
Making the gradients legal: superpose, represent, resample, and measure.

Once :mod:`~seqcraft.compiler.boundaries` has decided where the sequence is cut, each
interval may hold several gradient pieces on one axis -- a rewinder overlapping a prephaser, a
readout crossing a seam.  pulseq allows exactly one gradient per axis per block, so they have to
become one event, and this is the only place in the compiler that changes a waveform.

Everything is exact where it can be
-----------------------------------
Both operations -- split one gradient at a boundary, sum several that share an axis -- are exact
on a piecewise-linear representation, and **every pulseq gradient is piecewise linear**.  The sum
of PWL functions is PWL with knots at the union of theirs, so evaluating each piece on that union
and adding reproduces it with no error at all.  :func:`superpose` does that, in integer ticks so a
knot shared by two pieces stays one knot after a sequence that has been running for minutes.

Where it cannot be, it says so
------------------------------
The result then has to be carried by a pulseq event, and the two available shapes reach different
knot sets: an extended trapezoid's knots must lie on the gradient raster, while an arbitrary
gradient's samples sit at raster *centres* (``make_arbitrary_grad`` sets
``tt = (arange(n) + 0.5) * raster``).  :func:`axis_gradient` picks whichever holds the knots
untouched.  Only a waveform bending both on and off the raster -- a trapezoid summed with a
spiral -- fits neither, and that one is resampled by :func:`_resampled`, which measures the
error it introduced and reports it rather than moving the waveform quietly.  A silent 2.5 %
amplitude loss on a spiral is exactly how that requirement was discovered.

A resample is also bounded (ADR-008).  For gradients sampled at raster centres or on raster edges
within the slew limit, it moves k by at most ``r**2 * max_slew / 4`` per raster-centre stretch,
at every instant, whatever the stretch's duration; m1 then moves by at most twice the block's
length times that.  A resample that measures more is refused, which is what a gradient sampled
every half raster can do.

Limits are measured here too (:func:`check_limits`), on the summed waveform, because that is the
only place the truth is visible: two individually legal gradients on one axis can sum to an
illegal one, and no module can see that in isolation.  A per-axis violation **raises**: there is
no legal sequence to return, so returning one with a note attached would only be a way of not
noticing.
"""

from __future__ import annotations

import functools
import math
from types import SimpleNamespace
from typing import TYPE_CHECKING, cast

import numpy as np
import pypulseq as pp

from ..design import events as ev
from ..design import units
from ..design.events import GRADIENT_KINDS, POINT_KINDS
from ..design.logic import BARRIER
from ..design.timing import EPS, TICKS_PER_SECOND, to_ticks
from ..errors import format_error
from .errors import CompileError, HardwareLimitError
from .model import (
    EXCLUSIVE_KINDS,
    LegalizationResult,
    PlacedEvent,
    PulseqReadyBlock,
    in_block_delay,
)
from .verification import require_valid_contract, verify_ready_blocks

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping, Sequence

    from pypulseq.opts import Opts

    from ..design.timing import Raster

__all__ = [
    'axis_gradient',
    'check_limits',
    'common_path',
    'legalize_blocks',
    'required_duration',
    'superpose',
]

#: A superposed block's edge value below this fraction of ``max_grad`` is set to exactly zero (see
#: :func:`axis_gradient`).  The against-tree check models the change that makes.
_NEGLIGIBLE_EDGE = 1e-6

#: PyPulseq stores a raster-centre shape of four samples or fewer without its last half raster, and
#: registration refuses that (ADR-007).  A delayed or zero-padded run that short is resampled
#: instead of being emitted as an arbitrary gradient it would then be refused for.
_SHORTEST_STORED_RUN = 5

#: Relative float allowance on the resample bound: the measured k change is a sum of products of
#: float times and amplitudes, and a waveform built exactly at the slew limit reaches the bound.
_BOUND_ROUNDING = 1e-9


def superpose(
    pieces: Sequence[PlacedEvent],
    a: float,
    b: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return the exact superposition of `pieces` restricted to ``[a, b]``, as PWL knots.

    Both things the compiler does to a gradient -- split one at a block boundary, sum several
    that share an axis -- are *exact* on a piecewise-linear representation, and every pulseq
    gradient is piecewise linear.  The sum of PWL functions is PWL with knots at the union of
    theirs, so evaluating each piece on that union and adding reproduces it with no error.

    This replaces sampling onto a uniform raster grid, which is exact only when every knot
    happens to land on it.  An arbitrary gradient's samples sit at raster *centres*
    (``make_arbitrary_grad`` sets ``tt = (arange(n) + 0.5) * raster``), so they never do, and a
    uniform resample rounded the peaks off a spiral by 2.5 % while preserving its area exactly --
    which is why the m0 invariant could not see it.

    Times are carried as integer ticks so that a knot shared by two pieces is *the same* knot
    after a sequence that has been running for minutes, rather than two knots a femtosecond apart
    that then both survive into the shape.

    Returns
    -------
    rel_ticks, amps
        Knot times in integer ticks from `a`, and amplitudes in Hz/m.
    """
    a_ps, b_ps = to_ticks(a), to_ticks(b)
    marks = {a_ps, b_ps}
    knots: list[tuple[np.ndarray, np.ndarray]] = []
    for p in pieces:
        times, amps = ev.knots_of(p.event, p.node_t)
        t_ps = np.array([to_ticks(t) for t in times], dtype=np.int64)
        marks.update(int(t) for t in t_ps if a_ps < t < b_ps)
        knots.append((t_ps, amps))

    rel_ps = np.array(sorted(marks), dtype=np.int64) - a_ps
    # Interpolated in ticks from the block's start: differences of exact integers, so a value at a
    # cut or at another piece's knot is the same wherever the block sits.  At absolute float
    # seconds its error grew with the block's position in the sequence.
    grid = rel_ps.astype(float)
    total = np.zeros(len(rel_ps), dtype=float)
    for t_ps, amps in knots:
        # left/right zero: outside its own span a gradient contributes nothing.  At its first and
        # last knot np.interp returns the knot value, so `first` and `last` are carried exactly.
        total += np.interp(grid, (t_ps - a_ps).astype(float), amps, left=0.0, right=0.0)
    return rel_ps, total


def _as_arbitrary(
    rel_ps: np.ndarray,
    amps: np.ndarray,
    raster_ps: int,
) -> tuple[np.ndarray, float, float, int] | None:
    """
    Return ``(waveform, first, last, delay_ps)`` if these knots are exactly an arbitrary gradient's.

    A gradient sampled at raster centres has knots at ``d, d + r/2, ..., d + (n - 1/2) r, d + n r``
    -- the two edges carrying ``first`` and ``last``, the interior points being the samples.  Those
    times are **not** on the gradient raster, so :func:`pypulseq.make_extended_trapezoid` rejects
    them; recognising the pattern is what lets a split spiral stay a spiral instead of being
    resampled onto raster edges.

    Splitting works because a block boundary is on the raster: cutting there leaves each piece's
    samples at the centres of *its own* raster intervals, with the seam amplitude becoming one
    piece's ``last`` and the other's ``first``.

    The run normally starts at the block's start (``d = 0``).  It may also start later, at a raster
    edge where it and everything before it are zero, which is a raster-centre gradient that begins
    inside its block and is then cut.  It may end before the block's end, where it and everything
    after are zero.  Both are the same event with a delay, so they are held exactly instead of
    resampled.  A padded run of four samples or fewer is left to the resample, since PyPulseq
    would store it without its last half raster.
    """
    if raster_ps % 2:
        return None
    half = raster_ps // 2
    centres = np.flatnonzero(rel_ps % raster_ps == half)
    n = centres.size
    if n < 1 or centres[-1] - centres[0] + 1 != n:
        return None
    lo, hi = int(centres[0]), int(centres[-1])
    start = int(rel_ps[lo]) - half
    want = start + (np.arange(n, dtype=np.int64) * 2 + 1) * half
    if not np.array_equal(rel_ps[lo:hi + 1], want):
        return None
    end = start + n * raster_ps
    lead, trail = rel_ps[:lo], rel_ps[hi + 1:]
    if start == 0:
        if lead.size != 1:
            return None
    elif lead.size != 2 or int(lead[1]) != start or amps[0] != 0.0 or amps[1] != 0.0:
        return None
    if trail.size == 0 or int(trail[0]) != end:
        return None
    if trail.size > 2 or (trail.size == 2 and (amps[hi + 1] != 0.0 or amps[hi + 2] != 0.0)):
        return None
    if (start > 0 or trail.size == 2) and n < _SHORTEST_STORED_RUN:
        return None
    return amps[lo:hi + 1], float(amps[lo - 1]), float(amps[hi + 1]), start


def axis_gradient(
    axis: str,
    pieces: Sequence[PlacedEvent],
    a: float,
    b: float,
    opts: Opts,
    notes: dict[str, list[str]],
    block_index: int,
    *,
    resampled: list[tuple[int, str]] | None = None,
    moved: list[tuple[str, float, float, float]] | None = None,
) -> SimpleNamespace | None:
    """
    Return the single gradient event for `axis` over the interval ``[a, b)``, or ``None``.

    Fast path -- one gradient, entirely inside the interval -- passes through untouched apart
    from its in-block delay.  No new shape, so a spiral readout or a long diffusion lobe costs
    exactly what it did before compilation.

    Otherwise the pieces are superposed exactly (:func:`superpose`) and the result is emitted in
    whichever pulseq representation can hold those knots without moving them: an extended
    trapezoid when they are all on the gradient raster, an arbitrary gradient when they are the
    raster-centre pattern, delayed or zero-padded within the block.  Only a waveform needing
    *both* -- a trapezoid summed with a spiral -- can be held by neither, and that one is resampled
    within a bound and reported rather than done quietly: in the warning notes, as
    ``(block_index, axis)`` in `resampled`, and as ``(axis, sup |Δk|, Δm0, Δm1)`` in `moved`, when
    those are given.
    """
    if not pieces:
        return None

    if len(pieces) == 1:
        p = pieces[0]
        if p.start >= a - EPS and p.end <= b + EPS:
            return ev.derive(p.event, delay=in_block_delay(p, a, opts))

    if len(pieces) > 1:
        notes.setdefault('merge', []).append(
            '+'.join(sorted({p.where for p in pieces})) + f' (axis {axis})'
        )

    rel_ps, amps = superpose(pieces, a, b)
    if not np.any(amps):
        return None

    # Snap negligible endpoints to exactly zero.  Two pieces that cancel at a seam leave ~1e-12 of
    # the amplifier's range rather than 0.0, and pypulseq tests ``first != 0`` *exactly* before
    # demanding the previous block continue it -- so without this a rounding artifact becomes a
    # continuity error hundreds of blocks away from anything that looks wrong.
    negligible = _NEGLIGIBLE_EDGE * float(opts.max_grad)
    for edge in (0, -1):
        if abs(amps[edge]) < negligible:
            amps[edge] = 0.0

    raster_ps = to_ticks(float(opts.grad_raster_time))
    if np.all(rel_ps % raster_ps == 0):
        # The interval edges are block boundaries, so a non-zero value at either end has to be
        # carried by first/last with delay 0: pypulseq rejects a delayed gradient that starts away
        # from zero, and checks the join against the neighbouring block.
        grad = pp.make_extended_trapezoid(
            channel=axis,
            times=rel_ps / TICKS_PER_SECOND,
            amplitudes=amps,
            system=opts,
            max_grad=math.inf,
            max_slew=math.inf,
            skip_check=True,
        )
        grad.delay = 0.0
        return cast('SimpleNamespace', grad)

    centre = _as_arbitrary(rel_ps, amps, raster_ps)
    if centre is not None:
        waveform, first, last, delay_ps = centre
        # Limits are measured on the finished block by check_limits.  Letting make_arbitrary_grad
        # raise here would report the violation before the waveform is complete, on a piece rather
        # than on the sum -- which is the one place the truth is visible.
        grad = pp.make_arbitrary_grad(
            channel=axis,
            waveform=waveform,
            first=first,
            last=last,
            delay=delay_ps / TICKS_PER_SECOND,
            max_grad=math.inf,
            max_slew=math.inf,
            system=opts,
        )
        return cast('SimpleNamespace', grad)

    if resampled is not None:
        resampled.append((block_index, axis))
    return _resampled(axis, rel_ps, amps, opts, notes, block_index, pieces, a, moved)


def _resampled(
    axis: str,
    rel_ps: np.ndarray,
    amps: np.ndarray,
    opts: Opts,
    notes: dict[str, list[str]],
    block_index: int,
    pieces: Sequence[PlacedEvent],
    at: float,
    moved: list[tuple[str, float, float, float]] | None = None,
) -> SimpleNamespace:
    """
    Force knots onto the gradient raster, refuse it beyond the bound, and say what moved.

    The one case pulseq cannot represent exactly: a trapezoid's corners are on raster edges, an
    arbitrary gradient's samples are at raster centres, and their sum bends at both.  Neither
    event type has room for that, so something has to give -- but silently giving is how a 2.5 %
    amplitude error reaches a scanner, so the error is measured and reported.

    The amplitude bound is exact rather than estimated: two PWL functions differ most at a knot of
    one of them, so comparing at the union of both knot sets *is* the supremum.

    **The moment bound (ADR-008).**  Keeping only the raster-edge values replaces, in each raster
    interval, a raster-centre sample's tent by its chord.  Over a stretch of samples ``c_i`` the
    area changes telescope to ``(r/8)·(D_end - D_start)``, ``D`` being the sample step across each
    end, so a stretch within the slew limit moves the running k by at most ``r²·max_slew/4`` at
    every instant, whatever its duration.  Edge knots are kept exactly and add nothing.  m1 needs
    no bound of its own: by parts, ``Δm1 = L·Δk(L) - ∫Δk dt``, so ``|Δm1| <= 2·L·sup|Δk|`` for a
    block of length ``L``.  The bound is therefore on the running k, not on the block's m0, which
    can cancel while k strays.  Half-raster input has no telescoping, so its resample is
    measured against the same bound, never assumed to meet it.
    """
    raster_ps = to_ticks(float(opts.grad_raster_time))
    n = int(rel_ps[-1] // raster_ps)
    grid_ps = np.arange(n + 1, dtype=np.int64) * raster_ps
    grid_amps = np.interp(grid_ps / TICKS_PER_SECOND, rel_ps / TICKS_PER_SECOND, amps)

    union = np.union1d(rel_ps, grid_ps) / TICKS_PER_SECOND
    before = np.interp(union, rel_ps / TICKS_PER_SECOND, amps)
    after = np.interp(union, grid_ps / TICKS_PER_SECOND, grid_amps)
    error = float(np.max(np.abs(before - after)))
    k_moved, dm0, dm1 = _resample_moments(union, after - before)

    sites = '+'.join(sorted({p.where for p in pieces}))
    r = float(opts.grad_raster_time)
    stretches = sum(1 for p in pieces if not _on_raster_edges(p, raster_ps))
    bound = stretches * r * r * float(opts.max_slew) / 4.0
    if k_moved > bound * (1.0 + _BOUND_ROUNDING):
        msg = format_error(
            f'block {block_index} cannot be resampled onto the raster on axis {axis} within the '
            f'resample bound: it would move k by {k_moved:.4g} 1/m.',
            {
                'from': sites,
                'at': f'{at * 1e3:.3f} ms',
                'k moved by': f'{k_moved:.4g} 1/m',
                'bound': f'{bound:.4g} 1/m = {stretches} x raster^2 x max_slew / 4',
                'block m0 / m1 change': f'{dm0:+.4g} 1/m / {dm1:+.4g} s/m',
            },
            [
                'the bound holds for gradients sampled at raster centres or on raster edges within '
                'the slew limit; a gradient sampled every half raster (oversampling=True) has '
                'none, so give it no same-axis overlap, or a barrier at each end so it keeps its '
                'own block',
                'a piece built past the slew limit has no bound either: build it against the '
                'system limits',
            ],
        )
        raise CompileError(msg)
    if moved is not None:
        moved.append((axis, k_moved, dm0, dm1))
    notes.setdefault('resample', []).append(
        f'{sites} (axis {axis}) at {at * 1e3:.3f} ms, '
        f'by at most {error / float(opts.max_grad) * 100:.2f} % of max_grad, '
        f'k by up to {k_moved:.3g} 1/m of a {bound:.3g} 1/m bound, '
        f'block m0 {dm0:+.3g} 1/m and m1 {dm1:+.3g} s/m'
    )
    grad = pp.make_extended_trapezoid(
        channel=axis,
        times=grid_ps / TICKS_PER_SECOND,
        amplitudes=grid_amps,
        system=opts,
        max_grad=math.inf,
        max_slew=math.inf,
        skip_check=True,
    )
    grad.delay = 0.0
    return cast('SimpleNamespace', grad)


def _resample_moments(t: np.ndarray, d: np.ndarray) -> tuple[float, float, float]:
    """
    A resample's change as ``(sup |Δk|, Δm0, Δm1)``, about the block's start.

    `d` is the resampled waveform minus the exact one at `t`, the knots of both, so it is linear
    between them and the running integral is exact there.  Inside a segment the running integral
    peaks only where `d` crosses zero, so those points are checked as well.
    """
    h = np.diff(t)
    k = np.concatenate(([0.0], np.cumsum(h * (d[:-1] + d[1:]) / 2.0)))
    peak = float(np.max(np.abs(k)))
    cross = np.flatnonzero(d[:-1] * d[1:] < 0.0)
    if cross.size:
        reach = h[cross] * d[cross] / (d[cross] - d[cross + 1])
        peak = max(peak, float(np.max(np.abs(k[cross] + reach * d[cross] / 2.0))))
    return peak, float(k[-1]), ev.pwl_moment(t, d, 1)


def _on_raster_edges(piece: PlacedEvent, raster_ps: int) -> bool:
    """Whether every knot of `piece` is on a raster edge, so that a resample keeps it exactly."""
    times, _ = ev.knots_of(piece.event, piece.node_t)
    return all(to_ticks(float(t)) % raster_ps == 0 for t in times)


def check_limits(
    events: Iterable[SimpleNamespace],
    opts: Opts,
    block_index: int,
    origin: str,
    notes: dict[str, list[str]],
    *,
    start: float = 0.0,
) -> None:
    """
    Measure amplitude and slew on a compiled block; raise on a per-axis violation.

    Raised rather than reported, because there is no legal sequence to hand back.  A returned
    object carrying a note is a way of not noticing: it writes a ``.seq`` whose problem surfaces an
    hour later, and the report explaining why is on an object nobody looked at.

    The vector norm across simultaneous axes is the exception, and stays a **warning**: two axes
    ramping together reach ``sqrt(2)`` times the per-axis slew in vector magnitude, which real
    amplifiers permit.  It becomes a hard limit only once a rotation can concentrate the whole
    vector onto one physical axis -- so it is measured, said out loud, and not fatal.

    Parameters
    ----------
    events
        The finished events of one compiled block.
    opts
        The scanner whose ``max_grad`` and ``max_slew`` are the ceiling.
    block_index, origin, start
        Where this block is, for the message: its index, the tag paths that fed it, and its
        absolute start time.  The time is what a caller can act on -- a block index means
        nothing until the sequence is written.
    notes
        Appended to under ``'norm'`` with the vector-norm findings, which are informational.

    Raises
    ------
    HardwareLimitError
        On the first per-axis amplitude or slew violation.
    """
    raster = float(opts.grad_raster_time)
    gamma = float(opts.gamma)
    for kind, where, achieved, limit in ev.check_limits(events, opts, raster):
        pulseq_unit, unit = ('Hz/m', 'mT/m') if kind.startswith('grad') else ('Hz/m/s', 'T/m/s')
        to_display = functools.partial(
            units.convert, from_unit=pulseq_unit, to_unit=unit, gamma=gamma
        )
        if kind.endswith('_norm'):
            notes.setdefault('norm', []).append(
                f'{origin} ({kind.replace("_norm", "")} {achieved / limit * 100:.0f}% at '
                f'{start * 1e3:.3f} ms)'
            )
            continue
        what = 'gradient' if kind == 'grad' else 'slew'
        knob = 'grad' if kind == 'grad' else 'slew'
        # Floored rather than rounded, so designing against the suggested derating actually
        # clears the limit instead of landing on 100.4 % of it.
        headroom = max(0.01, math.floor(limit / achieved * 100.0) / 100.0)
        remedy = (
            'lengthen the lobe, or lower the readout bandwidth'
            if kind == 'slew'
            else 'lower the amplitude, or lengthen the lobe'
        )
        msg = format_error(
            f'{what} {achieved / limit * 100:.0f}% of the {to_display(limit):.0f} {unit} '
            f'limit on axis {where}.',
            {
                'from': origin,
                'at': f'{start * 1e3:.3f} ms (block {block_index})',
                'reached': f'{to_display(achieved):.1f} {unit}',
            },
            [
                remedy,
                f'or design that part against sc.opts.derate(opts, {knob}={headroom:.2f}) -- '
                f'the finished sequence is still compiled against the full opts',
            ],
        )
        raise HardwareLimitError(msg)


def required_duration(events: Sequence[SimpleNamespace], opts: Opts) -> float:
    """Return the shortest block that can hold `events`, by PyPulseq's own rules."""
    out = 0.0
    for event in events:
        kind = getattr(event, 'type', None)
        delay = float(getattr(event, 'delay', 0.0) or 0.0)
        if kind == 'rf':
            ring = float(getattr(event, 'ringdown_time', opts.rf_ringdown_time) or 0.0)
            out = max(out, delay + float(event.shape_dur) + ring)
        elif kind == 'adc':
            dead = float(getattr(event, 'dead_time', opts.adc_dead_time) or 0.0)
            out = max(out, delay + float(event.num_samples) * float(event.dwell) + dead)
        elif kind != BARRIER:
            out = max(out, float(pp.calc_duration(event)))
    return out


def _adc_conflict(edge: float, adcs: Sequence[PlacedEvent]) -> PlacedEvent | None:
    """Return an ADC whose sampling window contains the proposed boundary `edge`."""
    for adc in adcs:
        if adc.start + EPS < edge < adc.end - EPS:
            return adc
    return None


def common_path(paths: Sequence[tuple[str, ...]]) -> tuple[str, ...]:
    """Return the longest tag path common to all non-empty paths."""
    real = [path for path in paths if path]
    if not real:
        return ()
    common = list(real[0])
    for path in real[1:]:
        keep = 0
        for left, right in zip(common, path):
            if left != right:
                break
            keep += 1
        common = common[:keep]
        if not common:
            break
    return tuple(common)


def legalize_blocks(  # noqa: C901, PLR0912
    edges: Sequence[float],
    placed: Sequence[PlacedEvent],
    targets: Mapping[int, float],
    opts: Opts,
    raster: Raster,
) -> LegalizationResult:
    """
    Assemble and validate every ready block for the fixed boundary and label decisions.

    This is the last policy-bearing stage. It schedules point events, splits and superposes
    gradients, checks the finished waveforms and returns immutable blocks plus explicit notes.
    PyPulseq sequence mutation belongs to :func:`seqcraft.compiler.emission.emit_blocks`.
    """
    notes: dict[str, list[str]] = {}
    blocks: list[PulseqReadyBlock] = []
    resampled: list[tuple[int, str]] = []
    moved: list[tuple[str, float, float, float]] = []

    grads: dict[str, list[PlacedEvent]] = {}
    for event in placed:
        if event.kind in GRADIENT_KINDS:
            grads.setdefault(event.event.channel, []).append(event)
    for pieces in grads.values():
        pieces.sort(key=lambda event: event.start)
    cursor = dict.fromkeys(grads, 0)
    active: dict[str, list[PlacedEvent]] = {axis: [] for axis in grads}

    # Labels are assigned at their target ADC time. The event time remains the secondary key so
    # diagnostics stay intuitive, although label_order_conflict has already rejected groups whose
    # meaning would depend on intra-block insertion order.
    schedule = sorted(
        (
            (targets.get(index, event.res_start), event.res_start, event)
            for index, event in enumerate(placed)
            if event.kind in EXCLUSIVE_KINDS or event.kind in POINT_KINDS
        ),
        key=lambda item: (item[0], item[1]),
    )
    assign_at = [time for time, _, _ in schedule]
    singles = [event for _, _, event in schedule]
    single_cursor = 0
    adcs = [event for event in placed if event.kind == 'adc']

    for index, (start, end) in enumerate(zip(edges, edges[1:])):
        events: list[SimpleNamespace] = []
        paths: list[tuple[str, ...]] = []

        while single_cursor < len(singles) and assign_at[single_cursor] < end - EPS:
            event = singles[single_cursor]
            single_cursor += 1
            if event.kind in POINT_KINDS and not hasattr(event.event, 'delay'):
                events.append(event.event)
            else:
                events.append(ev.derive(event.event, delay=in_block_delay(event, start, opts)))
            paths.append(event.path)

        for axis, pieces in grads.items():
            while cursor[axis] < len(pieces) and pieces[cursor[axis]].start < end - EPS:
                active[axis].append(pieces[cursor[axis]])
                cursor[axis] += 1
            active[axis] = [event for event in active[axis] if event.end > start + EPS]
            here = [event for event in active[axis] if event.start < end - EPS]
            if not here:
                continue
            crossing = [
                event for event in here
                if event.start < start - EPS or event.end > end + EPS
            ]
            for edge in (start, end) if crossing else ():
                blocked = _adc_conflict(edge, adcs)
                if blocked is None:
                    continue
                msg = format_error(
                    f'a block boundary at {edge * 1e6:.1f} us falls inside a gradient on axis '
                    f'{axis} that an ADC is sampling.',
                    {
                        'gradient from': ', '.join(sorted({event.where for event in crossing})),
                        'adc from': blocked.where,
                        'adc window': (
                            f'{blocked.start * 1e6:.1f} .. {blocked.end * 1e6:.1f} us'
                        ),
                    },
                    [
                        'a readout gradient must stay one event -- ramp sampling and vendor '
                        'gridding both depend on that, so the compiler will not split it',
                        'the boundary comes from an explicit barrier, or this is a compiler bug; '
                        'please report it with the tree that produced it',
                    ],
                )
                raise CompileError(msg)
            gradient = axis_gradient(axis, here, start, end, opts, notes, index,
                                     resampled=resampled, moved=moved)
            if gradient is not None:
                events.append(gradient)
                paths.extend(event.path for event in here)

        duration = raster.nearest(end - start)
        ready = PulseqReadyBlock(
            index=index,
            start=start,
            end=end,
            duration=duration,
            events=tuple(events),
            source_paths=tuple(paths),
            origin=common_path(paths),
        )
        origin = ', '.join(sorted({'.'.join(path) for path in ready.source_paths if path})) or '?'
        if ready.events:
            check_limits(ready.events, opts, index, origin, notes, start=start)
            needed = required_duration(ready.events, opts)
            if needed > ready.duration + EPS:
                msg = format_error(
                    f'block {index} spans {ready.duration * 1e6:.1f} us but its events need '
                    f'{needed * 1e6:.1f} us.',
                    {'from': origin, 'shortfall_us': round((needed - ready.duration) * 1e6, 3)},
                    [
                        'usually an ADC whose trailing dead time, or an RF whose ringdown, extends '
                        'past the interval -- the module should report a longer duration',
                        'this is a compiler bug if the module\'s duration property is correct; '
                        'please report it with the tree that produced it',
                    ],
                )
                raise CompileError(msg)
        blocks.append(ready)

    # Each resample is bounded on its own; several on one axis add up.  The net m0 is what a
    # balanced design sees, and the summed magnitudes are what cancellation would hide.  m1 has no
    # signed total: block m1 values share no origin, and about t = 0 the total would bring back the
    # absolute-time lever arm ADR-007 retired.
    for axis in sorted({axis for axis, _, _, _ in moved}):
        rows = [row for row in moved if row[0] == axis]
        if len(rows) > 1:
            notes.setdefault('resample_total', []).append(
                f'axis {axis} over {len(rows)} blocks: m0 {sum(r[2] for r in rows):+.3g} 1/m net; '
                f'k {sum(r[1] for r in rows):.3g} 1/m and m1 {sum(abs(r[3]) for r in rows):.3g} '
                f's/m summed in magnitude'
            )

    ready_blocks = tuple(blocks)
    require_valid_contract(
        'ready-block',
        verify_ready_blocks(
            ready_blocks,
            expected_start=edges[0] if edges else None,
        ),
    )
    frozen_notes = tuple(
        (category, tuple(entries))
        for category, entries in sorted(notes.items())
        if entries
    )
    return LegalizationResult(ready_blocks, frozen_notes, tuple(resampled))
