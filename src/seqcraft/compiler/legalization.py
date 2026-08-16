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

Limits are measured here too (:func:`check_limits`), on the summed waveform, because that is the
only place the truth is visible: two individually legal gradients on one axis can sum to an
illegal one, and no module can see that in isolation.  A per-axis violation **raises**: there is
no legal sequence to return, so returning one with a note attached would only be a way of not
noticing.
"""

from __future__ import annotations

import functools
import math
from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ..design import events as ev
from ..design import units
from ..design.timing import EPS, TICKS_PER_SECOND, to_ticks
from ..errors import format_error
from .errors import HardwareLimitError
from .model import PlacedEvent, in_block_delay

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from types import SimpleNamespace

    from pypulseq.opts import Opts

__all__ = ['axis_gradient', 'check_limits', 'superpose']


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

    grid_ps = np.array(sorted(marks), dtype=np.int64)
    grid = grid_ps / TICKS_PER_SECOND
    total = np.zeros(len(grid_ps), dtype=float)
    for t_ps, amps in knots:
        # left/right zero: outside its own span a gradient contributes nothing.  At its first and
        # last knot np.interp returns the knot value, so `first` and `last` are carried exactly.
        total += np.interp(grid, t_ps / TICKS_PER_SECOND, amps, left=0.0, right=0.0)
    return grid_ps - a_ps, total


def _as_arbitrary(
    rel_ps: np.ndarray,
    amps: np.ndarray,
    raster_ps: int,
) -> tuple[np.ndarray, float, float] | None:
    """
    Return ``(waveform, first, last)`` if these knots are exactly an arbitrary gradient's.

    A gradient sampled at raster centres has knots at ``0, r/2, 3r/2, ..., (n - 1/2) r, n r`` --
    the two edges carrying ``first`` and ``last``, the interior points being the samples.  Those
    times are **not** on the gradient raster, so :func:`pypulseq.make_extended_trapezoid` rejects
    them; recognising the pattern is what lets a split spiral stay a spiral instead of being
    resampled onto raster edges.

    Splitting works because a block boundary is on the raster: cutting there leaves each piece's
    samples at the centres of *its own* raster intervals, with the seam amplitude becoming one
    piece's ``last`` and the other's ``first``.
    """
    n = len(rel_ps) - 2
    if n < 1 or raster_ps % 2:
        return None
    if int(rel_ps[0]) != 0 or int(rel_ps[-1]) != n * raster_ps:
        return None
    want = (np.arange(n, dtype=np.int64) * 2 + 1) * (raster_ps // 2)
    if not np.array_equal(rel_ps[1:-1], want):
        return None
    return amps[1:-1], float(amps[0]), float(amps[-1])


def axis_gradient(
    axis: str,
    pieces: Sequence[PlacedEvent],
    a: float,
    b: float,
    opts: Opts,
    notes: dict[str, list[str]],
    block_index: int,
) -> SimpleNamespace | None:
    """
    Return the single gradient event for `axis` over the interval ``[a, b)``, or ``None``.

    Fast path -- one gradient, entirely inside the interval -- passes through untouched apart
    from its in-block delay.  No new shape, so a spiral readout or a long diffusion lobe costs
    exactly what it did before compilation.

    Otherwise the pieces are superposed exactly (:func:`superpose`) and the result is emitted in
    whichever pulseq representation can hold those knots without moving them: an extended
    trapezoid when they are all on the gradient raster, an arbitrary gradient when they are the
    raster-centre pattern.  Only a waveform needing *both* -- a trapezoid summed with a spiral --
    can be held by neither, and that one is resampled and reported rather than done quietly.
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
    negligible = 1e-6 * float(opts.max_grad)
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
        return grad

    centre = _as_arbitrary(rel_ps, amps, raster_ps)
    if centre is not None:
        waveform, first, last = centre
        # Limits are measured on the finished block by check_limits.  Letting make_arbitrary_grad
        # raise here would report the violation before the waveform is complete, on a piece rather
        # than on the sum -- which is the one place the truth is visible.
        return pp.make_arbitrary_grad(
            channel=axis,
            waveform=waveform,
            first=first,
            last=last,
            delay=0.0,
            max_grad=math.inf,
            max_slew=math.inf,
            system=opts,
        )

    return _resampled(axis, rel_ps, amps, opts, notes, block_index, pieces, a)


def _resampled(
    axis: str,
    rel_ps: np.ndarray,
    amps: np.ndarray,
    opts: Opts,
    notes: dict[str, list[str]],
    block_index: int,
    pieces: Sequence[PlacedEvent],
    at: float,
) -> SimpleNamespace:
    """
    Force knots onto the gradient raster, and say by how much the waveform moved.

    The one case pulseq cannot represent exactly: a trapezoid's corners are on raster edges, an
    arbitrary gradient's samples are at raster centres, and their sum bends at both.  Neither
    event type has room for that, so something has to give -- but silently giving is how a 2.5 %
    amplitude error reaches a scanner, so the error is measured and reported.

    The bound is exact rather than estimated: two PWL functions differ most at a knot of one of
    them, so comparing at the union of both knot sets *is* the supremum.
    """
    raster_ps = to_ticks(float(opts.grad_raster_time))
    n = int(rel_ps[-1] // raster_ps)
    grid_ps = np.arange(n + 1, dtype=np.int64) * raster_ps
    grid_amps = np.interp(grid_ps / TICKS_PER_SECOND, rel_ps / TICKS_PER_SECOND, amps)

    union = np.union1d(rel_ps, grid_ps) / TICKS_PER_SECOND
    before = np.interp(union, rel_ps / TICKS_PER_SECOND, amps)
    after = np.interp(union, grid_ps / TICKS_PER_SECOND, grid_amps)
    error = float(np.max(np.abs(before - after)))
    notes.setdefault('resample', []).append(
        f'{"+".join(sorted({p.where for p in pieces}))} (axis {axis}) at {at * 1e3:.3f} ms, '
        f'by at most {error / float(opts.max_grad) * 100:.2f} % of max_grad'
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
    return grad


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
    object carrying a note is a way of not noticing: it writes a ``.seq`` the console refuses an
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
