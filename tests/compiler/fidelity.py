"""
An exact oracle for "does the compiled sequence play what the tree said?".

Every pulseq gradient is a **piecewise-linear** function of time, so two of them are equal
everywhere if and only if they agree at the union of their knots.  That makes this a proof
rather than a sampled approximation: no tolerance on the time axis, and no risk of a
distortion hiding between sample points.

Deliberately independent of :mod:`seqcraft.compiler` and
:mod:`seqcraft.design.events` -- it re-derives the pulseq semantics from the raw event fields.
An oracle that shared code with the thing it checks would hide exactly the bugs it exists to
find.

The pulseq semantics being asserted, which is the whole subtlety
-----------------------------------------------------------------
``trap``
    Knots at ``delay``, ``+rise_time``, ``+flat_time``, ``+fall_time`` with amplitudes
    ``0, A, A, 0``.  A triangle (``flat_time == 0``) collapses the middle pair.

``grad``
    Knots at ``delay + tt`` with ``waveform``.  Plus the two edge values pulseq keeps
    *outside* the sample array: ``first`` at ``delay`` and ``last`` at
    ``delay + shape_dur``.  For an extended trapezoid ``tt[0] == 0`` and
    ``tt[-1] == shape_dur`` so those coincide with existing knots and are dropped; for an
    arbitrary gradient ``tt`` sits at raster *centres* -- ``(k + 0.5) * raster`` -- so
    ``first`` and ``last`` are genuinely extra knots half a raster outside the samples.
    Getting this wrong is what makes an arbitrary gradient appear to lose amplitude.

Outside its own span a gradient contributes exactly zero.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from seqcraft.design.logic import BARRIER, LogicBlock, flatten

if TYPE_CHECKING:
    from types import SimpleNamespace

    import pypulseq as pp

__all__ = ['axis_knots', 'compare', 'compiled_knots', 'describe_mismatch', 'tree_knots']

_AXES = ('x', 'y', 'z')
#: Time quantum for knot identity, seconds.  One picosecond: finer than any pulseq raster by
#: five orders of magnitude, and exactly representable as an integer in float64.
_PS = 1e-12


def _q(t: float) -> int:
    """Quantise an absolute time to integer picoseconds, so knots compare exactly."""
    return round(float(t) / _PS)


def event_knots(event: SimpleNamespace, t0: float) -> tuple[list[int], list[float]]:
    """
    Return the exact PWL knots of one gradient event, in absolute time.

    Parameters
    ----------
    event
        A ``trap`` or ``grad`` pulseq event.
    t0
        Absolute time of the event's node, *before* its own ``delay``.

    Returns
    -------
    times_ps, amps
        Knot times in integer picoseconds and amplitudes in Hz/m, strictly ascending.
    """
    kind = getattr(event, 'type', None)
    delay = float(getattr(event, 'delay', 0.0) or 0.0)
    start = t0 + delay

    if kind == 'trap':
        rise = float(event.rise_time)
        flat = float(event.flat_time)
        fall = float(event.fall_time)
        amp = float(event.amplitude)
        if flat > 0:
            ts = [start, start + rise, start + rise + flat, start + rise + flat + fall]
            gs = [0.0, amp, amp, 0.0]
        else:
            ts = [start, start + rise, start + rise + fall]
            gs = [0.0, amp, 0.0]
    elif kind == 'grad':
        tt = np.asarray(event.tt, dtype=float)
        wf = np.asarray(event.waveform, dtype=float)
        shape_dur = float(getattr(event, 'shape_dur', tt[-1]))
        ts = list(start + tt)
        gs = list(wf)
        # The two edge values pulseq stores outside the sample array.  Only real knots when
        # the samples do not already reach the edges, which is the raster-centre case.
        first = getattr(event, 'first', None)
        last = getattr(event, 'last', None)
        if first is not None and tt[0] > _PS:
            ts.insert(0, start)
            gs.insert(0, float(first))
        if last is not None and shape_dur - tt[-1] > _PS:
            ts.append(start + shape_dur)
            gs.append(float(last))
    else:
        msg = f'not a gradient event: type={kind!r}'
        raise ValueError(msg)

    return [_q(t) for t in ts], [float(g) for g in gs]


def _eval(times_ps: list[int], amps: list[float], at: np.ndarray) -> np.ndarray:
    """Evaluate one PWL segment set at `at` (integer ps), zero outside its own span."""
    xs = np.asarray(times_ps, dtype=np.float64)
    ys = np.asarray(amps, dtype=np.float64)
    out = np.interp(at.astype(np.float64), xs, ys, left=0.0, right=0.0)
    out[(at < times_ps[0]) | (at > times_ps[-1])] = 0.0
    return out


def axis_knots(pieces: list[tuple[list[int], list[float]]]) -> np.ndarray:
    """Return the sorted union of knot times over `pieces`, as integer ps."""
    if not pieces:
        return np.zeros(0, dtype=np.int64)
    return np.unique(np.concatenate([np.asarray(t, dtype=np.int64) for t, _ in pieces]))


def tree_knots(root: LogicBlock) -> dict[str, list[tuple[list[int], list[float]]]]:
    """Collect the PWL pieces of every gradient in the tree, keyed by axis."""
    out: dict[str, list[tuple[list[int], list[float]]]] = {ax: [] for ax in _AXES}
    for t0, event, _path in flatten(root):
        kind = getattr(event, 'type', None)
        if kind in ('trap', 'grad') and kind != BARRIER:
            out[str(event.channel)].append(event_knots(event, t0))
    return out


def compiled_knots(seq: pp.Sequence) -> dict[str, tuple[list[int], list[float]]]:
    """
    Reconstruct the single PWL waveform each axis actually plays, over the whole sequence.

    **Blocks are concatenated, not summed.**  A pulseq block holds at most one gradient per
    axis and blocks play in sequence, so the played waveform is one continuous function -- the
    blocks' contributions laid end to end.  Summing them instead double-counts every seam by
    the amplitude at that seam, which looks exactly like a real waveform error.

    Where a block carries no gradient on an axis, or the gradient does not fill the block, the
    axis is zero for the remainder -- and pulseq guarantees a gradient that stops short of a
    block edge reaches zero there.
    """
    out: dict[str, tuple[list[int], list[float]]] = {}
    per_axis: dict[str, list[tuple[int, float]]] = {ax: [] for ax in _AXES}

    t = 0.0
    for index in sorted(seq.block_events):
        block = seq.get_block(index)
        dur = float(seq.block_durations[index])
        for ax in _AXES:
            grad = getattr(block, f'g{ax}', None)
            acc = per_axis[ax]
            if grad is None:
                acc.append((_q(t), 0.0))
                acc.append((_q(t + dur), 0.0))
                continue
            ts, gs = event_knots(grad, t)
            if ts[0] > _q(t):
                acc.append((_q(t), 0.0))
            acc.extend(zip(ts, gs, strict=True))
            if ts[-1] < _q(t + dur):
                acc.append((_q(t + dur), 0.0))
        t += dur

    for ax in _AXES:
        acc = per_axis[ax]
        if not acc:
            out[ax] = ([], [])
            continue
        # Collapse duplicate knot times at block seams.  The two sides must agree there --
        # pulseq's own continuity rule -- so keeping the first is lossless; a disagreement is
        # a real defect and is reported by compare().
        times: list[int] = []
        amps: list[float] = []
        for tq, g in acc:
            if times and tq == times[-1]:
                continue
            times.append(tq)
            amps.append(g)
        out[ax] = (times, amps)
    return out


def seam_discontinuities(
    seq: pp.Sequence,
    tol: float = 1e-6,
) -> list[tuple[float, str, float, float]]:
    """
    Return every block seam where an axis's amplitude jumps.

    pypulseq enforces this when building, but checking it independently is cheap and it is the
    invariant a bad split would break.
    """
    bad: list[tuple[float, str, float, float]] = []
    t = 0.0
    tail: dict[str, float] = dict.fromkeys(_AXES, 0.0)
    for index in sorted(seq.block_events):
        block = seq.get_block(index)
        dur = float(seq.block_durations[index])
        for ax in _AXES:
            grad = getattr(block, f'g{ax}', None)
            if grad is None:
                head, new_tail = 0.0, 0.0
            else:
                ts, gs = event_knots(grad, t)
                head = gs[0] if ts[0] <= _q(t) else 0.0
                new_tail = gs[-1] if ts[-1] >= _q(t + dur) else 0.0
            if abs(head - tail[ax]) > tol:
                bad.append((t, ax, tail[ax], head))
            tail[ax] = new_tail
        t += dur
    # The sequence must also end at zero on every axis.
    bad.extend((t, ax, tail[ax], 0.0) for ax in _AXES if abs(tail[ax]) > tol)
    return bad


def compare(
    root: LogicBlock,
    compiled: pp.Sequence,
    *,
    atol: float = 1e-6,
    rtol: float = 1e-9,
) -> dict[str, dict[str, float]]:
    """
    Compare the tree's gradient waveform with the compiled one, exactly.

    Both sides are summed per axis and evaluated at the **union of every knot on either
    side**, which for piecewise-linear functions is a complete test: agreement there implies
    agreement everywhere.

    Parameters
    ----------
    root
        The tree that was compiled.
    compiled
        What :func:`seqcraft.compile` returned for it.
    atol
        Absolute amplitude tolerance in Hz/m.  1e-6 Hz/m is about 1e-13 of a 40 mT/m
        amplifier's range -- float noise, not physics.
    rtol
        Relative tolerance, scaled by the peak amplitude on that axis.

    Returns
    -------
    dict
        Per axis: ``max_abs_error`` (Hz/m), ``peak`` (Hz/m), ``n_points`` compared, and
        ``t_worst`` (seconds) where the largest error occurred.
    """
    want_all = tree_knots(root)
    got_all = compiled_knots(compiled)
    report: dict[str, dict[str, float]] = {}

    for ax in _AXES:
        want = want_all[ax]
        got_t, got_g = got_all.get(ax, ([], []))
        if not want and not got_t:
            continue
        at = np.unique(
            np.concatenate([axis_knots(want), np.asarray(got_t, dtype=np.int64)])
        )
        # Also probe the midpoint of every knot interval.  Redundant for correct PWL sums but
        # it catches an implementation that agrees at knots while inserting a spurious extra
        # one between them.
        if len(at) > 1:
            mids = ((at[:-1] + at[1:]) // 2).astype(np.int64)
            at = np.unique(np.concatenate([at, mids]))

        # The tree superposes: overlapping events on one axis genuinely add.
        w = np.zeros(len(at))
        for t_ps, amps in want:
            w += _eval(t_ps, amps, at)
        # The compiled side is one concatenated function, already assembled.
        g = _eval(got_t, got_g, at) if got_t else np.zeros(len(at))

        err = np.abs(g - w)
        worst = int(np.argmax(err)) if len(err) else 0
        peak = float(np.max(np.abs(w))) if len(w) else 0.0
        report[ax] = {
            'max_abs_error': float(err[worst]) if len(err) else 0.0,
            'peak': peak,
            'n_points': float(len(at)),
            't_worst': float(at[worst]) * _PS if len(at) else 0.0,
            'tolerance': atol + rtol * peak,
        }
    return report


def describe_mismatch(report: dict[str, dict[str, float]]) -> str:
    """Render :func:`compare`'s output as a one-line-per-axis failure message."""
    lines = []
    for ax, r in sorted(report.items()):
        ok = r['max_abs_error'] <= r['tolerance']
        lines.append(
            f'  axis {ax}: max error {r["max_abs_error"]:.6g} Hz/m at '
            f'{r["t_worst"] * 1e6:.3f} us  (peak {r["peak"]:.6g}, tol {r["tolerance"]:.3g}, '
            f'{int(r["n_points"])} points)  {"OK" if ok else "<-- MISMATCH"}'
        )
    return '\n'.join(lines)


def assert_matches(
    root: LogicBlock,
    compiled: pp.Sequence,
    *,
    atol: float = 1e-6,
    rtol: float = 1e-9,
) -> dict[str, dict[str, float]]:
    """Assert the compiled gradient waveform equals the tree's, on every axis."""
    report = compare(root, compiled, atol=atol, rtol=rtol)
    bad = {ax: r for ax, r in report.items() if r['max_abs_error'] > r['tolerance']}
    if bad:
        raise AssertionError(
            'compiled gradient waveform differs from the logic block:\n'
            + describe_mismatch(report)
        )
    return report
