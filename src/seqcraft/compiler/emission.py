"""
Ready blocks to a :class:`pypulseq.Sequence`.

The last mechanical stage.  By the time it runs, every decision has been taken: the boundaries are
fixed, the label assignments are fixed, and the gradients are legal.  What is left is to walk the
intervals, collect what falls in each, hand the result to ``add_block``, and record where it came
from.

Three things it checks on the way, all of them shaped by pulseq being permissive where it should
not be:

**A block must be long enough for its contents.**  ``set_block`` takes a block's duration as the
max over its events, so an interval shorter than what it holds silently produces an off-raster
block instead of an error.  :func:`required_duration` recomputes the floor by pulseq's own rules
and the caller stops there, where the boundary that caused it can still be named.

**No boundary may split an ADC's sampling window.**  :func:`_adc_conflict` is a safety net rather
than an expected condition -- boundaries are chosen outside every reservation, and a reservation
contains its ADC's window -- but a silently split readout reaching a scanner is bad enough to be
worth one comparison per crossing gradient.

**A block's origin is the longest tag path common to everything in it.**  :func:`common_path`
gives a block built from one module that module's full path, and a block where three modules
overlap their shared ancestor, which is the honest answer to "where did this come from".  Picking
one arbitrarily would be misleading.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import pypulseq as pp

from ..design import events as ev
from ..design.events import GRADIENT_KINDS, POINT_KINDS
from ..design.logic import BARRIER
from ..design.timing import EPS
from ..errors import CompileError, format_error
from .legalization import axis_gradient, check_limits
from .model import EXCLUSIVE_KINDS, PlacedEvent, PulseqReadyBlock, in_block_delay
from .verification import require_valid_contract, verify_ready_blocks

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from types import SimpleNamespace

    from pypulseq.opts import Opts

    from ..design.timing import Raster

__all__ = ['common_path', 'emit_blocks', 'required_duration']


def required_duration(events: Sequence[SimpleNamespace], opts: Opts) -> float:
    """
    Return the shortest block that can hold `events`, by pulseq's own rules.

    Mirrors what ``set_block`` computes: an RF needs its delay, shape and ringdown; an ADC needs
    its delay, window and *trailing* dead time; a gradient needs its delay and extent.
    """
    out = 0.0
    for e in events:
        kind = getattr(e, 'type', None)
        delay = float(getattr(e, 'delay', 0.0) or 0.0)
        if kind == 'rf':
            ring = float(getattr(e, 'ringdown_time', opts.rf_ringdown_time) or 0.0)
            out = max(out, delay + float(e.shape_dur) + ring)
        elif kind == 'adc':
            dead = float(getattr(e, 'dead_time', opts.adc_dead_time) or 0.0)
            out = max(out, delay + float(e.num_samples) * float(e.dwell) + dead)
        elif kind == BARRIER:
            continue
        else:
            out = max(out, float(pp.calc_duration(e)))
    return out


def _adc_conflict(edge: float, adcs: Sequence[PlacedEvent]) -> PlacedEvent | None:
    """
    Return an ADC whose sampling window the boundary `edge` falls strictly inside.

    A safety net rather than an expected condition: :func:`find_boundaries` only accepts times
    outside every reservation, and a reservation contains its ADC's window.  If this ever fires
    it is a compiler bug, and it fires with the information needed to find it rather than letting
    a silently-split readout reach the scanner.
    """
    for adc in adcs:
        if adc.start + EPS < edge < adc.end - EPS:
            return adc
    return None


def common_path(paths: Sequence[tuple[str, ...]]) -> tuple[str, ...]:
    """
    Return the longest tag path common to everything in a block.

    A block built from one module gets that module's full path; a block where three modules
    overlap gets their shared ancestor, which is the honest answer to "where did this come
    from" -- the alternative, picking one arbitrarily, would be misleading.
    """
    real = [p for p in paths if p]
    if not real:
        return ()
    common = list(real[0])
    for path in real[1:]:
        keep = 0
        for x, y in zip(common, path):
            if x != y:
                break
            keep += 1
        common = common[:keep]
        if not common:
            break
    return tuple(common)


def emit_blocks(  # noqa: C901, PLR0912
    seq: Any,
    edges: Sequence[float],
    placed: Sequence[PlacedEvent],
    targets: Mapping[int, float],
    opts: Opts,
    raster: Raster,
    notes: dict[str, list[str]],
) -> list[tuple[str, ...]]:
    """
    Add one pulseq block per interval of `edges` to `seq`, and return each block's origin path.

    Parameters
    ----------
    seq
        The :class:`pypulseq.Sequence` under construction.  Mutated: this is the stage whose whole
        purpose is the side effect.
    edges
        The boundaries from :func:`~seqcraft.compiler.boundaries.find_boundaries`.  Blocks
        are the consecutive pairs, so ``len(edges) - 1`` of them are added.
    placed, targets
        The placed events, and the label assignment times chosen for them.
    opts, raster
        The scanner, and the block-duration raster.
    notes
        ``category -> entries``, appended to: same-axis merges, resamplings and vector-norm
        findings made while emitting.  Aggregated into one warning per category at the end of
        the compile.

    Raises
    ------
    CompileError
        If a boundary would split an ADC's sampling window, if an interval is shorter than the
        events it holds, or if pypulseq refuses the block.
    HardwareLimitError
        If the *summed* waveform in a block exceeds ``max_grad`` or ``max_slew`` on an axis.

    Notes
    -----
    Both sweeps are single-pass.  Gradients are indexed per axis and advanced by a cursor, and
    the singles are pre-sorted by assignment time -- because 1700 TRs is roughly half a million
    events, and rescanning them per interval would be quadratic.
    """
    origins: list[tuple[str, ...]] = []
    previous_block_end: float | None = None

    grads: dict[str, list[PlacedEvent]] = {}
    for p in placed:
        if p.kind in GRADIENT_KINDS:
            grads.setdefault(p.event.channel, []).append(p)
    for pieces in grads.values():
        pieces.sort(key=lambda p: p.start)
    cursor = dict.fromkeys(grads, 0)
    active: dict[str, list[PlacedEvent]] = {ax: [] for ax in grads}

    # Labels are assigned at their *target ADC's* time rather than their own, so a boundary
    # landing between a label and the readout it addresses cannot change which ADC sees it.
    # Everything else is assigned where it sits.
    #
    # The secondary sort key is the event's own time, which keeps the emitted order intuitive --
    # but it is presentation only: pypulseq sorts a block's extensions by library id, so
    # intra-block label order carries no meaning.  `label_order_conflict` rejects the groups where
    # that would matter, so nothing here depends on it.
    schedule = sorted(
        (
            (targets.get(i, p.res_start), p.res_start, p)
            for i, p in enumerate(placed)
            if p.kind in EXCLUSIVE_KINDS or p.kind in POINT_KINDS
        ),
        key=lambda item: (item[0], item[1]),
    )
    assign_at = [t for t, _, _ in schedule]
    singles = [p for _, _, p in schedule]
    single_cursor = 0
    adcs = [p for p in placed if p.kind == 'adc']

    for index, (a, b) in enumerate(zip(edges, edges[1:])):
        block: list[SimpleNamespace] = []
        paths: list[tuple[str, ...]] = []

        while single_cursor < len(singles) and assign_at[single_cursor] < b - EPS:
            p = singles[single_cursor]
            single_cursor += 1
            if p.kind in POINT_KINDS and not hasattr(p.event, 'delay'):
                block.append(p.event)
            else:
                block.append(ev.derive(p.event, delay=in_block_delay(p, a, opts)))
            paths.append(p.path)

        for axis, pieces in grads.items():
            while cursor[axis] < len(pieces) and pieces[cursor[axis]].start < b - EPS:
                active[axis].append(pieces[cursor[axis]])
                cursor[axis] += 1
            active[axis] = [p for p in active[axis] if p.end > a + EPS]
            here = [p for p in active[axis] if p.start < b - EPS]
            if not here:
                continue
            crossing = [p for p in here if p.start < a - EPS or p.end > b + EPS]
            for edge in (a, b) if crossing else ():
                blocked = _adc_conflict(edge, adcs)
                if blocked is not None:
                    msg = format_error(
                        f'a block boundary at {edge * 1e6:.1f} us falls inside a gradient on axis '
                        f'{axis} that an ADC is sampling.',
                        {
                            'gradient from': ', '.join(sorted({p.where for p in crossing})),
                            'adc from': blocked.where,
                            'adc window': f'{blocked.start * 1e6:.1f} .. {blocked.end * 1e6:.1f} us',
                        },
                        [
                            'a readout gradient must stay one event -- ramp sampling and vendor '
                            'gridding both depend on that, so the compiler will not split it',
                            'the boundary comes from an explicit barrier, or this is a compiler '
                            'bug; please report it with the tree that produced it',
                        ],
                    )
                    raise CompileError(msg)
            grad = axis_gradient(axis, here, a, b, opts, notes, index)
            if grad is not None:
                block.append(grad)
                paths.extend(p.path for p in here)

        duration = raster.nearest(b - a)
        ready = PulseqReadyBlock(
            index=index,
            start=a,
            end=b,
            duration=duration,
            events=tuple(block),
            source_paths=tuple(paths),
            origin=common_path(paths),
        )
        require_valid_contract(
            'ready-block',
            verify_ready_blocks(
                (ready,),
                expected_first_index=index,
                expected_start=previous_block_end,
            ),
        )
        previous_block_end = ready.end
        origins.append(ready.origin)
        if not ready.events:
            seq.add_block(pp.make_delay(ready.duration))
            continue

        origin = ', '.join(sorted({'.'.join(p) for p in ready.source_paths if p})) or '?'
        check_limits(ready.events, opts, index, origin, notes, start=a)
        # pypulseq takes a block's duration as the max over its events, so an interval shorter
        # than its contents silently produces an off-raster block instead of an error.  Catch it
        # here, where the boundary that caused it can still be named.
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
        # add_block() takes no duration argument; a delay event is how pulseq states an explicit
        # block length, and `set_block` takes the max of it and the event extents.
        try:
            seq.add_block(*ready.events, pp.make_delay(ready.duration))
        except (ValueError, RuntimeError) as err:
            msg = format_error(
                f'pypulseq rejected block {index} at {a * 1e6:.1f} us: {err}',
                {
                    'duration_us': f'{ready.duration * 1e6:.1f}',
                    'events': ', '.join(getattr(e, 'type', '?') for e in ready.events),
                    'from': origin,
                },
                ['this is a compiler bug unless the tree contains raw events built against '
                 'a different Opts -- please report it with the tree that produced it'],
            )
            raise CompileError(msg) from err

    return origins
