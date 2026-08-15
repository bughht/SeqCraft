"""
The compiler: a logic-block tree in, legal pulseq blocks out.

This is a scheduling problem.  A logic block lets anything overlap anything; a pulseq block
holds **at most one RF, one ADC and one gradient per axis**, must last a whole number of block
rasters, must leave the RF dead time before a pulse and the ringdown after it, and must have
gradients that join continuously across its boundaries.  The compiler's job is to find block
boundaries satisfying all of that, and to combine whatever lands inside each one.

How boundaries are chosen
-------------------------
From two independent properties of an event, which the reservation model keeps apart.

**Indivisible** -- no boundary may fall strictly inside it.  An RF's dead time and ringdown, an
ADC's window and trailing dead time, a trigger's pulse: each is one hardware action rather than a
waveform, so splitting it is meaningless.

**Exclusive** -- a block may hold at most one, so a boundary is *required* between two.  RF and
ADC only.  Triggers are ``TRIGGERS`` extensions and pulseq accepts several per block, so treating
them as exclusive would invent a constraint the hardware does not have.

Putting a boundary somewhere in the gap between each pair of consecutive *exclusive* reservations
guarantees at most one RF and one ADC per block **by construction**, with nothing left to check.
Gradients constrain nothing -- they follow wherever the boundaries turn out to be, split and
summed as needed -- and labels are retimed onto the block holding the ADC they address, so a
boundary cannot change which readout sees them.

If two exclusive reservations overlap, the sequence is physically impossible and compilation
stops.  If nothing in the gap between two of them can be cut, that is impossible too, and the
error names what is in the way.

Three rules on overlap, each chosen deliberately
------------------------------------------------
**Gradients on different axes: silent.**  A phase-encode blip beside a slice rewinder beside a
readout prephaser is the normal way to build a sequence, not a problem.  Warning about it would
only teach people to ignore warnings.

**Gradients on the same axis: warned, then summed.**  Summing is almost always what was meant --
a rewinder sharing time with a prephaser -- so it happens automatically.  But it is the one
place the compiler changes a waveform, so it says so.

**Two RF or two ADC events overlapping: an error.**  You cannot transmit twice at once, sample
twice at once, or transmit and receive at once.

Limits are checked *after* summing
----------------------------------
Two individually legal gradients on one axis can sum to an illegal one, and no module can see
that in isolation: adding an area-100 and an area-200 trapezoid on a 40 mT/m, 150 T/m/s system
reaches 93 % of the amplitude limit but **189 %** of the slew limit.  So amplitude and slew are
measured on the compiled waveform, which is the only place the truth is visible.

Examples
--------
>>> import pypulseq as pp
>>> import seqcraft as sc
>>> system = sc.System.preset('generic_3t')
>>> o = system.default
>>> gentle = {'area': 100.0, 'duration': 2e-3, 'rise_time': 200e-6, 'system': o}
>>> lb = sc.LogicBlock('demo')
>>> _ = lb.add(0.0, pp.make_trapezoid('x', **gentle))
>>> _ = lb.add(0.0, pp.make_trapezoid('y', **gentle))   # different axis: nothing to report
>>> out = sc.compile(lb, system)
>>> out.report.ok, len(out.report.warnings)
(True, 0)
>>> out.n_blocks
1
"""

from __future__ import annotations

import bisect
import functools
import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np
import pypulseq as pp

from . import events as ev
from . import units
from ._compiler.model import (
    ADDRESS_KEYS,
    AXES,
    EXCLUSIVE_KINDS,
    GRADIENT_KINDS,
    HANDLED_KINDS,
    INDIVISIBLE_KINDS,
    LABEL_KINDS,
    POINT_KINDS,
    PlacedEvent,
    PulseqReadyBlock,
)
from ._compiler.verification import (
    require_valid_contract,
    verify_placed_events,
    verify_ready_blocks,
)
from .errors import CompileError, format_error
from .logic import BARRIER, flatten
from .report import Issue, Report
from .timing import EPS, TICKS_PER_SECOND, Raster, exact_diff, to_ticks
from .validate import merge_definitions

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Sequence
    from pathlib import Path
    from types import SimpleNamespace

    from pypulseq.opts import Opts

    from .geometry import Geometry
    from .logic import LogicBlock
    from .system import System

__all__ = ['CompiledSequence', 'WriteResult', 'compile_sequence']

# Temporary aliases keep the established compiler implementation readable while each owning stage
# moves behind the private package.  `_Placed` names the same class; there is no parallel IR.
_Placed = PlacedEvent
_GRAD = GRADIENT_KINDS
_POINT = POINT_KINDS
_INDIVISIBLE = INDIVISIBLE_KINDS
_EXCLUSIVE = EXCLUSIVE_KINDS
_LABEL = LABEL_KINDS
_HANDLED = HANDLED_KINDS
_AXES = AXES
_ADDRESS_KEYS = ADDRESS_KEYS

#: Types pypulseq 1.5 can produce that seqcraft recognises but does not emit, each with what it
#: is and the way round.  Kept apart from "never heard of it" because the remedy differs
#: completely -- these have one, an unknown type is a bug or a version skew.
_UNSUPPORTED: dict[str, tuple[str, tuple[str, ...]]] = {
    'rot3D': (
        'the pulseq rotation extension',
        (
            'seqcraft bakes rotations into the waveform instead of emitting the extension, so '
            'that limits, moments and k-space all describe what actually plays',
            'rotate at build time, the way SpiralVDS rotates its interleaves and the diffusion '
            'modules resolve a direction',
            'for a one-off, pypulseq.rotate() / rotate_3d() return rotated gradient events that '
            'can be added to a LogicBlock directly',
        ),
    ),
    'soft_delay': (
        'the pulseq soft-delay extension',
        (
            'a soft delay must occupy a block containing no other events, which the compiler '
            'does not yet arrange',
            'use a plain delay event (pp.make_delay) if the duration is fixed at compile time',
        ),
    ),
    'rf_shim': (
        'the pulseq RF-shim extension',
        ('not yet emitted; leave it out, or open an issue with the sequence that needs it',),
    ),
}
# --------------------------------------------------------------------------------- placing
def _intrinsic_duration(event: SimpleNamespace) -> float:
    """Return how long an event is active, excluding its own leading delay."""
    kind = getattr(event, 'type', None)
    if kind == 'trap':
        return float(event.rise_time) + float(event.flat_time) + float(event.fall_time)
    if kind == 'grad':
        shape = getattr(event, 'shape_dur', None)
        return float(shape) if shape is not None else float(event.tt[-1])
    if kind == 'rf':
        return float(event.shape_dur)
    if kind == 'adc':
        return float(event.num_samples) * float(event.dwell)
    if kind in ('trigger', 'output'):
        # pypulseq builds scanner inputs as 'trigger' and scanner outputs as 'output', and both
        # hold their block open for the length of the pulse.
        return float(getattr(event, 'duration', 0.0) or 0.0)
    return 0.0


def _unsupported(kind: str | None, path: tuple[str, ...], node_t: float) -> CompileError:
    """Build the error for an event type the compiler will not emit."""
    where = '.'.join(path) if path else '(untagged)'
    if kind in _UNSUPPORTED:
        what, hints = _UNSUPPORTED[kind]
        return CompileError(format_error(
            f'{what} ({kind!r}) is not supported by the compiler.',
            {'from': where, 'at': f'{node_t * 1e6:.1f} us'},
            list(hints),
        ))
    return CompileError(format_error(
        f'unknown event type {kind!r}.',
        {'from': where, 'at': f'{node_t * 1e6:.1f} us', 'handled': ', '.join(sorted(_HANDLED))},
        [
            'LogicBlock.add() accepts any object with a .type attribute, so a typo or a '
            'hand-built namespace reaches the compiler unchecked',
            'if this is a real pulseq event from a newer pypulseq, seqcraft needs updating -- '
            'please open an issue naming the type',
        ],
    ))


def _place(root: LogicBlock, opts: Opts) -> list[_Placed]:
    """
    Flatten `root` and resolve every event to absolute times and reservations.

    A ``delay`` event is the one special case: pypulseq stores its length in ``delay`` and it
    has no waveform, so it occupies ``[t, t + delay]`` rather than starting after its own
    delay.  It exists to hold a block open -- a b=0 diffusion volume that must fill the same
    slot as an encoded one.

    Raises
    ------
    CompileError
        If any leaf carries a type the compiler does not emit.  Rejecting here is the point:
        every branch below is a positive match, so an unrecognised type would otherwise be
        placed with a zero-width reservation, collected by nothing, and silently lost.
    """
    out: list[_Placed] = []
    for node_t, event, path in flatten(root):
        kind = getattr(event, 'type', None)
        if kind not in _HANDLED:
            raise _unsupported(kind, path, node_t)
        if kind == BARRIER:
            out.append(_Placed(node_t, node_t, node_t, node_t, node_t, event, path))
            continue
        delay = float(getattr(event, 'delay', 0.0) or 0.0)
        if kind == 'delay':
            out.append(_Placed(node_t, node_t, node_t + delay, node_t, node_t + delay, event, path))
            continue

        start = node_t + delay
        end = start + _intrinsic_duration(event)
        if kind == 'rf':
            ring = float(getattr(event, 'ringdown_time', opts.rf_ringdown_time) or 0.0)
            res_start, res_end = node_t, end + ring
        elif kind == 'adc':
            dead = float(getattr(event, 'dead_time', opts.adc_dead_time) or 0.0)
            res_start, res_end = node_t, end + dead
        else:
            res_start, res_end = start, end
        out.append(_Placed(node_t, start, end, res_start, res_end, event, path))
    return out


def _check_exclusive(placed: Sequence[_Placed]) -> None:
    """
    Raise if two RF or ADC reservations overlap.

    Checked before any scheduling, because it is the one class of conflict no choice of block
    boundaries can fix: the sequence is asking the hardware to do two things it physically
    cannot do at once.  The message names both tag paths and the overlap -- information
    pypulseq cannot give, since it would fail much later, on a block index, with
    ``Multiple ADC events were specified``.
    """
    exclusive = sorted(
        (p for p in placed if p.kind in ('rf', 'adc')),
        key=lambda p: (p.res_start, p.res_end),
    )
    for a, b in zip(exclusive, exclusive[1:]):
        if b.res_start < a.res_end - EPS:
            overlap = min(a.res_end, b.res_end) - b.res_start
            pair = f'{a.kind.upper()} and {b.kind.upper()}'
            reason = {
                'RF and RF': 'the transmitter cannot play two pulses at once',
                'ADC and ADC': 'the receiver cannot open two sampling windows at once',
            }.get(pair, 'the coil cannot transmit and receive at the same time')
            msg = format_error(
                f'{pair} overlap by {overlap * 1e6:.1f} us -- {reason}.',
                {
                    f'{a.kind} at {a.start * 1e6:.1f} us': a.where,
                    'reserves': f'{a.res_start * 1e6:.1f} .. {a.res_end * 1e6:.1f} us',
                    f'{b.kind} at {b.start * 1e6:.1f} us': b.where,
                    'reserves ': f'{b.res_start * 1e6:.1f} .. {b.res_end * 1e6:.1f} us',
                },
                [
                    f'move one of them by at least {overlap * 1e6:.1f} us',
                    'a reservation includes the RF dead time and ringdown, and the ADC dead '
                    'time before and after sampling -- so two events can conflict while their '
                    'waveforms do not visibly touch',
                ],
            )
            raise CompileError(msg)


# ---------------------------------------------------------------------------- boundaries
class _Spans:
    """
    A sorted set of intervals, answering "is this time strictly inside one of them?" in log time.

    Built once per compile.  With 1700 TRs there are half a million spans and a boundary
    candidate for each, so a linear test per candidate would be quadratic.
    """

    def __init__(self, spans: Iterable[tuple[float, float]]) -> None:
        pairs = [(lo, hi) for lo, hi in spans if hi > lo + EPS]
        self._starts = sorted(lo for lo, _ in pairs)
        self._ends = sorted(hi for _, hi in pairs)

    def contains_strictly(self, t: float) -> bool:
        """``True`` if some span ``(lo, hi)`` satisfies ``lo < t < hi``."""
        started = bisect.bisect_left(self._starts, t - EPS)
        finished = bisect.bisect_right(self._ends, t + EPS)
        return started > finished


def _label_targets(placed: Sequence[_Placed]) -> dict[int, float]:
    """
    Return, per label, the time it should be *assigned* at -- its target ADC's reservation start.

    A pulseq label is a running register, not an instant.  The interpreter applies a block's
    labels when it reaches that block, and an ADC in the same block then samples the state:
    pypulseq's ``evaluate_labels`` applies the labels first and records afterwards.  So a label
    may live in **any** block strictly after the previous ADC's and at or before its target
    ADC's, with identical results.

    That freedom means the compiler has to *choose*, and choosing by containment is wrong.  A
    boundary pushed later than the ADC's reservation end -- which is what the gap-midpoint
    fallback does whenever a gradient covers the natural candidates -- puts a label placed
    comfortably after one readout into that readout's own block, silently overwriting its
    k-space address.  Measured: two ADCs, a label between them, both ending up with the same
    LIN.

    Assigning instead to the block holding the **first ADC at or after the label's own time**
    makes the outcome independent of where boundaries land, which is the property worth having.
    A label with no ADC after it keeps containment assignment -- there is nothing for it to
    address, so there is nothing to get wrong.

    A **barrier is never crossed**.  A barrier is an explicit request for a seam at that instant,
    so moving a label over one would override the one instruction the user gave about block
    structure.  It is also what makes a barrier a working remedy for the order-dependent case in
    :func:`_label_order_conflict`: separated into different blocks, two labels are applied in
    block order, which pulseq does define.

    Returns
    -------
    dict
        Index into `placed` -> assignment time.  Labels absent from the mapping keep their own
        ``res_start``: either there is no ADC after them, or a barrier stands in the way.
    """
    adc_starts = sorted(p.res_start for p in placed if p.kind == 'adc')
    barriers = sorted(p.node_t for p in placed if p.kind == BARRIER)
    out: dict[int, float] = {}
    if not adc_starts:
        return out
    for i, p in enumerate(placed):
        if p.kind not in _LABEL:
            continue
        j = bisect.bisect_left(adc_starts, p.res_start - EPS)
        if j >= len(adc_starts):
            continue
        target = adc_starts[j]
        k = bisect.bisect_right(barriers, p.res_start + EPS)
        if k < len(barriers) and barriers[k] < target - EPS:
            continue
        out[i] = target
    return out


def _label_order_conflict(
    placed: Sequence[_Placed],
    targets: dict[int, float],
) -> CompileError | None:
    """
    Return an error if two labels would share a block but their order matters.

    **pypulseq discards the order labels were added in.**  ``Sequence/block.py`` sorts a block's
    extensions by their library reference id -- "we rely on the sorting of the extension IDs" --
    so ``add_block(set, inc)`` and ``add_block(inc, set)`` produce the same block and the same
    result.  Verified: both give ``LIN = 10``, never 13.

    Only order-*independent* groups are therefore expressible within one block:

    * labels on **different** keys -- they do not interact;
    * several ``labelinc`` on one key -- addition commutes.

    A ``labelset`` sharing a key with anything else does not commute, so pulseq cannot express
    the intent and guessing would silently pick one of two different k-space addressings.  This
    is not new with target assignment -- containment could already put two such labels in one
    block -- but grouping by target makes it reachable on purpose, so it has to be detected.
    """
    groups: dict[tuple[float, str], list[_Placed]] = {}
    for i, p in enumerate(placed):
        if p.kind not in _LABEL:
            continue
        key = (targets.get(i, p.res_start), str(getattr(p.event, 'label', '?')))
        groups.setdefault(key, []).append(p)

    for (_, label_key), members in sorted(groups.items()):
        if len(members) < 2 or all(m.kind == 'labelinc' for m in members):
            continue
        kinds = ', '.join(f'{m.kind}={getattr(m.event, "value", "?")} from {m.where}'
                          for m in members)
        return CompileError(format_error(
            f'{len(members)} labels set or increment {label_key!r} for the same readout, and '
            f'their order cannot be expressed.',
            {
                'labels': kinds,
                'times': ', '.join(f'{m.res_start * 1e6:.1f} us' for m in members),
            },
            [
                'pypulseq sorts a block\'s extensions by library id, so the order they are '
                'added in is discarded -- only labels on different keys, or several labelinc '
                'on one key, are order-independent',
                'put a barrier between them so they land in separate blocks, where block order '
                'does determine the order',
                'or combine them into the single label the readout should actually see',
            ],
        ))
    return None




def _expected_addresses(
    placed: Sequence[_Placed],
    targets: dict[int, float],
) -> list[dict[str, int]]:
    """
    Fold the tree's labels the way the interpreter will, giving the address of each ADC.

    Mirrors ``Sequence.evaluate_labels(evolution='adc')``: labels are applied in assignment
    order, and the running state is recorded at every ADC.  Comparing the two is what makes
    label placement checkable rather than merely intended -- and it catches the failure the
    duplicate-address check cannot see, an addressing shifted by one but still unique.
    """
    adcs = sorted((p for p in placed if p.kind == 'adc'), key=lambda p: p.res_start)
    labels = sorted(
        (
            (targets.get(i, p.res_start), p.res_start, p)
            for i, p in enumerate(placed)
            if p.kind in _LABEL
        ),
        key=lambda item: (item[0], item[1]),
    )
    state: dict[str, int] = {}
    out: list[dict[str, int]] = []
    cursor = 0
    for adc in adcs:
        while cursor < len(labels) and labels[cursor][0] <= adc.res_start + EPS:
            _, _, p = labels[cursor]
            cursor += 1
            key = str(getattr(p.event, 'label', '?'))
            value = int(getattr(p.event, 'value', 0))
            state[key] = state.get(key, 0) + value if p.kind == 'labelinc' else value
        out.append(dict(state))
    return out


def _orphan_label_issues(
    placed: Sequence[_Placed],
    targets: dict[int, float],
) -> list[Issue]:
    """
    Report labels with no ADC after them.

    Such a label addresses nothing, so it is placed by containment -- and containment can put it
    in the *last* ADC's block, changing that readout's address.  Rather than guess, say so: the
    fix is always to move the label before the readout it was meant for.
    """
    out: list[Issue] = []
    for i, p in enumerate(placed):
        if p.kind in _LABEL and i not in targets:
            out.append(Issue(
                'label',
                p.where,
                f'label {getattr(p.event, "label", "?")} at {p.res_start * 1e6:.1f} us has no '
                f'ADC after it, so it addresses no readout; it may land in the preceding '
                f'readout\'s block and change that address',
                'warning',
            ))
    return out


def _covering(placed: Sequence[_Placed], t: float) -> _Placed | None:
    """Return an indivisible event whose reservation strictly contains `t`, for error messages."""
    for p in placed:
        if p.kind in _INDIVISIBLE and p.res_start + EPS < t < p.res_end - EPS:
            return p
    return None


def _barrier_conflict(
    barrier: _Placed,
    at: float,
    placed: Sequence[_Placed],
) -> CompileError:
    """Build the error for a barrier that would split an indivisible event."""
    blocker = _covering(placed, at)
    detail = {'barrier': barrier.where, 'at': f'{at * 1e6:.1f} us'}
    if blocker is not None:
        detail[f'would split {blocker.kind}'] = blocker.where
        detail['which reserves'] = (
            f'{blocker.res_start * 1e6:.1f} .. {blocker.res_end * 1e6:.1f} us'
        )
    return CompileError(format_error(
        'a barrier falls inside an event that cannot be split.',
        detail,
        [
            'an RF reserves its dead time and ringdown, an ADC its window and trailing dead '
            'time, and a trigger its whole pulse -- a boundary inside any of them has no '
            'meaning, since each is one hardware action rather than a waveform',
            'move the barrier outside that span, or drop it and let the compiler choose',
        ],
    ))


def _gap_boundary(
    a_end: float,
    b_start: float,
    raster: Raster,
    cuttable: Callable[[float], bool],
) -> float | None:
    """
    Return a legal boundary strictly between two exclusive reservations, or ``None``.

    The midpoint is tried first because for an EPI train it is the middle of the blip, which is
    the natural place to split the readout gradient.  The two edges are the fallbacks for a gap
    that something indivisible partly covers.
    """
    options = (
        raster.ceil(0.5 * (a_end + b_start)),
        raster.ceil(a_end),
        raster.floor(b_start),
    )
    for t in options:
        if a_end - EPS <= t <= b_start + EPS and cuttable(t):
            return t
    return None


def _gap_blocked(a_end: float, b_start: float, placed: Sequence[_Placed]) -> CompileError:
    """Build the error for a mandatory gap that nothing may be cut in."""
    mid = 0.5 * (a_end + b_start)
    blocker = _covering(placed, mid)
    detail = {'gap': f'{a_end * 1e6:.1f} .. {b_start * 1e6:.1f} us'}
    if blocker is not None:
        detail[f'covered by {blocker.kind}'] = blocker.where
        detail['which reserves'] = (
            f'{blocker.res_start * 1e6:.1f} .. {blocker.res_end * 1e6:.1f} us'
        )
    return CompileError(format_error(
        'two RF/ADC events need a block boundary between them, but nothing in the gap can be '
        'cut.',
        detail,
        [
            'pulseq holds one RF and one ADC per block, so consecutive ones must be separated -- '
            'and the only times available are outside every indivisible span',
            'usually a trigger stretched across the gap: shorten it, or move it so it sits '
            'wholly within one of the two blocks',
        ],
    ))


def _boundaries(
    placed: Sequence[_Placed],
    total: float,
    raster: Raster,
    max_block: float,
) -> list[float]:
    """
    Return the sorted block boundaries.

    Two independent properties decide where a boundary may and must go, and keeping them apart
    is what makes triggers work.  An event is

    **indivisible** when no boundary may fall strictly inside its span -- an RF's dead time and
    ringdown, an ADC's window and trailing dead time, a trigger's pulse.  Splitting any of them
    is meaningless: they are single hardware actions, not waveforms.

    **exclusive** when a block may hold at most one -- RF and ADC, because pulseq stores one of
    each per block.  Triggers are *not* exclusive: they are TRIGGERS extensions and pulseq
    accepts several in one block, so demanding a boundary between two of them would invent a
    constraint the hardware does not have.

    From those, three rules produce boundaries:

    **Mandatory.**  One boundary somewhere in the gap between each pair of consecutive
    *exclusive* reservations, which guarantees at most one RF and one ADC per block by
    construction.  Zero and the total duration are boundaries too, as are explicit barriers.

    **Opportunistic.**  Every gradient edge and reservation edge is a *candidate*, accepted only
    if it falls strictly inside neither a gradient nor an indivisible span.  That is what keeps a
    trapezoid a trapezoid: a slice rephaser on z overlapping a phase blip on y stays two clean
    events in one block, rather than being cut where the blip happens to end.  It is also why a
    readout gradient survives whole -- its own edges are candidates, and nothing else's edge
    lands inside it.

    **Forced.**  Intervals longer than `max_block` are subdivided, since pulseq stores a block
    duration in a fixed-width field.

    Raises
    ------
    CompileError
        If a barrier sits inside an indivisible span, or if a mandatory gap contains no legal
        boundary at all because something indivisible covers the whole of it.
    """
    exclusive = sorted(
        ((p.res_start, p.res_end) for p in placed if p.kind in _EXCLUSIVE),
        key=lambda s: s[0],
    )
    grad_spans = _Spans((p.start, p.end) for p in placed if p.kind in _GRAD)
    indivisible = _Spans(
        (p.res_start, p.res_end) for p in placed if p.kind in _INDIVISIBLE
    )

    def cuttable(t: float) -> bool:
        """Hard requirement on *every* boundary: it may not split a single hardware action."""
        return not indivisible.contains_strictly(t)

    def acceptable(t: float) -> bool:
        """Additionally leaves every gradient whole -- a preference, not a requirement."""
        return cuttable(t) and not grad_spans.contains_strictly(t)

    marks = {0.0, total}
    for p in placed:
        if p.kind != BARRIER:
            continue
        at = raster.nearest(p.node_t)
        if not cuttable(at):
            raise _barrier_conflict(p, at, placed)
        marks.add(at)

    # Snap the start of anything **down** and the end of anything **up**, so a block only ever
    # grows to fit its contents.  Nearest-rounding an end is a real bug: an ADC whose window plus
    # trailing dead time lands 4 us past a raster edge would get a block 6 us too short, which
    # pypulseq reports as an unaligned duration rather than as the missing microseconds it is.
    candidates: set[float] = set()
    for p in placed:
        if p.kind in _GRAD:
            candidates.add(raster.floor(p.start))
            candidates.add(raster.ceil(p.end))
        elif p.kind in _INDIVISIBLE:
            candidates.add(raster.floor(p.res_start))
            candidates.add(raster.ceil(p.res_end))
    marks.update(t for t in candidates if acceptable(t))

    # Guarantee one boundary per reservation gap.  Anything in the gap will do -- the gap is by
    # definition outside every reservation -- so prefer a point already chosen, then a gradient
    # edge inside the gap, and fall back to the midpoint, which for an EPI train is the middle
    # of the blip and therefore the natural place to split the readout gradient.
    #
    # Gaps are disjoint and ordered, so `cursor` only ever moves forward and one scan answers
    # every gap.  Re-sorting `marks` after each added midpoint instead is quadratic, and it fires
    # on precisely the case this fallback exists for: in an EPI train a continuous readout
    # gradient makes both natural candidates unacceptable, so *every* echo lands here.  A midpoint
    # can only fall inside its own gap, so no later gap can need to see it.
    #
    # The chosen point must be `cuttable`, not merely inside the gap.  Being outside every
    # *exclusive* reservation is guaranteed by construction, but something indivisible and
    # non-exclusive -- a trigger -- can still cover part or all of a gap, and cutting it would be
    # meaningless.  Cutting a *gradient* here is fine and expected: for an EPI train the midpoint
    # is the middle of the blip, which is the natural place to split the readout.
    chosen = sorted(marks)
    cursor = 0
    extra: list[float] = []
    for (_, a_end), (b_start, _) in zip(exclusive, exclusive[1:]):
        while cursor < len(chosen) and chosen[cursor] < a_end - EPS:
            cursor += 1
        if cursor < len(chosen) and chosen[cursor] <= b_start + EPS:
            continue
        at = _gap_boundary(a_end, b_start, raster, cuttable)
        if at is None:
            raise _gap_blocked(a_end, b_start, placed)
        extra.append(at)
    marks.update(extra)

    snapped = sorted({raster.nearest(min(max(0.0, t), total)) for t in marks})

    out: list[float] = []
    for a, b in zip(snapped, snapped[1:]):
        out.append(a)
        n_extra = int(math.ceil((b - a) / max_block)) - 1
        if n_extra > 0:
            step = raster.nearest((b - a) / (n_extra + 1))
            out.extend(a + step * (k + 1) for k in range(n_extra) if a + step * (k + 1) < b - EPS)
    out.append(snapped[-1])
    return sorted(set(out))


# ------------------------------------------------------------------------------ gradients
def _superpose(
    pieces: Sequence[_Placed],
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


def _axis_gradient(
    axis: str,
    pieces: Sequence[_Placed],
    a: float,
    b: float,
    opts: Opts,
    issues: list[Issue],
    block_index: int,
) -> SimpleNamespace | None:
    """
    Return the single gradient event for `axis` over the interval ``[a, b)``, or ``None``.

    Fast path -- one gradient, entirely inside the interval -- passes through untouched apart
    from its in-block delay.  No new shape, so a spiral readout or a long diffusion lobe costs
    exactly what it did before compilation.

    Otherwise the pieces are superposed exactly (:func:`_superpose`) and the result is emitted in
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
            return ev.derive(p.event, delay=_in_block_delay(p, a, opts))

    if len(pieces) > 1:
        issues.append(
            Issue(
                'grad_merge',
                f'block {block_index}',
                f'{len(pieces)} gradients overlap on axis {axis} over '
                f'{a * 1e6:.1f}..{b * 1e6:.1f} us and were summed: '
                + ', '.join(sorted({p.where for p in pieces})),
                'warning',
            )
        )

    rel_ps, amps = _superpose(pieces, a, b)
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
        # Limits are measured on the compiled block by _limit_issues, which reports rather than
        # raises; letting make_arbitrary_grad raise here would turn a reportable violation into a
        # failed compile, and would do it before the waveform is complete.
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

    return _resampled(axis, rel_ps, amps, opts, issues, block_index, pieces)


def _resampled(
    axis: str,
    rel_ps: np.ndarray,
    amps: np.ndarray,
    opts: Opts,
    issues: list[Issue],
    block_index: int,
    pieces: Sequence[_Placed],
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
    issues.append(
        Issue(
            'grad_resample',
            f'block {block_index}',
            f'axis {axis}: a trapezoid summed with a raster-centre waveform bends both on and '
            f'off the gradient raster, which no pulseq gradient event can hold; resampled onto '
            f'the raster, moving the waveform by at most '
            f'{error / float(opts.max_grad) * 100:.2f} % of max_grad. From '
            + ', '.join(sorted({p.where for p in pieces})),
            'warning',
        )
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


def _in_block_delay(p: _Placed, block_start: float, opts: Opts) -> float:
    """
    Return an event's delay within its block, quantised onto that event's own raster.

    Two reasons this cannot be a plain subtraction.  The absolute times come from arithmetic over
    a sequence that may run for minutes, so by the last TR the float resolution is coarser than a
    picosecond and ``p.start - block_start`` drifts -- pypulseq then reports an RF delay of
    ``129.9999999986us`` and rejects the block.  And pulseq requires each event's delay to sit on
    its own raster -- 1 us for RF, 100 ns for ADC and 10 us for gradients on Siemens, whatever
    the scanner reports elsewhere.  Subtracting in integer ticks and then snapping satisfies both.
    """
    dt = max(0.0, exact_diff(p.start, block_start))
    raster = {
        'rf': Raster(float(opts.rf_raster_time), 'RF'),
        'adc': Raster(float(opts.adc_raster_time), 'ADC'),
    }.get(p.kind, Raster(float(opts.grad_raster_time), 'gradient'))
    return raster.nearest(dt)


def _required_duration(events: Sequence[SimpleNamespace], opts: Opts) -> float:
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


def _adc_conflict(edge: float, adcs: Sequence[_Placed]) -> _Placed | None:
    """
    Return an ADC whose sampling window the boundary `edge` falls strictly inside.

    A safety net rather than an expected condition: :func:`_boundaries` only accepts times
    outside every reservation, and a reservation contains its ADC's window.  If this ever fires
    it is a compiler bug, and it fires with the information needed to find it rather than letting
    a silently-split readout reach the scanner.
    """
    for adc in adcs:
        if adc.start + EPS < edge < adc.end - EPS:
            return adc
    return None


# --------------------------------------------------------------------------------- limits
def _limit_issues(
    events: Iterable[SimpleNamespace],
    opts: Opts,
    block_index: int,
    origin: str,
) -> list[Issue]:
    """
    Measure amplitude and slew on a compiled block and describe any violation.

    Reported rather than raised so a run surfaces every violation at once; one mistimed module
    usually produces a run of them, and stopping at the first hides the pattern.

    The vector norm across simultaneous axes is a warning, not an error: two axes ramping
    together reach ``sqrt(2)`` times the per-axis slew in vector magnitude, which real
    amplifiers permit.  It becomes a hard limit only once a rotation can concentrate the whole
    vector onto one physical axis.
    """
    out: list[Issue] = []
    raster = float(opts.grad_raster_time)
    gamma = float(opts.gamma)
    for kind, where, achieved, limit in ev.check_limits(events, opts, raster):
        pulseq_unit, unit = ('Hz/m', 'mT/m') if kind.startswith('grad') else ('Hz/m/s', 'T/m/s')
        to_display = functools.partial(
            units.convert, from_unit=pulseq_unit, to_unit=unit, gamma=gamma
        )
        out.append(
            Issue(
                f'{kind}_limit',
                f'block {block_index}',
                f'{kind.replace("_", " ")} on {where} reaches {to_display(achieved):.1f} {unit}, '
                f'limit {to_display(limit):.1f} {unit} '
                f'({achieved / limit * 100:.0f}%); from {origin}',
                'warning' if kind.endswith('_norm') else 'error',
            )
        )
    return out


# ------------------------------------------------------------------------------- the pass
def compile_sequence(  # noqa: C901, PLR0912, PLR0915
    root: LogicBlock,
    system: System,
    *,
    geometry: Geometry | None = None,
    name: str = '',
    regime: str = 'default',
    definitions: dict[str, Any] | None = None,
) -> CompiledSequence:
    """
    Turn a logic-block tree into a :class:`pypulseq.Sequence`.

    Parameters
    ----------
    root
        The tree to compile.  Usually built by adding module outputs to one
        :class:`~seqcraft.core.logic.LogicBlock`.
    system
        The scanner: rasters, dead times and limits.
    geometry
        Optional; contributes ``[DEFINITIONS]`` to the written file.
    name
        Sequence name for the definitions.  Defaults to ``root.tag``.
    regime
        Which of `system`'s named limit regimes to validate against.  Modules may design
        against a derated regime individually; this is the ceiling the *combined* waveform
        has to respect.
    definitions
        Extra ``[DEFINITIONS]`` entries, merged with a collision check so two sources claiming
        the same key with different values raises rather than one silently winning.

    Returns
    -------
    CompiledSequence
        The pypulseq sequence, the compile report, and per-block provenance.

    Raises
    ------
    CompileError
        Two RF or ADC events overlap in time, an absolute start is negative, or a block
        boundary would have to fall inside a gradient an ADC is sampling.

    Notes
    -----
    Amplitude and slew violations are *reported*, not raised.  Call
    :meth:`CompiledSequence.check` then
    :meth:`~seqcraft.core.report.Report.raise_if_failed` to stop on them.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> lb = sc.LogicBlock('t').add(0.0, pp.make_delay(1e-3))
    >>> out = compile_sequence(lb, system)
    >>> round(out.duration_s * 1e3, 1)
    1.0
    """
    opts = system.limits(regime)
    raster = system.block_raster
    placed = _place(root, opts)
    require_valid_contract('placed-event', verify_placed_events(placed))

    negative = [p for p in placed if p.res_start < -EPS]
    if negative:
        worst = min(negative, key=lambda p: p.res_start)
        msg = format_error(
            f'event starts at {worst.res_start * 1e6:.1f} us, before the start of the sequence.',
            {'event': worst.kind, 'from': worst.where},
            [
                'a negative start usually means a module was placed at '
                '"te - module.time_to_echo" with a TE shorter than the module needs',
                'increase TE, or shift the whole tree later',
            ],
        )
        raise CompileError(msg)

    # A gradient off the gradient raster used to be snapped silently by _in_block_delay, moving it
    # by up to half a raster -- a gradient asked for at 5 us played at 10 us, and m0 could not see
    # it because area does not depend on when a lobe plays.  There is no correct snap to make: the
    # tree asked for something the hardware cannot do, and only the caller knows which way to round.
    grad_raster = system.grad_raster
    off = [p for p in placed if p.kind in _GRAD and not grad_raster.holds(p.start)]
    if off:
        worst = off[0]
        msg = format_error(
            f'gradient starts at {worst.start * 1e6:.3f} us, which is not a multiple of the '
            f'{grad_raster.dt * 1e6:.0f} us gradient raster.',
            {
                'from': worst.where,
                'nearest below': f'{grad_raster.floor(worst.start) * 1e6:.1f} us',
                'nearest above': f'{grad_raster.ceil(worst.start) * 1e6:.1f} us',
                'others': len(off) - 1,
            },
            [
                'round the start time where it is computed, with system.grad_raster.ceil(t)',
                'a start derived from "te - module.time_to_echo" lands off the raster whenever TE '
                'does, so quantise TE first',
            ],
        )
        raise CompileError(msg)

    _check_exclusive(placed)

    # Ceiling, not nearest: rounding the total *down* truncates the final block below what its
    # events need, which surfaces later as an unaligned block duration rather than as the four
    # microseconds it actually is.
    total = raster.ceil(max((p.res_end for p in placed), default=0.0))
    if total <= 0.0:
        msg = format_error(
            'nothing to compile: the tree contains no events that occupy time.',
            {'tag': root.tag or '(untagged)', 'nodes': len(root)},
        )
        raise CompileError(msg)

    issues: list[Issue] = [
        Issue('raster', p.where, f'start {p.res_start * 1e6:.4f} us snapped to '
              f'{raster.nearest(p.res_start) * 1e6:.1f} us', 'warning')
        for p in placed
        if p.kind in ('rf', 'adc') and not raster.holds(p.res_start)
    ]

    max_block = raster.dt * 2**24
    edges = _boundaries(placed, total, raster, max_block)

    seq = pp.Sequence(system=opts)
    origins: list[tuple[str, ...]] = []
    previous_block_end: float | None = None

    # Sort once per axis so each interval is served by a sweep rather than a rescan: 1700 TRs
    # is ~500k events, and rescanning per interval would be quadratic.
    grads: dict[str, list[_Placed]] = {}
    for p in placed:
        if p.kind in _GRAD:
            grads.setdefault(p.event.channel, []).append(p)
    for pieces in grads.values():
        pieces.sort(key=lambda p: p.start)
    cursor = dict.fromkeys(grads, 0)
    active: dict[str, list[_Placed]] = {ax: [] for ax in grads}

    # Labels are assigned at their *target ADC's* time rather than their own, so a boundary
    # landing between a label and the readout it addresses cannot change which ADC sees it.
    # Everything else is assigned where it sits.
    #
    # The secondary sort key is the event's own time, which keeps the emitted order intuitive --
    # but it is presentation only: pypulseq sorts a block's extensions by library id, so
    # intra-block label order carries no meaning.  _label_order_conflict rejects the groups where
    # that would matter, so nothing here depends on it.
    targets = _label_targets(placed)
    conflict = _label_order_conflict(placed, targets)
    if conflict is not None:
        raise conflict
    issues.extend(_orphan_label_issues(placed, targets))
    schedule = sorted(
        (
            (targets.get(i, p.res_start), p.res_start, p)
            for i, p in enumerate(placed)
            if p.kind in _EXCLUSIVE or p.kind in _POINT
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
            if p.kind in _POINT and not hasattr(p.event, 'delay'):
                block.append(p.event)
            else:
                block.append(ev.derive(p.event, delay=_in_block_delay(p, a, opts)))
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
            grad = _axis_gradient(axis, here, a, b, opts, issues, index)
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
            origin=_common_path(paths),
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
        issues.extend(_limit_issues(ready.events, opts, index, origin))
        # pypulseq takes a block's duration as the max over its events, so an interval shorter
        # than its contents silently produces an off-raster block instead of an error.  Catch it
        # here, where the boundary that caused it can still be named.
        needed = _required_duration(ready.events, opts)
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
                 'a different System -- please report it with the tree that produced it'],
            )
            raise CompileError(msg) from err

    # Named sources, so a conflict says *who* claimed the key twice rather than "already"/"also".
    defs = merge_definitions({
        'the sequence name': {'Name': name or root.tag or 'seqcraft'},
        'the geometry': geometry.definitions() if geometry is not None else {},
        'the definitions= argument': definitions or {},
    })

    out = CompiledSequence(
        seq=seq,
        system=system,
        regime=regime,
        report=Report(tuple(issues), subject=defs['Name']),
        origins=tuple(origins),
        definitions=defs,
        tree_duration_s=total,
        geometry=geometry,
    )
    out._verify(placed, targets)
    return out


def _common_path(paths: Sequence[tuple[str, ...]]) -> tuple[str, ...]:
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


# ------------------------------------------------------------------------------ the result
@dataclass(frozen=True)
class WriteResult:
    """What :meth:`CompiledSequence.write` produced."""

    path: Path
    sha256: str
    n_blocks: int
    duration_s: float
    sidecar: Path | None = None


@dataclass
class CompiledSequence:
    """
    A compiled sequence: the pypulseq object, the compile report, and provenance.

    Attributes
    ----------
    seq
        The :class:`pypulseq.Sequence`.  Yours to use directly; seqcraft never hides it.
    report
        Everything the compile found: every same-axis merge, every limit violation, every
        snapped time.
    origins
        One tag path per compiled block, so a block index traces back to the module that
        produced it.
    """

    seq: Any
    system: System
    regime: str
    report: Report
    origins: tuple[tuple[str, ...], ...]
    definitions: dict[str, Any]
    tree_duration_s: float
    geometry: Geometry | None = None
    _checked: Report | None = field(default=None, repr=False)

    # ------------------------------------------------------------------------ properties
    @property
    def n_blocks(self) -> int:
        """Number of pulseq blocks."""
        return len(self.seq.block_events)

    @property
    def duration_s(self) -> float:
        """Total duration, seconds.  ``Sequence.duration()`` returns a tuple; this does not."""
        return float(self.seq.duration()[0])

    def origin(self, block_index: int) -> tuple[str, ...]:
        """
        Return the tag path of the module that produced block `block_index`.

        Examples
        --------
        >>> import pypulseq as pp
        >>> import seqcraft as sc
        >>> system = sc.System.preset('generic_3t')
        >>> inner = sc.LogicBlock('spoiler')
        >>> _ = inner.add(0.0, pp.make_trapezoid('z', area=500.0, system=system.default))
        >>> out = sc.compile(sc.LogicBlock('tr').add(0.0, inner), system)
        >>> out.origin(0)
        ('tr', 'spoiler')
        """
        return self.origins[block_index]

    def __repr__(self) -> str:
        """One line: name, block count, duration, error and warning counts."""
        return (
            f'CompiledSequence({self.definitions.get("Name", "?")}, {self.n_blocks} blocks, '
            f'{self.duration_s:.3f} s, {len(self.report.errors)} errors, '
            f'{len(self.report.warnings)} warnings)'
        )

    # -------------------------------------------------------------------------- checking
    def _verify(
        self,
        placed: Sequence[_Placed],
        targets: dict[int, float] | None = None,
    ) -> None:
        """
        Assert the compile preserved what the tree meant.

        Four invariants, all cheap, each catching a class of compiler bug no individual test case
        would:

        **Total duration** must match the tree's.  Guards a boundary merge that dropped time
        rather than dropping a cut.

        **Zeroth moment (m0)** per axis must match the sum over the flattened events.  A split
        that lost a tail, or a merge that dropped a piece, changes it.

        **First moment (m1)** per axis, referenced to the start of the sequence.  m0 is exactly
        the quantity that survives a *time shift* -- an event moved by a whole raster leaves it
        untouched -- so m0 alone cannot see a gradient that plays at the wrong moment.  m1 can,
        because a piece of area ``A`` displaced by ``dt`` changes it by ``A * dt``.

        **Label addresses** must match the fold of the tree's labels.  The duplicate-address
        check in :meth:`check` only fires when two addresses *collide*; an addressing shifted by
        one but still unique passes it, which is exactly how mis-retimed labels used to escape.
        Skipped when a label has no ADC after it, since its placement is then reported as a
        warning rather than defined.
        """
        extra: list[Issue] = []
        got = self.duration_s
        # Tolerance scales with the total, because it has to: float64 resolves about 4 ns at
        # 20 000 s, so demanding nanosecond agreement on a long acquisition would flag every
        # sequence over an hour -- and print two identical-looking numbers while doing it.
        tolerance = max(EPS, 1e-12 * self.tree_duration_s)
        if abs(got - self.tree_duration_s) > tolerance:
            extra.append(
                Issue(
                    'duration',
                    'sequence',
                    f'compiled duration {got * 1e6:.1f} us differs from the tree total '
                    f'{self.tree_duration_s * 1e6:.1f} us',
                    'error',
                )
            )

        want: dict[str, float] = {}
        # Tolerance is scaled by the total area *traversed*, not by the net.  A readout and its
        # prephaser very nearly cancel, so a relative tolerance on the net would demand exactness
        # that float summation over thousands of pieces cannot deliver -- while a tolerance on the
        # total still catches a whole lost lobe, which is what this is for.
        magnitude: dict[str, float] = {}
        for p in placed:
            if p.kind in _GRAD:
                channel = p.event.channel
                area = ev.pwl_moment(*ev.knots_of(p.event, p.node_t), 0)
                want[channel] = want.get(channel, 0.0) + area
                magnitude[channel] = magnitude.get(channel, 0.0) + abs(area)
        # m1 referenced to the start of the sequence, computed **exactly** on both sides.
        #
        # Both are sums of per-piece moments, which is legitimate because a moment is linear in
        # the waveform: the moment of a sum is the sum of the moments, whether the pieces overlap
        # (the tree, where they superpose) or abut (the compiled blocks, where they concatenate).
        # No union of knots needs building.
        want_m1: dict[str, float] = {}
        for p in placed:
            if p.kind in _GRAD:
                want_m1[p.event.channel] = want_m1.get(p.event.channel, 0.0) + ev.pwl_moment(
                    *ev.knots_of(p.event, p.node_t), 1
                )
        actual_m1 = self.moments(order=1)

        if want:
            actual = self.moments()
            for axis, expected in want.items():
                scale = max(magnitude[axis], 1.0)
                if abs(actual.get(axis, 0.0) - expected) > 1e-6 * scale:
                    extra.append(
                        Issue(
                            'moment',
                            f'axis {axis}',
                            f'compiled m0 {actual.get(axis, 0.0):.6g} 1/m differs from the tree '
                            f'sum {expected:.6g} 1/m -- a split or a merge lost area',
                            'error',
                        )
                    )
                # Both sides are exact closed forms, so the only error is float summation over the
                # pieces.  Scaled by area *traversed* times the horizon, for the same reason m0's
                # is scaled by area traversed: a readout and its prephaser nearly cancel, so a
                # tolerance on the net would demand more than float64 can carry.  A lobe displaced
                # by one raster changes m1 by area * 10 us, comfortably above this.
                m1_scale = max(scale * max(self.tree_duration_s, 1e-3), 1.0)
                if abs(actual_m1.get(axis, 0.0) - want_m1.get(axis, 0.0)) > 1e-9 * m1_scale:
                    extra.append(
                        Issue(
                            'moment',
                            f'axis {axis}',
                            f'compiled m1 {actual_m1.get(axis, 0.0):.6g} s/m differs from the '
                            f'tree sum {want_m1.get(axis, 0.0):.6g} s/m -- a gradient plays at '
                            f'the wrong time, which m0 cannot see',
                            'error',
                        )
                    )

        extra.extend(self._address_issues(placed, targets))
        if extra:
            object.__setattr__(self, 'report', Report(
                (*self.report.issues, *extra), subject=self.report.subject
            ))

    def _address_issues(
        self,
        placed: Sequence[_Placed],
        targets: dict[int, float] | None,
    ) -> list[Issue]:
        """
        Check every ADC's compiled label state against the fold of the tree's labels.

        The check the duplicate-address test cannot do.  A labelling shifted by one readout stays
        unique, so ``_label_issues`` sees nothing wrong -- and that is precisely the shape of the
        bug retiming exists to prevent, so it needs an invariant of its own rather than trust.
        """
        if targets is None:
            return []
        n_labels = sum(1 for p in placed if p.kind in _LABEL)
        if not n_labels:
            return []
        # An orphan label's placement is reported as a warning rather than defined, so its effect
        # is not predictable from a time fold and this check would be comparing against a guess.
        if any(i for i, p in enumerate(placed) if p.kind in _LABEL and i not in targets):
            return []

        expected = _expected_addresses(placed, targets)
        try:
            got = self.seq.evaluate_labels(evolution='adc')
        except (AttributeError, ValueError, IndexError):  # pragma: no cover - older pypulseq
            return []
        if not got or not expected:
            return []

        out: list[Issue] = []
        for key in sorted({k for state in expected for k in state}):
            series = np.atleast_1d(np.asarray(got.get(key, 0)))
            series = np.resize(series, len(expected)) if series.size else np.zeros(len(expected))
            for index, state in enumerate(expected):
                if int(series[index]) != int(state.get(key, 0)):
                    out.append(Issue(
                        'address',
                        f'adc {index}',
                        f'label {key} is {int(series[index])} on readout {index} but the tree '
                        f'implies {int(state.get(key, 0))} -- a label reached the wrong readout',
                        'error',
                    ))
                    break       # one report per key; a shift corrupts every later readout too
        return out

    def moments(self, order: int = 0) -> dict[str, float]:
        """
        Return the whole-sequence gradient moment per axis, integrated from the compiled blocks.

        Parameters
        ----------
        order
            ``0`` for area in 1/m, ``1`` for s/m, ``2`` for s^2/m.  Referenced to the start of
            the sequence.

        Notes
        -----
        Integrated from each block's exact knots, not from raster samples.  The difference is
        not cosmetic: sampling is exact for ``order == 0`` and only then, and an arbitrary
        gradient's samples sit at raster *centres*, so a raster-sampled m0 quietly matched a
        raster-sampled tree even when the compiled waveform had moved.
        """
        out: dict[str, float] = dict.fromkeys(_AXES, 0.0)
        t = 0.0
        # Block IDs are 1-based, and block_durations is a dict keyed by them, not a list.
        for index in sorted(self.seq.block_events):
            block = self.seq.get_block(index)
            for axis in _AXES:
                grad = getattr(block, f'g{axis}', None)
                if grad is not None:
                    out[axis] += ev.pwl_moment(*ev.knots_of(grad, t), order)
            t += float(self.seq.block_durations[index])
        return out

    def check(self, *, allow_timing: Sequence[str] = ('TotalDuration',)) -> Report:
        """
        Run every post-compile check and return one report.

        Combines the compile report with ``Sequence.check_timing`` and label-address
        uniqueness.

        Parameters
        ----------
        allow_timing
            Substrings of ``check_timing`` messages to downgrade to information.  Defaults to
            the ``TotalDuration`` float-equality artifact, which pypulseq emits even on
            pulseq's own approved reference files.

        Examples
        --------
        >>> import pypulseq as pp
        >>> import seqcraft as sc
        >>> system = sc.System.preset('generic_3t')
        >>> lb = sc.LogicBlock('t')
        >>> _ = lb.add(0.0, pp.make_trapezoid('x', area=100.0, system=system.default))
        >>> sc.compile(lb, system).check().ok
        True
        """
        if self._checked is not None:
            return self._checked
        issues = list(self.report.issues)
        ok, errors = self.seq.check_timing()
        if not ok:
            for line in errors:
                text = str(line).strip()
                allowed = any(token in text for token in allow_timing)
                issues.append(Issue('timing', 'sequence', text, 'info' if allowed else 'error'))
        issues.extend(self._label_issues())
        issues.extend(self._event_size_issues())
        out = Report(tuple(issues), subject=self.report.subject, values={
            'n_blocks': self.n_blocks,
            'duration_s': self.duration_s,
        })
        object.__setattr__(self, '_checked', out)
        return out

    def _event_size_issues(self) -> list[Issue]:
        """
        Check every ADC and RF event against the interpreter's per-event sample limits.

        These limits live in ``Opts`` as ``adc_samples_limit`` and ``rf_samples_limit`` and default
        to ``0``, pypulseq's "no limit".  Nothing checked them until a 67 388-sample spiral readout
        reached a scanner, which refused the block with ``fRTEBFinish() failed for block type:
        ArbX ArbY ADC`` -- a message that names the block type and says nothing about samples.

        The limit is the vendor interpreter's, not the amplifier's, so it has to be set from the
        installation: ``System.preset`` puts 8192 on the Siemens entries, which is the common value.
        A readout longer than one event's worth has to be split into several ADCs, at the cost of
        ``adc_dead_time`` between them.
        """
        out: list[Issue] = []
        limits = self.system.limits(self.regime)
        for kind, attribute, count in (
            ('adc', 'adc_samples_limit', lambda block: int(block.adc.num_samples)),
            ('rf', 'rf_samples_limit', lambda block: int(np.size(block.rf.signal))),
        ):
            limit = int(getattr(limits, attribute, 0) or 0)
            if limit <= 0:
                continue
            worst = 0
            where = ''
            for index in sorted(self.seq.block_events):
                block = self.seq.get_block(index)
                if getattr(block, kind, None) is None:
                    continue
                samples = count(block)
                if samples > worst:
                    worst, where = samples, f'block {index} ({self.origin(index)})'
            if worst > limit:
                out.append(Issue(
                    f'{kind}_samples_limit',
                    where,
                    f'{worst} {kind.upper()} samples in one event, above the {limit} the '
                    f'interpreter accepts; split it into '
                    f'{-(-worst // limit)} events or lengthen the dwell',
                    'error',
                ))
        return out

    def _label_issues(self) -> list[Issue]:
        """
        Check that no two imaging ADCs write the same k-space address.

        The highest-value check available on a finished sequence: a duplicate address means two
        readouts landing in the same place, which catches a wrong slice order, an off-by-one
        partial-Fourier start and a mis-nested loop from one assertion.
        """
        try:
            labels = self.seq.evaluate_labels(evolution='adc')
        except (AttributeError, ValueError, IndexError):  # pragma: no cover - older pypulseq
            return []
        keys = [k for k in _ADDRESS_KEYS if k in labels]
        if not keys:
            return []
        arrays = [np.atleast_1d(np.asarray(labels[k])) for k in keys]
        n = max(len(a) for a in arrays)
        arrays = [np.resize(a, n) for a in arrays]
        skip = np.zeros(n, dtype=bool)
        for flag in ('NOISE', 'REF', 'NAV'):
            if flag in labels:
                skip |= np.resize(np.atleast_1d(np.asarray(labels[flag])).astype(bool), n)
        addresses = [tuple(int(a[i]) for a in arrays) for i in range(n) if not skip[i]]
        duplicates = len(addresses) - len(set(addresses))
        if not duplicates:
            return []
        seen: set[tuple[int, ...]] = set()
        first: tuple[int, ...] = ()
        for address in addresses:
            if address in seen:
                first = address
                break
            seen.add(address)
        return [
            Issue(
                'label',
                'sequence',
                f'{duplicates} imaging ADC(s) repeat a k-space address; first repeat is '
                f'{dict(zip(keys, first))} -- two readouts are writing the same location',
                'error',
            )
        ]

    # ---------------------------------------------------------------------------- output
    def kspace(self) -> dict[str, np.ndarray]:
        """
        Return the k-space trajectory, in 1/m.

        Returns
        -------
        dict
            ``k_adc`` (3 x n_samples, at the ADC sample times), ``t_adc``, ``k`` (dense),
            ``t_k``, ``t_excitation``, ``t_refocusing``.

        Notes
        -----
        ``calculate_kspacePP`` returns its tuple in a different order from
        ``calculate_kspace``; getting that wrong silently swaps the trajectory for its
        timebase.
        """
        k_adc, t_adc, k, t_k, t_exc, t_refoc = self.seq.calculate_kspacePP()[:6]
        return {
            'k_adc': np.asarray(k_adc),
            't_adc': np.asarray(t_adc),
            'k': np.asarray(k),
            't_k': np.asarray(t_k),
            't_excitation': np.asarray(t_exc),
            't_refocusing': np.asarray(t_refoc),
        }

    def pns(self, hardware: SimpleNamespace | None = None) -> Report:
        """
        Predict peripheral nerve stimulation against a gradient hardware model.

        Parameters
        ----------
        hardware
            Defaults to the model attached to the :class:`~seqcraft.core.system.System`.
            :func:`~seqcraft.core.system.synthetic_hardware` provides a vendor-free stand-in
            for CI; it is **not** a real scanner and must never be used to clear a sequence for
            human scanning.
        """
        model = hardware if hardware is not None else self.system.hardware
        if model is None:
            return Report((
                Issue(
                    'pns',
                    'sequence',
                    'no gradient hardware model attached; call System.with_hardware() with '
                    'load_hardware() or synthetic_hardware()',
                    'info',
                ),
            ))
        ok, pns_norm, _components, _t = self.seq.calculate_pns(model, do_plots=False)
        peak = float(np.max(pns_norm)) if np.size(pns_norm) else 0.0
        note = (
            ' (synthetic model -- not valid for clearing a human scan)'
            if getattr(model, 'is_synthetic', False)
            else ''
        )
        return Report(
            (
                Issue(
                    'pns',
                    'sequence',
                    f'peak PNS {peak * 100:.0f}% of the stimulation limit{note}',
                    'info' if ok else 'error',
                ),
            ),
            values={'peak_pns_fraction': peak},
        )

    def write(self, path: str | Path, *, sidecar: bool = True) -> WriteResult:
        """
        Write the ``.seq`` file, and by default a JSON provenance sidecar beside it.

        Takes no geometry, matrix or FOV arguments: everything written comes from what was
        compiled, so the file's metadata cannot disagree with what it plays.

        Parameters
        ----------
        path
            Destination ``.seq`` path.
        sidecar
            Also write ``<path>.json`` recording versions, git state, the definitions, the
            achieved duration and the file's sha256.
        """
        import hashlib  # noqa: PLC0415
        from pathlib import Path as _Path  # noqa: PLC0415

        target = _Path(path)
        for key, value in self.definitions.items():
            self.seq.set_definition(key, value)
        self.seq.set_definition('TotalDuration', self.duration_s)
        self.seq.write(str(target))
        digest = hashlib.sha256(target.read_bytes()).hexdigest()

        side: _Path | None = None
        if sidecar:
            from ..provenance import write_sidecar  # noqa: PLC0415

            side = write_sidecar(
                target,
                {
                    'definitions': {k: _jsonable(v) for k, v in self.definitions.items()},
                    'system': self.system.params(),
                    'regime': self.regime,
                    'n_blocks': self.n_blocks,
                    'duration_s': self.duration_s,
                    'sha256': digest,
                    'issues': [
                        {'kind': i.kind, 'where': i.where, 'message': i.message,
                         'severity': i.severity}
                        for i in self.report.issues
                    ],
                },
            )
        return WriteResult(
            path=target,
            sha256=digest,
            n_blocks=self.n_blocks,
            duration_s=self.duration_s,
            sidecar=side,
        )


def _jsonable(value: Any) -> Any:
    """Convert numpy scalars and arrays to plain Python, for the sidecar."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value
