"""
Where the blocks are cut.

A pulseq block holds at most one RF and one ADC, must last a whole number of block rasters, and
must join its neighbours continuously.  This stage decides, from the placed events alone, the
sorted list of times at which the sequence is cut -- before a single pulseq object exists.

Two independent properties of an event decide where a boundary may and must go, and keeping them
apart is what makes triggers work.

**Indivisible** -- no boundary may fall strictly inside its reservation.  An RF's dead time and
ringdown, an ADC's window and trailing dead time, a trigger's pulse: each is one hardware action
rather than a waveform, so cutting it is meaningless.

**Exclusive** -- a block holds at most one, so a boundary is *required* between two.  RF and ADC
only.  Triggers are ``TRIGGERS`` extensions and pulseq accepts several per block, so demanding a
boundary between two of them would invent a constraint the hardware does not have.

Putting one boundary in the gap between each pair of consecutive exclusive reservations guarantees
at most one RF and one ADC per block **by construction**, with nothing left to check afterwards.
Gradients constrain nothing: they follow wherever the boundaries land, split and summed by
:mod:`~seqcraft.compiler.legalization`.

Labels are the other half of this module, and they are here rather than in emission because the
answer they need -- which readout does this label address? -- has to be independent of where the
boundaries turn out to be.  :func:`_label_targets` assigns each label to its target ADC's time, and
:func:`_label_order_conflict` rejects the groups pulseq cannot express.
"""

from __future__ import annotations

import bisect
import math
from typing import TYPE_CHECKING

from ..design.events import GRADIENT_KINDS, LABEL_KINDS
from ..design.logic import BARRIER
from ..design.timing import EPS
from ..errors import CompileError, format_error
from ..report import Issue
from .model import EXCLUSIVE_KINDS, INDIVISIBLE_KINDS, PlacedEvent

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Sequence

    from ..design.timing import Raster

__all__ = [
    'check_exclusive',
    'find_boundaries',
    'label_order_conflict',
    'label_targets',
    'orphan_label_issues',
]


def check_exclusive(placed: Sequence[PlacedEvent]) -> None:
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


def label_targets(placed: Sequence[PlacedEvent]) -> dict[int, float]:
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
    :func:`label_order_conflict`: separated into different blocks, two labels are applied in
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
        if p.kind not in LABEL_KINDS:
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


def label_order_conflict(
    placed: Sequence[PlacedEvent],
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
    groups: dict[tuple[float, str], list[PlacedEvent]] = {}
    for i, p in enumerate(placed):
        if p.kind not in LABEL_KINDS:
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




def orphan_label_issues(
    placed: Sequence[PlacedEvent],
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
        if p.kind in LABEL_KINDS and i not in targets:
            out.append(Issue(
                'label',
                p.where,
                f'label {getattr(p.event, "label", "?")} at {p.res_start * 1e6:.1f} us has no '
                f'ADC after it, so it addresses no readout; it may land in the preceding '
                f'readout\'s block and change that address',
                'warning',
            ))
    return out


def _covering(placed: Sequence[PlacedEvent], t: float) -> PlacedEvent | None:
    """Return an indivisible event whose reservation strictly contains `t`, for error messages."""
    for p in placed:
        if p.kind in INDIVISIBLE_KINDS and p.res_start + EPS < t < p.res_end - EPS:
            return p
    return None


def _barrier_conflict(
    barrier: PlacedEvent,
    at: float,
    placed: Sequence[PlacedEvent],
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


def _gap_blocked(a_end: float, b_start: float, placed: Sequence[PlacedEvent]) -> CompileError:
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


def find_boundaries(
    placed: Sequence[PlacedEvent],
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
        ((p.res_start, p.res_end) for p in placed if p.kind in EXCLUSIVE_KINDS),
        key=lambda s: s[0],
    )
    grad_spans = _Spans((p.start, p.end) for p in placed if p.kind in GRADIENT_KINDS)
    indivisible = _Spans(
        (p.res_start, p.res_end) for p in placed if p.kind in INDIVISIBLE_KINDS
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
        if p.kind in GRADIENT_KINDS:
            candidates.add(raster.floor(p.start))
            candidates.add(raster.ceil(p.end))
        elif p.kind in INDIVISIBLE_KINDS:
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
        # A separate name from the barrier loop's `at` above: that one is where a barrier snapped
        # to and is always a time, this one is a search that can come back empty.
        cut = _gap_boundary(a_end, b_start, raster, cuttable)
        if cut is None:
            raise _gap_blocked(a_end, b_start, placed)
        extra.append(cut)
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
