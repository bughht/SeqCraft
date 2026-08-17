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

What the compiler needs from the scanner
----------------------------------------
One :class:`pypulseq.Opts`, and eight fields of it: the two amplitude limits, the four rasters,
the ADC dead time and the RF ringdown.  It takes that object directly rather than any seqcraft
wrapper, so the same thing that configures ``pp.make_trapezoid`` configures the compile.

Examples
--------
>>> import pypulseq as pp
>>> import seqcraft as sc
>>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
>>> gentle = {'area': 100.0, 'duration': 2e-3, 'rise_time': 200e-6, 'system': opts}
>>> lb = sc.LogicBlock('demo')
>>> _ = lb.add(0.0, pp.make_trapezoid('x', **gentle))
>>> _ = lb.add(0.0, pp.make_trapezoid('y', **gentle))   # different axis: nothing to report
>>> seq = sc.compile(lb, opts)
>>> type(seq).__module__.startswith('pypulseq')
True
>>> len(seq.block_events)
1

What a compile returns
----------------------
A :class:`pypulseq.Sequence`, and nothing else.  Not a wrapper, not a result object, not a pair of
a sequence and a report.  ``seq.write(...)``, ``seq.plot()``, ``seq.block_events`` and
``seq.definitions`` are pypulseq's own and are what the rest of the ecosystem already reads.

That is only tenable because **every legality failure raises**.  A returned object carrying a
report is a way of not noticing: it writes a ``.seq`` the console refuses an hour later, and the
explanation is on an object nobody looked at.  What the compile *did* rather than refused -- summed
two gradients on an axis, resampled one onto the raster -- is a
:class:`~seqcraft.errors.SeqCraftWarning`, one per category, so the standard ``warnings`` machinery
decides what happens to it.
"""

from __future__ import annotations

import collections
import warnings
from typing import TYPE_CHECKING, Any

import pypulseq as pp

from ..design.events import GRADIENT_KINDS
from ..design.timing import EPS, Raster
from ..errors import SeqCraftWarning, format_error
from .boundaries import (
    check_exclusive,
    find_boundaries,
    label_order_conflict,
    label_targets,
    orphan_label_notes,
)
from .emission import emit_blocks
from .errors import CompileError, DefinitionConflict
from .legalization import legalize_blocks
from .placement import place_events
from .verification import (
    _sequence_moments,
    check_event_sizes,
    check_label_addresses,
    require_valid_contract,
    verify_against_tree,
    verify_placed_events,
)

if TYPE_CHECKING:
    from collections.abc import Mapping

    from pypulseq.opts import Opts

    from ..design.logic import LogicBlock

__all__ = ['compile_sequence']

#: Substrings of a ``check_timing`` complaint that stay informational rather than failing the
#: compile.  ``TotalDuration`` is a float-equality artifact pypulseq emits even on pulseq's own
#: approved reference files, so treating it as an error would fail every sequence.
#:
#: Deliberately not a parameter.  It was ``check(allow_timing=...)`` and no caller ever passed it,
#: which makes it a knob whose only effect is to let a real timing failure through.
_ALLOWED_TIMING = ('TotalDuration',)

#: What each note category is called in its warning, singular and plural.  Every category the
#: compiler can produce appears here, so adding one without naming it is a visible omission rather
#: than a bare key -- and both forms are spelled out because the plural is not always on the last
#: word ("blocks over the ... limit").
_WARNING_TEXT = {
    'merge': ('same-axis gradient merge', 'same-axis gradient merges'),
    'norm': ('block over the vector-norm limit (legal on real amplifiers)',
             'blocks over the vector-norm limit (legal on real amplifiers)'),
    'orphan_label': ('label with no ADC after it', 'labels with no ADC after them'),
    'resample': ('gradient resampled onto the raster', 'gradients resampled onto the raster'),
    'snap': ('reservation snapped to the block raster',
             'reservations snapped to the block raster'),
}

#: How many distinct sites a warning names before it starts counting.
_WARNING_SITES = 5


def _warn(notes: dict[str, list[str]]) -> None:
    """
    Emit one aggregated :class:`~seqcraft.errors.SeqCraftWarning` per non-empty category.

    **Aggregation is not cosmetic.**  Python's default filter shows a warning once per unique
    ``(message, category, module, lineno)``, so one ``warn`` per merge would print the first and
    silently swallow the other eleven -- strictly worse than the count a report gave.  One warning
    per category carries the whole picture in one line.

    Identical entries are counted rather than repeated, because the interesting number is how many
    *sites* merged, not how many TRs the loop ran for: a 64-line acquisition merges the same two
    modules 64 times and that is one fact, not sixty-four.

    ``stacklevel=2`` so the warning points at the caller's ``sc.compile(...)`` rather than at the
    compiler's guts.
    """
    for category, entries in sorted(notes.items()):
        if not entries:
            continue
        counted = collections.Counter(entries)
        shown = ', '.join(
            site if n == 1 else f'{site} x{n}'
            for site, n in counted.most_common(_WARNING_SITES)
        )
        extra = len(counted) - _WARNING_SITES
        names = _WARNING_TEXT.get(category, (category, category))
        warnings.warn(
            f'{len(entries)} {names[len(entries) != 1]}: {shown}'
            + (f' (+{extra} more sites)' if extra > 0 else ''),
            SeqCraftWarning,
            stacklevel=3,        # _warn -> compile_sequence -> the caller
        )


def compile_sequence(  # noqa: C901, PLR0912, PLR0915
    root: LogicBlock,
    opts: Opts,
    *,
    name: str = '',
    definitions: Mapping[str, Any] | None = None,
) -> pp.Sequence:
    """
    Turn a logic-block tree into a :class:`pypulseq.Sequence`.

    Parameters
    ----------
    root
        The tree to compile.  Usually built by adding module outputs to one
        :class:`~seqcraft.design.logic.LogicBlock`.
    opts
        The scanner: rasters, dead times and limits.  A plain :class:`pypulseq.Opts`, the same
        object passed to pypulseq's ``make_*`` functions.  This is the ceiling the *combined*
        waveform has to respect; a part designed against a derated limit is designed against a
        second ``Opts`` (:func:`seqcraft.opts.derate`), which does not change what is validated
        here.
    name
        Sequence name for the definitions.  Defaults to ``root.tag``.
    definitions
        ``[DEFINITIONS]`` entries for the written file -- pulseq's own vocabulary, exactly as it
        will be written: ``FOV``, ``SliceThickness``, ``kSpaceCenterLine``, ``TE``, ``TR``, and
        anything else your acquisition wants recorded.  Merged with the sequence name under a
        collision check, so two sources claiming the same key with different values raises rather
        than one silently winning.

    Returns
    -------
    pypulseq.Sequence
        Yours to use directly.  The definitions you passed are already set on it, along with
        ``Name`` and ``TotalDuration``, so ``seq.write(path)`` needs nothing else.

    Raises
    ------
    CompileError
        Two RF or ADC events overlap in time, an absolute start is negative, a gradient starts off
        the gradient raster, no boundary can be cut in the gap between two exclusive events, a
        block boundary would fall inside a gradient an ADC is sampling, two ADCs write the same
        k-space address, or ``check_timing`` fails.
    HardwareLimitError
        The *summed* waveform exceeds ``max_grad`` or ``max_slew`` on an axis, or an ADC or RF
        event exceeds the interpreter's per-event sample limit.
    DefinitionConflict
        ``name=`` and ``definitions['Name']`` disagree.
    CompilerContractError
        The compiled sequence does not match the tree.  A compiler bug; report it with the tree.

    Notes
    -----
    **There is no geometry argument.**  There was, and it existed to call ``definitions()`` on a
    ``Geometry`` and merge the eight keys that came back -- which is what ``definitions=`` already
    does, for any source, with no dataclass in between.  FOV, matrix and slice order are decisions
    about the scan you are running; the compiler turns a tree into legal pulseq blocks and is
    indifferent to why the tree looks the way it does.  Hand it the keys::

        sc.compile(tree, opts, definitions=my_geometry.definitions())

    ``salvage/geometry.py`` holds the dataclass that used to be here, standalone, for when a module
    library wants one again.

    **Every legality failure raises**, including amplitude and slew: there is no legal sequence to
    hand back, and a returned object carrying a report is a way of not noticing.  What the compile
    *did* -- summed two gradients on an axis, resampled one onto the raster -- is a
    :class:`~seqcraft.errors.SeqCraftWarning` instead.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
    ...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> lb = sc.LogicBlock('t').add(0.0, pp.make_delay(1e-3))
    >>> seq = compile_sequence(lb, opts, definitions={'FOV': [0.25, 0.25, 0.005]})
    >>> round(seq.duration()[0] * 1e3, 1), seq.definitions['FOV']
    (1.0, [0.25, 0.25, 0.005])
    """
    raster = Raster(float(opts.block_duration_raster), 'block')
    placed = place_events(root, opts)
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
    grad_raster = Raster(float(opts.grad_raster_time), 'gradient')
    off = [p for p in placed if p.kind in GRADIENT_KINDS and not grad_raster.holds(p.start)]
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
                'round the start time where it is computed, with '
                'seqcraft.Raster(opts.grad_raster_time).ceil(t)',
                'a start derived from "te - module.time_to_echo" lands off the raster whenever TE '
                'does, so quantise TE first',
            ],
        )
        raise CompileError(msg)

    check_exclusive(placed)

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

    notes: dict[str, list[str]] = {}
    snapped = [
        f'{p.where} {p.res_start * 1e6:.4f} -> {raster.nearest(p.res_start) * 1e6:.1f} us'
        for p in placed
        if p.kind in ('rf', 'adc') and not raster.holds(p.res_start)
    ]
    if snapped:
        notes['snap'] = snapped

    max_block = raster.dt * 2**24
    edges = find_boundaries(placed, total, raster, max_block)

    # Labels are assigned at their *target ADC's* time rather than their own, so a boundary
    # landing between a label and the readout it addresses cannot change which ADC sees it.
    targets = label_targets(placed)
    conflict = label_order_conflict(placed, targets)
    if conflict is not None:
        raise conflict
    orphans = orphan_label_notes(placed, targets)
    if orphans:
        notes['orphan_label'] = orphans

    legalized = legalize_blocks(edges, placed, targets, opts, raster)
    for category, entries in legalized.notes:
        notes.setdefault(category, []).extend(entries)
    origins = [block.origin for block in legalized.blocks]

    # Legalization returns a complete immutable tuple so its cross-block contract can be checked.
    # Drain a private queue during emission so events already copied into PyPulseq are not retained
    # by that tuple for the rest of a long sequence build.
    ready = collections.deque(legalized.blocks)
    del legalized

    seq = pp.Sequence(system=opts)
    emit_blocks(seq, (ready.popleft() for _ in range(len(ready))))

    # The two questions only a built sequence can answer: does any one event exceed the
    # interpreter's sample limit, and do two imaging ADCs write the same k-space address.
    check_event_sizes(seq, opts, origins)
    check_label_addresses(seq)

    # pypulseq's own timing audit, run here rather than offered as a method: a `.seq` that fails
    # it is one the console will refuse, so there is nothing to hand back.
    timing_ok, complaints = seq.check_timing()
    if not timing_ok:
        fatal = [
            text for text in (str(line).strip() for line in complaints)
            if not any(token in text for token in _ALLOWED_TIMING)
        ]
        if fatal:
            msg = format_error(
                f'pypulseq\'s timing check found {len(fatal)} problem(s) in the compiled '
                f'sequence.',
                dict(enumerate(fatal[:5], start=1)),
                ['this is a compiler bug unless the tree contains raw events built against a '
                 'different Opts -- please report it with the tree that produced it'],
            )
            raise CompileError(msg)

    # Two sources, so a collision check is two lines rather than a merge algorithm.  It was a
    # named-source merge while `Geometry` was a third source; that dataclass is in salvage/ now,
    # and FOV and matrix arrive through `definitions=` like everything else.
    defs = dict(definitions or {})
    wanted = name or root.tag or 'seqcraft'
    if defs.get('Name', wanted) != wanted:
        raise DefinitionConflict(format_error(
            f'two sources set Name: {wanted!r} and {defs["Name"]!r}.',
            {'from name=/root.tag': wanted, 'from definitions=': defs['Name']},
            ['pass one or the other, not both'],
        ))
    defs['Name'] = wanted

    # This is what makes the returned pp.Sequence self-sufficient: the definitions are on it, so
    # nothing has to survive until write time to put them there.
    for key, value in defs.items():
        seq.set_definition(key, value)
    seq.set_definition('TotalDuration', float(seq.duration()[0]))

    # The last compile stage: check the produced sequence against the tree it came from.  It reads
    # the private IR, so it is a compile stage that happens to run last rather than an accessor
    # somebody has to remember to call.
    verify_against_tree(
        placed,
        targets,
        duration_s=float(seq.duration()[0]),
        tree_duration_s=total,
        moments=lambda order: _sequence_moments(seq, order),
        label_states=lambda: seq.evaluate_labels(evolution='adc'),
    )
    # Last, so that nothing warns about a compile that then failed.
    _warn(notes)
    return seq
