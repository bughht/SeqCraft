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
>>> out = sc.compile(lb, opts)
>>> out.report.ok, len(out.report.warnings)
(True, 0)
>>> out.n_blocks
1
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import pypulseq as pp

from ..design.events import GRADIENT_KINDS
from ..design.timing import EPS, Raster
from ..errors import CompileError, format_error
from ..report import Issue, Report
from ..result import CompiledSequence
from .boundaries import (
    check_exclusive,
    find_boundaries,
    label_order_conflict,
    label_targets,
    orphan_label_issues,
)
from .definitions import merge_definitions
from .emission import emit_blocks
from .placement import place_events
from .verification import (
    require_valid_contract,
    verify_against_tree,
    verify_placed_events,
)

if TYPE_CHECKING:
    from collections.abc import Mapping

    from pypulseq.opts import Opts

    from ..design.logic import LogicBlock

__all__ = ['compile_sequence']


def compile_sequence(  # noqa: C901, PLR0912, PLR0915
    root: LogicBlock,
    opts: Opts,
    *,
    name: str = '',
    definitions: Mapping[str, Any] | None = None,
) -> CompiledSequence:
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
    CompiledSequence
        The pypulseq sequence, the compile report, and per-block provenance.

    Raises
    ------
    CompileError
        Two RF or ADC events overlap in time, an absolute start is negative, or a block
        boundary would have to fall inside a gradient an ADC is sampling.

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

    Amplitude and slew violations are *reported*, not raised.  Call
    :meth:`CompiledSequence.check` then
    :meth:`~seqcraft.report.Report.raise_if_failed` to stop on them.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
    ...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> lb = sc.LogicBlock('t').add(0.0, pp.make_delay(1e-3))
    >>> out = compile_sequence(lb, opts, definitions={'FOV': [0.25, 0.25, 0.005]})
    >>> round(out.duration_s * 1e3, 1), out.definitions['FOV']
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

    issues: list[Issue] = [
        Issue('raster', p.where, f'start {p.res_start * 1e6:.4f} us snapped to '
              f'{raster.nearest(p.res_start) * 1e6:.1f} us', 'warning')
        for p in placed
        if p.kind in ('rf', 'adc') and not raster.holds(p.res_start)
    ]

    max_block = raster.dt * 2**24
    edges = find_boundaries(placed, total, raster, max_block)

    # Labels are assigned at their *target ADC's* time rather than their own, so a boundary
    # landing between a label and the readout it addresses cannot change which ADC sees it.
    targets = label_targets(placed)
    conflict = label_order_conflict(placed, targets)
    if conflict is not None:
        raise conflict
    issues.extend(orphan_label_issues(placed, targets))

    seq = pp.Sequence(system=opts)
    origins = emit_blocks(seq, edges, placed, targets, opts, raster, issues)

    # Named sources, so a conflict says *who* claimed the key twice rather than "already"/"also".
    defs = merge_definitions({
        'the sequence name': {'Name': name or root.tag or 'seqcraft'},
        'the definitions= argument': dict(definitions or {}),
    })

    out = CompiledSequence(
        seq=seq,
        opts=opts,
        report=Report(tuple(issues), subject=defs['Name']),
        origins=tuple(origins),
        definitions=defs,
        tree_duration_s=total,
    )
    # The last compile stage: check the produced sequence against the tree it came from.  It runs
    # here, on the object just built, rather than as a method on it -- it reads the private IR,
    # which is the compiler's, not the result's.
    verify_against_tree(
        placed,
        targets,
        duration_s=out.duration_s,
        tree_duration_s=total,
        moments=out.moments,
        label_states=lambda: seq.evaluate_labels(evolution='adc'),
    )
    return out
