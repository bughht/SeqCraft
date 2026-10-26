"""Verification: the compiler checking its own work.

Two layers, and the difference is worth stating.

**Contracts** (:func:`verify_placed_events`, :func:`verify_ready_blocks`) protect *stage
boundaries*.  They ask whether one stage handed the next a well-formed intermediate
representation -- finite times, intervals the right way round, blocks that abut, at most one RF
and one ADC.  A failure is a compiler bug by definition, so it raises
:class:`CompilerContractError` rather than becoming a report.

**Semantics** (:func:`verify_against_tree`) asks the harder question: did the compile preserve
what the *tree* meant?  Duration, zeroth and first gradient moment per axis, and the label
address each ADC ends up with, all compared against the placed events the compile started from.
A failure here is a compiler bug too -- the tree was legal, since every legality check has already
passed -- so it raises :class:`CompilerContractError` as well.  All four invariants are measured
before anything is raised, because seeing every affected axis at once beats stopping at the first.

**Finished-sequence checks** (:func:`check_event_sizes`, :func:`check_rf_amplitude`,
:func:`check_label_addresses`) ask what only a built :class:`pypulseq.Sequence` can answer: does
any one event exceed the interpreter's sample limit, does any RF event's peak exceed the transmit
limit, and do two imaging ADCs write the same k-space address?  None is an internal invariant --
all three are constraints that can only be evaluated on the built sequence, so seqcraft raises on
them before returning it rather than handing back something it knows to be unusable: the first two
:class:`~seqcraft.errors.HardwareLimitError` and the third
:class:`~seqcraft.errors.CompileError`.

They live here rather than on what a compile returns because a check nobody has to call is a
check nobody calls.  ``pSeq_Base``'s ``get_report()`` printed its result and returned
``None``; its successor returned an object with a ``check()`` method, which is the same failure
one indirection later.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from ..design.events import (
    ADDRESS_KEYS,
    AXES,
    GRADIENT_KINDS,
    HANDLED_KINDS,
    LABEL_KINDS,
    knots_of,
    pwl_moment,
)
from ..design.timing import EPS
from ..errors import format_error
from .errors import CompileError, HardwareLimitError
from .model import (
    EXCLUSIVE_KINDS,
    PlacedEvent,
    PulseqReadyBlock,
    interval_duration,
    time_at_or_before,
    time_equal,
)

if TYPE_CHECKING:
    from collections.abc import Callable


@dataclass(frozen=True)
class ContractViolation:
    """One internal contract failure with a stable location and explanation."""

    location: str
    message: str

    def __str__(self) -> str:
        """Format the violation for an internal compiler exception."""
        return f'{self.location}: {self.message}'


class CompilerContractError(RuntimeError):
    """Raised when one private compiler stage breaks an internal IR contract."""

    def __init__(self, name: str, violations: Sequence[ContractViolation]) -> None:
        details = '\n'.join(f'- {violation}' for violation in violations)
        super().__init__(f'internal {name} contract violated:\n{details}')


def verify_placed_events(events: Sequence[PlacedEvent]) -> tuple[ContractViolation, ...]:
    """Return structural violations in an ordered placed-event sequence."""
    violations: list[ContractViolation] = []
    for index, placed in enumerate(events):
        location = f'placed event {index} ({placed.where})'
        times = (placed.node_t, placed.start, placed.end, placed.res_start, placed.res_end)
        if not all(math.isfinite(value) for value in times):
            violations.append(ContractViolation(location, 'all time fields must be finite'))
        if placed.kind not in HANDLED_KINDS:
            violations.append(ContractViolation(location, f'unclassified kind {placed.kind!r}'))
        if not time_at_or_before(placed.start, placed.end):
            violations.append(ContractViolation(location, 'active interval is reversed'))
        if not time_at_or_before(placed.res_start, placed.start):
            violations.append(ContractViolation(location, 'reservation starts after activity'))
        if not time_at_or_before(placed.end, placed.res_end):
            violations.append(ContractViolation(location, 'reservation ends before activity'))
        if not isinstance(placed.path, tuple) or not all(isinstance(part, str) for part in placed.path):
            violations.append(ContractViolation(location, 'source path must be a tuple of strings'))
    return tuple(violations)


def verify_ready_blocks(
    blocks: Sequence[PulseqReadyBlock],
    *,
    expected_first_index: int = 0,
    expected_start: float | None = None,
) -> tuple[ContractViolation, ...]:
    """Return structural violations in one complete sequence or a streamed block window."""
    violations: list[ContractViolation] = []
    previous_end = expected_start
    for offset, block in enumerate(blocks):
        expected_index = expected_first_index + offset
        location = f'ready block {block.index}'
        if block.index != expected_index:
            violations.append(ContractViolation(location, f'expected contiguous index {expected_index}'))
        if not all(math.isfinite(value) for value in (block.start, block.end, block.duration)):
            violations.append(ContractViolation(location, 'all time fields must be finite'))
        if not time_at_or_before(block.start, block.end):
            violations.append(ContractViolation(location, 'block interval is reversed'))
        if not time_equal(block.duration, interval_duration(block.start, block.end)):
            violations.append(ContractViolation(location, 'duration does not match start/end'))
        if previous_end is not None and not time_equal(previous_end, block.start):
            violations.append(ContractViolation(location, 'block does not abut its predecessor'))
        previous_end = block.end

        kinds = block.kinds
        for kind in EXCLUSIVE_KINDS:
            if kinds.count(kind) > 1:
                violations.append(ContractViolation(location, f'contains more than one {kind.upper()} event'))
        if not isinstance(block.events, tuple):
            violations.append(ContractViolation(location, 'events must be stored as a tuple'))
        paths = (*block.source_paths, block.origin)
        if not all(isinstance(path, tuple) for path in paths):
            violations.append(ContractViolation(location, 'provenance paths must be tuples'))
    return tuple(violations)


def require_valid_contract(
    name: str,
    violations: Sequence[ContractViolation],
) -> None:
    """Raise an internal error when a compiler stage violates its declared contract."""
    if not violations:
        return
    raise CompilerContractError(name, violations)


# ------------------------------------------------------------------ what only the sequence knows
def _sequence_moments(seq: Any, order: int = 0) -> dict[str, float]:
    """
    Return the whole-sequence gradient moment per axis, integrated from the compiled blocks.

    The compiled side of :func:`verify_against_tree`.  Deliberately *not* shared with
    :func:`seqcraft.analysis.moments`, which walks the tree: these are different iterations over
    different structures, and having the self-check call the same code as the thing it checks
    would make it compare a number with itself.

    Parameters
    ----------
    seq
        The built :class:`pypulseq.Sequence`.
    order
        ``0`` for area in 1/m, ``1`` for s/m, ``2`` for s^2/m, referenced to the sequence start.

    Notes
    -----
    Integrated from each block's exact knots, not from raster samples.  The difference is not
    cosmetic: sampling is exact for ``order == 0`` and only then, and an arbitrary gradient's
    samples sit at raster *centres*, so a raster-sampled m0 quietly matched a raster-sampled tree
    even when the compiled waveform had moved.
    """
    out: dict[str, float] = dict.fromkeys(AXES, 0.0)
    t = 0.0
    # Block IDs are 1-based, and block_durations is a dict keyed by them, not a list.
    for index in sorted(seq.block_events):
        block = seq.get_block(index)
        for axis in AXES:
            grad = getattr(block, f'g{axis}', None)
            if grad is not None:
                out[axis] += pwl_moment(*knots_of(grad, t), order)
        t += float(seq.block_durations[index])
    return out


def check_event_sizes(seq: Any, opts: Any, origins: Sequence[tuple[str, ...]] = ()) -> None:
    """
    Check every ADC and RF event against the interpreter's per-event sample limits.

    These limits live in ``Opts`` as ``adc_samples_limit`` and ``rf_samples_limit`` and default to
    ``0``, pypulseq's "no limit".  Nothing checked them until a 67 388-sample spiral readout was
    run, and in that run the block was refused with ``fRTEBFinish() failed for block type:
    ArbX ArbY ADC`` -- a message that names the block type and says nothing about samples.

    The limit is the vendor interpreter's, not the amplifier's, so it has to be set from the
    installation: :func:`seqcraft.scanner.opts.from_scanner` takes ``adc_samples_limit=`` for
    exactly this, and the examples in this repository use 8192.  A readout longer than one
    event's worth
    has to be split into several ADCs, at the cost of ``adc_dead_time`` between them.

    Raises
    ------
    HardwareLimitError
        Naming the worst event, where it came from, and how many events it needs splitting into.
    """
    for kind, attribute, count in (
        ('adc', 'adc_samples_limit', lambda block: int(block.adc.num_samples)),
        ('rf', 'rf_samples_limit', lambda block: int(np.size(block.rf.signal))),
    ):
        limit = int(getattr(opts, attribute, 0) or 0)
        if limit <= 0:
            continue
        worst = 0
        where = ''
        for index in sorted(seq.block_events):
            block = seq.get_block(index)
            if getattr(block, kind, None) is None:
                continue
            samples = count(block)
            if samples > worst:
                # Block ids are 1-based; `origins` is a list in emission order, so it is not.
                path = origins[index - 1] if 0 < index <= len(origins) else ()
                worst = samples
                where = f'block {index} ({".".join(path) or "?"})'
        if worst > limit:
            msg = format_error(
                f'{worst} {kind.upper()} samples in one event, above the {limit} the '
                f'interpreter accepts.',
                {'from': where, 'limit': f'{kind}_samples_limit = {limit}'},
                [
                    f'split it into {-(-worst // limit)} events, at the cost of '
                    f'adc_dead_time between them',
                    'or lengthen the dwell, which lowers the readout bandwidth',
                ],
            )
            raise HardwareLimitError(msg)


def check_rf_amplitude(seq: Any, opts: Any, origins: Sequence[tuple[str, ...]] = ()) -> None:
    """
    Check every emitted RF event's peak amplitude against ``opts.max_b1``.

    The **backstop** for the peak-B1 contract.  Every RF-producing module checks its own pulse at
    design time, where it can say what to change; this one measures the waveform that is actually
    being emitted, so it also catches a raw ``pp.make_sinc_pulse`` added straight to a
    :class:`~seqcraft.LogicBlock`, a future module that forgets its own check, and any path that
    reaches a block without passing one.

    The comparison is literal, so ``max_b1 = +inf`` is how a caller designs without a transmit
    limit -- and it is the only value that does.  Zero, a negative, ``NaN`` and an absent limit
    are each refused as unusable rather than treated as unlimited; see
    :func:`~seqcraft.modules._support.check_peak_b1`, which the modules use and which this agrees
    with.  ``pp.Opts()`` **defaults it to 851.52 Hz (20 uT)**, so the limit is live unless a
    caller has deliberately raised it.

    Raises
    ------
    HardwareLimitError
        Naming the worst event, where it came from, and the peak as a percentage of the limit.
    """
    configured = getattr(opts, 'max_b1', None)
    limit = float('nan') if configured is None else float(configured)
    # **Positive** infinity, and only that, means "design without a transmit limit".  NaN, -inf
    # and an absent limit are not spellings of it: every comparison against NaN is False, so
    # treating it as unlimited would let any pulse through silently, which is the failure this
    # check exists to prevent.
    if np.isposinf(limit):
        return
    worst = 0.0
    where = ''
    saw_rf = False
    for index in sorted(seq.block_events):
        block = seq.get_block(index)
        rf = getattr(block, 'rf', None)
        if rf is None:
            continue
        saw_rf = True
        peak = float(np.abs(np.asarray(rf.signal)).max())
        if peak > worst:
            # Block ids are 1-based; `origins` is a list in emission order, so it is not.
            path = origins[index - 1] if 0 < index <= len(origins) else ()
            worst, where = peak, f'block {index} ({".".join(path) or "?"})'
    if not saw_rf:
        return
    usable = limit > 0.0                              # False for 0, negatives and NaN alike
    if usable and worst <= limit:
        return
    stated = 'max_b1 is not set' if configured is None else f'max_b1 is {limit:g}'
    headline = (f'an RF event peaks at {worst / limit * 100:.0f} % of max_b1.' if usable
                else f'{stated}, which is not a transmit limit an RF event can meet.')
    fixes = ([
        'lengthen the pulse or lower its flip angle -- peak B1 scales with both for a fixed shape',
        'an adiabatic pulse instead needs a narrower frequency sweep',
        'or raise max_b1, if the scanner really does deliver it',
    ] if usable else [
        'set a real max_b1 -- sc.convert(20, "uT", "Hz") is a common value',
        'or use max_b1=inf to design without a transmit limit',
    ])
    raise HardwareLimitError(format_error(
        headline, {'from': where, 'peak_b1_hz': round(worst, 2), 'max_b1_hz': limit}, fixes))


def check_label_addresses(seq: Any) -> None:
    """
    Check that no two imaging ADCs write the same k-space address.

    The highest-value check available on a finished sequence: a duplicate address means two
    readouts landing in the same place, which catches a wrong slice order, an off-by-one
    partial-Fourier start and a mis-nested loop from one assertion.

    Raises
    ------
    CompileError
        Naming how many readouts repeat and the first address that does.
    """
    try:
        labels = seq.evaluate_labels(evolution='adc')
    except (AttributeError, ValueError, IndexError):  # pragma: no cover - older pypulseq
        return
    keys = [k for k in ADDRESS_KEYS if k in labels]
    if not keys:
        return
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
        return
    seen: set[tuple[int, ...]] = set()
    first: tuple[int, ...] = ()
    for address in addresses:
        if address in seen:
            first = address
            break
        seen.add(address)
    msg = format_error(
        f'{duplicates} imaging ADC(s) repeat a k-space address -- two readouts are writing the '
        f'same location.',
        {'first repeat': dict(zip(keys, first)), 'readouts': len(addresses)},
        [
            'check the slice order, the phase-encode table, and that every loop index reaches '
            'a label',
            'a readout that is deliberately not imaged should carry a NOISE, REF or NAV label, '
            'which excludes it from this check',
        ],
    )
    raise CompileError(msg)


# --------------------------------------------------------------- did the compile mean the tree
def expected_addresses(
    placed: Sequence[PlacedEvent],
    targets: Mapping[int, float],
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
            if p.kind in LABEL_KINDS
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


def _address_violations(
    placed: Sequence[PlacedEvent],
    targets: Mapping[int, float] | None,
    label_states: Callable[[], Mapping[str, Any]],
) -> list[ContractViolation]:
    """
    Check every ADC's compiled label state against the fold of the tree's labels.

    The check the duplicate-address test cannot do.  A labelling shifted by one readout stays
    unique, so a uniqueness test sees nothing wrong -- and that is precisely the shape of the bug
    retiming exists to prevent, so it needs an invariant of its own rather than trust.
    """
    if targets is None:
        return []
    n_labels = sum(1 for p in placed if p.kind in LABEL_KINDS)
    if not n_labels:
        return []
    # An orphan label's placement is reported as a warning rather than defined, so its effect
    # is not predictable from a time fold and this check would be comparing against a guess.
    if any(i for i, p in enumerate(placed) if p.kind in LABEL_KINDS and i not in targets):
        return []

    expected = expected_addresses(placed, targets)
    try:
        got = label_states()
    except (AttributeError, ValueError, IndexError):  # pragma: no cover - older pypulseq
        return []
    if not got or not expected:
        return []

    out: list[ContractViolation] = []
    for key in sorted({k for state in expected for k in state}):
        series = np.atleast_1d(np.asarray(got.get(key, 0)))
        series = np.resize(series, len(expected)) if series.size else np.zeros(len(expected))
        for index, state in enumerate(expected):
            if int(series[index]) != int(state.get(key, 0)):
                out.append(ContractViolation(
                    f'adc {index}',
                    f'label {key} is {int(series[index])} on readout {index} but the tree '
                    f'implies {int(state.get(key, 0))} -- a label reached the wrong readout',
                ))
                break       # one report per key; a shift corrupts every later readout too
    return out


def _traversed_area(times: np.ndarray, amps: np.ndarray) -> float:
    """
    Return the integral of ``|g|`` over piecewise-linear knots, in 1/m.

    Exact: a segment whose ends have opposite signs is split at its zero crossing, where a
    trapezoidal rule on ``|amps|`` would overstate it.  Additive over a partition at knots, so
    dividing one waveform into consecutive events leaves the total unchanged.
    """
    a, b, dt = amps[:-1], amps[1:], np.diff(times)
    out = 0.5 * np.abs(a + b) * dt
    crossing = a * b < 0
    if np.any(crossing):
        ac, bc, dc = a[crossing], b[crossing], dt[crossing]
        out[crossing] = dc * (ac * ac + bc * bc) / (2.0 * (np.abs(ac) + np.abs(bc)))
    return float(np.sum(out))


# ------------------------------------------------- m0 and m1 about each block's own start (trial)
#: PyPulseq's compress_shape rounds each peak-normalised sample to the nearest 1e-7.
_QUANTUM = 1e-7
#: PyPulseq keys every library entry by its values printed with '%.6g'; two events sharing a key
#: share the first one's stored values, which differ by less than 1e-5 relative.
_KEY_SHARING = 1e-5
#: Float, on both sides, and in PyPulseq's cumulative-sum decompression.
_FLOAT = 1e-12
#: axis_gradient sets a superposed block edge below this fraction of max_grad to exactly 0.
_SNAP = 1e-6
_ROUNDOFF = 2.0 ** -53
_TPS = 1e12


def _ticks(times: Any) -> np.ndarray:
    return np.rint(np.asarray(times, dtype=float) * _TPS).astype(np.int64)


def ready_gradients(blocks: Sequence[PulseqReadyBlock]) -> tuple[np.ndarray, list]:
    """Block edges in ticks, and each ready block's gradients as absolute-tick knots."""
    edges = [int(round(blocks[0].start * _TPS))] if blocks else [0]
    pieces = []
    for block in blocks:
        start = int(round(block.start * _TPS))
        edges.append(int(round(block.end * _TPS)))
        for event in block.events:
            if getattr(event, 'type', None) in GRADIENT_KINDS:
                rel, amps = knots_of(event, 0.0)
                pieces.append((event.channel, start + _ticks(rel), np.asarray(amps, float)))
    return np.asarray(edges, np.int64), pieces


def _registered_gradients(seq: Any) -> tuple[np.ndarray, list]:
    """The same for the built sequence: edges from its cumulative durations, summed in ticks."""
    ids = sorted(seq.block_events)
    durations = [int(round(float(seq.block_durations[i]) * _TPS)) for i in ids]
    edges = np.concatenate([[0], np.cumsum(durations)]).astype(np.int64)
    pieces = []
    for k, i in enumerate(ids):
        block = seq.get_block(i)
        for axis in AXES:
            grad = getattr(block, f'g{axis}', None)
            if grad is not None:
                rel, amps = knots_of(grad, 0.0)
                pieces.append((axis, edges[k] + _ticks(rel), np.asarray(amps, float)))
    return edges, pieces


def _tree_gradients(placed: Sequence[PlacedEvent]) -> list:
    out = []
    for p in placed:
        if p.kind in GRADIENT_KINDS:
            rel, amps = knots_of(p.event, 0.0)
            out.append((p.event.channel, int(round(p.node_t * _TPS)) + _ticks(rel),
                        np.asarray(amps, float)))
    return out


def _segments(t: np.ndarray, g: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, ...]:
    """Split one piece at the partition edges strictly inside it, exactly (it is linear)."""
    i0 = int(np.searchsorted(edges, t[0], side='right'))
    i1 = int(np.searchsorted(edges, t[-1], side='left'))
    if i1 > i0:
        tt = np.union1d(t, edges[i0:i1])
        gg = np.interp(tt.astype(float), t.astype(float), g)
    else:
        tt, gg = t, g
    return np.searchsorted(edges, tt[:-1], side='right') - 1, tt[:-1], tt[1:], gg[:-1], gg[1:]


def _interval_moments(pieces: list, edges: np.ndarray) -> tuple[dict, dict]:
    """
    Per axis and interval [a, b): m0, m1 about a, and the tolerance ingredients.

    T0 and T1 integrate |g| (weighted by t - a for T1), S0 and S1 the support, P the summed peak
    and slope the summed peak slope.  Content outside the partition is returned separately.
    """
    n = edges.size - 1
    out: dict = {}
    outside: dict = {}
    for axis, t, g in pieces:
        if t.size < 2:
            continue
        idx, t0, t1, g0, g1 = _segments(t, g, edges)
        inside = (idx >= 0) & (idx < n)
        if not inside.all():
            outside[axis] = outside.get(axis, 0.0) + float(np.sum(
                0.5 * (t1[~inside] - t0[~inside]) / _TPS * np.abs(g0[~inside] + g1[~inside])))
        idx, t0, t1, g0, g1 = idx[inside], t0[inside], t1[inside], g0[inside], g1[inside]
        if idx.size == 0:
            continue
        a = edges[idx]
        h = (t1 - t0) / _TPS
        u0, u1 = (t0 - a) / _TPS, (t1 - a) / _TPS
        b0, b1 = np.abs(g0), np.abs(g1)
        support = ((g0 != 0) | (g1 != 0)).astype(float)
        values = {
            'm0': 0.5 * h * (g0 + g1),
            'm1': h * (g0 * (2 * u0 + u1) + g1 * (u0 + 2 * u1)) / 6.0,
            'T0': 0.5 * h * (b0 + b1),
            'T1': h * (b0 * (2 * u0 + u1) + b1 * (u0 + 2 * u1)) / 6.0,
            'S0': h * support,
            'S1': 0.5 * h * (u0 + u1) * support,
        }
        starts = np.flatnonzero(np.r_[True, idx[1:] != idx[:-1]])
        keys = idx[starts]
        acc = out.setdefault(axis, {f: np.zeros(n) for f in (*values, 'P', 'slope')})
        for field, v in values.items():
            acc[field][keys] += np.add.reduceat(v, starts)
        acc['P'][keys] += np.maximum.reduceat(np.maximum(b0, b1), starts)
        acc['slope'][keys] += np.maximum.reduceat(np.abs(g1 - g0) / h, starts)
    return out, outside


def _merge_resampled(edges: np.ndarray, resampled: frozenset[int]) -> np.ndarray:
    """Drop the edge between two consecutive resampled blocks: a resample keeps the run's moments."""
    keep = [0] + [k for k in range(1, edges.size - 1)
                  if not (k - 1 in resampled and k in resampled)] + [edges.size - 1]
    return edges[keep]


def _snap_bounds(pieces: list, blocks: np.ndarray, edges: np.ndarray, max_grad: float) -> dict:
    """
    What the endpoint snap can change, from the tree: at each inner block edge where the tree's
    summed value v is under the threshold but not zero, a triangle of height |v| reaching the
    nearest knot on each side.
    """
    n = edges.size - 1
    length = np.diff(edges) / _TPS
    out: dict = {}
    for e in blocks[1:-1]:
        for left in (True, False):
            value: dict = {}
            reach: dict = {}
            for axis, t, g in pieces:
                if not (t[0] < e <= t[-1] if left else t[0] <= e < t[-1]):
                    continue
                value[axis] = value.get(axis, 0.0) + float(np.interp(float(e), t.astype(float), g))
                h = int(e - t[t < e][-1]) if left else int(t[t > e][0] - e)
                reach[axis] = min(reach.get(axis, h), h)
            for axis, v in value.items():
                if not 0.0 < abs(v) < _SNAP * max_grad:
                    continue
                lo, hi = (e - reach[axis], e) if left else (e, e + reach[axis])
                i0 = max(int(np.searchsorted(edges, lo, side='right')) - 1, 0)
                i1 = min(int(np.searchsorted(edges, hi, side='left')), n)
                acc = out.setdefault(axis, (np.zeros(n), np.zeros(n)))
                for i in range(i0, i1):
                    overlap = (min(hi, edges[i + 1]) - max(lo, edges[i])) / _TPS
                    if overlap > 0:
                        acc[0][i] += abs(v) * overlap
                        acc[1][i] += abs(v) * overlap * length[i]
    return out


def _local_violations(
    placed: Sequence[PlacedEvent],
    seq: Any,
    ready: tuple[np.ndarray, list],
    resampled: frozenset[int],
    max_grad: float,
) -> list[ContractViolation]:
    """
    m0 and m1 per compiled block, each about the block's own start, in two stages.

    **Stage A, the compile:** the tree's events restricted exactly to each block, against the
    ready blocks SeqCraft handed to PyPulseq.  Only float, the endpoint snap and the float of
    interpolating at absolute times separate them.  **Stage B, registration:** the ready blocks
    against what PyPulseq stored, within its storage precision -- compression to 1e-7 of the peak
    and the six-digit key it deduplicates by.  Neither stage sees an absolute time, so a correct
    waveform compares the same wherever it is placed.
    """
    blocks, registered = _registered_gradients(seq)
    ready_edges, ready_pieces = ready
    out: list[ContractViolation] = []
    if not np.array_equal(blocks, ready_edges):
        return [ContractViolation('block grid', 'the sequence holds different block edges from '
                                  'the ones emission was handed')]
    edges = _merge_resampled(blocks, resampled)
    length = np.diff(edges) / _TPS
    tree, tree_out = _interval_moments(_tree_gradients(placed), edges)
    mine, mine_out = _interval_moments(ready_pieces, edges)
    stored, stored_out = _interval_moments(registered, edges)
    snap = _snap_bounds(_tree_gradients(placed), blocks, edges, max_grad)

    for axis in sorted(set(tree_out) | set(mine_out) | set(stored_out)):
        out.append(ContractViolation(f'axis {axis}', 'gradient lies outside the sequence'))

    zero = np.zeros(edges.size - 1)
    for axis in sorted(set(tree) | set(mine) | set(stored)):
        t, r, s = tree.get(axis), mine.get(axis), stored.get(axis)
        p = t['P'] if t else zero
        s0, s1 = snap.get(axis, (zero, zero))
        placement = 2.0 * _ROUNDOFF * (edges[1:] / _TPS) * (t['slope'] if t else zero) * length
        tolerance_a = {
            'm0': (_FLOAT * (t['T0'] + p * length) + s0 + placement) if t else zero,
            'm1': (_FLOAT * (t['T0'] * length + p * length ** 2 / 2) + s1 + placement * length)
            if t else zero,
        }
        rp = r['P'] if r else zero
        tolerance_b = {
            'm0': ((_QUANTUM / 2 + _KEY_SHARING) * rp * r['S0']
                   + _FLOAT * (r['T0'] + rp * length)) if r else zero,
            'm1': ((_QUANTUM / 2 + _KEY_SHARING) * rp * r['S1']
                   + _FLOAT * (r['T0'] * length + rp * length ** 2 / 2)) if r else zero,
        }
        for stage, got, want, tolerances in (('compile', r, t, tolerance_a),
                                             ('registration', s, r, tolerance_b)):
            for moment in ('m0', 'm1'):
                g = got[moment] if got else zero
                w = want[moment] if want else zero
                difference = g - w
                bad = np.flatnonzero(np.abs(difference) > tolerances[moment])
                if not bad.size:
                    continue
                k = int(bad[np.argmax(np.abs(difference[bad]) / np.maximum(tolerances[moment][bad], 1e-300))])
                unit = '1/m' if moment == 'm0' else 's/m'
                what = ('the tree asked for' if stage == 'compile' else 'emission handed PyPulseq')
                out.append(ContractViolation(
                    f'axis {axis} {moment}, block interval at {edges[k] / _TPS * 1e3:.3f} ms',
                    f'compiled {moment} {g[k]:.6g} {unit} about the block start differs from what '
                    f'{what}, {w[k]:.6g} {unit}, by {difference[k]:.3g}, beyond the {stage} '
                    f'tolerance {tolerances[moment][k]:.3g}'
                    + (f' ({bad.size} intervals)' if bad.size > 1 else ''),
                ))
    return out


def verify_against_tree(
    placed: Sequence[PlacedEvent],
    targets: Mapping[int, float] | None,
    *,
    duration_s: float,
    tree_duration_s: float,
    moments: Callable[[int], dict[str, float]],
    label_states: Callable[[], Mapping[str, Any]],
    sequence: Any = None,
    ready: tuple[np.ndarray, list] | None = None,
    resampled: frozenset[int] = frozenset(),
    max_grad: float = 0.0,
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

    **Label addresses** must match the fold of the tree's labels.  A duplicate-address check
    only fires when two addresses *collide*; an addressing shifted by one but still unique
    passes it, which is exactly how mis-retimed labels used to escape.  Skipped when a label has
    no ADC after it, since its placement is then reported as a warning rather than defined.

    Parameters
    ----------
    placed, targets
        The compiler's own IR: the placed events, and the label-assignment times chosen for them.
    duration_s, tree_duration_s
        What the compiled sequence measures, and what the tree asked for.
    moments
        ``order -> {axis: moment}`` over the *compiled* blocks.  A callable rather than two
        dicts, so a tree with no gradients pays for neither.
    label_states
        ``Sequence.evaluate_labels(evolution='adc')``, deferred for the same reason.

    Raises
    ------
    CompilerContractError
        If any invariant fails.  Every one is measured first, so the message carries all of them:
        one mis-timed stage usually breaks several axes, and the pattern is the diagnosis.
    """
    violations: list[ContractViolation] = []
    # Tolerance scales with the total, because it has to: float64 resolves about 4 ns at
    # 20 000 s, so demanding nanosecond agreement on a long acquisition would flag every
    # sequence over an hour -- and print two identical-looking numbers while doing it.
    tolerance = max(EPS, 1e-12 * tree_duration_s)
    if abs(duration_s - tree_duration_s) > tolerance:
        violations.append(
            ContractViolation(
                'sequence duration',
                f'compiled duration {duration_s * 1e6:.1f} us differs from the tree total '
                f'{tree_duration_s * 1e6:.1f} us',
            )
        )

    # Trial: with the ready blocks in hand, m0 and m1 are compared per block, about each block's
    # start, and the global sums below are not checked at all -- they are implied.
    if sequence is not None and ready is not None:
        violations.extend(_local_violations(placed, sequence, ready, resampled, max_grad))
        violations.extend(_address_violations(placed, targets, label_states))
        require_valid_contract('against-tree', violations)
        return

    want: dict[str, float] = {}
    # m0's tolerance is scaled by the area each event *traverses*, the integral of |g|, not by its
    # net.  A readout and its prephaser very nearly cancel, so a relative tolerance on the net
    # would demand exactness that float summation over thousands of pieces cannot deliver -- while
    # a tolerance on the total still catches a whole lost lobe, which is what this is for.
    #
    # Traversed area is additive over a division of one waveform in time, so cutting it into
    # consecutive events leaves the scale unchanged; overlapping pieces that cancel can only raise
    # it, because the sum of their integrals of |g| is at least the integral of |their sum|.  It is
    # therefore never below the waveform's own traversed area -- which is not the same as being
    # independent of how events are packaged.  Summing each event's |net| collapsed to the 1.0
    # floor for any one event that changes sign and nets zero, such as a closed trajectory's axis.
    # At that floor PyPulseq's registration was refused: it rounds a compressible shape's
    # peak-normalised samples to 1e-7, moving the area by at most 5e-8 * peak * duration, and in
    # practice by a few percent of that.  Under that model 1e-6 of the traversed area is
    # sufficient to cover the worst case whenever a shape's mean |g| is at least 5 % of its peak;
    # it is not necessary, and lower-duty shapes are not thereby refused.  Writing the file
    # rounds again, separately, and this check never sees that.
    traversed: dict[str, float] = {}
    # Each event's |net area|, summed: m1's relative term, unchanged -- see the m1 check below.
    magnitude: dict[str, float] = {}
    # Peak |g| per axis.  m1's floor tracks this and nothing else -- see the m1 check below.
    peak: dict[str, float] = {}
    # m1 referenced to the start of the sequence, computed **exactly** on both sides.
    #
    # Both are sums of per-piece moments, which is legitimate because a moment is linear in
    # the waveform: the moment of a sum is the sum of the moments, whether the pieces overlap
    # (the tree, where they superpose) or abut (the compiled blocks, where they concatenate).
    # No union of knots needs building.
    want_m1: dict[str, float] = {}
    for p in placed:
        if p.kind in GRADIENT_KINDS:
            channel = p.event.channel
            knots = knots_of(p.event, p.node_t)
            area = pwl_moment(*knots, 0)
            want[channel] = want.get(channel, 0.0) + area
            traversed[channel] = traversed.get(channel, 0.0) + _traversed_area(*knots)
            magnitude[channel] = magnitude.get(channel, 0.0) + abs(area)
            peak[channel] = max(peak.get(channel, 0.0), float(np.max(np.abs(knots[1]))))
            want_m1[channel] = want_m1.get(channel, 0.0) + pwl_moment(*knots, 1)

    if want:
        actual = moments(0)
        actual_m1 = moments(1)
        for axis, expected in want.items():
            difference = actual.get(axis, 0.0) - expected
            m0_tolerance = 1e-6 * max(traversed[axis], 1.0)
            if abs(difference) > m0_tolerance:
                violations.append(
                    ContractViolation(
                        f'axis {axis} m0',
                        f'compiled m0 {actual.get(axis, 0.0):.6g} 1/m differs from the tree '
                        f'sum {expected:.6g} 1/m by {difference:.3g}, beyond the tolerance '
                        f'{m0_tolerance:.3g} (1e-6 of the {traversed[axis]:.6g} 1/m the axis '
                        f'traverses) -- the emitted blocks do not carry the area the tree asked for',
                    )
                )
            # Both sides are exact closed forms, so the only error is float summation over the
            # pieces.  Scaled by each event's |net area|, summed, times the horizon: a readout
            # and its prephaser nearly cancel, so a tolerance on the total net would demand more
            # than float64 can carry.  Deliberately left on |net| when m0 moved to traversed
            # area -- m1 has its own floor for representation noise below, and a wider relative
            # term would only loosen the time-shift check this exists for.  A lobe displaced by
            # one raster changes m1 by area * 10 us, comfortably above this.
            #
            # **Not covered: registration error far from the sequence start.**  m1 is referenced
            # to t = 0, so a shape whose stored area moved by delta_m0 at time t moves m1 by
            # t * delta_m0 as well, while this tolerance does not grow with t.  A compressible shape
            # that compiles at the start of a sequence can therefore be refused here when the same
            # shape is placed a second or more later -- on main before the m0 change as after it.
            # No tolerance alone separates that from a mistimed lobe late in a long sequence;
            # comparing m1 per block would.
            #
            # **Plus a floor proportional to peak |g|**, because the relative term alone sits
            # inside the noise for a whole class of legal input.  An arbitrary gradient that the
            # compiler has to split and rebuild does not round-trip exactly, and the residual is
            # ``~8.4e-14 * peak`` -- measured constant across knot count, amplitude and
            # ``max_grad``, in ``tests/compiler/test_m1_floor.py``.  The relative term scales as
            # traversed x horizon, which grows faster than that, so the two cross: short or weak
            # waveforms were refused while long strong ones passed, *for the same physics*.  Six
            # of ten legal sign-changing gradients were refused before this floor existed, and
            # which six was not predictable from the sequence.
            #
            # 1e-11 clears the measured floor by ~120 and stays three to four orders of magnitude
            # below the signal: one raster of displacement moves m1 by ``traversed * raster``,
            # which is 1e-2 .. 1e0 s/m on the same waveforms.  Both bounds are properties of the
            # representation and of the failure being caught, not of any one module.
            scale = max(magnitude[axis], 1.0)
            m1_scale = max(scale * max(tree_duration_s, 1e-3), 1.0)
            m1_tolerance = max(1e-9 * m1_scale, 1e-11 * peak.get(axis, 0.0))
            if abs(actual_m1.get(axis, 0.0) - want_m1.get(axis, 0.0)) > m1_tolerance:
                violations.append(
                    ContractViolation(
                        f'axis {axis} m1',
                        f'compiled m1 {actual_m1.get(axis, 0.0):.6g} s/m differs from the '
                        f'tree sum {want_m1.get(axis, 0.0):.6g} s/m -- a gradient plays at '
                        f'the wrong time, which m0 cannot see',
                    )
                )

    violations.extend(_address_violations(placed, targets, label_states))
    require_valid_contract('against-tree', violations)


__all__ = [
    'CompilerContractError',
    'ContractViolation',
    'check_event_sizes',
    'check_label_addresses',
    'check_rf_amplitude',
    'expected_addresses',
    'require_valid_contract',
    'verify_against_tree',
    'verify_placed_events',
    'verify_ready_blocks',
]
