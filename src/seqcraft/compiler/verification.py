"""Verification: the compiler checking its own work.

Two layers, and the difference is worth stating.

**Contracts** (:func:`verify_placed_events`, :func:`verify_ready_blocks`) protect *stage
boundaries*.  They ask whether one stage handed the next a well-formed intermediate
representation -- finite times, intervals the right way round, blocks that abut, at most one RF
and one ADC.  A failure is a compiler bug by definition, so it raises
:class:`CompilerContractError` rather than becoming a report.

**Semantics** (:func:`verify_against_tree`) asks the harder question: did the compile preserve
what the *tree* meant?  Duration, the zeroth and first gradient moments of every compiled block
about its own start, and the label address each ADC ends up with, all compared against the placed
events the compile started from.  The moments are compared in two stages: the tree against the
ready blocks emission hands PyPulseq, and those against what PyPulseq stored.  A failure here is a
compiler bug too -- the tree was legal, since every legality check has already passed -- so it
raises :class:`CompilerContractError` as well.  All four invariants are measured before anything is
raised, because seeing every affected axis at once beats stopping at the first.

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
from ..design.timing import EPS, TICKS_PER_SECOND, to_ticks
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

    At absolute times, so its float grows with placement.  :func:`verify_against_tree` compares
    block by block instead (:func:`stored_gradients`); this remains for the compiler baseline, which
    digests these numbers.  Deliberately *not* shared with :func:`seqcraft.analysis.moments`, which
    walks the tree: these are different iterations over different structures.

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


# ----------------------------------------------- m0 and m1 per compiled block, about its own start
#: The unit roundoff of float64.
_U = 2.0**-53
_TPS = float(TICKS_PER_SECOND)
#: At most this many roundings enter one segment's term before it is summed.  A sum of N such terms,
#: in any order, is then within (N + _TERM) * u * sum|term| of the exact one.
_TERM = 10
#: PyPulseq's ``compress_shape`` stores the differences of a peak-normalised shape, each rounded to
#: this quantum.
_QUANTUM = 1e-7
#: What that leaves in a stored sample, as a fraction of the peak.  Half a quantum from rounding --
#: plus one more quantum in every sample after the first when the first lies exactly on a half
#: quantum, because rounding half away from zero then carries an offset into the cumulative sum.
#: Ordinary shapes see the half quantum; the bound has to be the edge case.
_COMPRESSION = 1.5 * _QUANTUM
#: PyPulseq keys every library entry by its values printed with ``'%.6g'``, and stores an event that
#: prints like an earlier one as that one.  Two numbers with one such key differ by less than a unit
#: in their sixth significant digit: under 1e-5 of the larger.
_KEY = 1e-5
#: A normalised shape difference this large prints more coarsely than one quantum, so two shapes can
#: share a key while differing there.  Smaller differences print every quantum.
_KEY_COARSE = 0.1
#: Absolute times reach the compiler as float seconds.  The raster check puts a gradient's start
#: within half a tick of its tick, and turning a knot's time into ticks takes a few more roundings
#: of about u * t each: adding the knot's offset to the start, converting, and a block edge's own
#: round trip through seconds.  So a knot can land a tick off on either side, and further once
#: u * t approaches a tick, past about 1100 s.
_TIME_OPERATIONS = 4


@dataclass(frozen=True)
class BlockMoments:
    """
    One side's gradients per compiled block, as m0 and m1 about each block's own start.

    Rows are blocks and columns the axes ``x, y, z``.  ``edges`` are the block boundaries in ticks,
    the cumulative block durations from zero, so no lever arm involves an absolute float time.
    ``error0`` and ``error1`` bound the float error of ``m0`` and ``m1`` as computed here.
    ``overrun`` names each ``(block, axis)`` whose gradient reaches outside its own block.
    """

    edges: np.ndarray
    m0: np.ndarray
    m1: np.ndarray
    error0: np.ndarray
    error1: np.ndarray
    overrun: tuple[tuple[int, str], ...]


@dataclass(frozen=True)
class HandedGradients(BlockMoments):
    """
    The ready blocks as emission hands them to PyPulseq, measured before registration.

    ``storage0`` and ``storage1`` bound how far PyPulseq's storage precision may move ``m0`` and
    ``m1``: compression rounding, and the six-digit key it deduplicates events by.
    """

    storage0: np.ndarray
    storage1: np.ndarray


def _knots(event: Any) -> tuple[np.ndarray, np.ndarray]:
    """Return a gradient event's knots in ticks from its own reference, and their amplitudes."""
    times, amps = knots_of(event, 0.0)
    return np.rint(np.asarray(times) * _TPS).astype(np.int64), np.asarray(amps, dtype=float)


def _segments(
    pieces: Sequence[tuple[np.ndarray, np.ndarray]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return every piece's consecutive knots as segments: ``piece, t0, t1, g0, g1``."""
    if not pieces:
        empty = np.zeros(0, dtype=np.int64)
        return empty, empty, empty, np.zeros(0), np.zeros(0)
    sizes = np.array([times.size for times, _ in pieces])
    times = np.concatenate([times for times, _ in pieces])
    amps = np.concatenate([amps for _, amps in pieces])
    opens = np.ones(times.size, dtype=bool)
    opens[np.cumsum(sizes) - 1] = False            # a piece's last knot opens no segment
    i = np.flatnonzero(opens)
    piece = np.repeat(np.arange(len(pieces)), sizes)[i]
    return piece, times[i], times[i + 1], amps[i], amps[i + 1]


def _sums(
    axis: np.ndarray,
    block: np.ndarray,
    t0: np.ndarray,
    t1: np.ndarray,
    g0: np.ndarray,
    g1: np.ndarray,
    edges: np.ndarray,
) -> dict[str, np.ndarray]:
    """
    Per block and axis: m0 and m1 about the block start, the ``|g|`` integrals that bound their
    float error (``T0``, and ``T1`` weighted by ``t - a``), and the number of segments summed.

    Exact per linear segment; every lever arm is a difference of ticks.
    """
    n = edges.size - 1
    a = edges[block]
    h = (t1 - t0) / _TPS
    u0 = (t0 - a) / _TPS
    u1 = (t1 - a) / _TPS
    b0, b1 = np.abs(g0), np.abs(g1)
    key = block * 3 + axis

    def total(weights: np.ndarray | None = None) -> np.ndarray:
        return np.bincount(key, weights=weights, minlength=3 * n).reshape(n, 3).astype(float)

    return {
        'm0': total(0.5 * h * (g0 + g1)),
        'm1': total(h * (g0 * (2 * u0 + u1) + g1 * (u0 + 2 * u1)) / 6.0),
        'T0': total(0.5 * h * (b0 + b1)),
        'T1': total(h * (b0 * (2 * u0 + u1) + b1 * (u0 + 2 * u1)) / 6.0),
        'N': total(),
    }


def _float_error(sums: Mapping[str, np.ndarray], order: int) -> np.ndarray:
    """Bound the float error of one side's m0 (order 0) or m1 (order 1), as :func:`_sums` sums it."""
    return np.asarray((sums['N'] + _TERM) * _U * sums[f'T{order}'], dtype=float)


def _block_side(
    pieces: Sequence[tuple[int, int, np.ndarray, np.ndarray]],
    edges: np.ndarray,
) -> tuple[dict[str, np.ndarray], tuple[tuple[int, str], ...]]:
    """Sum pieces that each belong to one block, ``(block, axis, ticks, amps)``, per block and axis."""
    piece, t0, t1, g0, g1 = _segments([(times, amps) for _, _, times, amps in pieces])
    block = np.array([k for k, _, _, _ in pieces], dtype=np.int64)[piece]
    axis = np.array([a for _, a, _, _ in pieces], dtype=np.int64)[piece]
    late = (t0 < edges[block]) | (t1 > edges[block + 1])
    overrun = tuple(sorted({(int(k), AXES[int(a)]) for k, a in zip(block[late], axis[late])}))
    return _sums(axis, block, t0, t1, g0, g1, edges), overrun


def handed_gradients(blocks: Sequence[PulseqReadyBlock]) -> HandedGradients:
    """
    Measure the ready blocks' gradients as emission hands them to PyPulseq.

    Run before emission, while legalization's tuple still holds every block, so nothing PyPulseq
    does when it registers them can change what is measured.  Only numbers are kept: about 170 bytes
    per block, rather than the events.
    """
    durations = [to_ticks(block.duration) for block in blocks]
    edges = np.concatenate([[0], np.cumsum(durations, dtype=np.int64)]).astype(np.int64)
    pieces: list[tuple[int, int, np.ndarray, np.ndarray]] = []
    shapes: list[tuple[int, int, float, float, float, float]] = []
    for k, block in enumerate(blocks):
        for event in block.events:
            if getattr(event, 'type', None) not in GRADIENT_KINDS:
                continue
            ticks, amps = _knots(event)
            if ticks.size < 2:
                continue
            axis = AXES.index(event.channel)
            pieces.append((k, axis, edges[k] + ticks, amps))
            shapes.append((k, axis, *_storage_terms(event, ticks, amps)))
    sums, overrun = _block_side(pieces, edges)

    # PyPulseq's storage precision.  Compression moves every stored sample of a shape by at most
    # _COMPRESSION of its peak, anywhere in the shape's span; PyPulseq's own float, decompressing by
    # a cumulative sum, adds a roundoff per sample.  A shared six-digit key moves each value by
    # under _KEY of itself -- pointwise, so the |g| integrals bound it -- and each coarse shape
    # difference by under _KEY of that difference, carried into every later sample.
    n = edges.size - 1
    compression = np.zeros((2, n, 3))
    coarse = np.zeros((2, n, 3))
    for k, axis, quantum, steps, length, lever in shapes:
        compression[:, k, axis] += (quantum * length, quantum * lever)
        coarse[:, k, axis] += (steps * length, steps * lever)
    key = _KEY * (1.0 + _KEY)
    return HandedGradients(
        edges=edges,
        m0=sums['m0'],
        m1=sums['m1'],
        error0=_float_error(sums, 0),
        error1=_float_error(sums, 1),
        overrun=overrun,
        storage0=compression[0] + key * (sums['T0'] + coarse[0]),
        storage1=compression[1] + key * (sums['T1'] + coarse[1]),
    )


def _storage_terms(
    event: Any, ticks: np.ndarray, amps: np.ndarray,
) -> tuple[float, float, float, float]:
    """
    One ready gradient's storage-precision ingredients: the compression bound per sample in Hz/m,
    the sum of its coarse shape differences in Hz/m, and its span's length and first moment about
    the block start.

    A trapezoid is stored as its parameters, so only the key can move it.  A shape is stored
    normalised by its peak sample, which is what the coarse differences are measured against.
    """
    length = (ticks[-1] - ticks[0]) / _TPS
    lever = length * 0.5 * (ticks[0] + ticks[-1]) / _TPS
    if getattr(event, 'type', None) != 'grad':
        return 0.0, 0.0, length, lever
    samples = np.asarray(event.waveform, dtype=float)
    peak = float(np.max(np.abs(samples))) if samples.size else 0.0
    steps = np.abs(np.diff(samples, prepend=0.0))
    coarse = float(np.sum(steps[steps >= (_KEY_COARSE - 2.0 * _QUANTUM) * peak])) if peak else 0.0
    quantum = (_COMPRESSION + (samples.size + _TERM) * _U) * float(np.max(np.abs(amps)))
    return quantum, coarse, length, lever


def stored_gradients(seq: Any) -> BlockMoments:
    """
    Measure what PyPulseq stored, block by block, as :func:`handed_gradients` measures the blocks
    it was handed.  The edges are the stored block durations, summed in ticks.
    """
    ids = sorted(seq.block_events)
    durations = [to_ticks(float(seq.block_durations[i])) for i in ids]
    edges = np.concatenate([[0], np.cumsum(durations, dtype=np.int64)]).astype(np.int64)
    pieces: list[tuple[int, int, np.ndarray, np.ndarray]] = []
    for k, index in enumerate(ids):
        block = seq.get_block(index)
        for axis_index, axis in enumerate(AXES):
            grad = getattr(block, f'g{axis}', None)
            if grad is None:
                continue
            ticks, amps = _knots(grad)
            if ticks.size >= 2:
                pieces.append((k, axis_index, edges[k] + ticks, amps))
    sums, overrun = _block_side(pieces, edges)
    return BlockMoments(edges, sums['m0'], sums['m1'], _float_error(sums, 0), _float_error(sums, 1),
                        overrun)


def _split(
    t0: np.ndarray, t1: np.ndarray, g0: np.ndarray, g1: np.ndarray, edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Cut every segment at the block edges strictly inside it.

    Exact apart from the interpolated value at each cut, since a segment is linear; the order of
    the segments is kept.  Returns each part's original segment, its ends and its end values.
    """
    lo = np.searchsorted(edges, t0, side='right')
    cuts = np.maximum(np.searchsorted(edges, t1, side='left') - lo, 0)
    if not cuts.any():
        return np.arange(t0.size), t0, t1, g0, g1
    parts = cuts + 1
    owner = np.repeat(np.arange(t0.size), parts)
    q = np.arange(owner.size) - np.repeat(np.cumsum(parts) - parts, parts)
    first, last = q == 0, q == cuts[owner]
    top = edges.size - 1
    left = np.where(first, t0[owner], edges[np.clip(lo[owner] + q - 1, 0, top)])
    right = np.where(last, t1[owner], edges[np.clip(lo[owner] + q, 0, top)])
    a, ga, gb = t0[owner], g0[owner], g1[owner]
    with np.errstate(divide='ignore', invalid='ignore'):
        slope = (gb - ga) / (t1[owner] - a)
    return (owner, left, right, np.where(first, ga, ga + slope * (left - a)),
            np.where(last, gb, ga + slope * (right - a)))


def _tree_side(
    placed: Sequence[PlacedEvent],
    edges: np.ndarray,
    resampled: Sequence[tuple[int, str]],
    raster: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, float], np.ndarray]:
    """
    What the tree asked for per block and axis, and how far the compile may move it.

    Every placed gradient is its start in ticks plus its own knots, restricted exactly to each block
    it overlaps.  A block-axis legalization reports as resampled is replaced by the tree
    resampled the same way, onto the gradient raster from the block's start.

    Returns m0, m1, and the two tolerances: this side's float error plus what the compile itself
    may change -- the float of superposing the pieces, which also bounds what ``superpose``'s
    zeroing of an end value can move.  Also the area lying outside the sequence, per axis, and
    which block-axes were modelled as resampled.
    """
    n = edges.size - 1
    pieces: list[tuple[int, np.ndarray, np.ndarray]] = []
    for p in placed:
        if p.kind in GRADIENT_KINDS:
            ticks, amps = _knots(p.event)
            if ticks.size >= 2:
                # Anchored at the start the compiler checked onto the gradient raster.
                offset = to_ticks(p.start) - to_ticks(float(getattr(p.event, 'delay', 0.0) or 0.0))
                pieces.append((AXES.index(p.event.channel), offset + ticks, amps))
    piece, t0, t1, g0, g1 = _segments([(times, amps) for _, times, amps in pieces])
    owner, left, right, gl, gr = _split(t0, t1, g0, g1, edges)
    piece = piece[owner]
    # A value at a cut is interpolated from the whole segment it lies on, whose other end may lie
    # outside the block and be larger than anything in it: its float error scales with those ends.
    bracket = np.maximum(np.abs(g0[owner]), np.abs(g1[owner]))
    axis = np.array([a for a, _, _ in pieces], dtype=np.int64)[piece]
    block = np.searchsorted(edges, left, side='right') - 1
    inside = (block >= 0) & (block < n)
    outside: dict[str, float] = {}
    for a in np.unique(axis[~inside]):
        out = ~inside & (axis == a)
        outside[AXES[int(a)]] = float(np.sum(0.5 * (right[out] - left[out]) / _TPS
                                              * np.abs(gl[out] + gr[out])))
    piece, axis, block = piece[inside], axis[inside], block[inside]
    left, right, gl, gr = left[inside], right[inside], gl[inside], gr[inside]
    bracket = bracket[inside]
    sums = _sums(axis, block, left, right, gl, gr, edges)

    # Per block and axis, summed over the pieces there: their peak |g|, the largest end value of a
    # segment they are interpolated from, and their peak slope; and how many there are.  These bound
    # the float of superposing them.
    flat = 3 * n
    key = block * 3 + axis
    length = right - left

    def per_bin(where: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
        return np.bincount(where, weights=weights, minlength=flat).reshape(n, 3).astype(float)

    peak = ends = steep = count = np.zeros((n, 3))
    if key.size:
        start = np.flatnonzero(np.r_[True, (piece[1:] != piece[:-1]) | (block[1:] != block[:-1])])
        with np.errstate(divide='ignore', invalid='ignore'):
            slope = np.where(length > 0, np.abs(gr - gl) / (length / _TPS), 0.0)
        peak = per_bin(key[start], np.maximum.reduceat(np.maximum(np.abs(gl), np.abs(gr)), start))
        ends = per_bin(key[start], np.maximum.reduceat(bracket, start))
        steep = per_bin(key[start], np.maximum.reduceat(slope, start))
        count = per_bin(key[start])
    variation = per_bin(key, np.abs(gr - gl))

    span = (np.diff(edges) / _TPS)[:, None]
    modelled = np.zeros((n, 3), dtype=bool)
    if resampled:
        modelled, model = _resampled_tree(pieces, edges, resampled, raster)
        sums = {name: np.where(modelled, model[name], given) for name, given in sums.items()}

    # superpose zeroes an end value only when it is within `value` of zero, and that changes the
    # waveform by less than this term already allows.
    value = (count + 8.0) * _U * ends                         # per value: k pieces and their sum
    lever = 4.0 * _U * span * steep                           # the resample, in seconds from the start
    noise = value + lever

    # A knot, on either side, or a block edge can land `shift` from where the other side puts it
    # (_TIME_OPERATIONS).  Moving every knot of a piece by that changes its area by at most `shift`
    # times its variation and its two end values, and its m1 by that times the lever arm, plus the
    # area times the shift.
    shift = (np.ceil(0.5 + _TIME_OPERATIONS * _U * edges[1:]) / _TPS)[:, None]
    moved = shift * (variation + 2.0 * peak)

    interpolated = 4.0 * _U * ends                            # this side's own interpolated values
    tolerance0 = _float_error(sums, 0) + (interpolated + noise) * span + moved
    tolerance1 = (_float_error(sums, 1) + (interpolated + noise) * span**2 / 2.0
                  + moved * span + shift * sums['T0'])
    return sums['m0'], sums['m1'], tolerance0, tolerance1, outside, modelled


def _resampled_tree(
    pieces: Sequence[tuple[int, np.ndarray, np.ndarray]],
    edges: np.ndarray,
    resampled: Sequence[tuple[int, str]],
    raster: int,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """
    The tree interpolated onto the gradient raster from each reported block's start, as
    ``_resampled`` does to the superposed waveform, summed per block and axis like :func:`_sums`.

    Pieces count when they overlap the open block, which is what legalization superposes there.
    One sweep per axis over pieces sorted by start, so the cost is linear in pieces and blocks.
    Returns which block-axes it covers, and the sums.
    """
    modelled = np.zeros((edges.size - 1, 3), dtype=bool)
    wanted: dict[int, list[int]] = {}
    for k, name in sorted(set(resampled)):
        wanted.setdefault(AXES.index(name), []).append(k)
    seg_axis, seg_block, seg_t0, seg_t1, seg_g0, seg_g1 = [], [], [], [], [], []
    for axis, blocks in wanted.items():
        order = sorted((i for i, (a, _, _) in enumerate(pieces) if a == axis),
                       key=lambda i: int(pieces[i][1][0]))
        cursor = 0
        active: list[int] = []
        as_float: dict[int, np.ndarray] = {}
        for k in blocks:
            a, b = int(edges[k]), int(edges[k + 1])
            while cursor < len(order) and int(pieces[order[cursor]][1][0]) < b:
                active.append(order[cursor])
                cursor += 1
            active = [i for i in active if int(pieces[i][1][-1]) > a]
            grid = a + raster * np.arange((b - a) // raster + 1, dtype=np.int64)
            x = grid.astype(float)
            value = np.zeros(grid.size)
            for i in active:
                if i not in as_float:
                    as_float[i] = pieces[i][1].astype(float)
                value += np.interp(x, as_float[i], pieces[i][2], left=0.0, right=0.0)
            modelled[k, axis] = True
            seg_axis.append(np.full(grid.size - 1, axis))
            seg_block.append(np.full(grid.size - 1, k))
            seg_t0.append(grid[:-1])
            seg_t1.append(grid[1:])
            seg_g0.append(value[:-1])
            seg_g1.append(value[1:])
    return modelled, _sums(np.concatenate(seg_axis), np.concatenate(seg_block),
                           np.concatenate(seg_t0), np.concatenate(seg_t1),
                           np.concatenate(seg_g0), np.concatenate(seg_g1), edges)


def _moment_violations(
    stage: str,
    got: tuple[np.ndarray, np.ndarray],
    want: tuple[np.ndarray, np.ndarray],
    tolerance: tuple[np.ndarray, np.ndarray],
    edges: np.ndarray,
    measured: str,
    expected: str,
) -> list[ContractViolation]:
    """One violation per axis and moment: the block furthest beyond its tolerance, and a count."""
    out: list[ContractViolation] = []
    for order, unit in ((0, '1/m'), (1, 's/m')):
        difference = got[order] - want[order]
        bad = (np.abs(difference) > tolerance[order]) | ~np.isfinite(difference)
        for axis_index, axis in enumerate(AXES):
            rows = np.flatnonzero(bad[:, axis_index])
            if not rows.size:
                continue
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = np.abs(difference[rows, axis_index]) / tolerance[order][rows, axis_index]
            k = int(rows[int(np.argmax(np.nan_to_num(ratio, nan=np.inf)))])
            more = f'; {rows.size} blocks differ on this axis' if rows.size > 1 else ''
            out.append(ContractViolation(
                f'axis {axis} m{order}, PyPulseq block {k + 1} at {edges[k] / _TPS * 1e3:.3f} ms',
                f'{measured} m{order} {got[order][k, axis_index]:.6g} {unit} about the block start '
                f'differs from {expected}, {want[order][k, axis_index]:.6g} {unit}, by '
                f'{difference[k, axis_index]:.3g}, beyond the {stage} tolerance '
                f'{tolerance[order][k, axis_index]:.3g}{more}',
            ))
    return out


def verify_against_tree(
    placed: Sequence[PlacedEvent],
    targets: Mapping[int, float] | None,
    *,
    duration_s: float,
    tree_duration_s: float,
    handed: HandedGradients,
    stored: Callable[[], BlockMoments],
    resampled: Sequence[tuple[int, str]],
    grad_raster: float,
    label_states: Callable[[], Mapping[str, Any]],
) -> None:
    """
    Assert the compile preserved what the tree meant.

    Four invariants, each catching a class of compiler bug no individual test case would:

    **Total duration** must match the tree's.  Guards a boundary merge that dropped time
    rather than dropping a cut.

    **m0 and m1 per compiled block and axis, about each block's own start.**  m0 catches a split
    that lost a tail or a merge that dropped a piece.  m1 catches a gradient that plays at the
    wrong moment, which m0 cannot see: a piece of area ``A`` displaced by ``dt`` moves m1 by
    ``A * dt``.  Each block is its own reference, so a lever arm is at most a block long and a
    correct waveform compares the same wherever it sits in the sequence.  The sequence's whole m0
    and m1 are sums of these, ``M0 = sum m0_k`` and ``M1 = sum (m1_k + a_k m0_k)``, so they add
    nothing a block-wise check misses.  Compared in two stages, each against its own input:

    - **compile** -- the tree's placed events, restricted exactly to each block, against the ready
      blocks emission hands to PyPulseq (`handed`).  The tolerance is SeqCraft's own arithmetic:
      float, and a knot landing a tick off.  A block-axis legalization reports as resampled is
      compared with the tree resampled onto the raster.
    - **registration** -- the ready blocks against what PyPulseq stored (`stored`).  The tolerance
      is PyPulseq's storage precision: compression rounding, and its six-digit deduplication key.

    So a refusal names which of the two moved the waveform, and stored against tree is bounded by
    the two tolerances together.

    **Label addresses** must match the fold of the tree's labels.  A duplicate-address check
    only fires when two addresses *collide*; an addressing shifted by one but still unique
    passes it, which is exactly how mis-retimed labels used to escape.  Skipped when a label has
    no ADC after it, since its placement is then reported as a warning rather than defined.

    m0 and m1 are two linear functionals per block, **not waveform fidelity**.  They cannot see a
    change that keeps every block's area and centroid -- a distortion orthogonal to ``1`` and
    ``t - a``, a zero-net piece moved inside its block, two equal pieces swapped -- nor anything
    below the tolerance.  RF and ADC timing, limits and seam continuity are checked elsewhere, and
    the written file is never seen.

    Parameters
    ----------
    placed, targets
        The compiler's own IR: the placed events, and the label-assignment times chosen for them.
        `placed` is the expected side, never read from the output.
    duration_s, tree_duration_s
        What the compiled sequence measures, and what the tree asked for.
    handed
        :func:`handed_gradients` of the ready blocks, measured before emission.
    stored
        :func:`stored_gradients` of the built sequence, deferred so a tree with no gradients pays
        nothing.
    resampled
        The ``(block, axis)`` pairs legalization resampled onto the raster.
    grad_raster
        The gradient raster in seconds, which a resample interpolates onto.
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

    if (any(p.kind in GRADIENT_KINDS for p in placed) or bool(np.any(handed.error0))
            or handed.overrun):
        edges = handed.edges
        tree0, tree1, tolerance0, tolerance1, outside, modelled = _tree_side(
            placed, edges, resampled, to_ticks(grad_raster))
        for axis, area in sorted(outside.items()):
            violations.append(ContractViolation(
                f'axis {axis}', f'the tree puts {area:.3g} 1/m of gradient outside the sequence'))
        for k, axis in handed.overrun:
            violations.append(ContractViolation(
                f'ready block {k}', f'its gradient on axis {axis} reaches outside the block'))
        expected = 'what the tree asked for'
        if modelled.any():
            expected += ' (resampled onto the raster where legalization did)'
        violations.extend(_moment_violations(
            'compile', (handed.m0, handed.m1), (tree0, tree1),
            (handed.error0 + tolerance0, handed.error1 + tolerance1), edges,
            'compiled', expected))

        registered = stored()
        if registered.edges.shape != edges.shape or np.any(registered.edges != edges):
            violations.append(ContractViolation(
                'block grid', f'PyPulseq holds {registered.edges.size - 1} blocks lasting '
                f'{registered.edges[-1] / _TPS * 1e6:.1f} us, but emission handed it '
                f'{edges.size - 1} lasting {edges[-1] / _TPS * 1e6:.1f} us'))
        else:
            for k, axis in registered.overrun:
                violations.append(ContractViolation(
                    f'PyPulseq block {k + 1}', f'its gradient on axis {axis} reaches outside the '
                    f'block'))
            violations.extend(_moment_violations(
                'registration', (registered.m0, registered.m1), (handed.m0, handed.m1),
                (handed.storage0 + handed.error0 + registered.error0,
                 handed.storage1 + handed.error1 + registered.error1), edges,
                'registered', 'what emission handed PyPulseq'))

    violations.extend(_address_violations(placed, targets, label_states))
    require_valid_contract('against-tree', violations)


__all__ = [
    'BlockMoments',
    'CompilerContractError',
    'ContractViolation',
    'HandedGradients',
    'check_event_sizes',
    'check_label_addresses',
    'check_rf_amplitude',
    'expected_addresses',
    'handed_gradients',
    'require_valid_contract',
    'stored_gradients',
    'verify_against_tree',
    'verify_placed_events',
    'verify_ready_blocks',
]
