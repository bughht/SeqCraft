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

Why the semantic check lives in the compiler rather than on the result it checks: it takes
``Sequence[PlacedEvent]`` -- the compiler's private IR, which the result type never sees -- and it
runs once, on the object :func:`~seqcraft.compiler.compile_sequence` has just built.  It is a
compile stage that happens to run last, not an accessor.  Keeping it here is what lets the result
types stay free of any compiler import.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from ..design.events import GRADIENT_KINDS, HANDLED_KINDS, LABEL_KINDS, knots_of, pwl_moment
from ..design.timing import EPS
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


def verify_against_tree(
    placed: Sequence[PlacedEvent],
    targets: Mapping[int, float] | None,
    *,
    duration_s: float,
    tree_duration_s: float,
    moments: Callable[[int], dict[str, float]],
    label_states: Callable[[], Mapping[str, Any]],
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

    want: dict[str, float] = {}
    # Tolerance is scaled by the total area *traversed*, not by the net.  A readout and its
    # prephaser very nearly cancel, so a relative tolerance on the net would demand exactness
    # that float summation over thousands of pieces cannot deliver -- while a tolerance on the
    # total still catches a whole lost lobe, which is what this is for.
    magnitude: dict[str, float] = {}
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
            magnitude[channel] = magnitude.get(channel, 0.0) + abs(area)
            want_m1[channel] = want_m1.get(channel, 0.0) + pwl_moment(*knots, 1)

    if want:
        actual = moments(0)
        actual_m1 = moments(1)
        for axis, expected in want.items():
            scale = max(magnitude[axis], 1.0)
            if abs(actual.get(axis, 0.0) - expected) > 1e-6 * scale:
                violations.append(
                    ContractViolation(
                        f'axis {axis} m0',
                        f'compiled m0 {actual.get(axis, 0.0):.6g} 1/m differs from the tree '
                        f'sum {expected:.6g} 1/m -- a split or a merge lost area',
                    )
                )
            # Both sides are exact closed forms, so the only error is float summation over the
            # pieces.  Scaled by area *traversed* times the horizon, for the same reason m0's
            # is scaled by area traversed: a readout and its prephaser nearly cancel, so a
            # tolerance on the net would demand more than float64 can carry.  A lobe displaced
            # by one raster changes m1 by area * 10 us, comfortably above this.
            m1_scale = max(scale * max(tree_duration_s, 1e-3), 1.0)
            if abs(actual_m1.get(axis, 0.0) - want_m1.get(axis, 0.0)) > 1e-9 * m1_scale:
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
    'expected_addresses',
    'require_valid_contract',
    'verify_against_tree',
    'verify_placed_events',
    'verify_ready_blocks',
]
