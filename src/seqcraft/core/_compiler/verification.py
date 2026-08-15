"""Structural verification for compiler intermediate representations.

These checks protect stage boundaries rather than scanner semantics.  The existing compiler
invariants remain authoritative for waveform moments, labels, limits, and emitted duration.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from .model import (
    EXCLUSIVE_KINDS,
    HANDLED_KINDS,
    PlacedEvent,
    PulseqReadyBlock,
    interval_duration,
    time_at_or_before,
    time_equal,
)


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


__all__ = [
    'CompilerContractError',
    'ContractViolation',
    'require_valid_contract',
    'verify_placed_events',
    'verify_ready_blocks',
]
