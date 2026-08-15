"""Immutable compiler contracts shared by placement, legalization, and emission.

The contracts deliberately hold references to PyPulseq events instead of copying them.  SeqCraft
treats those event objects as read-only during compilation; tuples freeze the compiler-owned
structure without pretending that third-party ``SimpleNamespace`` instances are deeply immutable.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from ..logic import BARRIER
from ..timing import EPS, exact_diff

# Event groups encode block-format constraints, not Python class relationships.  Keeping them in
# one module makes every compiler stage use the same classification vocabulary.
GRADIENT_KINDS = frozenset({'trap', 'grad'})
POINT_KINDS = frozenset({'labelset', 'labelinc', 'trigger', 'output'})
INDIVISIBLE_KINDS = frozenset({'rf', 'adc', 'trigger', 'output'})
EXCLUSIVE_KINDS = frozenset({'rf', 'adc'})
LABEL_KINDS = frozenset({'labelset', 'labelinc'})
HANDLED_KINDS = frozenset(
    {
        BARRIER,
        'adc',
        'delay',
        'grad',
        'labelinc',
        'labelset',
        'output',
        'rf',
        'trap',
        'trigger',
    }
)
AXES = ('x', 'y', 'z')
ADDRESS_KEYS = ('SLC', 'LIN', 'PAR', 'AVG', 'REP', 'SEG', 'ECO', 'SET')


def time_equal(left: float, right: float) -> bool:
    """Return whether two times are equal under the compiler tolerance."""
    return abs(left - right) <= EPS


def time_before(left: float, right: float) -> bool:
    """Return whether ``left`` is materially earlier than ``right``."""
    return left < right - EPS


def time_at_or_before(left: float, right: float) -> bool:
    """Return whether ``left`` is earlier than or tolerance-equal to ``right``."""
    return left <= right + EPS


def time_strictly_between(value: float, start: float, end: float) -> bool:
    """Return whether ``value`` lies materially inside the open interval ``(start, end)``."""
    return start + EPS < value < end - EPS


def interval_duration(start: float, end: float) -> float:
    """Return an interval duration using the compiler's integer-tick arithmetic."""
    return exact_diff(end, start)


@dataclass(frozen=True)
class PlacedEvent:
    """One source-tree leaf resolved to absolute active and reservation intervals."""

    node_t: float
    start: float
    end: float
    res_start: float
    res_end: float
    event: Any
    path: tuple[str, ...]

    @property
    def kind(self) -> str:
        """Return the Pulseq event type, or ``'?'`` for a malformed event."""
        return str(getattr(self.event, 'type', '?'))

    @property
    def where(self) -> str:
        """Return the source tag path, falling back to the event kind."""
        return '.'.join(self.path) if self.path else self.kind

    @property
    def source_path(self) -> tuple[str, ...]:
        """Return the immutable source path using the shared provenance vocabulary."""
        return self.path

    @property
    def duration(self) -> float:
        """Return the active duration using exact tick subtraction."""
        return interval_duration(self.start, self.end)

    @property
    def reservation_duration(self) -> float:
        """Return the duration reserved against indivisible or exclusive events."""
        return interval_duration(self.res_start, self.res_end)

    def summary(self) -> str:
        """Return a concise, stable description suitable for diagnostics."""
        return (
            f'{self.kind} {self.start * 1e6:.1f}..{self.end * 1e6:.1f} us '
            f'(reserved {self.res_start * 1e6:.1f}..{self.res_end * 1e6:.1f} us) '
            f'from {self.where}'
        )

    def __repr__(self) -> str:
        """Avoid expanding waveform arrays embedded in the PyPulseq event."""
        return f'PlacedEvent({self.summary()})'


@dataclass(frozen=True)
class PulseqReadyBlock:
    """An immutable block contract ready for mechanical PyPulseq emission."""

    index: int
    start: float
    end: float
    duration: float
    events: tuple[Any, ...]
    source_paths: tuple[tuple[str, ...], ...]
    origin: tuple[str, ...]

    @property
    def kinds(self) -> tuple[str, ...]:
        """Return event types in their emission order."""
        return tuple(str(getattr(event, 'type', '?')) for event in self.events)

    def summary(self) -> str:
        """Return a concise, stable description suitable for diagnostics."""
        kinds = ','.join(self.kinds) or 'delay'
        where = '.'.join(self.origin) or '?'
        return (
            f'block {self.index} {self.start * 1e6:.1f}..{self.end * 1e6:.1f} us ' f'[{kinds}] from {where}'
        )

    def __repr__(self) -> str:
        """Avoid expanding event payloads while keeping timing and provenance visible."""
        return f'PulseqReadyBlock({self.summary()})'


__all__ = [
    'ADDRESS_KEYS',
    'AXES',
    'EXCLUSIVE_KINDS',
    'GRADIENT_KINDS',
    'HANDLED_KINDS',
    'INDIVISIBLE_KINDS',
    'LABEL_KINDS',
    'POINT_KINDS',
    'PlacedEvent',
    'PulseqReadyBlock',
    'interval_duration',
    'time_at_or_before',
    'time_before',
    'time_equal',
    'time_strictly_between',
]
