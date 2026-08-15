"""Immutable compiler contracts shared by placement, legalization, and emission.

The contracts deliberately hold references to PyPulseq events instead of copying them.  SeqCraft
treats those event objects as read-only during compilation; tuples freeze the compiler-owned
structure without pretending that third-party ``SimpleNamespace`` instances are deeply immutable.

**Block-format policy only.**  What an event *is* -- which types carry a gradient, which are
instants, which are labels, which the compiler handles at all -- is event identity rather than a
compiler decision, and lives in :mod:`seqcraft.design.events` beside the functions that read pulseq's
``type`` field.  The two constants here are the ones the compiler genuinely owns: they answer *may
a boundary fall inside this* and *may two of these share a block*, which are properties of the
pulseq block format, not of the event.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from ..design.timing import EPS, Raster, exact_diff

if TYPE_CHECKING:
    from pypulseq.opts import Opts

#: No block boundary may fall strictly inside one of these: each is a single hardware action --
#: an RF's dead time and ringdown, an ADC's window and trailing dead time, a trigger's pulse --
#: rather than a waveform, so a cut inside it has no meaning.
INDIVISIBLE_KINDS = frozenset({'rf', 'adc', 'trigger', 'output'})

#: A block holds at most one of these, so a boundary is *required* between two.  RF and ADC only:
#: triggers are ``TRIGGERS`` extensions and pulseq accepts several per block, so treating them as
#: exclusive would invent a constraint the hardware does not have.
EXCLUSIVE_KINDS = frozenset({'rf', 'adc'})


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


def in_block_delay(p: PlacedEvent, block_start: float, opts: Opts) -> float:
    """
    Return an event's delay within its block, quantised onto that event's own raster.

    Two reasons this cannot be a plain subtraction.  The absolute times come from arithmetic over
    a sequence that may run for minutes, so by the last TR the float resolution is coarser than a
    picosecond and ``p.start - block_start`` drifts -- pypulseq then reports an RF delay of
    ``129.9999999986us`` and rejects the block.  And pulseq requires each event's delay to sit on
    its own raster -- 1 us for RF, 100 ns for ADC and 10 us for gradients on Siemens, whatever
    the scanner reports elsewhere.  Subtracting in integer ticks and then snapping satisfies both.

    Time policy rather than a stage: legalization needs it to delay a gradient that passes
    through untouched, and emission needs it for every RF, ADC and label it places.  Owning it
    here is what keeps those two stages in one dependency direction.
    """
    dt = max(0.0, exact_diff(p.start, block_start))
    raster = {
        'rf': Raster(float(opts.rf_raster_time), 'RF'),
        'adc': Raster(float(opts.adc_raster_time), 'ADC'),
    }.get(p.kind, Raster(float(opts.grad_raster_time), 'gradient'))
    return raster.nearest(dt)


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
    'EXCLUSIVE_KINDS',
    'INDIVISIBLE_KINDS',
    'PlacedEvent',
    'in_block_delay',
    'PulseqReadyBlock',
    'interval_duration',
    'time_at_or_before',
    'time_before',
    'time_equal',
    'time_strictly_between',
]
