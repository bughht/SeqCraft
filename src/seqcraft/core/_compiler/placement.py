"""Resolve a :class:`LogicBlock` tree into immutable absolute-time events.

Placement owns tree traversal, event duration interpretation, hardware reservation spans, and
unsupported-kind rejection.  It does not select block boundaries, inspect downstream block state,
or construct a PyPulseq sequence.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ..errors import CompileError, format_error
from ..logic import BARRIER, flatten
from .model import HANDLED_KINDS, PlacedEvent

if TYPE_CHECKING:
    from types import SimpleNamespace

    from pypulseq.opts import Opts

    from ..logic import LogicBlock

# Recognized PyPulseq kinds that require explicit rejection rather than an unknown-kind fallback.
UNSUPPORTED_KINDS: dict[str, tuple[str, tuple[str, ...]]] = {
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
        # Scanner inputs and outputs are pulses whose duration holds their containing block open.
        return float(getattr(event, 'duration', 0.0) or 0.0)
    return 0.0


def _with_context(error: CompileError, path: tuple[str, ...]) -> CompileError:
    """Attach machine-readable stage and source context without changing the frozen message."""
    setattr(error, 'stage', 'placement')
    setattr(error, 'source_path', path)
    return error


def _unsupported(kind: str | None, path: tuple[str, ...], node_t: float) -> CompileError:
    """Build the placement error for an event kind the compiler will not emit."""
    where = '.'.join(path) if path else '(untagged)'
    if kind in UNSUPPORTED_KINDS:
        what, hints = UNSUPPORTED_KINDS[kind]
        return _with_context(
            CompileError(
                format_error(
                    f'{what} ({kind!r}) is not supported by the compiler.',
                    {'from': where, 'at': f'{node_t * 1e6:.1f} us'},
                    list(hints),
                )
            ),
            path,
        )
    return _with_context(
        CompileError(
            format_error(
                f'unknown event type {kind!r}.',
                {'from': where, 'at': f'{node_t * 1e6:.1f} us', 'handled': ', '.join(sorted(HANDLED_KINDS))},
                [
                    'LogicBlock.add() accepts any object with a .type attribute, so a typo or a '
                    'hand-built namespace reaches the compiler unchecked',
                    'if this is a real pulseq event from a newer pypulseq, seqcraft needs '
                    'updating -- please open an issue naming the type',
                ],
            )
        ),
        path,
    )


def place_events(root: LogicBlock, opts: Opts) -> tuple[PlacedEvent, ...]:
    """Flatten ``root`` into ordered absolute active and reservation intervals.

    A delay event occupies ``[node_t, node_t + delay]`` because its stored delay is its duration.
    Other events start after their own leading delay. RF and ADC reservations additionally include
    their hardware dead-time/ringdown spans.

    Raises
    ------
    CompileError
        If a leaf has an unsupported or unknown event kind. The exception keeps its frozen text and
        carries private ``stage='placement'`` and ``source_path`` attributes.
    """
    placed: list[PlacedEvent] = []
    for node_t, event, path in flatten(root):
        kind = getattr(event, 'type', None)
        if kind not in HANDLED_KINDS:
            raise _unsupported(kind, path, node_t)
        if kind == BARRIER:
            placed.append(PlacedEvent(node_t, node_t, node_t, node_t, node_t, event, path))
            continue
        delay = float(getattr(event, 'delay', 0.0) or 0.0)
        if kind == 'delay':
            placed.append(PlacedEvent(node_t, node_t, node_t + delay, node_t, node_t + delay, event, path))
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
        placed.append(PlacedEvent(node_t, start, end, res_start, res_end, event, path))
    return tuple(placed)


__all__ = ['UNSUPPORTED_KINDS', 'place_events']
