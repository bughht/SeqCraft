"""Phase 1 compiler contracts and their adapters into the authoritative compile path."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from types import SimpleNamespace

import pypulseq as pp
import pytest

import seqcraft as sc
import seqcraft.compiler as compiler
from seqcraft.compiler import emission
from seqcraft.compiler.model import (
    EXCLUSIVE_KINDS,
    PlacedEvent,
    PulseqReadyBlock,
    interval_duration,
    time_at_or_before,
    time_before,
    time_equal,
    time_strictly_between,
)
from seqcraft.compiler.verification import (
    verify_placed_events,
    verify_ready_blocks,
)
from seqcraft.design.events import GRADIENT_KINDS, HANDLED_KINDS


def _event(kind: str) -> SimpleNamespace:
    """Build the smallest event payload needed by the contract tests."""
    return SimpleNamespace(type=kind)


def test_placed_event_is_immutable_and_has_concise_provenance() -> None:
    """The first IR freezes compiler-owned fields without expanding event payloads in repr."""
    placed = PlacedEvent(
        node_t=1e-3,
        start=1.1e-3,
        end=1.4e-3,
        res_start=1e-3,
        res_end=1.5e-3,
        event=_event('rf'),
        path=('tr', 'excitation'),
    )

    assert placed.kind == 'rf'
    assert placed.where == 'tr.excitation'
    assert placed.source_path == ('tr', 'excitation')
    assert placed.duration == pytest.approx(300e-6)
    assert placed.reservation_duration == pytest.approx(500e-6)
    assert repr(placed) == (
        'PlacedEvent(rf 1100.0..1400.0 us ' '(reserved 1000.0..1500.0 us) from tr.excitation)'
    )
    with pytest.raises(FrozenInstanceError):
        placed.start = 0.0  # type: ignore[misc]


def test_ready_block_is_immutable_and_summarises_events() -> None:
    """The second IR carries stable timing and source paths immediately before emission."""
    ready = PulseqReadyBlock(
        index=2,
        start=1e-3,
        end=2e-3,
        duration=1e-3,
        events=(_event('adc'), _event('grad')),
        source_paths=(('tr', 'readout'), ('tr', 'readout')),
        origin=('tr', 'readout'),
    )

    assert ready.kinds == ('adc', 'grad')
    assert repr(ready) == ('PulseqReadyBlock(block 2 1000.0..2000.0 us [adc,grad] from tr.readout)')
    with pytest.raises(FrozenInstanceError):
        ready.duration = 0.0  # type: ignore[misc]


def test_event_classification_is_central_and_explicit() -> None:
    """Every stage sees the same hardware-oriented kind groups."""
    assert GRADIENT_KINDS == {'trap', 'grad'}
    assert EXCLUSIVE_KINDS == {'rf', 'adc'}
    assert {'rf', 'adc', 'trap', 'grad', 'trigger', 'output'} <= HANDLED_KINDS
    assert 'trigger' not in EXCLUSIVE_KINDS


def test_time_policy_distinguishes_tolerance_from_exact_tick_arithmetic() -> None:
    """Comparisons use EPS while interval arithmetic uses integer ticks."""
    assert time_equal(1e-3, 1e-3 + 0.5e-9)
    assert not time_before(1e-3, 1e-3 + 0.5e-9)
    assert time_before(1e-3, 1e-3 + 2e-9)
    assert time_at_or_before(1e-3 + 0.5e-9, 1e-3)
    assert time_strictly_between(2e-3, 1e-3, 3e-3)
    assert interval_duration(8e-3, 12.3e-3) == 4.3e-3


def test_placed_verifier_reports_reversed_intervals() -> None:
    """Contract verification returns structured findings instead of mutating the IR."""
    placed = PlacedEvent(0.0, 2e-3, 1e-3, 0.0, 3e-3, _event('adc'), ('readout',))
    violations = verify_placed_events((placed,))
    assert [violation.message for violation in violations] == ['active interval is reversed']


def test_ready_verifier_reports_timing_and_exclusivity() -> None:
    """The skeleton catches stage bugs before PyPulseq receives a malformed block."""
    ready = PulseqReadyBlock(
        index=0,
        start=0.0,
        end=1e-3,
        duration=2e-3,
        events=(_event('adc'), _event('adc')),
        source_paths=(('a',), ('b',)),
        origin=(),
    )
    messages = [violation.message for violation in verify_ready_blocks((ready,))]
    assert 'duration does not match start/end' in messages
    assert 'contains more than one ADC event' in messages


def test_authoritative_compile_path_produces_both_contracts(monkeypatch, opts) -> None:
    """
    The compile path really does produce both IRs, rather than running a second algorithm beside
    them.

    Patched in two places because the two contracts are checked by two stages: the pass verifies
    the placed events it receives, and emission verifies each ready block as it builds it.
    """
    seen_placed: list[PlacedEvent] = []
    seen_ready: list[PulseqReadyBlock] = []
    real_placed = compiler.verify_placed_events
    real_ready = emission.verify_ready_blocks

    def capture_placed(events):
        seen_placed.extend(events)
        return real_placed(events)

    def capture_ready(blocks, **kwargs):
        seen_ready.extend(blocks)
        return real_ready(blocks, **kwargs)

    monkeypatch.setattr(compiler, 'verify_placed_events', capture_placed)
    monkeypatch.setattr(emission, 'verify_ready_blocks', capture_ready)

    grad = pp.make_trapezoid('x', area=100.0, duration=1e-3, system=opts)
    tree = sc.LogicBlock('tr').add(0.0, grad).add(2e-3, pp.make_delay(1e-3))
    out = sc.compile(tree, opts)

    assert seen_placed and all(isinstance(event, PlacedEvent) for event in seen_placed)
    assert seen_ready and all(isinstance(block, PulseqReadyBlock) for block in seen_ready)
    assert len(out.block_events) == len(seen_ready)
    # Provenance is no longer returned -- it is used where it is produced, to name the source in
    # an error message -- so what is checkable from here is that emission computed one per block.
    assert all(isinstance(block.origin, tuple) for block in seen_ready)
