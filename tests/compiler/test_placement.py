"""Phase 2 tests for the independent absolute-time placement stage."""

from __future__ import annotations

from types import SimpleNamespace

import pypulseq as pp
import pytest

import seqcraft as sc
import seqcraft.core.compiler as compiler
from seqcraft.core._compiler.model import PlacedEvent
from seqcraft.core._compiler.placement import place_events
from seqcraft.core.events import content_hash


def test_nested_offsets_keep_insertion_order_and_provenance(opts) -> None:
    """Traversal resolves absolute times without sorting equal or out-of-order siblings."""
    first = pp.make_delay(100e-6)
    second = pp.make_delay(200e-6)
    tail = pp.make_delay(300e-6)
    inner = sc.LogicBlock('inner').add(2e-3, first).add(1e-3, second)
    root = sc.LogicBlock('tr').add(3e-3, inner).add(0.0, tail)

    placed = place_events(root, opts)

    assert isinstance(placed, tuple)
    assert [event.event for event in placed] == [first, second, tail]
    assert [event.node_t for event in placed] == pytest.approx([5e-3, 4e-3, 0.0])
    assert [event.path for event in placed] == [
        ('tr', 'inner'),
        ('tr', 'inner'),
        ('tr',),
    ]


def test_active_and_reservation_intervals_match_each_event_kind(opts) -> None:
    """Delay, gradient, RF, ADC, and trigger timing retain the frozen interpretation."""
    delay = pp.make_delay(300e-6)
    trap = pp.make_trapezoid('x', area=100.0, duration=1e-3, delay=40e-6, system=opts)
    rf = pp.make_block_pulse(
        flip_angle=0.2,
        duration=400e-6,
        delay=200e-6,
        system=opts,
        use='excitation',
    )
    adc = pp.make_adc(num_samples=8, dwell=10e-6, delay=120e-6, system=opts)
    trigger = pp.make_trigger('physio1', duration=250e-6, delay=30e-6, system=opts)
    root = sc.LogicBlock('timing')
    for index, event in enumerate((delay, trap, rf, adc, trigger)):
        root.add(index * 2e-3, event)

    by_kind = {event.kind: event for event in place_events(root, opts)}

    placed_delay = by_kind['delay']
    assert placed_delay.start == placed_delay.node_t
    assert placed_delay.end == pytest.approx(placed_delay.node_t + float(delay.delay))
    assert (placed_delay.res_start, placed_delay.res_end) == (
        placed_delay.start,
        placed_delay.end,
    )

    placed_trap = by_kind['trap']
    trap_duration = float(trap.rise_time) + float(trap.flat_time) + float(trap.fall_time)
    assert placed_trap.start == pytest.approx(placed_trap.node_t + float(trap.delay))
    assert placed_trap.end == pytest.approx(placed_trap.start + trap_duration)
    assert (placed_trap.res_start, placed_trap.res_end) == (
        placed_trap.start,
        placed_trap.end,
    )

    placed_rf = by_kind['rf']
    assert placed_rf.start == pytest.approx(placed_rf.node_t + float(rf.delay))
    assert placed_rf.end == pytest.approx(placed_rf.start + float(rf.shape_dur))
    assert placed_rf.res_start == placed_rf.node_t
    assert placed_rf.res_end == pytest.approx(placed_rf.end + float(rf.ringdown_time))

    placed_adc = by_kind['adc']
    assert placed_adc.start == pytest.approx(placed_adc.node_t + float(adc.delay))
    assert placed_adc.end == pytest.approx(placed_adc.start + float(adc.num_samples) * float(adc.dwell))
    assert placed_adc.res_start == placed_adc.node_t
    assert placed_adc.res_end == pytest.approx(placed_adc.end + float(adc.dead_time))

    placed_trigger = by_kind['trigger']
    assert placed_trigger.start == pytest.approx(placed_trigger.node_t + float(trigger.delay))
    assert placed_trigger.end == pytest.approx(placed_trigger.start + float(trigger.duration))
    assert (placed_trigger.res_start, placed_trigger.res_end) == (
        placed_trigger.start,
        placed_trigger.end,
    )


def test_placement_does_not_mutate_source_events(opts) -> None:
    """The stage stores read-only event references and leaves their numeric content untouched."""
    events = (
        pp.make_trapezoid('x', area=100.0, system=opts),
        pp.make_adc(num_samples=16, dwell=10e-6, system=opts),
    )
    root = sc.LogicBlock('readout').add(0.0, *events)
    before = tuple(content_hash(event) for event in events)

    placed = place_events(root, opts)

    assert tuple(content_hash(event) for event in events) == before
    assert [event.event for event in placed] == list(events)


def test_placement_does_not_construct_a_pulseq_sequence(monkeypatch, opts) -> None:
    """The stage is independently callable and has no dependency on block emission."""

    def forbidden_sequence(*args, **kwargs):
        raise AssertionError

    monkeypatch.setattr(pp, 'Sequence', forbidden_sequence)
    root = sc.LogicBlock('t').add(0.0, pp.make_delay(1e-3))
    assert place_events(root, opts)[0].kind == 'delay'


def test_placement_error_carries_stage_and_nested_source_without_changing_text(opts) -> None:
    """Machine-readable context supplements the frozen human-facing error message."""
    bogus = SimpleNamespace(type='trapezoid', channel='x', delay=0.0)
    inner = sc.LogicBlock('spoiler').add(2e-3, bogus)
    root = sc.LogicBlock('tr').add(3e-3, inner)

    with pytest.raises(sc.CompileError) as caught:
        place_events(root, opts)

    assert str(caught.value).startswith("unknown event type 'trapezoid'.")
    assert 'tr.spoiler' in str(caught.value)
    assert getattr(caught.value, 'stage') == 'placement'
    assert getattr(caught.value, 'source_path') == ('tr', 'spoiler')


def test_compiler_facade_uses_the_single_placement_implementation(opts) -> None:
    """The migration alias points at the extracted function instead of retaining old traversal."""
    assert compiler._place is place_events
    placed = compiler._place(sc.LogicBlock('t').add(0.0, pp.make_delay(1e-3)), opts)
    assert isinstance(placed, tuple)
    assert all(isinstance(event, PlacedEvent) for event in placed)
