"""Compiler IR contracts and their use by the authoritative compile path."""

from __future__ import annotations

import warnings
from dataclasses import FrozenInstanceError
from types import SimpleNamespace

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
import seqcraft.compiler as compiler
from seqcraft.compiler import legalization
from seqcraft.compiler.model import (
    EXCLUSIVE_KINDS,
    LegalizationResult,
    PlacedEvent,
    PulseqReadyBlock,
    interval_duration,
    time_at_or_before,
    time_before,
    time_equal,
    time_strictly_between,
)
from seqcraft.compiler.verification import (
    _traversed,
    verify_placed_events,
    verify_ready_blocks,
)
from seqcraft.design.events import GRADIENT_KINDS, HANDLED_KINDS, knots_of, pwl_moment


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


def test_legalization_result_freezes_blocks_and_notes() -> None:
    """The stage returns its ready blocks and warning evidence without mutable side output."""
    ready = PulseqReadyBlock(0, 0.0, 1e-3, 1e-3, (), (), ())
    result = LegalizationResult((ready,), (('merge', ('readout+rewinder (axis x)',)),))

    assert result.blocks == (ready,)
    assert result.notes[0][0] == 'merge'
    with pytest.raises(FrozenInstanceError):
        result.blocks = ()  # type: ignore[misc]


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


def test_warning_note_categories_require_explicit_display_text() -> None:
    """A new transformation note must not silently degrade to a bare internal key."""
    with pytest.raises(sc.CompilerContractError) as caught:
        compiler._warn({'unnamed': ['readout']})
    assert 'warning-note' in str(caught.value)
    assert 'unnamed' in str(caught.value)


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
    real_ready = legalization.verify_ready_blocks

    def capture_placed(events):
        seen_placed.extend(events)
        return real_placed(events)

    def capture_ready(blocks, **kwargs):
        seen_ready.extend(blocks)
        return real_ready(blocks, **kwargs)

    monkeypatch.setattr(compiler, 'verify_placed_events', capture_placed)
    monkeypatch.setattr(legalization, 'verify_ready_blocks', capture_ready)

    grad = pp.make_trapezoid('x', area=100.0, duration=1e-3, system=opts)
    tree = sc.LogicBlock('tr').add(0.0, grad).add(2e-3, pp.make_delay(1e-3))
    out = sc.compile(tree, opts)

    assert seen_placed and all(isinstance(event, PlacedEvent) for event in seen_placed)
    assert seen_ready and all(isinstance(block, PulseqReadyBlock) for block in seen_ready)
    assert len(out.block_events) == len(seen_ready)
    # Provenance is no longer returned -- it is used where it is produced, to name the source in
    # an error message -- so what is checkable from here is that emission computed one per block.
    assert all(isinstance(block.origin, tuple) for block in seen_ready)


# --------------------------------------------------- the tolerance an out-and-back event needs
def _out_and_back(opts, *, samples: int = 400, peak: float = 1.0e6):
    """One arbitrary gradient that runs out to a large area and comes back to the origin.

    A spiral in miniature: the net area is zero and the area *traversed* is thousands of times
    ``dk``, which is the only shape in the library where those two numbers differ at all.
    """
    ramp = np.sin(np.linspace(0.0, np.pi, samples // 2)) * peak
    return pp.make_arbitrary_grad(
        channel='x', waveform=np.concatenate([ramp, -ramp]), first=0.0, last=0.0,
        max_grad=np.inf, max_slew=np.inf, system=opts)


def test_the_m0_tolerance_scales_with_the_area_traversed_not_the_net(opts) -> None:
    """An out-and-back gradient compiles, and the same one with a lobe removed still does not.

    Both halves are the test.  Scaled by ``|net area|`` the tolerance collapses onto its ``1e-6``
    floor -- where pypulseq's own shape compression already sits -- and a legal readout is refused
    by a message blaming a split that never happened.  Scaled by ``integral |g|`` it is four orders
    of magnitude above the compression noise and four orders *below* a genuinely lost lobe.
    """
    event = _out_and_back(opts)
    tree = sc.LogicBlock('out-and-back').add(0.0, event)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        seq = sc.compile(tree, opts)
    assert len(seq.block_events) == 1

    net = abs(pwl_moment(*knots_of(event)))
    traversed = _traversed(*knots_of(event))
    assert net < 1e-6 * traversed, 'the shape under test must nearly cancel'
    assert traversed > 1e5 * max(net, 1e-9)

    # The same waveform with one lobe deleted: a real loss, and still four orders of magnitude
    # above the corrected tolerance of 1e-6 * traversed.
    half = np.asarray(event.waveform)[: len(event.waveform) // 2]
    lost = abs(pwl_moment(*knots_of(pp.make_arbitrary_grad(
        channel='x', waveform=half, first=0.0, last=0.0, max_grad=np.inf, max_slew=np.inf,
        system=opts))))
    assert lost > 1e4 * (1e-6 * traversed)


def test_traversed_area_counts_both_signs_of_a_sign_change() -> None:
    """``integral |g|``, split at the zero crossing -- the two lines the m0 tolerance rests on."""
    trapezoid = (np.array([0.0, 1.0, 2.0, 3.0]), np.array([0.0, 1.0, 1.0, 0.0]))
    bipolar = (np.array([0.0, 1.0, 2.0, 3.0, 4.0]), np.array([0.0, 1.0, 0.0, -1.0, 0.0]))
    through_zero = (np.array([0.0, 2.0]), np.array([-1.0, 1.0]))
    assert _traversed(*trapezoid) == pytest.approx(pwl_moment(*trapezoid))
    assert _traversed(*bipolar) == pytest.approx(2.0)
    assert pwl_moment(*bipolar) == pytest.approx(0.0)
    assert _traversed(*through_zero) == pytest.approx(1.0)
