"""
Triggers: indivisible, but not exclusive.

A trigger holds its block open for the length of its pulse, so a boundary inside one is
meaningless -- it is a single hardware action, not a waveform.  But pulseq stores triggers as
``TRIGGERS`` extensions and accepts several in one block, so requiring a boundary *between* two
of them would invent a constraint the hardware does not have.

Before the reservation model separated those two properties, a trigger contributed nothing to
boundary selection at all: it held its block open while having no say in where the block ended,
so a physio trigger overlapping a readout failed with a shortfall error that blamed the module.
"""

from __future__ import annotations

import pypulseq as pp
import pytest

import seqcraft as sc


def _trigger(opts, duration: float, channel: str = 'physio1'):
    return pp.make_trigger(channel, duration=duration, system=opts)


def _triggers_of(block) -> list:
    """``get_block().trig`` is a list when present -- pulseq allows several per block."""
    trig = getattr(block, 'trig', None)
    if trig is None:
        return []
    return list(trig) if isinstance(trig, (list, tuple)) else [trig]


def _block_starts(compiled) -> dict[int, float]:
    """Absolute start time of each block, by 1-based block index."""
    out, t = {}, 0.0
    for i in sorted(compiled.seq.block_events):
        out[i] = t
        t += float(compiled.seq.block_durations[i])
    return out


def test_a_trigger_may_overlap_a_readout(opts) -> None:
    """
    The case that used to fail outright.

    A physio trigger spanning the gap between two ADCs is an ordinary cardiac-gated acquisition,
    and it used to raise 'block 0 spans 660.0 us but its events need 2500.0 us'.
    """
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = (
        sc.LogicBlock('t')
        .add(0.0, adc)
        .add(0.5e-3, _trigger(opts, 2e-3))
        .add(3e-3, adc)
    )
    out = sc.compile(tree, opts)
    assert out.check().ok, out.check()


def test_a_trigger_is_never_split(opts) -> None:
    """Its whole pulse must sit inside one block, so exactly one block carries it."""
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    trig = _trigger(opts, 2e-3)
    tree = sc.LogicBlock('t').add(0.0, adc).add(0.5e-3, trig).add(3e-3, adc)
    out = sc.compile(tree, opts)
    starts = _block_starts(out)
    carrying = [i for i in sorted(out.seq.block_events) if _triggers_of(out.seq.get_block(i))]
    assert len(carrying) == 1, f'trigger appears in {len(carrying)} blocks, must be exactly one'

    index = carrying[0]
    start, dur = starts[index], float(out.seq.block_durations[index])
    trig = _triggers_of(out.seq.get_block(index))[0]
    t0 = start + float(trig.delay)
    assert t0 >= start - 1e-9, 'the trigger must not begin before its block'
    assert t0 + float(trig.duration) <= start + dur + 1e-9, (
        f'the trigger runs to {(t0 + float(trig.duration)) * 1e6:.1f} us but its block ends at '
        f'{(start + dur) * 1e6:.1f} us'
    )
    assert t0 == pytest.approx(0.5e-3, abs=1e-9), 'and it must play when the tree said'


def test_two_triggers_need_no_boundary_between_them(opts) -> None:
    """
    Not exclusive: pulseq accepts several triggers per block, so no boundary is *required*.

    Asserted as "one block holds both" rather than "the sequence has one block" -- a trailing
    empty block from the trigger's own reservation edge is a separate matter (D7 / W6), and
    conflating the two would make this test fail for an unrelated reason.
    """
    tree = (
        sc.LogicBlock('t')
        .add(0.0, _trigger(opts, 1e-3, 'physio1'))
        .add(0.0, _trigger(opts, 1e-3, 'physio2'))
        .add(0.0, pp.make_delay(2e-3))
    )
    out = sc.compile(tree, opts)
    counts = [len(_triggers_of(out.seq.get_block(i))) for i in sorted(out.seq.block_events)]
    assert max(counts) == 2, f'both triggers must share one block, got {counts} per block'
    assert sum(counts) == 2, 'and neither may be duplicated'
    assert out.check().ok


def test_a_trigger_alongside_an_rf_and_a_gradient(opts) -> None:
    rf = pp.make_sinc_pulse(flip_angle=0.5, duration=1e-3, system=opts, use='excitation',
                            slice_thickness=5e-3, apodization=0.5, time_bw_product=4)
    tree = (
        sc.LogicBlock('t')
        .add(0.0, _trigger(opts, 200e-6))
        .add(0.0, rf)
        .add(1.5e-3, pp.make_trapezoid('z', area=200.0, system=opts))
    )
    out = sc.compile(tree, opts)
    assert out.check().ok, out.check()


def test_a_digital_output_is_treated_the_same(opts) -> None:
    """``output`` is the other half of the TRIGGERS extension and follows identical rules."""
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    out_pulse = pp.make_digital_output_pulse('osc0', duration=2e-3, system=opts)
    tree = sc.LogicBlock('t').add(0.0, adc).add(0.5e-3, out_pulse).add(3e-3, adc)
    compiled = sc.compile(tree, opts)
    assert compiled.check().ok, compiled.check()


def test_a_trigger_covering_a_whole_gap_is_a_clear_error(opts) -> None:
    """
    The genuinely impossible case, reported as such.

    Two ADCs must be in different blocks; a trigger spanning the entire gap between them cannot
    be cut and cannot be in both.  That is infeasible in pulseq, so the error names the trigger
    rather than surfacing later as a duration shortfall.
    """
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    # adc reservations are 0..660 us and 2000..2660 us; cover 600..2100 us entirely.
    tree = (
        sc.LogicBlock('t')
        .add(0.0, adc)
        .add(600e-6, _trigger(opts, 1500e-6))
        .add(2e-3, adc)
    )
    with pytest.raises(sc.CompileError) as err:
        sc.compile(tree, opts)
    text = str(err.value)
    assert 'nothing in the gap can be cut' in text
    assert 'trigger' in text, 'the message must name what is blocking the gap'


def test_a_barrier_inside_a_trigger_is_rejected(opts) -> None:
    """A barrier is an explicit request for a boundary; inside a trigger there cannot be one."""
    tree = (
        sc.LogicBlock('t')
        .add(0.0, _trigger(opts, 2e-3))
        .add(1e-3, sc.barrier('mid'))
        .add(0.0, pp.make_delay(3e-3))
    )
    with pytest.raises(sc.CompileError, match='cannot be split'):
        sc.compile(tree, opts)


def test_a_barrier_inside_an_rf_is_rejected(opts) -> None:
    """Same rule, and previously this surfaced as an unaligned-duration complaint instead."""
    rf = pp.make_sinc_pulse(flip_angle=0.5, duration=1e-3, system=opts, use='excitation',
                            slice_thickness=5e-3, apodization=0.5, time_bw_product=4)
    tree = sc.LogicBlock('t').add(0.0, rf).add(500e-6, sc.barrier('mid'))
    with pytest.raises(sc.CompileError) as err:
        sc.compile(tree, opts)
    assert 'cannot be split' in str(err.value)
    assert 'rf' in str(err.value)


def test_a_barrier_inside_an_adc_is_rejected(opts) -> None:
    adc = pp.make_adc(num_samples=256, dwell=10e-6, system=opts)
    tree = sc.LogicBlock('t').add(0.0, adc).add(1e-3, sc.barrier('mid'))
    with pytest.raises(sc.CompileError, match='cannot be split'):
        sc.compile(tree, opts)


def test_a_barrier_at_a_trigger_edge_is_fine(opts) -> None:
    """The rule is *strictly* inside; touching an edge splits nothing."""
    tree = (
        sc.LogicBlock('t')
        .add(0.0, _trigger(opts, 2e-3))
        .add(2e-3, sc.barrier('after'))
        .add(0.0, pp.make_delay(4e-3))
    )
    out = sc.compile(tree, opts)
    assert out.n_blocks == 2
    assert out.check().ok
