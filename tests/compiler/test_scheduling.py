"""
The compiler: the case table from the plan, one test each.

This is the heart of the suite.  The data model is small enough to be obviously correct; the
compiler is where every decision about overlap, splitting, merging and limits is made, so this is
where the tests go.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest
from pypulseq.opts import Opts

import seqcraft as sc
from seqcraft.design.events import content_hash


def compile_one(opts: Opts, *nodes: tuple[float, object]) -> sc.CompiledSequence:
    """Compile a flat tree of ``(start, event)`` pairs."""
    lb = sc.LogicBlock('t')
    for start, event in nodes:
        lb.add(start, event)
    return sc.compile(lb, opts)


def merges(out: sc.CompiledSequence) -> int:
    """How many same-axis merges the compile reported."""
    return len(out.report.of_kind('grad_merge'))


def limit_errors(out: sc.CompiledSequence) -> list[sc.Issue]:
    """Per-axis amplitude and slew violations (errors, not the vector-norm warnings)."""
    return [i for i in out.report.errors if i.kind.endswith('_limit')]


# ------------------------------------------------------------------------- overlap: the 3 rules
def test_different_axes_at_once_is_one_block_and_silent(opts) -> None:
    """
    The normal way to build a sequence: a slice rewinder, a phase blip and a readout prephaser.

    Warning about it would train people to ignore warnings, so nothing about the *overlap* may be
    reported.  `rise_time` is given explicitly because ``make_trapezoid`` with only a duration uses
    the shortest legal ramp, which puts each axis near the slew limit on its own -- three of those
    together do earn the vector-norm warning, and that is a different claim with its own test.
    """
    gentle = {'area': 100.0, 'duration': 2e-3, 'rise_time': 200e-6, 'system': opts}
    out = compile_one(
        opts,
        (0.0, pp.make_trapezoid('x', **gentle)),
        (0.0, pp.make_trapezoid('y', **gentle)),
        (0.0, pp.make_trapezoid('z', **gentle)),
    )
    assert out.n_blocks == 1
    assert merges(out) == 0
    assert out.report.issues == (), f'expected silence, got {out.report}'
    assert out.check().ok


def test_same_axis_at_once_merges_with_one_warning(opts) -> None:
    """Summing is what was meant, so it happens -- but it is the one waveform change, so it says so."""
    out = compile_one(
        opts,
        (0.0, pp.make_trapezoid('x', area=100.0, system=opts)),
        (0.0, pp.make_trapezoid('x', area=200.0, system=opts)),
    )
    assert out.n_blocks == 1
    assert merges(out) == 1
    assert out.moments()['x'] == pytest.approx(300.0, abs=1e-6)


def test_the_merge_warning_names_both_sources(opts) -> None:
    """A merge you did not expect has to be traceable to the two modules that caused it."""
    a = sc.LogicBlock('rewinder').add(0.0, pp.make_trapezoid('x', area=100.0, system=opts))
    b = sc.LogicBlock('prephaser').add(0.0, pp.make_trapezoid('x', area=-50.0, system=opts))
    out = sc.compile(sc.LogicBlock('tr').add(0.0, a).add(0.0, b), opts)
    message = out.report.of_kind('grad_merge')[0].message
    assert 'tr.rewinder' in message
    assert 'tr.prephaser' in message


def test_two_legal_gradients_can_sum_to_an_illegal_one(opts) -> None:
    """
    The reason limits are checked after merging rather than per module.

    Two 60 %-amplitude gradients on one axis are each perfectly legal and together are not, and no
    module can see that in isolation.
    """
    big = pp.make_trapezoid('x', amplitude=0.6 * opts.max_grad, duration=2e-3, system=opts)
    out = compile_one(opts, (0.0, big), (0.0, big))
    errors = limit_errors(out)
    assert errors, 'a 120% merged amplitude must be an error'
    assert any(i.kind == 'grad_limit' for i in errors)
    assert '120%' in ' '.join(i.message for i in errors)
    assert not out.check().ok


def test_merging_can_break_the_slew_limit_too(opts) -> None:
    """Amplitude is not the only thing a merge can break; on my first attempt slew hit 189 %."""
    out = compile_one(
        opts,
        (0.0, pp.make_trapezoid('x', area=100.0, system=opts)),
        (0.0, pp.make_trapezoid('x', area=200.0, system=opts)),
    )
    assert any(i.kind == 'slew_limit' for i in limit_errors(out))


def test_vector_norm_is_a_warning_not_an_error(opts) -> None:
    """
    Two axes ramping together reach sqrt(2) times the per-axis slew, which real amplifiers allow.

    Making this an error would reject the ordinary three-way winder overlap that the whole design
    exists to support.
    """
    strong = {'amplitude': 0.9 * opts.max_grad, 'duration': 1e-3, 'system': opts}
    out = compile_one(
        opts,
        (0.0, pp.make_trapezoid('x', **strong)),
        (0.0, pp.make_trapezoid('y', **strong)),
    )
    assert any(i.kind.endswith('_norm_limit') for i in out.report.warnings)
    assert not limit_errors(out)
    assert out.check().ok


# ---------------------------------------------------------------------------- RF/ADC exclusivity
def test_two_overlapping_rf_pulses_are_an_error(opts) -> None:
    rf = pp.make_sinc_pulse(flip_angle=1.57, duration=1e-3, system=opts, use='excitation')
    with pytest.raises(sc.CompileError, match='cannot play two pulses at once'):
        compile_one(opts, (0.0, rf), (0.5e-3, rf))


def test_two_overlapping_adcs_are_an_error(opts) -> None:
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    with pytest.raises(sc.CompileError, match='cannot open two sampling windows'):
        compile_one(opts, (0.0, adc), (0.5e-3, adc))


def test_transmitting_while_receiving_is_an_error(opts) -> None:
    rf = pp.make_sinc_pulse(flip_angle=1.57, duration=1e-3, system=opts, use='excitation')
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    with pytest.raises(sc.CompileError, match='cannot transmit and receive'):
        compile_one(opts, (0.0, adc), (0.3e-3, rf))


def test_the_exclusivity_error_names_both_paths_and_the_overlap(opts) -> None:
    """pypulseq would fail later, on a block index, with 'Multiple ADC events were specified'."""
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    a = sc.LogicBlock('readout').add(0.0, adc)
    b = sc.LogicBlock('navigator').add(0.5e-3, adc)
    with pytest.raises(sc.CompileError) as err:
        sc.compile(sc.LogicBlock('tr').add(0.0, a).add(0.0, b), opts)
    text = str(err.value)
    assert 'tr.readout' in text
    assert 'tr.navigator' in text
    assert 'overlap by' in text


def test_dead_time_conflicts_are_caught_even_when_waveforms_do_not_touch(opts) -> None:
    """
    An RF ending 5 us before an ADC starts looks fine and is not.

    The RF's ringdown and the ADC's dead time both sit inside their reservations, which is why the
    reservation rather than the waveform is what gets compared.
    """
    rf = pp.make_sinc_pulse(flip_angle=1.57, duration=1e-3, system=opts, use='excitation')
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    end_of_rf = float(rf.delay) + float(rf.shape_dur)
    with pytest.raises(sc.CompileError, match='transmit and receive'):
        compile_one(opts, (0.0, rf), (end_of_rf + 5e-6, adc))


# ------------------------------------------------------------------------ boundaries and splits
def test_a_gradient_spanning_an_rf_stays_one_event(opts) -> None:
    """
    pulseq allows one RF and one gradient per block, so nothing needs splitting here.

    Splitting anyway would turn a trapezoid into two arbitrary waveforms and cost two shapes.
    """
    rf = pp.make_sinc_pulse(flip_angle=1.57, duration=1e-3, system=opts, use='excitation')
    long_g = pp.make_trapezoid('x', area=2000.0, duration=4e-3, system=opts)
    out = compile_one(opts, (0.0, long_g), (1.5e-3, rf))
    assert out.n_blocks == 1
    assert out.seq.get_block(1).gx.type == 'trap'
    assert out.moments()['x'] == pytest.approx(2000.0, rel=1e-9)
    assert out.check().ok


def test_a_boundary_never_falls_inside_an_adc_window(opts) -> None:
    """Ramp sampling and vendor gridding both need the sampled gradient to be one event."""
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    long_g = pp.make_trapezoid('x', area=2000.0, duration=4e-3, system=opts)
    out = compile_one(opts, (0.0, long_g), (1.5e-3, adc))
    assert out.n_blocks == 1
    assert out.seq.get_block(1).gx.type == 'trap'
    assert out.check().ok


def test_a_split_preserves_area_and_continuity(opts) -> None:
    """When a boundary is forced, the two halves must join and their areas must add up."""
    long_g = pp.make_trapezoid('x', area=2000.0, duration=4e-3, system=opts)
    out = sc.compile(
        sc.LogicBlock('t').add(0.0, long_g).add(2e-3, sc.barrier('mid')), opts
    )
    assert out.n_blocks == 2
    assert out.moments()['x'] == pytest.approx(2000.0, rel=1e-9)
    first, second = out.seq.get_block(1).gx, out.seq.get_block(2).gx
    assert float(first.last) == pytest.approx(float(second.first), rel=1e-9)
    assert out.check().ok


def test_a_split_mid_ramp_keeps_the_slew(opts) -> None:
    """The seam must not become a kink."""
    ramp = pp.make_trapezoid('x', amplitude=0.5 * opts.max_grad, duration=2e-3, system=opts)
    out = sc.compile(
        sc.LogicBlock('t').add(0.0, ramp).add(float(ramp.rise_time) / 2.0, sc.barrier()), opts
    )
    assert out.n_blocks == 2
    assert not limit_errors(out)
    assert out.check().ok


def test_barrier_forces_a_boundary_and_costs_no_time(opts) -> None:
    g = pp.make_trapezoid('x', area=500.0, duration=2e-3, system=opts)
    plain = compile_one(opts, (0.0, g))
    split = sc.compile(sc.LogicBlock('t').add(0.0, g).add(1e-3, sc.barrier()), opts)
    assert split.n_blocks == plain.n_blocks + 1
    assert split.duration_s == pytest.approx(plain.duration_s)


def test_a_delay_only_block_compiles(opts) -> None:
    """The b=0 diffusion volume: correct duration, no events."""
    out = compile_one(opts, (0.0, pp.make_delay(4.2e-3)))
    assert out.duration_s == pytest.approx(4.2e-3)
    assert out.n_blocks == 1
    assert out.check().ok


def test_an_over_long_interval_is_subdivided(opts) -> None:
    """pulseq stores a block duration in a fixed-width field."""
    out = compile_one(opts, (0.0, pp.make_delay(1.0)))
    assert out.duration_s == pytest.approx(1.0)
    assert out.check().ok


# ------------------------------------------------------------------------------ error reporting
def test_a_negative_start_is_an_error_naming_the_likely_cause(opts) -> None:
    """The usual cause is a TE shorter than the readout needs, and the message says so."""
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    with pytest.raises(sc.CompileError, match='before the start of the sequence'):
        compile_one(opts, (-1e-3, g))


def test_an_empty_tree_is_an_error(opts) -> None:
    with pytest.raises(sc.CompileError, match='nothing to compile'):
        sc.compile(sc.LogicBlock('empty'), opts)


@pytest.mark.parametrize('point_kind', ['barrier', 'label'])
def test_a_tree_with_only_zero_duration_points_is_an_error(opts, point_kind) -> None:
    """Point events cannot create a physical block without an event that occupies time."""
    point = sc.barrier('only') if point_kind == 'barrier' else pp.make_label('LIN', 'SET', 0)
    tree = sc.LogicBlock('point_only').add(0.0, point)
    with pytest.raises(sc.CompileError, match='nothing to compile'):
        sc.compile(tree, opts)


def test_compile_does_not_mutate_the_input_tree_or_events(opts) -> None:
    """Compilation may register derived events, but source nodes and numeric content stay intact."""
    gradient = pp.make_trapezoid('x', area=100.0, system=opts)
    adc = pp.make_adc(num_samples=64, dwell=4e-6, system=opts)
    label = pp.make_label('LIN', 'SET', 7)
    inner = sc.LogicBlock('inner').add(0.0, gradient).add(1e-3, adc, label)
    tree = sc.LogicBlock('root').add(2e-3, inner)

    before = [
        (start, path, id(event), content_hash(event))
        for start, event, path in sc.flatten(tree)
    ]
    out = sc.compile(tree, opts)
    after = [
        (start, path, id(event), content_hash(event))
        for start, event, path in sc.flatten(tree)
    ]

    assert after == before
    assert out.check().ok


# ---------------------------------------------------------------------------------- invariants
def test_compiled_duration_equals_the_tree_duration(opts) -> None:
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g).add(5e-3, g).add(20e-3, pp.make_delay(1e-3))
    out = sc.compile(tree, opts)
    assert out.duration_s == pytest.approx(21e-3)
    assert not out.report.of_kind('duration')


def test_per_axis_m0_survives_compilation(opts) -> None:
    """
    The invariant that catches a split which lost a tail or a merge which dropped a piece.

    It runs on every compile, so no individual test case has to think of the combination that
    would break.
    """
    tree = sc.LogicBlock('t')
    tree.add(0.0, pp.make_trapezoid('x', area=100.0, system=opts))
    tree.add(0.0, pp.make_trapezoid('y', area=-250.0, system=opts))
    tree.add(1e-3, pp.make_trapezoid('x', area=-40.0, system=opts))
    tree.add(2e-3, sc.barrier())
    out = sc.compile(tree, opts)
    assert out.moments()['x'] == pytest.approx(60.0, abs=1e-6)
    assert out.moments()['y'] == pytest.approx(-250.0, abs=1e-6)
    assert not out.report.of_kind('moment')


def test_block_durations_land_on_the_block_raster(opts) -> None:
    """
    Every block duration must be an exact multiple of the block raster.

    The failure mode this guards is subtle: nearest-rounding a reservation *end* can leave a block
    a few microseconds short of what its ADC needs, and pypulseq then reports an unaligned duration
    rather than the missing microseconds it actually is.
    """
    rf = pp.make_sinc_pulse(flip_angle=1.57, duration=1e-3, system=opts, use='excitation')
    adc = pp.make_adc(num_samples=2372, dwell=2e-6, system=opts)
    g = pp.make_trapezoid('x', area=400.0, system=opts)
    out = sc.compile(
        sc.LogicBlock('t').add(0.0, rf).add(2e-3, g).add(2e-3, adc), opts
    )
    raster = sc.Raster(opts.block_duration_raster)
    for index, duration in out.seq.block_durations.items():
        assert raster.holds(duration), f'block {index} is {duration * 1e6} us'
    assert out.check().ok


def test_event_delays_land_on_their_own_raster(opts) -> None:
    """
    RF on 1 us, ADC on 100 ns, gradients on 10 us -- and computed in integer picoseconds.

    A plain subtraction of absolute times drifts: at 39 s into a sequence it produced an RF delay
    of 129.9999999986 us, which pypulseq rejects.
    """
    rf = pp.make_sinc_pulse(flip_angle=1.57, duration=1e-3, system=opts, use='excitation')
    tree = sc.LogicBlock('t')
    for i in range(40):
        tree.add(i * 1.997e-3 + 3.0, rf)
    out = sc.compile(tree, opts)
    for index in out.seq.block_events:
        block = out.seq.get_block(index)
        if getattr(block, 'rf', None) is not None:
            assert sc.Raster(opts.rf_raster_time).holds(float(block.rf.delay))
    assert out.check().ok


# ---------------------------------------------------------------------------------- provenance
def test_origin_traces_a_block_back_to_its_module(opts) -> None:
    inner = sc.LogicBlock('spoiler').add(0.0, pp.make_trapezoid('z', area=500.0, system=opts))
    out = sc.compile(sc.LogicBlock('tr').add(0.0, inner), opts)
    assert out.origin(0) == ('tr', 'spoiler')


def test_origin_of_a_shared_block_is_the_common_ancestor(opts) -> None:
    """
    Three modules in one block have no single origin, and saying otherwise would mislead.

    The honest answer is the path they share.
    """
    tree = sc.LogicBlock('tr')
    tree.add(0.0, sc.LogicBlock('rephaser').add(0.0, pp.make_trapezoid('z', area=100.0, system=opts)))
    tree.add(0.0, sc.LogicBlock('blip').add(0.0, pp.make_trapezoid('y', area=100.0, system=opts)))
    out = sc.compile(tree, opts)
    assert out.origin(0) == ('tr',)


# -------------------------------------------------------------------------------------- labels
def test_duplicate_kspace_addresses_are_an_error(opts) -> None:
    """
    Two readouts writing the same k-space location is the single highest-value finished-sequence
    check: it catches a wrong slice order, an off-by-one partial-Fourier start, and a mis-nested
    loop, all from one assertion.
    """
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = sc.LogicBlock('t')
    for i in range(3):
        tree.add(i * 5e-3, adc, pp.make_label('LIN', 'SET', 7))
    out = sc.compile(tree, opts)
    report = out.check()
    assert not report.ok
    assert any(i.kind == 'label' for i in report.errors)


def test_unique_kspace_addresses_pass(opts) -> None:
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = sc.LogicBlock('t')
    for i in range(3):
        tree.add(i * 5e-3, adc, pp.make_label('LIN', 'SET', i))
    assert sc.compile(tree, opts).check().ok


def test_a_label_attaches_to_the_block_containing_its_start(opts) -> None:
    """A readout's LIN must land on the block holding the ADC, which is what the recon reads."""
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    out = sc.compile(
        sc.LogicBlock('t').add(0.0, g).add(1e-3, adc, pp.make_label('LIN', 'SET', 5)), opts
    )
    labels = out.seq.evaluate_labels(evolution='adc')
    assert int(np.atleast_1d(labels['LIN'])[0]) == 5


# ------------------------------------------------------------------------------- shape economy
def test_merging_two_trapezoids_stays_a_short_shape(opts) -> None:
    """
    Reducing collinear points first is what keeps a merged trapezoid cheap.

    Without it a summed trapezoid would be hundreds of raster samples, and a long sequence would
    grow a shape library out of proportion to what it plays.
    """
    out = compile_one(
        opts,
        (0.0, pp.make_trapezoid('x', area=100.0, duration=2e-3, system=opts)),
        (0.0, pp.make_trapezoid('x', area=200.0, duration=2e-3, system=opts)),
    )
    grad = out.seq.get_block(1).gx
    if grad.type == 'grad':
        assert len(grad.tt) <= 8, f'merged shape has {len(grad.tt)} points'


def test_a_lone_arbitrary_gradient_passes_through_untouched(opts) -> None:
    """
    A spiral must not be resampled: it takes the fast path when it sits alone in its interval.

    Resampling would be lossy exactly where the diffusion measurement is most sensitive.
    """
    wave = np.sin(np.linspace(0.0, np.pi, 500)) * 0.5 * opts.max_grad
    g = pp.make_arbitrary_grad(channel='x', waveform=wave, system=opts)
    out = compile_one(opts, (0.0, g))
    assert out.seq.get_block(1).gx.type == 'grad'
    assert len(out.seq.get_block(1).gx.waveform) == len(wave)
