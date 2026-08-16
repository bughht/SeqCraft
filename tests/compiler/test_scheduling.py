"""
The compiler: the case table from the plan, one test each.

This is the heart of the suite.  The data model is small enough to be obviously correct; the
compiler is where every decision about overlap, splitting, merging and limits is made, so this is
where the tests go.
"""

from __future__ import annotations

import warnings

import numpy as np
import pypulseq as pp
import pytest
from fidelity import compiled_knots
from pypulseq.opts import Opts

import seqcraft as sc
from seqcraft.compiler.emission import common_path
from seqcraft.design.events import content_hash, trapz


def compile_one(opts: Opts, *nodes: tuple[float, object]):
    """Compile a flat tree of ``(start, event)`` pairs."""
    lb = sc.LogicBlock('t')
    for start, event in nodes:
        lb.add(start, event)
    return sc.compile(lb, opts)


def compiled_m0(out, axis: str) -> float:
    """
    The area an axis actually plays, in 1/m, read off the compiled blocks.

    Via the independent oracle in ``fidelity.py`` rather than the compiler's own arithmetic:
    a check that shares code with what it checks compares a number with itself.  m0 is the
    integral of a piecewise-linear function, so the trapezoidal rule is exact here.
    """
    times, amps = compiled_knots(out)[axis]
    return float(trapz(amps, np.asarray(times, dtype=float) * 1e-12))


def warned(kind: str, make) -> list[str]:
    """
    Return the messages of every :class:`sc.SeqCraftWarning` matching `kind` that `make` emitted.

    A recorder rather than ``pytest.warns``, because half the claims here are *negative* -- no
    merge warning for three axes at once -- and ``pytest.warns`` has no clean negative form.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        make()
    return [
        str(w.message) for w in caught
        if issubclass(w.category, sc.SeqCraftWarning) and kind in str(w.message)
    ]


#: A same-axis pair that sums to a legal waveform.  ``rise_time`` is stated because
#: ``make_trapezoid`` with only a duration uses the shortest legal ramp, which puts each lobe near
#: the slew limit on its own -- so their sum is 158 % of it and the compile raises.
GENTLE = {'duration': 2e-3, 'rise_time': 400e-6}


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

    def build():
        return compile_one(
            opts,
            (0.0, pp.make_trapezoid('x', **gentle)),
            (0.0, pp.make_trapezoid('y', **gentle)),
            (0.0, pp.make_trapezoid('z', **gentle)),
        )

    assert len(build().block_events) == 1
    assert warned('', build) == [], 'expected silence about three axes at one instant'


def test_same_axis_at_once_merges_with_one_warning(opts) -> None:
    """Summing is what was meant, so it happens -- but it is the one waveform change, so it says so."""
    def build():
        return compile_one(
            opts,
            (0.0, pp.make_trapezoid('x', area=100.0, system=opts, **GENTLE)),
            (0.0, pp.make_trapezoid('x', area=200.0, system=opts, **GENTLE)),
        )

    out = build()
    assert len(out.block_events) == 1
    assert compiled_m0(out, 'x') == pytest.approx(300.0, abs=1e-6)
    with pytest.warns(sc.SeqCraftWarning, match='1 same-axis gradient merge:'):
        build()


def test_a_run_of_merges_is_one_warning_not_one_each(opts) -> None:
    """
    Python's default filter shows a warning once per source line.

    One ``warn`` per merge would print the first and swallow the other eleven, which is strictly
    worse than the count the report used to give.  So the compile emits exactly one, carrying the
    count and the sites.
    """
    tree = sc.LogicBlock('tr')
    for i in range(8):
        tree.add(
            i * 5e-3,
            sc.LogicBlock('rewinder').add(
                0.0, pp.make_trapezoid('x', area=100.0, system=opts, **GENTLE)),
            sc.LogicBlock('prephaser').add(
                0.0, pp.make_trapezoid('x', area=-50.0, system=opts, **GENTLE)),
        )
    messages = warned('merge', lambda: sc.compile(tree, opts))
    assert len(messages) == 1, f'expected exactly one warning, got {messages}'
    assert messages[0].startswith('8 same-axis gradient merges'), messages[0]


def test_the_merge_warning_names_both_sources(opts) -> None:
    """A merge you did not expect has to be traceable to the two modules that caused it."""
    a = sc.LogicBlock('rewinder').add(
        0.0, pp.make_trapezoid('x', area=100.0, system=opts, **GENTLE))
    b = sc.LogicBlock('prephaser').add(
        0.0, pp.make_trapezoid('x', area=-50.0, system=opts, **GENTLE))
    tree = sc.LogicBlock('tr').add(0.0, a).add(0.0, b)
    message = warned('merge', lambda: sc.compile(tree, opts))[0]
    assert 'tr.rewinder' in message
    assert 'tr.prephaser' in message


def test_two_legal_gradients_can_sum_to_an_illegal_one(opts) -> None:
    """
    The reason limits are checked after merging rather than per module.

    Two 60 %-amplitude gradients on one axis are each perfectly legal and together are not, and no
    module can see that in isolation.  There is no legal sequence to hand back, so the compile
    raises rather than returning one with a note attached.
    """
    big = pp.make_trapezoid('x', amplitude=0.6 * opts.max_grad, duration=2e-3, system=opts)
    with pytest.raises(sc.HardwareLimitError) as err:
        compile_one(opts, (0.0, big), (0.0, big))
    text = str(err.value)
    assert text.startswith('gradient 120% of the 40 mT/m limit on axis x.'), text
    assert 'derate' in text, 'the message must offer a way forward'


def test_merging_can_break_the_slew_limit_too(opts) -> None:
    """Amplitude is not the only thing a merge can break; on my first attempt slew hit 189 %."""
    with pytest.raises(sc.HardwareLimitError, match='slew 189% of the 150 T/m/s limit on axis x'):
        compile_one(
            opts,
            (0.0, pp.make_trapezoid('x', area=100.0, system=opts)),
            (0.0, pp.make_trapezoid('x', area=200.0, system=opts)),
        )


def test_the_limit_error_names_the_time_and_the_source(opts) -> None:
    """A block index means nothing until the sequence is written; a time can be acted on."""
    big = pp.make_trapezoid('x', amplitude=0.6 * opts.max_grad, duration=2e-3, system=opts)
    inner = sc.LogicBlock('prephaser').add(0.0, big).add(0.0, big)
    with pytest.raises(sc.HardwareLimitError) as err:
        sc.compile(sc.LogicBlock('tr').add(0.0, pp.make_delay(1e-3)).add(1e-3, inner), opts)
    text = str(err.value)
    assert 'tr.prephaser' in text
    assert '1.000 ms' in text
    assert 'reached:' in text


def test_vector_norm_is_a_warning_not_an_error(opts) -> None:
    """
    Two axes ramping together reach sqrt(2) times the per-axis slew, which real amplifiers allow.

    Making this an error would reject the ordinary three-way winder overlap that the whole design
    exists to support -- so it is measured, said out loud, and not fatal.
    """
    strong = {'amplitude': 0.9 * opts.max_grad, 'duration': 1e-3, 'system': opts}

    def build():
        return compile_one(
            opts,
            (0.0, pp.make_trapezoid('x', **strong)),
            (0.0, pp.make_trapezoid('y', **strong)),
        )

    build()                      # a per-axis violation would have raised
    assert warned('vector-norm', build), 'the norm exceedance must not be silent either'


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
    assert len(out.block_events) == 1
    assert out.get_block(1).gx.type == 'trap'
    assert compiled_m0(out, 'x') == pytest.approx(2000.0, rel=1e-9)


def test_a_boundary_never_falls_inside_an_adc_window(opts) -> None:
    """Ramp sampling and vendor gridding both need the sampled gradient to be one event."""
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    long_g = pp.make_trapezoid('x', area=2000.0, duration=4e-3, system=opts)
    out = compile_one(opts, (0.0, long_g), (1.5e-3, adc))
    assert len(out.block_events) == 1
    assert out.get_block(1).gx.type == 'trap'


def test_a_split_preserves_area_and_continuity(opts) -> None:
    """When a boundary is forced, the two halves must join and their areas must add up."""
    long_g = pp.make_trapezoid('x', area=2000.0, duration=4e-3, system=opts)
    out = sc.compile(
        sc.LogicBlock('t').add(0.0, long_g).add(2e-3, sc.barrier('mid')), opts
    )
    assert len(out.block_events) == 2
    assert compiled_m0(out, 'x') == pytest.approx(2000.0, rel=1e-9)
    first, second = out.get_block(1).gx, out.get_block(2).gx
    assert float(first.last) == pytest.approx(float(second.first), rel=1e-9)


def test_a_split_mid_ramp_keeps_the_slew(opts) -> None:
    """The seam must not become a kink."""
    ramp = pp.make_trapezoid('x', amplitude=0.5 * opts.max_grad, duration=2e-3, system=opts)
    out = sc.compile(
        sc.LogicBlock('t').add(0.0, ramp).add(float(ramp.rise_time) / 2.0, sc.barrier()), opts
    )
    assert len(out.block_events) == 2


def test_barrier_forces_a_boundary_and_costs_no_time(opts) -> None:
    g = pp.make_trapezoid('x', area=500.0, duration=2e-3, system=opts)
    plain = compile_one(opts, (0.0, g))
    split = sc.compile(sc.LogicBlock('t').add(0.0, g).add(1e-3, sc.barrier()), opts)
    assert len(split.block_events) == len(plain.block_events) + 1
    assert split.duration()[0] == pytest.approx(plain.duration()[0])


def test_a_delay_only_block_compiles(opts) -> None:
    """The b=0 diffusion volume: correct duration, no events."""
    out = compile_one(opts, (0.0, pp.make_delay(4.2e-3)))
    assert out.duration()[0] == pytest.approx(4.2e-3)
    assert len(out.block_events) == 1


def test_an_over_long_interval_is_subdivided(opts) -> None:
    """pulseq stores a block duration in a fixed-width field."""
    out = compile_one(opts, (0.0, pp.make_delay(1.0)))
    assert out.duration()[0] == pytest.approx(1.0)


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
    sc.compile(tree, opts)
    after = [
        (start, path, id(event), content_hash(event))
        for start, event, path in sc.flatten(tree)
    ]

    assert after == before


# --------------------------------------------------------------------------------- definitions
def test_the_definitions_reach_the_sequence(opts) -> None:
    """
    Set during the compile, not at write time.

    That is what makes the returned ``pp.Sequence`` self-sufficient: nothing has to survive until
    ``write()`` to put them there, so there is no wrapper to carry them.
    """
    tree = sc.LogicBlock('gre').add(0.0, pp.make_delay(1e-3))
    out = sc.compile(tree, opts, definitions={'FOV': [0.25, 0.25, 0.005], 'TE': 8e-3})
    assert out.definitions['FOV'] == [0.25, 0.25, 0.005]
    assert out.definitions['Name'] == 'gre', 'the tag names the sequence by default'
    assert out.definitions['TotalDuration'] == pytest.approx(1e-3)


def test_two_sources_claiming_Name_is_an_error(opts) -> None:
    """
    Last-writer-wins is how a file came to say kSpaceCenterLine = 73 while its own navigator
    used 36.5: neither source was wrong about itself, and nothing compared them.
    """
    tree = sc.LogicBlock('gre').add(0.0, pp.make_delay(1e-3))
    with pytest.raises(sc.DefinitionConflict) as err:
        sc.compile(tree, opts, name='gre', definitions={'Name': 'something_else'})
    text = str(err.value)
    assert 'gre' in text and 'something_else' in text, 'the message must name both claimants'


def test_a_Name_that_agrees_is_not_a_conflict(opts) -> None:
    """Saying the same thing twice is not a disagreement, and stopping for it would be noise."""
    tree = sc.LogicBlock('gre').add(0.0, pp.make_delay(1e-3))
    out = sc.compile(tree, opts, name='gre', definitions={'Name': 'gre'})
    assert out.definitions['Name'] == 'gre'


# ---------------------------------------------------------------------------------- invariants
def test_compiled_duration_equals_the_tree_duration(opts) -> None:
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g).add(5e-3, g).add(20e-3, pp.make_delay(1e-3))
    out = sc.compile(tree, opts)
    # The duration invariant runs on every compile and raises, so returning at all is half the
    # assertion; the number itself is the other half.
    assert out.duration()[0] == pytest.approx(21e-3)


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
    assert compiled_m0(out, 'x') == pytest.approx(60.0, abs=1e-6)
    assert compiled_m0(out, 'y') == pytest.approx(-250.0, abs=1e-6)
    assert sc.moments(tree) == pytest.approx({'x': 60.0, 'y': -250.0}, abs=1e-6), (
        'and the tree-side measurement must agree, which is what the m0 invariant asserts'
    )


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
    for index, duration in out.block_durations.items():
        assert raster.holds(duration), f'block {index} is {duration * 1e6} us'


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
    for index in out.block_events:
        block = out.get_block(index)
        if getattr(block, 'rf', None) is not None:
            assert sc.Raster(opts.rf_raster_time).holds(float(block.rf.delay))


# ---------------------------------------------------------------------------------- provenance
def test_a_block_built_from_one_module_carries_that_modules_path(opts) -> None:
    """
    Provenance is no longer returned -- it is used where it is produced, to name the source in an
    error message -- so it is checked at the stage that computes it.
    """
    inner = sc.LogicBlock('spoiler').add(0.0, pp.make_trapezoid('z', area=500.0, system=opts))
    assert common_path([('tr', 'spoiler'), ('tr', 'spoiler')]) == ('tr', 'spoiler')

    # And end to end: an error about that gradient names where it came from.
    tree = sc.LogicBlock('tr').add(3e-6, inner)          # off the gradient raster
    with pytest.raises(sc.CompileError, match=r'tr\.spoiler'):
        sc.compile(tree, opts)


def test_the_origin_of_a_shared_block_is_the_common_ancestor(opts) -> None:
    """
    Three modules in one block have no single origin, and saying otherwise would mislead.

    The honest answer is the path they share; picking one arbitrarily would name a module that
    is only half responsible.
    """
    assert common_path([('tr', 'rephaser'), ('tr', 'blip')]) == ('tr',)
    assert common_path([('tr', 'ro'), ('other', 'ro')]) == ()
    assert common_path([(), ('tr', 'blip')]) == ('tr', 'blip'), 'an untagged node abstains'


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
    with pytest.raises(sc.CompileError) as err:
        sc.compile(tree, opts)
    text = str(err.value)
    assert 'repeat a k-space address' in text
    assert "'LIN': 7" in text, 'the message must name the address that repeats'


def test_unique_kspace_addresses_pass(opts) -> None:
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = sc.LogicBlock('t')
    for i in range(3):
        tree.add(i * 5e-3, adc, pp.make_label('LIN', 'SET', i))


def test_a_label_attaches_to_the_block_containing_its_start(opts) -> None:
    """A readout's LIN must land on the block holding the ADC, which is what the recon reads."""
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    out = sc.compile(
        sc.LogicBlock('t').add(0.0, g).add(1e-3, adc, pp.make_label('LIN', 'SET', 5)), opts
    )
    labels = out.evaluate_labels(evolution='adc')
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
        (0.0, pp.make_trapezoid('x', area=100.0, system=opts, **GENTLE)),
        (0.0, pp.make_trapezoid('x', area=200.0, system=opts, **GENTLE)),
    )
    grad = out.get_block(1).gx
    if grad.type == 'grad':
        assert len(grad.tt) <= 8, f'merged shape has {len(grad.tt)} points'


def test_a_lone_arbitrary_gradient_passes_through_untouched(opts) -> None:
    """
    A spiral must not be resampled: it takes the fast path when it sits alone in its interval.

    Resampling would be lossy exactly where the diffusion measurement is most sensitive.

    ``first`` and ``last`` are stated: without them ``make_arbitrary_grad`` extrapolates from the
    end samples, which puts the waveform at -2680 Hz/m half a raster before it starts.  pypulseq's
    own timing check calls that a step from zero, and the compile now stops on it.
    """
    wave = np.sin(np.linspace(0.0, np.pi, 500)) * 0.5 * opts.max_grad
    g = pp.make_arbitrary_grad(channel='x', waveform=wave, first=0.0, last=0.0, system=opts)
    out = compile_one(opts, (0.0, g))
    assert out.get_block(1).gx.type == 'grad'
    assert len(out.get_block(1).gx.waveform) == len(wave)
