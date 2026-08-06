"""
The waveform contract: what the compiler emits must equal what the tree said, exactly.

One oracle (:mod:`tests.compiler.fidelity`), many trees.  These are the cases where a
scheduler can plausibly go wrong: a boundary landing inside a waveform, two gradients summed,
a gradient sharing a block with an RF or an ADC, an arbitrary waveform that must not be
resampled, and the long trains where float error accumulates.

Marked ``xfail`` where the *current* compiler is known to lose fidelity; each such mark names
the defect from PLAN_COMPILER_V2.md and is removed by the work item that fixes it.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest
from fidelity import assert_matches, compare, compiled_knots, seam_discontinuities

import seqcraft as sc


# ------------------------------------------------------------------------------- one gradient
def test_single_trapezoid(system, opts) -> None:
    tree = sc.LogicBlock('t').add(0.0, pp.make_trapezoid('x', area=100.0, system=opts))
    assert_matches(tree, sc.compile(tree, system))


def test_triangle_no_flat_top(system, opts) -> None:
    g = pp.make_trapezoid('x', area=5.0, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g)
    assert_matches(tree, sc.compile(tree, system))


def test_three_axes_at_once(system, opts) -> None:
    gentle = {'area': 100.0, 'duration': 2e-3, 'rise_time': 200e-6, 'system': opts}
    tree = sc.LogicBlock('t')
    for ax in ('x', 'y', 'z'):
        tree.add(0.0, pp.make_trapezoid(ax, **gentle))
    assert_matches(tree, sc.compile(tree, system))


def test_sequential_gradients_same_axis(system, opts) -> None:
    """Disjoint lobes on one axis: nothing to sum, nothing to split."""
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g).add(1e-3, g).add(2e-3, g)
    assert_matches(tree, sc.compile(tree, system))


# ------------------------------------------------------------------------------------ merging
def test_two_trapezoids_summed(system, opts) -> None:
    tree = (
        sc.LogicBlock('t')
        .add(0.0, pp.make_trapezoid('x', area=100.0, duration=2e-3, system=opts))
        .add(0.0, pp.make_trapezoid('x', area=200.0, duration=2e-3, system=opts))
    )
    assert_matches(tree, sc.compile(tree, system))


def test_partially_overlapping_trapezoids(system, opts) -> None:
    """The hard merge: the sum has corners neither input has."""
    a = pp.make_trapezoid('x', area=100.0, duration=1e-3, rise_time=100e-6, system=opts)
    b = pp.make_trapezoid('x', area=-60.0, duration=1e-3, rise_time=100e-6, system=opts)
    tree = sc.LogicBlock('t').add(0.0, a).add(400e-6, b)
    assert_matches(tree, sc.compile(tree, system))


def test_three_overlapping_on_one_axis(system, opts) -> None:
    a = pp.make_trapezoid('x', area=80.0, duration=1e-3, rise_time=100e-6, system=opts)
    tree = sc.LogicBlock('t').add(0.0, a).add(200e-6, a).add(500e-6, a)
    assert_matches(tree, sc.compile(tree, system))


# ----------------------------------------------------------------------------------- splitting
def test_barrier_split_preserves_the_waveform(system, opts) -> None:
    g = pp.make_trapezoid('x', area=2000.0, duration=4e-3, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g).add(2e-3, sc.barrier('mid'))
    assert_matches(tree, sc.compile(tree, system))


def test_barrier_split_mid_ramp(system, opts) -> None:
    """The seam lands part-way up a ramp, so both halves start and end off zero."""
    g = pp.make_trapezoid('x', amplitude=0.5 * opts.max_grad, duration=2e-3, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g).add(float(g.rise_time) / 2.0, sc.barrier())
    assert_matches(tree, sc.compile(tree, system))


def test_many_barriers_across_one_gradient(system, opts) -> None:
    g = pp.make_trapezoid('x', area=2000.0, duration=4e-3, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g)
    for k in range(1, 8):
        tree.add(k * 500e-6, sc.barrier(f'b{k}'))
    assert_matches(tree, sc.compile(tree, system))


# ------------------------------------------------------------------- gradients beside rf / adc
def test_gradient_spanning_an_rf(system, opts) -> None:
    rf = pp.make_sinc_pulse(flip_angle=1.57, duration=1e-3, system=opts, use='excitation')
    g = pp.make_trapezoid('x', area=2000.0, duration=4e-3, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g).add(1.5e-3, rf)
    assert_matches(tree, sc.compile(tree, system))


def test_gradient_spanning_an_adc(system, opts) -> None:
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    g = pp.make_trapezoid('x', area=2000.0, duration=4e-3, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g).add(1.5e-3, adc)
    assert_matches(tree, sc.compile(tree, system))


def test_readout_with_prephaser(system, opts) -> None:
    ro = pp.make_trapezoid('x', amplitude=0.25 * opts.max_grad, duration=3.4e-3, system=opts)
    adc = pp.make_adc(num_samples=256, dwell=12e-6, delay=float(ro.rise_time), system=opts)
    pre = pp.make_trapezoid('x', area=-0.5 * float(ro.area), system=opts)
    tree = sc.LogicBlock('t').add(0.0, pre).add(float(pp.calc_duration(pre)), ro, adc)
    assert_matches(tree, sc.compile(tree, system))


# ------------------------------------------------------------------------ arbitrary waveforms
def _spiral_like(opts, n=200, scale=0.3, turns=2.0):
    """
    A smooth two-axis arbitrary waveform pair, zero at both ends and within the slew limit.

    `turns` is kept low deliberately: a fast oscillation on the gradient raster violates slew
    long before it stresses the compiler, and ``make_arbitrary_grad`` would refuse it.
    """
    s = np.linspace(0.0, 1.0, n)
    env = np.sin(np.pi * s)
    wx = env * np.cos(turns * np.pi * s) * scale * opts.max_grad
    wy = env * np.sin(turns * np.pi * s) * scale * opts.max_grad
    return wx, wy


def test_lone_arbitrary_gradient_untouched(system, opts) -> None:
    """The fast path: a spiral alone in its interval must survive bit-for-bit."""
    wx, _ = _spiral_like(opts)
    g = pp.make_arbitrary_grad('x', waveform=wx, first=0.0, last=0.0, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g)
    assert_matches(tree, sc.compile(tree, system))


def test_two_axis_arbitrary_pair(system, opts) -> None:
    wx, wy = _spiral_like(opts)
    gx = pp.make_arbitrary_grad('x', waveform=wx, first=0.0, last=0.0, system=opts)
    gy = pp.make_arbitrary_grad('y', waveform=wy, first=0.0, last=0.0, system=opts)
    tree = sc.LogicBlock('t').add(0.0, gx).add(0.0, gy)
    assert_matches(tree, sc.compile(tree, system))


def test_arbitrary_merged_with_a_trapezoid_is_reported_not_silent(system, opts) -> None:
    """
    The one case pulseq's two gradient representations cannot both be held.

    An arbitrary gradient is sampled at raster *centres*; a trapezoid's corners are on raster
    *edges*.  Their sum bends at both, and no pulseq gradient event has room for that -- so this
    is the one place the compiler cannot be exact.  What it must not do is be inexact quietly:
    the claim under test is that the resample is *reported*, with a bound, and that the bound is
    honest.
    """
    wx, _ = _spiral_like(opts, n=60)
    g = pp.make_arbitrary_grad('x', waveform=wx, first=0.0, last=0.0, system=opts)
    trap = pp.make_trapezoid('x', area=20.0, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g).add(0.0, trap)
    out = sc.compile(tree, system)

    reported = out.report.of_kind('grad_resample')
    assert reported, 'a resample that moves the waveform must never be silent'
    assert 'axis x' in reported[0].message

    # The reported bound must actually bound the measured deviation.
    claimed = float(reported[0].message.split('at most ')[1].split(' %')[0]) / 100.0
    worst = max(
        v['max_abs_error'] for v in compare(tree, out).values()
    ) / float(opts.max_grad)
    assert worst <= claimed + 1e-9, f'measured {worst:.4%} exceeds reported bound {claimed:.4%}'
    assert worst < 0.05, f'resample moved the waveform by {worst:.2%} of max_grad'

    # Area is what a resample preserves, so it is not evidence of fidelity -- but losing it would
    # mean something worse than a resample happened.
    assert not out.report.of_kind('moment')


def test_arbitrary_gradient_split_by_a_barrier(system, opts) -> None:
    wx, _ = _spiral_like(opts, n=60)
    g = pp.make_arbitrary_grad('x', waveform=wx, first=0.0, last=0.0, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g).add(300e-6, sc.barrier('mid'))
    assert_matches(tree, sc.compile(tree, system))


def test_extended_trapezoid_passes_through(system, opts) -> None:
    g = pp.make_extended_trapezoid(
        'x',
        times=np.array([0.0, 200e-6, 600e-6, 700e-6]),
        amplitudes=np.array([0.0, 1e4, 1e4, 0.0]),
        system=opts,
    )
    tree = sc.LogicBlock('t').add(0.0, g)
    assert_matches(tree, sc.compile(tree, system))


# -------------------------------------------------------------------------------- off raster
def test_gradient_started_off_the_gradient_raster_is_an_error(system, opts) -> None:
    """
    There is no correct snap, so the compiler must not pick one.

    This used to be silent: ``_in_block_delay`` rounded the start onto the raster, so a gradient
    asked for at 5 us played at 10 us.  m0 could not see it -- area does not depend on when a
    lobe plays -- and only m1 moved.  The tree asked for something the hardware cannot do, and
    which way to round is the caller's decision, not the compiler's.
    """
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    tree = sc.LogicBlock('t').add(5e-6, g)
    with pytest.raises(sc.CompileError, match='not a multiple of the 10 us gradient raster'):
        sc.compile(tree, system)


def test_the_off_raster_error_names_both_neighbouring_rasters(system, opts) -> None:
    """A message that only says "wrong" leaves the fix to be guessed."""
    g = pp.make_trapezoid('y', area=100.0, system=opts)
    tree = sc.LogicBlock('tr').add(0.0, sc.LogicBlock('spoiler').add(3e-6, g))
    with pytest.raises(sc.CompileError) as err:
        sc.compile(tree, system)
    text = str(err.value)
    assert 'tr.spoiler' in text
    assert '0.0 us' in text and '10.0 us' in text


def test_a_gradient_on_the_raster_is_unaffected(system, opts) -> None:
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    tree = sc.LogicBlock('t').add(20e-6, g)
    assert_matches(tree, sc.compile(tree, system))


# ---------------------------------------------------------------------------------- long runs
def test_epi_like_train(system, opts) -> None:
    """One readout gradient over many ADCs: every boundary lands inside the gradient."""
    adc = pp.make_adc(num_samples=64, dwell=4e-6, system=opts)
    n, period = 24, 500e-6
    g = pp.make_trapezoid('x', amplitude=0.2 * opts.max_grad,
                          duration=n * period, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g)
    for i in range(n):
        tree.add(100e-6 + i * period, adc)
    assert_matches(tree, sc.compile(tree, system))


def test_repeated_tr_accumulates_no_drift(system, opts) -> None:
    """40 TRs at a non-raster-friendly spacing: float drift would show as a shifted lobe."""
    rf = pp.make_sinc_pulse(flip_angle=1.57, duration=1e-3, system=opts, use='excitation')
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    inner = sc.LogicBlock('tr').add(0.0, rf).add(1.5e-3, g)
    tree = sc.LogicBlock('t')
    for i in range(40):
        tree.add(i * 5e-3, inner)
    assert_matches(tree, sc.compile(tree, system))


# ------------------------------------------------------------------------------- real modules
def _gre_tree(system, n_tr: int = 3):
    """A spoiled GRE, built the way tests/integration/conftest.py builds it."""
    exc = sc.modules.SincExcitation(
        system, flip_deg=15.0, duration_us=1000, slice_thickness_mm=5.0, rephase=False)
    ro = sc.modules.CartesianLine(
        system, fov_ro_mm=250.0, matrix_ro=32, readout_duration_us=3200, prephase=False)
    pe = sc.modules.PhaseEncode(system, fov_pe_mm=250.0, matrix_pe=32)
    spoil = sc.modules.Spoiler(system, twists=4, voxel_mm=5.0)

    raster = system.block_raster_s
    t_winders = exc.duration
    t_readout = sc.raster.ceil_to(exc.isodelay + 8e-3 - ro.time_to_echo, raster)
    t_spoil = sc.raster.ceil_to(t_readout + ro.duration, raster)

    tree = sc.LogicBlock('gre')
    for index in range(n_tr):
        t0 = index * 20e-3
        tree.add(t0, exc.build())
        tree.add(t0 + t_winders, exc.rephaser())
        tree.add(t0 + t_winders, pe.build(line=index - n_tr // 2))
        tree.add(t0 + t_winders, ro.prephaser_block())
        tree.add(t0 + t_readout, ro.build())
        tree.add(t0 + t_spoil, spoil.build())
    return tree


def test_gre_tr_from_modules(system) -> None:
    """
    A whole TR built the way a user would: three winders at one time on three axes.

    The case the design exists to support, and the one where a boundary chosen for the RF or the
    ADC has to leave three unrelated gradients intact.
    """
    tree = _gre_tree(system)
    assert_matches(tree, sc.compile(tree, system))


# ------------------------------------------------------------------------ seam continuity
def _seam_trees(system, opts):
    """Trees whose block seams are the interesting part."""
    long_g = pp.make_trapezoid('x', area=2000.0, duration=4e-3, system=opts)
    ramp = pp.make_trapezoid('x', amplitude=0.5 * opts.max_grad, duration=2e-3, system=opts)
    adc = pp.make_adc(num_samples=64, dwell=4e-6, system=opts)
    n, period = 12, 500e-6
    epi_g = pp.make_trapezoid('x', amplitude=0.2 * opts.max_grad,
                              duration=n * period, system=opts)
    epi = sc.LogicBlock('epi').add(0.0, epi_g)
    for i in range(n):
        epi.add(100e-6 + i * period, adc)
    return {
        'barrier_flat_top': sc.LogicBlock('t').add(0.0, long_g).add(2e-3, sc.barrier()),
        'barrier_mid_ramp': sc.LogicBlock('t').add(
            0.0, ramp).add(float(ramp.rise_time) / 2.0, sc.barrier()),
        'epi_train': epi,
        'gre': _gre_tree(system, n_tr=2),
    }


def test_no_amplitude_jumps_at_block_seams(system, opts) -> None:
    """
    A split must not become a kink: the amplitude leaving one block must equal the amplitude
    entering the next, and the sequence must start and end at zero.

    pypulseq enforces this when building, so this is a second opinion -- but it is the exact
    invariant a bad split breaks, and it is cheap.
    """
    for name, tree in _seam_trees(system, opts).items():
        bad = seam_discontinuities(sc.compile(tree, system))
        assert not bad, (
            f'{name}: amplitude jumps at block seams: '
            + '; '.join(f'{t * 1e6:.1f} us on {ax}: {a:.6g} -> {b:.6g}' for t, ax, a, b in bad)
        )


# --------------------------------------------------------------------------- the oracle itself
def test_the_oracle_catches_a_corrupted_waveform(system, opts) -> None:
    """
    An oracle that cannot fail proves nothing.

    Scale one compiled block's gradient by 1.01 and the comparison must notice.  This guards
    against the whole fidelity suite silently degrading into a no-op.
    """
    g = pp.make_trapezoid('x', area=2000.0, duration=4e-3, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g).add(2e-3, sc.barrier('mid'))
    out = sc.compile(tree, system)
    assert_matches(tree, out)                                   # clean to begin with

    block = out.seq.get_block(1)
    block.gx.waveform = np.asarray(block.gx.waveform) * 1.01    # 1 % too strong
    # get_block() rebuilds from the library, so patch the library entry the block came from.
    out.seq.get_block = lambda i, _b=block, _o=out.seq.get_block: (  # type: ignore[method-assign]
        _b if i == 1 else _o(i)
    )
    report = compare(tree, out)
    assert report['x']['max_abs_error'] > report['x']['tolerance'], (
        'the oracle did not notice a 1 % amplitude error'
    )
    with pytest.raises(AssertionError, match='differs from the logic block'):
        assert_matches(tree, out)


def test_the_oracle_reconstructs_a_continuous_waveform(system, opts) -> None:
    """
    Blocks are concatenated, not summed.

    Summing them double-counts every seam by the amplitude there, which is indistinguishable
    from a real waveform error -- so this pins the reconstruction the oracle depends on.
    """
    g = pp.make_trapezoid('x', amplitude=0.4 * opts.max_grad, duration=2e-3, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g).add(1e-3, sc.barrier())
    out = sc.compile(tree, system)
    times, amps = compiled_knots(out)['x']
    assert times == sorted(times), 'knot times must be strictly increasing'
    assert len(times) == len(set(times)), 'seam knots must be collapsed, not duplicated'
    assert abs(amps[0]) < 1e-9, 'must start at zero'
    assert abs(amps[-1]) < 1e-9, 'must end at zero'
    assert max(abs(a) for a in amps) == pytest.approx(float(g.amplitude), rel=1e-9), (
        'the reconstructed peak must be the trapezoid amplitude, not twice it'
    )


# ------------------------------------------------------------------------- the real modules
def _roundtrip_floor(event, opts) -> float:
    """
    Return pypulseq's own shape-library round-trip error for `event`, in Hz/m.

    A ``grad`` waveform is stored compressed -- run-length encoded derivatives -- and rebuilt by
    ``cumsum``, so what comes back out of ``get_block`` is not bit-identical to what went in.
    That is a floor no compiler can get under, and measuring it here rather than assuming a
    tolerance is what keeps the assertions below honest: they compare against pypulseq's floor,
    not against a number chosen to make them pass.
    """
    seq = pp.Sequence(system=opts)
    seq.add_block(event, pp.make_delay(float(pp.calc_duration(event))))
    before = np.asarray(event.waveform, dtype=float)
    after = np.asarray(seq.get_block(1).gx.waveform, dtype=float)
    return float(np.max(np.abs(before - after)))


def test_spiral_module_survives_compilation(system, opts) -> None:
    """
    The real SpiralVDS, not a synthetic stand-in: its readout must reach the .seq unmoved.

    The residual is pypulseq's shape compression, so the tolerance is *derived* from that floor
    rather than picked.  If seqcraft ever starts moving the waveform itself, this fails even
    though the absolute number is tiny.
    """
    ro = sc.modules.SpiralVDS(system, fov_mm=240, matrix=96, n_interleaves=8,
                              density=0.6, regime='default')
    tree = sc.LogicBlock('t').add(0.0, ro.build(interleaf=3))
    out = sc.compile(tree, system)

    floor = _roundtrip_floor(ro.gx, opts)
    report = compare(tree, out, atol=10.0 * max(floor, 1e-9), rtol=0.0)
    for ax, r in report.items():
        assert r['max_abs_error'] <= r['tolerance'], (
            f'axis {ax}: {r["max_abs_error"]:.4g} Hz/m exceeds 10x pypulseq\'s own '
            f'{floor:.4g} Hz/m round-trip floor -- seqcraft moved the waveform'
        )
    assert not out.report.of_kind('grad_resample'), 'a spiral alone must never be resampled'


def test_every_interleaf_of_a_spiral_is_faithful(system, opts) -> None:
    """Rotation is baked into the waveform, so each interleaf is a different shape to preserve."""
    ro = sc.modules.SpiralVDS(system, fov_mm=240, matrix=64, n_interleaves=4,
                              density=0.6, regime='default')
    floor = 10.0 * max(_roundtrip_floor(ro.gx, opts), 1e-9)
    for k in range(4):
        tree = sc.LogicBlock('t').add(0.0, ro.build(interleaf=k))
        out = sc.compile(tree, system)
        for ax, r in compare(tree, out, atol=floor, rtol=0.0).items():
            assert r['max_abs_error'] <= r['tolerance'], f'interleaf {k}, axis {ax}'


def test_a_diffusion_prepared_spiral_tr_is_faithful(system, opts) -> None:
    """
    The shape the DTI notebooks actually build: diffusion lobes, an excitation, a spiral readout.

    Traps and a raster-centre waveform on the same axes, in one tree -- the combination every
    individual test above isolates one part of.
    """
    rf = pp.make_sinc_pulse(flip_angle=1.57, duration=2e-3, system=opts, use='excitation',
                            slice_thickness=5e-3, apodization=0.5, time_bw_product=4)
    diff = sc.modules.MonopolarDiffusion(system, b_value_s_per_mm2=1000.0,
                                         refocus_duration_us=3000.0, regime='default')
    ro = sc.modules.SpiralVDS(system, fov_mm=240, matrix=96, n_interleaves=8,
                              density=0.6, regime='default')
    tree = sc.LogicBlock('tr')
    tree.add(0.0, rf)
    tree.add(3e-3, diff.build(part='pre', direction=(1.0, 1.0, 0.0)))
    tree.add(3e-3 + diff.duration + 1e-3, ro.build(interleaf=2))
    out = sc.compile(tree, system)

    floor = 10.0 * max(_roundtrip_floor(ro.gx, opts), 1e-9)
    for ax, r in compare(tree, out, atol=floor, rtol=0.0).items():
        assert r['max_abs_error'] <= r['tolerance'], (
            f'axis {ax}: {r["max_abs_error"]:.4g} Hz/m > {r["tolerance"]:.4g}'
        )
    assert not seam_discontinuities(out), 'a gradient must join across every block boundary'
