"""
:class:`Spiral2D`, checked against the oracles that share no code with it.

The order is the plan's, and it is the order the risk is in:

1. **the trajectory** against ``sc.kspace``, which is an independent integration of the compiled
   file.  Everything below leans on it, so it is first;
2. **the path and the traversal**, which import no pypulseq at all -- that is the structural claim
   of the module, and a test that reaches them without a compile is what asserts it;
3. **the ADC**, whose guard/divisor/limit arithmetic is the one hardware constraint a spiral meets
   that nothing else in this library does;
4. **the variants and the joins**, where ``v = 0`` at both ends and a join measured in context are
   the two things whose failure is silent for the two variants anyone builds first.
"""

from __future__ import annotations

import warnings

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.modules._support import oblique_trapezoid, sample_quantum, segment_samples
from seqcraft.modules.readout.spiral_2d import _design_path, _fov_curve, _traverse

FOV_MM = 220.0
MATRIX = 128
DWELL_S = 2.5e-6
DK_PER_M = 1e3 / FOV_MM
KMAX_PER_M = MATRIX / (2 * FOV_MM * 1e-3)

#: pypulseq's shape compression round-trips a thousand-sample waveform with a relative error of
#: 5e-8, which over a spiral's whole gradient is 1.3-4.7e-5 1/m of ``k`` -- and it says so, twice
#: per block.  This is that floor, not a design tolerance: tightening it is a day spent on
#: something no module can improve.
KSPACE_TOLERANCE_PER_M = 1e-4


@pytest.fixture(scope='module')
def opts() -> pp.Opts:
    """The reference protocol every number here is measured against."""
    return pp.Opts(
        max_grad=40, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s', B0=3.0,
        rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
        adc_samples_limit=8192,
    )


def spiral(opts: pp.Opts, **kwargs) -> sc.modules.Spiral2D:
    """One readout at the reference protocol, with whatever the test is varying.

    A keyword given as ``None`` is *removed* rather than passed on, so a refusal test can ask
    for the argument to be absent rather than empty -- which for `dwell_s` are two different
    questions and the module answers them the same way.
    """
    base = {'opts': opts, 'fov_mm': FOV_MM, 'matrix': MATRIX, 'shots': 4, 'dwell_s': DWELL_S}
    merged = {**base, **kwargs}
    return sc.modules.Spiral2D(**{k: v for k, v in merged.items()
                                  if v is not None or k not in base})


def compiled(tree: sc.LogicBlock, opts: pp.Opts) -> tuple[object, list[str]]:
    """Compile, keeping the warnings that are not the vector-norm one every module produces."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', sc.SeqCraftWarning)
        seq = sc.compile(tree, opts)
    return seq, [str(w.message) for w in caught
                 if issubclass(w.category, sc.SeqCraftWarning)
                 and 'vector-norm' not in str(w.message)]


# ------------------------------------------------------------------ 1. the trajectory that plays
@pytest.mark.parametrize('shot', range(4))
def test_k_per_m_agrees_with_kspace_on_every_interleaf(opts: pp.Opts, shot: int) -> None:
    """The module's array is where the samples are, and ``sc.kspace`` is the independent oracle.

    Per interleaf rather than once, so a rotation applied twice or not at all cannot pass.
    """
    readout = spiral(opts)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        measured = sc.kspace(sc.LogicBlock('ro').add(0.0, readout(shot=shot)), opts)['k_adc']
    assert np.abs(measured[:2] - readout.k_per_m(shot=shot)).max() < KSPACE_TOLERANCE_PER_M


@pytest.mark.parametrize(('variant', 'echoes'), [
    ('out', 1), ('out', 3), ('in', 1), ('in', 2), ('in-out', 1), ('in-out', 2), ('out-in', 2),
])
def test_every_variant_agrees_with_kspace_and_stays_an_arbitrary_gradient(
    opts: pp.Opts, variant: str, echoes: int,
) -> None:
    """A split spiral stays a spiral: every compiled block's gradient is still a ``grad``."""
    readout = spiral(opts, variant=variant, echoes=echoes)
    tree = sc.LogicBlock('ro').add(0.0, readout(shot=1))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        measured = sc.kspace(tree, opts)['k_adc']
    seq, notes = compiled(tree, opts)
    kinds = {getattr(getattr(seq.get_block(i), axis, None), 'type', None)
             for i in sorted(seq.block_events) for axis in ('gx', 'gy')} - {None}
    assert np.abs(measured[:2] - readout.k_per_m(shot=1)).max() < KSPACE_TOLERANCE_PER_M
    assert kinds == {'grad'}
    assert not notes, notes


def test_the_reported_kmax_is_the_built_one_and_not_the_designed_one(opts: pp.Opts) -> None:
    """The design pass's ``k`` is discarded, and this is what says so without reading the source.

    A trajectory read off the design would reach ``matrix / (2 fov)`` **exactly**, because that is
    the number the path was integrated to.  Rasterising it loses 0.03 % of that, which is a 0.03 %
    scaling of the resolution and nothing else -- and it is the fingerprint of a measurement.
    """
    readout = spiral(opts)
    built = float(np.hypot(*readout.k_per_m(shot=0)).max())
    assert built != readout.kmax_per_m
    assert abs(built / readout.kmax_per_m - 1.0) < 5e-3


# ------------------------------------------------- 2. the path and the traversal, without pypulseq
def test_the_path_satisfies_the_nyquist_condition_and_nothing_else() -> None:
    """Turn spacing is ``shots / FOV``, measured off the curve rather than read off the formula."""
    path, _ = _design_path(fov_m=FOV_MM / 1e3, matrix=MATRIX, shots=1, curve=_fov_curve((1.0,)))
    angle = np.unwrap(np.angle(path[1:]))
    radius = np.abs(path[1:])
    spacing = radius - np.interp(angle - 2 * np.pi, angle, radius, left=np.nan)
    inside = (radius > 20.0) & (radius < 0.9 * KMAX_PER_M)
    assert abs(float(np.nanmedian(spacing[inside])) / DK_PER_M - 1.0) < 2e-3
    assert abs(angle[-1] / (2 * np.pi) - MATRIX / 2) < 1e-2


@pytest.mark.parametrize('density', [(1.0,), (1.0, 0.5), (1.0, 1.0, 0.7, 0.33)])
def test_arc_length_matches_the_closed_form_and_scales_with_shots(density: tuple) -> None:
    """``integral 2 pi FOV(kr) kr dkr``, and ``shots`` divides one curve rather than shortening it."""
    curve = _fov_curve(density)
    u = np.linspace(0.0, 1.0, 200_001)
    predicted = float(np.trapezoid(
        2 * np.pi * (FOV_MM / 1e3) * curve(u) * u * KMAX_PER_M, u * KMAX_PER_M))
    lengths = [float(_design_path(fov_m=FOV_MM / 1e3, matrix=MATRIX, shots=n, curve=curve)[1][-1])
               for n in (1, 2, 4)]
    assert abs(lengths[0] / predicted - 1.0) < 5e-3
    for n, length in zip((1, 2, 4), lengths):
        assert abs(n * length / lengths[0] - 1.0) < 2e-3


def test_the_traversal_starts_and_ends_at_rest_and_beats_its_own_lower_bound() -> None:
    """``v = 0`` at both ends is the design decision, and the closed form is a *bound*."""
    max_grad, max_slew = 1.70304e6, 7.66368e9
    curve = _fov_curve((1.0,))
    path, arclength = _design_path(fov_m=FOV_MM / 1e3, matrix=MATRIX, shots=1, curve=curve)
    _, speed, times, _ = _traverse(path, arclength, max_grad=max_grad, max_slew=max_slew)
    u = np.linspace(0.0, 1.0, 200_001)
    bound = float(2 * np.pi / np.sqrt(max_slew) * np.trapezoid(
        (FOV_MM / 1e3) * curve(u) * np.sqrt(u * KMAX_PER_M), u * KMAX_PER_M))
    assert speed[0] == 0.0
    assert speed[-1] == 0.0
    assert bound <= float(times[-1]) <= 1.05 * bound
    assert float(speed.max()) <= max_grad * (1 + 1e-9)


def test_density_is_sampled_values_and_not_polynomial_coefficients(opts: pp.Opts) -> None:
    """``(1.0, 0.5)`` halves the FOV at the edge; under the other convention it would raise it."""
    tapered = spiral(opts, density=(1.0, 0.5))
    uniform = spiral(opts)
    assert tapered.edge_undersampling == pytest.approx(2.0)
    assert float(tapered.fov_mm_of(1.0)) == pytest.approx(FOV_MM / 2)
    assert tapered.arm_duration_s < uniform.arm_duration_s
    assert abs(tapered.worst_turn_gap_dk - 2.0) < 0.05
    assert tapered.kmax_per_m == uniform.kmax_per_m       # a taper does not move the resolution


def test_a_callable_density_is_data_too(opts: pp.Opts) -> None:
    """Nyquist inside a radius and a smooth taper past it, which a sequence cannot say."""
    def dual(u, u0=0.4, r_edge=3.0):
        t = np.clip((u - u0) / (1 - u0), 0.0, 1.0)
        return 1.0 / (1.0 + (r_edge - 1.0) * t ** 2)

    readout = spiral(opts, density=dual)
    assert readout.edge_undersampling == pytest.approx(3.0)
    assert float(readout.fov_mm_of(0.4)) == pytest.approx(FOV_MM)
    assert readout.peak_slew_hz_m_s <= float(opts.max_slew) * (1 + 1e-9)


# -------------------------------------------------------- 3. the ADC: guard, divisor, limit, split
def test_the_sample_split_meets_the_limit_and_the_divisor(opts: pp.Opts) -> None:
    """Measured: a single-shot arm needs three events, and four shots need one."""
    counts = {n: spiral(opts, shots=n).segment_samples for n in (1, 2, 4, 8)}
    assert [len(c) for c in counts.values()] == [3, 2, 1, 1]
    assert counts[1] == (8192, 8192, 4580)
    divisor = int(getattr(opts, 'adc_samples_divisor', 1) or 1)
    for split in counts.values():
        assert all(count % divisor == 0 for count in split)
        assert all(count <= int(opts.adc_samples_limit) for count in split)


def test_a_split_readout_passes_pypulseqs_own_timing_check(opts: pp.Opts) -> None:
    """`adc_samples_divisor` is pypulseq's, and it is what a scanner refuses a block over."""
    seq, notes = compiled(sc.LogicBlock('ro').add(0.0, spiral(opts, shots=1)(shot=0)), opts)
    ok, error = seq.check_timing()[:2]
    assert ok, error
    assert not notes, notes


def test_adc_segments_of_one_refuses_and_names_a_shots_that_fits(opts: pp.Opts) -> None:
    """A split costs a hole in the trajectory, so it is a recovery and `shots` is the answer."""
    with pytest.raises(sc.ConfigurationError, match='adc_segments=1'):
        spiral(opts, shots=1, adc_segments=1)
    try:
        spiral(opts, shots=1, adc_segments=1)
    except sc.ConfigurationError as error:
        named = int(str(error).split('shots=')[1].split(' ')[0])
    assert spiral(opts, shots=named, adc_segments=1).adc_segments == 1


def test_an_acquisition_region_is_not_an_arm(opts: pp.Opts) -> None:
    """One ADC spans an in-out pair, which is what puts ``k = 0`` inside a sampling window.

    Built per arm instead, the sample nearest the seam would be a whole guard away -- and the
    contrast is the assertion: a guard is four dwells here, and the measured offset is half of one.
    """
    readout = spiral(opts, variant='in-out', prephase=False)
    offset_s = abs(readout.time_to_echo(0) - readout.arm_duration_s - readout.guard_s)
    k = readout.k_per_m(shot=0)
    assert offset_s <= readout.dwell_s
    assert np.hypot(*k[:, readout.echo_sample(0)]) < 0.05 * DK_PER_M


def test_the_dwell_warns_rather_than_refusing_and_names_a_legal_one(opts: pp.Opts) -> None:
    """Deliberate along-arm undersampling is a real choice; not knowing about it is not."""
    with pytest.warns(sc.SeqCraftWarning, match='undersampled along its own path'):
        coarse = spiral(opts, dwell_s=6e-6)
    assert coarse.worst_sample_gap_dk > 1.0
    assert spiral(opts).worst_sample_gap_dk < 1.0


# --------------------------------------------------------- 4. the variants, the echoes, the joins
@pytest.mark.parametrize('variant', ['out', 'in', 'in-out', 'out-in'])
@pytest.mark.parametrize('echoes', [1, 2, 3])
def test_every_join_is_continuous_and_nothing_is_resampled(
    opts: pp.Opts, variant: str, echoes: int,
) -> None:
    """``v = 0`` at both ends of the arm, so every seam meets at ``g = 0`` with nothing between.

    Silent for ``'out'`` and ``'in'``, which build either way -- which is why the two-arm variants
    are in the same test rather than a later one.
    """
    if variant == 'out-in' and echoes == 1:
        with pytest.raises(sc.ConfigurationError, match='at least two echoes'):
            spiral(opts, variant=variant, echoes=echoes)
        return
    readout = spiral(opts, variant=variant, echoes=echoes)
    _, notes = compiled(sc.LogicBlock('ro').add(0.0, readout(shot=2)), opts)
    assert not notes, notes
    assert readout.peak_slew_hz_m_s <= float(opts.max_slew) * (1 + 1e-6)
    assert readout.peak_grad_hz_m <= float(opts.max_grad) * (1 + 1e-6)
    assert len(readout.te_s) == echoes


@pytest.mark.parametrize('shots', [1, 2, 4, 8])
@pytest.mark.parametrize('variant', ['out', 'in', 'in-out', 'out-in'])
def test_k_returns_to_the_origin_at_the_end_of_every_block(
    opts: pp.Opts, shots: int, variant: str,
) -> None:
    """A residual here is a phase ramp across the *next* excitation's image.

    Uncorrected, the same waveform leaves 0.108 1/m -- 0.024 ``dk``, invisible in every plot.
    """
    echoes = 2 if variant == 'out-in' else 1
    readout = spiral(opts, shots=shots, variant=variant, echoes=echoes,
                     prephase=None if variant in ('out', 'out-in') else True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        k = sc.kspace(sc.LogicBlock('ro').add(0.0, readout(shot=0)), opts)['k']
    assert np.abs(k[:2, -1]).max() < 1e-3
    assert np.abs(readout.k_end_per_m(0)).max() < 1e-3


def test_out_in_retraces_rather_than_covering_twice_as_much(opts: pp.Opts) -> None:
    """The return arm adds no k-space coverage; its value is the **second echo**.

    Written as a test so the \"out-in covers twice as much\" reading dies here rather than in
    a notebook.  Reflected in *time* rather than in sample index, because the ADC seam puts a
    guard in the middle of the acquisition and the indices are no longer symmetric across it.
    """
    readout = spiral(opts, variant='out-in', echoes=2)
    k, t = readout.k_per_m(shot=1), readout.t_adc_s
    mirrored = t[0] + t[-1] - t
    apart = np.hypot(np.interp(mirrored, t, k[0]) - k[0],
                     np.interp(mirrored, t, k[1]) - k[1])
    assert float(np.median(apart)) < 1e-6
    assert float(apart.max()) < 0.2 * DK_PER_M     # the worst is the seam, interpolated across
    assert len(readout.te_s) == 2


def test_te_s_is_read_off_the_measured_trajectory(opts: pp.Opts) -> None:
    """``k = 0`` is a sample, not an instant, so `echo_sample` is ``argmin |k|`` and not a solve."""
    for variant, echoes in (('out', 3), ('in', 2), ('in-out', 2), ('out-in', 3)):
        readout = spiral(opts, variant=variant, echoes=echoes)
        k = readout.k_per_m(shot=0)
        assert len(readout.te_s) == echoes
        for echo in range(echoes):
            sample = readout.echo_sample(echo)
            assert readout.te_s[echo] == readout.t_adc_s[sample]
            assert np.hypot(*k[:, sample]) < 0.05 * DK_PER_M
        assert readout.echo_spacing_s == pytest.approx(readout.te_s[1] - readout.te_s[0])


def test_the_labels_are_read_off_the_compiled_sequence(opts: pp.Opts) -> None:
    """`LIN` for the interleaf, `ECO` above one echo, `SEG`/`REF`/`IMA` on request -- and no `REV`."""
    readout = spiral(opts, variant='out', echoes=2)
    tree = sc.LogicBlock('ro').add(0.0, readout(shot=2, segment=5, reference=True))
    seq, _ = compiled(tree, opts)
    labels = seq.evaluate_labels(evolution='adc')
    assert list(np.atleast_1d(labels['LIN'])) == [2, 2]
    assert list(np.atleast_1d(labels['ECO'])) == [0, 1]
    assert list(np.atleast_1d(labels['SEG'])) == [5, 5]
    assert list(np.atleast_1d(labels['REF'])) == [1, 1]
    assert list(np.atleast_1d(labels['IMA'])) == [0, 0]
    assert 'REV' not in labels

    plain, _ = compiled(sc.LogicBlock('ro').add(0.0, spiral(opts)(shot=0)), opts)
    assert 'ECO' not in plain.evaluate_labels(evolution='adc')


def test_a_dummy_shot_plays_the_same_gradients_with_no_adc(opts: pp.Opts) -> None:
    """A dummy has to load the amplifier exactly as a real shot does, or it establishes nothing."""
    readout = spiral(opts)
    real, dummy = readout(shot=1), readout(shot=1, acquire=False)
    assert real.duration == dummy.duration
    assert sc.moments(real, 0) == pytest.approx(sc.moments(dummy, 0))
    seq, _ = compiled(sc.LogicBlock('ro').add(0.0, dummy), opts)
    assert not any(getattr(seq.get_block(i), 'adc', None) for i in sorted(seq.block_events))


# ------------------------------------------------------------------------------- 5. the refusals
@pytest.mark.parametrize(('kwargs', 'match'), [
    ({'dwell_s': None}, 'exactly one of dwell_s and bandwidth_hz'),
    ({'dwell_s': DWELL_S, 'bandwidth_hz': 4e5}, 'exactly one of dwell_s and bandwidth_hz'),
    ({'dwell_s': DWELL_S, 'variant': 'inward'}, 'variant must be one of'),
    ({'dwell_s': DWELL_S, 'variant': 'out-in', 'echoes': 1}, 'at least two echoes'),
    ({'dwell_s': DWELL_S, 'prephase': True}, 'nothing for prephase=True to do'),
    ({'dwell_s': DWELL_S, 'variant': 'in', 'rewind': True}, 'nothing for rewind=True to do'),
    ({'dwell_s': DWELL_S, 'density': (1.0, 0.0)}, 'zero or negative'),
    ({'dwell_s': DWELL_S, 'density': (1.0, -0.5)}, 'zero or negative'),
    ({'dwell_s': DWELL_S, 'axes': ('x', 'x')}, 'two different channels'),
    ({'dwell_s': DWELL_S, 'echo_spacing_s': 20e-3}, 'at least two echoes to be the spacing'),
    ({'dwell_s': DWELL_S, 'echoes': 2, 'echo_spacing_s': 1e-3}, 'shorter than the'),
])
def test_what_cannot_take_effect_raises(opts: pp.Opts, kwargs: dict, match: str) -> None:
    """Every refusal names what to change, which is the difference from a silent default."""
    with pytest.raises(sc.ConfigurationError, match=match):
        spiral(opts, **kwargs)


def test_build_requires_exactly_one_of_shot_and_angle_rad(opts: pp.Opts) -> None:
    """A four-shot acquisition that silently acquires one interleaf four times is a plausible bug."""
    readout = spiral(opts)
    with pytest.raises(sc.ConfigurationError, match='exactly one of shot and angle_rad'):
        readout()
    with pytest.raises(sc.ConfigurationError, match='exactly one of shot and angle_rad'):
        readout(shot=0, angle_rad=0.5)
    with pytest.raises(sc.ConfigurationError, match=r'shot must be in range\(4\)'):
        readout(shot=4)
    assert readout(angle_rad=0.5).duration == readout(shot=0).duration


def test_echo_spacing_lengthens_the_joins_to_the_number_asked_for(opts: pp.Opts) -> None:
    """A `g = 0` join takes dead time for nothing but time, which is what a chosen dTE costs."""
    minimum = spiral(opts, echoes=3)
    stretched = spiral(opts, echoes=3, echo_spacing_s=20e-3)
    assert stretched.echo_spacing_s == pytest.approx(20e-3, abs=float(opts.grad_raster_time))
    assert stretched.echo_spacing_s > minimum.echo_spacing_s
    with pytest.raises(sc.ConfigurationError, match='gap between two echoes'):
        _ = spiral(opts).echo_spacing_s


# ----------------------------------------------------------------------------- 6. the helpers
def test_oblique_trapezoid_holds_the_vector_limit_where_a_per_axis_pair_does_not(
    opts: pp.Opts,
) -> None:
    """Two trapezoids stretched to a common duration ask for ``sqrt(2)`` times the vector slew."""
    shapes = set()
    for angle in (0.0, np.pi / 4, 1.0):
        area = KMAX_PER_M * np.array([np.cos(angle), np.sin(angle)])
        pair = oblique_trapezoid(area_per_m=area, opts=opts)
        peak = float(np.hypot(*[float(g.amplitude) for g in pair]))
        shapes.add((round(float(pp.calc_duration(*pair)) * 1e9), round(peak)))
        assert peak / float(pair[0].rise_time) <= float(opts.max_slew) * (1 + 1e-9)
        assert np.hypot(float(pair[0].area), float(pair[1].area)) == pytest.approx(
            float(np.hypot(*area)))
    assert len(shapes) == 1, 'the timing and the vector peak must not depend on the angle'

    naive = [pp.make_trapezoid(axis, area=KMAX_PER_M / np.sqrt(2), duration=340e-6, system=opts)
             for axis in 'xy']
    naive_slew = np.hypot(*[float(g.amplitude) for g in naive]) / float(naive[0].rise_time)
    assert naive_slew > 1.3 * float(opts.max_slew)


def test_sample_quantum_and_the_split_land_every_seam_on_the_gradient_raster(
    opts: pp.Opts,
) -> None:
    """Solved in integer ticks rather than searched, so the boundary is legal by construction."""
    raster_ticks = sc.timing.to_ticks(float(opts.grad_raster_time))
    for dwell_s in (2.5e-6, 2e-6, 3e-6, 5e-6):
        quantum = sample_quantum(dwell_s=dwell_s, opts=opts)
        assert quantum % int(getattr(opts, 'adc_samples_divisor', 1) or 1) == 0
        assert (quantum * sc.timing.to_ticks(dwell_s)) % raster_ticks == 0
        split = segment_samples(20_000, dwell_s=dwell_s, limit=8192, opts=opts)
        assert all(count % quantum == 0 for count in split)
        assert all(count <= 8192 for count in split)
        assert split[-1] <= split[0], 'the short segment goes last'
    assert segment_samples(20_000, dwell_s=2.5e-6, limit=0, opts=opts) == (20_000,)
