"""
EPI readout physics.

The contract assertions come from :mod:`seqcraft.testing`; what is here is the numbers an
independent calculation gives -- the sampled k-space extent against the lobe's own area, where the
echo lands, and that the blip never touches an ADC window.  Every one of them was a bug first.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import seqcraft as sc
from seqcraft.core.errors import ConfigurationError

BIG = {'fov_ro_mm': 240.0, 'matrix_ro': 128, 'fov_pe_mm': 240.0, 'matrix_pe': 128}
SMALL = {'fov_ro_mm': 240.0, 'matrix_ro': 64, 'fov_pe_mm': 240.0, 'matrix_pe': 64}


@pytest.fixture
def system() -> sc.System:
    """A Cima.X with an EPI regime at full slew, which is where the echo spacing wants to be."""
    return sc.System.preset('cima_x').derate('epi', grad=0.85, slew=1.0)


def readout(system: sc.System, **kwargs: object) -> sc.modules.EPIReadout:
    """An EPI readout on the small matrix unless the test says otherwise."""
    merged = {**SMALL, 'regime': 'epi'}
    merged.update(kwargs)
    return sc.modules.EPIReadout(system, **merged)  # type: ignore[arg-type]


def spin_echo(system: sc.System, ro: sc.modules.EPIReadout) -> sc.CompiledSequence:
    """
    Place `ro` inside a spin echo, which is the only way ``kspace()`` can be asked where k is.

    ``calculate_kspacePP`` needs the refocusing pulse to apply its sign flip, and a readout on its
    own has no echo to be referenced to.
    """
    exc = sc.modules.SincExcitation(
        system, flip_deg=90, duration_us=2000, slice_thickness_mm=4, rephase=True, regime='epi')
    refoc = sc.modules.SincRefocusing(
        system, flip_deg=180, duration_us=4000, slice_thickness_mm=4,
        crusher_twists=4, crusher_voxel_mm=4, regime='epi')
    raster = system.block_raster
    first_half = raster.ceil((exc.duration - exc.isodelay) + refoc.isodelay)
    second_half = raster.ceil((refoc.duration - refoc.isodelay) + ro.time_to_echo)
    te = 2 * max(first_half, second_half)
    t_refoc = raster.ceil(exc.isodelay + te / 2 - refoc.isodelay)
    t_read = raster.ceil(exc.isodelay + te - ro.time_to_echo)

    tree = sc.LogicBlock('se_epi')
    tree.add(0.0, exc.build())
    tree.add(t_refoc, refoc.build())
    tree.add(t_read, ro.build())
    return sc.compile(tree, system, regime='epi', name='se_epi')


# ------------------------------------------------------------------------------- the contract
@pytest.mark.parametrize(
    ('kwargs', 'build_args'),
    [
        ({}, {}),
        ({'partial_fourier_pe': 0.75}, {}),
        ({'partial_echo': 0.75}, {}),
        ({'n_shots': 2}, {'shot': 1}),
        ({'n_shots': 4, 'partial_fourier_pe': 0.75}, {'shot': 3}),
        ({}, {'pe_polarity': -1}),
        ({'labels': False}, {}),
        ({'flat_time_us': 200.0}, {}),
    ],
    ids=['plain', 'pf', 'partial_echo', 'two_shot', 'four_shot', 'blip_down', 'no_labels', 'flat'],
)
def test_contract(system, kwargs, build_args) -> None:
    """Purity, determinism, an honest duration, per-axis limits, raster, and a clean solo compile."""
    sc.testing.assert_all(readout(system, **kwargs), **build_args)


def test_the_train_compiles_without_merge_or_raster_warnings(system) -> None:
    """
    One gradient event per axis, an ADC node on the block raster: nothing to report.

    Built per-echo instead, the tail of echo *n* and the head of *n+1* share a block and the
    compiler warns once per echo -- 96 warnings on the example sequence, which is the state in
    which a warning stops meaning anything.  And an ADC node off the block raster is *snapped*,
    which moves sampling against the gradient by up to half a raster: 10 1/m of kx here.
    """
    out = spin_echo(system, readout(system, **BIG, partial_fourier_pe=0.75))
    issues = out.check().issues
    assert not [i for i in issues if i.kind == 'grad_merge'], 'same-axis merges in the train'
    assert not [i for i in issues if i.kind == 'raster'], 'an ADC or RF was snapped'
    assert out.check().ok


def test_sampled_extent_is_k_max_not_the_lobe_area(system) -> None:
    """
    What has to cover 2 k_max is the **sampled** area, not the lobe's own.

    The ADC skips the first and last ``blip / 2``, so the lobe must overshoot by the two unsampled
    corners.  On the example sequence that is 540.5 against 533.3 1/m -- 1.7 dk that a prephaser of
    minus k_max would displace every shot by, and a k-space offset is a linear phase ramp across the
    image: the magnitude image looks perfect and every phase-derived quantity is wrong.
    """
    ro = readout(system, **BIG)
    lobe_area = ro.peak_amplitude_Hz_per_m * (ro.ramp_time + ro.flat_time)
    k = ro.trajectory(0)
    sampled = float(k[0].max() - k[0].min())
    wanted = (1.0 + ro.partial_echo) * ro.k_max_ro_per_m

    assert lobe_area > wanted * 1.005, 'the lobe should overshoot the sampled extent'
    # Samples sit half a dwell inside the window edges, so the measured span falls just short of
    # the designed one -- by exactly that half dwell of gradient, and no more.
    slack = ro.peak_amplitude_Hz_per_m * ro.dwell_s
    assert wanted - slack <= sampled <= wanted, (
        f'sampled extent {sampled:.3f} 1/m against a designed {wanted:.3f}; '
        f'the lobe carries {lobe_area:.3f}'
    )


@pytest.mark.parametrize('partial_echo', [1.0, 0.75, 0.6])
def test_the_prephaser_places_kx_exactly(system, partial_echo) -> None:
    """
    kx must enter the first echo at ``-partial_echo * k_max`` and leave it at ``+k_max``.

    Checked against pypulseq's own trajectory rather than the module's, so a prephaser and a
    ``trajectory()`` that agree with each other but not with the waveform cannot pass.
    """
    ro = readout(system, partial_echo=partial_echo)
    out = spin_echo(system, ro)
    k_adc = out.kspace()['k_adc']
    first_echo = k_adc[0][: ro.samples_per_echo]
    slack = ro.peak_amplitude_Hz_per_m * ro.dwell_s

    assert first_echo[0] == pytest.approx(-partial_echo * ro.k_max_ro_per_m, abs=slack)
    assert first_echo[-1] == pytest.approx(ro.k_max_ro_per_m, abs=slack)


def test_ky_is_constant_while_a_line_is_read(system) -> None:
    """
    The blip is centred on the junction, so it lies wholly outside every ADC window.

    If it ever leaked in, ky would drift within a line and the trajectory would stop being a set of
    horizontal lines -- which the reconstruction can handle, but the phase-encode bandwidth and the
    distortion figure would both become fiction.
    """
    ro = readout(system, **BIG, partial_fourier_pe=0.75)
    out = spin_echo(system, ro)
    per_echo = out.kspace()['k_adc'][1][: ro.n_samples].reshape(ro.n_echoes, ro.samples_per_echo)
    assert float(np.ptp(per_echo, axis=1).max()) == 0.0

    # And the sampling window really is inside the junctions, with the ADC dead time fitting.
    assert ro.sample_offset >= ro.blip_duration / 2.0
    tail = ro.echo_spacing - (ro.sample_offset + ro.samples_per_echo * ro.dwell_s)
    assert tail >= float(ro.opts.adc_dead_time)


def test_k_is_zero_at_the_echo_on_every_axis(system) -> None:
    """
    The assertion that catches a missing slice rewinder, and a mis-timed readout.

    ``time_to_echo`` is what TE is measured from, so if it points anywhere but k = 0 the diffusion
    weighting is referenced to the wrong instant.
    """
    ro = readout(system, **BIG, partial_fourier_pe=0.75)
    out = spin_echo(system, ro)
    k = out.kspace()
    at_echo = int(np.argmin(np.abs(k['t_adc'] - (k['t_excitation'][0] + _te_of(out)))))
    slack = ro.peak_amplitude_Hz_per_m * ro.dwell_s
    for axis, value in zip('xyz', k['k_adc'][:, at_echo]):
        assert abs(float(value)) < slack, f'k_{axis} = {float(value):+.2f} 1/m at the echo'


def _te_of(out: sc.CompiledSequence) -> float:
    """TE as the compiled sequence measures it: twice excitation-to-refocusing."""
    k = out.kspace()
    return 2.0 * float(k['t_refocusing'][0] - k['t_excitation'][0])


@pytest.mark.parametrize('partial_fourier_pe', [1.0, 0.75, 0.6])
@pytest.mark.parametrize('n_shots', [1, 2, 4])
def test_time_to_echo_is_shot_independent_and_points_at_ky_zero(
    system, partial_fourier_pe, n_shots
) -> None:
    """
    Every shot has identical timing, and the echo it points at is the ky = 0 one.

    Shot-independence is required twice over: a build argument may not change a timing property,
    and segmented EPI plays every shot at the same TE so they carry the same T2 weighting.  Only
    one shot actually samples the centre line; the others pass within ``(n_shots - 1) * dk`` of it,
    which is what ``p // n_shots`` delivers.
    """
    ro = readout(system, partial_fourier_pe=partial_fourier_pe, n_shots=n_shots)
    durations = {ro.build(shot=s).duration for s in range(n_shots)}
    assert len(durations) == 1, 'shots differ in duration'

    holders = [s for s in range(n_shots) if ro.centre_line in ro.lines(s)]
    assert len(holders) == 1
    assert ro.lines(holders[0])[ro.echo_index] == ro.centre_line
    for shot in range(n_shots):
        offset = abs(ro.signed_lines(shot)[ro.echo_index]) * ro.dk_pe_per_m
        assert offset <= (n_shots - 1) * ro.dk_pe_per_m + 1e-12


def test_partial_fourier_brings_the_echo_forward(system) -> None:
    """
    The entire reason to use it.  Nothing else about the readout moves.

    Taking the train's midpoint instead would put TE late by half the lines that were dropped --
    7.7 ms on the example sequence, silently.
    """
    full = readout(system, **BIG)
    short = readout(system, **BIG, partial_fourier_pe=0.75)
    assert short.n_echoes < full.n_echoes
    assert short.time_to_echo < full.time_to_echo
    assert short.echo_spacing == full.echo_spacing
    assert short.resolution_mm == full.resolution_mm
    assert short.time_to_echo == pytest.approx(
        full.time_to_echo - (full.echo_index - short.echo_index) * full.echo_spacing, abs=1e-9
    )
    # Halfway through the train is the wrong answer, and wrong by a knowable amount.
    midpoint = short.prephase_duration + 0.5 * short.train_duration
    assert midpoint - short.time_to_echo > 7e-3


def test_partial_fourier_is_snapped_up_to_divide_the_shots(system) -> None:
    """
    The line count must divide `n_shots`, or the shots would differ in length and so in timing.

    Snapped by acquiring **more** lines, never fewer, and the attribute reports what was achieved
    rather than what was asked -- the convention ``CartesianLine.partial_echo`` already follows.
    """
    ro = readout(system, partial_fourier_pe=0.7, n_shots=4)
    assert ro.n_echoes * 4 == sum(len(ro.lines(s)) for s in range(4))
    assert ro.partial_fourier_pe >= 0.7
    assert ro.partial_fourier_pe * ro.matrix_pe == pytest.approx(
        round(ro.partial_fourier_pe * ro.matrix_pe)
    )

    with pytest.raises(ConfigurationError, match='cannot be split across'):
        sc.modules.EPIReadout(
            system, fov_ro_mm=240, matrix_ro=64, fov_pe_mm=240, matrix_pe=66, n_shots=4,
            regime='epi',
        )


# ------------------------------------------------------------------ what sets the echo spacing
def test_echo_spacing_is_slew_limited_not_amplitude_limited(system) -> None:
    """
    Doubling the slew shortens the echo by ``sqrt(2)``; doubling the amplitude changes nothing.

    The second half is the interesting one.  For a lobe of fixed area the minimum-time trapezoid is
    a triangle at ``G = sqrt(A S)``, and at 240 mm on a 128 matrix that is 49 mT/m against a
    170 mT/m limit -- so the amplitude limit is nowhere near binding, and a caller reaching for a
    readout duration is reaching for the wrong knob.
    """
    base = readout(system, **BIG)
    faster = sc.System.preset('cima_x').derate('epi', grad=0.85, slew=2.0)
    stronger = sc.System.preset('cima_x').derate('epi', grad=1.7, slew=1.0)

    quick = sc.modules.EPIReadout(faster, **BIG, regime='epi')
    # The echo spacing is two raster-rounded ramps, so it is quantised in steps of two rasters;
    # comparing to a real-valued sqrt(2) needs that granularity allowed for, not a relative epsilon.
    assert abs(quick.echo_spacing - base.echo_spacing / math.sqrt(2.0)) <= (
        2.0 * system.grad_raster.dt
    )
    assert quick.flat_time == 0.0 and base.flat_time == 0.0

    same = sc.modules.EPIReadout(stronger, **BIG, regime='epi')
    assert same.echo_spacing == base.echo_spacing
    assert same.peak_amplitude_Hz_per_m == pytest.approx(base.peak_amplitude_Hz_per_m, rel=1e-9)


def test_a_flat_top_appears_only_when_the_amplitude_limit_binds(system) -> None:
    """
    On a weak-gradient system the lobe stops being a triangle, which is the classic trapezoid.

    The solve has to switch from lengthening the ramp -- which cannot help once the amplitude is
    capped -- to growing a flat top, and that is a different branch rather than a longer search.
    """
    weak = sc.System.from_limits(sc.core.system.Limits(8.0, 400.0), name='weak gradients')
    ro = sc.modules.EPIReadout(weak, fov_ro_mm=240, matrix_ro=128, fov_pe_mm=240, matrix_pe=128)
    assert ro.flat_time > 0.0
    # Running at the cap, but not exactly on it: the flat top grows a whole raster at a time, so the
    # first one that fits overshoots slightly and the amplitude lands just under.
    ceiling = float(ro.opts.max_grad)
    assert 0.95 * ceiling < ro.peak_amplitude_Hz_per_m <= ceiling
    assert ro.ramp_time == pytest.approx(
        weak.grad_raster.ceil(ceiling / float(ro.opts.max_slew)), abs=1e-12)


def test_a_forced_flat_top_lowers_the_amplitude(system) -> None:
    """A deliberately derated readout: longer echo, lower gradient, quieter and easier on PNS."""
    base = readout(system, **BIG)
    slow = readout(system, **BIG, flat_time_us=400.0)
    assert slow.flat_time == pytest.approx(400e-6)
    assert slow.echo_spacing > base.echo_spacing
    assert slow.peak_amplitude_Hz_per_m < base.peak_amplitude_Hz_per_m


# --------------------------------------------------------------------------------- the sampling
def test_the_default_dwell_is_the_longest_that_satisfies_nyquist(system) -> None:
    """
    A finer dwell buys nothing and costs samples; a coarser one aliases along the readout.

    The constraint is circular -- a longer dwell leaves a bigger unsampled remainder, which pushes
    sampling into the ramp, which raises the amplitude the extent needs, which tightens Nyquist --
    so the answer is found by shrinking until legal and then growing back.
    """
    ro = readout(system, **BIG)
    fov_ro_m = ro.fov_ro_mm / 1e3
    spacing = ro.peak_amplitude_Hz_per_m * ro.dwell_s
    assert spacing <= 1.0 / fov_ro_m, 'the chosen dwell aliases along the readout'

    one_raster_longer = ro.dwell_s + system.adc_raster.dt
    with pytest.raises(ConfigurationError, match='aliases along the readout'):
        readout(system, **BIG, adc_dwell_ns=one_raster_longer * 1e9)


def test_the_sample_count_respects_the_divisor_and_the_event_limit(system) -> None:
    """Siemens needs a multiple of four samples, and no more than 8192 in one ADC event."""
    ro = readout(system, **BIG)
    assert ro.samples_per_echo % system.adc_samples_divisor == 0
    assert ro.samples_per_echo <= system.limits('epi').adc_samples_limit

    # A long flat top plus the finest dwell the raster offers: 10 660 samples in one ADC event,
    # which the interpreter refuses.  The default dwell can never reach it -- Nyquist keeps it
    # coarse -- so this is only reachable by asking for both.
    with pytest.raises(ConfigurationError, match='samples per echo, above'):
        readout(system, fov_ro_mm=240, matrix_ro=128, fov_pe_mm=240, matrix_pe=128,
                flat_time_us=1000.0, adc_dwell_ns=100.0)


def test_a_blip_shorter_than_its_own_floors_is_refused(system) -> None:
    """
    Three floors, and the message names which one binds.

    On this system the blip's own **slew** limit asks for 60 us while the ADC dead time asks for
    only 20; on a lower-slew system, or with the larger blips a segmented acquisition needs, the
    ordering swaps -- so the error reports all three rather than guessing.
    """
    ro = readout(system, **BIG)
    assert max(ro.blip_floors, key=lambda k: ro.blip_floors[k]) == 'slew'
    assert ro.blip_duration == pytest.approx(max(ro.blip_floors.values()))
    assert system.grad_raster.count(ro.blip_duration) % 2 == 0, 'blip/2 must land on the raster'

    with pytest.raises(ConfigurationError, match='shorter than the'):
        readout(system, **BIG, blip_duration_us=20.0)


def test_a_bigger_blip_is_needed_for_more_shots(system) -> None:
    """Its area is ``n_shots / FOV``, so the floor rises with the segmentation."""
    one = readout(system, **BIG)
    four = readout(system, **BIG, n_shots=4)
    assert four.blip_duration > one.blip_duration
    assert four.n_echoes * 4 == one.n_echoes


# ------------------------------------------------------------------------------- the trajectory
@pytest.mark.parametrize('n_shots', [1, 2, 4])
@pytest.mark.parametrize('partial_echo', [1.0, 0.8])
def test_the_trajectory_matches_the_compiled_sequence(system, n_shots, partial_echo) -> None:
    """
    The module's own trajectory against pypulseq's, which is the only check a sidecar can get.

    A sidecar that agrees with the module but not with the waveform describes a different sequence,
    and nothing downstream would notice.

    Parametrised over `n_shots` because a single shot did not catch the asymmetric-window bug: with
    one shot the leftover between the window and the samples happened to split evenly, and with two
    it did not.
    """
    ro = readout(system, **BIG, partial_fourier_pe=0.75, n_shots=n_shots, partial_echo=partial_echo)
    out = spin_echo(system, ro)
    mine = ro.trajectory(0)
    theirs = out.kspace()['k_adc'][:2, : ro.n_samples]
    worst = float(np.abs(mine - theirs).max())
    assert worst < 0.01 * ro.dk_pe_per_m, f'trajectory disagrees by {worst:.4f} 1/m'


@pytest.mark.parametrize('n_shots', [1, 2, 4, 8, 16])
def test_the_sampled_window_is_centred_on_the_lobe(system, n_shots) -> None:
    """
    Equal unsampled time before and after the samples, exactly.

    Otherwise the leading and trailing unsampled corners of the lobe differ, and they no longer
    cancel between a forward echo and the reversed one after it: the reversed echoes' ``kx`` drifts
    by the difference.  Measured at 800 ns of asymmetry on a two-shot design -- **0.28 1/m, 6.8 % of
    dk** -- while a single shot was exact, which is why this is parametrised.

    A k-space offset that alternates echo to echo is the worst kind: it is a phase ramp of
    alternating sign, which is a Nyquist ghost, and it looks exactly like an uncalibrated gradient
    delay rather than like an arithmetic mistake.
    """
    ro = readout(system, **BIG, partial_fourier_pe=0.75, n_shots=n_shots)
    lead = ro.sample_offset
    trail = ro.echo_spacing - (ro.sample_offset + ro.samples_per_echo * ro.dwell_s)
    assert lead == pytest.approx(trail, abs=1e-15), (
        f'window asymmetric by {(lead - trail) * 1e9:.0f} ns, which is '
        f'{abs(lead - trail) * ro.peak_amplitude_Hz_per_m / ro.ramp_time * lead / ro.dk_pe_per_m:.2f}'
        f' dk of alternating kx offset'
    )
    assert lead >= ro.blip_duration / 2.0
    assert ro.system.rf_raster.holds(ro.sample_offset), 'pypulseq needs the ADC delay on the RF raster'

    # And the consequence: every echo covers the same kx interval, forward or reversed.
    kx = ro.trajectory(0)[0].reshape(ro.n_echoes, ro.samples_per_echo)
    lows, highs = kx.min(axis=1), kx.max(axis=1)
    assert float(np.ptp(lows)) < 1e-9
    assert float(np.ptp(highs)) < 1e-9


def test_reversing_the_polarity_mirrors_ky_and_nothing_else(system) -> None:
    """The blip-down half of a distortion-correction pair: same timing, mirrored phase encode."""
    ro = readout(system, **BIG, partial_fourier_pe=0.75)
    up, down = ro.trajectory(0), ro.trajectory(0, pe_polarity=-1)
    assert np.allclose(up[0], down[0])
    assert np.allclose(up[1], -down[1])
    assert ro.build().duration == ro.build(pe_polarity=-1).duration

    out = spin_echo(system, ro)
    flipped = sc.compile(
        sc.LogicBlock('down').add(0.0, ro.build(pe_polarity=-1)), system, regime='epi')
    played = flipped.kspace()['k_adc'][1][: ro.n_samples]
    assert np.allclose(played, down[1], atol=0.01 * ro.dk_pe_per_m)
    assert out.check().ok


def test_sample_times_are_referenced_to_the_echo(system) -> None:
    """
    Negative before the echo, which is the whole difference from a spiral.

    A reconstruction that rebased these on the first sample would put ``2 pi df * time_to_echo`` of
    *spatially varying* phase on everything -- not a constant, since df varies with position.
    """
    ro = readout(system, **BIG, partial_fourier_pe=0.75)
    t = ro.sample_times()
    assert t.size == ro.n_samples
    assert np.all(np.diff(t) > 0.0)
    assert t.min() < 0.0 < t.max()
    # The zero crossing lands in the echo `time_to_echo` points at.
    assert int(np.argmin(np.abs(t))) // ro.samples_per_echo == ro.echo_index
    # First sample of the first echo to last sample of the last, which is one echo spacing and one
    # dwell short of the train: sampling starts inside the first lobe and stops inside the last.
    assert t.max() - t.min() == pytest.approx(
        (ro.n_echoes - 1) * ro.echo_spacing + (ro.samples_per_echo - 1) * ro.dwell_s, rel=1e-9
    )


def test_the_distortion_figure_follows_the_train_length(system) -> None:
    """
    ``pe_bandwidth_per_pixel_hz`` is what predicts the displacement, and nothing else does.

    Two shots halve the train, so they double the bandwidth and halve the distortion -- which is the
    whole trade the segmented variant exists to make.
    """
    one = readout(system, **BIG, partial_fourier_pe=0.75)
    two = readout(system, **BIG, partial_fourier_pe=0.75, n_shots=2)
    assert one.pe_bandwidth_per_pixel_hz == pytest.approx(
        1.0 / (one.n_echoes * one.echo_spacing), rel=1e-12)
    assert two.pe_bandwidth_per_pixel_hz == pytest.approx(
        2.0 * one.pe_bandwidth_per_pixel_hz, rel=0.05)
    assert one.total_readout_time_s == pytest.approx(
        (one.n_echoes - 1) * one.echo_spacing, rel=1e-12)


# ------------------------------------------------------------------------------------ labels
def test_labels_address_the_lines_the_table_says(system) -> None:
    """
    ``LIN`` per echo, ``SEG`` per shot, ``REV`` on the reversed ones.

    Label addressing is where an EPI train is most easily wrong: a boundary pushed later by a
    gradient used to put a label in the *previous* readout's block and overwrite that k-space
    address.  The compiler attaches each to the first ADC at or after its own time, and this
    asserts the result rather than the mechanism.
    """
    ro = readout(system, n_shots=2, partial_fourier_pe=0.75)
    out = spin_echo(system, ro)
    assert out.check().ok

    seen: list[int] = []
    reversed_flags: list[bool] = []
    for index in sorted(out.seq.block_events):
        block = out.seq.get_block(index)
        labels = getattr(block, 'label', None)
        if labels is None or getattr(block, 'adc', None) is None:
            continue
        for label in (labels if isinstance(labels, (list, tuple)) else [labels]):
            if label.label == 'LIN':
                seen.append(int(label.value))
            elif label.label == 'REV':
                reversed_flags.append(bool(label.value))
    assert tuple(seen) == ro.lines(0)
    assert reversed_flags == [echo % 2 == 1 for echo in range(ro.n_echoes)]


def test_labels_can_be_turned_off(system) -> None:
    """A sequence with its own labelling scheme should not have to unpick the module's."""
    plain = readout(system, labels=False).build()
    assert not [n for n in plain if getattr(n.item, 'type', None) == 'labelset']


# ------------------------------------------------------------------------------- refusals
def test_the_two_axes_must_differ(system) -> None:
    """Readout and phase encode on one axis is a typo, not a sequence."""
    with pytest.raises(ConfigurationError, match='axes must differ'):
        readout(system, axes=('x', 'x'))


def test_prephasers_can_be_placed_separately(system) -> None:
    """
    The reason ``prephase=False`` exists: three winders on three axes, one block, no coordination.

    The readout prephaser on x, the phase-encode prephaser on y and a slice rewinder on z coincide,
    and the compiler is silent because they are on different axes.
    """
    ro = readout(system, prephase=False)
    assert ro.prephase_duration == 0.0
    assert ro.time_to_echo < readout(system).time_to_echo

    exc = sc.modules.SincExcitation(
        system, flip_deg=90, duration_us=2000, slice_thickness_mm=4, rephase=False, regime='epi')
    tree = sc.LogicBlock('overlap')
    tree.add(0.0, exc.build())
    tree.add(exc.duration, exc.rephaser())
    tree.add(exc.duration, ro.prephaser_block())
    out = sc.compile(tree, system, regime='epi')
    assert not [i for i in out.check().issues if i.kind == 'grad_merge']
    assert out.check().ok


def test_the_readout_train_is_two_gradient_events(system) -> None:
    """
    One per axis, whatever the echo count.  This is the decision the compile hygiene rests on.

    Stated as a test because the obvious implementation -- one trapezoid per echo -- also works and
    is what produces 96 warnings a shot.
    """
    ro = readout(system, **BIG, partial_fourier_pe=0.75, labels=False)
    events = [event for _, event, _ in sc.flatten(ro.build())]
    gradients = [e for e in events if getattr(e, 'type', None) in ('trap', 'grad')]
    # the two prephasers plus the two train events, and no more however many echoes there are
    assert len(gradients) == 4
    assert ro.n_echoes == 96
    assert len([e for e in events if getattr(e, 'type', None) == 'adc']) == ro.n_echoes
