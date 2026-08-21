"""
``EPI2D`` -- the symmetry that stops an EPI ghosting, measured rather than asserted about.

The module's whole claim is that the two readout polarities sample **one** grid.  It is checked
here against ``calculate_kspacePP``, which integrates the *compiled* waveform at the true ADC
sample times and knows nothing about the arithmetic that produced it -- the same independent
oracle ``test_cartesian_line.py`` uses, and for a sharper reason: a half-sample asymmetry between
even and odd echoes is an N/2 ghost, it is invisible in the block structure, and against a
symmetric phantom the image still looks like an image.

Every geometric assertion runs over **six modes** -- ramp sampling, flat top, two partial echoes,
a multi-shot blip and navigators -- because each of them moves the guard, and the guard is what
the centring is built out of.  The protocol is the reference one: 220 mm, 64 x 64, 3.5 us dwell at
``oversampling=2``, 40 mT/m, 180 T/m/s, 10 us receive dead time.

The other half of the module is what the *compiler* does with it.  A blip centred on the seam
straddles the only boundary the compiler can cut between two ADCs, so the barrier that pins the
boundary to the seam is load-bearing, and its failure is silent: the readout lobes come out as
arbitrary gradients and only a merge warning says so.  Hence two tests on gradient *types* and one
that asserts the merge warning rather than suppressing it.
"""

from __future__ import annotations

import warnings

import numpy as np
import pypulseq as pp
import pytest
from pypulseq.opts import Opts

import seqcraft as sc

FOV_MM = 220.0
MATRIX = (64, 64)
DWELL_S = 3.5e-6                    # the reference protocol: 2232 Hz per reconstructed pixel


@pytest.fixture(scope='module')
def epi_opts() -> Opts:
    """A scanner an EPI can actually be built on: the shared fixture with a faster slew rate."""
    return Opts(
        max_grad=40, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s', B0=3.0,
        rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
    )


def _build(opts: Opts, **kwargs) -> sc.modules.EPI2D:
    """The reference readout, with anything the caller wants to vary replacing its default."""
    base = dict(opts=opts, fov_mm=FOV_MM, matrix=MATRIX, dwell_s=DWELL_S)
    if 'bandwidth_hz_px' in kwargs:
        base.pop('dwell_s')
    return sc.modules.EPI2D(**{**base, **kwargs})


def _gradient_types(seq) -> dict[str, set[str]]:
    """Which pulseq representation each axis came out as -- a split lobe stops being a ``trap``."""
    types: dict[str, set[str]] = {}
    for index in range(1, len(seq.block_events) + 1):
        block = seq.get_block(index)
        for axis in ('gx', 'gy', 'gz'):
            event = getattr(block, axis, None)
            if event is not None:
                types.setdefault(axis, set()).add(getattr(event, 'type', '?'))
    return types


def _lines(mode: dict) -> list[int]:
    """Eight echoes this mode's blip can actually step between."""
    step = mode.get('blip_lines', 1)
    return list(range(0, 8 * step, step))


def _labels(block) -> dict[str, list[int]]:
    """Every ``SET`` label the block emits, in the order it emits them."""
    found: dict[str, list[int]] = {}
    for _, event, _ in sc.flatten(block):
        if getattr(event, 'type', '') == 'labelset':
            found.setdefault(event.label, []).append(int(event.value))
    return found


def _k_at_samples(epi: sc.modules.EPI2D, opts: Opts, lines) -> tuple[np.ndarray, np.ndarray]:
    """``(kx, ky)`` at every ADC sample of an imaging train, from pypulseq's own integration."""
    table = list(lines)
    tree = sc.LogicBlock('probe').add(0.0, epi(lines=table))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        k = sc.kspace(tree, opts)
    total = epi.navigator_echoes + len(table)
    shaped = k['k_adc'].reshape(3, total, epi.num_samples)
    return shaped[0][epi.navigator_echoes:], shaped[1][epi.navigator_echoes:]


# ------------------------------------------------------------------- the centred window
MODES = [
    pytest.param({}, id='ramp_sampling'),
    pytest.param({'ramp_sampling': False}, id='flat_top'),
    pytest.param({'partial_fourier': 0.75}, id='partial_echo'),
    pytest.param({'partial_fourier': 0.625}, id='partial_echo_hard'),
    pytest.param({'blip_lines': 4}, id='multi_shot_blip'),
    pytest.param({'navigator_echoes': 3}, id='navigators'),
]


@pytest.mark.parametrize('mode', MODES)
def test_the_sampling_window_is_exactly_centred_in_its_lobe(epi_opts, mode) -> None:
    """
    ``T == 2*guard + N*dwell``, exactly, in every mode.

    This is the whole design.  A lobe whose flat top was rounded up onto the gradient raster with
    the slack left at one end satisfies ``T >= 2*guard + N*dwell`` and fails this, and the
    failure is an N/2 ghost of the sequence's own making.
    """
    epi = _build(epi_opts, **mode)

    window = 2 * epi.guard_s + epi.num_samples * epi.dwell_s

    assert window == pytest.approx(epi.echo_spacing_s, abs=1e-12)


@pytest.mark.parametrize('mode', MODES)
def test_the_two_polarities_sample_the_same_k_grid(epi_opts, mode) -> None:
    """
    A reverse echo's samples, reversed, land on the forward echo's positions.

    Measured on the compiled waveform rather than on ``k_read_per_m``, because the array and the
    trajectory agreeing with each other proves nothing: both come from this module.
    """
    epi = _build(epi_opts, **mode)
    lines = _lines(mode)

    kx, _ = _k_at_samples(epi, epi_opts, lines)

    forward = kx[0] if epi.polarity(0) > 0 else kx[0][::-1]
    for echo in range(len(lines)):
        here = kx[echo] if epi.polarity(echo) > 0 else kx[echo][::-1]
        assert np.abs(here - forward).max() < 1e-4 * epi.dk_read_per_m


@pytest.mark.parametrize('mode', MODES)
def test_k_is_zero_at_the_echo_sample_of_every_echo(epi_opts, mode) -> None:
    """
    ``kx`` is zero at :attr:`pre_echo_samples` on a forward echo and at its mirror on a reverse one.

    The mirror is the assertion that matters: ``num_samples - 1 - pre_echo_samples`` is only the
    echo of a reverse lobe if the window is centred, so this fails for the same reason the test
    above does and says a different thing about it.
    """
    epi = _build(epi_opts, **mode)
    lines = _lines(mode)

    kx, _ = _k_at_samples(epi, epi_opts, lines)

    for echo in range(len(lines)):
        assert kx[echo, epi.echo_sample(echo)] == pytest.approx(
            0.0, abs=1e-4 * epi.dk_read_per_m)


@pytest.mark.parametrize('mode', MODES)
def test_the_declared_trajectory_is_the_one_that_plays(epi_opts, mode) -> None:
    """``k_read_per_m`` is what a reconstruction grids from, so it is checked against the file."""
    epi = _build(epi_opts, **mode)
    lines = _lines(mode)

    kx, _ = _k_at_samples(epi, epi_opts, lines)

    for echo in range(len(lines)):
        declared = epi.k_read_per_m if epi.polarity(echo) > 0 else epi.k_read_per_m[::-1]
        assert np.abs(declared - kx[echo]).max() < 1e-4 * epi.dk_read_per_m


# ----------------------------------------------------------------- what each mode means
def test_ramp_sampling_shortens_the_lobe_and_raises_the_peak(epi_opts) -> None:
    """
    The trade, in one test: less time per echo, more amplitude, the same sampled extent.

    The extent is the part worth pinning.  A ramp-sampled lobe that traversed a different amount
    of k would change the field of view rather than the echo spacing, and the image would come
    back the right shape at the wrong scale.
    """
    ramp, flat = _build(epi_opts), _build(epi_opts, ramp_sampling=False)

    assert ramp.echo_spacing_s < flat.echo_spacing_s
    assert ramp.readout_amplitude_hz_m > flat.readout_amplitude_hz_m
    assert np.ptp(np.diff(flat.k_read_per_m)) == pytest.approx(0.0, abs=1e-9)
    assert np.ptp(np.diff(ramp.k_read_per_m)) > 0.3 * ramp.dk_read_per_m
    # Both traverse one dk per *reconstructed* sample over the window, which is what fixes the
    # field of view.  Oversampling divides the spacing and leaves this alone.
    for epi in (ramp, flat):
        span = np.ptp(epi.k_read_per_m)
        assert span == pytest.approx(
            (epi.num_samples / epi.oversampling - 1) * epi.dk_read_per_m, rel=2e-2)


def test_partial_echo_moves_the_echo_and_shortens_the_train(epi_opts) -> None:
    """Fewer pre-echo samples is a shorter lobe, and therefore an earlier echo."""
    full, partial = _build(epi_opts), _build(epi_opts, partial_fourier=0.75)

    assert partial.num_samples == 96                    # 0.75 * 64, even, times oversampling
    assert partial.pre_echo_samples == 32
    assert full.pre_echo_samples == 64
    assert partial.echo_spacing_s < full.echo_spacing_s
    assert partial.time_to_echo(0) < full.time_to_echo(0)


def test_the_sample_count_is_oversampling_times_an_even_number(epi_opts) -> None:
    """Odd counts leave the guard half an ADC raster short, so they are rounded away."""
    for fraction in (0.61, 0.64, 0.67, 0.72, 0.79):
        epi = _build(epi_opts, partial_fourier=fraction)
        assert epi.num_samples % (2 * epi.oversampling) == 0


def test_oversampling_costs_no_echo_spacing_at_all(epi_opts) -> None:
    """
    Twice the samples at half the dwell is the same sampling duration, so the lobe does not move.

    One assertion, because that identity is the entire argument for defaulting it to 2: the
    readout is unchanged and only the sample *spacing* improves.  The bandwidth is quoted per
    *reconstructed* pixel, which is why it can be held fixed across the three.
    """
    built = [_build(epi_opts, bandwidth_hz_px=1953.125, oversampling=n) for n in (1, 2, 4)]

    for epi in built[1:]:
        assert (epi.echo_spacing_s, epi.guard_s, epi.readout_amplitude_hz_m) == (
            built[0].echo_spacing_s, built[0].guard_s, built[0].readout_amplitude_hz_m)
        assert epi.bandwidth_hz_px == built[0].bandwidth_hz_px
    assert [epi.num_samples for epi in built] == [64, 128, 256]
    assert [epi.dwell_s for epi in built] == [8e-6, 4e-6, 2e-6]


def test_oversampling_is_what_brings_the_worst_gap_below_one_dk(epi_opts) -> None:
    """
    The sampling-theory reason for the default, and it comes before any question of interpolation.

    A ramp-sampled lobe is dense on the ramps and sparse on the flat top.  At ``oversampling=1``
    the sparse part is **wider than one reconstruction dk**, so the readout does not sample the
    grid it is reconstructed onto and no interpolator can recover what was never measured.
    """
    one, two = (_build(epi_opts, bandwidth_hz_px=2232.142857142857, oversampling=n)
                for n in (1, 2))

    gaps = [np.diff(epi.k_read_per_m).max() / epi.dk_read_per_m for epi in (one, two)]

    assert gaps[0] == pytest.approx(1.119, abs=5e-3)
    assert gaps[1] == pytest.approx(0.559, abs=5e-3)


def test_the_dwell_is_snapped_up_onto_a_quantum_coarser_than_the_adc_raster(epi_opts) -> None:
    """
    The guard is the ADC's delay, and pypulseq measures an ADC delay in **RF** rasters.

    So ``num_samples * dwell_s`` has to be an even number of RF rasters, which at 128 samples
    makes the dwell a multiple of 500 ns rather than of the 100 ns ADC raster.  Snapped **up**,
    because rounding down raises the bandwidth, shortens the lobe, and turns a feasible protocol
    into a refusal.
    """
    epi = _build(epi_opts, dwell_s=3.1e-6)

    assert epi.dwell_s == 3.5e-6
    assert epi.bandwidth_hz_px == pytest.approx(1.0 / (3.5e-6 * MATRIX[0] * epi.oversampling))
    # And the point of the quantum: the guard lands on the RF raster, exactly.
    assert round(epi.guard_s / epi_opts.rf_raster_time) == pytest.approx(
        epi.guard_s / epi_opts.rf_raster_time)


def test_a_bigger_blip_costs_echo_spacing(epi_opts) -> None:
    """
    `blip_lines` prices an ordering at design time, which is the reason it is an argument.

    The price is asserted rather than only its direction: a centric single shot steps 63 lines at
    one point, every blip in the train is designed for that worst step, and the echo spacing it
    costs is 1.7x the linear ordering's.  No table of integers carries that on its own.
    """
    spacings = [round(_build(epi_opts, blip_lines=n).echo_spacing_s * 1e6)
                for n in (1, 2, 4, 8, 63)]

    assert spacings == [510, 530, 550, 590, 850]


def test_the_guard_holds_half_a_blip_and_the_receiver_dead_time(epi_opts) -> None:
    """
    Three claims on one interval, and which of the two terms binds is a property of the scanner.

    At the reference protocol the blip binds -- 60 us of blip against 10 us of dead time.  On a
    receiver with a 60 us dead time the dead time binds instead, and the echo spacing goes from
    510 us to 570.  Both are the same ``max``, so both are checked at once.
    """
    blip_bound = _build(epi_opts)
    slow_receiver = Opts(
        max_grad=40, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s', B0=3.0,
        rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=60e-6,
    )
    dead_time_bound = _build(slow_receiver)

    for epi, opts in ((blip_bound, epi_opts), (dead_time_bound, slow_receiver)):
        assert epi.guard_s >= max(epi.blip_duration_s / 2, opts.adc_dead_time) - 1e-12
    assert (blip_bound.guard_s, round(blip_bound.echo_spacing_s * 1e6)) == (31e-6, 510)
    assert (dead_time_bound.guard_s, round(dead_time_bound.echo_spacing_s * 1e6)) == (61e-6, 570)


def test_a_guard_exactly_equal_to_the_dead_time_still_compiles(epi_opts) -> None:
    """
    The degenerate case, which needs no margin and no special case -- so it needs a test.

    When ``guard == adc_dead_time`` the interval the compiler has to cut a boundary in collapses
    to a single instant, and that instant is the seam the barrier already occupies.
    """
    exact = Opts(
        max_grad=40, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s', B0=3.0,
        rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=31e-6,
    )
    epi = _build(exact)
    assert epi.guard_s == exact.adc_dead_time

    tree = sc.LogicBlock('probe').add(0.0, epi(lines=range(8)))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        seq, k = sc.compile(tree, exact), sc.kspace(tree, exact)

    assert len(seq.block_events) == 9
    assert _gradient_types(seq)['gx'] == {'trap'}
    kx = k['k_adc'].reshape(3, 8, epi.num_samples)[0]
    assert max(abs(kx[n, epi.echo_sample(n)]) for n in range(8)) < 1e-4 * epi.dk_read_per_m


def test_the_blip_lands_on_an_even_number_of_gradient_rasters(epi_opts) -> None:
    """
    A centred blip starts at ``seam - blip/2``, and a gradient may not start off the raster.

    So the duration is rounded onto *two* rasters, not one.  An odd raster count is the one way
    to make this design fail at compile time rather than in the image, and a request between two
    even counts is the direct way to ask for one.
    """
    epi = _build(epi_opts, blip_duration_s=70e-6)

    assert epi.blip_duration_s == pytest.approx(80e-6)
    assert round(epi.blip_duration_s / epi_opts.grad_raster_time) % 2 == 0


def test_every_blip_takes_the_same_time_whatever_it_steps(epi_opts) -> None:
    """The echo spacing is a constant of the instance, not a property of the table."""
    epi = _build(epi_opts, blip_lines=4)

    tables = ([0, 4, 8, 12], [0, 1, 2, 3], [0, 4, 5, 9], [12, 8, 4, 0])
    durations = {round(epi(lines=table).duration, 12) for table in tables}

    assert len(durations) == 1


def test_navigators_sample_the_centre_of_the_blip_axis(epi_opts) -> None:
    """
    Three lobes before the blip-axis prephaser, at ``ky = 0``.

    At the edge of k-space a navigator measures the same odd/even phase multiplied by whatever
    signal is there, which on a real object is close to none -- so where the prephaser goes is
    the whole point of the argument.
    """
    epi = _build(epi_opts, navigator_echoes=3)
    lines = list(range(8))

    tree = sc.LogicBlock('probe').add(0.0, epi(lines=lines))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        ky = sc.kspace(tree, epi_opts)['k_adc'].reshape(3, 3 + len(lines), epi.num_samples)[1]

    assert np.abs(ky[:3]).max() < 1e-4 * epi.dk_blip_per_m
    assert ky[3].mean() == pytest.approx(epi.k_blip_per_m(lines[0]), abs=1e-4 * epi.dk_blip_per_m)
    # An odd number of navigators starts the imaging train on a reverse lobe, and REV says so.
    assert epi.polarity(0) == -1


def test_a_navigator_is_a_negative_echo_index(epi_opts) -> None:
    """
    The whole train is addressable, so a reconstruction can write ``range(-navigators, len(lines))``
    and get the file's readouts in order.

    ``time_to_echo`` is the reason this is a method and not ``echo * echo_spacing_s``: the
    blip-axis prephaser plays *between* the last navigator and imaging echo zero, so a navigator is
    further back than its index suggests.  Measured against the true ADC sample times: 0 ps.
    """
    epi = _build(epi_opts, navigator_echoes=3)
    lines = list(range(8))
    nav = epi.navigator_echoes

    tree = sc.LogicBlock('probe').add(0.0, epi(lines=lines))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        k = sc.kspace(tree, epi_opts)
    rows = nav + len(lines)
    kx = k['k_adc'].reshape(3, rows, epi.num_samples)[0]
    t_adc = k['t_adc'].reshape(rows, epi.num_samples)

    assert [epi.polarity(e) for e in range(-nav, 3)] == [1, -1, 1, -1, 1, -1]
    for row, echo in enumerate(range(-nav, len(lines))):
        sample = epi.echo_sample(echo)
        assert abs(kx[row, sample]) < 1e-4 * epi.dk_read_per_m, f'echo {echo}'
        assert t_adc[row, sample] == pytest.approx(epi.time_to_echo(echo), abs=1e-12)

    # A navigator that was never built is still refused, and the message names the range.
    with pytest.raises(sc.ConfigurationError, match='before the start of the train'):
        epi.polarity(-nav - 1)
    with pytest.raises(sc.ConfigurationError, match='builds no navigator echoes'):
        _build(epi_opts).polarity(-1)


def test_the_labels_carry_the_line_and_the_polarity(epi_opts) -> None:
    """A reconstruction reads the file, not the module, so the file has to say which is which."""
    epi = _build(epi_opts, navigator_echoes=2)
    lines = [3, 4, 5]

    labels = _labels(epi(lines=lines, segment=7, reference=True))

    assert labels['LIN'] == lines
    assert labels['SEG'] == [7]
    assert labels['NAV'] == [1, 0]
    assert labels['REV'] == [0, 1, 0, 1, 0]              # two navigators, then three echoes
    # REF and IMA are complementary and say what the *shot* is for.  A reconstruction has to tell
    # a calibration shot from an imaging one, and counting shots is how it comes to guess wrong.
    assert (labels['REF'], labels['IMA']) == ([1], [0])
    assert _labels(epi(lines=lines, reference=False))['REF'] == [0]
    assert 'REF' not in _labels(epi(lines=lines))


def test_a_dummy_shot_plays_the_same_gradients_and_no_labels(epi_opts) -> None:
    """The gradients are what establish the state a real shot then samples."""
    epi = _build(epi_opts)

    real, dummy = epi(lines=range(8)), epi(lines=range(8), acquire=False)
    assert dummy.duration == pytest.approx(real.duration, abs=1e-12)
    assert sc.moments(dummy, order=0) == pytest.approx(sc.moments(real, order=0))
    assert not [e for _, e, _ in sc.flatten(dummy) if getattr(e, 'type', '') in ('adc', 'labelset')]
    assert len([e for _, e, _ in sc.flatten(real) if getattr(e, 'type', '') == 'adc']) == 8


def test_the_naive_lobe_is_the_bug_this_module_exists_to_avoid(epi_opts) -> None:
    """
    The failure mode as a test, because a test that only asserts the fix passes on a new bug.

    ``CartesianLine`` rounds the flat top up onto the gradient raster and leaves the slack after
    the last sample, which is correct for a *line*.  Do that in a *train* and the forward and
    reverse lobes sample two grids up to half a gradient raster apart -- an N/2 ghost the sequence
    made itself, invisible in the block structure and indistinguishable from the real artefact.
    """
    epi = _build(epi_opts, ramp_sampling=False)
    raster = epi_opts.grad_raster_time
    sample_time = epi.num_samples * epi.dwell_s
    slack = float(np.ceil(sample_time / raster) * raster - sample_time)
    assert slack == pytest.approx(2e-6)            # 448 us of samples in a 450 us flat top

    offsets = (np.arange(epi.num_samples) + 0.5) * epi.dwell_s
    centred = epi.guard_s + offsets                            # the module: slack split in two
    naive = (epi.guard_s - slack / 2) + offsets                # the line: slack left at the end

    # ``t_i + t_{N-1-i} == T`` is the whole design, and it is what makes a reverse lobe sample the
    # forward lobe's grid.  The naive window misses it by the slack, every sample, every echo.
    assert np.abs(centred + centred[::-1] - epi.echo_spacing_s).max() < 1e-12
    assert np.abs(naive + naive[::-1] - epi.echo_spacing_s).max() == pytest.approx(slack)


# --------------------------------------------------------------------------- the timing
def test_time_to_echo_names_the_sample_the_echo_actually_is(epi_opts) -> None:
    """The one question a tree cannot answer, checked against the true ADC sample times."""
    epi = _build(epi_opts)
    lines = list(range(12))

    tree = sc.LogicBlock('probe').add(0.0, epi(lines=lines))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        t_adc = sc.kspace(tree, epi_opts)['t_adc'].reshape(len(lines), epi.num_samples)

    for echo in range(len(lines)):
        assert t_adc[echo, epi.echo_sample(echo)] == pytest.approx(
            epi.time_to_echo(echo), abs=1e-12)


def test_the_echo_time_anchor_moves_with_the_table(epi_opts) -> None:
    """Which echo carries ``k = 0`` is the ordering's decision, so this is a method."""
    epi = _build(epi_opts)

    linear = epi.time_to_center_line(range(64))
    reverse = epi.time_to_center_line(list(range(64))[::-1])

    assert linear == epi.time_to_echo(32)
    assert reverse == epi.time_to_echo(31)
    # One echo apart in the table and one echo spacing *plus one dwell* apart in time.  The two
    # echoes have opposite polarity, so their k = 0 samples sit on opposite sides of their lobes'
    # midpoints -- and the offset is half a dwell each way, because k = 0 belongs at sample
    # N // 2 rather than at the middle of the window.  A query returning n * echo_spacing would
    # be exactly one dwell wrong, every other echo.
    assert linear - reverse == pytest.approx(epi.echo_spacing_s + epi.dwell_s, abs=1e-12)


def test_the_prephaser_cancels_exactly_what_the_first_lobe_accumulates(epi_opts) -> None:
    """``area_to_echo_per_m`` is the physics number, and ``prephase=False`` leaves it available."""
    with_prephaser, without = _build(epi_opts), _build(epi_opts, prephase=False)

    assert without.area_to_echo_per_m == pytest.approx(with_prephaser.area_to_echo_per_m)
    prephaser = next(n.item for n in with_prephaser(lines=[0])
                     if getattr(n.item, 'type', '') == 'trap' and n.item.channel == 'x')
    assert float(prephaser.area) == pytest.approx(-with_prephaser.area_to_echo_per_m, rel=1e-9)


def test_prephase_false_starts_the_block_at_the_first_lobe(epi_opts) -> None:
    """The caller places the dephasers, so the module's own block begins with the readout."""
    with_prephaser, without = _build(epi_opts), _build(epi_opts, prephase=False)

    gap = with_prephaser.time_to_echo(0) - without.time_to_echo(0)

    assert gap == pytest.approx(with_prephaser.prephaser_duration_s, abs=1e-12)
    assert not [e for _, e, _ in sc.flatten(without(lines=[0]))
                if getattr(e, 'type', '') == 'trap' and e.channel == 'y']


def test_a_dephaser_placed_before_a_refocusing_pulse_still_lands_on_k_zero(epi_opts) -> None:
    """
    The spin-echo placement, end to end: the conjugation flips both dephasers.

    A sign error here is a whole readout's worth of ``k``, every waveform is legal, and against a
    symmetric phantom the image still looks like an image -- so it is asserted through
    ``calculate_kspacePP``, which is the only thing that models the conjugation.
    """
    epi = _build(epi_opts, prephase=False)
    exc = sc.modules.Excitation(opts=epi_opts, flip_deg=90.0, thickness_mm=5.0, duration_s=2e-3)
    refoc = sc.modules.Refocusing(opts=epi_opts, thickness_mm=6.25, duration_s=3e-3,
                                  crush_cycles_per_voxel=4.0, crush_voxel_mm=5.0)
    line = 20
    # Opposite sign on both, because a refocusing pulse maps k to -k.
    deph_x = pp.make_trapezoid('x', area=epi.area_to_echo_per_m, duration=600e-6, system=epi_opts)
    deph_y = pp.make_trapezoid('y', area=-epi.k_blip_per_m(line), duration=600e-6, system=epi_opts)

    tree = sc.LogicBlock('probe').add([
        [0.0, exc()],
        [exc().duration, deph_x, deph_y],
        [exc().duration + 600e-6, refoc()],
        [exc().duration + 600e-6 + refoc().duration, epi(lines=[line])],
    ])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        k = sc.kspace(tree, epi_opts)['k_adc']

    assert k[0, epi.echo_sample(0)] == pytest.approx(0.0, abs=1e-4 * epi.dk_read_per_m)
    assert k[1, epi.echo_sample(0)] == pytest.approx(epi.k_blip_per_m(line),
                                                     abs=1e-4 * epi.dk_blip_per_m)


# ------------------------------------------------------------------------ the compiler
@pytest.mark.parametrize('mode', MODES)
def test_every_readout_lobe_survives_as_a_trapezoid(epi_opts, mode) -> None:
    """
    The barrier's test, and the one that would have caught the silent split.

    The blip is centred on the seam, so it *straddles* the only boundary the compiler can cut
    between two ADCs.  Left to choose, the compiler takes the midpoint of the gap -- which is
    inside the fall ramp -- and every readout lobe in the train comes out as two arbitrary
    gradients instead of one trapezoid.  That compiles, passes every k-space check here and
    simulates correctly; the only report of it is a ``merge`` warning naming the readout axis.

    ``sc.barrier`` at the seam is the fix, and this is what says it worked.
    """
    epi = _build(epi_opts, **mode)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        seq = sc.compile(sc.LogicBlock('probe').add(0.0, epi(lines=_lines(mode))), epi_opts)

    types = _gradient_types(seq)
    assert types['gx'] == {'trap'}
    # The blip axis is where the split goes, and where it is *meant* to go: each block carries the
    # trailing half of one blip and the leading half of the next.
    assert 'grad' in types['gy']


def test_the_only_merge_the_compiler_reports_is_on_the_blip_axis(epi_opts) -> None:
    """
    Asserted rather than suppressed, so that a merge naming the **readout** axis is a failure.

    Summing the two blip halves back together at every boundary is what a centred blip costs, and
    it is reported instead of hidden.  ``FSE2D`` already produces 71 of the same kind.  A merge
    on the readout axis would mean the barrier stopped working, and nothing else would say so.
    """
    epi = _build(epi_opts)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        sc.compile(sc.LogicBlock('probe').add(0.0, epi(lines=range(64))), epi_opts)

    merges = [str(w.message) for w in caught if 'gradient merge' in str(w.message)]

    assert len(merges) == 1
    assert '62 same-axis gradient merges' in merges[0]
    assert f'axis {epi.blip_axis}' in merges[0]
    assert f'axis {epi.axis}' not in merges[0]


def test_the_train_compiles_on_its_own(epi_opts) -> None:
    """A component that only works beside something else is not reusable."""
    epi = _build(epi_opts)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        seq = sc.compile(sc.LogicBlock('probe').add(0.0, epi(lines=range(64))), epi_opts)

    # One winder, then one per echo: the barrier adds no block, because the seam needed a
    # boundary between two ADCs anyway.  It only decides *where*.
    assert len(seq.block_events) == 65


def test_two_calls_produce_identical_events(epi_opts) -> None:
    """Purity: the canonical bug is an involution, so it is called twice."""
    epi = _build(epi_opts, navigator_echoes=3)

    hashes = [
        [sc.events.content_hash(e) for _, e, _ in sc.flatten(epi(lines=range(16)))
         if getattr(e, 'type', '') in ('trap', 'grad', 'adc')]
        for _ in range(3)
    ]

    assert hashes[0] == hashes[1] == hashes[2]


# ------------------------------------------------------------------------ the refusals
def test_a_step_beyond_the_designed_blip_is_refused(epi_opts) -> None:
    """Lengthening one blip silently would move every echo after it."""
    epi = _build(epi_opts, blip_lines=2)

    with pytest.raises(sc.ConfigurationError, match='exceed the 2-line blip'):
        epi(lines=[0, 4])


def test_an_infeasible_readout_names_two_remedies_and_both_work(epi_opts) -> None:
    """
    A message carrying a remedy that does not work is worse than one carrying none.

    So both are checked by *building at them*.  The second is the one nobody guesses: a **longer**
    blip is a longer guard and a longer lobe, and the lobe grows faster than the ramp corners the
    receiver then misses -- so lengthening the blip fixes an infeasible readout as surely as
    lowering the bandwidth does.
    """
    with pytest.raises(sc.ConfigurationError, match='cannot reach') as caught:
        _build(epi_opts, dwell_s=2.5e-6)

    message = str(caught.value)
    suggested = float(message.split('lower bandwidth_hz_px to at most ')[1].split(',')[0])

    assert 'lengthen* blip_duration_s' in message
    assert _build(epi_opts, bandwidth_hz_px=suggested).echo_spacing_s > 0
    assert _build(epi_opts, dwell_s=2.5e-6, blip_duration_s=140e-6).echo_spacing_s > 0
    assert _build(epi_opts, dwell_s=2.5e-6, ramp_sampling=False).echo_spacing_s > 0


def test_an_over_limit_amplitude_is_this_module_s_refusal_and_not_pypulseq_s(epi_opts) -> None:
    """
    The readout lobe is built twice, and **both** builds go through the amplitude check.

    ``pp.make_trapezoid`` refuses an over-limit amplitude with ``ValueError: Amplitude violation
    (133%)`` and nothing else -- no bandwidth, no field of view, no note that ramp sampling puts
    the peak above the nominal ``dk/dwell``, and no indication of which of this module's several
    trapezoids raised it.  Checking the amplitude only *after* the event was built let that escape,
    which is a refusal the caller cannot act on arriving as somebody else's exception.

    A small field of view is the shortest way there: ``dk`` goes up, and with a flat top the
    amplitude *is* ``dk / dwell`` with nothing to trade against it.
    """
    with pytest.raises(sc.ConfigurationError, match='above the') as caught:
        _build(epi_opts, fov_mm=60.0, ramp_sampling=False)

    message = str(caught.value)
    assert 'lower bandwidth_hz_px, or widen fov_mm' in message
    assert 'ramp_sampling' in message
    # And the remedy works: a wider field of view is a smaller dk and a smaller amplitude.
    widened = _build(epi_opts, fov_mm=220.0, ramp_sampling=False)
    assert widened.readout_amplitude_hz_m <= float(epi_opts.max_grad)


@pytest.mark.parametrize(('kwargs', 'match'), [
    ({'blip_axis': 'x'}, 'both'),
    ({'blip_lines': 0}, 'at least 1'),
    ({'oversampling': 0}, 'at least 1'),
    ({'navigator_echoes': -1}, 'not be negative'),
    ({'partial_fourier': 0.4}, 'no pre-echo samples'),
    ({'blip_duration_s': 1e-6}, 'shorter than the minimum'),
    ({'prephaser_duration_s': 1e-6}, 'shorter than the minimum'),
    ({'prephase': False, 'prephaser_duration_s': 1e-3}, 'cannot take effect'),
    ({'dwell_s': None}, 'exactly one of'),
])
def test_the_construction_refusals(epi_opts, kwargs, match) -> None:
    """Each one refused once, so that a message that stops naming its fix is noticed."""
    with pytest.raises(sc.ConfigurationError, match=match):
        _build(epi_opts, **kwargs)


@pytest.mark.parametrize(('lines', 'match'), [
    ([], 'nothing to acquire'),
    ([0, 64], 'outside'),
    ([3, 3], 'more than once'),
])
def test_the_table_refusals(epi_opts, lines, match) -> None:
    with pytest.raises(sc.ConfigurationError, match=match):
        _build(epi_opts)(lines=lines)


def test_a_table_without_the_centre_has_no_echo_time(epi_opts) -> None:
    """Under segmentation only one shot has a meaningful TE, and this is what says which."""
    epi = _build(epi_opts)

    with pytest.raises(sc.ConfigurationError, match='omits the centre'):
        epi.time_to_center_line(range(0, 20))


def test_a_prephaser_question_with_no_prephaser_points_at_the_physics_number(epi_opts) -> None:
    epi = _build(epi_opts, prephase=False)

    with pytest.raises(sc.ConfigurationError, match='area_to_echo_per_m'):
        _ = epi.prephaser_duration_s


def test_the_receiver_takes_the_carrier_phase_its_excitation_was_given(epi_opts) -> None:
    """A spoiled sequence needs `phase_deg`, and it has to reach **every** ADC of the train.

    The transmitter runs a schedule and the receiver is phase-locked to it.  Applying the phase to
    only some echoes -- the imaging ones but not the navigators, say -- is worse than applying it
    to none: it is a constant phase step partway through one shot, which is a discontinuity in
    ``ky`` where none of the physics has one.
    """
    epi = _build(epi_opts, navigator_echoes=2)
    wanted = float(np.deg2rad(117.0))

    for phase_deg, expected in ((117.0, wanted), (0.0, 0.0)):
        adcs = [e for _, e, _ in sc.flatten(epi(lines=range(8), phase_deg=phase_deg))
                if getattr(e, 'type', '') == 'adc']
        assert len(adcs) == 10                      # two navigators and eight imaging echoes
        assert [float(a.phase_offset) for a in adcs] == pytest.approx([expected] * 10, abs=1e-12)


def test_the_carrier_phase_moves_nothing_but_the_phase(epi_opts) -> None:
    """It is a receiver setting, so the gradients and the sample times must not notice it."""
    epi = _build(epi_opts)
    plain = sc.kspace(sc.LogicBlock('p').add(0.0, epi(lines=range(8))), epi_opts)
    spoiled = sc.kspace(sc.LogicBlock('s').add(0.0, epi(lines=range(8), phase_deg=117.0)),
                        epi_opts)

    assert spoiled['k_adc'] == pytest.approx(plain['k_adc'], abs=1e-9)
    assert spoiled['t_adc'] == pytest.approx(plain['t_adc'], abs=1e-12)
