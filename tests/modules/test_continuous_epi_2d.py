"""
``ContinuousEPI2D`` -- one ADC window across a blipped echo-planar train, measured from outside.

Every claim about where a sample lands is checked against an integration of the **compiled**
gradients written in this file, which shares no code with the module, and against
``calculate_kspacePP``.  In memory the two agree with the module to well below 1e-9 of a k-space
step.  A written file stores each gradient amplitude to six significant digits, so read-back
checks use the bound that rounding allows: 5e-6 of the area each event has traversed.

Two designs run through the geometric tests, because where ``k = 0`` falls depends on the parity
of the samples per echo: 5 us at ``oversampling=1.5`` gives 118 samples per echo, and 4 us gives
135.  The scanner is 32 mT/m and 130 T/m/s with a 10 us receiver dead time and an 8192-sample
per-event limit.
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
DK = 1e3 / FOV_MM
EVEN = {'dwell_s': 5e-6, 'oversampling': 1.5}          # 118 samples per echo
ODD = {'dwell_s': 4e-6, 'oversampling': 1.5}           # 135 samples per echo
#: Each written gradient amplitude has six significant digits, so it is within this fraction.
SIX_DIGITS = 5e-6 + 5e-9


def _scanner(**overrides) -> Opts:
    base = dict(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s',
                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
                adc_samples_limit=8192)
    return Opts(**{**base, **overrides})


@pytest.fixture(scope='module')
def cepi_opts() -> Opts:
    return _scanner()


def _build(opts: Opts, **kwargs) -> sc.modules.ContinuousEPI2D:
    """The reference readout, with anything the caller wants to vary replacing its default."""
    base = dict(opts=opts, fov_mm=FOV_MM, matrix=MATRIX, **EVEN)
    return sc.modules.ContinuousEPI2D(**{**base, **kwargs})


def _compiled(block: sc.LogicBlock, opts: Opts):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        return sc.compile(sc.LogicBlock('probe').add(0.0, block), opts)


def _read_back(seq, tmp_path) -> pp.Sequence:
    path = tmp_path / 'train.seq'
    seq.write(str(path))
    back = pp.Sequence()
    back.read(str(path))
    return back


# ------------------------------------------- an integration of the compiled gradients, from here
def _event_knots(event, t0: float) -> tuple[np.ndarray, np.ndarray]:
    """Knots of one ``trap`` or ``grad`` event, read here rather than through SeqCraft."""
    start = t0 + float(event.delay)
    if event.type == 'trap':
        rise, flat, fall = float(event.rise_time), float(event.flat_time), float(event.fall_time)
        amp = float(event.amplitude)
        return (start + np.array([0.0, rise, rise + flat, rise + flat + fall]),
                np.array([0.0, amp, amp, 0.0]))
    tt, wf = np.asarray(event.tt, float), np.asarray(event.waveform, float)
    end = float(event.shape_dur) if hasattr(event, 'shape_dur') else tt[-1]
    if tt[0] > 1e-12:
        tt, wf = np.concatenate([[0.0], tt]), np.concatenate([[float(event.first)], wf])
    if end - tt[-1] > 1e-12:
        tt, wf = np.concatenate([tt, [end]]), np.concatenate([wf, [float(event.last)]])
    return start + tt, wf


def _axis_knots(seq, axis: str) -> tuple[np.ndarray, np.ndarray]:
    t, times, amps = 0.0, [], []
    for index in sorted(seq.block_events):
        grad = getattr(seq.get_block(index), f'g{axis}', None)
        if grad is not None:
            kt, ka = _event_knots(grad, t)
            times.append(kt)
            amps.append(ka)
        t += float(seq.block_durations[index])
    return np.concatenate(times), np.concatenate(amps)


def _integral(knots: tuple[np.ndarray, np.ndarray], at) -> np.ndarray:
    """The exact integral from zero of piecewise-linear knots, at every time in `at`."""
    times, amps = knots
    at = np.atleast_1d(np.asarray(at, float))
    seg = np.diff(times)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (amps[1:] + amps[:-1]) * seg)])
    idx = np.clip(np.searchsorted(times, at, side='right') - 1, 0, len(times) - 2)
    into = np.clip(at - times[idx], 0.0, seg[idx])
    safe = np.where(seg[idx] > 0, seg[idx], 1.0)
    slope = np.where(seg[idx] > 0, (amps[idx + 1] - amps[idx]) / safe, 0.0)
    return cum[idx] + into * (amps[idx] + 0.5 * slope * into)


def _adc_events(seq) -> list[dict]:
    t, found = 0.0, []
    for index in sorted(seq.block_events):
        adc = seq.get_block(index).adc
        if adc is not None:
            found.append({'block_start': t, 'samples': int(adc.num_samples),
                          'dwell': float(adc.dwell), 'delay': float(adc.delay),
                          'block_duration': float(seq.block_durations[index])})
        t += float(seq.block_durations[index])
    return found


def _sample_times(seq) -> np.ndarray:
    return np.concatenate([e['block_start'] + e['delay'] + (np.arange(e['samples']) + 0.5)
                           * e['dwell'] for e in _adc_events(seq)])


def _measured_k(seq) -> tuple[np.ndarray, np.ndarray]:
    """``(2, n)`` k at the sequence's own ADC sample times, and those times."""
    times = _sample_times(seq)
    return np.vstack([_integral(_axis_knots(seq, axis), times) for axis in 'xy']), times


def _six_digit_bound(epi, lines) -> np.ndarray:
    """Per axis, the most a written file's amplitude rounding can move k, in 1/m."""
    k = epi.k_per_m(lines)
    start = np.array([-epi.area_to_echo_per_m, epi.k_blip_per_m(list(lines)[0])])
    swing = np.abs(k - start[:, None]).max(axis=1)
    return SIX_DIGITS * (np.abs(start) + swing) + 1e-9 * DK


# ----------------------------------------------------------------- the window and the samples
def test_one_adc_event_spans_exactly_the_train(cepi_opts, tmp_path) -> None:
    """
    One event per train, whatever its length, opening on the first lobe and closing on the last.

    Measured on the read-back file: the window against the readout gradient's own support.
    """
    epi = _build(cepi_opts)
    back = _read_back(_compiled(epi(lines=range(64)), cepi_opts), tmp_path)

    events = _adc_events(back)
    assert len(events) == 1
    (adc,) = events
    assert adc['samples'] == 64 * epi.samples_per_echo
    assert adc['dwell'] == pytest.approx(epi.dwell_s, abs=1e-15)

    times, amps = _axis_knots(back, 'x')
    inside = times >= adc['block_start']                 # the prephaser block comes first
    lobes = times[inside][np.flatnonzero(np.abs(amps[inside]) > 0)]
    opens = adc['block_start'] + adc['delay']
    closes = opens + adc['samples'] * adc['dwell']
    # The first nonzero knot is one ramp after the lobe starts, the last one ramp before it ends.
    ramp = float(epi.gx.rise_time)
    assert opens == pytest.approx(lobes[0] - ramp, abs=1e-9)
    assert closes == pytest.approx(lobes[-1] + ramp, abs=1e-9)
    assert closes - opens == pytest.approx(64 * epi.echo_spacing_s, abs=1e-12)


@pytest.mark.parametrize('design', [pytest.param(EVEN, id='even'), pytest.param(ODD, id='odd')])
def test_every_reported_sample_is_where_the_emitted_gradients_put_it(cepi_opts, tmp_path,
                                                                     design) -> None:
    """
    ``k_per_m`` and ``sample_times_s`` against the compiled sequence, at two dwells.

    Two independent integrations of the compiled gradients -- this file's and
    ``calculate_kspacePP`` -- agree with the module in memory.  The written file agrees within its
    six-digit amplitude rounding, and the two integrations of it agree with each other.
    """
    epi = _build(cepi_opts, **design)
    lines = list(range(64)) if epi.samples_per_echo % 2 == 0 else list(range(4, 60))
    block = epi(lines=lines)
    seq = _compiled(block, cepi_opts)
    reported = epi.k_per_m(lines)

    measured, times = _measured_k(seq)
    assert times == pytest.approx(epi.sample_times_s(lines), abs=1e-12)
    assert np.abs(measured - reported).max() < 1e-9 * DK
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        k = sc.kspace(sc.LogicBlock('probe').add(0.0, block), cepi_opts)
    assert np.abs(k['k_adc'][:2] - reported).max() < 1e-9 * DK

    back = _read_back(seq, tmp_path)
    from_file, file_times = _measured_k(back)
    k_pp, t_pp = back.calculate_kspacePP()[:2]
    assert file_times == pytest.approx(epi.sample_times_s(lines), abs=1e-12)
    assert np.asarray(t_pp) == pytest.approx(file_times, abs=1e-12)
    assert np.abs(np.asarray(k_pp)[:2] - from_file).max() < 1e-9 * DK
    error = np.abs(from_file - reported).max(axis=1)
    assert np.all(error <= _six_digit_bound(epi, lines))


def test_flat_top_samples_sit_on_their_line_and_only_blip_samples_leave_it(cepi_opts) -> None:
    """
    A sample is off its ky line exactly when a blip is playing: the set is the blips' support.

    Off-line is measured as ky against the lobe's own midpoint; the blip support is read from the
    compiled blip-axis gradient, not from the design's ``blip_duration_s``.
    """
    epi = _build(cepi_opts)
    lines = list(range(64))
    seq = _compiled(epi(lines=lines), cepi_opts)
    k, times = _measured_k(seq)
    per_echo = epi.samples_per_echo

    ky = k[1].reshape(64, per_echo)
    mid = _integral(_axis_knots(seq, 'y'),
                    np.asarray([epi.time_to_echo(n) for n in range(64)]))
    off_line = np.abs(ky - mid[:, None]) > 1e-9 * DK

    y_times, y_amps = _axis_knots(seq, 'y')
    playing = np.zeros(times.size, bool)
    starts = y_times[:-1][(np.abs(y_amps[:-1]) > 0) | (np.abs(y_amps[1:]) > 0)]
    ends = y_times[1:][(np.abs(y_amps[:-1]) > 0) | (np.abs(y_amps[1:]) > 0)]
    for a, b in zip(starts, ends):
        if b > epi.sample_times_s(lines)[0]:               # the blip-axis prephaser is not a blip
            playing |= (times > a + 1e-12) & (times < b - 1e-12)
    assert np.array_equal(off_line.ravel(), playing)
    # The flat top is on its line, at the k the table asks for.
    on_line = ky[~off_line]
    wanted = np.broadcast_to(((np.asarray(lines) - epi.center_line) * DK)[:, None], ky.shape)
    assert np.abs(on_line - wanted[~off_line]).max() < 1e-9 * DK
    # And a blip moves ky one way only.
    assert np.all(np.diff(k[1]) >= -1e-12)


@pytest.mark.parametrize('design', [pytest.param(EVEN, id='even'), pytest.param(ODD, id='odd')])
def test_forward_and_reverse_lobes_sample_one_kx_grid(cepi_opts, design) -> None:
    """A reverse lobe's samples, reversed, land on the forward lobe's kx, every lobe."""
    epi = _build(cepi_opts, **design)
    lines = list(range(28, 36))
    k, _ = _measured_k(_compiled(epi(lines=lines), cepi_opts))

    kx = k[0].reshape(len(lines), epi.samples_per_echo)
    for n in range(len(lines) - 1):
        assert np.abs(kx[n] - kx[n + 1][::-1]).max() < 1e-9 * DK
    # The sampled extent is the lobe's area less two half-dwell ramp corners.
    corner = epi.readout_amplitude_hz_m / float(epi.gx.rise_time) * epi.dwell_s ** 2 / 8
    assert np.ptp(kx[0]) == pytest.approx(MATRIX[0] * DK - 2 * corner, abs=1e-9 * DK)


@pytest.mark.parametrize('delta', [0.25, 0.5])
def test_a_window_opened_part_of_a_dwell_late_breaks_coincidence_in_proportion(cepi_opts,
                                                                              delta) -> None:
    """
    The control for the test above: shift the module's own window and the grids separate.

    The module cannot open its window late, so its ADC is replaced here by one whose delay is
    ``delta`` of a dwell longer.  Forward and reverse lobes then disagree by ``2 G dwell delta``.
    """
    epi = _build(cepi_opts, **ODD)
    lines = list(range(28, 36))
    block = epi(lines=lines)
    for node in block.nodes:
        if getattr(node.item, 'type', '') == 'adc':
            node.item = sc.events.derive(node.item,
                                         delay=float(node.item.delay) + delta * epi.dwell_s)
    k, _ = _measured_k(_compiled(block, cepi_opts))

    kx = k[0].reshape(len(lines), epi.samples_per_echo)
    error = max(np.abs(kx[n] - kx[n + 1][::-1]).max() for n in range(len(lines) - 1))
    assert error == pytest.approx(2 * epi.readout_amplitude_hz_m * epi.dwell_s * delta,
                                  rel=1e-9)


@pytest.mark.parametrize('design', [pytest.param(EVEN, id='even'), pytest.param(ODD, id='odd')])
def test_each_echo_instant_is_its_lobe_midpoint_inside_the_window(cepi_opts, tmp_path,
                                                                 design) -> None:
    """
    Every lobe's ``kx = 0``, found by bisection on the compiled gradient, is ``time_to_echo``.

    It sits on a sample when the samples per echo are odd and halfway between two when they are
    even, and consecutive instants are exactly one echo spacing apart.  Read back, the instant can
    move by the k rounding over the gradient at the crossing.
    """
    epi = _build(cepi_opts, **design)
    lines = list(range(28, 36))
    seq = _compiled(epi(lines=lines), cepi_opts)
    sample_times = epi.sample_times_s(lines).reshape(len(lines), epi.samples_per_echo)

    def crossing(knots, n) -> float:
        lo = epi.time_to_echo(n) - epi.echo_spacing_s / 4
        hi = epi.time_to_echo(n) + epi.echo_spacing_s / 4
        f_lo = _integral(knots, lo)[0]
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if (_integral(knots, mid)[0] > 0) == (f_lo > 0):
                lo, f_lo = mid, _integral(knots, mid)[0]
            else:
                hi = mid
        return 0.5 * (lo + hi)

    knots = _axis_knots(seq, 'x')
    found = np.array([crossing(knots, n) for n in range(len(lines))])
    assert np.abs(found - [epi.time_to_echo(n) for n in range(len(lines))]).max() < 1e-9
    assert np.diff(found) == pytest.approx([epi.echo_spacing_s] * (len(lines) - 1), abs=1e-9)
    nearest = {round(float(np.abs(sample_times[n] - found[n]).min() / epi.dwell_s), 9)
               for n in range(len(lines))}
    assert nearest == ({0.0} if epi.samples_per_echo % 2 else {0.5})

    back_knots = _axis_knots(_read_back(seq, tmp_path), 'x')
    allowed = _six_digit_bound(epi, lines)[0] / epi.readout_amplitude_hz_m
    back = np.array([crossing(back_knots, n) for n in range(len(lines))])
    assert np.abs(back - found).max() <= allowed


def test_the_origin_crossing_has_ky_zero(cepi_opts) -> None:
    """At ``time_to_center_line`` both kx and ky are zero: the effective echo time is one instant."""
    epi = _build(cepi_opts)
    for table in (list(range(64)), list(range(63, -1, -1)), list(range(20, 50))):
        seq = _compiled(epi(lines=table), cepi_opts)
        at = epi.time_to_center_line(table)
        origin = [_integral(_axis_knots(seq, axis), at)[0] for axis in 'xy']
        assert np.abs(origin).max() < 1e-9 * DK


def test_the_largest_kx_step_is_the_oversampling_bound(cepi_opts) -> None:
    """
    Every kx step between consecutive samples is at most ``dk / oversampling``.

    The largest is on the flat top, ``G * dwell``, and the design meets the bound there.
    """
    for oversampling in (1.0, 1.5, 2.0):
        epi = _build(cepi_opts, oversampling=oversampling)
        k, _ = _measured_k(_compiled(epi(lines=range(30, 34)), cepi_opts))
        steps = np.abs(np.diff(k[0]))
        assert steps.max() <= DK / oversampling * (1 + 1e-9)
        assert steps.max() == pytest.approx(epi.readout_amplitude_hz_m * epi.dwell_s, rel=1e-9)


def test_the_lobe_is_the_shortest_that_meets_every_limit(cepi_opts) -> None:
    """
    No shorter echo spacing on the dwell's lattice admits a legal lobe, checked by brute force.

    For each candidate spacing every ramp on the gradient raster is tried against the four
    limits: amplitude, slew, the sampling bound and half a blip on each ramp.
    """
    epi = _build(cepi_opts)
    assert (epi.samples_per_echo, round(epi.echo_spacing_s * 1e6)) == (118, 590)
    assert (round(float(epi.gx.rise_time) * 1e6), round(epi.readout_amplitude_hz_m)) == (110, 606061)

    area, raster = MATRIX[0] * DK, 10e-6
    lattice = 10e-6                                      # lcm(5 us, 10 us)

    def legal(spacing: float, ramp: float) -> bool:
        g = area / (spacing - ramp)
        return (g <= cepi_opts.max_grad and g / ramp <= cepi_opts.max_slew * (1 + 1e-9)
                and g * epi.dwell_s <= DK / epi.oversampling * (1 + 1e-12)
                and ramp >= epi.blip_duration_s / 2 - 1e-12 and 2 * ramp <= spacing + 1e-12)

    for steps in range(1, round(epi.echo_spacing_s / lattice)):
        spacing = steps * lattice
        assert not any(legal(spacing, r * raster) for r in range(1, steps // 2 + 1)), spacing
    ramps = [r * raster for r in range(1, 30) if legal(epi.echo_spacing_s, r * raster)]
    assert ramps[0] == pytest.approx(float(epi.gx.rise_time))


def test_each_half_blip_lies_inside_a_readout_ramp(cepi_opts) -> None:
    """
    While a blip plays the readout is never on its flat top, including with a lengthened blip.

    The second case is where the ramp has to grow to hold the blip, which lengthens the lobe.
    """
    for blip_s in (None, 300e-6):
        epi = _build(cepi_opts, blip_duration_s=blip_s)
        assert float(epi.gx.rise_time) >= epi.blip_duration_s / 2 - 1e-12
        seq = _compiled(epi(lines=range(20, 28)), cepi_opts)
        y_times, y_amps = _axis_knots(seq, 'y')
        busy = y_times[np.abs(y_amps) > 0]
        after = busy[busy > epi.sample_times_s(range(20, 28))[0]]
        gx_at = np.interp(after, *_axis_knots(seq, 'x'))
        assert np.abs(gx_at).max() < epi.readout_amplitude_hz_m * (1 - 1e-6)
    assert _build(cepi_opts, blip_duration_s=300e-6).echo_spacing_s > _build(
        cepi_opts).echo_spacing_s


# ------------------------------------------------------------------- the block and its timing
def test_sample_times_include_the_adc_delay_and_the_window_ends_inside_its_block(
        cepi_opts, tmp_path) -> None:
    """
    The block reserves the receiver's trailing dead time, so whatever follows it starts after.

    The successor varies -- nothing, a spoiler on the readout axis, a second train -- because a
    reservation that only ever met one neighbour is untested against the others.
    """
    epi = _build(cepi_opts)
    lines = list(range(24, 40))
    train = epi(lines=lines)
    spoil = sc.modules.spoiler(cepi_opts, cycles_per_voxel=4.0, voxel_mm=3.4, axis='x')
    for successor in (None, spoil, epi(lines=lines[::-1])):
        tree = sc.LogicBlock('probe').add(0.0, train)
        if successor is not None:
            tree.add(train.duration, successor)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', sc.SeqCraftWarning)
            back = _read_back(sc.compile(tree, cepi_opts), tmp_path)
        first = _adc_events(back)[0]
        assert _sample_times(back)[:first['samples']] == pytest.approx(
            epi.sample_times_s(lines), abs=1e-12)
        dead_time = float(cepi_opts.adc_dead_time)        # the file does not carry it
        window_end = first['block_start'] + first['delay'] + first['samples'] * first['dwell']
        assert first['block_duration'] >= first['delay'] + first['samples'] * first['dwell'] \
            + dead_time - 1e-12
        assert train.duration == pytest.approx(window_end + dead_time, abs=1e-12)
        if successor is not None:
            assert first['block_start'] + first['block_duration'] <= train.duration + 1e-12


def test_a_dummy_shot_occupies_the_same_span_and_plays_the_same_gradients(cepi_opts) -> None:
    """A dummy loads the gradients exactly as a real shot does, over exactly the same time."""
    epi = _build(cepi_opts, **{**EVEN, 'shots': 2, 'blip_lines': 2})
    table = list(range(1, 64, 2))
    for segment in (None, 1):
        real = epi(lines=table, segment=segment)
        dummy = epi(lines=table, segment=segment, acquire=False)
        assert dummy.duration == pytest.approx(real.duration, abs=1e-15)

        def gradients(block):
            return [(round(t, 12), sc.events.content_hash(e)) for t, e, _ in sc.flatten(block)
                    if e.type in ('trap', 'grad')]

        assert gradients(dummy) == gradients(real)
        assert not [e for _, e, _ in sc.flatten(dummy) if e.type in ('adc', 'labelset')]


def test_prephase_false_opens_the_block_with_the_receiver_dead_time(cepi_opts) -> None:
    """
    No prephaser, and the train starts one dead time in, where the window's lead ends.

    ``area_to_echo_per_m`` and ``k_blip_per_m`` stay available, because they are what a caller
    dephasing before a refocusing pulse has to supply.
    """
    with_prephaser, without = _build(cepi_opts), _build(cepi_opts, prephase=False)
    block = without(lines=range(64))

    adc = next(n for n in block.nodes if getattr(n.item, 'type', '') == 'adc')
    assert adc.start == 0.0
    assert float(adc.item.delay) == pytest.approx(cepi_opts.adc_dead_time)
    assert without.sample_times_s(range(64))[0] == pytest.approx(
        cepi_opts.adc_dead_time + without.dwell_s / 2)
    assert not [n for n in block.nodes if getattr(n.item, 'channel', '') == 'y'
                and n.start < cepi_opts.adc_dead_time]
    assert without.area_to_echo_per_m == pytest.approx(with_prephaser.area_to_echo_per_m)
    assert with_prephaser.time_to_echo(0) - without.time_to_echo(0) == pytest.approx(
        with_prephaser.sample_times_s([0, 1])[0] - without.sample_times_s([0, 1])[0], abs=1e-15)


def test_a_spin_echo_composition_with_prephase_false_lands_on_the_nominal_trajectory(
        cepi_opts) -> None:
    """
    The nominal trajectory is what a correct composition delivers -- checked on one composition.

    The dephasers go before the refocusing pulse with the signs of ``area_to_echo_per_m`` and
    ``-k_blip_per_m``, because the pulse conjugates k.  ``calculate_kspacePP`` models that
    conjugation; with the dephasers dropped, the measurement leaves the nominal trajectory by
    exactly what they would have supplied.
    """
    opts = _scanner(max_grad=40, max_slew=150)
    epi = _build(opts, prephase=False)
    exc = sc.modules.Excitation(opts=opts, flip_deg=90.0, thickness_mm=5.0, duration_s=2e-3)
    refoc = sc.modules.Refocusing(opts=opts, thickness_mm=6.25, duration_s=3e-3,
                                  crush_cycles_per_voxel=4.0, crush_voxel_mm=5.0)
    lines = list(range(64))
    deph_x = pp.make_trapezoid('x', area=epi.area_to_echo_per_m, duration=800e-6, system=opts)
    deph_y = pp.make_trapezoid('y', area=-epi.k_blip_per_m(lines[0]), duration=800e-6,
                               system=opts)

    def acquired(*dephasers) -> np.ndarray:
        start = exc().duration + 800e-6 + refoc().duration
        tree = sc.LogicBlock('se').add(0.0, exc())
        if dephasers:
            tree.add(exc().duration, *dephasers)
        tree.add(exc().duration + 800e-6, refoc()).add(start, epi(lines=lines))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', sc.SeqCraftWarning)
            return sc.kspace(tree, opts)['k_adc'][:2]

    nominal = epi.k_per_m(lines)
    assert np.abs(acquired(deph_x, deph_y) - nominal).max() < 1e-6 * DK
    missing = acquired() - nominal
    assert missing[0] == pytest.approx(np.full(missing.shape[1], epi.area_to_echo_per_m),
                                       abs=1e-6 * DK)
    assert missing[1] == pytest.approx(np.full(missing.shape[1], -epi.k_blip_per_m(0)),
                                       abs=1e-6 * DK)


def test_shot_s_waits_exactly_s_esp_over_shots_and_ky_time_is_affine(cepi_opts) -> None:
    """
    The interleave delay is a whole number of block rasters, and t(ky) is a straight line.

    With the echo at each lobe's midpoint the line is exact: no alternating half-dwell residual.
    """
    shots = 4
    epi = _build(cepi_opts, blip_lines=shots, shots=shots)
    raster = float(cepi_opts.block_duration_raster)
    assert epi.shot_delay_s == pytest.approx(epi.echo_spacing_s / shots, abs=1e-15)

    times = np.empty(MATRIX[1])
    spans = set()
    for segment in range(shots):
        offset = epi.shot_offset_s(segment)
        assert offset / raster == pytest.approx(round(offset / raster), abs=1e-9)
        table = list(range(segment, MATRIX[1], shots))
        for echo, line in enumerate(table):
            times[line] = epi.time_to_echo(echo, segment=segment)
        spans.add(round(epi(lines=table, segment=segment).duration, 12))
    ky = np.arange(MATRIX[1])
    residual = times - np.polyval(np.polyfit(ky, times, 1), ky)
    assert np.abs(residual).max() < 1e-12
    assert np.diff(times) == pytest.approx([epi.shot_delay_s] * (MATRIX[1] - 1), abs=1e-12)
    assert len(spans) == 1


def test_the_receiver_phase_reaches_the_adc_and_moves_nothing_else(cepi_opts) -> None:
    epi = _build(cepi_opts)
    plain, spoiled = epi(lines=range(8)), epi(lines=range(8), phase_deg=117.0)

    adcs = [e for _, e, _ in sc.flatten(spoiled) if e.type == 'adc']
    assert [float(a.phase_offset) for a in adcs] == pytest.approx([np.deg2rad(117.0)])
    others = [[(round(t, 12), sc.events.content_hash(e)) for t, e, _ in sc.flatten(b)
               if e.type != 'adc'] for b in (plain, spoiled)]
    assert others[0] == others[1]


def test_the_only_label_is_seg_and_only_when_asked(cepi_opts) -> None:
    """One ADC event cannot carry a line or a polarity per lobe, so the file carries neither."""
    epi = _build(cepi_opts)
    labels = [(e.label, int(e.value)) for _, e, _ in sc.flatten(epi(lines=range(8), segment=3))
              if e.type == 'labelset']
    assert labels == [('SEG', 3)]
    assert not [e for _, e, _ in sc.flatten(epi(lines=range(8))) if e.type == 'labelset']


# --------------------------------------------------------------------------- the compiler
def test_the_train_compiles_on_its_own_as_one_adc_block(cepi_opts) -> None:
    """
    Two blocks: the prephasers, then the window with each axis's train summed into one gradient.

    The summing is the compiler's, and it reports it once per axis.
    """
    epi = _build(cepi_opts)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        seq = sc.compile(sc.LogicBlock('probe').add(0.0, epi(lines=range(64))), cepi_opts)

    assert len(seq.block_events) == 2
    window = seq.get_block(2)
    assert int(window.adc.num_samples) == 64 * epi.samples_per_echo
    assert (window.gx.type, window.gy.type) == ('grad', 'grad')
    merges = [str(w.message) for w in caught if 'gradient merge' in str(w.message)]
    assert len(merges) == 1
    assert '2 same-axis gradient merges' in merges[0]
    assert '(axis x)' in merges[0] and '(axis y)' in merges[0]


@pytest.mark.parametrize('kwargs', [
    pytest.param({}, id='reference'),
    pytest.param(ODD, id='odd'),
    pytest.param({'blip_lines': 4, 'shots': 4}, id='four-shot'),
    pytest.param({'blip_duration_s': 300e-6}, id='long-blip'),
    pytest.param({'prephase': False}, id='no-prephaser'),
])
def test_the_emitted_gradients_stay_within_per_axis_limits(cepi_opts, kwargs) -> None:
    """Per-axis amplitude and slew on the compiled knots; the vector norm is not a limit here."""
    epi = _build(cepi_opts, **kwargs)
    table = list(range(2, 62, kwargs.get('blip_lines', 1)))
    seq = _compiled(epi(lines=table, segment=0 if kwargs.get('shots') else None), cepi_opts)
    for axis in 'xy':
        times, amps = _axis_knots(seq, axis)
        assert np.abs(amps).max() <= cepi_opts.max_grad * (1 + 1e-9)
        span = np.diff(times)
        slews = np.abs(np.diff(amps)[span > 1e-12] / span[span > 1e-12])
        assert slews.max() <= cepi_opts.max_slew * (1 + 1e-6)


def test_two_calls_produce_identical_events(cepi_opts, component_checks) -> None:
    """Purity, determinism and compiling alone, for a reversed and an interleaved table."""
    component_checks.all(_build(cepi_opts), lines=list(range(63, -1, -1)))
    component_checks.all(_build(cepi_opts, blip_lines=4, shots=4),
                         lines=list(range(2, 64, 4)), segment=2)


# ------------------------------------------------------------------- the event and its limits
def test_a_train_at_the_event_limit_builds_and_one_sample_over_is_refused() -> None:
    """
    At a fixed design and table, a limit of exactly the train's count builds one event of it.

    One sample less refuses at build, naming the count and the limit, and offers no split.  At
    8192 the longest table that fits builds and the next refuses.
    """
    design = {'dwell_s': 3.5e-6, 'oversampling': 1.0}                  # 140 samples per echo
    lines = list(range(24, 40))
    count = len(lines) * 140

    exact = _build(_scanner(adc_samples_limit=count), **design)
    assert exact.samples_per_echo == 140
    seq = _compiled(exact(lines=lines), exact.opts)
    assert [e['samples'] for e in _adc_events(seq)] == [count]

    with pytest.raises(sc.ConfigurationError) as caught:
        _build(_scanner(adc_samples_limit=count - 1), **design)(lines=lines)
    message = str(caught.value)
    assert message.splitlines()[0] == (
        f"{count} ADC samples in this train's one ADC event (16 lines x 140 per echo), above "
        f'the {count - 1} that opts.adc_samples_limit allows.')
    assert 'split it into' not in message
    assert 'never splits its window' in message

    at_8192 = _build(_scanner(), **design)
    fits = 8192 // 140
    assert [e['samples'] for e in _adc_events(_compiled(at_8192(lines=range(fits)),
                                                        at_8192.opts))] == [fits * 140]
    with pytest.raises(sc.ConfigurationError, match='above the 8192'):
        at_8192(lines=range(fits + 1))


def test_the_limit_refusal_offers_fixes_that_build(cepi_opts) -> None:
    """The largest table it names builds, and so does the oversampling it names."""
    epi = _build(_scanner(adc_samples_limit=4000), oversampling=2.0)
    with pytest.raises(sc.ConfigurationError) as caught:
        epi(lines=range(64))
    message = str(caught.value)

    longest = int(message.split('a table of at most ')[1].split()[0])
    _compiled(epi(lines=range(longest)), epi.opts)
    for longer in (longest + 1, longest + 2):       # by the divisor, then by the limit
        with pytest.raises(sc.ConfigurationError):
            epi(lines=range(longer))
    offered = message.split('oversampling=1 gives ')[1]
    per_echo = int(offered.split()[0])
    builds = int(offered.split('a table of up to ')[1].split()[0])
    lower = _build(_scanner(adc_samples_limit=4000), oversampling=1.0)
    assert lower.samples_per_echo == per_echo
    assert builds > longest
    _compiled(lower(lines=range(builds)), lower.opts)


def test_no_configured_limit_is_refused_at_construction() -> None:
    """``adc_samples_limit = 0`` is PyPulseq's "no limit", and there is then nothing to fit."""
    for opts in (_scanner(adc_samples_limit=0),
                 Opts(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s')):
        with pytest.raises(sc.ConfigurationError, match='adc_samples_limit is 0'):
            _build(opts)


def test_a_count_the_divisor_does_not_divide_is_refused_by_name(cepi_opts) -> None:
    """
    118 samples per echo over 63 lines is 7434, which 4 does not divide.  62 and 64 lines do.

    The tables and the dwell the refusal offers are each built and compiled, and each says what
    it changes about the acquisition.
    """
    epi = _build(cepi_opts)
    assert epi.samples_per_echo % 4 == 2
    with pytest.raises(sc.ConfigurationError) as caught:
        epi(lines=range(63))
    message = str(caught.value)
    assert message.splitlines()[0] == (
        '7434 ADC samples (63 lines x 118 per echo) is not a multiple of '
        "opts.adc_samples_divisor = 4, which PyPulseq's timing check requires of every ADC event.")
    assert 'a table of 62 or 64 lines divides' in message
    assert 'a different set of lines' in message
    for length in (62, 64):
        _compiled(epi(lines=range(length)), cepi_opts)

    dwell = float(message.split('dwell_s=')[1].split(':')[0])
    other = _build(cepi_opts, dwell_s=dwell)
    assert other.samples_per_echo % 4 == 0
    assert f'{other.echo_spacing_s * 1e6:.1f} us' in message
    _compiled(other(lines=range(63)), cepi_opts)


def test_a_divisor_of_one_constrains_nothing() -> None:
    epi = _build(_scanner(adc_samples_divisor=1))
    _compiled(epi(lines=range(63)), epi.opts)


# ------------------------------------------------------------------ shared with EPI2D, and its own
@pytest.mark.parametrize(('lines', 'kwargs'), [
    pytest.param([], {}, id='empty'),
    pytest.param([0, 64], {}, id='outside'),
    pytest.param([3, 3], {}, id='repeated'),
    pytest.param([0, 4], {'blip_lines': 2}, id='step-beyond-blip'),
])
def test_both_epi_readouts_refuse_the_same_tables_for_the_same_reasons(cepi_opts, lines,
                                                                       kwargs) -> None:
    """One owner for the table rules: the same table is refused by both, with the same words."""
    epi2d = sc.modules.EPI2D(opts=cepi_opts, fov_mm=FOV_MM, matrix=MATRIX, dwell_s=3.5e-6,
                             **kwargs)
    texts = []
    for readout in (epi2d, _build(cepi_opts, **kwargs)):
        with pytest.raises(sc.ConfigurationError) as caught:
            readout(lines=lines)
        texts.append(str(caught.value).splitlines())
    assert texts[0][0] == texts[1][0]
    if 'blip_lines' in kwargs:
        assert 'the readout ramps that must each hold half of it' in '\n'.join(texts[1])
    else:
        assert texts[0] == texts[1]


@pytest.mark.parametrize(('kwargs', 'match'), [
    ({'dwell_s': 4.05e-6}, 'not a whole number of adc_raster_time'),
    ({'oversampling': 0.9}, 'undersamples the field of view'),
    ({'blip_axis': 'x'}, 'both'),
    ({'blip_lines': 0}, 'at least 1'),
    ({'shots': 0}, 'at least 1'),
    ({'matrix': (1, 64)}, 'at least 2'),
    ({'blip_duration_s': 1e-6}, 'shorter than the minimum'),
])
def test_the_construction_refusals(cepi_opts, kwargs, match) -> None:
    """Each one refused once, so that a message that stops naming its fix is noticed."""
    with pytest.raises(sc.ConfigurationError, match=match):
        _build(cepi_opts, **kwargs)


def test_a_dwell_off_the_adc_raster_names_the_two_nearest_and_both_build(cepi_opts) -> None:
    with pytest.raises(sc.ConfigurationError) as caught:
        _build(cepi_opts, dwell_s=4.05e-6)
    named = str(caught.value).split('the nearest legal dwells are ')[1].split(' and ')
    for dwell in named:
        assert _build(cepi_opts, dwell_s=float(dwell)).samples_per_echo > 0


def test_every_query_refuses_what_build_refuses(cepi_opts) -> None:
    """A table that cannot be played has no sample times, no trajectory and no echo time."""
    epi = _build(cepi_opts)
    for query in (epi.sample_times_s, epi.k_per_m, epi.time_to_center_line):
        with pytest.raises(sc.ConfigurationError, match='not a multiple'):
            query(range(63))
        with pytest.raises(sc.ConfigurationError, match='more than once'):
            query([3, 3])


def test_a_table_without_the_centre_has_no_echo_time(cepi_opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='omits the centre'):
        _build(cepi_opts).time_to_center_line(range(0, 20))


def test_a_segment_outside_the_interleave_is_refused(cepi_opts) -> None:
    epi = _build(cepi_opts, blip_lines=4, shots=4)
    with pytest.raises(sc.ConfigurationError, match='not one of the 4 interleaved shots'):
        epi(lines=range(0, 64, 4), segment=4)
    assert epi.shot_offset_s(None) == 0.0
    with pytest.raises(sc.ConfigurationError, match='before the start of the train'):
        epi.polarity(-1)


# ----------------------------------------------------------------- what moves with the table
def test_reversing_the_table_reverses_ky_and_leaves_kx(cepi_opts) -> None:
    """kx depends on the lobe index only; ky follows the table; the echo time moves with it."""
    epi = _build(cepi_opts)
    forward, reverse = list(range(64)), list(range(63, -1, -1))
    k_f, k_r = epi.k_per_m(forward), epi.k_per_m(reverse)
    per_echo = epi.samples_per_echo

    assert np.array_equal(k_f[0], k_r[0])
    mid = [n * per_echo + per_echo // 2 for n in range(64)]
    assert k_r[1][mid] == pytest.approx(k_f[1][mid][::-1], abs=1e-9 * DK)
    assert epi.time_to_center_line(reverse) - epi.time_to_center_line(forward) == pytest.approx(
        -epi.echo_spacing_s, abs=1e-12)


def test_an_interleaved_table_steps_ky_by_its_blip(cepi_opts) -> None:
    """With ``blip_lines = 2`` every mid-lobe step is exactly two lines."""
    epi = _build(cepi_opts, blip_lines=2)
    table = list(range(1, 64, 2))
    k = epi.k_per_m(table)
    mid = k[1][[n * epi.samples_per_echo + epi.samples_per_echo // 2 for n in range(len(table))]]
    assert np.diff(mid) == pytest.approx([2 * DK] * (len(table) - 1), abs=1e-9 * DK)
