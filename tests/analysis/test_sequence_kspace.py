"""
``sc.sequence_kspace`` against an independent exact integration of the emitted gradients.

For an ideal sequence the nominal k-space position at an ADC sample has one value: the integral of
the piecewise-linear gradient from the last excitation's RF centre to the sample time, conjugated
at each refocusing pulse.  ``Exact`` below integrates the block events' own knots with NumPy,
independently of PyPulseq's k-space code, and ``sequence_kspace`` has to match it at the ADC
samples, in memory and after a write and read.

The sequences are built with PyPulseq directly, so that each gradient shape is the one named:
ramps, trapezoids, arbitrary gradients, Wave readouts with odd and even cycle counts, and a gradient
running through the RF centre.  Each is sampled on the raster edges, half a raster off them, and at
an incommensurate dwell.  The SeqCraft-built case checks that ``sc.kspace`` of a tree and
``sequence_kspace`` of its written file report one trajectory.
"""

import numpy as np
import pypulseq as pp
import pytest
from pypulseq.opts import Opts

import seqcraft as sc

SYSTEM = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
                 rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=0.0)
D = SYSTEM.grad_raster_time
#: The backend pads waveform edges by 1 ps; measured agreement is <= 3.4e-7 1/m.
TOLERANCE = 1e-6


# ------------------------------------------------------------------ independent exact integration
def _knots(event, t0):
    """The piecewise-linear knots of a gradient event, in absolute time."""
    if event.type == 'trap':
        t = np.cumsum([0.0, event.rise_time, event.flat_time, event.fall_time]) + event.delay + t0
        return t, np.array([0.0, event.amplitude, event.amplitude, 0.0])
    tt = np.asarray(event.tt, dtype=float)
    waveform = np.asarray(event.waveform, dtype=float)
    if tt[0] > 1e-12:  # samples at raster centres: the edge values live in first/last
        tt = np.concatenate([[0.0], tt, [event.shape_dur]])
        waveform = np.concatenate([[event.first], waveform, [event.last]])
    return tt + event.delay + t0, waveform


class Exact:
    """Exact k at arbitrary times, from the block events of `seq` alone."""

    def __init__(self, seq):
        pieces = {axis: [] for axis in 'xyz'}
        self.t_adc, self.excitations, self.refocusings = [], [], []
        t0 = 0.0
        for index in seq.block_events:
            block = seq.get_block(index)
            if block.rf is not None:
                centre = t0 + block.rf.delay + pp.calc_rf_center(block.rf)[0]
                (self.refocusings if block.rf.use == 'refocusing' else self.excitations).append(centre)
            for axis in 'xyz':
                event = getattr(block, 'g' + axis, None)
                if event is not None:
                    pieces[axis].append(_knots(event, t0))
            if block.adc is not None:
                adc = block.adc
                self.t_adc.append(t0 + adc.delay + (np.arange(adc.num_samples) + 0.5) * adc.dwell)
            t0 += seq.block_durations[index]
        self.t_adc = np.concatenate(self.t_adc) if self.t_adc else np.zeros(0)
        self.knots = {}
        for axis in 'xyz':
            t = np.concatenate([[0.0]] + [p[0] for p in pieces[axis]] + [[t0]])
            g = np.concatenate([[0.0]] + [p[1] for p in pieces[axis]] + [[0.0]])
            keep = np.concatenate([[True], np.diff(t) > 1e-15])
            t, g = t[keep], g[keep]
            area = np.concatenate([[0.0], np.cumsum(np.diff(t) * (g[:-1] + g[1:]) / 2)])
            self.knots[axis] = (t, g, area)

    def moment(self, axis, t):
        """Exact integral of the gradient on `axis` from 0 to `t` (constant after the end)."""
        tk, gk, area = self.knots[axis]
        t = min(max(t, tk[0]), tk[-1])
        i = int(np.clip(np.searchsorted(tk, t, side='right') - 1, 0, tk.size - 2))
        tau, width = t - tk[i], tk[i + 1] - tk[i]
        return area[i] + gk[i] * tau + (gk[i + 1] - gk[i]) * tau**2 / (2 * width)

    def gradient(self, axis, t):
        tk, gk, _ = self.knots[axis]
        return np.interp(t, tk, gk)

    def k(self, times=None):
        """``(3, n)`` k at `times` (the ADC samples by default)."""
        times = self.t_adc if times is None else np.asarray(times, dtype=float)
        out = np.zeros((3, times.size))
        for n, t in enumerate(times):
            start = max([e for e in self.excitations if e <= t], default=0.0)
            flips = sorted(r for r in self.refocusings if start < r <= t)
            for row, axis in enumerate('xyz'):
                k, previous = 0.0, start
                for r in flips:
                    k = -(k + self.moment(axis, r) - self.moment(axis, previous))
                    previous = r
                out[row, n] = k + self.moment(axis, t) - self.moment(axis, previous)
        return out


# ------------------------------------------------------------------ the sequences
SAMPLINGS = {  # (dwell, ADC delay) placing the samples on raster edges, half a raster off, or neither
    'on raster edges': (10e-6, 5e-6),
    'half a raster off': (10e-6, 0.0),
    'incommensurate dwell': (6.5e-6, 0.0),
}


def _adc(sampling, window, extra_delay=0.0):
    dwell, delay = SAMPLINGS[sampling]
    return pp.make_adc(num_samples=int((window - delay) // dwell), dwell=dwell,
                       delay=extra_delay + delay, system=SYSTEM)


def ramp(sampling):
    """No RF: one ramp and plateau on x, sampled throughout."""
    seq = pp.Sequence(system=SYSTEM)
    g = 20e-3 * SYSTEM.gamma
    gx = pp.make_extended_trapezoid('x', times=np.array([0, 1e-3, 2e-3, 2.2e-3]),
                                    amplitudes=np.array([0, g, g, 0]), system=SYSTEM)
    seq.add_block(gx, _adc(sampling, 2e-3))
    return seq


def _excitation(thickness=5e-3):
    return pp.make_sinc_pulse(flip_angle=np.deg2rad(15), duration=1e-3, slice_thickness=thickness,
                              apodization=0.5, time_bw_product=4, delay=SYSTEM.rf_dead_time,
                              system=SYSTEM, return_gz=True, use='excitation')


def trapezoid(sampling):
    """A slice-selective excitation (gradient through the RF centre), then a trapezoid readout."""
    seq = pp.Sequence(system=SYSTEM)
    rf, gz, gzr = _excitation()
    dwell = SAMPLINGS[sampling][0]
    flat = np.ceil(round(128 * dwell / D, 6)) * D + D
    gx = pp.make_trapezoid('x', flat_area=128 / 0.256, flat_time=flat, system=SYSTEM)
    adc = _adc(sampling, 128 * dwell + SAMPLINGS[sampling][1], extra_delay=gx.rise_time)
    seq.add_block(rf, gz)
    seq.add_block(gzr, pp.make_trapezoid('x', area=-gx.area / 2, duration=1e-3, system=SYSTEM))
    seq.add_block(gx, adc)
    return seq


def arbitrary(sampling):
    """Arbitrary piecewise-linear gradients: irregular on-raster knots, and raster-centre samples."""
    seq = pp.Sequence(system=SYSTEM)
    g = 10e-3 * SYSTEM.gamma
    times = np.array([0, 30, 110, 140, 260, 300, 420, 500]) * 1e-5
    gx = pp.make_extended_trapezoid(
        'x', times=times, amplitudes=np.array([0, 0.3, 0.9, -0.4, 0.6, 0.6, -0.2, 0]) * g,
        system=SYSTEM)
    # A raster-centre shape whose edge values PyPulseq's reader restores to the piecewise-linear
    # interpolation of its samples (`restore_additional_shape_samples`), so its knots are known.
    t = (np.arange(200) + 0.5) * D
    shape = np.sin(np.pi * t / t[-1] * 0.999) ** 2 * np.sin(2 * np.pi * 2.3e3 * t)
    gy = pp.make_arbitrary_grad('y', 8e-3 * SYSTEM.gamma * shape, first=0, last=0, system=SYSTEM)
    seq.add_block(gx, gy, _adc(sampling, 5e-3))
    return seq


def wave(sampling, cycles):
    """A slab-selective excitation, then a Wave readout: sine on y, cosine on z, knots on edges."""
    seq = pp.Sequence(system=SYSTEM)
    rf, gz, gzr = _excitation(thickness=0.12)
    window, g0, ramp_time = 2.56e-3, 3.5e-3 * SYSTEM.gamma, 4e-5
    edges = np.arange(0, round(window / D) + 1) * D
    phase = 2 * np.pi * cycles * edges / window
    sine = g0 * np.sin(phase)
    sine[[0, -1]] = 0.0  # whole cycles: the sine starts and ends at rest
    gy = pp.make_extended_trapezoid('y', times=ramp_time + edges, amplitudes=sine, system=SYSTEM)
    gz_wave = pp.make_extended_trapezoid(
        'z', times=np.concatenate([[0.0], ramp_time + edges, [2 * ramp_time + window]]),
        amplitudes=np.concatenate([[0.0], g0 * np.cos(phase), [0.0]]), system=SYSTEM)
    gx = pp.make_trapezoid('x', flat_area=32 / 0.2, flat_time=window, system=SYSTEM)
    gx.delay = ramp_time - gx.rise_time
    adc = _adc(sampling, window - D, extra_delay=ramp_time)
    seq.add_block(rf, gz)
    seq.add_block(gzr)
    seq.add_block(gx, gy, gz_wave, adc)
    return seq


def spin_echo(sampling):
    """Excitation, a lobe, a refocusing pulse, the same lobe again: k conjugates at the pulse."""
    seq = pp.Sequence(system=SYSTEM)
    ninety = pp.make_block_pulse(flip_angle=np.pi / 2, duration=1e-3, delay=1e-4, use='excitation',
                                 system=SYSTEM)
    one_eighty = pp.make_block_pulse(flip_angle=np.pi, duration=2e-3, delay=1e-4, use='refocusing',
                                     system=SYSTEM)
    seq.add_block(ninety)
    seq.add_block(pp.make_trapezoid('x', area=300.0, duration=1e-3, system=SYSTEM))
    seq.add_block(one_eighty)
    seq.add_block(pp.make_trapezoid('x', area=300.0, duration=1e-3, system=SYSTEM),
                  _adc(sampling, 1e-3))
    return seq


CASES = {
    **{f'ramp, {s}': (lambda s=s: ramp(s)) for s in SAMPLINGS},
    **{f'trapezoid after excitation, {s}': (lambda s=s: trapezoid(s)) for s in SAMPLINGS},
    **{f'arbitrary, {s}': (lambda s=s: arbitrary(s)) for s in SAMPLINGS},
    **{f'Wave {c} cycles, {s}': (lambda s=s, c=c: wave(s, c)) for s in SAMPLINGS for c in (5, 6)},
    **{f'spin echo, {s}': (lambda s=s: spin_echo(s)) for s in SAMPLINGS},
}


def _read_back(seq, path):
    seq.write(str(path))
    back = pp.Sequence(system=SYSTEM)
    back.read(str(path))
    return back


# ------------------------------------------------------------------ the contract
# The raster-centre shape in `arbitrary` is read back without edge restoration, which PyPulseq
# announces; the knots it then uses are the piecewise-linear ones `_knots` integrates.
@pytest.mark.filterwarnings('ignore:.*Last restored point differs too much:UserWarning')
@pytest.mark.parametrize('read_back', [False, True], ids=['in memory', 'read back'])
@pytest.mark.parametrize('case', list(CASES))
def test_sequence_kspace_is_the_exact_trajectory_at_the_adc_samples(case, read_back, tmp_path):
    seq = CASES[case]()
    if read_back:
        seq = _read_back(seq, tmp_path / 'case.seq')
    exact = Exact(seq)
    trajectory = sc.sequence_kspace(seq)
    assert trajectory['t_adc'] == pytest.approx(exact.t_adc, abs=1e-12)
    assert np.max(np.abs(trajectory['k_adc'] - exact.k())) < TOLERANCE


@pytest.mark.parametrize('case', ['trapezoid after excitation, incommensurate dwell',
                                  'Wave 5 cycles, half a raster off',
                                  'spin echo, on raster edges'])
def test_the_dense_trajectory_is_exact_at_its_own_times(case):
    """
    ``k`` on ``t_k``, away from the RF centres, where k restarts or conjugates.

    The dense trajectory marks each restart with NaN just before its excitation, so that a plot
    breaks there rather than drawing a line back to the origin.
    """
    seq = CASES[case]()
    exact = Exact(seq)
    trajectory = sc.sequence_kspace(seq)
    k, t_k = trajectory['k'], trajectory['t_k']
    assert k.shape == (3, t_k.size)
    rf = np.concatenate([exact.excitations, exact.refocusings])
    marked = np.isnan(k).any(axis=0)
    restarts = t_k[marked][:, None] - np.asarray(exact.excitations)[None, :]
    assert np.all(np.any((restarts <= 0) & (restarts > -D), axis=1))
    away = ~marked & (np.min(np.abs(t_k[:, None] - rf[None, :]), axis=1) > 1e-9)
    assert np.max(np.abs(k[:, away] - exact.k(t_k[away]))) < TOLERANCE


def test_the_rf_centres_are_named_for_what_they_do():
    seq = CASES['spin echo, on raster edges']()
    exact = Exact(seq)
    trajectory = sc.sequence_kspace(seq)
    assert trajectory['t_excitation'] == pytest.approx(exact.excitations, abs=1e-12)
    assert trajectory['t_refocusing'] == pytest.approx(exact.refocusings, abs=1e-12)
    assert trajectory['k_adc'].shape == (3, trajectory['t_adc'].size)


# ------------------------------------------------------------------ one trajectory, every route
OPTS = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s', B0=3.0,
            rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)


@pytest.mark.parametrize('dwell_s', [10e-6, 6.5e-6], ids=['raster-commensurate', 'incommensurate'])
def test_a_tree_and_its_written_file_report_one_trajectory(dwell_s, tmp_path):
    """
    ``sc.kspace`` of the tree and ``sc.sequence_kspace`` of the file it writes agree.

    A slice-selective excitation puts a gradient through the RF centre; the echo sample of the
    readout sits at kx = 0, and k on the slice axis carries no restart constant.
    """
    excitation = sc.modules.Excitation(opts=OPTS, flip_deg=15.0, thickness_mm=5.0, axis='z',
                                       duration_s=1e-3)
    line = sc.modules.CartesianLine(opts=OPTS, fov_mm=256.0, matrix=64, dwell_s=dwell_s, axis='x')
    first = excitation()
    tree = sc.LogicBlock('probe').add(0.0, first).add(first.duration, line())

    from_tree = sc.kspace(tree, OPTS)
    back = pp.Sequence(system=OPTS)
    path = tmp_path / 'probe.seq'
    sc.compile(tree, OPTS).write(str(path))
    back.read(str(path))
    from_file = sc.sequence_kspace(back)

    assert from_file['t_adc'] == pytest.approx(from_tree['t_adc'], abs=1e-12)
    assert np.max(np.abs(from_file['k_adc'] - from_tree['k_adc'])) < 1e-3 * line.dk_per_m
    assert np.max(np.abs(from_tree['k_adc'] - Exact(back).k())) < 1e-3 * line.dk_per_m
    echo = from_tree['k_adc'][:, line.echo_sample(0)]
    assert abs(echo[0]) < 1e-6 * line.dk_per_m
    assert abs(echo[2]) < 1e-6


# ------------------------------------------------------------------ why not calculate_kspace()
def test_the_pinned_calculate_kspace_is_not_the_exact_trajectory():
    """
    Why SeqCraft does not use the pinned PyPulseq's ``Sequence.calculate_kspace()``.

    It samples the gradient on raster edges and sums rectangles, so on a ramp it reads ahead by
    half a raster of gradient (and has no value before the first raster edge).  A fixed ``trajectory_delay = -D/2`` does not repair it: the
    rectangle sum restarted at the excitation keeps a constant from a slice gradient running
    through the RF centre.  If this fails, the pinned dependency changed;
    ``docs/pypulseq_dependency.md`` says what to re-check.
    """
    seq = CASES['ramp, half a raster off']()
    exact = Exact(seq)
    bias = seq.calculate_kspace()[0][0] - exact.k()[0]
    assert np.nanmax(np.abs(bias - D / 2 * exact.gradient('x', exact.t_adc))) < 0.02
    assert np.nanmax(np.abs(bias)) > 4.0

    seq = CASES['trapezoid after excitation, on raster edges']()
    exact = Exact(seq)
    k_z = seq.calculate_kspace(trajectory_delay=-D / 2)[0][2]
    slice_gradient = float(exact.gradient('z', exact.excitations[0]))
    assert np.nanmax(np.abs(k_z - exact.k()[2])) == pytest.approx(D / 2 * slice_gradient, rel=1e-3)
