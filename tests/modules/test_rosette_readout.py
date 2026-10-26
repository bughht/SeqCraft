"""
``RosetteReadout`` -- the emitted trajectory against the analytic curve, and the report against
the emitted trajectory.

Every measurement here is made by code in this file, which shares nothing with the module: the
knots of the events ``build`` returns, of the compiled sequence in memory, or of the written and
read-back file are integrated here, or read through ``calculate_kspacePP``.  The analytic curve is
written here from the formula, from the readout's public parameters.

The raster departure is checked in three stages, each against its own allowance:

    designed knots (the events build returns)   <= B
    registered (the compiled sequence)          <= B + registration allowance
    read back (the written file)                <= B + registration + serialisation allowance

B is the closed-form bound.  The two allowances are derived from what PyPulseq does to a gradient
shape: registration rounds a compressible shape's peak-normalised samples to 1e-7, and the file
stores the amplitude to six significant digits and the normalised shape to nine.
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc

GOLDEN = (math.sqrt(5.0) - 1.0) / 2.0
TICK = 1e-12


def scanner(grad: float = 32.0, slew: float = 120.0, *, raster: float = 10e-6,
            limit: int = 8192) -> pp.Opts:
    """A scanner with a configured per-event sample limit and divisor."""
    return pp.Opts(max_grad=grad, grad_unit='mT/m', max_slew=slew, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
                   grad_raster_time=raster, block_duration_raster=raster,
                   adc_samples_limit=limit, adc_samples_divisor=4)


def rosette(opts: pp.Opts, **overrides: object) -> sc.modules.RosetteReadout:
    """256 mm, matrix 128, 13 petals at ratio 7/13 and a 2 us dwell, with `overrides` applied."""
    kwargs: dict = {'fov_mm': 256.0, 'matrix': 128, 'petals': 13, 'omega_ratio': 7 / 13,
                    'dwell_s': 2e-6}
    kwargs.update(overrides)
    return sc.modules.RosetteReadout(opts=opts, **kwargs)


#: Pairwise over ratio, petals, dwell, geometry, hardware, raster and angle.  ``record-7-810us``
#: is a window that half a dwell of padding would leave uncovered; ``slow-petal`` puts the first
#: sample exactly on the first crossing, and ``resting-34us`` samples after the gradient has ended.
#: Hardware is (mT/m, T/m/s, gradient raster).
CASES = {
    'record-13': ((32, 120, 10e-6), {}, 0.7),
    'spike-1': ((20, 80, 10e-6), {'omega_ratio': 0.0, 'petals': 1}, 0.0),
    'narrow-2-20us': ((32, 120, 20e-6),
                      {'fov_mm': 220.0, 'matrix': 64, 'omega_ratio': 1 / 13, 'petals': 2,
                       'dwell_s': 6.7e-6}, 0.7),
    'golden-5': ((80, 200, 10e-6), {'fov_mm': 400.0, 'omega_ratio': GOLDEN, 'petals': 5}, 0.0),
    'circles-7': ((32, 120, 10e-6), {'omega_ratio': 1.0, 'petals': 7, 'dwell_s': 4e-6}, 0.7),
    'wide-1-200ns': ((80, 200, 10e-6),
                     {'fov_mm': 220.0, 'matrix': 64, 'omega_ratio': 1.5, 'petals': 1,
                      'dwell_s': 200e-9}, 0.0),
    'wide-13-20us': ((20, 80, 20e-6),
                     {'fov_mm': 400.0, 'omega_ratio': 1.5, 'petals': 13, 'dwell_s': 5.6e-6}, 0.7),
    'slow-petal': ((32, 120, 10e-6),
                   {'fov_mm': 220.0, 'matrix': 64, 'petals': 1, 'dwell_s': 40e-6,
                    'petal_duration_s': 10e-3}, 0.0),
    'record-7-810us': ((32, 120, 10e-6), {'petals': 7, 'petal_duration_s': 810e-6}, 0.0),
    'resting-34us': ((32, 120, 10e-6),
                     {'fov_mm': 220.0, 'matrix': 64, 'petals': 1, 'dwell_s': 34e-6,
                      'petal_duration_s': 10e-3}, 0.7),
}


def case(name: str) -> tuple[pp.Opts, sc.modules.RosetteReadout, float]:
    """The scanner, the readout and the angle of one named case."""
    (grad, slew, raster), overrides, angle = CASES[name]
    opts = scanner(grad, slew, raster=raster)
    return opts, rosette(opts, **overrides), angle


# ---------------------------------------------------------------------------- the oracle
def departure_bound(k_max: float, w1: float, w2: float, dt: float, core: float) -> float:
    """B: Euler--Maclaurin at the raster points, growth within a raster, the remainder."""
    def m(n: int) -> float:
        return 0.5 * k_max * ((w1 + w2) ** n + abs(w1 - w2) ** n)
    return dt**2 / 6 * m(2) + dt**3 / 12 * m(3) + core * dt**4 / 720 * m(5)


class Analytic:
    """The curve a readout claims, written from the formula and its public parameters."""

    def __init__(self, ro: sc.modules.RosetteReadout, opts: pp.Opts) -> None:
        fov_m = ro.fov_mm / 1e3
        self.k_max, self.dk = ro.matrix / (2.0 * fov_m), 1.0 / fov_m
        self.w1 = math.pi / ro.petal_duration_s
        self.w2 = ro.omega_ratio * self.w1
        self.petal, self.petals = ro.petal_duration_s, ro.petals
        self.T = ro.petals * ro.petal_duration_s
        self.dt = float(opts.grad_raster_time)
        self.rasters = round(self.T / self.dt)

    def k(self, tau: np.ndarray) -> np.ndarray:
        return self.k_max * np.sin(self.w1 * tau) * np.exp(1j * self.w2 * tau)

    def g(self, tau: np.ndarray) -> np.ndarray:
        w1, w2 = self.w1, self.w2
        return self.k_max * (w1 * np.cos(w1 * tau) + 1j * w2 * np.sin(w1 * tau)) * np.exp(1j * w2 * tau)

    def m(self, n: int) -> float:
        """A bound on the n-th derivative of k: k is a sum of two phasors."""
        return 0.5 * self.k_max * ((self.w1 + self.w2) ** n + abs(self.w1 - self.w2) ** n)

    def bound(self) -> float:
        return departure_bound(self.k_max, self.w1, self.w2, self.dt, self.T)


# ---------------------------------------------------------------------- knots and integrals
class Knots:
    """One in-plane gradient as knots: absolute times, and x + i y."""

    def __init__(self, times: np.ndarray, values: np.ndarray) -> None:
        self.times, self.values = times, values

    @classmethod
    def of_events(cls, events: list, start: float = 0.0) -> Knots:
        """From extended-trapezoid events on x and y, placed at `start`."""
        axes = {}
        for event in events:
            assert event.type == 'grad'
            axes[event.channel] = (start + float(event.delay) + np.asarray(event.tt, float),
                                   np.asarray(event.waveform, float))
        return cls.combine(axes)

    @classmethod
    def of_sequence(cls, seq: pp.Sequence) -> Knots:
        """From every block of a sequence, junction knots merged; a jump fails the merge."""
        pieces: dict[str, tuple[list, list]] = {'x': ([], []), 'y': ([], [])}
        t = 0.0
        for n in sorted(seq.block_events):
            block = seq.get_block(n)
            for axis in 'xy':
                grad = getattr(block, f'g{axis}', None)
                if grad is not None:
                    assert grad.type == 'grad'
                    pieces[axis][0].append(t + float(grad.delay) + np.asarray(grad.tt, float))
                    pieces[axis][1].append(np.asarray(grad.waveform, float))
            t += float(seq.block_durations[n])
        axes = {}
        for axis, (ts, vs) in pieces.items():
            if ts:
                tt, vv = np.concatenate(ts), np.concatenate(vs)
                same = np.flatnonzero(np.abs(np.diff(tt)) < 1e-13)
                assert np.all(vv[same] == vv[same + 1]), f'gradient jump on {axis}'
                keep = np.ones(tt.size, bool)
                keep[same + 1] = False
                axes[axis] = (tt[keep], vv[keep])
        return cls.combine(axes)

    @classmethod
    def combine(cls, axes: dict) -> Knots:
        times = np.unique(np.concatenate([t for t, _ in axes.values()]))
        values = np.zeros(times.size, complex)
        for axis, (t, v) in axes.items():
            values += (1.0 if axis == 'x' else 1j) * np.interp(times, t, v, left=0.0, right=0.0)
        return cls(times, values)

    def k(self, at: np.ndarray) -> np.ndarray:
        """The exact integral from the first knot, held beyond the last."""
        t, g = self.times, self.values
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(t))])
        i = np.clip(np.searchsorted(t, at, side='right') - 1, 0, t.size - 2)
        span = t[i + 1] - t[i]
        s = np.clip(at - t[i], 0.0, span)
        return cum[i] + g[i] * s + 0.5 * (g[i + 1] - g[i]) / span * s * s

    def k_peak(self, part) -> float:
        """max |k| of one axis, exactly: at a knot or where that axis's gradient is zero."""
        g = part(self.values)
        t = self.times
        cross = np.flatnonzero(g[:-1] * g[1:] < 0)
        zeros = t[cross] - g[cross] * (t[cross + 1] - t[cross]) / (g[cross + 1] - g[cross])
        k = part(self.k(np.concatenate([t, zeros])))
        return float(np.max(np.abs(k)))

    def allowances(self) -> tuple[float, float]:
        """
        (registration, serialisation), as vectors in 1/m.

        Registration: a compressed shape's samples move by at most half a 1e-7 quantum of the
        peak.  Serialisation: the amplitude is written to six significant digits, a scale error
        of 5e-6 on k, and each normalised sample to nine.  Both carry a float floor.
        """
        duration = float(self.times[-1] - self.times[0])
        reg, ser = [], []
        for part in (np.real, np.imag):
            peak = float(np.max(np.abs(part(self.values))))
            reg.append(5e-8 * peak * duration)
            ser.append(5e-6 * self.k_peak(part) + 5e-9 * peak * duration)
        floor = 1e-12 * float(np.max(np.abs(self.values))) * duration
        return math.hypot(*reg) + floor, math.hypot(*ser) + floor


def read_back(seq: pp.Sequence) -> pp.Sequence:
    """Write `seq` to a file and read it again."""
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'rosette.seq'
        seq.write(str(path))
        back = pp.Sequence()
        back.read(str(path))
    return back


def compiled(ro: sc.modules.RosetteReadout, opts: pp.Opts, angle: float):
    """The designed, registered and read-back knots of one readout compiled alone."""
    block = ro(angle_rad=angle)
    seq = sc.compile(block, opts, name='rosette')
    designed = Knots.of_events([n.item for n in block.nodes if n.item.type == 'grad'])
    return block, seq, read_back(seq), designed


def core_start(knots: Knots, oracle: Analytic, angle: float) -> float:
    """
    Where the knots begin to play the analytic core: the first knot from which the next
    ``rasters`` knots are the analytic gradient on raster edges.
    """
    expected = oracle.g(np.arange(oracle.rasters + 1) * oracle.dt) * np.exp(1j * angle)
    scale = float(np.max(np.abs(expected)))
    for i in range(knots.times.size - oracle.rasters):
        t = knots.times[i:i + oracle.rasters + 1]
        if (np.allclose(np.diff(t), oracle.dt, rtol=0, atol=1e-12)
                and np.max(np.abs(knots.values[i:i + oracle.rasters + 1] - expected)) < 1e-5 * scale):
            return float(t[0])
    pytest.fail('the analytic core was not found in the knots')


def departure(knots: Knots, oracle: Analytic, t0: float, angle: float) -> float:
    """
    An upper bound on max |k - k_analytic| over the core, with its error controlled.

    The largest departure on 32 points per raster, plus half their spacing times the largest
    slope the departure can have: the knots' distance from the analytic gradient at the raster
    edges, plus how far a linear interpolant of it can stray, dt^2 M3 / 8.
    """
    points = 32
    tau = np.linspace(0.0, oracle.T, points * oracle.rasters + 1)
    analytic = oracle.k(tau) * np.exp(1j * angle)
    grid = float(np.max(np.abs(knots.k(t0 + tau) - analytic)))
    edges = np.arange(oracle.rasters + 1) * oracle.dt
    on_edges = np.interp(t0 + edges, knots.times, knots.values.real) \
        + 1j * np.interp(t0 + edges, knots.times, knots.values.imag)
    slope = float(np.max(np.abs(on_edges - oracle.g(edges) * np.exp(1j * angle)))) \
        + oracle.dt**2 / 8 * oracle.m(3)
    return grid + 0.5 * oracle.dt / points * slope


def crossings(knots: Knots, oracle: Analytic, t0: float) -> list[tuple[float, float]]:
    """
    Every local minimum of |k| below dk/2 near the core, refined: ``(time, |k|)``.

    Scanned at 64 points per raster from one raster before the core to one after it, never at
    the module's stride, then refined by golden section within a scan step either side.
    """
    step = oracle.dt / 64
    t = np.arange(t0 - oracle.dt, t0 + oracle.T + oracle.dt + step / 2, step)
    r = np.abs(knots.k(t))
    found = []
    for i in np.flatnonzero((r[1:-1] <= r[:-2]) & (r[1:-1] < r[2:]) & (r[1:-1] < oracle.dk / 2)) + 1:
        a, b = t[i - 1], t[i + 1]
        phi = (math.sqrt(5.0) - 1.0) / 2.0
        for _ in range(80):
            c, d = b - phi * (b - a), a + phi * (b - a)
            fc, fd = np.abs(knots.k(np.array([c, d])))
            a, b = (a, d) if fc <= fd else (c, b)
        tm = 0.5 * (a + b)
        found.append((float(tm), float(abs(knots.k(np.array([tm]))[0]))))
    return found


def adc_of(seq: pp.Sequence):
    """The one ADC of a sequence, with its block's start time."""
    t = 0.0
    for n in sorted(seq.block_events):
        block = seq.get_block(n)
        if block.adc is not None:
            return t, block.adc, float(seq.block_durations[n])
        t += float(seq.block_durations[n])
    pytest.fail('no ADC')


def sample_ticks(start: float, adc) -> np.ndarray:
    """Sample times of an ADC in integer picoseconds: delay + (j + 1/2) dwell."""
    dwell = round(float(adc.dwell) / TICK)
    first = round((start + float(adc.delay)) / TICK)
    return first + (2 * np.arange(int(adc.num_samples), dtype=np.int64) + 1) * dwell // 2


# --------------------------------------------------------------------------- the trajectory
@pytest.mark.parametrize('name', CASES)
def test_the_core_plays_the_analytic_gradient_on_the_raster(name: str) -> None:
    """
    Trajectory agreement with the analytic rosette (C1): the time law, measured.

    The core's knots are the analytic gradient on raster edges, and the compiled and read-back
    trajectories are its exact integral, within the allowance of each stage.
    """
    opts, ro, angle = case(name)
    oracle = Analytic(ro, opts)
    _, seq, back, designed = compiled(ro, opts, angle)
    t0 = core_start(designed, oracle, angle)
    edges = np.arange(oracle.rasters + 1) * oracle.dt
    raster_curve = Knots(t0 + edges, oracle.g(edges) * np.exp(1j * angle))
    tau = np.linspace(0.0, oracle.T, 16 * oracle.rasters + 1)
    expected = raster_curve.k(t0 + tau)

    reg, ser = designed.allowances()
    for knots, allowance in ((designed, 1e-12 * oracle.k_max), (Knots.of_sequence(seq), reg),
                             (Knots.of_sequence(back), reg + ser)):
        assert core_start(knots, oracle, angle) == pytest.approx(t0, abs=1e-12)
        error = float(np.max(np.abs(knots.k(t0 + tau) - expected)))
        assert error <= allowance, f'{error:.3g} 1/m against {allowance:.3g}'


@pytest.mark.parametrize('name', CASES)
def test_the_raster_departure_stays_within_its_allowance_at_every_stage(name: str) -> None:
    """
    The budget applies to the vector departure over the whole core, stage by stage (D7).

    ``departure`` is an error-controlled upper bound, not a grid maximum, so the first assertion
    is a full-interval claim.  The registered and read-back stages are allowed their own
    rounding on top of ``B``, never a wider ``B``.
    """
    opts, ro, angle = case(name)
    oracle = Analytic(ro, opts)
    _, seq, back, designed = compiled(ro, opts, angle)
    t0 = core_start(designed, oracle, angle)
    bound = oracle.bound()
    reg, ser = designed.allowances()
    assert bound <= 0.05 * oracle.dk

    stages = {'designed': (designed, bound), 'registered': (Knots.of_sequence(seq), bound + reg),
              'read back': (Knots.of_sequence(back), bound + reg + ser)}
    for stage, (knots, allowed) in stages.items():
        measured = departure(knots, oracle, t0, angle)
        assert measured <= allowed, f'{stage}: {measured / oracle.dk:.5g} dk > {allowed / oracle.dk:.5g}'

    # The reported departure is the module's own upper bound on the same quantity: never below
    # the largest value sampled here, and above B by no more than its own correction.
    tau = np.linspace(0.0, oracle.T, 32 * oracle.rasters + 1)
    sampled = float(np.max(np.abs(designed.k(t0 + tau) - oracle.k(tau) * np.exp(1j * angle))))
    assert sampled <= ro.departure_per_m <= bound + oracle.dt**3 * oracle.m(3) / 256


# ----------------------------------------------------------------------- the crossing report
@pytest.mark.parametrize('name', CASES)
def test_every_crossing_is_reported_as_ideal_closest_and_nearest_sample(name: str) -> None:
    """
    The crossing schedule (C2): P + 1 crossings, three facts each, never merged.

    Found on the read-back file by a search that does not use the module's instants, then
    compared: ideal instants exactly, closest approaches within what the file's rounding can
    move them, nearest samples by integer time with ties going to the earlier sample.
    """
    opts, ro, angle = case(name)
    oracle = Analytic(ro, opts)
    _, _, back, designed = compiled(ro, opts, angle)
    t0 = core_start(designed, oracle, angle)
    reg, ser = designed.allowances()
    rounding = reg + ser
    found = crossings(Knots.of_sequence(back), oracle, t0)

    assert len(found) == ro.petals + 1
    for field in ('origin_crossing_ideal_s', 'origin_crossing_closest_s',
                  'origin_crossing_residual_per_m', 'origin_crossing_samples'):
        assert isinstance(getattr(ro, field), tuple)
        assert len(getattr(ro, field)) == ro.petals + 1

    ideal_ticks = [round((t0 + p * oracle.petal) / TICK) for p in range(ro.petals + 1)]
    assert [round(t / TICK) for t in ro.origin_crossing_ideal_s] == ideal_ticks

    for p, (when, residual) in enumerate(found):
        assert abs(residual - ro.origin_crossing_residual_per_m[p]) <= rounding + 1e-9 * oracle.dk
        assert residual <= oracle.bound() + rounding
        reach = 2.0 * math.sqrt(rounding * (residual + rounding)) / (0.9 * oracle.k_max * oracle.w1)
        assert abs(when - ro.origin_crossing_closest_s[p]) <= reach + 1e-9

    start, adc, _ = adc_of(back)
    samples = sample_ticks(start, adc)
    ties = 0
    for p, at in enumerate(ideal_ticks):
        distance = np.abs(samples - at)
        nearest = np.flatnonzero(distance == distance.min())
        ties += nearest.size > 1
        assert ro.origin_crossing_samples[p] == int(nearest[0])
    if name == 'record-13':
        # The dwell divides the petal, so every crossing sits half a dwell from two samples.
        assert ties == ro.petals + 1


@pytest.mark.parametrize('name', CASES)
def test_the_samples_cover_every_crossing_and_its_closest_approach(name: str) -> None:
    """
    Coverage, on the read-back file: the first sample is at or before the first crossing and
    its closest approach, the last at or after the last crossing and its closest approach.
    """
    opts, ro, angle = case(name)
    oracle = Analytic(ro, opts)
    _, _, back, designed = compiled(ro, opts, angle)
    t0 = core_start(designed, oracle, angle)
    found = crossings(Knots.of_sequence(back), oracle, t0)
    start, adc, _ = adc_of(back)
    times = sample_ticks(start, adc) * TICK
    assert times[0] <= min(t0, found[0][0]) + 1e-12
    assert times[-1] >= max(t0 + oracle.T, found[-1][0]) - 1e-12


def test_coverage_is_checked_not_assumed() -> None:
    """
    At a 200 ns dwell, one dwell of padding put on its lattice misses the last closest
    approach, and the window grows until it covers it.
    """
    opts, ro, _ = case('wide-1-200ns')
    oracle = Analytic(ro, opts)
    _, _, back, designed = compiled(ro, opts, 0.0)
    t0 = core_start(designed, oracle, 0.0)
    last_closest = crossings(Knots.of_sequence(back), oracle, t0)[-1][0]

    dwell, rf = round(ro.dwell_s / TICK), round(float(opts.rf_raster_time) / TICK)
    t0_ticks, end_ticks = round(t0 / TICK), round((t0 + oracle.T) / TICK)
    delay = (t0_ticks - dwell) // rf * rf
    count = -(-(end_ticks - delay) // dwell) + 1
    count += (-count) % 4
    padded_last = (delay + (2 * count - 1) * dwell // 2) * TICK
    assert padded_last < last_closest

    assert ro.num_samples > count
    assert ro.sample_times_s()[-1] >= last_closest


def test_a_window_that_cannot_open_early_enough_is_refused() -> None:
    """
    ``pp.make_adc`` holds a delay to at least ``adc_dead_time``.  At 40 us the first sample
    lands exactly on the first crossing, and at 41 us after it.
    """
    opts = scanner()
    common = {'fov_mm': 220.0, 'matrix': 64, 'petals': 1, 'petal_duration_s': 10e-3}
    ro = rosette(opts, dwell_s=40e-6, **common)
    assert ro.sample_times_s()[0] == pytest.approx(ro.origin_crossing_ideal_s[0], abs=1e-12)
    with pytest.raises(sc.ConfigurationError, match=r'first sample falls 0\.5 us after'):
        rosette(opts, dwell_s=41e-6, **common)


# --------------------------------------------------------------------------------- the ADC
@pytest.mark.parametrize('name', CASES)
def test_the_emitted_adc_is_legal_and_its_times_are_pulseqs(name: str) -> None:
    """
    The emitted ADC, read back: on its lattice, within its limit, with its dead times inside
    the block, its sample times PyPulseq's own, and padding samples exactly those outside the
    core.
    """
    opts, ro, angle = case(name)
    oracle = Analytic(ro, opts)
    block, _, back, designed = compiled(ro, opts, angle)
    start, adc, span = adc_of(back)
    delay = float(adc.delay)
    assert start == 0.0
    assert abs(delay / float(opts.rf_raster_time) - round(delay / float(opts.rf_raster_time))) < 1e-9
    assert abs(float(adc.dwell) / float(opts.adc_raster_time)
               - round(float(adc.dwell) / float(opts.adc_raster_time))) < 1e-9
    assert int(adc.num_samples) == ro.num_samples
    assert ro.num_samples % 4 == 0 and ro.num_samples <= opts.adc_samples_limit
    assert delay >= float(opts.adc_dead_time) - 1e-12
    assert span >= delay + ro.num_samples * float(adc.dwell) + float(opts.adc_dead_time) - 1e-12

    k_pp, t_pp = back.calculate_kspacePP()[:2]
    assert np.max(np.abs(ro.sample_times_s() - np.asarray(t_pp))) < 1e-12
    reg, ser = designed.allowances()
    reported = ro.k_per_m(angle_rad=angle)
    assert np.max(np.hypot(*(np.asarray(k_pp)[:2] - reported))) <= reg + ser

    t0 = core_start(designed, oracle, angle)
    times = sample_ticks(start, adc)
    inside = (times >= round(t0 / TICK)) & (times <= round((t0 + oracle.T) / TICK))
    assert np.array_equal(np.flatnonzero(inside), np.arange(ro.num_samples)[ro.core_samples])
    assert block.duration == pytest.approx(span, abs=1e-12)


# ------------------------------------------------------------------------------ boundaries
@pytest.mark.parametrize('name', CASES)
def test_the_readout_starts_and_ends_at_rest_with_no_net_moment(name: str) -> None:
    """
    At rest at both ends, continuous at every junction, and the end point of the core cancelled,
    in memory and in the file (C4).
    """
    opts, ro, angle = case(name)
    _, seq, back, designed = compiled(ro, opts, angle)
    reg, ser = designed.allowances()
    duration = float(designed.times[-1] - designed.times[0])
    net_floor = math.hypot(*(5e-9 * float(np.max(np.abs(part(designed.values)))) * duration
                             for part in (np.real, np.imag)))
    for knots, allowed in ((Knots.of_sequence(seq), reg), (Knots.of_sequence(back), net_floor)):
        assert knots.values[0] == 0.0 and knots.values[-1] == 0.0
        assert knots.times[0] == 0.0
        assert np.all(np.diff(knots.times) > 0.0)
        net = abs(knots.k(np.array([knots.times[-1]]))[0])
        assert net <= allowed + 1e-12 * float(np.max(np.abs(knots.values))) * duration


# -------------------------------------------------------------------------------- rotation
def test_every_in_plane_rotation_is_legal() -> None:
    """
    The exact bound at every rotation is read off the angle-0 vectors (C6): |g| at the knots,
    where a piecewise-linear |g| peaks, and |dg/dt| on each segment.  ``limits()`` is that
    measurement, and the worst angles compile.
    """
    opts = scanner()
    ro = rosette(opts, petals=5)
    knots = Knots.of_events([n.item for n in ro().nodes if n.item.type == 'grad'])
    slew = np.diff(knots.values) / np.diff(knots.times)
    peak, steepest = float(np.max(np.abs(knots.values))), float(np.max(np.abs(slew)))
    assert peak <= float(opts.max_grad)
    assert steepest <= float(opts.max_slew)
    assert ro.limits() == pytest.approx((peak, steepest), rel=1e-12)
    assert steepest > 0.99 * float(opts.max_slew)

    worst = [-float(np.angle(knots.values[np.argmax(np.abs(knots.values))])),
             -float(np.angle(slew[np.argmax(np.abs(slew))]))]
    for angle in (*worst, 0.0, 0.7, *(2 * math.pi * i / 5 + 0.1 for i in range(5))):
        sc.compile(ro(angle_rad=angle), opts, name='rotated')


def test_rotation_is_a_rotation() -> None:
    """Rotation equivariance (C7), in the report and in the file."""
    opts = scanner()
    ro = rosette(opts, petals=5)
    turn = np.exp(1j * 0.7)
    base, turned = ro.k_per_m(), ro.k_per_m(angle_rad=0.7)
    assert np.max(np.abs((base[0] + 1j * base[1]) * turn - (turned[0] + 1j * turned[1]))) \
        < 1e-12 * ro.k_max_per_m

    files = {}
    for angle in (0.0, 0.7):
        _, _, back, designed = compiled(ro, opts, angle)
        k = np.asarray(back.calculate_kspacePP()[0])
        files[angle] = (k[0] + 1j * k[1], sum(designed.allowances()))
    rounding = files[0.0][1] + files[0.7][1]
    assert np.max(np.abs(files[0.0][0] * turn - files[0.7][0])) <= rounding


# --------------------------------------------------------------------------- time reversal
@pytest.mark.parametrize('petals', [1, 2, 7])
def test_time_reversal_holds_up_to_the_closure_residual(petals: int) -> None:
    """
    On the raster the core reverses as ``k(T - t) - k(T) = c conj(k(t))``, with
    ``c = (-1)^(P+1) exp(i w2 T)`` (C8).  It checks the crossing schedule's symmetry, not the
    time law: ``test_a_retimed_traversal_passes_reversal_and_fails_the_time_law``.
    """
    opts = scanner()
    ro = rosette(opts, petals=petals, petal_duration_s=810e-6)
    oracle = Analytic(ro, opts)
    seq = sc.compile(ro(), opts, name='rosette')
    knots = Knots.of_sequence(seq)
    t0 = core_start(knots, oracle, 0.0)
    tau = np.linspace(0.0, oracle.T, 8 * oracle.rasters + 1)
    k = knots.k(t0 + tau) - knots.k(np.array([t0]))[0]
    c = (-1) ** (petals + 1) * np.exp(1j * oracle.w2 * oracle.T)
    assert np.max(np.abs(k[::-1] - k[-1] - c * np.conj(k))) < 1e-9 * oracle.dk


@pytest.mark.parametrize('epsilon', [0.05, 0.2])
def test_a_retimed_traversal_passes_reversal_and_fails_the_time_law(epsilon: float) -> None:
    """
    The negative control.  ``u(t) = t + epsilon / (2 w1) sin(2 w1 t)`` is monotone and fixes
    every crossing, so the re-timed curve passes the reversal check -- and misses the analytic
    curve by more than a dk, which only the trajectory comparison sees.
    """
    opts = scanner()
    oracle = Analytic(rosette(opts, petals=7, petal_duration_s=810e-6), opts)
    w1 = oracle.w1
    edges = np.arange(oracle.rasters + 1) * oracle.dt
    u = edges + epsilon / (2 * w1) * np.sin(2 * w1 * edges)
    retimed = Knots(edges, oracle.g(u) * (1 + epsilon * np.cos(2 * w1 * edges)))
    tau = np.linspace(0.0, oracle.T, 8 * oracle.rasters + 1)
    k = retimed.k(tau)
    c = (-1) ** (oracle.petals + 1) * np.exp(1j * oracle.w2 * oracle.T)
    assert np.max(np.abs(k[::-1] - k[-1] - c * np.conj(k))) < 1e-9 * oracle.dk
    assert np.max(np.abs(k - oracle.k(tau))) > oracle.dk


# ----------------------------------------------------------------------------- one reading
@pytest.mark.parametrize('name', ['record-13', 'golden-5', 'wide-13-20us'])
def test_the_module_the_compiler_and_the_file_read_one_trajectory(name: str) -> None:
    """
    The reported k is what the compiler holds and what the file plays (C9): extended
    trapezoids with knots on the raster, one per axis, read alike by every path.
    """
    opts, ro, angle = case(name)
    block, seq, back, designed = compiled(ro, opts, angle)
    reg, ser = designed.allowances()
    reported = ro.k_per_m(angle_rad=angle)
    in_memory = sc.kspace(block, opts)['k_adc'][:2]
    assert np.max(np.hypot(*(in_memory - reported))) <= reg
    k_pp = np.asarray(back.calculate_kspacePP()[0])[:2]
    assert np.max(np.hypot(*(k_pp - reported))) <= reg + ser

    raster = float(opts.grad_raster_time)
    for source in (seq, back):
        assert len(source.block_events) == 1
        events = [getattr(source.get_block(1), f'g{axis}') for axis in 'xy']
        for grad in events:
            assert grad.type == 'grad'
            tt = np.asarray(grad.tt, float)
            assert np.allclose(tt / raster, np.round(tt / raster), atol=1e-6)
            assert float(grad.first) == 0.0 and float(grad.last) == 0.0


# ---------------------------------------------------------------- sampling and demands
@pytest.mark.parametrize('name', CASES)
def test_no_two_samples_are_more_than_dk_apart(name: str) -> None:
    """Along-path sampling (C10), on the read-back file."""
    opts, ro, angle = case(name)
    _, _, back, _ = compiled(ro, opts, angle)
    k = np.asarray(back.calculate_kspacePP()[0])
    assert np.max(np.hypot(*np.diff(k[:2], axis=1))) <= ro.dk_per_m


def test_the_along_path_limit_is_refused_one_raster_above_it() -> None:
    """At a 790 us petal the fastest point of the curve allows 3.9 us, and 4.0 us is refused."""
    opts = scanner()
    rosette(opts, petals=1, dwell_s=3.9e-6)
    with pytest.raises(sc.ConfigurationError, match=r'moves k by 1\.02 dk.*\n(.*\n)*.*3\.9 us'):
        rosette(opts, petals=1, dwell_s=4e-6)


@pytest.mark.parametrize('ratio', [0.5385, 1.5, 3.0])
@pytest.mark.parametrize('shift', [0, 1])
def test_the_peak_demands_are_the_closed_forms(ratio: float, shift: int) -> None:
    """
    The core's peak |g| and slew against ``k_max max(w1, w2)`` and ``k_max (w1^2 + w2^2)`` (C5).

    Sampled on raster edges, the |g| peak is exact unless ``w2 > w1`` puts it mid-raster, which
    an odd raster count per petal does; then it is exactly the value half a raster away.  The
    sampled slew is a secant, at most the peak and at least ``cos^2(w1 dt / 2)`` of it, less the
    midpoint error ``dt^2 M4 / 24``.
    """
    opts = scanner(80, 200)
    default = rosette(opts, fov_mm=256.0, petals=2, omega_ratio=ratio)
    rasters = round(default.petal_duration_s / 10e-6) + shift
    ro = rosette(opts, fov_mm=256.0, petals=2, omega_ratio=ratio,
                 petal_duration_s=rasters * 10e-6)
    oracle = Analytic(ro, opts)
    knots = Knots.of_events([n.item for n in ro().nodes if n.item.type == 'grad'])
    t0 = core_start(knots, oracle, 0.0)
    core = (knots.times >= t0 - 1e-12) & (knots.times <= t0 + oracle.T + 1e-12)
    g = knots.values[core]
    peak = float(np.max(np.abs(g)))
    w1, w2, dt = oracle.w1, oracle.w2, oracle.dt
    if w2 > w1 and rasters % 2:
        a, b = w1 + w2, w1 - w2
        expected = 0.5 * oracle.k_max * math.sqrt(a * a + b * b - 2 * a * b * math.cos(w1 * dt))
    else:
        expected = oracle.k_max * max(w1, w2)
    assert peak == pytest.approx(expected, rel=1e-9)

    slew = float(np.max(np.abs(np.diff(g)))) / dt
    top = oracle.k_max * (w1**2 + w2**2)
    assert top * math.cos(w1 * dt / 2) ** 2 - dt**2 / 24 * oracle.m(4) <= slew <= top


# -------------------------------------------------------------------------- the block clock
def test_reported_times_are_on_the_block_clock() -> None:
    """
    The block clock (C3): sample times are what the compiler reads, and placed after an
    excitation the crossings fall where the report and the placement put them.
    """
    opts = scanner()
    ro = rosette(opts, petals=5)
    assert np.max(np.abs(sc.kspace(ro(), opts)['t_adc'] - ro.sample_times_s())) < 1e-12

    exc = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0, duration_s=1e-3)
    first = exc()
    tree = sc.LogicBlock('gre').add(0.0, first).add(first.duration, ro())
    back = read_back(sc.compile(tree, opts, name='gre'))
    oracle = Analytic(ro, opts)
    knots = Knots.of_sequence(back)
    t0 = core_start(knots, oracle, 0.0)
    assert round(t0 / TICK) == round((first.duration + ro.origin_crossing_ideal_s[0]) / TICK)

    block = back.get_block(sorted(back.block_events)[0])
    rf_centre = float(block.rf.delay) + float(pp.calc_rf_center(block.rf)[0])
    assert rf_centre == pytest.approx(exc.time_to_center(), abs=1e-9)
    found = crossings(knots, oracle, t0)
    for p, (when, _) in enumerate(found):
        assert abs(when - first.duration - ro.origin_crossing_closest_s[p]) < 0.5 * 100e-9


# -------------------------------------------------------------------------------- the domain
REFUSALS = [
    ({'omega_ratio': -0.5}, r'omega_ratio must be a finite number >= 0, got -0\.5'),
    ({'omega_ratio': math.inf}, r'omega_ratio must be a finite number >= 0'),
    ({'petal_duration_s': 805e-6}, r'not a whole grad_raster_time(.*\n)*.*800 us or 810 us'),
    ({'dwell_s': 2.05e-6}, r'not a whole adc_raster_time(.*\n)*.*2 us or 2\.1 us'),
    ({'petal_duration_s': 780e-6}, r'780 us petal needs 122\.9 T/m/s(.*\n)*.*>= 790 us'),
    ({'petals': 21, 'petal_duration_s': 810e-6},
     r'21 petals need 8508 samples(.*\n)*.*at most 20(.*\n)*.*does not split'),
]


@pytest.mark.parametrize(('overrides', 'match'), REFUSALS)
def test_each_refusal_names_its_way_out(overrides: dict, match: str) -> None:
    """Each edge of the domain at its measured boundary, refused with the value that fits."""
    with pytest.raises(sc.ConfigurationError, match=match):
        rosette(scanner(), **{'petals': 1, **overrides})


@pytest.mark.parametrize(('hardware', 'refused', 'accepted', 'match'), [
    ((20, 80), {'omega_ratio': 0.0, 'petal_duration_s': 920e-6},
     {'omega_ratio': 0.0, 'petal_duration_s': 930e-6}, r'needs 20\.05 mT/m(.*\n)*.*>= 930 us'),
    ((80, 200), {'fov_mm': 400.0, 'petal_duration_s': 490e-6},
     {'fov_mm': 400.0, 'petal_duration_s': 540e-6},
     r'departure bound for a 490 us petal is 0\.0592\d dk(.*\n)*.*not a physical limit(.*\n)*.*>= 540 us'),
    ((80, 200), {'fov_mm': 400.0, 'omega_ratio': 0.0, 'petal_duration_s': 440e-6},
     {'fov_mm': 400.0, 'omega_ratio': 0.0, 'petal_duration_s': 470e-6},
     r'0\.0563\d dk(.*\n)*.*>= 470 us'),
    ((32, 120), {'petal_duration_s': 780e-6}, {'petal_duration_s': 790e-6}, r'122\.9 T/m/s'),
    ((32, 120), {'petals': 21, 'petal_duration_s': 810e-6},
     {'petals': 20, 'petal_duration_s': 810e-6}, r'8508 samples'),
])
def test_the_accepted_side_of_each_boundary_compiles(hardware, refused, accepted, match) -> None:
    """The gradient, budget, slew and sample-limit boundaries, one raster or one petal apart."""
    opts = scanner(*hardware)
    with pytest.raises(sc.ConfigurationError, match=match):
        rosette(opts, **{'petals': 1, **refused})
    sc.compile(rosette(opts, **{'petals': 1, **accepted})(), opts, name='accepted')


def test_a_readout_needs_a_configured_event_limit() -> None:
    """``adc_samples_limit = 0`` sets no limit, and this readout refuses to have none."""
    opts = scanner(limit=0)
    with pytest.raises(sc.ConfigurationError, match='adc_samples_limit is 0'):
        rosette(opts, petals=1)


@pytest.mark.parametrize(('hardware', 'raster', 'overrides'), [
    ((32, 120), 10e-6, {'omega_ratio': 7 / 13}),
    ((20, 80), 10e-6, {'omega_ratio': 0.0}),
    ((80, 200), 10e-6, {'fov_mm': 400.0, 'omega_ratio': 7 / 13}),
    ((32, 120), 20e-6, {'omega_ratio': 1.5}),
])
def test_the_default_petal_is_the_shortest_that_fits(hardware, raster, overrides) -> None:
    """
    The default meets the vector gradient and slew limits and the budget, one raster less fails
    one of them, and the dwell plays no part.
    """
    opts = scanner(*hardware, raster=raster)
    ro = rosette(opts, petals=2, **overrides)
    assert rosette(opts, petals=2, dwell_s=1e-6, **overrides).petal_duration_s \
        == ro.petal_duration_s

    k_max, dk = ro.k_max_per_m, ro.dk_per_m

    def fits(petal: float) -> bool:
        w1 = math.pi / petal
        w2 = ro.omega_ratio * w1
        return (k_max * max(w1, w2) <= float(opts.max_grad)
                and k_max * (w1**2 + w2**2) <= float(opts.max_slew)
                and departure_bound(k_max, w1, w2, raster, 2 * petal) <= 0.05 * dk)

    assert fits(ro.petal_duration_s)
    assert not fits(ro.petal_duration_s - raster)


# ----------------------------------------------------------------------- the module conventions
def test_it_compiles_alone_and_a_call_is_pure(component_checks) -> None:
    """Twice-called, unchanged, deterministic, and legal on its own."""
    ro = rosette(scanner(), petals=2)
    component_checks.all(ro, angle_rad=0.7, phase_deg=30.0)


def test_a_dummy_shot_plays_the_same_gradients_over_the_same_span() -> None:
    """``acquire=False`` drops the ADC and nothing else: same gradients, same block length."""
    opts = scanner()
    ro = rosette(opts, petals=2)
    loud, quiet = ro(angle_rad=0.7), ro(angle_rad=0.7, acquire=False)
    assert quiet.duration == loud.duration
    assert not any(n.item.type == 'adc' for n in quiet.nodes)
    grads = [[sc.events.content_hash(n.item) for n in b.nodes if n.item.type == 'grad']
             for b in (loud, quiet)]
    assert grads[0] == grads[1]
    seq = sc.compile(quiet, opts, name='dummy')
    assert all(seq.get_block(n).adc is None for n in seq.block_events)


def test_a_window_that_outlasts_the_gradient_holds_the_block_open_at_rest() -> None:
    """
    A slow petal with a long dwell: the ADC's trailing dead time ends after the gradient.

    The block lasts until the later of the two, with or without the ADC; the gradient is zero
    after its last knot; and the samples taken there sit at the end point, which is the origin.
    """
    opts, ro, angle = case('resting-34us')
    block, seq, back, designed = compiled(ro, opts, angle)
    gradient_end = float(designed.times[-1])
    start, adc, span = adc_of(back)
    adc_end = float(adc.delay) + ro.num_samples * float(adc.dwell) + float(opts.adc_dead_time)
    assert adc_end > gradient_end
    assert block.duration == pytest.approx(math.ceil(adc_end / 10e-6 - 1e-9) * 10e-6, abs=1e-12)
    assert ro(angle_rad=angle, acquire=False).duration == block.duration
    assert sc.compile(ro(angle_rad=angle, acquire=False), opts).duration()[0] == pytest.approx(span)

    knots = Knots.of_sequence(back)
    assert knots.times[-1] == pytest.approx(gradient_end, abs=1e-12)
    assert knots.values[-1] == 0.0
    k_pp, t_pp = back.calculate_kspacePP()[:2]
    resting = np.asarray(t_pp) > gradient_end
    assert resting.any()
    reg, ser = designed.allowances()
    assert np.max(np.hypot(*np.asarray(k_pp)[:2, resting])) <= reg + ser


# --------------------------------------------------------------------------- composition
def test_a_small_acquisition_reads_back_to_the_reported_trajectory() -> None:
    """
    A dummy shot and four shots: excitation, rosette, spoiler, fixed TR, RF spoiling.

    Read back, each shot's samples are where the readout says at that angle, the slice axis is
    rephased at every crossing, and the dummy shot plays the same gradients without sampling.
    """
    opts = scanner()
    exc = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=5.0, duration_s=1e-3)
    ro = rosette(opts, petals=5)
    spoil = pp.make_trapezoid(channel='z', area=4.0 / 5e-3, system=opts)
    first = exc()
    readout_at = first.duration
    z_rounding = 5e-6 * sum(abs(float(event.area)) for _, event, _ in sc.flatten(first)
                            if getattr(event, 'channel', None) == 'z')
    assert z_rounding > 0.0
    tr = math.ceil((readout_at + ro().duration + pp.calc_duration(spoil)) / 1e-4) * 1e-4 + 1e-3
    angles = [2 * math.pi * s / 4 for s in range(4)]

    scan = sc.LogicBlock('rosette_gre')
    for shot in range(-1, 4):
        phase = 117.0 * (shot + 1) * (shot + 2) / 2 % 360.0
        angle = angles[shot] if shot >= 0 else 0.0
        at = (shot + 1) * tr
        readout = ro(angle_rad=angle, phase_deg=phase, acquire=shot >= 0)
        scan.add(at, exc(phase_deg=phase), pp.make_delay(tr))
        scan.add(at + readout_at, readout)
        scan.add(at + readout_at + readout.duration, spoil)
    back = read_back(sc.compile(scan, opts, name='rosette_gre'))

    k_pp, t_pp = back.calculate_kspacePP()[:2]
    k_pp, t_pp = np.asarray(k_pp), np.asarray(t_pp)
    assert k_pp.shape[1] == 4 * ro.num_samples
    for shot in range(4):
        designed = Knots.of_events([n.item for n in ro(angle_rad=angles[shot]).nodes
                                    if n.item.type == 'grad'])
        rounding = sum(designed.allowances())
        rows = slice(shot * ro.num_samples, (shot + 1) * ro.num_samples)
        expected_t = (shot + 1) * tr + readout_at + ro.sample_times_s()
        assert np.max(np.abs(t_pp[rows] - expected_t)) < 1e-9
        assert np.max(np.hypot(*(k_pp[:2, rows] - ro.k_per_m(angle_rad=angles[shot])))) <= rounding
        # The readout plays nothing on z, so kz is the excitation's: zero as designed, and in
        # the file within the six-digit rounding of its own lobes.
        crossing_kz = k_pp[2, rows][list(ro.origin_crossing_samples)]
        assert np.ptp(crossing_kz) < 1e-9
        assert np.max(np.abs(crossing_kz)) <= z_rounding

    knots = Knots.of_sequence(back)
    window = np.linspace(0.0, ro().duration, 4001)

    def played(at: float) -> np.ndarray:
        return np.interp(at + window, knots.times, knots.values.real) \
            + 1j * np.interp(at + window, knots.times, knots.values.imag)

    dummy, shot_0 = played(readout_at), played(tr + readout_at)
    assert np.max(np.abs(dummy - shot_0)) < 1e-12 * np.max(np.abs(shot_0))
    assert back.duration()[0] == pytest.approx(5 * tr)
