r"""
:class:`RosetteReadout` -- an in-plane rosette readout with the analytic time law.

``readout/`` because it contains an ADC and no RF.

The trajectory
--------------
The analytic core plays

.. math:: k(t) = k_{max} \sin(\omega_1 t)\, e^{i \omega_2 t}, \qquad 0 \le t \le T = P\pi/\omega_1

in the x-y plane, x the real part and y the imaginary part.  ``matrix`` and ``fov_mm`` set
:math:`k_{max} = \mathrm{matrix} / (2\,\mathrm{FOV})`, one petal lasts :math:`\pi/\omega_1`, and
``omega_ratio`` is :math:`\omega_2/\omega_1`: ``0`` is a spike through the centre, below 1 gives
narrow petals, ``1`` circles through the origin, and above 1 wide petals.  The core passes
through the origin at :math:`t_p = p\pi/\omega_1` for ``p = 0 ... P``, which is ``P + 1``
crossings counting both ends.  A shot need not close on itself.

At ``angle_rad = 0`` the gradient at the first crossing points along +x and the petals advance
counter-clockwise.

The gradient is the analytic :math:`dk/dt` sampled on the gradient raster's edges and played
piecewise linearly, one extended trapezoid per in-plane axis.  The block holds, in order::

    lead-in      along g(0): from rest down to -a, then up to |g(0)|, with zero area
    core         the analytic gradient on raster edges, t = 0 ... T
    lead-out     along the core's last gradient g(T): from |g(T)| down to -b, then to rest,
                 cancelling the core's end point along g(T)
    correction   a triangle perpendicular to g(T), cancelling the rest of the end point

so the readout starts and ends at rest, consecutive parts share a knot, and the in-plane moment
of the whole block is zero.  Each two-piece ramp takes the fewest gradient rasters the slew limit
allows, and the correction the fewest for its area.  There is no prephaser and no rewinder.

The departure budget
--------------------
On the raster the core departs from the analytic curve.  Its largest vector departure over the
core is at most

.. math::

    B = \frac{\Delta t^2}{6} M_2 + \frac{\Delta t^3}{12} M_3 + \frac{T \Delta t^4}{720} M_5,
    \qquad M_n = \frac{k_{max}}{2}\left((\omega_1 + \omega_2)^n + |\omega_1 - \omega_2|^n\right)

with :math:`\Delta t` the gradient raster: the Euler--Maclaurin error of the trapezoidal rule at
the raster points, the most the error can grow inside one raster, and the remainder.  This readout
supports protocols with :math:`B \le 0.05\,\Delta k`, :math:`\Delta k = 1/\mathrm{FOV}`.  That is
a fixed engineering budget, not a scanner or image-quality guarantee.  The bound is conservative,
so a petal it refuses may in fact depart by less: a refusal marks the supported domain, not a
physical limit.

The gradient and slew limits are applied to the in-plane **vector**, so every rotation is legal
on every axis.  The analytic peaks are :math:`k_{max}\max(\omega_1, \omega_2)` and
:math:`k_{max}(\omega_1^2 + \omega_2^2)`; the sampled core stays below both.

What it reports, and from what
------------------------------
Every time is on the **block clock**, from the start of the block :meth:`~RosetteReadout.build`
returns.  Every trajectory report is computed from the gradient knots this module designs,
before any compiler registers them, and is the same at every rotation except ``k_per_m``.  A
compiled sequence holds those knots as pypulseq registered them, which may round a compressible
shape's samples to 1e-7 of its peak, and a written file stores each amplitude to six significant
digits.  Neither is exactly what this module reports.

- **Ideal crossings** (``origin_crossing_ideal_s``) are the formula's instants: the core start plus
  ``p`` petals.  They are on the gradient raster.
- **Closest approaches** (``origin_crossing_closest_s``, ``origin_crossing_residual_per_m``) are
  where the designed trajectory comes nearest the origin around each ideal instant, and how near.
  They differ from the ideal instants because the raster departs from the curve, and the end
  crossing's may fall just inside the lead-out.
- **Nearest samples** (``origin_crossing_samples``) are the ADC samples whose times are closest to
  each ideal instant, compared in integer picoseconds.  When two are equally close, the earlier is
  reported.

There is no scalar echo time: which crossing a protocol calls TE is the composing kernel's choice.

The ADC window
--------------
One ADC event, opened one dwell before the core starts and sized so its last sample lies at least
half a dwell after the core ends.  The delay is floored onto ``rf_raster_time``, and the sample
count rounded up to ``adc_samples_divisor``, extra samples going after the core.  The delay is
read back from the built event, because ``pp.make_adc`` raises a delay below ``adc_dead_time`` to
it.  The window then has to contain every ideal crossing and every closest approach: when the last
one falls after the last sample, one divisor of samples is added, and when the window cannot open
early enough the readout is refused.  Samples outside the core are real data, taken on the lead-in
before it and on the lead-out, the correction or the gradient at rest after it; ``core_samples``
separates them, and ``k_per_m`` reports where every sample is.

``opts.adc_samples_limit`` must be positive.  A readout needing more samples than one event may
hold is refused: this readout does not split its window.

The block lasts until both the gradient and the ADC's trailing dead time have ended, the same
with or without the ADC, so a dummy shot occupies the same time.  When the ADC outlasts the
gradient the gradient is zero for the rest of the block and the trajectory rests at its end point.

Out of scope
------------
A minimum-time traversal, ADC segmentation, negative ratios and 3D.  The excitation, the
shot-angle schedule, RF spoiling, dummy shots, TR, flip angle and contrast belong to whatever
composes the readout.

Examples
--------
>>> import pypulseq as pp
>>> import seqcraft as sc
>>> opts = pp.Opts(max_grad=32, grad_unit='mT/m', max_slew=120, slew_unit='T/m/s',
...                adc_dead_time=10e-6, adc_samples_limit=8192)
>>> ro = sc.modules.RosetteReadout(opts=opts, fov_mm=256.0, matrix=128, petals=13,
...                                omega_ratio=7 / 13, dwell_s=2e-6)
>>> round(ro.petal_duration_s * 1e6)
790
>>> len(ro.origin_crossing_ideal_s), ro.num_samples
(14, 5140)
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ...design.events import derive
from ...design.logic import LogicBlock
from ...design.module import Module
from ...design.timing import Raster, from_ticks, to_ticks
from ...errors import ConfigurationError, format_error
from .._support import require_count, require_positive

if TYPE_CHECKING:
    from pypulseq.opts import Opts

__all__ = ['RosetteReadout']

#: The supported raster departure, in units of dk = 1/FOV.
_BUDGET_DK = 0.05

#: Departure samples per raster interval; the bound on what lies between them is added.
_DEPARTURE_POINTS = 16

#: Relative margin kept below the gradient and slew limits, so that float rounding in a rotation
#: cannot carry an emitted value over a limit pypulseq checks exactly.
_GUARD = 1e-9


# ------------------------------------------------------------------ the analytic rosette
def _k(t: np.ndarray, k_max: float, w1: float, w2: float) -> np.ndarray:
    """The analytic trajectory, complex, 1/m."""
    return k_max * np.sin(w1 * t) * np.exp(1j * w2 * t)


def _g(t: np.ndarray, k_max: float, w1: float, w2: float) -> np.ndarray:
    """Its time derivative, complex, Hz/m."""
    return k_max * (w1 * np.cos(w1 * t) + 1j * w2 * np.sin(w1 * t)) * np.exp(1j * w2 * t)


def _derivative_bound(k_max: float, w1: float, w2: float, order: int) -> float:
    """``M_n``: a bound on the n-th derivative of k over all time."""
    return 0.5 * k_max * ((w1 + w2) ** order + abs(w1 - w2) ** order)


def _departure_bound(k_max: float, w1: float, w2: float, raster: float, core: float) -> float:
    """``B``: the core's largest vector departure from the curve, at most."""
    def m(n: int) -> float:
        return _derivative_bound(k_max, w1, w2, n)
    return raster**2 / 6 * m(2) + raster**3 / 12 * m(3) + core * raster**4 / 720 * m(5)


# ------------------------------------------------------- the ramps around the core
def _ramp_to_rest(start: float, area: float, slew: float, raster: float) -> tuple[int, int, float]:
    """
    Return ``(m1, m2, b)`` for ``start -> -b -> 0`` over ``m1`` then ``m2`` rasters with `area`.

    The fewest rasters in total, then the fewest in the first piece.  For a fixed total, a longer
    first piece only eases its own slew and only tightens the second's, so the first piece that
    passes its own two conditions is the only candidate worth checking.
    """
    alpha = 2.0 * area / raster
    total = 2
    while True:
        excess = slew * raster * total - start
        if excess > 0.0:
            floor = max((start * total - alpha) / excess, alpha / start, 1.0)
            m1 = max(1, math.ceil(floor) - 1)
            while m1 < total:
                b = (m1 * start - alpha) / total
                if b >= 0.0 and (start + b) / (m1 * raster) <= slew:
                    break
                m1 += 1
            if m1 < total and (m1 * start - alpha) / total / ((total - m1) * raster) <= slew:
                return m1, total - m1, (m1 * start - alpha) / total
        total += 1


def _lead_in(start: float, slew: float, raster: float) -> tuple[int, int, float]:
    """
    Return ``(n1, n2, a)`` for ``0 -> -a -> start`` over ``n1`` then ``n2`` rasters, zero area.

    The fewest rasters in total, then the fewest in the first piece.  ``a = start n2 / (n1 + n2)``
    makes the area exactly zero; a longer first piece eases its own slew and tightens the second's.
    """
    total = 2
    while True:
        n1 = max(1, math.ceil(start * total / (start + slew * raster * total)) - 1)
        while n1 < total and start * (total - n1) / total / (n1 * raster) > slew:
            n1 += 1
        if n1 < total:
            a = start * (total - n1) / total
            if (start + a) / ((total - n1) * raster) <= slew:
                return n1, total - n1, a
        total += 1


def _triangle(area: float, slew: float, raster: float) -> tuple[int, float]:
    """Return ``(q, h)`` for ``0 -> h -> 0`` over ``2 q`` rasters with `area`, ``q`` smallest."""
    q = max(1, math.ceil(math.sqrt(abs(area) / (slew * raster**2))) - 1)
    while abs(area) / (q * raster) / (q * raster) > slew:
        q += 1
    return q, area / (q * raster)


def _closest(plan: dict, window: tuple[float, float]) -> tuple[float, float]:
    """
    Return ``(time, |k|)`` where the planned trajectory is nearest the origin within `window`.

    Exact on each knot interval: ``k`` is quadratic there, so ``|k|^2`` is a quartic whose
    stationary points are the real roots of a cubic.  The earliest of equal minima wins.
    """
    lo, hi = window
    times, values, areas = plan['times'], plan['values'], plan['areas']
    best_t, best_r = lo, math.inf
    for i in range(max(int(np.searchsorted(times, lo, side='right')) - 1, 0), len(times) - 1):
        if times[i] > hi:
            break
        span = float(times[i + 1] - times[i])
        a, b = complex(areas[i]), complex(values[i]) * span
        c = complex(values[i + 1] - values[i]) * span / 2.0
        u_lo, u_hi = max(0.0, (lo - times[i]) / span), min(1.0, (hi - times[i]) / span)
        if u_hi < u_lo:
            continue
        cubic = [2.0 * abs(c) ** 2, 3.0 * (c * b.conjugate()).real,
                 abs(b) ** 2 + 2.0 * (c * a.conjugate()).real, (b * a.conjugate()).real]
        roots = [r.real for r in np.roots(cubic) if abs(r.imag) < 1e-9 and u_lo < r.real < u_hi]
        for u in sorted([u_lo, *roots, u_hi]):
            r = abs(a + b * u + c * u * u)
            if r < best_r:
                best_t, best_r = float(times[i]) + u * span, r
    return best_t, best_r


class RosetteReadout(Module):
    """
    One rosette readout: ``petals`` whole petals of the analytic curve, at rest at both ends.

    Parameters
    ----------
    opts
        The scanner.  ``max_grad`` and ``max_slew`` bound the in-plane vector; the gradient, RF
        and ADC rasters, ``adc_dead_time``, ``adc_samples_divisor`` and a positive
        ``adc_samples_limit`` shape the ADC window.
    fov_mm, matrix
        The image.  ``matrix`` sets ``k_max = matrix / (2 fov)``, so it is a resolution rather
        than a sample count; ``dk = 1 / fov``.
    petals
        Whole petals, ``P``.  The core lasts ``P`` petal durations and returns to the origin.
    omega_ratio
        ``omega2 / omega1``, finite and at least zero.
    dwell_s
        ADC dwell, a whole ``adc_raster_time``.  No sample may move more than ``dk`` along the
        path: ``k_max max(omega1, omega2) dwell_s <= dk``.
    petal_duration_s
        ``pi / omega1``, a whole ``grad_raster_time``; ``None`` for the shortest petal on the
        raster that the gradient and slew limits and the departure budget allow.  The dwell does
        not enter that choice.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    k_max_per_m, dk_per_m : float
        ``matrix / (2 fov)`` and ``1 / fov``, in 1/m.
    petal_duration_s : float
        The petal played, given or derived.
    num_samples : int
        Samples in the one ADC event.
    core_samples : slice
        The samples whose times lie from the first ideal crossing to the last, inclusive.
    origin_crossing_ideal_s : tuple[float, ...]
        ``P + 1`` ideal crossing instants, seconds.
    origin_crossing_closest_s, origin_crossing_residual_per_m : tuple[float, ...]
        The designed trajectory's closest approach to the origin around each, and its distance
        from the origin there.
    origin_crossing_samples : tuple[int, ...]
        The sample nearest each ideal instant; the earlier of two equally near.
    departure_per_m : float
        An upper bound on the core's largest vector departure from the analytic curve, in 1/m,
        on the designed knots: the largest departure at 16 points per raster plus the most it can
        change between them.  It does not exceed ``B`` by more than that addition.
    """

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: float,
        matrix: int,
        petals: int,
        omega_ratio: float,
        dwell_s: float,
        petal_duration_s: float | None = None,
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.fov_mm = require_positive(fov_mm, 'fov_mm')
        self.matrix = require_count(matrix, 'matrix', low=2)
        self.petals = require_count(petals, 'petals', low=1)
        self.omega_ratio = self._check_ratio(omega_ratio)
        self.dwell_s = require_positive(dwell_s, 'dwell_s')
        self._limit = self._check_event_limit()
        self._on_raster(self.dwell_s, 'dwell_s', 'adc_raster_time')

        fov_m = self.fov_mm / 1e3
        self.dk_per_m = 1.0 / fov_m
        self.k_max_per_m = self.matrix / (2.0 * fov_m)
        self._raster = float(opts.grad_raster_time)
        self._gmax = float(opts.max_grad) * (1.0 - _GUARD)
        self._smax = float(opts.max_slew) * (1.0 - _GUARD)

        rasters = self._petal_rasters(petal_duration_s)
        self.petal_duration_s = Raster(self._raster).at(rasters)
        self._w1 = math.pi / self.petal_duration_s
        self._w2 = self.omega_ratio * self._w1
        self._check_sampling()

        plan = self._plan(self.petals, rasters)
        self._times, self._values, self._areas = plan['times'], plan['values'], plan['areas']
        self._core_start_s, self._core_end_s = plan['core']
        self.departure_per_m = self._departure(self.petals * rasters)
        self.origin_crossing_ideal_s = tuple(
            from_ticks(to_ticks(self._core_start_s) + p * to_ticks(self.petal_duration_s))
            for p in range(self.petals + 1)
        )
        closest = [_closest(plan, window) for window in plan['windows']]
        self.origin_crossing_closest_s = tuple(t for t, _ in closest)
        self.origin_crossing_residual_per_m = tuple(r for _, r in closest)

        self._adc = self._window(self._core_end_s, self.origin_crossing_closest_s[-1])
        self.num_samples = int(self._adc.num_samples)
        if self.num_samples > self._limit:
            self._refuse_samples(rasters)
        times = self._sample_ticks()
        ideal = [to_ticks(t) for t in self.origin_crossing_ideal_s]
        inside = np.flatnonzero((times >= ideal[0]) & (times <= ideal[-1]))
        self.core_samples = slice(int(inside[0]), int(inside[-1]) + 1)
        self.origin_crossing_samples = tuple(self._nearest(times, t) for t in ideal)

        end = max(to_ticks(float(self._times[-1])), to_ticks(float(pp.calc_duration(self._adc))))
        self._occupancy_s = Raster(self._raster).ceil(from_ticks(end))

    # ------------------------------------------------------------------ the design
    def _plan(self, petals: int, rasters: int) -> dict:
        """
        Return the knots of the whole readout at angle 0, and where its crossings are searched.

        Times are seconds from the block start, on raster edges; values are complex Hz/m, x the
        real part.  ``areas`` is the integral from the block start to each knot, in 1/m.
        """
        dt, k_max, w1, w2 = self._raster, self.k_max_per_m, self._w1, self._w2
        n_core = petals * rasters
        core_g = _g(np.arange(n_core + 1) * dt, k_max, w1, w2)
        g0, g_end = float(abs(core_g[0])), float(abs(core_g[-1]))
        u0, u_end = core_g[0] / g0, core_g[-1] / g_end
        residual = complex(np.sum(0.5 * (core_g[1:] + core_g[:-1])) * dt)

        n1, n2, a = _lead_in(g0, self._smax, dt)
        along = -float(np.real(residual * np.conj(u_end)))
        across = -float(np.real(residual * np.conj(1j * u_end)))
        m1, m2, b = _ramp_to_rest(g_end, along, self._smax, dt)
        q, h = _triangle(across, self._smax, dt)

        lead = n1 + n2
        out = lead + n_core
        index = np.concatenate([
            [0, n1], lead + np.arange(n_core + 1),
            out + np.array([m1, m1 + m2]), out + m1 + m2 + np.array([q, 2 * q]),
        ])
        values = np.concatenate([
            [0.0, -a * u0], core_g, [-b * u_end, 0.0], [h * 1j * u_end, 0.0],
        ])
        times = index * to_ticks(dt) / 1e12
        areas = np.concatenate([[0.0], np.cumsum(0.5 * (values[1:] + values[:-1]) * np.diff(times))])

        # Around an interior crossing the core is within B of the curve, whose |k| grows as
        # k_max |sin(w1 tau)|, so nothing nearer the origin than the crossing's own point lies
        # beyond |tau| = asin(2B / k_max) / w1.  Outside the core the search runs over the
        # monotone pieces: the lead-in after its gradient last passes zero, the lead-out before
        # it first does.
        bound = _departure_bound(k_max, w1, w2, dt, petals * self.petal_duration_s)
        reach = min(math.asin(min(1.0, 2.0 * bound / k_max)) / w1 + dt,
                    0.5 * self.petal_duration_s)
        t0, t_end = float(times[2]), float(times[2 + n_core])
        before = float(times[1]) + n2 * dt * a / (a + g0)
        after = t_end + m1 * dt * g_end / (g_end + b)
        crossings = [t0 + p * self.petal_duration_s for p in range(petals + 1)]
        windows = [(t - reach, t + reach) for t in crossings]
        windows[0] = (before, windows[0][1])
        windows[-1] = (windows[-1][0], after)
        return {'times': times, 'values': values, 'areas': areas, 'core': (t0, t_end),
                'windows': windows}

    def _departure(self, intervals: int) -> float:
        """
        Return an upper bound on ``max |k_r - k|`` over the core, from the designed knots.

        The largest departure on ``_DEPARTURE_POINTS`` points per raster, plus half their spacing
        times the largest slope the departure can have, ``|g_r - g| <= dt^2 M3 / 8``.
        """
        dt, k_max, w1, w2 = self._raster, self.k_max_per_m, self._w1, self._w2
        s = np.linspace(0.0, dt, _DEPARTURE_POINTS + 1)
        worst = 0.0
        for lo in range(0, intervals, 4096):
            j = np.arange(lo, min(lo + 4096, intervals))[:, None]
            g_a, g_b = self._values[2 + j], self._values[3 + j]
            k_r = self._areas[2 + j] + g_a * s + (g_b - g_a) * s * s / (2.0 * dt)
            worst = max(worst, float(np.max(np.abs(k_r - _k(j * dt + s, k_max, w1, w2)))))
        slope = dt**2 / 8.0 * _derivative_bound(k_max, w1, w2, 3)
        return worst + 0.5 * dt / _DEPARTURE_POINTS * slope

    # ------------------------------------------------------------------ the ADC
    def _window(self, core_end_s: float, closest_end_s: float) -> object:
        """
        Return the ADC: one dwell of padding each side, put on its lattice, then checked.

        Times are integer ticks.  Sample ``j`` is at ``delay + (j + 1/2) dwell``, so doubled ticks
        keep half a dwell exact.
        """
        opts = self.opts
        dwell, rf = to_ticks(self.dwell_s), to_ticks(float(opts.rf_raster_time))
        divisor = int(getattr(opts, 'adc_samples_divisor', 1) or 1)
        t0, t_end = to_ticks(self._core_start_s), to_ticks(core_end_s)
        start = 2 * min(t0, to_ticks(self.origin_crossing_closest_s[0]))
        finish = 2 * max(t_end, to_ticks(closest_end_s))

        def count(delay: int) -> int:
            n = -(-(t_end - delay) // dwell) + 1
            return n + (-n) % divisor

        delay = (t0 - dwell) // rf * rf
        n = count(delay)
        while True:
            adc = pp.make_adc(num_samples=n, dwell=self.dwell_s, delay=from_ticks(delay),
                              system=opts)
            stored = to_ticks(float(adc.delay))
            if 2 * stored + dwell > start:
                if stored > delay:
                    self._refuse_late_window(stored)
                delay = (delay - dwell) // rf * rf
                n = max(n, count(delay))
            elif count(stored) > n:
                n = count(stored)
            elif 2 * stored + (2 * n - 1) * dwell < finish:
                n += divisor
            else:
                return adc

    def _sample_ticks(self) -> np.ndarray:
        """Each sample's time on the block clock, in integer ticks."""
        delay, dwell = to_ticks(float(self._adc.delay)), to_ticks(self.dwell_s)
        return delay + (2 * np.arange(self.num_samples, dtype=np.int64) + 1) * dwell // 2

    def _nearest(self, times: np.ndarray, at: int) -> int:
        """The sample nearest the instant `at`, in ticks; the earlier of two equally near."""
        dist = np.abs(times - at)
        return int(np.flatnonzero(dist == dist.min())[0])

    # ------------------------------------------------------------------ what it reports
    def sample_times_s(self) -> np.ndarray:
        """Return each ADC sample's time, in seconds from the start of the block."""
        return self._sample_ticks() / 1e12

    def k_per_m(self, *, angle_rad: float = 0.0) -> np.ndarray:
        """
        Return ``(2, num_samples)`` k-space positions of the samples, in 1/m.

        Integrated exactly from the designed gradient knots at :meth:`sample_times_s`, and
        rotated by `angle_rad` as :meth:`build` rotates the events.  Samples after the gradient
        has ended sit at its end point.
        """
        k = self._k_at(self.sample_times_s()) * np.exp(1j * float(angle_rad))
        return np.vstack([k.real, k.imag])

    def _k_at(self, t: np.ndarray) -> np.ndarray:
        """Designed k at block times `t`, complex: exact on the piecewise-linear knots."""
        times, values, areas = self._times, self._values, self._areas
        i = np.clip(np.searchsorted(times, t, side='right') - 1, 0, len(times) - 2)
        span = times[i + 1] - times[i]
        s = np.clip(t - times[i], 0.0, span)
        return areas[i] + values[i] * s + (values[i + 1] - values[i]) * s * s / (2.0 * span)

    def limits(self) -> tuple[float, float]:
        """
        Return the peak in-plane vector gradient and slew this readout emits, Hz/m and Hz/m/s.

        Over the knots and the segments between them, lead-in and lead-out included.  Both are
        the bound for each axis at every rotation.
        """
        slew = np.diff(self._values) / np.diff(self._times)
        return float(np.max(np.abs(self._values))), float(np.max(np.abs(slew)))

    # ------------------------------------------------------------------ assembly
    def build(self, *, angle_rad: float = 0.0, acquire: bool = True,
              phase_deg: float = 0.0) -> LogicBlock:
        """
        Return the readout: one extended trapezoid per in-plane axis, the ADC, and its span.

        Parameters
        ----------
        angle_rad
            In-plane rotation, applied to the emitted gradients, so the block is the rotated
            readout.  An axis whose gradient is identically zero is omitted.
        acquire
            ``False`` plays the gradients without the ADC, for a dummy shot.  The block lasts as
            long either way.
        phase_deg
            ADC phase offset, for the RF spoiling a caller schedules.
        """
        out = LogicBlock()
        values = self._values * np.exp(1j * float(angle_rad))
        for axis, amplitudes in (('x', values.real), ('y', values.imag)):
            if np.any(amplitudes != 0.0):
                out.add(0.0, pp.make_extended_trapezoid(axis, times=self._times,
                                                        amplitudes=amplitudes, system=self.opts))
        if acquire:
            phase = math.radians(float(phase_deg))
            out.add(0.0, derive(self._adc, phase_offset=phase) if phase else self._adc)
        out.add(0.0, pp.make_delay(self._occupancy_s))
        return out

    # ------------------------------------------------------------------ the domain
    def _check_ratio(self, ratio: float) -> float:
        """Return `ratio` as a float, having checked it is finite and at least zero."""
        value = float(ratio)
        if not math.isfinite(value) or value < 0.0:
            msg = format_error(
                f'omega_ratio must be a finite number >= 0, got {value:g}.',
                {'omega_ratio': value},
                ['0 is a spike through the centre, below 1 narrow petals, 1 circles through the '
                 'origin, above 1 wide petals',
                 'a negative ratio is the mirror image of a positive one, and it is not supported'],
            )
            raise ConfigurationError(msg)
        return value

    def _check_event_limit(self) -> int:
        """Return ``opts.adc_samples_limit``, refusing a limit that bounds nothing."""
        limit = int(getattr(self.opts, 'adc_samples_limit', 0) or 0)
        if limit <= 0:
            msg = format_error(
                f'opts.adc_samples_limit is {limit}, so there is no per-event sample limit for '
                f"this readout's one ADC event to fit.",
                {'adc_samples_limit': limit},
                ['set it to the interpreter\'s per-event limit, e.g. pp.Opts(..., '
                 'adc_samples_limit=8192)',
                 'this readout never splits its window, so the limit is what bounds it'],
            )
            raise ConfigurationError(msg)
        return limit

    def _on_raster(self, value: float, name: str, raster_name: str) -> None:
        """Refuse `value` unless it is a whole number of ``opts.<raster_name>``."""
        raster = Raster(float(getattr(self.opts, raster_name)))
        if raster.holds(value):
            return
        below, above = raster.floor(value), raster.ceil(value)
        msg = format_error(
            f'{name} = {value * 1e6:g} us is not a whole {raster_name} ({raster.dt * 1e6:g} us).',
            {name: value, raster_name: raster.dt},
            [f'pass {below * 1e6:g} us or {above * 1e6:g} us'],
        )
        raise ConfigurationError(msg)

    def _peaks(self, rasters: int) -> tuple[float, float, float]:
        """Return the analytic peak |g|, the peak slew and ``B`` for a petal of `rasters`."""
        w1 = math.pi / (rasters * self._raster)
        w2 = self.omega_ratio * w1
        bound = _departure_bound(self.k_max_per_m, w1, w2, self._raster,
                                 self.petals * rasters * self._raster)
        return (self.k_max_per_m * w1 * max(1.0, self.omega_ratio),
                self.k_max_per_m * w1**2 * (1.0 + self.omega_ratio**2), bound)

    def _fits(self, rasters: int) -> bool:
        """Whether a petal of `rasters` meets the gradient and slew limits and the budget."""
        grad, slew, bound = self._peaks(rasters)
        return grad <= self._gmax and slew <= self._smax and bound <= _BUDGET_DK * self.dk_per_m

    def _shortest(self) -> int:
        """The fewest rasters per petal that fit.  Fitting is monotone in the petal."""
        hi = 1
        while not self._fits(hi):
            hi *= 2
        lo = hi // 2
        while hi - lo > 1:
            mid = (lo + hi) // 2
            lo, hi = (lo, mid) if self._fits(mid) else (mid, hi)
        return hi

    def _petal_rasters(self, petal_duration_s: float | None) -> int:
        """Return rasters per petal: the shortest that fits, or the caller's, checked."""
        shortest = self._shortest()
        if petal_duration_s is None:
            return shortest
        petal = require_positive(petal_duration_s, 'petal_duration_s')
        self._on_raster(petal, 'petal_duration_s', 'grad_raster_time')
        rasters = Raster(self._raster).count(petal)
        grad, slew, bound = self._peaks(rasters)
        way_out = (f'pass petal_duration_s >= {shortest * self._raster * 1e6:g} us, or omit it '
                   f'for the shortest petal that fits')
        fields = {'petal_duration_s': petal, 'k_max_per_m': self.k_max_per_m,
                  'omega_ratio': round(self.omega_ratio, 4)}
        gamma = float(self.opts.gamma)
        vector = 'the peak is an in-plane vector, and some rotation puts all of it on one axis'
        geometry = 'or lower matrix, or widen fov_mm (k_max = matrix / (2 fov))'
        if grad > self._gmax:
            msg = format_error(
                f'a {petal * 1e6:g} us petal needs {grad / gamma * 1e3:.4g} mT/m at its peak, '
                f'above max_grad = {float(self.opts.max_grad) / gamma * 1e3:.4g} mT/m; {vector}.',
                fields,
                [way_out, geometry, *(['or lower omega_ratio'] if self.omega_ratio > 1 else [])],
            )
            raise ConfigurationError(msg)
        if slew > self._smax:
            msg = format_error(
                f'a {petal * 1e6:g} us petal needs {slew / gamma:.4g} T/m/s at its peak, above '
                f'max_slew = {float(self.opts.max_slew) / gamma:.4g} T/m/s; {vector}.',
                fields,
                [way_out, geometry, *(['or lower omega_ratio'] if self.omega_ratio > 0 else [])],
            )
            raise ConfigurationError(msg)
        if bound > _BUDGET_DK * self.dk_per_m:
            msg = format_error(
                f'the raster departure bound for a {petal * 1e6:g} us petal is '
                f'{bound / self.dk_per_m:.4g} dk, above the {_BUDGET_DK:g} dk this readout '
                f'supports.  The bound is conservative, so the petal may depart by less: this is '
                f'the edge of the supported domain, not a physical limit.',
                {**fields, 'grad_raster_time': self._raster, 'matrix': self.matrix},
                [way_out, 'or lower matrix'],
            )
            raise ConfigurationError(msg)
        return rasters

    def _check_sampling(self) -> None:
        """Refuse a dwell that moves k by more than dk where the trajectory is fastest."""
        speed = self.k_max_per_m * self._w1 * max(1.0, self.omega_ratio)
        step = speed * self.dwell_s
        if step <= self.dk_per_m * (1.0 + _GUARD):
            return
        longest = Raster(float(self.opts.adc_raster_time)).floor(self.dk_per_m / speed)
        fixes = [f'pass dwell_s <= {longest * 1e6:g} us'] if longest > 0 else []
        fixes.append('or pass a longer petal_duration_s, which slows the trajectory')
        msg = format_error(
            f'dwell_s = {self.dwell_s * 1e6:g} us moves k by {step / self.dk_per_m:.3g} dk per '
            f'sample where the trajectory is fastest; at most 1 dk is allowed.',
            {'dwell_s': self.dwell_s, 'petal_duration_s': self.petal_duration_s},
            fixes,
        )
        raise ConfigurationError(msg)

    def _refuse_late_window(self, stored: int) -> None:
        """Refuse a window that cannot open before the first crossing."""
        dead = float(self.opts.adc_dead_time)
        longest = Raster(float(self.opts.adc_raster_time)).floor(
            2.0 * (self._core_start_s - dead))
        late = from_ticks(stored) + 0.5 * self.dwell_s - self._core_start_s
        fixes = ([f'pass dwell_s <= {longest * 1e6:g} us'] if longest > 0 else
                 ['no dwell fits: the lead-in is shorter than adc_dead_time'])
        msg = format_error(
            f'with dwell_s = {self.dwell_s * 1e6:g} us the first sample falls '
            f'{late * 1e6:.3g} us after the first crossing: pp.make_adc holds an ADC delay to at '
            f'least adc_dead_time ({dead * 1e6:g} us), and the readout starts with the lead-in.',
            {'dwell_s': self.dwell_s, 'adc_dead_time': dead,
             'lead-in': f'{self._core_start_s * 1e6:g} us'},
            fixes,
        )
        raise ConfigurationError(msg)

    def _refuse_samples(self, rasters: int) -> None:
        """Refuse a window above the per-event limit, naming the most petals that fit."""
        # Each petal adds petal / dwell samples, so start from that estimate and step down until
        # the window that the same rule builds fits.
        per_petal = self.petal_duration_s / self.dwell_s
        fewer = min(self.petals - 1,
                    self.petals + math.floor((self._limit - self.num_samples) / per_petal) + 1)
        while fewer >= 1:
            plan = self._plan(fewer, rasters)
            end = _closest(plan, plan['windows'][-1])[0]
            if int(self._window(plan['core'][1], end).num_samples) <= self._limit:
                break
            fewer -= 1
        speed = self.k_max_per_m * self._w1 * max(1.0, self.omega_ratio)
        longest = Raster(float(self.opts.adc_raster_time)).floor(self.dk_per_m / speed)
        fixes = [f'lower petals to at most {fewer}'] if fewer >= 1 else []
        if longest > self.dwell_s:
            fixes.append(f'or lengthen dwell_s, up to the {longest * 1e6:g} us along-path limit')
        fixes.append('this readout acquires in one ADC event and does not split it')
        msg = format_error(
            f'{self.petals} petals need {self.num_samples} samples in one ADC event, above '
            f'adc_samples_limit = {self._limit}.',
            {'petals': self.petals, 'dwell_s': self.dwell_s, 'adc_samples_limit': self._limit},
            fixes,
        )
        raise ConfigurationError(msg)
