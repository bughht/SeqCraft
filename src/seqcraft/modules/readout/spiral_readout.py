r"""
:class:`SpiralReadout` -- one reversible spiral arm, and the four ways to traverse it.

``readout/`` because it contains ADCs and no RF.

What makes this one module rather than four
-------------------------------------------
Everything below rests on a single design decision: **an arm begins and ends at rest.**

.. math:: g(\text{start}) = 0 \qquad g(\text{end}) = 0

An arm at rest at both ends can be played backwards, so ``'in'`` is ``'out'`` reversed; two arms
can be laid end to end with no connector, because both are at :math:`g = 0` where they meet.  The
four variants are then a choice of *how many arms and in which order*, not four implementations
behind a string flag.

The independent reference -- ``pulseq/pulseq``'s ``writeSpiral.m`` -- does the opposite: it ends
its spiral-out at full gradient and ramps down inside the spoiler.  That is coherent for a family
of one, and it is why that reference has **no spiral-in**: an arm ending at full gradient
time-reverses into one starting there, which no block can begin with.

**What the policy costs, measured** (``tools/module_mining/candidates/spiral/run_endpoint_cost.py``,
seven protocols across two hardware regimes): 0.097--0.137 ms, which is 0.20--0.78 % of readout
duration.  The penalty is a *braking distance* and so is fixed in absolute terms -- it shrinks as
readouts lengthen, and its worst case is the shortest readout.  k-space extent is identical and
peak gradient is unchanged or lower, because braking happens where a spiral is fastest.

Three passes, and only the middle one knows the scanner
-------------------------------------------------------
::

    path         the Nyquist condition and nothing else:  dtheta/dr = 2 pi FOV(r) / shots
      |          no gradient limit, no slew limit, no raster
      v
    traversal    time-optimal under |g| <= max_grad and |dg/dt| <= max_slew,
      |          swept forward from rest and backward to rest
      v
    realization  raster waveform, ADC segmentation, prephaser, rewinder

That separation is **internal**.  It is not two public modules: a path is meaningful without a
scanner, but no traversal is meaningful without its path, and nothing wants to swap traversal
policies on a fixed path.  ``writeSpiral.m`` arrives at the same separation independently, which
is the strongest single piece of corroboration this module has.

What it reports, and why the plural matters
-------------------------------------------
``'out-in'`` crosses the origin **twice** and every other variant once.  A module whose semantic
echo were a scalar would work for three variants and change meaning on the fourth, so
:attr:`origin_crossing_samples` is a tuple from the first line of code -- even for ``'out'``,
whose tuple always has length one.

An **origin crossing is a fact about the trajectory**, not about magnetisation.  This module owns
where and when the trajectory passes through :math:`k = 0`; the kernel above it owns whether the
physical echo is aligned with the intended crossing.  A readout that claimed to own TE would hide
exactly that error.

What is deliberately not here
-----------------------------
``echoes > 1`` refuses, naming the contract: multi-echo is the last step of the recorded
implementation order and needs the transition rules -- a fly-back for one-arm variants, continuous
traversal for two-arm ones -- which nothing yet exercises.  ``'out-in'`` already produces two
crossings and a derived spacing at ``echoes=1``, so the plural machinery is exercised, not merely
declared.

No 3D, no anisotropic field of view, no shipped density presets, and no ordering tables:
``density`` is data and ``angle_rad`` is an argument.

Examples
--------
>>> import pypulseq as pp
>>> import seqcraft as sc
>>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
>>> spiral = sc.modules.SpiralReadout(opts=opts, fov_mm=220.0, matrix=64, shots=4,
...                                   dwell_s=4e-6)
>>> spiral.origin_crossing_samples
(0,)
>>> out_in = sc.modules.SpiralReadout(opts=opts, fov_mm=220.0, matrix=64, shots=4,
...                                   dwell_s=4e-6, variant='out-in')
>>> len(out_in.origin_crossing_samples)
2
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ...design.events import derive
from ...design.logic import LogicBlock, barrier
from ...design.module import Module
from ...errors import ConfigurationError, format_error
from .._support import ceil_raster, require_count, require_positive

if TYPE_CHECKING:
    from pypulseq.opts import Opts

__all__ = ['SpiralReadout']

#: How many arms each variant plays, and in which direction.  ``+1`` is outward, ``-1`` inward.
#: This table *is* the variant contract: everything else about a variant follows from it.
_VARIANTS: dict[str, tuple[int, ...]] = {
    'out': (+1,),
    'in': (-1,),
    'in-out': (-1, +1),
    'out-in': (+1, -1),
}

#: Where path sampling starts.  It has to stay finer than the gradient raster the traversal is
#: resampled onto -- otherwise the resampling *upsamples*, interpolating a coarse path into kinks
#: that read as slew violations no margin can fix.  A long arm needs more, so :meth:`_design`
#: refines against the raster count it actually produced rather than guessing from the protocol.
_PATH_STEPS = 6000

#: Path knots per gradient raster point after refinement.
_PATH_OVERSAMPLING = 24

#: Where the design margin starts.  A time-optimal sweep sits exactly on the constraint by
#: construction, and resampling onto the gradient raster rounds some segments over it --
#: ``writeSpiral.m`` carries the same margin for the same reason, calling it "otherwise we just
#: about violate the slew rate due to the rounding errors".
#:
#: But a fixed margin is a guess.  The contract is about the **emitted** waveform, so this is a
#: starting point for :meth:`SpiralReadout._design`, which measures the rastered result and
#: tightens until it actually fits.  Inferring the slew from the algorithm is what §9 of the
#: acceptance plan forbids.
_MARGIN_START = 0.97

#: Multiplied into the margin on each retry, and the number of retries before giving up.
_MARGIN_STEP = 0.93
_MARGIN_TRIES = 12

#: Below this fraction of `dk`, a sample counts as the origin.  The crossing is found by
#: interpolation, not by nearest-sample search, for the reason `RadialReadout` records: a float
#: tie on `argmin|k|` is how a comparator gets the echo wrong.
_ORIGIN_FRACTION = 0.5

#: Closer than this fraction of `dk` to the origin counts as being at it, for deciding whether a
#: prephaser or a rewinder is required at all.
_AT_ORIGIN = 1e-3


def _split(total: int, segments: int, divisor: int) -> tuple[int, ...]:
    """Divide `total` samples into `segments` events, each a whole multiple of `divisor`."""
    per = total // segments
    per -= per % divisor
    spare = total - per * segments
    counts = [per] * segments
    counts[-1] += spare - spare % divisor
    return tuple(counts)


class SpiralReadout(Module):
    """
    A variable-density spiral arm, traversed one of four ways.

    Parameters
    ----------
    opts
        The scanner.  ``max_grad`` and ``max_slew`` shape the traversal; nothing else does.
    fov_mm, matrix
        The image.  ``matrix`` sets ``k_max = matrix / (2 * fov)``, so it is a resolution rather
        than a sample count -- a spiral has no grid to count against.
    shots
        Interleaves.  **This enters the path**, not the schedule: the Nyquist condition is
        ``dtheta/dr = 2 pi FOV(r) / shots``, so one shot of a four-shot spiral is a different
        curve from one shot of a one-shot spiral, not the same curve rotated.
    dwell_s
        ADC dwell.
    density
        Sampled field-of-view multipliers across the radius, ``(1.0,)`` for uniform and
        ``(1.0, 0.5)`` for half the field of view at the edge.  **Sampled values, not polynomial
        coefficients**: under the coefficient convention that `vds.m` uses, ``(1.0, 0.5)`` means
        the field of view *grows* by half -- an oversampled outer ring and a longer readout, with
        nothing raising.  Every number written here is the thing it controls.
    variant
        ``'out'``, ``'in'``, ``'in-out'`` or ``'out-in'``.
    echoes
        Reserved; only ``1`` is implemented.  See the module docstring.
    adc_segments
        Force the acquisition into this many ADC events, or ``None`` to use the fewest that fit
        ``opts.adc_samples_limit``.  ``1`` refuses rather than splitting, and names the smallest
        ``shots`` that would fit.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    origin_crossing_samples : tuple[int, ...]
    origin_crossing_times : tuple[float, ...]
    k_max_per_m, dk_per_m, num_samples, duration_s : float
    """

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: float,
        matrix: int,
        shots: int,
        dwell_s: float,
        density: tuple[float, ...] = (1.0,),
        variant: str = 'out',
        echoes: int = 1,
        adc_segments: int | None = None,
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.fov_mm = require_positive(fov_mm, 'fov_mm')
        self.matrix = require_count(matrix, 'matrix', low=2)
        self.shots = require_count(shots, 'shots', low=1)
        self.dwell_s = require_positive(dwell_s, 'dwell_s')
        self.density = self._check_density(density)
        self.variant = self._check_variant(variant)
        self.echoes = self._check_echoes(echoes)

        fov_m = self.fov_mm / 1e3
        self.dk_per_m = 1.0 / fov_m
        self.k_max_per_m = self.matrix / (2.0 * fov_m)

        # --- passes 1 to 3.  The design loop validates the COMPOSED waveform, because that is
        # what plays: a two-arm variant has a seam the single arm does not, and checking the arm
        # alone would certify something nobody emits.
        self._design()
        self._knot_cache: tuple | None = None

        self.duration_s = float(self._waveform.shape[1] - 1) * float(opts.grad_raster_time)
        self.num_samples, self._segments = self._plan_adc(adc_segments)
        # Build the ADC events now, because pypulseq gives an ADC a LEADING delay of one
        # `adc_dead_time` as well as a trailing one.  Assuming the first sample sits at
        # `dwell/2` puts every reported sample time -- and therefore the whole reported
        # trajectory -- one dead time early, which is a blur and a wrong field map that no
        # k-space check would show.  Read it off the event instead of deriving it.
        self._adcs = tuple(
            pp.make_adc(num_samples=count, dwell=self.dwell_s, system=opts)
            for count in self._segments
        )
        self._adc_starts = self._place_adcs()
        self.origin_crossing_samples, self.origin_crossing_times = self._crossings()

    # ------------------------------------------------------------------ pass 1: the path
    def _path(self, steps: int = _PATH_STEPS) -> np.ndarray:
        r"""
        Return ``(2, n)`` k-space points from the origin to ``k_max``, in 1/m.

        The Nyquist condition and nothing else.  With ``FOV(r)`` sampled from `density`,

        .. math:: \frac{d\theta}{dr} = \frac{2 \pi\, \mathrm{FOV}(r)}{\text{shots}}

        integrated by cumulative trapezoid.  A uniform density recovers the Archimedean spiral
        ``writeSpiral.m`` builds.
        """
        radius = np.linspace(0.0, self.k_max_per_m, steps)
        u = radius / self.k_max_per_m
        fov_m = (self.fov_mm / 1e3) * np.interp(
            u, np.linspace(0.0, 1.0, len(self.density)), np.asarray(self.density, dtype=float),
        )
        dtheta_dr = 2.0 * np.pi * fov_m / self.shots
        theta = np.concatenate([[0.0], np.cumsum(
            0.5 * (dtheta_dr[1:] + dtheta_dr[:-1]) * np.diff(radius),
        )])
        k = np.vstack([radius * np.cos(theta), radius * np.sin(theta)])

        # Resample uniformly in ARCLENGTH, not in radius.  A spiral's outer turns are far longer
        # than its inner ones, so a radius-uniform grid starves exactly the region where the
        # traversal brakes -- and the resampling onto the gradient raster then interpolates a
        # coarse curve into kinks that read as slew violations no design margin can fix.
        step = np.hypot(*np.diff(k, axis=1))
        arclength = np.concatenate([[0.0], np.cumsum(step)])
        even = np.linspace(0.0, arclength[-1], steps)
        return np.vstack([np.interp(even, arclength, k[0]), np.interp(even, arclength, k[1])])

    # -------------------------------------------------------------- pass 2: the traversal
    def _design(self) -> None:
        """
        Solve the traversal and compose the acquisition, **measuring** that what plays fits.

        The traversal solves against a margin, and resampling onto the raster can still round a
        segment over -- so the composed result is checked and the margin tightened until it fits.
        Verifying the thing that plays rather than trusting the solver is the point; it is also
        what makes :meth:`limits` a measurement rather than a restatement.
        """
        # Refine the path until it is comfortably finer than the raster the traversal lands on.
        probe = self._traverse(self._path(), _MARGIN_START)
        path = self._path(max(_PATH_STEPS, _PATH_OVERSAMPLING * probe.shape[1]))

        margin = _MARGIN_START
        for _ in range(_MARGIN_TRIES):
            self._arm_k = self._traverse(path, margin)
            self._waveform = self._compose()
            gradient, slew = self._limits_of(self._waveform)
            if gradient <= float(self.opts.max_grad) and slew <= float(self.opts.max_slew):
                self._margin = margin
                return
            margin *= _MARGIN_STEP
        msg = format_error(
            'no traversal of this path fits the scanner after tightening the design margin.',
            {'fov_mm': self.fov_mm, 'matrix': self.matrix, 'shots': self.shots,
             'variant': self.variant, 'max_grad': float(self.opts.max_grad),
             'max_slew': float(self.opts.max_slew)},
            ['raise shots, which shortens and straightens each arm',
             'or lower matrix, or widen fov_mm'],
        )
        raise ConfigurationError(msg)

    def _limits_of(self, k: np.ndarray) -> tuple[float, float]:
        """
        Return the peak vector gradient and slew of the waveform `k` implies, **as emitted**.

        On the emitted knots, not on raster-spaced differences.  An arbitrary gradient's samples
        sit at raster *centres* while its enforced zeros sit at the block edges, so the first and
        last intervals are **half a raster** -- and the slew across them is twice what a
        full-raster difference reports.  Measuring the convenient lattice instead of the emitted
        one is how a design passes its own check and is then refused by the factory that builds it.
        """
        raster = float(self.opts.grad_raster_time)
        waveform = np.diff(k, axis=1) / raster
        if waveform.shape[1] < 1:
            return 0.0, 0.0
        centres = (np.arange(waveform.shape[1]) + 0.5) * raster
        times = np.concatenate([[0.0], centres, [waveform.shape[1] * raster]])
        amplitudes = np.hstack([np.zeros((2, 1)), waveform, np.zeros((2, 1))])
        slew = np.diff(amplitudes, axis=1) / np.diff(times)
        return float(np.max(np.hypot(*waveform))), float(np.max(np.hypot(*slew)))

    def limits(self) -> tuple[float, float]:
        """
        Return the peak vector gradient and slew this readout actually emits, in Hz/m and Hz/m/s.

        Measured off the composed waveform, not asserted from the traversal's constraints.
        """
        return self._limits_of(self._waveform)

    def _traverse(self, k: np.ndarray, margin: float) -> np.ndarray:
        """
        Return ``(times, k)`` on the gradient raster, at rest at both ends.

        Time-optimal under the two limits, by the standard forward--backward sweep over
        arclength: the curvature ceiling first -- a bend spends slew budget before any
        acceleration does -- then what is left for tangential acceleration.  Both sweeps start
        from ``v = 0``, and *that* is the endpoint policy.
        """
        d1 = np.gradient(k, axis=1)
        step = np.hypot(d1[0], d1[1])
        arclength = np.concatenate([[0.0], np.cumsum(step[1:])])
        d2 = np.gradient(d1, axis=1)
        curvature = np.nan_to_num(
            np.abs(d1[0] * d2[1] - d1[1] * d2[0]) / np.maximum(step, 1e-30) ** 3,
            nan=0.0, posinf=0.0,
        )

        g_max = float(self.opts.max_grad) * margin
        s_max = float(self.opts.max_slew) * margin
        ceiling = np.minimum(g_max, np.sqrt(s_max / np.maximum(curvature, 1e-30)))
        ds = np.diff(arclength)

        def sweep(reverse: bool) -> np.ndarray:
            """One monotone pass from rest, in the given direction."""
            v = ceiling.copy()
            order = range(len(ds) - 1, -1, -1) if reverse else range(len(ds))
            v[-1 if reverse else 0] = 0.0
            for i in order:
                here, nxt = (i + 1, i) if reverse else (i, i + 1)
                lateral = (v[here] ** 2) * curvature[here]
                tangential = np.sqrt(max(s_max**2 - lateral**2, 0.0))
                v[nxt] = min(ceiling[nxt],
                             np.sqrt(max(v[here] ** 2 + 2.0 * tangential * ds[i], 0.0)))
            return v

        speed = np.minimum(sweep(reverse=False), sweep(reverse=True))
        dt = ds / np.maximum(0.5 * (speed[:-1] + speed[1:]), 1e-12)
        t = np.concatenate([[0.0], np.cumsum(dt)])

        raster = float(self.opts.grad_raster_time)
        # Round UP and add one knot past the end.  Rounding down truncates the braking ramp, so
        # the waveform would stop at whatever speed it had reached and `last=0` would then force
        # a one-raster jump to zero -- a slew violation manufactured entirely by the resampling,
        # at the exact place the endpoint policy is supposed to be arriving gently.  `np.interp`
        # holds the final value past `t[-1]`, so the extra knot carries zero gradient.
        n = int(np.ceil(t[-1] / raster)) + 2
        grid = np.arange(n) * raster
        return np.vstack([np.interp(grid, t, k[0]), np.interp(grid, t, k[1])])

    # ----------------------------------------------------------- pass 3: the realization
    def _compose(self) -> np.ndarray:
        """
        Return the whole acquisition's k, by laying arms end to end in the variant's order.

        Every arm is the same arm; ``-1`` plays it backwards.  Consecutive arms meet at
        :math:`g = 0`, so the concatenation is continuous with no connector -- which is the
        property the endpoint policy exists to provide.
        """
        pieces: list[np.ndarray] = []
        for index, direction in enumerate(_VARIANTS[self.variant]):
            arm = self._arm_k if direction > 0 else self._arm_k[:, ::-1]
            pieces.append(arm if index == 0 else arm[:, 1:])   # drop the shared knot
        return np.hstack(pieces)

    def _plan_adc(self, forced: int | None) -> tuple[int, tuple[int, ...]]:
        """
        Return the total sample count and the per-segment split.

        One acquisition region per contiguous non-zero gradient -- which for every variant here is
        the whole readout, because the arms join at rest rather than stopping.  The split is the
        interpreter's per-event sample limit and nothing else, so a segment boundary may fall
        anywhere, **including on an origin crossing**; the contract requires only that the
        crossing be inside a sampling window.

        **The acquisition must finish before the gradient does.**  Each ADC event spends a dead
        time before and after sampling, so a segmented readout fits fewer samples than its
        duration divided by the dwell -- and an acquisition that outlasts its gradient makes the
        compiler hold the last block open, padding the waveform with area nobody designed.  So
        the budget is solved for the segment count rather than guessed.
        """
        divisor = int(getattr(self.opts, 'adc_samples_divisor', 1) or 1)
        limit = int(getattr(self.opts, 'adc_samples_limit', 0) or 0)
        raster = float(self.opts.grad_raster_time)

        probe = pp.make_adc(num_samples=divisor, dwell=self.dwell_s, system=self.opts)
        lead, trail = float(probe.delay), float(probe.dead_time)

        def placed_end(split: tuple[int, ...]) -> float:
            """Where the last ADC event finishes once every span is on the gradient raster."""
            return float(sum(ceil_raster(lead + count * self.dwell_s + trail, raster)
                             for count in split))

        def budget(count: int) -> int:
            """
            Samples that fit **inside** the gradient when split into `count` events.

            Sized against the placement, not against an estimate of it.  Each event spends a lead
            delay and a trailing dead time *and* is rounded up to the gradient raster, so the
            overhead per segment is up to ``lead + trail + raster`` rather than the ``lead +
            trail`` a subtraction would guess.  Under-estimating it by one raster is enough: the
            acquisition then outlasts the gradient, the compiler holds the last block open, and
            the waveform is padded with area nobody designed -- which is exactly the m0 the
            contract check refuses.

            The estimate is deliberately conservative, so the search has to run **both ways**.
            Shrinking alone leaves whole samples unacquired at the tail, and for ``'in'`` the tail
            is where the origin is: an arm that brakes to rest at k=0 puts its last dk in its last
            few microseconds, so discarding 20 us of acquisition can discard the echo.
            `placed_end` is non-decreasing in the sample count, so growing while it still fits
            terminates, and lands on the largest acquisition the raster admits.
            """
            samples = int((self.duration_s - count * (lead + trail + raster)) / self.dwell_s)
            samples -= samples % divisor
            while samples > 0 and placed_end(_split(samples, count, divisor)) > self.duration_s:
                samples -= divisor
            while placed_end(_split(samples + divisor, count, divisor)) <= self.duration_s:
                samples += divisor
            return samples

        segments = 1 if forced is None else require_count(forced, 'adc_segments', low=1)
        if forced is None and limit > 0:
            while budget(segments) > segments * limit and segments < 1024:
                segments += 1

        total = budget(segments)
        if total <= 0:
            msg = format_error(
                'the readout is too short to hold one ADC sample block.',
                {'duration_ms': self.duration_s * 1e3, 'dwell_us': self.dwell_s * 1e6,
                 'adc_segments': segments, 'adc_samples_divisor': divisor},
                ['lower dwell_s', 'or raise matrix, which lengthens the arm'],
            )
            raise ConfigurationError(msg)

        if limit > 0 and total > segments * limit:
            needed = segments
            while budget(needed) > needed * limit and needed < 1024:
                needed += 1
            msg = format_error(
                f'{total} samples will not fit in {segments} ADC event(s) of at most {limit}.',
                {'samples': total, 'adc_segments': segments, 'adc_samples_limit': limit,
                 'shots': self.shots},
                [f'pass adc_segments={needed}',
                 f'or raise shots to about {self.shots * needed}, which shortens each arm'],
            )
            raise ConfigurationError(msg)

        split = _split(total, segments, divisor)
        if min(split) <= 0:
            msg = format_error(
                f'{segments} ADC events cannot each hold a multiple of {divisor} samples.',
                {'samples': total, 'adc_segments': segments, 'adc_samples_divisor': divisor},
                ['lower adc_segments', 'or adjust dwell_s so the total divides evenly'],
            )
            raise ConfigurationError(msg)
        return sum(split), split

    # ------------------------------------------------------------------ what it knows
    def _crossings(self) -> tuple[tuple[int, ...], tuple[float, ...]]:
        """
        Return the ADC samples and times at which the trajectory passes through the origin.

        Interpolated, never a nearest-sample search.  A tuple **always** -- ``'out-in'`` has two,
        and a scalar here would change meaning when it arrived.
        """
        times = self._arm_sample_times()
        k = self._k_at(times)
        radius = np.hypot(k[0], k[1])
        near = radius < _ORIGIN_FRACTION * self.dk_per_m
        if not near.any():
            return (), ()
        # Group contiguous runs of near-origin samples; each run is one crossing.
        edges = np.flatnonzero(np.diff(near.astype(int)) != 0) + 1
        runs = np.split(np.flatnonzero(near), np.searchsorted(np.flatnonzero(near), edges))
        samples = tuple(int(run[int(np.argmin(radius[run]))]) for run in runs if len(run))
        # Reported on the block's clock, which is what a composing kernel places against.
        lead = self.prephaser_duration_s
        return samples, tuple(float(times[s]) + lead for s in samples)

    def _place_adcs(self) -> tuple[float, ...]:
        """
        Return where each ADC event starts, from the readout block's origin.

        An ADC reserves a dead time **before** sampling as well as after, so consecutive events
        cannot simply abut: laying them end to end by sampling duration alone overlaps them by two
        dead times, and the compiler refuses it -- correctly, since a receiver cannot open two
        windows at once.

        The gap is real and is carried into :meth:`sample_times_s`, so samples either side of a
        segment boundary are **not** uniformly spaced.  Reporting them as if they were would put
        every later sample's k in the wrong place, which is the same class of error as reporting
        the design trajectory.
        """
        raster = float(self.opts.grad_raster_time)
        starts: list[float] = []
        cursor = 0.0
        for adc in self._adcs:
            starts.append(cursor)
            span = float(adc.delay) + float(adc.num_samples) * self.dwell_s + float(adc.dead_time)
            cursor += float(ceil_raster(span, raster))
        return tuple(starts)

    def _arm_sample_times(self) -> np.ndarray:
        """
        Return each ADC sample's time from the start of the **arm**.

        The internal clock.  `_knots` starts at the arm, so every k integration uses this one;
        :meth:`sample_times_s` is the same thing on the block's clock, and the two differ by a
        prephaser.

        Per segment, because the segments are not contiguous -- and offset by the ADC's own
        leading delay, which pypulseq gives it and which is one dead time long.
        """
        times: list[np.ndarray] = []
        for start, adc in zip(self._adc_starts, self._adcs, strict=True):
            count = int(adc.num_samples)
            times.append(start + float(adc.delay) + (np.arange(count) + 0.5) * self.dwell_s)
        return np.concatenate(times)

    def sample_times_s(self) -> np.ndarray:
        """
        Return each ADC sample's time from the start of the readout block :meth:`build` returns.

        **Which includes the prephaser**, and that is the whole reason this is not
        :meth:`_arm_sample_times`.  ``'in'`` and ``'in-out'`` begin at ``k_max``, so their block
        opens with a dephaser and the arm does not start at zero -- reporting arm time here makes
        every sample, the origin crossing and therefore ``TE`` early by that dephaser, in a
        sequence that compiles and whose trajectory is exactly right.  `RadialReadout.
        time_to_center` has always been block-relative and this now matches it.

        With ``build(prephase=False)`` a caller owns the dephaser, and the block then starts at
        the arm: subtract :attr:`prephaser_duration_s`.
        """
        return self._arm_sample_times() + self.prephaser_duration_s

    def _k_at(self, times: np.ndarray) -> np.ndarray:
        """
        Return ``(2, n)`` k at `times`, **integrated from the emitted events' own knots**.

        Not interpolated from the design pass's k.  An arbitrary gradient's samples sit at raster
        *centres* and it is piecewise linear through them, so its integral is piecewise quadratic
        -- linearly interpolating the design knots instead disagrees by about 0.03 `dk`, which is
        invisible in a trajectory plot and is a blur and a shading in an image.

        Exact: cumulative trapezoid over the knots, plus the partial trapezoid into the interval
        each sample time falls in.
        """
        out = np.empty((2, len(times)), dtype=float)
        for axis in (0, 1):
            knot_t, knot_g = self._knots[axis]
            cumulative = np.concatenate([[0.0], np.cumsum(
                0.5 * (knot_g[1:] + knot_g[:-1]) * np.diff(knot_t),
            )])
            index = np.clip(np.searchsorted(knot_t, times, side='right') - 1,
                            0, len(knot_t) - 2)
            span = knot_t[index + 1] - knot_t[index]
            into = np.clip(times - knot_t[index], 0.0, span)
            # g is linear across the interval, so the partial area is a trapezoid of its own.
            slope = (knot_g[index + 1] - knot_g[index]) / np.maximum(span, 1e-30)
            # The integral starts wherever the acquisition starts.  `'in'` and `'in-out'` begin
            # at k_max, delivered by the prephaser, and reporting their k relative to zero would
            # describe a trajectory twice the requested extent.
            out[axis] = (self._waveform[axis, 0] + cumulative[index]
                         + into * (knot_g[index] + 0.5 * slope * into))
        return out

    @property
    def _knots(self) -> tuple[tuple[np.ndarray, np.ndarray], ...]:
        """Return ``((t, g), (t, g))`` for x and y, as the emitted arbitrary gradients hold them."""
        if getattr(self, '_knot_cache', None) is None:
            raster = float(self.opts.grad_raster_time)
            waveform = np.diff(self._waveform, axis=1) / raster
            centres = (np.arange(waveform.shape[1]) + 0.5) * raster
            # `first=0` and `last=0` are the endpoint policy, and they are knots like any other.
            t = np.concatenate([[0.0], centres, [waveform.shape[1] * raster]])
            self._knot_cache = tuple(
                (t, np.concatenate([[0.0], waveform[axis], [0.0]])) for axis in (0, 1)
            )
        return self._knot_cache

    def k_per_m(self, *, angle_rad: float = 0.0) -> np.ndarray:
        """
        Return ``(2, num_samples)`` sample positions in 1/m, **measured off the built events**.

        Integrated from the emitted gradients' own knots, not from the design pass -- whose ``k``
        is discarded and held on no attribute.  A NUFFT is *told* where the samples are, so a
        trajectory reported one raster from where it played is a blur, a shading and a wrong field
        map, and no k-space check would show it.
        """
        k = self._k_at(self._arm_sample_times())
        if not angle_rad:
            return k
        cos, sin = np.cos(angle_rad), np.sin(angle_rad)
        return np.vstack([cos * k[0] - sin * k[1], sin * k[0] + cos * k[1]])

    def echo_sample(self, index: int = 0) -> int:
        """Return the ADC sample at origin crossing `index`."""
        return self.origin_crossing_samples[self._check_index(index)]

    def time_to_echo(self, index: int = 0) -> float:
        """
        Return the time from the readout block's start to origin crossing `index`.

        **This is a fact about the trajectory, not about magnetisation.**  Whether the physical
        echo lands here is the surrounding kernel's to establish; a spin echo placed against the
        wrong crossing is a legal sequence whose echo time is wrong by about one arm.
        """
        return self.origin_crossing_times[self._check_index(index)]

    @property
    def echo_spacing_s(self) -> float | None:
        """
        Seconds between consecutive origin crossings, or ``None`` when there is only one.

        **Derived from the realised trajectory**, never stored: it is fixed by the arm's traversal
        and the variant, and a caller cannot lengthen it without changing `density`, `shots` or
        the scanner.
        """
        times = self.origin_crossing_times
        return None if len(times) < 2 else float(times[1] - times[0])

    @property
    def start_k_per_m(self) -> float:
        """Radius the acquisition starts at: ``0`` for ``'out'`` and ``'out-in'``."""
        return float(np.hypot(*self._waveform[:, 0]))

    @property
    def end_k_per_m(self) -> float:
        """Radius the acquisition ends at.  Non-zero means a rewinder is required."""
        return float(np.hypot(*self._emitted_end))

    @property
    def needs_prephase(self) -> bool:
        """Whether the acquisition starts away from the origin, so a moment is required first."""
        return self.start_k_per_m > _AT_ORIGIN * self.dk_per_m

    @property
    def prephaser_duration_s(self) -> float:
        """
        How long the dephaser at the front of :meth:`build`'s block lasts, or ``0``.

        Public because it is the offset between this module's two clocks, and a caller who passes
        ``prephase=False`` needs it to put the reported times back on the block it assembled.
        The sibling `RadialReadout.prephaser_duration_s` is the same quantity.
        """
        if not self.needs_prephase:
            return 0.0
        if getattr(self, '_prephaser_cache', None) is None:
            self._prephaser_cache = float(self._ramp(self._waveform[:, 0], 0.0).duration)
        return self._prephaser_cache

    @property
    def needs_rewind(self) -> bool:
        """
        Whether the acquisition ends away from the origin, so a moment is required after.

        Thresholded rather than compared with zero: the end point is *integrated* from the
        emitted knots, so a variant that finishes at the origin finishes there to a float, not
        exactly.
        """
        return self.end_k_per_m > _AT_ORIGIN * self.dk_per_m

    @property
    def _emitted_end(self) -> np.ndarray:
        """
        Where k actually finishes, integrated from the **emitted** knots.

        Not ``_waveform[:, -1]``, which is where the design pass meant to finish.  The emitted
        gradient is piecewise linear through raster-*centre* samples with zeros forced at both
        block edges, so its integral differs from the design by a fraction of a step -- and a
        rewinder sized against the design leaves that fraction behind.

        PR23 measured the same thing at 0.108 1/m and called it invisible in every plot and a
        phase ramp across the next excitation's image.  Correcting on the assembled waveform is
        the only place it can be correct, because the residual belongs to the assembly.
        """
        end = float(self._waveform.shape[1] - 1) * float(self.opts.grad_raster_time)
        return self._k_at(np.array([end]))[:, 0]

    # ----------------------------------------------------------------------- assembly
    def build(self, *, angle_rad: float = 0.0, acquire: bool = True,
              phase_deg: float = 0.0, prephase: bool = True, rewind: bool = True) -> LogicBlock:
        """
        Return the readout: an optional prephaser, the arms, the ADCs, an optional rewinder.

        Parameters
        ----------
        angle_rad
            Interleaf rotation.  Applied to the emitted events, so the block **is** the rotated
            acquisition -- the module's reported trajectory and its waveform have one owner.
        acquire
            ``False`` plays the gradients without ADCs, for a dummy shot.
        prephase, rewind
            Whether this module *realises* the moment it requires at each end.  It always
            **states** them -- :attr:`start_k_per_m` and :attr:`end_k_per_m` -- and ``False``
            hands the realisation to a caller who is folding it into something else.  The same
            split ``CartesianLine(prephase=False)`` and ``Excitation.build(rephase=False)`` use.
        """
        out = LogicBlock()
        start = 0.0
        if prephase and self.needs_prephase:
            lead = self._ramp(self._waveform[:, 0], angle_rad)
            out.add(0.0, lead)
            # From the property, so the reported clock and the emitted block cannot drift: a
            # rotation changes which axes carry the dephaser, never how long it takes.
            start = self.prephaser_duration_s
            # A boundary, for the same reason the rewinder gets one: superposed into the arm's
            # block, a trapezoid's raster-edge knots and an arbitrary gradient's raster-centre
            # samples are summed on lattices that disagree, and the emitted area stops matching
            # the tree's -- silently, and by about 0.013 `dk`.
            out.add(start, barrier())

        for axis_index, axis in enumerate('xy'):
            event = self._rotated(axis_index, axis, angle_rad)
            if event is not None:
                out.add(start, event)

        if acquire:
            phase = float(np.deg2rad(phase_deg))
            for offset, adc in zip(self._adc_starts, self._adcs, strict=True):
                out.add(start + offset, derive(adc, phase_offset=phase) if phase else adc)

        if rewind and self.needs_rewind:
            # After the acquisition *ends*, not after the gradient does.  An ADC reserves its
            # window plus a trailing dead time, which can outlast the arm, and a rewinder placed
            # at the arm's end would land inside that reservation.
            #
            # `sc.barrier()` then forces a block boundary, and it is load-bearing: superposed in
            # one block, the arm and the rewinder are summed on lattices that disagree -- an
            # arbitrary gradient's samples sit at raster centres and a trapezoid's knots at raster
            # edges -- and the emitted area silently stops matching the tree's.  Measured at 0.058
            # 1/m on a 64-matrix arm, which is 0.013 `dk`: invisible in any plot, and a phase ramp
            # across the next excitation's image.
            tail = start + self.acquisition_end_s
            out.add(tail, barrier())
            out.add(tail, self._ramp(-self._emitted_end, angle_rad))
        return out

    @property
    def acquisition_end_s(self) -> float:
        """
        Seconds from the **arm's** start until nothing is reserved any more.

        The ADC's own reservation, read off the event: a leading delay, the sampling window, and
        a trailing dead time.  It can outlast the gradient, and a rewinder placed at the arm's
        end would then land inside it.

        On the **arm's** clock, deliberately -- it is compared against :attr:`duration_s`, which
        is the arm's gradient, and "the acquisition finishes inside its gradient" is a statement
        about those two and nothing else.  :meth:`sample_times_s` is on the block's clock because
        a composing kernel places against the block; see :attr:`prephaser_duration_s`.
        """
        last = self._adcs[-1]
        adc = (self._adc_starts[-1] + float(last.delay)
               + float(last.num_samples) * self.dwell_s + float(last.dead_time))
        # Onto the gradient raster: what follows is a gradient, and a block boundary off the
        # raster is refused where it is computed rather than where it is used.
        return float(ceil_raster(max(self.duration_s, adc), float(self.opts.grad_raster_time)))

    def _rotated(self, index: int, axis: str, angle_rad: float):
        """Return this axis's gradient for an acquisition rotated by `angle_rad`."""
        raster = float(self.opts.grad_raster_time)
        waveform = np.diff(self._waveform, axis=1) / raster
        if angle_rad:
            cos, sin = np.cos(angle_rad), np.sin(angle_rad)
            waveform = np.vstack([cos * waveform[0] - sin * waveform[1],
                                  sin * waveform[0] + cos * waveform[1]])
        if not np.any(np.abs(waveform[index]) > 0.0):
            return None
        return pp.make_arbitrary_grad(channel=axis, waveform=waveform[index], system=self.opts,
                                      delay=0.0, first=0.0, last=0.0)

    def _ramp(self, k_target: np.ndarray, angle_rad: float) -> LogicBlock:
        """
        Return a two-axis trapezoid delivering `k_target`, rotated.

        Private on purpose.  An oblique dephaser split across two axes by direction cosines is
        plausibly shareable, but nothing else in the package needs it today, and one consumer does
        not justify promotion to ``_support``.
        """
        cos, sin = np.cos(angle_rad), np.sin(angle_rad)
        target = np.array([cos * k_target[0] - sin * k_target[1],
                           sin * k_target[0] + cos * k_target[1]])
        total = float(np.hypot(*target))
        probe = pp.make_trapezoid(channel='x', area=total, system=self.opts)
        duration = float(pp.calc_duration(probe))
        block = LogicBlock()
        for index, axis in enumerate('xy'):
            if abs(target[index]) > 0.0:
                block.add(0.0, pp.make_trapezoid(channel=axis, area=float(target[index]),
                                                 duration=duration, system=self.opts))
        return block

    # -------------------------------------------------------------------- the refusals
    def _check_variant(self, variant: str) -> str:
        """Return `variant` having checked it names one of the four."""
        if variant not in _VARIANTS:
            listed = ', '.join(repr(name) for name in _VARIANTS)
            msg = format_error(
                f'variant must be one of {listed}, got {variant!r}.', {'variant': variant},
                ['all four are one arm traversed differently; the arm is at rest at both ends'],
            )
            raise ConfigurationError(msg)
        return variant

    def _check_echoes(self, echoes: int) -> int:
        """Return `echoes`, refusing the multi-echo case this increment does not implement."""
        value = require_count(echoes, 'echoes', low=1)
        if value != 1:
            msg = format_error(
                f'echoes={value} is not implemented yet; only echoes=1 is.',
                {'echoes': value, 'variant': self.variant},
                ['multi-echo is the last step of the recorded implementation order: it needs '
                 'the transition rules -- a fly-back for one-arm variants, continuous traversal '
                 'for two-arm ones -- and nothing exercises them yet',
                 "variant='out-in' already gives two origin crossings at echoes=1"],
            )
            raise ConfigurationError(msg)
        return value

    def _check_density(self, density: tuple[float, ...]) -> tuple[float, ...]:
        """Return `density` as sampled field-of-view multipliers, all positive."""
        values = tuple(float(v) for v in np.atleast_1d(density))
        if not values or any(v <= 0.0 for v in values):
            msg = format_error(
                'density must be one or more positive field-of-view multipliers.',
                {'density': density},
                ['(1.0,) is uniform', '(1.0, 0.5) halves the field of view at the edge',
                 'these are SAMPLED VALUES, not polynomial coefficients'],
            )
            raise ConfigurationError(msg)
        return values

    def _check_index(self, index: int) -> int:
        """Return `index` having checked this variant has that many origin crossings."""
        count = len(self.origin_crossing_samples)
        if not 0 <= index < count:
            msg = format_error(
                f'origin crossing {index} does not exist; {self.variant!r} has {count}.',
                {'index': index, 'variant': self.variant, 'crossings': count},
                ["'out-in' has two; the other three have one"],
            )
            raise ConfigurationError(msg)
        return index
