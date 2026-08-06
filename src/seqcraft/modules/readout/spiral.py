"""
Variable-density spiral readout.

The trajectory is the classic Glover dual-rate design: a slew-limited spiral near the origin,
becoming amplitude-limited further out, with a density that can taper from fully sampled at the
centre to undersampled at the edge.  It is generated analytically and then integrated to a
gradient waveform, which is what makes the module honest about its own limits -- the waveform is
checked against the amplifier after generation, not assumed to fit.

Why a spiral for diffusion
--------------------------
A spiral collects a whole 2D plane per shot from a single echo, like EPI, but its k-space
sampling starts at the origin.  That puts the highest-SNR samples at k=0 where the diffusion
attenuation is measured, and it makes the readout immune to the phase-encode-direction
distortion that makes EPI-based DWI hard to register.  The cost is that off-resonance blurs
rather than shifts, so the readout must be kept short -- which is why `readout_duration_us` is a
first-class parameter here rather than something derived.

Examples
--------
>>> import seqcraft as sc
>>> system = sc.System.preset('generic_3t').derate('spiral', grad=0.9, slew=0.7)
>>> ro = SpiralVDS(system, fov_mm=240, matrix=96, n_interleaves=8, regime='spiral')
>>> ro.n_interleaves
8
>>> round(ro.k_max_per_m, 1)                         # matrix / (2 * FOV)
200.0
>>> ro.readout_duration_us > 0
True
>>> ro.build(interleaf=0)
LogicBlock(spiral, 5 nodes, 6.12 ms)
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import pypulseq as pp

from ...core import events as ev
from ...core.errors import ConfigurationError, format_error
from ...core.logic import LogicBlock
from ...core.module import Module
from ...core.units import convert
from ...core.validate import Range, require_in_range, require_int_in, require_positive

if TYPE_CHECKING:
    from types import SimpleNamespace

    from ...core.system import System

__all__ = ['Spiral', 'SpiralVDS', 'vds_trajectory']

_FOV_RANGE = Range(0.5, 2000.0, 'mm', ((1e3, 'm'), (10.0, 'cm')))


class Spiral(NamedTuple):
    """
    A generated spiral: its k-space trajectory and the gradient waveform that produces it.

    Attributes
    ----------
    kx, ky
        k-space coordinates in 1/m, starting at the origin.
    gx, gy
        Gradient amplitudes in Hz/m on the gradient raster.  Both start and end at exactly zero,
        so the waveform can be followed by any block without a discontinuity.
    """

    kx: np.ndarray
    ky: np.ndarray
    gx: np.ndarray
    gy: np.ndarray


def vds_trajectory(
    *,
    fov_m: float,
    k_max_per_m: float,
    n_interleaves: int,
    max_grad: float,
    max_slew: float,
    raster_s: float,
    alpha: float = 1.0,
    density: float = 1.0,
    max_samples: int = 200_000,
) -> Spiral:
    """
    Generate a variable-density spiral, traversed as fast as the amplifier allows.

    Parameters
    ----------
    fov_m
        Field of view, metres.
    k_max_per_m
        Outer k-space radius, 1/m.  ``matrix / (2 * FOV)``.
    n_interleaves
        Number of rotated shots covering k-space together.
    max_grad, max_slew
        Amplifier limits in Hz/m and Hz/m/s.
    raster_s
        Gradient raster, seconds.
    alpha
        Variable-density exponent.  ``1`` is a uniform-density (Archimedean) spiral; larger values
        sample the centre more densely and the edge more sparsely.
    density
        Turns relative to Nyquist.  ``1`` fully samples k-space for `n_interleaves`; ``0.6`` winds
        60 % of the turns, so the readout is 40 % shorter and the outer k-space is undersampled by
        1/0.6.  That trade is what makes a **single-shot** spiral possible at a useful resolution:
        aliasing from an undersampled outer ring is far more benign than the T2 decay and
        off-resonance blur of a readout twice as long.
    max_samples
        Safety bound on the integration, so a bad parameter set fails rather than hangs.

    Returns
    -------
    Spiral

    Notes
    -----
    Writes the spiral as ``k(theta) = r(theta) * exp(i*theta)`` with ``r = c * theta**(1/alpha)``
    and walks ``theta`` forward one raster at a time.  Each step is the largest that satisfies the
    constraints the hardware actually applies, which are on the **sampled** waveform:

    .. code-block:: text

        g[n] = (k(theta[n]) - k(theta[n-1])) / raster     exact, so cumsum(g) * raster == k

        |g[n]|          <= G_max
        |g[n] - g[n-1]| <= S_max * raster

    capped additionally at the speed the curvature ahead can sustain,

    .. code-block:: text

        |dk/dtheta|    = sqrt(r'^2 + r^2)                 tangential scale
        |d2k/dtheta^2| = sqrt((r'' - r)^2 + (2 r')^2)     curvature scale

        theta' <= min( G_max / |dk/dtheta|,  sqrt(S_max / |d2k/dtheta^2|) )

    Neither half is redundant.  Without the discrete constraints, a continuous-time model of the
    slew has to be paid for twice -- once in its own error near the origin, where one raster is not
    a small angle, and again in the global derating needed to bring the realised waveform back
    inside the limits.  Without the cap, a greedy forward pass accelerates to a speed the next turn
    cannot hold, and since shedding speed also costs slew there is then no legal move at all: the
    integration stalls a few radians out.  With both, the result runs at 100 % of peak slew and
    **99.8 % of it on average**, which is what "slew-limited spiral" ought to mean.

    Near the origin the curvature cap binds, because the trajectory is turning hard while its
    radius is still small; further out the amplitude limit binds, because ``|dk/dtheta|`` grows
    with the radius.  That crossover is the classic dual-rate spiral.

    The waveform is then ramped down to zero at the slew limit, so it can be followed directly by
    a rewinder with no gradient discontinuity -- pypulseq rejects a block whose gradient starts
    away from the previous block's end value.  That ramp carries ``|k|`` about 12 % past ``k_max``;
    those samples are ordinary k-space points a little beyond the nominal resolution, and density
    compensation takes them in its stride.  Scaling the finished waveform down to bring ``max|k|``
    back to exactly ``k_max`` is the tempting alternative and is a trap: it scales the **turn
    spacing** with it, so a spiral asked for Nyquist sampling silently delivers 16 % coarser.

    Examples
    --------
    >>> spiral = vds_trajectory(fov_m=0.24, k_max_per_m=133.3, n_interleaves=8,
    ...                         max_grad=1.45e6, max_slew=4.15e9, raster_s=1e-5)
    >>> k_max = float(np.hypot(spiral.kx, spiral.ky).max())
    >>> 1.0 <= k_max / 133.3 <= 1.2   # reaches k_max; the ramp-down goes a bit past
    True
    >>> float(spiral.gx[0]), float(spiral.gx[-1])                  # starts and ends at rest
    (0.0, 0.0)
    >>> bool(np.hypot(spiral.gx, spiral.gy).max() <= 1.45e6 * 1.001)
    True
    >>> slew = np.hypot(np.diff(spiral.gx), np.diff(spiral.gy)) / 1e-5
    >>> bool(slew.max() <= 4.15e9 * 1.001)
    True
    """
    if k_max_per_m <= 0 or n_interleaves < 1 or density <= 0 or alpha <= 0:
        msg = format_error(
            'spiral parameters must be positive.',
            {
                'k_max_per_m': k_max_per_m,
                'n_interleaves': n_interleaves,
                'density': density,
                'alpha': alpha,
            },
        )
        raise ConfigurationError(msg)

    # One interleaf must reach k_max while adjacent turns stay 1/FOV apart, with the winding shared
    # between n_interleaves shots.  `density` then scales the turn count directly: 0.6 winds 60 % of
    # them, which is 40 % less readout and 1/0.6 undersampling at the edge.
    #
    # It multiplies rather than divides.  Dividing -- which is what this did at first -- makes
    # `density=0.6` wind *more* turns and run 66 % longer, the exact opposite of what the name says,
    # and it silently turns the one lever a single-shot spiral has into a penalty.
    theta_max = 2.0 * math.pi * k_max_per_m * fov_m * density / n_interleaves
    if theta_max <= 0:
        msg = format_error('spiral has no extent; check fov_mm and matrix.', {})
        raise ConfigurationError(msg)

    # The rate-limited forward integration is a very good approximation, not an exact solution --
    # a true time-optimal traversal needs a backward pass as well.  Rather than trust it, measure
    # the finished waveform and tighten the budget until it genuinely fits, so this function's
    # contract is "within the limits given" rather than "ought to be".
    budget = 1.0
    for _attempt in range(24):
        gx, gy = _integrate_spiral(
            theta_max=theta_max,
            k_max_per_m=k_max_per_m,
            power=1.0 / alpha,
            max_grad=max_grad * budget,
            max_slew=max_slew * budget,
            raster_s=raster_s,
            max_samples=max_samples,
        )
        gx, gy = _ramp_down(gx, gy, max_slew=max_slew, raster_s=raster_s)
        peak_grad = float(np.hypot(gx, gy).max())
        peak_slew = (
            float(np.hypot(np.diff(gx), np.diff(gy)).max() / raster_s) if len(gx) > 1 else 0.0
        )
        if peak_grad <= max_grad * 1.0005 and peak_slew <= max_slew * 1.0005:
            return Spiral(kx=np.cumsum(gx) * raster_s, ky=np.cumsum(gy) * raster_s, gx=gx, gy=gy)
        # Amplitude scales with the budget, slew roughly with its square.
        over = max(peak_grad / max_grad, math.sqrt(peak_slew / max_slew))
        budget = budget / max(over, 1.005) * 0.999

    msg = format_error(  # pragma: no cover - needs contradictory limits
        'could not find a spiral within the given gradient limits.',
        {'k_max_per_m': k_max_per_m, 'n_interleaves': n_interleaves, 'alpha': alpha},
        ['increase n_interleaves', 'reduce matrix', 'increase density'],
    )
    raise ConfigurationError(msg)


def _integrate_spiral(
    *,
    theta_max: float,
    k_max_per_m: float,
    power: float,
    max_grad: float,
    max_slew: float,
    raster_s: float,
    max_samples: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Walk theta forward, taking the largest step the **discrete** limits allow at each raster.

    The waveform returned is the exact finite difference of the analytic trajectory,
    ``g[n] = (k(theta[n]) - k(theta[n-1])) / raster``, so ``cumsum(g) * raster`` reproduces
    ``k(theta[n])`` to the last bit.  Both quantities that matter are then exact by construction and
    need no correction afterwards: the spiral ends at ``k_max``, and adjacent turns are
    ``1 / (FOV * density)`` apart.

    That is worth more than it sounds.  Differentiating the trajectory analytically instead --
    ``g = k'(theta) * theta'`` -- leaves the cumulative sum of the *sampled* gradient short of the
    analytic radius, because near the origin the spiral turns by an eighth of a radian per raster
    and the polygon cuts the corners.  Correcting that by scaling the finished waveform to land on
    ``k_max`` scales the **turn spacing** with it, so a spiral asked for Nyquist sampling delivers
    16 % coarser, and the aliasing that follows looks like undersampling nobody asked for.

    The step is chosen against the discrete constraints themselves,

    .. code-block:: text

        |g[n]|             <= G_max
        |g[n] - g[n-1]|    <= S_max * raster

    rather than against a continuous-time model of them.  A continuous model has to be paid for
    twice: once in the approximation near the origin, where a raster is not a small angle, and again
    in the global derating needed to bring the realised waveform back inside the limits.  Solving
    the constraint that the hardware actually applies costs a bisection per raster step and leaves
    the amplifier running at its limit for the whole readout -- which for a slew-limited spiral is
    the difference between a 20 ms readout and a 29 ms one.

    Feasibility is an *interval* in the step size, not a half-line: too large a step outruns the
    slew limit, and too small a one cannot decelerate from the previous gradient inside it either.
    So the bracket starts from the previous step, which sits near the middle of that interval, and
    expands upward.
    """
    c = k_max_per_m / theta_max**power
    slew_step = max_slew * raster_s
    theta = 0.0
    position_x = position_y = 0.0
    previous_gx = previous_gy = 0.0
    gx: list[float] = [0.0]
    gy: list[float] = [0.0]
    # Seed the bracket below any real step, so the first iteration expands up to the answer.
    delta = theta_max * 1e-9

    def gradient_after(step: float) -> tuple[float, float, float]:
        """Return the gradient produced by advancing theta, and the angle actually reached."""
        reached = min(theta + step, theta_max)
        radius = c * reached**power
        return (
            (radius * math.cos(reached) - position_x) / raster_s,
            (radius * math.sin(reached) - position_y) / raster_s,
            reached,
        )

    def fits(step: float) -> bool:
        if step <= 0.0:
            return False
        candidate_x, candidate_y, _ = gradient_after(step)
        return (
            math.hypot(candidate_x, candidate_y) <= max_grad
            and math.hypot(candidate_x - previous_gx, candidate_y - previous_gy) <= slew_step
        )

    for _ in range(max_samples):
        if theta >= theta_max:
            break
        # The sustainable rate. Taking the largest step the slew limit allows and nothing more is
        # not enough: that step can leave theta' above the speed the curvature ahead can hold, and
        # since shedding speed also costs slew the trajectory then has no legal move at all. The
        # spiral stalls a few radians in. Capping at the steady-state ceiling is what a backward
        # pass would otherwise be needed for.
        guard = max(theta, theta_max * 1e-6)
        radius = c * guard**power
        first = c * power * guard ** (power - 1.0)
        second = c * power * (power - 1.0) * guard ** (power - 2.0)
        tangential = max(math.hypot(first, radius), 1e-12)              # |dk/dtheta|
        curvature = max(math.hypot(second - radius, 2.0 * first), 1e-12)  # |d2k/dtheta^2|
        ceiling = min(max_grad / tangential, math.sqrt(max_slew / curvature)) * raster_s

        if fits(ceiling):
            delta = ceiling            # cruising: the curvature binds, not the slew
        else:
            # The slew binds, so bisect down from the ceiling -- but the bisection needs a step that
            # fits to start from, and the window of steps that do is narrow. Its width is about one
            # step's own turn angle, a few percent, because a step much smaller than the last one
            # has to shed most of the current gradient in one raster and that costs slew too. So the
            # ladder is fine near 1.0: cruising, the answer sits a fraction of a percent below the
            # previous step, and a coarse ladder jumps clean over the window and reports a spiral
            # that cannot continue.
            base = min(delta, ceiling)
            low = 0.0
            for factor in (1.0, 0.9998, 0.9995, 0.999, 0.998, 0.995, 0.99, 0.98, 0.95, 0.9, 0.7, 0.4):
                if fits(base * factor):
                    low = base * factor
                    break
            if low == 0.0:
                if theta_max - theta <= 2.0 * delta:
                    # The last hop onto theta_max exactly is shorter than a full step, so the
                    # gradient it needs is a step change the amplifier cannot make. Stop here and
                    # let the ramp-down close the gap: what is given up is a fraction of one turn
                    # out of hundreds, which moves the edge of k-space by under one part in 10^3.
                    break
                msg = format_error(  # pragma: no cover - would need contradictory limits
                    'the spiral cannot continue within the slew limit.',
                    {'theta': theta, 'theta_max': theta_max},
                    ['increase n_interleaves', 'reduce matrix', 'reduce density'],
                )
                raise ConfigurationError(msg)
            high = ceiling
            for _ in range(40):
                middle = 0.5 * (low + high)
                if fits(middle):
                    low = middle
                else:
                    high = middle
            delta = low

        new_gx, new_gy, theta = gradient_after(delta)
        position_x += new_gx * raster_s
        position_y += new_gy * raster_s
        previous_gx, previous_gy = new_gx, new_gy
        gx.append(new_gx)
        gy.append(new_gy)
    else:  # pragma: no cover - only on absurd parameters
        msg = format_error(
            f'spiral did not reach k_max within {max_samples} raster steps.',
            {'theta_max': theta_max, 'reached': theta},
            ['reduce matrix, increase n_interleaves, or increase density'],
        )
        raise ConfigurationError(msg)

    return np.asarray(gx), np.asarray(gy)


def _ramp_down(
    gx: np.ndarray,
    gy: np.ndarray,
    *,
    max_slew: float,
    raster_s: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Append a slew-limited ramp to zero, then pad one zero sample at each end.

    A spiral naturally finishes at full amplitude.  Leaving it there means the next block starts
    away from the previous block's end value, which pypulseq rejects -- and physically it is a step
    change in gradient, which no amplifier can produce.  Ramping the vector down along its own
    direction is the cheapest way to land at rest.

    The **padding** is subtler and matters just as much.  ``make_arbitrary_grad`` treats waveform
    values as samples at raster *centres*, so ``first`` and ``last`` are the extrapolated values at
    the raster *edges* -- ``waveform[0] - (waveform[1] - waveform[0]) / 2`` and its mirror.  Ending
    the waveform with a single zero sample does not make that extrapolation zero: the ramp's slope
    carries it past.  Declaring ``last=0`` anyway makes pypulseq's shape round-trip disagree with
    the recorded value by a few tenths of a mT/m, which it reports as *"Last restored point differs
    too much from the recorded last; skipping shape restoration"* -- and the restored waveform is
    then not the one that was designed.  A zero sample at each end makes the extrapolation genuinely
    zero, so the declaration is true rather than asserted.  It costs one raster period per end.
    """
    end = math.hypot(float(gx[-1]), float(gy[-1]))
    if end > 0.0:
        n = max(1, int(math.ceil(end / (max_slew * raster_s))))
        taper = np.linspace(1.0, 0.0, n + 1)[1:]
        gx = np.concatenate([gx, float(gx[-1]) * taper])
        gy = np.concatenate([gy, float(gy[-1]) * taper])
    pad = np.zeros(1)
    return (
        np.concatenate([pad, gx, pad]),
        np.concatenate([pad, gy, pad]),
    )


class SpiralVDS(Module):
    """
    Variable-density spiral readout with its own rewinder.

    Parameters
    ----------
    system
        The scanner.
    fov_mm
        In-plane field of view, millimetres.
    matrix
        In-plane matrix size; sets ``k_max = matrix / (2 * FOV)``.
    n_interleaves
        Number of rotated shots.  More interleaves means a shorter readout each, so less
        off-resonance blur, at the cost of more shots.
    density
        Sampling density relative to Nyquist for `n_interleaves`.  ``1`` is fully sampled;
        ``0.6`` undersamples the outer k-space, which is normal for diffusion where SNR is the
        binding constraint.
    alpha
        Variable-density exponent; ``1`` is uniform density.
    adc_dwell_ns
        ADC dwell time, nanoseconds.  Must be a multiple of the ADC raster (100 ns).
    rewind
        Append a rewinder returning both axes to k=0, so the next TR starts from the origin.
    regime
        Which limit regime to design the waveform against.  A spiral is usually run derated,
        because the continuous slewing is the worst case for peripheral nerve stimulation and
        for acoustic noise.

    Properties
    ----------
    k_max_per_m
        Outer k-space radius, 1/m.
    readout_duration_us
        Duration of the spiral gradient itself.
    time_to_echo
        Seconds from the start of the built block to k=0 -- **which for a spiral is the first
        sample**, not the middle.  This is what makes TE placement differ from a Cartesian
        readout, and getting it wrong shifts the diffusion weighting rather than the image.
    duration
        Seconds occupied by the whole built block.
    resolution_mm
        ``FOV / matrix``.

    Build arguments
    ---------------
    interleaf : int, default 0
        Which rotated shot to play.  Rotation is applied to the waveform, so all interleaves
        share one design.

    `rewind` is deliberately **not** a build argument: dropping the rewinder changes the block's
    duration, which would silently invalidate :attr:`duration`.  Set it on the constructor.

    Notes
    -----
    **k=0 is the first sample.**  A Cartesian readout puts the echo in the middle of its ADC
    window; a spiral starts at the origin, so ``time_to_echo`` is the ADC delay alone.  Placing a
    spiral at ``te - time_to_echo`` therefore puts the *start* of sampling at TE, which is what
    the diffusion b-value is referenced to.

    The rewinder is part of this module rather than something the caller adds, because its area
    is exactly the negative of the spiral's own and nothing else can know that.

    Examples
    --------
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t').derate('spiral', grad=0.9, slew=0.7)
    >>> ro = SpiralVDS(system, fov_mm=240, matrix=96, n_interleaves=8, density=0.6,
    ...                regime='spiral')
    >>> round(ro.resolution_mm, 2)
    2.5
    >>> round(ro.time_to_echo * 1e6) == round(ro.adc.delay * 1e6)
    True
    >>> len(ro.build(interleaf=3).nodes)             # gx, gy, adc, and the rewinder pair
    5
    """

    def __init__(
        self,
        system: System,
        *,
        fov_mm: float,
        matrix: int,
        n_interleaves: int = 1,
        density: float = 1.0,
        alpha: float = 1.0,
        adc_dwell_ns: float = 2000.0,
        rewind: bool = True,
        axes: tuple[str, str] = ('x', 'y'),
        regime: str = 'default',
    ) -> None:
        super().__init__(system, regime=regime)
        self.fov_mm = float(fov_mm)
        self.matrix = int(matrix)
        self.n_interleaves = int(n_interleaves)
        self.density = float(density)
        self.alpha = float(alpha)
        self.adc_dwell_ns = float(adc_dwell_ns)
        self.rewind = bool(rewind)
        self.axes = tuple(axes)

        require_positive(self, 'fov_mm', 'density', 'alpha', 'adc_dwell_ns')
        require_in_range(self, 'fov_mm', _FOV_RANGE.lo, _FOV_RANGE.hi, unit='mm')
        require_int_in(self, 'matrix', lo=8, hi=2048)
        require_int_in(self, 'n_interleaves', lo=1, hi=1024)

        raster = self.system.grad_raster.dt
        spiral = vds_trajectory(
            fov_m=self.fov_mm / 1e3,
            k_max_per_m=self.k_max_per_m,
            n_interleaves=self.n_interleaves,
            max_grad=float(self.opts.max_grad),
            max_slew=float(self.opts.max_slew),
            raster_s=raster,
            alpha=self.alpha,
            density=self.density,
        )
        self._kx, self._ky = spiral.kx, spiral.ky
        # The gradient is the exact finite difference of the analytic trajectory, so the k-space
        # positions the reconstruction is told about are the ones the waveform really produces.
        gx_wave, gy_wave = spiral.gx, spiral.gy
        self._check_waveform(gx_wave, gy_wave, raster)

        # first/last stated rather than extrapolated: the trajectory generator ramps the waveform
        # to exactly zero at both ends, and saying so is what lets the block before and after carry
        # no gradient at all without tripping pypulseq's continuity check.
        self.gx = pp.make_arbitrary_grad(
            channel=self.axes[0], waveform=gx_wave, first=0.0, last=0.0,
            system=self.opts, delay=0.0,
        )
        self.gy = pp.make_arbitrary_grad(
            channel=self.axes[1], waveform=gy_wave, first=0.0, last=0.0,
            system=self.opts, delay=0.0,
        )

        dwell = self.adc_dwell_ns / 1e9
        divisor = self.system.adc_samples_divisor
        n_samples = int(len(gx_wave) * raster / dwell) // divisor * divisor
        if n_samples < divisor:
            msg = format_error(
                'the spiral is shorter than one ADC sample group.',
                {
                    'spiral_us': convert(len(gx_wave) * raster, 's', 'us'),
                    'dwell_ns': self.adc_dwell_ns,
                },
                ['reduce adc_dwell_ns, or increase matrix'],
            )
            raise ConfigurationError(msg)
        self.adc = pp.make_adc(
            num_samples=n_samples,
            dwell=dwell,
            delay=0.0,
            system=self.opts,
        )

        self.rewinder = self._design_rewinder()

    def _check_waveform(self, gx: np.ndarray, gy: np.ndarray, raster: float) -> None:
        """
        Reject the design if the generated waveform exceeds the regime it was designed for.

        The trajectory generator respects the limits it was given by construction, so a
        violation here means the finite-difference differentiation at the very start of the
        spiral overshot -- worth catching, because a spiral that clips at the origin distorts
        exactly the samples the diffusion measurement depends on.
        """
        peak_grad = float(max(np.max(np.abs(gx)), np.max(np.abs(gy))))
        slew = float(
            max(
                np.max(np.abs(np.diff(gx) / raster)) if len(gx) > 1 else 0.0,
                np.max(np.abs(np.diff(gy) / raster)) if len(gy) > 1 else 0.0,
            )
        )
        as_mT_m = self.system.convert
        if peak_grad > float(self.opts.max_grad) * 1.001:
            msg = format_error(
                f'the generated spiral reaches '
                f'{as_mT_m(peak_grad, "Hz/m", "mT/m"):.1f} mT/m, above the '
                f'{as_mT_m(float(self.opts.max_grad), "Hz/m", "mT/m"):.1f} mT/m limit of regime '
                f'{self.regime!r}.',
                {'matrix': self.matrix, 'fov_mm': self.fov_mm, 'n_interleaves': self.n_interleaves},
                ['increase n_interleaves', 'reduce matrix', 'design against a less derated regime'],
            )
            raise ConfigurationError(msg)
        if slew > float(self.opts.max_slew) * 1.05:
            msg = format_error(
                f'the generated spiral reaches '
                f'{as_mT_m(slew, "Hz/m/s", "T/m/s"):.0f} T/m/s, above the '
                f'{as_mT_m(float(self.opts.max_slew), "Hz/m/s", "T/m/s"):.0f} T/m/s limit of '
                f'regime {self.regime!r}.',
                {'matrix': self.matrix, 'n_interleaves': self.n_interleaves},
                ['increase n_interleaves', 'design against a less derated regime'],
            )
            raise ConfigurationError(msg)

    def _design_rewinder(self) -> tuple[SimpleNamespace, SimpleNamespace] | None:
        """
        Build the pair of trapezoids returning both axes from the spiral's end to k=0.

        Designed for the **magnitude** of the end k-vector on each axis, not for its unrotated
        components.  Rotating an interleaf moves area between x and y, so an axis that needed
        ``kx_end`` unrotated can need up to ``hypot(kx_end, ky_end)`` at some other angle.
        Designing for the components and scaling up would exceed the slew limit on exactly those
        interleaves -- which is the sort of bug that shows up as a handful of corrupted shots
        rather than as an obvious failure.
        """
        if not self.rewind:
            return None
        worst = math.hypot(float(self._kx[-1]), float(self._ky[-1]))
        if worst == 0.0:
            return None
        reference = pp.make_trapezoid(channel=self.axes[0], area=worst, system=self.opts)
        duration = self.system.grad_raster.ceil(float(pp.calc_duration(reference)))
        return tuple(  # type: ignore[return-value]
            pp.make_trapezoid(channel=axis, area=worst, duration=duration, system=self.opts)
            for axis in self.axes
        )

    # -------------------------------------------------------------------------- properties
    @property
    def k_max_per_m(self) -> float:
        """Outer k-space radius, ``matrix / (2 * FOV)``, in 1/m."""
        return self.matrix / (2.0 * convert(self.fov_mm, 'mm', 'm'))

    @property
    def resolution_mm(self) -> float:
        """In-plane resolution, ``FOV / matrix``, in millimetres."""
        return self.fov_mm / self.matrix

    @property
    def n_samples(self) -> int:
        """Number of ADC samples per interleaf."""
        return int(self.adc.num_samples)

    @property
    def readout_duration_us(self) -> float:
        """Duration of the spiral gradient itself, microseconds."""
        return convert(self.system.grad_raster.at(len(self._kx)), 's', 'us')

    @property
    def time_to_echo(self) -> float:
        """
        Seconds from the start of the built block to k=0.

        For a spiral that is the ADC delay alone, because sampling **starts** at the origin.
        A Cartesian readout puts k=0 in the middle of its window; conflating the two shifts the
        diffusion weighting.
        """
        return float(self.adc.delay)

    @property
    def adc_duration(self) -> float:
        """Length of the ADC window, seconds.  An ADC event has no ``duration`` attribute."""
        return float(self.adc.num_samples) * float(self.adc.dwell)

    @property
    def sampling_block_duration(self) -> float:
        """
        Seconds the spiral's own block must occupy.

        The longer of the gradient waveform and the ADC window plus its **trailing** dead time.
        The second is easy to overlook and is usually the binding one: the ADC runs almost to the
        end of the gradient, and pulseq requires ``adc.delay + samples*dwell + adc_dead_time`` to
        fit inside the block.  Starting the rewinder before that is not merely illegal -- the
        compiler would have nowhere to put a block boundary and would sum the rewinder into the
        spiral, producing a step discontinuity in the middle of the readout.
        """
        gradient = float(pp.calc_duration(self.gx, self.gy))
        sampling = float(self.adc.delay) + self.adc_duration + float(self.opts.adc_dead_time)
        return self.system.block_raster.ceil(max(gradient, sampling))

    @property
    def duration(self) -> float:
        """Seconds occupied by the whole built block, rewinder included."""
        return self.sampling_block_duration + self.rewind_duration

    @property
    def rewind_duration(self) -> float:
        """Seconds occupied by the rewinder, or zero when there is none."""
        if self.rewinder is None:
            return 0.0
        return self.system.block_raster.ceil(float(pp.calc_duration(*self.rewinder)))

    def trajectory(self, interleaf: int = 0) -> tuple[np.ndarray, np.ndarray]:
        """
        Return the k-space trajectory of one interleaf, in 1/m.

        Examples
        --------
        >>> import numpy as np
        >>> import seqcraft as sc
        >>> ro = SpiralVDS(sc.System.preset('generic_3t'), fov_mm=240, matrix=64,
        ...                n_interleaves=4)
        >>> kx0, ky0 = ro.trajectory(0)
        >>> kx1, ky1 = ro.trajectory(1)
        >>> round(float(np.arctan2(ky1[-1], kx1[-1]) - np.arctan2(ky0[-1], kx0[-1])), 4)
        1.5708
        """
        angle = 2.0 * math.pi * interleaf / self.n_interleaves
        cos, sin = math.cos(angle), math.sin(angle)
        return self._kx * cos - self._ky * sin, self._kx * sin + self._ky * cos

    # ------------------------------------------------------------------------------- build
    def build(self, *, interleaf: int = 0) -> LogicBlock:
        """
        Return the spiral gradients, the ADC, and the rewinder if this readout owns one.

        Parameters
        ----------
        interleaf
            Which rotated shot.  The rotation is applied to the stored waveform, so all
            interleaves share one design.

        Examples
        --------
        >>> import seqcraft as sc
        >>> ro = SpiralVDS(sc.System.preset('generic_3t'), fov_mm=240, matrix=64,
        ...                n_interleaves=4, rewind=False)
        >>> ro.build()                              # gx, gy, adc -- no rewinder
        LogicBlock(spiral, 3 nodes, 5.08 ms)
        """
        angle = 2.0 * math.pi * (interleaf % self.n_interleaves) / self.n_interleaves
        cos, sin = math.cos(angle), math.sin(angle)
        wx = np.asarray(self.gx.waveform)
        wy = np.asarray(self.gy.waveform)
        # Rotating a waveform that starts and ends at zero leaves it starting and ending at zero,
        # so first/last stay exactly zero rather than being recomputed from the rotated samples.
        gx = ev.derive(self.gx, waveform=wx * cos - wy * sin, first=0.0, last=0.0)
        gy = ev.derive(self.gy, waveform=wx * sin + wy * cos, first=0.0, last=0.0)

        out = LogicBlock('spiral')
        out.add(0.0, gx, gy, self.adc)

        if self.rewinder is not None:
            # The rewinder must undo this interleaf's own end point, which is the unrotated end
            # point turned through the interleaf angle.
            end_kx = float(self._kx[-1]) * cos - float(self._ky[-1]) * sin
            end_ky = float(self._kx[-1]) * sin + float(self._ky[-1]) * cos
            rx, ry = self.rewinder
            out.add(
                self.sampling_block_duration,
                _retarget(rx, -end_kx),
                _retarget(ry, -end_ky),
            )
        return out


def _retarget(reference: SimpleNamespace, area: float) -> SimpleNamespace:
    """
    Scale a trapezoid to a new area, keeping its duration.

    Keeping the duration is the point: every interleaf's rewinder occupies the same time, so the
    next TR starts at the same instant regardless of which shot was played.
    """
    factor = 0.0 if reference.area == 0 else area / float(reference.area)
    return ev.derive(
        reference,
        amplitude=float(reference.amplitude) * factor,
        area=float(reference.area) * factor,
        flat_area=float(getattr(reference, 'flat_area', 0.0)) * factor,
    )
