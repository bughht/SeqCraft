"""
Variable-density spiral trajectory design.

Lifted from the deleted ``seqcraft.modules.readout.spiral``.  Depends on numpy and nothing else --
the seqcraft error types have been replaced by ``ValueError`` so this file stands alone.

The classic Glover dual-rate design: slew-limited near the origin, amplitude-limited further out,
with a density that can taper from fully sampled at the centre to undersampled at the edge.  It is
generated analytically and then integrated to a gradient waveform, and the waveform is **measured**
against the limits afterwards rather than assumed to fit.

Why this was worth keeping
--------------------------
Three traps are solved here, each of which silently degrades an image rather than raising:

1. **The waveform is the exact finite difference of the trajectory**, ``g[n] = (k[n] - k[n-1])/dt``,
   so ``cumsum(g) * dt`` reproduces ``k`` to the last bit.  Differentiating analytically instead
   leaves the sampled gradient's cumulative sum short of the analytic radius, because near the
   origin the spiral turns by an eighth of a radian per raster and the polygon cuts the corners.
2. **Never rescale the finished waveform to land on ``k_max``.**  It scales the *turn spacing* with
   it, so a spiral asked for Nyquist sampling delivers 16 % coarser and the aliasing looks like
   undersampling nobody asked for.
3. **The step is chosen against the discrete constraints**, not a continuous-time model of them.
   The continuous model is paid for twice -- in its own error near the origin, and again in the
   global derating needed to bring the realised waveform back inside the limits.  Solving what the
   hardware actually applies leaves the amplifier at 100 % of peak slew and 99.8 % on average,
   which for a slew-limited spiral is a 20 ms readout instead of a 29 ms one.

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
"""

from __future__ import annotations

import math
from typing import NamedTuple

import numpy as np

__all__ = ['Spiral', 'vds_trajectory']


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

    Raises
    ------
    ValueError
        A parameter is non-positive, or no spiral fits the limits given.

    Notes
    -----
    Writes the spiral as ``k(theta) = r(theta) * exp(i*theta)`` with ``r = c * theta**(1/alpha)``
    and walks ``theta`` forward one raster at a time.  Each step is the largest that satisfies the
    constraints the hardware actually applies, which are on the **sampled** waveform::

        g[n] = (k(theta[n]) - k(theta[n-1])) / raster     exact, so cumsum(g) * raster == k

        |g[n]|          <= G_max
        |g[n] - g[n-1]| <= S_max * raster

    capped additionally at the speed the curvature ahead can sustain::

        |dk/dtheta|    = sqrt(r'^2 + r^2)                 tangential scale
        |d2k/dtheta^2| = sqrt((r'' - r)^2 + (2 r')^2)     curvature scale

        theta' <= min( G_max / |dk/dtheta|,  sqrt(S_max / |d2k/dtheta^2|) )

    Neither half is redundant.  Without the cap, a greedy forward pass accelerates to a speed the
    next turn cannot hold, and since shedding speed also costs slew there is then no legal move at
    all: the integration stalls a few radians out.

    Near the origin the curvature cap binds, because the trajectory is turning hard while its
    radius is still small; further out the amplitude limit binds, because ``|dk/dtheta|`` grows
    with the radius.  That crossover is the classic dual-rate spiral.

    The waveform is then ramped down to zero at the slew limit, so it can be followed directly by
    a rewinder with no gradient discontinuity.  That ramp carries ``|k|`` about 12 % past
    ``k_max``; those samples are ordinary k-space points a little beyond the nominal resolution,
    and density compensation takes them in its stride.
    """
    if k_max_per_m <= 0 or n_interleaves < 1 or density <= 0 or alpha <= 0:
        msg = (
            f'spiral parameters must be positive: k_max_per_m={k_max_per_m}, '
            f'n_interleaves={n_interleaves}, density={density}, alpha={alpha}'
        )
        raise ValueError(msg)

    # One interleaf must reach k_max while adjacent turns stay 1/FOV apart, with the winding shared
    # between n_interleaves shots.  `density` then scales the turn count directly: 0.6 winds 60 % of
    # them, which is 40 % less readout and 1/0.6 undersampling at the edge.
    #
    # It multiplies rather than divides.  Dividing -- which is what this did at first -- makes
    # `density=0.6` wind *more* turns and run 66 % longer, the exact opposite of what the name says,
    # and it silently turns the one lever a single-shot spiral has into a penalty.
    theta_max = 2.0 * math.pi * k_max_per_m * fov_m * density / n_interleaves
    if theta_max <= 0:
        msg = 'spiral has no extent; check fov_m and k_max_per_m'
        raise ValueError(msg)

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

    msg = (
        f'could not find a spiral within the given gradient limits '
        f'(k_max_per_m={k_max_per_m}, n_interleaves={n_interleaves}, alpha={alpha}); '
        f'increase n_interleaves, reduce k_max_per_m, or increase density'
    )
    raise ValueError(msg)


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
        tangential = max(math.hypot(first, radius), 1e-12)                 # |dk/dtheta|
        curvature = max(math.hypot(second - radius, 2.0 * first), 1e-12)   # |d2k/dtheta^2|
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
            for factor in (
                1.0, 0.9998, 0.9995, 0.999, 0.998, 0.995, 0.99, 0.98, 0.95, 0.9, 0.7, 0.4,
            ):
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
                msg = f'the spiral cannot continue within the slew limit at theta={theta}'
                raise ValueError(msg)
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
    else:
        msg = (
            f'spiral did not reach k_max within {max_samples} raster steps '
            f'(theta_max={theta_max}, reached={theta})'
        )
        raise ValueError(msg)

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
    change in gradient, which no amplifier can produce.

    The **padding** is subtler and matters just as much.  ``make_arbitrary_grad`` treats waveform
    values as samples at raster *centres*, so ``first`` and ``last`` are the extrapolated values at
    the raster *edges*.  Ending the waveform with a single zero sample does not make that
    extrapolation zero: the ramp's slope carries it past.  Declaring ``last=0`` anyway makes
    pypulseq's shape round-trip disagree with the recorded value by a few tenths of a mT/m, which
    it reports as *"Last restored point differs too much from the recorded last; skipping shape
    restoration"* -- and the restored waveform is then not the one that was designed.  A zero
    sample at each end makes the extrapolation genuinely zero, at a cost of one raster per end.
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
