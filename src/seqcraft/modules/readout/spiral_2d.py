"""
:class:`Spiral2D` -- a variable-density spiral readout, and the trajectory that plays.

``readout/`` for the reason :class:`~seqcraft.modules.EPI2D` is: it contains ADCs and it is a
*leaf*, holding no RF and knowing nothing about what excited the magnetisation it samples.  A
gradient echo and a spin echo differ only in what is played before it, which is what lets
``examples/gre_spiral_2d/`` and ``examples/se_spiral_2d/`` share one readout.

What a spiral asks that a Cartesian readout does not
----------------------------------------------------
A Cartesian readout has to land its samples **on a grid**, and :class:`~seqcraft.modules.EPI2D`
spends most of its arithmetic on that: a window exactly centred in its lobe, an even sample count,
a dwell snapped onto a quantum.  All of it exists because the reconstruction is an FFT and an FFT
assumes positions it was not told.

**A spiral has no grid, and a NUFFT is told the positions.**  So the requirement is different, and
stating it precisely is what keeps this module from importing arithmetic it does not need:

    The trajectory does not have to land anywhere in particular.  It has to be **known** -- the
    array handed back has to be where the samples actually are, to a small fraction of ``dk``.

Weaker in one direction and much stronger in the other.  ``kmax`` coming out 0.033 % low is a
0.033 % scaling of the resolution and nothing else; ``k`` reported one gradient raster from where
it played is a blur, a shading and a wrong field map.  So :meth:`Spiral2D.k_per_m` is **measured
off the built events** -- their own knots, integrated exactly -- and the design pass's ``k`` is
discarded.  ``sc.kspace(tree, opts)['k_adc']`` is the independent oracle, and the two agree to
1.3-4.7e-5 1/m, which is 3-10e-6 of one Nyquist step and is pypulseq's own shape-compression
floor rather than anything a module can improve.

The two things that *do* have to be exact are the **join areas**, because residual ``k`` at the end
of a TR is a phase ramp across the next image, and the **sample times**, because they are what the
off-resonance term multiplies.

The path and the timing are two questions
-----------------------------------------
Conflating them is the commonest way to get a spiral wrong, and the symptom is a readout that
cannot be reversed.  So the design is two passes, and neither imports pypulseq.

**The path** is the Nyquist condition and nothing else.  Adjacent turns of the *interleaved set*
must be no further apart than ``1 / FOV``, so with `shots` interleaves one arm satisfies
``dtheta / dkr = 2 pi FOV(kr) / shots``, integrated from the origin to ``kmax``.  No time appears
in that sentence, which is the point.

**The traversal** is a forward-backward sweep with ``v = 0`` at *both* ends.  Those boundary
conditions are the design decision, not the algorithm: they are what makes the arm reversible, and
reversibility is what makes four variants and any number of echoes one module rather than eight.
They also make the arm *placeable* -- a tree node that begins at 35 mT/m has to be continued by the
previous block.

Then the curve is resampled onto the gradient raster **on** ``k`` and differentiated for ``g``, not
the other way round, and the map from raster time back to arclength uses the exact
uniformly-accelerated law rather than interpolation.  Interpolating is accurate wherever the
arclength cells are short in *time* and wrong at the two ends, where ``v -> 0`` makes one cell span
two rasters: it put 268 T/m/s on a 180 T/m/s amplifier, in exactly two rasters out of 5241, at the
two instants a spiral most obviously ought to be fine.

One waveform, no connectors
---------------------------
Every join -- the prephaser, the fly-backs and the rewinder -- is built **inside the arm's own
arbitrary gradient**, as a triangle rather than a trapezoid.  An arbitrary gradient's samples sit
at raster *centres* and a trapezoid's knots at raster *edges*; the compiler holds either lattice
exactly and their sum on neither, so anything laid beside the arm would be resampled and warned
about.  Built in, there is nothing to superpose anywhere in the readout, and the only block
boundaries are the ones the ADC sample limit forces.

A join's area is also **not the join's alone**: a raster-centre waveform is piecewise linear
through its samples, so the segment running from the previous piece's last sample to this one's
first carries area that a join designed in isolation does not count.  Measured, that left ``k``
0.108 1/m from the origin after a rewinder exact to floating point on its own -- 0.024 ``dk``,
invisible in every plot, and a phase ramp across the next excitation's image.  So the joins are
measured on the *assembled* waveform and corrected by one rescale, which is exact because a
piecewise-linear integral is linear in the amplitudes.

At this protocol ``max_grad`` does nothing
------------------------------------------
The amplitude limit binds only beyond ``max_grad**2 / max_slew`` = 378.5 1/m, and ``kmax`` is
290.9 1/m, so the readout is slew-limited throughout: sweeping ``max_grad`` from 40 to 80 mT/m
moves nothing, and 30 mT/m does.  That is the opposite of the Cartesian intuition and it is why
`max_slew_hz_m_s` is the derating knob that matters -- a spiral slews **continuously at the limit**
for tens of milliseconds rather than in blips, which is the duty cycle peripheral nerve stimulation
responds to hardest.
"""

from __future__ import annotations

import math
import warnings
from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ...design.events import derive, knots_of
from ...design.logic import LogicBlock
from ...design.module import Module
from ...design.timing import Raster
from ...errors import ConfigurationError, SeqCraftWarning, format_error
from .._support import (
    require_axis,
    require_count,
    require_positive,
)
from .._support import segment_samples as _split_samples

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

    from pypulseq.opts import Opts

    from ...design.events import Event

__all__ = ['Spiral2D']

#: The four traversals of one designed arm.  Three are the same waveform played backwards and the
#: fourth is both, which is why they are one module -- see :func:`_plan`.
VARIANTS = ('out', 'in', 'in-out', 'out-in')

#: Points on the designed path, and on the arclength grid the traversal sweeps.  Both are far
#: finer than the gradient raster the answer lands on; the path is integrated once and the sweep
#: is a two-pass loop, so neither is on any hot path.
_PATH_POINTS = 400_000
_SWEEP_POINTS = 40_000


# ------------------------------------------------------------------ exact piecewise-linear parts
def _running_integral(times: np.ndarray, amps: np.ndarray) -> np.ndarray:
    """Exact cumulative integral of a piecewise-linear function, at its own knots."""
    return np.concatenate([[0.0], np.cumsum(0.5 * (amps[1:] + amps[:-1]) * np.diff(times))])


def _integral_at(times: np.ndarray, amps: np.ndarray, query: np.ndarray) -> np.ndarray:
    """
    Exact piecewise-linear integral at arbitrary query times.

    Exact, not sampled.  A PWL function's integral is piecewise *quadratic*, so evaluating it by
    interpolating a densely sampled running sum is wrong by an amount that accumulates over a
    readout with thousands of knots -- 12.5 1/m, 2.8 ``dk``, on the first draft of this arithmetic,
    which is a blur nothing else would have explained.
    """
    running = _running_integral(times, amps)
    idx = np.clip(np.searchsorted(times, query, side='right') - 1, 0, len(times) - 2)
    t0, a0 = times[idx], amps[idx]
    dt = times[idx + 1] - t0
    slope = np.where(dt > 0, (amps[idx + 1] - a0) / np.where(dt > 0, dt, 1.0), 0.0)
    step = query - t0
    return running[idx] + a0 * step + 0.5 * slope * step**2


def _k_of_pair(events: Sequence[Event], query: np.ndarray) -> np.ndarray:
    """``k`` at `query` times, from a built ``(gx, gy)`` pair's own knots.  Complex, in 1/m."""
    out = [_integral_at(*knots_of(event, 0.0), query) for event in events]
    return out[0] + 1j * out[1]


def _pwl_area(waveform: np.ndarray, raster: float) -> float | complex:
    """The exact area of a raster-centre waveform whose ``first`` and ``last`` are zero."""
    times = np.concatenate([[0.0], (np.arange(len(waveform)) + 0.5) * raster,
                            [len(waveform) * raster]])
    amps = np.concatenate([[0.0], waveform, [0.0]])
    return _running_integral(times, amps)[-1]


# --------------------------------------------------------------------- the path: density is data
def _fov_curve(density: Sequence[float] | Callable[[np.ndarray], np.ndarray]) -> Callable:
    """
    Return `density` as a callable of ``u = kr / kmax``, whatever form it arrived in.

    A sequence is **the FOV multiplier sampled at equally spaced radii** -- ``values[i]`` at
    ``u = i / (len(values) - 1)`` -- interpolated linearly, with a single value meaning a constant.
    So every number a caller writes *is* the thing it controls: ``(1.0, 0.5)`` is full FOV at the
    centre and half of it at the edge, and ``(1.0, 1.0, 0.7, 0.33)`` is Nyquist out to a third of
    the radius and then a taper to three times undersampled.

    **Not polynomial coefficients.**  Those are what ``vds.m`` takes and what the literature
    quotes, but they carry a sign the caller has to think about: the curve written ``(1.0, 0.5)``
    here is ``(1.0, -0.5)`` there, because a falling FOV needs negative higher-order terms.  And
    under that convention ``(1.0, 0.5)`` -- the obvious way to write "half at the edge" -- silently
    means a FOV that *grows* to 330 mm, which costs 60 % more readout and reads as a slow design
    rather than as a typo.  Sampled values cannot express that mistake.

    **Linear rather than a monotone spline, and that is measured rather than assumed.**  A kink in
    ``FOV(kr)`` is a jump in the curve's curvature and therefore in the traversal's speed ceiling,
    which sounds like it should cost slew.  It does not: across four sampled densities the peak
    raster slew before any correction is 179.67-179.90 T/m/s against a 180 limit, with linear
    interpolation and with ``scipy``'s ``PchipInterpolator`` alike, and the readouts differ by
    0-5 % -- linear being the *shorter*, because it cuts the corner the spline rounds.  A spline
    would buy a smoother slew profile, no headroom, and a scipy dependency this package does not
    have.  A caller who wants smoothness passes a callable.
    """
    if callable(density):
        return density
    values = np.asarray(density, dtype=float).ravel()
    if values.size == 0:
        msg = format_error('density must hold at least one value.', {'density': density})
        raise ConfigurationError(msg)
    if values.size == 1:
        return lambda u: np.full(np.shape(u), values[0], dtype=float)
    nodes = np.linspace(0.0, 1.0, values.size)
    return lambda u: np.interp(np.asarray(u, dtype=float), nodes, values)


def _check_density(curve: Callable, density: object) -> np.ndarray:
    """
    Return the multiplier on a dense ``u`` grid, refusing one that is not a field of view.

    A ``density`` that reaches zero or below anywhere inside ``u <= 1`` is a trajectory with
    infinite turn spacing, so it is refused naming the radius rather than integrated into a curve
    with a singularity in it.
    """
    grid = np.linspace(0.0, 1.0, 2049)
    values = np.asarray(curve(grid), dtype=float).ravel()
    if values.shape != grid.shape or not np.all(np.isfinite(values)):
        bad = grid[~np.isfinite(values)] if values.shape == grid.shape else grid[:1]
        where = f'u = {float(bad[0]):.3f}' if bad.size else 'the wrong shape'
        msg = format_error(
            f'density must return one finite value per radius; it did not at {where}.',
            {'density': density},
            ['a callable takes an array of u = kr/kmax and returns the FOV multiplier there'],
        )
        raise ConfigurationError(msg)
    if np.any(values <= 0.0):
        radius = float(grid[np.argmax(values <= 0.0)])
        msg = format_error(
            f'density reaches {float(values.min()):.3g} at u = {radius:.3f}, so the field of view '
            f'is zero or negative there.',
            {'density': density},
            ['density is the FOV multiplier, so every value must be above zero',
             '(1.0, 0.5) halves the FOV at the edge; (1.0, -0.5) is the polynomial convention '
             'and is not this one'],
        )
        raise ConfigurationError(msg)
    return values


def _design_path(*, fov_m: float, matrix: int, shots: int, curve: Callable,
                 points: int = _PATH_POINTS) -> tuple[np.ndarray, np.ndarray]:
    """
    Return ``(k, s)``: the spiral as complex 1/m, and its arclength, on a dense grid.

    No time, no gradient, no scanner -- and no pypulseq import in this function or the next, which
    is the structural claim of the whole module.  Written the other way round, as one loop
    producing gradients, the module works for ``'out'`` and cannot be reversed.

    The arc length has a closed form, ``integral 2 pi FOV(kr) kr dkr``, and it is **independent of
    `shots`**: interleaving divides the same curve rather than shortening it.
    """
    kmax = matrix / (2.0 * fov_m)
    kr = np.linspace(0.0, kmax, points)
    fov_of_kr = fov_m * np.asarray(curve(kr / kmax), dtype=float)
    integrand = 2.0 * np.pi * fov_of_kr / shots
    theta = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(kr))])
    k = kr * np.exp(1j * theta)
    return k, np.concatenate([[0.0], np.cumsum(np.abs(np.diff(k)))])


def _traverse(k: np.ndarray, s: np.ndarray, *, max_grad: float, max_slew: float,
              points: int = _SWEEP_POINTS) -> tuple[np.ndarray, ...]:
    """
    Return ``(grid, v, t, ds)``: how fast the curve may be taken, and when each point is reached.

    A forward-backward sweep on a uniform arclength grid, with ``v = 0`` at *both* ends.  The
    gradient is the velocity along the curve and the slew is the acceleration, so a corner of
    curvature ``kappa`` taken at speed ``v`` already spends ``kappa v**2`` of the slew budget on
    turning: the tangential acceleration left over is ``sqrt(max_slew**2 - (kappa v**2)**2)`` and
    the ceiling on speed is ``min(max_grad, sqrt(max_slew / kappa))``.

    **A curve cannot decelerate while it is on the curvature ceiling**, because there the turning
    term is the whole budget and the tangential term is exactly zero.  That is why the closed form
    ``(2 pi / sqrt(S)) integral FOV sqrt(kr) dkr``, which assumes the ceiling is ridden everywhere,
    is a genuine *lower bound* rather than an estimate -- tight to 0.4 % at the reference protocol,
    so a design that comes out below it has a bug.

    Returns the speed profile rather than a resampled curve, because the map from time back to
    arclength has to use the exact within-cell law and a function that hands back an interpolated
    curve has already thrown that away.
    """
    grid = np.linspace(0.0, float(s[-1]), points)
    curve = np.interp(grid, s, k.real) + 1j * np.interp(grid, s, k.imag)
    ds = float(grid[1] - grid[0])

    second = np.empty_like(curve)
    second[1:-1] = (curve[2:] - 2 * curve[1:-1] + curve[:-2]) / ds**2
    second[0], second[-1] = second[1], second[-2]
    kappa = np.maximum(np.abs(second), 1e-12)
    ceiling = np.minimum(max_grad, np.sqrt(max_slew / kappa))

    def sweep(limit: np.ndarray, curvature: np.ndarray) -> np.ndarray:
        v = np.empty_like(limit)
        v[0] = 0.0
        budget = max_slew**2
        for i in range(len(limit) - 1):
            turning = curvature[i] * v[i] ** 2
            tangential = math.sqrt(max(budget - turning * turning, 0.0))
            v[i + 1] = min(limit[i + 1], math.sqrt(v[i] ** 2 + 2.0 * tangential * ds))
        return v

    forward = sweep(ceiling, kappa)
    backward = sweep(ceiling[::-1], kappa[::-1])[::-1]   # the curvature reversed *with* the limit
    v = np.minimum(forward, backward)
    return grid, v, np.concatenate([[0.0], np.cumsum(2.0 * ds / (v[1:] + v[:-1]))]), ds


# ------------------------------------------------------- the raster: on k, not on g
def _k_on_raster(k, s, grid, v, t, ds, *, raster: float,
                 stretch: float = 1.0) -> tuple[np.ndarray, float]:
    """
    Return ``(k at the raster edges, the arm's duration)``.

    Resampled on ``k`` and differentiated for ``g``, not the other way round: ``k`` is what the
    reconstruction needs and what the check measures, and a uniform resample of ``g`` rounds a
    spiral's peaks off by 2.5 % while preserving its area exactly.

    **Time maps back to arclength by the exact uniformly-accelerated law.**  Within one arclength
    cell the sweep built ``v**2(sigma) = v_i**2 + 2 a sigma`` exactly, so ``sigma = v_i tau +
    a tau**2 / 2``.  Interpolating ``k`` against ``t`` instead is accurate everywhere the cells are
    short in *time* and wrong at the two ends, where ``v -> 0`` makes one cell span two rasters.

    `stretch` lengthens the whole arm: a time stretch of ``f`` divides the gradient by ``f`` and
    the slew by ``f**2``, so one rescale against a measured overshoot is exact.
    """
    count = int(math.ceil(t[-1] * stretch / raster))
    duration = count * raster
    times = t * (duration / t[-1])
    speeds = v * (t[-1] / duration)
    query = np.arange(count + 1) * raster

    cell = np.clip(np.searchsorted(times, query, side='right') - 1, 0, len(times) - 2)
    tau = query - times[cell]
    accel = (speeds[cell + 1] ** 2 - speeds[cell] ** 2) / (2.0 * ds)
    arc = np.clip(grid[cell] + speeds[cell] * tau + 0.5 * accel * tau**2, 0.0, float(s[-1]))
    return np.interp(arc, s, k.real) + 1j * np.interp(arc, s, k.imag), duration


def _rasterize(k, s, sweep, *, raster: float, max_grad: float,
               max_slew: float) -> tuple[np.ndarray, float]:
    """
    Return ``(g at the raster centres, the arm's duration)``, verified against the amplifier.

    Measured rather than trusted, and corrected by **one** rescale.  The sweep respects the limits
    on the continuous curve; the raster representation is a chord approximation of it and can come
    out a fraction of a per cent over.  Stretching time by ``sqrt(overshoot)`` fixes that exactly,
    so the loop below runs once in practice and is bounded rather than iterative in spirit.

    The two **half-raster edges** count too: an arbitrary gradient's ``first`` and ``last`` sit at
    the interval edges and its samples at the centres, so the outermost slope spans ``raster / 2``
    and a check that looks only at ``diff(g)`` misses it by a factor of two.
    """
    stretch, g, duration = 1.0, np.zeros(0, dtype=complex), 0.0
    for _ in range(8):
        k_edge, duration = _k_on_raster(k, s, *sweep, raster=raster, stretch=stretch)
        g = np.diff(k_edge) / raster
        edge_slew = 2.0 * max(abs(g[0]), abs(g[-1])) / raster
        over = max(float(np.abs(np.diff(g)).max() / raster / max_slew),
                   float(edge_slew / max_slew),
                   float(np.abs(g).max() / max_grad))
        if over <= 1.0:
            break
        stretch *= math.sqrt(over) * (1 + 1e-9)
    return g, duration


# ------------------------------------------------------------------- the joins, inside the arm
def _ramp_rasters(area: float, *, raster: float, max_grad: float, max_slew: float) -> int:
    """
    The fewest gradient rasters a **triangular** 2D join of `area` needs.

    A triangle, not a trapezoid: at fixed duration and area, lengthening the ramps at the expense
    of the flat top lowers the slew and raises the peak, and a spiral's joins are slew-limited
    rather than amplitude-limited.  ``n = 2 sqrt(A / (S r**2))`` is that solved, and for
    ``A = kmax`` it comes out at 39 rasters and 35.00 mT/m -- ``sqrt(max_slew * kmax)``, the same
    number as the spiral's own peak, and not by coincidence.
    """
    by_amplitude = 2.0 * area / (max_grad * raster)
    by_slew = 2.0 * math.sqrt(area / (max_slew * raster**2))
    return max(2, int(math.ceil(max(by_amplitude, by_slew))))


def _ramp_to(area: complex, *, rasters: int, raster: float) -> np.ndarray:
    """
    A 2D triangular join of exactly `area`, at raster centres, rescaled so the area is exact.

    **Sampled at the raster centres a raster-centre waveform means**, which is not the same as a
    ``linspace``: an arbitrary gradient's ``first`` sits at ``t = 0`` and its first sample at
    ``t = raster / 2``, so a shape whose first sample is already the full ramp step slews *twice*
    as hard over that half raster as it does anywhere else.  Built the obvious way, a join sized
    for ``kmax`` asked 390 T/m/s of a 180 T/m/s amplifier, in one half-raster at each end.

    Rescaled rather than trusted: a PWL waveform's area is exactly linear in its amplitude at fixed
    shape, so one rescale against the measured integral is exact.
    """
    centres = (np.arange(rasters) + 0.5) / rasters
    shape = 1.0 - np.abs(2.0 * centres - 1.0)
    achieved = _pwl_area(shape, raster)
    return area * shape / achieved if achieved else np.zeros(rasters, dtype=complex)


def _correct_joins(waveform: np.ndarray, joins: Sequence[tuple[int, int, int, complex]], *,
                   raster: float) -> np.ndarray:
    """
    Rescale each join so ``k`` reaches its target exactly, on the **assembled** waveform.

    A join's area is not the join's alone: a raster-centre waveform is piecewise linear through its
    samples, so the segment running from the previous piece's last sample to this one's first
    carries area too, and a join designed in isolation with ``first`` and ``last`` at zero does not
    count it.  Measured, a rewinder exact to floating point on its own left ``k`` 0.108 1/m from
    the origin -- 0.024 ``dk``, invisible in every plot, and a phase ramp across the next
    excitation's image.

    Each join carries a **checkpoint**: the sample index at which its target is measured.  For
    all but the last that is the join's own end, which is where the next arm begins.  For the
    last it is the end of the block, because a residual there is a phase ramp across the *next
    excitation's* image and there is nothing after it to absorb one -- and half a raster of the
    arm's own bracket, 0.02 ``dk``, is exactly what lands there otherwise.

    One rescale is enough, and that is not luck: a PWL integral is linear in the amplitudes, so
    scaling one join by ``(area + shortfall) / area`` moves the checkpoint by exactly the
    shortfall and nothing upstream of it moves at all.  Corrections are applied in time order
    for the same reason -- each is measured against a waveform whose earlier joins are right.
    """
    if not joins:
        return waveform
    times = np.concatenate([[0.0], (np.arange(len(waveform)) + 0.5) * raster,
                            [len(waveform) * raster]])
    for start, length, checkpoint, target in joins:
        amps = np.concatenate([[0.0], waveform, [0.0]])
        end = start + length
        actual = complex(_running_integral(times, amps.real)[checkpoint + 1],
                         _running_integral(times, amps.imag)[checkpoint + 1])
        area = complex(_pwl_area(waveform[start:end].real, raster),
                       _pwl_area(waveform[start:end].imag, raster))
        if area != 0:
            waveform[start:end] *= (area + (target - actual)) / area
    return waveform


def _plan(variant: str, echoes: int) -> tuple[list[int], str]:
    """
    Return ``(the direction of each arm in order, what goes between echoes)``.

    The two one-arm variants need a fly-back between echoes -- a rewinder for ``'out'``, a
    prephaser for ``'in'`` -- exactly as a **monopolar** multi-echo Cartesian train does.  The two
    two-arm variants need nothing at all, because both of their joins fall where ``g = 0``, which
    is the **bipolar** case.  That parallel is why `echoes` needs no `polarity` beside it:
    `variant` already carries it.

    ``'out-in'`` acquires ``k = 0`` at each *end* rather than in the middle, so it has an echo at
    the start and an echo at the finish and cannot have fewer than two.
    """
    if variant == 'out':
        return [+1] * echoes, 'flyback'
    if variant == 'in':
        return [-1] * echoes, 'flyback'
    if variant == 'in-out':
        return [-1, +1] * echoes, 'none'
    if echoes < 2:
        msg = format_error(
            "variant='out-in' acquires k = 0 at each end, so it has at least two echoes.",
            {'variant': variant, 'echoes': echoes},
            ["echoes=1 is variant='out'", "echoes=2 is the cheapest way to get exactly two"],
        )
        raise ConfigurationError(msg)
    return [+1, -1] * (echoes - 1), 'none'


class Spiral2D(Module):
    """
    A variable-density spiral readout: the path, the traversal, the ADCs, the joins, the labels.

    Parameters
    ----------
    opts
        The scanner.
    fov_mm
        The field of view **at ``k = 0``**, not an average.  With ``density=(1.0,)`` it is the FOV
        everywhere; any other density only ever reduces it.  Scalar rather than a pair, because a
        spiral's density is radial.
    matrix
        Square and in-plane.  Sets ``kmax = matrix / (2 fov_m)`` and therefore the resolution.  It
        does **not** set a sample count -- that is `dwell_s` and the traversal time.
    shots
        Interleaves.  ``build(shot=n)`` rotates the arm by ``2 pi n / shots``.  **Total readout
        time and total sample count barely depend on it** -- 52.45 ms at one shot against 53.68 ms
        at eight, a 2.3 % spread that is entirely the extra accelerate-and-stop at each arm's ends.
        What `shots` buys is a shorter *arm*, which is off-resonance blur, T2\\* decay and the ADC
        sample limit: three reasons, one number.
    density
        ``FOV(kr) / fov_mm``, **sampled at equally spaced radii and interpolated linearly**.
        ``(1.0,)`` is Archimedean, ``(1.0, 0.5)`` halves the FOV at the edge, and
        ``(1.0, 1.0, 0.7, 0.33)`` holds Nyquist to a third of the radius and then tapers.  Values,
        not polynomial coefficients -- :func:`_fov_curve` is why.  A callable is accepted for a
        density that came out of an optimiser, which is data too.
    variant
        ``'out'``, ``'in'``, ``'in-out'`` or ``'out-in'`` -- see :func:`_plan`.
    echoes
        The number of acquired ``k = 0`` crossings, which is the length of :attr:`te_s`, exactly as
        :class:`~seqcraft.modules.CartesianLine`'s `echoes` is the length of its.
    echo_spacing_s
        Lengthen the joins to reach a requested delta-TE.  For the two-arm variants that means
        inserting dead time at a ``g = 0`` join, which costs nothing but time.  ``None`` is the
        minimum.  Raises below it, naming it, and raises at ``echoes=1``, where there is no
        spacing.
    dwell_s, bandwidth_hz
        The sampling rate, either way round; **exactly one is required**.  `bandwidth_hz` is the
        whole receive bandwidth ``1 / dwell``, not a per-pixel figure: a spiral has no readout
        direction, so "Hz per pixel" would be a number with no referent.
    max_grad_hz_m, max_slew_hz_m_s
        Design against less than the amplifier can do, without building a second ``Opts``.  The
        usual reason is peripheral nerve stimulation.
    adc_segments
        The most ADC events one acquisition region may be split into.  ``None`` is however many the
        interpreter's sample limit needs; ``1`` refuses rather than splitting, naming the smallest
        `shots` that would fit.  A split costs a **hole in the trajectory** -- 10 us of dead time
        at 35 mT/m is 14.9 1/m, 3.28 ``dk`` -- so this is a *recovery* and `shots` is the answer.
    prephase, rewind
        Whether the module supplies the ramp that gets ``k`` to the start of the arm, and the one
        that brings it back to zero.  ``None`` is "whatever this variant needs".  ``False`` is the
        **spin-echo placement**: the dephaser goes before the refocusing pulse with the sign the
        conjugation will flip, where the interval already exists and it costs no echo time.  An
        explicit ``True`` on a variant that already starts (or ends) at ``k = 0`` raises, because
        there is nothing there to prephase or rewind.
    axes
        The two logical gradient channels.  They must differ.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    k_per_m : method
        ``(2, n)`` in 1/m -- **where the samples are**, measured off the built events rather than
        off the design.  A method because it depends on the shot; everything else here is a
        constant of the instance, which is the point.
    t_adc_s : numpy.ndarray
        ``(n,)`` sample times from the block start, in acquisition order.
    te_s : tuple of float
        One time per acquired ``k = 0``, from the block start.  Always a tuple, exactly like
        :attr:`CartesianLine.te_s`: a name that is a float for three variants and a tuple for one
        is a name every caller has to branch on.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s',
    ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
    ...          adc_samples_limit=8192, B0=3.0)
    >>> four = Spiral2D(opts=o, fov_mm=220.0, matrix=128, shots=4, dwell_s=2.5e-6)
    >>> round(four.kmax_per_m, 3), round(four.dk_per_m, 4), round(four.arm_duration_s * 1e3, 2)
    (290.909, 4.5455, 13.23)
    >>> four.num_samples, four.adc_segments, round(four.turns, 1)
    (5292, 1, 16.0)

    The readout is slew-limited throughout at this protocol, so ``max_grad`` does nothing:

    >>> harder = pp.Opts(max_grad=80, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s',
    ...                  rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> Spiral2D(opts=harder, fov_mm=220.0, matrix=128, shots=4,
    ...          dwell_s=2.5e-6).arm_duration_s == four.arm_duration_s
    True

    A taper is a shorter arc and therefore a shorter readout, at the same ``kmax``:

    >>> vds = Spiral2D(opts=o, fov_mm=220.0, matrix=128, shots=4, dwell_s=2.5e-6,
    ...                density=(1.0, 0.5))
    >>> round(vds.arm_duration_s / four.arm_duration_s, 3), vds.edge_undersampling
    (0.703, 2.0)

    ``echoes`` is the number of acquired ``k = 0`` crossings, and the variant decides what goes
    between them -- a fly-back for the one-arm variants, nothing at all for the two-arm ones:

    >>> me = Spiral2D(opts=o, fov_mm=220.0, matrix=128, shots=4, dwell_s=2.5e-6, echoes=3)
    >>> len(me.te_s), round(me.echo_spacing_s * 1e3, 2)
    (3, 13.62)
    >>> inout = Spiral2D(opts=o, fov_mm=220.0, matrix=128, shots=4, dwell_s=2.5e-6,
    ...                  variant='in-out', prephase=False)
    >>> len(inout.te_s), round(float(inout.k_start_per_m(0)[0]), 2)
    (1, 290.81)

    Notes
    -----
    **Exactly one of ``shot=`` and ``angle_rad=`` is required by** :meth:`build`.  A default of
    ``shot=0`` was rejected because a four-shot acquisition that silently acquires the same
    interleaf four times is a plausible bug with a plausible-looking image.

    **Nothing may overlap the spiral on either axis inside its own block.**  An arbitrary
    gradient's samples sit at raster centres and a trapezoid's knots at raster edges, and the
    compiler holds either lattice exactly and their sum on neither -- so it resamples and warns.
    Since the module builds one event per axis for the whole readout, there is no sum to take
    anywhere inside it; what a *caller* adds is the caller's to keep clear.
    """

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: float,
        matrix: int,
        shots: int = 1,
        density: Sequence[float] | Callable = (1.0,),
        variant: str = 'out',
        echoes: int = 1,
        echo_spacing_s: float | None = None,
        dwell_s: float | None = None,
        bandwidth_hz: float | None = None,
        max_grad_hz_m: float | None = None,
        max_slew_hz_m_s: float | None = None,
        adc_segments: int | None = None,
        prephase: bool | None = None,
        rewind: bool | None = None,
        axes: tuple[str, str] = ('x', 'y'),
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.fov_mm = require_positive(fov_mm, 'fov_mm')
        self.matrix = require_count(matrix, 'matrix', low=2)
        self.shots = require_count(shots, 'shots', hint='interleaves; 1 is single-shot')
        self.echoes = require_count(echoes, 'echoes', hint='acquired k = 0 crossings')
        self.variant = self._check_variant(variant)
        self.axes = self._check_axes(axes)
        self.density = density
        self._fov_of = _fov_curve(density)
        _check_density(self._fov_of, density)

        self.kmax_per_m = self.matrix / (2.0 * self.fov_mm * 1e-3)
        self.dk_per_m = 1e3 / self.fov_mm
        self.resolution_mm = self.fov_mm / self.matrix
        #: ``R`` at the edge of k-space, ``1 / density[-1]`` -- 1.0 at uniform density.
        self.edge_undersampling = 1.0 / float(self._fov_of(np.array(1.0)))

        self.raster_s = float(opts.grad_raster_time)
        self.dwell_s = self._resolve_dwell(dwell_s, bandwidth_hz)
        #: The whole receive bandwidth, ``1 / dwell_s``, read back after the raster snap.
        self.bandwidth_hz = 1.0 / self.dwell_s
        self.max_grad_hz_m_limit = min(float(opts.max_grad),
                                       float(max_grad_hz_m or math.inf))
        self.max_slew_hz_m_s_limit = min(float(opts.max_slew),
                                         float(max_slew_hz_m_s or math.inf))

        directions, join = _plan(self.variant, self.echoes)
        self._directions, self._join = directions, join
        self.prephase, self.rewind = self._resolve_ends(prephase, rewind)

        self._design_arm()
        self._assemble(echo_spacing_s, adc_segments)
        self._measure()

    # ------------------------------------------------------------------------- design
    def _design_arm(self) -> None:
        """Design one arm: the path, the traversal, the raster waveform, and what they measure."""
        path, arclength = _design_path(
            fov_m=self.fov_mm * 1e-3, matrix=self.matrix, shots=self.shots, curve=self._fov_of)
        #: Arc length of **one arm**, 1/m.  ``shots * arc_length_per_m`` is what is independent
        #: of `shots`: interleaving divides one curve rather than shortening it, which is why the
        #: total readout time barely moves with `shots` while the arm is inversely proportional.
        self.arc_length_per_m = float(arclength[-1])
        angle = np.unwrap(np.angle(path[1:]))
        #: Revolutions in one arm.
        self.turns = float(angle[-1] / (2 * np.pi))
        #: Widest radial gap between adjacent turns of the *interleaved set*, in ``dk``.
        radius = np.abs(path[1:])
        previous = np.interp(angle - 2 * np.pi / self.shots, angle, radius, left=np.nan)
        self.worst_turn_gap_dk = float(np.nanmax(radius - previous) / self.dk_per_m)

        sweep = _traverse(path, arclength, max_grad=self.max_grad_hz_m_limit,
                          max_slew=self.max_slew_hz_m_s_limit)
        self._arm_g, self.arm_duration_s = _rasterize(
            path, arclength, sweep, raster=self.raster_s,
            max_grad=self.max_grad_hz_m_limit, max_slew=self.max_slew_hz_m_s_limit)

    # ----------------------------------------------------------------------- assembly
    def _assemble(self, echo_spacing_s: float | None, adc_segments: int | None) -> None:
        """
        Lay out one whole shot at angle zero: joins, arms, pads, and where the ADCs open.

        **One continuous arbitrary gradient for the entire readout**, and one design for every
        shot: rotating a complex waveform *is* rotating the 2D gradient, and every step here --
        the triangular joins, the area correction, the k measurement -- is linear in it.  So the
        interleaf is a multiplication at build time rather than a second design, and an arm whose
        length changed with the interleaf (which would make every timing query below depend on the
        shot) is not expressible.
        """
        raster, directions = self.raster_s, self._directions
        #: Un-sampled time at each end of an acquisition region, on the gradient raster so that
        #: every ADC seam is an instant the compiler may cut a block at.
        self.guard_s = float(Raster(raster).ceil(float(self.opts.adc_dead_time)))
        guard_n = int(round(self.guard_s / raster))

        arm_end_k = complex(_pwl_area(self._arm_g.real, raster), _pwl_area(self._arm_g.imag, raster))
        join_n = _ramp_rasters(abs(arm_end_k), raster=raster, max_grad=self.max_grad_hz_m_limit,
                               max_slew=self.max_slew_hz_m_s_limit)
        pad_n = self._pad_rasters(echo_spacing_s, join_n)

        pieces: list[np.ndarray] = []
        joins: list[list] = []
        regions: list[tuple[int, int]] = []      # (first raster, rasters) of each acquisition
        echo_at: list[float] = []                # nominal k = 0 times, seconds from block start
        length = 0                               # rasters placed so far
        # Every target below is a k measured from wherever the *block* starts.  With
        # prephase=False -- where the caller supplies kmax and the module's own integral still
        # begins at zero -- they are offset by the k the caller is being asked for.
        self._caller_prephases = directions[0] < 0 and not self.prephase
        origin = arm_end_k if self._caller_prephases else 0.0j

        def place(piece: np.ndarray, *, target: complex | None = None) -> None:
            nonlocal length
            if target is not None:
                joins.append([length, len(piece), length + len(piece), target - origin])
            pieces.append(piece)
            length += len(piece)

        if directions[0] < 0 and self.prephase:
            place(_ramp_to(arm_end_k, rasters=join_n, raster=raster), target=arm_end_k)
        else:
            place(np.zeros(guard_n, dtype=complex))     # lead-in, so the ADC has its dead time

        region_start = length
        for index, direction in enumerate(directions):
            if direction > 0:
                echo_at.append(length * raster)
            place(self._arm_g if direction > 0 else -self._arm_g[::-1])
            if direction < 0:
                echo_at.append(length * raster)
            last = index + 1 == len(directions)
            if last or self._join == 'flyback':
                regions.append((region_start, length - region_start))
            if last:
                break
            if self._join == 'flyback':
                target = 0.0j if direction > 0 else arm_end_k
                source = arm_end_k if direction > 0 else 0.0j
                place(_ramp_to(target - source, rasters=join_n + pad_n, raster=raster),
                      target=target)
                region_start = length
            elif direction > 0 and pad_n:
                place(np.zeros(pad_n, dtype=complex))   # dead time at kmax, where g is zero

        if directions[-1] > 0 and self.rewind:
            place(_ramp_to(-arm_end_k, rasters=join_n, raster=raster), target=0.0 + 0.0j)
        else:
            place(np.zeros(guard_n, dtype=complex))     # tail, for the last ADC's dead time

        # Where the block has to leave k: the origin whenever the readout ends there or rewinds
        # itself, and kmax when the caller took the rewinder on.  The last join is retargeted
        # onto it, so the half-raster residual falls at the *start* of the readout -- a 0.007 %
        # scaling of kmax that `k_per_m` reports -- rather than at its end, where nothing could
        # absorb it and it would be a phase ramp across the next excitation's image.
        self._end_target = 0.0j if (directions[-1] < 0 or self.rewind) else arm_end_k
        if joins:
            joins[-1][2:] = [length, self._end_target - origin]
        self._waveform = _correct_joins(np.concatenate(pieces), joins, raster=raster)
        self.readout_duration_s = len(self._waveform) * raster
        self._nominal_echo_s = tuple(echo_at[i] for i in self._echo_arms())
        self._layout_adc(regions, adc_segments)

    def _echo_arms(self) -> tuple[int, ...]:
        """Which entries of the per-arm ``k = 0`` times are the acquired echoes."""
        if self.variant in ('out', 'in'):
            return tuple(range(self.echoes))
        if self.variant == 'in-out':
            return tuple(range(0, 2 * self.echoes, 2))          # each r|a seam
        return (0, *range(1, 2 * (self.echoes - 1), 2))

    def _pad_rasters(self, echo_spacing_s: float | None, join_n: int) -> int:
        """Extra rasters at each join, so a requested echo spacing is reached exactly."""
        arms = 1 if self._join == 'flyback' else 2
        minimum = arms * self.arm_duration_s + (join_n * self.raster_s if arms == 1 else 0.0)
        if echo_spacing_s is None:
            return 0
        if self.echoes < 2:
            msg = format_error(
                'echo_spacing_s needs at least two echoes to be the spacing between.',
                {'echoes': self.echoes, 'echo_spacing_s': float(echo_spacing_s)},
                ['raise echoes, or drop echo_spacing_s'],
            )
            raise ConfigurationError(msg)
        wanted = require_positive(echo_spacing_s, 'echo_spacing_s')
        if wanted < minimum - 1e-12:
            msg = format_error(
                f'echo_spacing_s = {wanted * 1e3:.3f} ms is shorter than the '
                f'{minimum * 1e3:.3f} ms this variant can reach.',
                {'echo_spacing_s': wanted, 'variant': self.variant},
                [f'the minimum is {minimum * 1e3:.3f} ms',
                 'raise shots or taper density to shorten the arm'],
            )
            raise ConfigurationError(msg)
        return int(round(Raster(self.raster_s).ceil(wanted - minimum) / self.raster_s))

    def _layout_adc(self, regions: Sequence[tuple[int, int]], adc_segments: int | None) -> None:
        """
        Decide how each acquisition region is sampled, and refuse a split the caller ruled out.

        **An acquisition region is not an arm.**  Arms that join at ``g = 0`` are one continuous
        gradient with nothing between them, so one ADC may span both -- and that is what puts an
        in-out readout's ``k = 0`` *inside* a sampling window instead of in the guard between two.
        Built per arm, the sample nearest the echo lands a whole guard away; built per region, it
        lands half a dwell from the seam, at ``|k|`` of 0.02 ``dk``.
        """
        span_s = regions[0][1] * self.raster_s
        counts = self._fit_region(span_s)
        if not counts:
            msg = format_error(
                f'no ADC samples fit in a {span_s * 1e6:.1f} us acquisition at a '
                f'{self.dwell_s * 1e9:.0f} ns dwell.',
                {'dwell_s': self.dwell_s, 'shots': self.shots},
                ['lower shots to lengthen the arm, or shorten dwell_s'],
            )
            raise ConfigurationError(msg)
        if adc_segments is not None and len(counts) > require_count(adc_segments, 'adc_segments'):
            self._refuse_segments(counts, int(adc_segments))

        #: One per ADC event of one region; they need not be equal.
        self.segment_samples = counts
        #: How many ADC events one acquisition region takes.
        self.adc_segments = len(counts)
        #: Samples in one acquisition region, summed over its ADC events.
        self.num_samples = sum(counts)

        # Aligned to the end of the region that carries a k = 0, so that echo is sampled exactly.
        # What is left over is at most one sample quantum -- zero at the reference protocol, where
        # four dwells make one gradient raster -- and it falls at the other end.
        spare = span_s - (self.num_samples * self.dwell_s
                          + (len(counts) - 1) * 2 * self.guard_s)
        offset = spare if self._directions[-1] < 0 else 0.0
        self._spans: list[tuple[float, int]] = []
        for start, _ in regions:
            when = start * self.raster_s + offset
            for count in counts:
                self._spans.append((when, count))
                when += count * self.dwell_s + 2 * self.guard_s

    def _fit_region(self, span_s: float) -> tuple[int, ...]:
        """
        How many samples fit in one acquisition of `span_s`, and how they split across events.

        Two costs, and the second is the one nobody predicts.  A **guard** of one
        ``adc_dead_time`` at each end, because the receiver needs it before it opens and after it
        closes and the gradient does not stop for either.  And a gap of **two** guards between
        consecutive ADC events of the same acquisition, for the same reason applied twice.  Both
        shorten the acquisition without shortening the gradient, so the sample count is not
        ``duration / dwell`` and never was.
        """
        limit = int(getattr(self.opts, 'adc_samples_limit', 0) or 0)
        capacity = span_s
        counts: tuple[int, ...] = ()
        for events in range(1, 8):
            usable = capacity - (events - 1) * 2 * self.guard_s
            counts = _split_samples(int(usable / self.dwell_s + 1e-9), dwell_s=self.dwell_s,
                                    limit=limit, opts=self.opts)
            if len(counts) <= events:
                return counts
        return counts

    def _refuse_segments(self, counts: tuple[int, ...], allowed: int) -> None:
        """Refuse a split the caller ruled out, naming the `shots` that would not need one."""
        limit = int(getattr(self.opts, 'adc_samples_limit', 0) or 0)
        per_event = limit * self.dwell_s + 2 * self.guard_s
        enough = max(self.shots + 1,
                     int(math.ceil(self.shots * self.arm_duration_s / per_event)))
        extra = ['a split costs one adc_dead_time of trajectory at each side of the seam -- '
                 f'{self.max_grad_hz_m_limit * float(self.opts.adc_dead_time) / self.dk_per_m:.1f}'
                 ' dk at this protocol'] if allowed == 1 else []
        msg = format_error(
            f'this acquisition needs {len(counts)} ADC events of {counts} samples against '
            f'adc_segments={allowed}; adc_samples_limit is {limit}.',
            {'adc_segments': allowed, 'shots': self.shots, 'adc_samples_limit': limit},
            [f'shots={enough} makes one arm fit in a single ADC event',
             'or raise adc_segments and accept the seam', *extra],
        )
        raise ConfigurationError(msg)

    # ---------------------------------------------------------------------- measurement
    def _measure(self) -> None:
        """Build the angle-zero events and read the trajectory, the echoes and the peaks off them."""
        pair = self._pair_at(0.0)
        # The whole block's area, exactly.  With prephase=False the caller's dephaser has to be
        # its negation, or the readout does not finish where it started -- a stronger and
        # simpler statement than 'supply kmax', and the only one that stays true once the joins
        # have been corrected against the assembled waveform.
        total = complex(_k_of_pair(pair, np.array([self.readout_duration_s]))[0])
        self._k_start = (self._end_target - total) if self._caller_prephases else 0.0j
        self._k_end = self._k_start + total

        self.t_adc_s = np.concatenate(
            [start + (np.arange(count) + 0.5) * self.dwell_s for start, count in self._spans])
        k = _k_of_pair(pair, self.t_adc_s) + self._k_start
        self._k0 = np.vstack([k.real, k.imag])

        magnitude = np.abs(k)
        half = self.arm_duration_s / 2
        samples = []
        for when in self._nominal_echo_s:
            window = np.flatnonzero(np.abs(self.t_adc_s - when) <= half)
            if window.size == 0:
                window = np.array([int(np.argmin(np.abs(self.t_adc_s - when)))])
            samples.append(int(window[np.argmin(magnitude[window])]))
        self._echo_samples = tuple(samples)
        #: One time per acquired ``k = 0``, from the block start.
        self.te_s = tuple(float(self.t_adc_s[n]) for n in samples)

        # On the events' own knots, which is where a PWL function's slope is piecewise constant
        # -- and which includes the two **half-raster edges**, where `first` and `last` sit and
        # where a check that looks only at the samples is out by a factor of two.
        times, amps = knots_of(pair[0], 0.0)
        _, amps_y = knots_of(pair[1], 0.0)
        #: The peak the amplifier is actually asked for, Hz/m -- joins included.
        self.peak_grad_hz_m = float(np.hypot(amps, amps_y).max())
        #: The peak **vector** slew, Hz/m/s -- not a per-axis one.
        self.peak_slew_hz_m_s = float(
            (np.abs(np.diff(amps) + 1j * np.diff(amps_y)) / np.diff(times)).max())

        gaps = []
        at = 0
        for _, count in self._spans:
            gaps.append(np.abs(np.diff(k[at:at + count])).max())
            at += count
        #: Widest ``|dk|`` between consecutive samples of one ADC event, in ``dk``.
        self.worst_sample_gap_dk = float(max(gaps) / self.dk_per_m)
        self._warn_if_undersampled()

        self._adc = {count: pp.make_adc(num_samples=count, dwell=self.dwell_s, system=self.opts)
                     for count in set(self.segment_samples)}

    def _warn_if_undersampled(self) -> None:
        """
        Warn -- not refuse -- when consecutive samples are further apart than one ``dk``.

        A warning because deliberate along-arm undersampling is a real choice and the module should
        not have an opinion about it: the along-arm direction is oversampled relative to the radial
        one on every spiral, and a caller trading it away is doing something coherent.  What is not
        a real choice is *not knowing*, so the message names the longest dwell that stays inside
        one ``dk``.
        """
        if self.worst_sample_gap_dk <= 1.0:
            return
        longest = self.dwell_s / self.worst_sample_gap_dk
        warnings.warn(
            f'consecutive samples are up to {self.worst_sample_gap_dk:.2f} dk apart at a '
            f'{self.dwell_s * 1e9:.0f} ns dwell, so the arm is undersampled along its own path; '
            f'{longest * 1e6:.3f} us is the longest dwell that stays inside one dk.',
            SeqCraftWarning, stacklevel=2,
        )

    # ------------------------------------------------------------------ what it knows
    @property
    def echo_spacing_s(self) -> float:
        """Seconds between consecutive acquired ``k = 0`` crossings.  Raises at one echo."""
        if len(self.te_s) < 2:
            msg = format_error(
                'echo_spacing_s is the gap between two echoes and this readout acquires one.',
                {'echoes': self.echoes, 'variant': self.variant},
                ['raise echoes', 'time_to_echo(0) is the one time this readout has'],
            )
            raise ConfigurationError(msg)
        return self.te_s[1] - self.te_s[0]

    def fov_mm_of(self, kr_norm: np.ndarray | float) -> np.ndarray:
        """
        Return the field of view supported at ``kr / kmax``, in millimetres.

        The density, evaluated -- and the reconstruction wants it, because the area one sample
        represents is the arc spacing times the *radial turn spacing*, and the turn spacing is
        exactly what `density` changed.  ``w ~ |dk/dt| / FOV(|k|)`` is the analytic density
        compensation that follows, and it is the one an approximation should be checked against.
        """
        return self.fov_mm * np.asarray(self._fov_of(np.asarray(kr_norm, dtype=float)), dtype=float)

    def echo_sample(self, echo: int = 0) -> int:
        """
        Which sample of :attr:`t_adc_s` carries ``k = 0`` for `echo`.

        Read off the measured trajectory rather than solved for.  ``k = 0`` is a *sample*, not an
        instant: solving for where ``k`` crosses zero and then finding the nearest sample is how
        the two come to disagree by half a dwell, and half a dwell of a spiral is a real position
        error in a real reconstruction.
        """
        return self._echo_samples[self._check_echo(echo)]

    def time_to_echo(self, echo: int = 0) -> float:
        """Seconds from the start of this module's block to `echo`'s ``k = 0`` sample."""
        return self.te_s[self._check_echo(echo)]

    def k_per_m(self, shot: int | None = None, *, angle_rad: float | None = None) -> np.ndarray:
        """
        Return ``(2, n)`` -- **where the samples are**, in 1/m, in acquisition order.

        Measured off the built events' own knots and integrated exactly, not read off the design
        pass, which is discarded and held on no attribute.  Rotating it for the interleaf is exact
        rather than approximate for the same reason the waveform is: the arm at angle ``theta`` is
        the angle-zero arm multiplied by ``exp(i theta)``, and every step between the two is linear.

        Against ``sc.kspace(tree, opts)['k_adc']`` -- which shares no code with this -- the two
        agree to 1.3-4.7e-5 1/m.  That residual is pypulseq's shape compression, not a design
        error, and it cannot be reduced: it round-trips a 1363-sample waveform with a relative
        error of 5e-8 and says so, twice per block.
        """
        _, angle = self._resolve_shot(shot, angle_rad)
        cos, sin = math.cos(angle), math.sin(angle)
        return np.vstack([cos * self._k0[0] - sin * self._k0[1],
                          sin * self._k0[0] + cos * self._k0[1]])

    def k_start_per_m(self, shot: int | None = None, *,
                      angle_rad: float | None = None) -> np.ndarray:
        """
        Return the ``(2,)`` k this readout expects to already be there when its block starts.

        Zero when the module prephases itself.  With ``prephase=False`` it is
        ``kmax * (cos theta, sin theta)``, and what the caller builds from it is a **separate
        event** -- so the vector slew limit becomes theirs, which is what
        :func:`~seqcraft.modules._support.oblique_trapezoid` exists for.  Two trapezoids designed
        independently and stretched to a common duration each reach ``max_slew`` and together ask
        for up to ``sqrt(2)`` times it -- measured at 250.3 T/m/s against a 180 T/m/s amplifier
        at 45 degrees, and 395.3 at 57.3 degrees, both compiling cleanly.
        """
        return self._rotate(self._k_start, shot, angle_rad)

    def k_end_per_m(self, shot: int | None = None, *,
                    angle_rad: float | None = None) -> np.ndarray:
        """Return the ``(2,)`` k a caller's rewinder has to cancel after this block."""
        return self._rotate(self._k_end, shot, angle_rad)

    # ----------------------------------------------------------------------- assembly
    def build(
        self,
        *,
        shot: int | None = None,
        angle_rad: float | None = None,
        acquire: bool = True,
        segment: int | None = None,
        reference: bool | None = None,
        phase_deg: float = 0.0,
    ) -> LogicBlock:
        """
        Return one interleaf: two arbitrary gradients, the ADCs the sample limit needs, the labels.

        Parameters
        ----------
        shot, angle_rad
            **Exactly one is required.**  ``shot=n`` is the uniform interleave at
            ``2 pi n / shots``; `angle_rad` is the escape hatch for an ordering the module should
            not know about -- a golden angle, a measured rotation, a re-ordered acquisition.
        acquire
            ``False`` plays the identical gradients with no ADCs and no labels, for a dummy shot.
            The duration is unchanged, which is the point: a dummy has to load the amplifier
            exactly as a real shot does, or the eddy-current state it is establishing is not the
            one that will be acquired.
        segment
            The shot's index, emitted as ``SEG`` on the first acquired ADC.  ``None`` emits none.
        reference
            Whether this is autocalibration rather than imaging data, emitted as the complementary
            ``REF``/``IMA`` pair.  ``None`` emits neither.
        phase_deg
            The receiver phase, and it is not optional in a spoiled sequence: the receiver is
            phase-locked to the transmitter, so every ADC of the readout an excitation opens has to
            demodulate at the carrier phase that excitation was given.

        Notes
        -----
        **There is no ``REV``.**  An earlier draft emitted it on every inward arm.  Nothing about a
        spiral has to be un-reversed -- :meth:`k_per_m` says where every sample is and
        :attr:`t_adc_s` says when, and a NUFFT needs nothing else -- so the only justification was
        a per-arm phase correction nobody has asked for.  An emitted label is a promise about the
        file, and this one had no reader.
        """
        index, angle = self._resolve_shot(shot, angle_rad)
        gx, gy = self._pair_at(angle)
        out = LogicBlock().add(0.0, gx).add(0.0, gy)
        if not acquire:
            return out

        phase_rad = float(np.deg2rad(phase_deg))
        per_region = self.adc_segments
        for event, (start, count) in enumerate(self._spans):
            adc = self._adc[count]
            if phase_rad:
                adc = derive(adc, phase_offset=phase_rad)
            # pypulseq's make_adc forces `delay >= adc_dead_time` when a system is passed, and
            # that delay is *inside* the event -- so a node placed at the wanted first-sample time
            # samples one dead time late.  Ten microseconds of it is 3.3 dk at this protocol, and
            # the trajectory is then simply wrong with nothing to look wrong.
            node = start - float(adc.delay)
            out.add(node, adc)
            if event == 0 and index is not None:
                out.add(node, pp.make_label(type='SET', label='LIN', value=int(index)))
            if event == 0:
                self._shot_labels(out, node, segment, reference)
            if self.echoes > 1 and event % per_region == 0:
                out.add(node, pp.make_label(type='SET', label='ECO',
                                            value=int(event // per_region)))
            if per_region > 1:
                out.add(node, pp.make_label(type='SET', label='SET',
                                            value=int(event % per_region)))
        return out

    def _shot_labels(self, out: LogicBlock, node: float, segment: int | None,
                     reference: bool | None) -> None:
        """Emit what describes the *shot* rather than the echo, on its first acquired ADC."""
        if segment is not None:
            out.add(node, pp.make_label(type='SET', label='SEG', value=int(segment)))
        if reference is not None:
            out.add(node, pp.make_label(type='SET', label='REF', value=int(bool(reference))))
            out.add(node, pp.make_label(type='SET', label='IMA', value=int(not reference)))

    def _pair_at(self, angle_rad: float) -> tuple[Event, Event]:
        """The rotated waveform as two pypulseq arbitrary gradients, ``first`` and ``last`` zero."""
        rotated = self._waveform * np.exp(1j * angle_rad)
        return tuple(
            derive(pp.make_arbitrary_grad(
                channel=axis, waveform=np.ascontiguousarray(component), first=0.0, last=0.0,
                max_grad=math.inf, max_slew=math.inf, system=self.opts))
            for axis, component in zip(self.axes, (rotated.real, rotated.imag))
        )

    # ------------------------------------------------------------------------ refusals
    def _resolve_dwell(self, dwell_s: float | None, bandwidth_hz: float | None) -> float:
        """Return the dwell, snapped **up** onto the ADC raster, from either spelling of it."""
        if (dwell_s is None) == (bandwidth_hz is None):
            msg = format_error(
                'exactly one of dwell_s and bandwidth_hz is required.',
                {'dwell_s': dwell_s, 'bandwidth_hz': bandwidth_hz},
                ['bandwidth_hz is the whole receive bandwidth 1/dwell, not a per-pixel figure -- '
                 'a spiral has no readout direction for a pixel to be measured along'],
            )
            raise ConfigurationError(msg)
        wanted = (require_positive(dwell_s, 'dwell_s') if dwell_s is not None
                  else 1.0 / require_positive(bandwidth_hz, 'bandwidth_hz'))
        # Up rather than down: rounding down raises the sample count, which is the direction that
        # turns a readout that fits the interpreter's sample limit into one that does not.
        return float(Raster(float(self.opts.adc_raster_time)).ceil(wanted))

    def _resolve_ends(self, prephase: bool | None, rewind: bool | None) -> tuple[bool, bool]:
        """
        Decide whether the module supplies each end ramp, refusing one that cannot take effect.

        ``None`` is "whatever this variant needs", which is what makes ``Spiral2D(variant='out')``
        legal without a caller thinking about it.  An explicit ``True`` where the arm already
        starts (or ends) at the origin raises rather than being ignored, because a caller who
        wrote it believes a ramp is being played and is placing an echo time against it.
        """
        needs_start, needs_end = self._directions[0] < 0, self._directions[-1] > 0
        for value, needed, name, where in ((prephase, needs_start, 'prephase', 'starts'),
                                           (rewind, needs_end, 'rewind', 'ends')):
            if value is True and not needed:
                msg = format_error(
                    f"variant={self.variant!r} {where} at k = 0, so there is nothing for "
                    f'{name}=True to do.',
                    {'variant': self.variant, name: value},
                    [f'drop {name}, which defaults to whatever the variant needs',
                     "variant='in' and 'in-out' are the ones that start away from the origin"],
                )
                raise ConfigurationError(msg)
        return (needs_start and prephase is not False, needs_end and rewind is not False)

    def _resolve_shot(self, shot: int | None,
                      angle_rad: float | None) -> tuple[int | None, float]:
        """Return ``(the interleaf index or None, the rotation in radians)``."""
        if (shot is None) == (angle_rad is None):
            msg = format_error(
                'exactly one of shot and angle_rad is required.',
                {'shot': shot, 'angle_rad': angle_rad, 'shots': self.shots},
                [f'shot=n is the uniform interleave at 2*pi*n/{self.shots}',
                 'angle_rad is for an ordering the module should not know about',
                 'there is no default, because acquiring one interleaf `shots` times is a '
                 'plausible bug with a plausible-looking image'],
            )
            raise ConfigurationError(msg)
        if angle_rad is not None:
            return None, float(angle_rad)
        index = int(shot)
        if not 0 <= index < self.shots:
            msg = format_error(
                f'shot must be in range({self.shots}), got {index}.',
                {'shot': index, 'shots': self.shots},
                ['build(angle_rad=...) is the way to rotate by anything else'],
            )
            raise ConfigurationError(msg)
        return index, 2.0 * math.pi * index / self.shots

    def _rotate(self, value: complex, shot: int | None, angle_rad: float | None) -> np.ndarray:
        """One complex k, rotated into the ``(2,)`` array a caller builds a gradient from."""
        _, angle = self._resolve_shot(shot, angle_rad)
        turned = value * np.exp(1j * angle)
        return np.array([turned.real, turned.imag], dtype=float)

    def _check_variant(self, variant: str) -> str:
        """Return `variant` having checked it names one of the four traversals."""
        if variant not in VARIANTS:
            listed = ', '.join(repr(name) for name in VARIANTS)
            msg = format_error(
                f'variant must be one of {listed}, got {variant!r}.',
                {'variant': variant},
                ["three of the four are one waveform played backwards and the fourth is both"],
            )
            raise ConfigurationError(msg)
        return variant

    def _check_axes(self, axes: tuple[str, str]) -> tuple[str, str]:
        """Return `axes` having checked they are two different logical channels."""
        first, second = (require_axis(axes[0], 'axes[0]'), require_axis(axes[1], 'axes[1]'))
        if first == second:
            msg = format_error(
                f'axes must be two different channels, got {axes!r}.', {'axes': axes},
                ['a spiral is two-dimensional; the pair is the plane it turns in'],
            )
            raise ConfigurationError(msg)
        return (first, second)

    def _check_echo(self, echo: int) -> int:
        """Return `echo` having checked this readout acquires it."""
        index = int(echo)
        if not 0 <= index < len(self.te_s):
            msg = format_error(
                f'echo must be in range({len(self.te_s)}), got {index}.',
                {'echo': index, 'echoes': self.echoes, 'variant': self.variant},
                ["variant='out-in' acquires k = 0 at each end, so it has echoes + 0 of them"],
            )
            raise ConfigurationError(msg)
        return index
