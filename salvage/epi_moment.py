"""
The unit-amplitude trapezoid moment integral, and its closed-form inverse.

Lifted from the deleted ``seqcraft.modules.readout.epi``.  Two functions, no dependencies.

Why this was worth keeping
--------------------------
Amplitude factors out of every question a ramp-sampled readout asks -- where the prephaser has to
leave k, where k crosses zero, what the sampled extent is -- so the shape is integrated **once at
unit amplitude** and scaled.  That single source is what keeps a prephaser and a ``time_to_echo``
from disagreeing; computing them separately is how a readout ends up prephased to one place and
claiming its echo is at another, which shifts the whole image and looks like a shim problem.

The inverse is solved in closed form on whichever of the three segments contains the value, not by
search, so the echo time is exact rather than converged.

Examples
--------
>>> unit_moment(2e-4, 1e-4, 0.0) == 1e-4            # the whole triangle: ramp * 1
True
>>> t = unit_moment_inverse(0.8e-4, 1e-4, 0.0)
>>> abs(unit_moment(t, 1e-4, 0.0) - 0.8e-4) < 1e-18
True
"""

from __future__ import annotations

import math

__all__ = ['unit_moment', 'unit_moment_inverse']


def unit_moment(t: float, ramp: float, flat: float) -> float:
    """
    Area under a **unit-amplitude** trapezoid from its own start to `t`.

    Parameters
    ----------
    t
        Time from the trapezoid's start, seconds.  Clamped: before the start the area is zero,
        after the end it is the whole area.
    ramp
        Ramp time of one edge, seconds.
    flat
        Flat-top duration, seconds.

    Examples
    --------
    >>> unit_moment(0.0, 1e-4, 0.0)
    0.0
    >>> unit_moment(1e-4, 1e-4, 0.0) == 0.5e-4          # half the rising ramp
    True
    >>> unit_moment(3e-4, 1e-4, 1e-4) == 2e-4           # with a flat top: ramp + flat
    True
    """
    if t <= 0.0:
        return 0.0
    if t < ramp:
        return 0.5 * t * t / ramp
    if t <= ramp + flat:
        return 0.5 * ramp + (t - ramp)
    past = min(t - ramp - flat, ramp)
    return 0.5 * ramp + flat + past * (1.0 - past / (2.0 * ramp))


def unit_moment_inverse(value: float, ramp: float, flat: float) -> float:
    """
    Return the time at which :func:`unit_moment` reaches `value`.

    Monotonic, so the inverse is well defined.

    Examples
    --------
    >>> round(unit_moment_inverse(0.5e-4, 1e-4, 0.0) * 1e6, 9)
    100.0
    """
    if value <= 0.0:
        return 0.0
    if value < 0.5 * ramp:
        return math.sqrt(2.0 * value * ramp)
    if value <= 0.5 * ramp + flat:
        return ramp + (value - 0.5 * ramp)
    rest = value - 0.5 * ramp - flat
    # On the falling ramp, m = ramp/2 + flat + p - p^2/(2 ramp); the root inside [0, ramp].
    disc = max(0.0, ramp * ramp - 2.0 * ramp * rest)
    return ramp + flat + (ramp - math.sqrt(disc))
