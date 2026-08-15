"""
Exact b-value of a trapezoidal diffusion encoding, and the two solvers around it.

Lifted from the deleted ``seqcraft.modules.encoding.diffusion`` as plain functions: they take
numbers in SI and hertz and return numbers, with no scanner object and no module in sight.

Everything here is in **s/m^2**, not the s/mm^2 a protocol is written in.  The conversion is
``1 s/mm^2 == 1e6 s/m^2``.

Why this was worth keeping
--------------------------
The commonly quoted ramp correction to ``b = (2 pi G)^2 delta^2 (Delta - delta/3)`` --
``+- eps^3/30 -+ delta eps^2/6`` -- is written with several incompatible conventions for whether
``delta`` includes the ramps.  Using one with the wrong convention overstates b by 2-4 % at the
short lobe durations diffusion sequences actually use, and a b-value wrong by 3 % biases every
reported diffusivity by 3 %.  :func:`b_of_monopolar` integrates ``|k(t)|^2`` directly instead, so
there is no convention to pick.

Examples
--------
The ramp-free limit is the textbook expression, which is the check worth making:

>>> import math
>>> G, delta, Delta = 1.0e6, 20e-3, 30e-3          # Hz/m, s, s
>>> exact = b_of_monopolar(G, delta, Delta, ramp=1e-9)
>>> textbook = (2 * math.pi * G) ** 2 * delta ** 2 * (Delta - delta / 3)
>>> abs(exact / textbook - 1) < 1e-5
True
"""

from __future__ import annotations

import math

__all__ = [
    'amplitude_for_b',
    'b_of_bipolar',
    'b_of_monopolar',
    'shortest_lobe_for_b',
]


def b_of_monopolar(amplitude: float, delta: float, big_delta: float, *, ramp: float) -> float:
    """
    Return the b-value of a trapezoidal **monopolar** pair, in s/m^2, exactly.

    Parameters
    ----------
    amplitude
        Gradient amplitude, Hz/m.
    delta
        Lobe duration including both ramps, seconds.
    big_delta
        Separation between the two lobes' starts, seconds.
    ramp
        Ramp time of one edge, seconds.

    Returns
    -------
    float
        b in s/m^2.

    Notes
    -----
    Derived by integrating ``b = (2 pi)^2 * integral |k(t)|^2 dt`` piecewise over the ramp, the
    flat top and the fall of each lobe, with the refocusing pulse's sign flip applied.  With
    ``A = G (delta - eps)`` the area of one lobe::

        I = G^2 [ eps^3/10 + ((delta - 1.5 eps)^3 - eps^3/8)/3
                  + (delta - eps)^2 eps - (delta - eps) eps^2/3 ]
        b = (2 pi)^2 [ 2 I + A^2 (Delta - delta) ]

    As ``eps -> 0`` this collapses to ``(2 pi G)^2 delta^2 (Delta - delta/3)``.
    """
    eps = ramp
    area = amplitude * (delta - eps)
    integral = amplitude**2 * (
        eps**3 / 10.0
        + ((delta - 1.5 * eps) ** 3 - eps**3 / 8.0) / 3.0
        + (delta - eps) ** 2 * eps
        - (delta - eps) * eps**2 / 3.0
    )
    return (2.0 * math.pi) ** 2 * (2.0 * integral + area**2 * (big_delta - delta))


def b_of_bipolar(amplitude: float, delta: float, big_delta: float, *, ramp: float) -> float:
    """
    Return the b-value of a **bipolar** (flow-compensated) pair, in s/m^2.

    Each lobe is two opposing half-length sub-lobes, so the pair behaves like a monopolar encoding
    of half the lobe duration whose accumulated phase reverses mid-lobe.  Integrating ``|k|^2`` over
    the four sub-lobes gives the expression below.

    `big_delta` is accepted for signature compatibility with :func:`b_of_monopolar` and does not
    enter: the net area of a bipolar lobe pair is zero, so the gap between them contributes
    nothing and all the weighting comes from the sub-lobes.  The factor of a quarter relative to a
    monopolar encoding of the same duration is the cost of that cancellation -- which is what buys
    insensitivity to bulk motion.
    """
    del big_delta
    half = delta / 2.0
    eps = ramp
    integral = amplitude**2 * (
        eps**3 / 10.0
        + ((half - 1.5 * eps) ** 3 - eps**3 / 8.0) / 3.0
        + (half - eps) ** 2 * eps
        - (half - eps) * eps**2 / 3.0
    )
    return (2.0 * math.pi) ** 2 * 4.0 * integral


def shortest_lobe_for_b(
    target_b: float,
    *,
    max_amplitude: float,
    ramp: float,
    gap: float,
    grad_raster: float,
    bipolar: bool = False,
    max_steps: int = 200_000,
) -> float:
    """
    Return the shortest lobe duration reaching `target_b` at full amplitude, in seconds.

    Parameters
    ----------
    target_b
        Wanted b-value, s/m^2.
    max_amplitude
        The amplitude to design at, Hz/m.
    ramp
        Ramp time of one edge, seconds.
    gap
        Time between the end of the first lobe and the start of the second -- in a spin echo, the
        refocusing pulse's duration.  ``Delta = delta + gap``.
    grad_raster
        Gradient raster, seconds.  The search steps on it, so the answer lands on it.
    bipolar
        Use :func:`b_of_bipolar` instead of :func:`b_of_monopolar`.

    Returns
    -------
    float
        The lobe duration, seconds.

    Raises
    ------
    ValueError
        The b-value is unreachable within `max_steps` rasters.

    Notes
    -----
    A closed form does not exist because ``Delta`` depends on ``delta``.  This steps up the raster
    from the minimum -- at most a few hundred iterations in practice, each cheap -- rather than
    running a Newton solve on a residual that is not monotonic everywhere.
    """
    b_of = b_of_bipolar if bipolar else b_of_monopolar
    shortest = math.ceil((2.0 * ramp + grad_raster) / grad_raster) * grad_raster
    if target_b <= 0.0:
        return shortest

    delta = shortest
    for _ in range(max_steps):
        if b_of(max_amplitude, delta, delta + gap, ramp=ramp) >= target_b:
            return delta
        delta += grad_raster

    msg = (
        f'cannot reach b = {target_b:g} s/m^2 at {max_amplitude:g} Hz/m '
        f'within {max_steps} raster steps'
    )
    raise ValueError(msg)


def amplitude_for_b(
    target_b: float,
    *,
    delta: float,
    big_delta: float,
    ramp: float,
    bipolar: bool = False,
) -> float:
    """
    Return the amplitude reaching `target_b` in a lobe of duration `delta`, in Hz/m.

    ``b`` scales as amplitude squared, so this is a square root rather than a search.  It does not
    check the result against an amplifier limit -- compare it with ``opts.max_grad`` at the call
    site, where the message can name the sequence parameter to change.

    Examples
    --------
    >>> g = amplitude_for_b(1e9, delta=20e-3, big_delta=30e-3, ramp=200e-6)
    >>> abs(b_of_monopolar(g, 20e-3, 30e-3, ramp=200e-6) / 1e9 - 1) < 1e-12
    True
    """
    if target_b <= 0.0:
        return 0.0
    b_of = b_of_bipolar if bipolar else b_of_monopolar
    unit_b = b_of(1.0, delta, big_delta, ramp=ramp)
    return math.sqrt(target_b / unit_b)
