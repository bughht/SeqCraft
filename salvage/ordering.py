"""
Loop-ordering utilities.  **Not part of the seqcraft package** -- see ``salvage/README.md``.

Moved out whole rather than deleted: four of the six functions here had never had a caller, and
ordering tables are sequence-programming choices rather than physics, so they belong in the
opinionated library rather than the core (``docs/adr/003-scanner-and-module-reform.md``).

Loops stay in user code -- a plain ``for``, which is the familiar idiom and keeps the
control flow visible.  What this module supplies is the handful of pieces that were
copy-pasted (and quietly diverged) between notebooks: slice interleaving, view ordering,
golden-angle increments, and the RF-spoiling phase.

**Every derived quantity is a closed form of the loop counter, never an accumulator.**
That is a hard rule with a concrete payoff: the phase of excitation *n* can be computed
without replaying excitations 0..n-1, so a single TR can be rebuilt in isolation for
debugging, and two runs of the same loop cannot drift apart.  The reference implementation
kept ``self.rf_spoil_idx`` and ``self.rf_spoil_phase`` on the module, which made the same
builder produce different sequences on a second call, and tracked readout polarity in
``self.ROpolarity`` and ``self.nav_sign_track`` for the same reason.

Examples
--------
>>> interleaved_slice_order(7)
(0, 2, 4, 6, 1, 3, 5)
>>> centric_order(8)
(4, 3, 5, 2, 6, 1, 7, 0)
>>> round(rf_spoil_phase(3) * 180 / 3.141592653589793, 3)
342.0
"""

from __future__ import annotations

import math

__all__ = [
    'GOLDEN_ANGLE_RAD',
    'bit_reversed_order',
    'centric_order',
    'golden_angle',
    'interleaved_slice_order',
    'linear_order',
    'rf_spoil_phase',
]

#: The golden angle in radians, ``pi * (3 - sqrt(5))`` -- 137.507...  degrees.
GOLDEN_ANGLE_RAD = math.pi * (3.0 - math.sqrt(5.0))


def linear_order(n: int) -> tuple[int, ...]:
    """
    ``0, 1, ..., n-1``.

    Examples
    --------
    >>> linear_order(4)
    (0, 1, 2, 3)
    """
    return tuple(range(n))


def interleaved_slice_order(n: int, step: int = 2) -> tuple[int, ...]:
    """
    Slice order that maximises the gap between consecutively excited slices.

    Reduces slice cross-talk from imperfect slice profiles: with ``step=2`` all even
    slices are excited, then all odd ones, so neighbours are separated by a full pass
    rather than one TR.

    Parameters
    ----------
    n
        Number of slices.
    step
        Interleave factor.  ``2`` is the usual choice.

    Returns
    -------
    tuple of int
        A permutation of ``range(n)``; every slice appears exactly once.

    Examples
    --------
    >>> interleaved_slice_order(6)
    (0, 2, 4, 1, 3, 5)
    >>> interleaved_slice_order(7)
    (0, 2, 4, 6, 1, 3, 5)
    >>> sorted(interleaved_slice_order(11)) == list(range(11))
    True
    >>> interleaved_slice_order(6, step=3)
    (0, 3, 1, 4, 2, 5)
    """
    if n <= 0:
        return ()
    out: list[int] = []
    for offset in range(step):
        out.extend(range(offset, n, step))
    return tuple(out)


def centric_order(n: int) -> tuple[int, ...]:
    """
    View order starting at the k-space centre and alternating outwards.

    Puts the highest-signal lines first, which matters whenever the magnetisation is
    decaying through the acquisition (single-shot, or after a magnetisation preparation).

    Parameters
    ----------
    n
        Number of lines.

    Returns
    -------
    tuple of int
        A permutation of ``range(n)`` beginning at ``n // 2``.

    Examples
    --------
    >>> centric_order(8)
    (4, 3, 5, 2, 6, 1, 7, 0)
    >>> centric_order(5)
    (2, 1, 3, 0, 4)
    >>> sorted(centric_order(9)) == list(range(9))
    True
    """
    if n <= 0:
        return ()
    centre = n // 2
    out = [centre]
    for d in range(1, n):
        for candidate in (centre - d, centre + d):
            if 0 <= candidate < n:
                out.append(candidate)
        if len(out) >= n:
            break
    return tuple(out[:n])


def bit_reversed_order(n: int) -> tuple[int, ...]:
    """
    Bit-reversed permutation of ``range(n)``, for `n` a power of two.

    Spreads acquisition order maximally in k-space, which decorrelates motion and
    system-drift artefacts from the phase-encode direction.

    Raises
    ------
    ValueError
        If `n` is not a positive power of two.

    Examples
    --------
    >>> bit_reversed_order(8)
    (0, 4, 2, 6, 1, 5, 3, 7)
    >>> sorted(bit_reversed_order(16)) == list(range(16))
    True
    """
    if n <= 0 or n & (n - 1):
        msg = f'bit_reversed_order needs a positive power of two, got {n}.'
        raise ValueError(msg)
    bits = n.bit_length() - 1
    return tuple(int(format(i, f'0{bits}b')[::-1], 2) for i in range(n))


def golden_angle(index: int, *, increment_rad: float = GOLDEN_ANGLE_RAD) -> float:
    """
    Rotation angle of spoke `index`, wrapped into ``[0, 2*pi)``.

    Closed form, so spoke 4711 can be produced without generating the previous 4710 --
    which is what makes a single TR reproducible in isolation.

    Parameters
    ----------
    index
        Spoke number, from zero.
    increment_rad
        Angular increment.  Defaults to the golden angle; use ``pi / n_spokes`` for a
        uniform radial acquisition.

    Returns
    -------
    float
        Angle in radians.

    Examples
    --------
    >>> round(golden_angle(0), 6)
    0.0
    >>> round(golden_angle(1) * 180 / math.pi, 3)
    137.508
    """
    return (index * increment_rad) % (2.0 * math.pi)


def rf_spoil_phase(index: int, increment_rad: float = math.pi * 117.0 / 180.0) -> float:
    """
    Quadratic RF-spoiling phase for excitation `index`, wrapped into ``[0, 2*pi)``.

    Uses the standard quadratic schedule ``phi_n = increment * n * (n + 1) / 2``, evaluated
    in closed form rather than accumulated.  117 degrees is Zur's widely used increment.

    Parameters
    ----------
    index
        Excitation number, from zero.  Count *excitations*, not lines: with multislice
        acquisition this is the global RF counter.
    increment_rad
        Phase increment.  Default 117 degrees.

    Returns
    -------
    float
        Phase in radians, to add to the RF carrier phase.

    Examples
    --------
    >>> round(rf_spoil_phase(0), 6)
    0.0
    >>> round(rf_spoil_phase(1) * 180 / math.pi, 3)
    117.0
    >>> round(rf_spoil_phase(2) * 180 / math.pi, 3)
    351.0
    >>> round(rf_spoil_phase(3) * 180 / math.pi, 3)   # 117*3*4/2 = 702, mod 360
    342.0
    """
    return (increment_rad * index * (index + 1) / 2.0) % (2.0 * math.pi)
