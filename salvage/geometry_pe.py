"""
The phase-encode table.  **Not part of the seqcraft package** -- see ``salvage/README.md``.

Lifted out of ``seqcraft.core.geometry.Geometry`` rather than deleted, because the sharing it
was written for is currently aspirational.  ``docs/architecture.md`` defends ``geometry`` as
holding *"one authoritative phase-encode index computation shared by ``kspace_center_line`` and
the ``LIN`` label values, so the two cannot disagree"* -- but nothing in the repository consumes
the table, because the module library that would have has been deleted
(``docs/adr/003-scanner-and-module-reform.md``).  ``kspace_center_line`` stays on ``Geometry``,
since ``definitions()`` writes it; the table that was meant to agree with it is here.

What is worth keeping is one line, and it is not obvious:

    skip += (c - skip) % r

Partial Fourier truncates the low end of k-space, and acceleration then skips forward from
whatever index that leaves.  Without the residue nudge, **an even skip and an odd centre miss
k = 0 entirely** -- the sequence acquires a full table, reconstructs, and produces an image with
no DC term, which looks like a windowing problem rather than a sampling bug.  Nudging the start
*up* by ``(c - skip) mod r`` puts the centre line back on the sampled lattice at the cost of at
most ``r - 1`` lines at the truncated end, which is the end that was being thrown away anyway.

Index conventions
-----------------
Phase-encode lines are numbered on the **reconstruction grid**: ``0 .. matrix_pe - 1``, with the
k-space centre at ``matrix_pe // 2`` (the Siemens convention).  That is the index space of the
``LIN`` label and of ``kSpaceCenterLine``.

These are plain functions taking numbers, not properties on a dataclass, so a readout module can
call them with its own matrix without constructing a ``Geometry`` first.  When the opinionated
module library exists, they move into it.

Examples
--------
A 146 matrix at partial Fourier 0.75 with twofold acceleration -- the reference
implementation's own numbers:

>>> center_index(146)
73
>>> first_index(146, partial_fourier=0.75, accel=2)
37
>>> lines = table(146, partial_fourier=0.75, accel=2)
>>> len(lines)
55
>>> center_index(146) in lines                  # the centre line is always sampled
True

The nudge, isolated.  An even skip with an odd centre would land on the wrong parity:

>>> first_index(146, partial_fourier=0.75, accel=2) % 2 == center_index(146) % 2
True

Interleaved shots take every ``n_shots``-th entry of the table:

>>> table(8, n_shots=2, shot=0), table(8, n_shots=2, shot=1)
((0, 2, 4, 6), (1, 3, 5, 7))
"""

from __future__ import annotations

import math

__all__ = ['center_index', 'first_index', 'round_half_up', 'table']


def round_half_up(x: float) -> int:
    """
    Round half away from zero, unlike Python's banker's rounding.

    ``round(109.5)`` is 110 in Python only by luck of the float representation; ``round(0.5)``
    is 0.  A partial-Fourier line count must not depend on that.

    Examples
    --------
    >>> round_half_up(0.5), round_half_up(1.5), round_half_up(109.5)
    (1, 2, 110)
    >>> round_half_up(-0.5)
    -1
    """
    return int(math.floor(x + 0.5)) if x >= 0 else -int(math.floor(-x + 0.5))


def center_index(matrix: int) -> int:
    """
    Recon-grid index of the k-space centre: ``matrix // 2``.

    The same expression ``Geometry.kspace_center_line`` uses, which is the point: this is what
    the ``LIN`` label and the ``kSpaceCenterLine`` definition must agree on.
    """
    return matrix // 2


def first_index(matrix: int, *, partial_fourier: float = 1.0, accel: int = 1) -> int:
    """
    Recon-grid index of the first acquired line.

    Parameters
    ----------
    matrix
        Number of lines on the reconstruction grid.
    partial_fourier
        Fraction of k-space sampled, in ``(0.5, 1.0]``.  Truncates the **low** end.
    accel
        Undersampling (skip) factor.

    Notes
    -----
    The start is nudged **up** by the residue ``(centre - skip) mod accel`` so the centre line is
    always sampled; see the module docstring for why that is not optional.
    """
    skip = matrix - round_half_up(partial_fourier * matrix)
    skip += (center_index(matrix) - skip) % accel
    return skip


def table(
    matrix: int,
    *,
    partial_fourier: float = 1.0,
    accel: int = 1,
    n_shots: int = 1,
    shot: int = 0,
) -> tuple[int, ...]:
    """
    Every acquired phase-encode line, in ascending recon-grid order.

    Parameters
    ----------
    matrix, partial_fourier, accel
        As for :func:`first_index`.
    n_shots
        Number of interleaved shots the table is split across.
    shot
        Which shot to return.  Shots interleave, so shot 0 takes entries
        ``0, n_shots, 2 * n_shots, ...`` of the full table.

    Raises
    ------
    ValueError
        If `shot` is outside ``[0, n_shots - 1]``, or if the parameters leave no line to acquire.
    """
    if not 0 <= shot < n_shots:
        msg = f'shot must be in [0, {n_shots - 1}], got {shot}.'
        raise ValueError(msg)
    lines = tuple(range(first_index(matrix, partial_fourier=partial_fourier, accel=accel),
                        matrix, accel))
    if not lines:
        msg = (
            f'the phase-encode table is empty for matrix={matrix}, '
            f'partial_fourier={partial_fourier}, accel={accel}.'
        )
        raise ValueError(msg)
    return lines[shot::n_shots]
