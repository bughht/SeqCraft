"""
DTI direction tables, and the number that says whether one is any good.

Lifted from the deleted ``seqcraft.modules.encoding.diffusion``.  Depends on numpy and nothing else
-- the seqcraft error types have been replaced by ``ValueError`` so this file stands alone.

Why this was worth keeping
--------------------------
The relaxation is the part that is easy to skip and expensive to skip.  A raw golden-angle spiral is
uniform in *area*, which is not the same as well separated: for 30 directions it leaves a pair only
14 degrees apart, and near-duplicate directions make the tensor fit ill-conditioned **without making
it look wrong**.  The scan completes, the fit converges, and the fractional anisotropy is biased.

:func:`direction_condition_number` is the number to compare schemes by -- not the raw count, and not
the minimum pairwise angle.

Examples
--------
>>> dirs = dti_directions(6)
>>> len(dirs)
6
>>> bool(direction_condition_number(dirs) < 2.5)
True
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Sequence

__all__ = ['direction_condition_number', 'dti_directions', 'normalise']


def dti_directions(n: int, *, iterations: int = 400) -> tuple[tuple[float, float, float], ...]:
    """
    Return `n` well-separated unit vectors for a DTI acquisition.

    Starts from a golden-angle spiral on the hemisphere and then relaxes it under **antipodal**
    electrostatic repulsion -- each direction is repelled by every other one *and by its
    reflection through the origin*.  That is the right symmetry: diffusion weighting depends on
    the square of the direction, so ``+d`` and ``-d`` are the same measurement, and a pair that
    is nearly antiparallel is just as redundant as a pair that is nearly parallel.

    Parameters
    ----------
    n
        Number of directions.  Six is the arithmetic minimum -- a diffusion tensor has six
        independent components -- 30 is a common clinical compromise, 60 or more for
        higher-order models.
    iterations
        Relaxation steps.  A few hundred is enough to converge for any practical `n`, and it
        takes milliseconds.

    Returns
    -------
    tuple of (x, y, z)
        Unit vectors, each on the upper hemisphere by convention.

    Examples
    --------
    >>> dirs = dti_directions(6)
    >>> all(abs(sum(c * c for c in d) - 1.0) < 1e-9 for d in dirs)
    True
    >>> all(d[2] >= -1e-12 for d in dirs)                     # upper hemisphere
    True
    >>> worst = max(abs(sum(a * b for a, b in zip(p, q)))
    ...             for i, p in enumerate(dirs) for q in dirs[i + 1:])
    >>> bool(worst < 0.55)                                    # no near-duplicate pair
    True
    """
    if n < 1:
        msg = f'dti_directions needs at least one direction, got {n}'
        raise ValueError(msg)
    if n == 1:
        return ((0.0, 0.0, 1.0),)

    golden = math.pi * (3.0 - math.sqrt(5.0))
    index = np.arange(n)
    z = 1.0 - (index + 0.5) / n
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    phi = index * golden
    points = np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)

    step = 0.5 / n
    for _ in range(iterations):
        # Repulsion from every other point and from its antipode, which is what makes the metric
        # the angle between *lines* rather than between vectors.
        force = np.zeros_like(points)
        for sign in (1.0, -1.0):
            delta = points[:, None, :] - sign * points[None, :, :]
            distance = np.linalg.norm(delta, axis=2)
            np.fill_diagonal(distance, np.inf)
            force += np.sum(delta / distance[:, :, None] ** 3, axis=1)
        points += step * force
        points /= np.linalg.norm(points, axis=1, keepdims=True)

    # Fold onto the upper hemisphere: a direction and its antipode are the same measurement, so
    # reporting them consistently makes the table easier to read and to compare against a vendor's.
    points[points[:, 2] < 0] *= -1.0
    return tuple((float(x), float(y), float(zz)) for x, y, zz in points)


def direction_condition_number(directions: Sequence[tuple[float, float, float]]) -> float:
    """
    Return the condition number of the diffusion-tensor design matrix for `directions`.

    Each direction contributes one row ``[dx^2, dy^2, dz^2, 2 dx dy, 2 dx dz, 2 dy dz]``, which is
    the linear system a tensor fit solves.  Its condition number says how much measurement noise
    the fit amplifies: the theoretical minimum for a well-spread set is around 1.6, six directions
    arranged along the icosahedral edges reach about 1.6, and anything above ~2.5 means some
    tensor components are much noisier than others.

    Examples
    --------
    >>> bool(direction_condition_number(dti_directions(6)) < 2.5)
    True
    >>> bool(direction_condition_number(dti_directions(30)) < 2.0)
    True
    """
    d = np.asarray(directions, dtype=float)
    design = np.stack(
        [
            d[:, 0] ** 2,
            d[:, 1] ** 2,
            d[:, 2] ** 2,
            2.0 * d[:, 0] * d[:, 1],
            2.0 * d[:, 0] * d[:, 2],
            2.0 * d[:, 1] * d[:, 2],
        ],
        axis=1,
    )
    return float(np.linalg.cond(design))


def normalise(direction: tuple[float, float, float]) -> tuple[float, float, float]:
    """
    Return `direction` scaled to unit length.

    Normalising is not a convenience -- it is a correctness requirement.  An unnormalised
    ``[1, 1, 1]`` asks for ``sqrt(3)`` times the intended vector amplitude, which both overstates
    the b-value and can exceed the amplifier on the vector norm while every individual axis looks
    legal.

    Examples
    --------
    >>> [round(c, 6) for c in normalise((1.0, 1.0, 0.0))]
    [0.707107, 0.707107, 0.0]
    """
    norm = math.sqrt(sum(c * c for c in direction))
    if norm == 0.0:
        msg = 'diffusion direction is the zero vector; use b = 0 for an unweighted volume'
        raise ValueError(msg)
    return (direction[0] / norm, direction[1] / norm, direction[2] / norm)
