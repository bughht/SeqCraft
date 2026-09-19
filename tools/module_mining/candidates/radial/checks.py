"""
Radial trajectory measurements -- the L4 layer for a non-Cartesian readout.

The Cartesian checks in :mod:`module_mining.fingerprint` ask where ``kx`` returns to zero and
whether areas balance about a refocusing pulse.  A spoke has neither question: what it has is an
**orientation**, an **extent**, a **sample spacing** and a **centre**, and the one property that
no single spoke can show is that the family of them is a rotation of one design.

Everything here is measured from the compiled trajectory, so it applies equally to a reference
script and to a candidate module, and neither has to say what it intended.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Sequence


def spoke_geometry(measured: dict[str, Any], k_adc: np.ndarray) -> list[dict[str, Any]]:
    """
    Per readout: its angle, its extent, its sample spacing and where its centre sample sits.

    The angle is taken from the **fitted direction** of the sampled points rather than from the
    first and last of them, so that a spoke which does not start at ``-kmax`` -- a centre-out one
    -- reports the same quantity as one that does.
    """
    rows, start = [], 0
    for readout in measured['readouts']:
        n = readout['n_samples']
        k = k_adc[:, start:start + n]
        start += n

        # Direction as the principal axis of the sampled points, sign-fixed by the sweep.
        centred = k[:2] - k[:2].mean(axis=1, keepdims=True)
        direction = np.linalg.svd(centred, full_matrices=False)[0][:, 0]
        if np.dot(k[:2, -1] - k[:2, 0], direction) < 0:
            direction = -direction

        radius = np.linalg.norm(k[:2], axis=0)
        along = direction @ k[:2]
        centre = int(np.argmin(radius))
        steps = np.diff(along)
        rows.append({
            'index': readout['index'],
            'angle_deg': float(np.degrees(np.arctan2(direction[1], direction[0]))),
            'centre_sample': centre,
            'k_at_centre_per_m': float(radius[centre]),
            'k_min_per_m': float(along.min()),
            'k_max_per_m': float(along.max()),
            'dk_per_m': float(np.mean(steps)),
            'dk_spread_per_m': float(np.ptp(steps)),
            # How far the samples stray from the straight line through them: a spoke that is not
            # straight is not a spoke, and nothing else here would notice.
            'straightness_per_m': float(
                np.abs(direction[0] * k[1] - direction[1] * k[0]).max()),
            'kz_max_per_m': float(np.abs(k[2]).max()),
        })
    return rows


def check_spokes(rows: Sequence[dict[str, Any]], *, expected_angles_deg: Sequence[float] | None = None,
                 k_tol_per_m: float = 1e-3, angle_tol_deg: float = 1e-6) -> list[dict[str, Any]]:
    """The radial invariant table, measured: centre, straightness, spacing, extent and angle."""
    def worst(values):
        return max(values, key=lambda pair: pair[0]) if values else (0.0, {})

    checks = [
        ('R1', 'a sample lands on k = 0',
         [(r['k_at_centre_per_m'], {'readout': r['index'], 'sample': r['centre_sample']})
          for r in rows], k_tol_per_m, '1/m'),
        ('R2', 'the samples lie on a straight line through the origin',
         [(r['straightness_per_m'], {'readout': r['index']}) for r in rows], k_tol_per_m, '1/m'),
        ('R3', 'the sample spacing is uniform along the spoke',
         [(r['dk_spread_per_m'], {'readout': r['index'], 'dk': r['dk_per_m']}) for r in rows],
         k_tol_per_m, '1/m'),
        ('R6', 'the spoke stays in plane',
         [(r['kz_max_per_m'], {'readout': r['index']}) for r in rows], k_tol_per_m, '1/m'),
    ]
    out = []
    for name, what, values, tol, unit in checks:
        value, where = worst(values)
        out.append({'invariant': name, 'what': what, 'worst': value, 'tolerance': tol,
                    'unit': unit, 'where': where, 'pass': value <= tol})

    extents = [r['k_max_per_m'] - r['k_min_per_m'] for r in rows]
    spacing = [r['dk_per_m'] for r in rows]
    out.append({'invariant': 'R4', 'what': 'every spoke has the same extent and spacing',
                'worst': float(max(np.ptp(extents), np.ptp(spacing))), 'tolerance': k_tol_per_m,
                'unit': '1/m', 'where': {'extent_per_m': extents[0], 'dk_per_m': spacing[0]},
                'pass': max(np.ptp(extents), np.ptp(spacing)) <= k_tol_per_m})

    if expected_angles_deg is not None:
        error = [(abs((r['angle_deg'] - want + 90) % 180 - 90), {'readout': r['index'],
                                                                 'want_deg': want,
                                                                 'got_deg': r['angle_deg']})
                 for r, want in zip(rows, expected_angles_deg)]
        value, where = worst(error)
        out.append({'invariant': 'R5', 'what': 'the spoke points where it was asked to',
                    'worst': value, 'tolerance': angle_tol_deg, 'unit': 'deg', 'where': where,
                    'pass': value <= angle_tol_deg})
    return out


def rotation_equivariance(k_zero: np.ndarray, k_phi: np.ndarray, angle_deg: float) -> dict[str, Any]:
    """
    The metamorphic test: is the spoke at ``phi`` the spoke at 0, rotated?

    Compares **k-space sample positions** rather than gradient waveforms, because a rotation is a
    statement about the trajectory and an implementation is free to reach it by rotating either.
    A design whose spokes are not rotations of one another is not a radial readout, however
    plausible each spoke looks on its own.
    """
    phi = np.deg2rad(angle_deg)
    rotation = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
    rotated = rotation @ k_zero[:2]
    error = float(np.abs(rotated - k_phi[:2]).max())
    return {
        'invariant': 'R7',
        'what': f'the spoke at {angle_deg:g} deg is the spoke at 0 deg, rotated',
        'worst': error, 'tolerance': 1e-3, 'unit': '1/m',
        'where': {'samples': int(k_zero.shape[1])},
        'pass': error <= 1e-3,
    }
