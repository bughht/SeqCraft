"""
The extracted ``RadialReadout`` against the invariant table, the metamorphic test and a sweep.

Three things, in the order that makes each of them mean something:

1. **the invariants the references were measured against** (R1-R7), re-run on the candidate;
2. **the rotation-equivariance test, unchanged** -- it was written and passed against the
   official reference before this module existed, which is what stops it encoding the module's
   own behaviour;
3. **a sweep**, because one protocol proves nothing about a reusable abstraction.

    python tools/module_mining/candidates/radial/run_validation.py
"""

from __future__ import annotations

import json
import logging
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from module_mining import fingerprint  # noqa: E402
from module_mining.candidates.radial import checks  # noqa: E402
from module_mining.candidates.radial.adapters import seqcraft_radial  # noqa: E402

ANGLES = [math.pi * n / 8 for n in range(8)]

#: Chosen so that each case can fail on its own account: the reference protocol, a second
#: matrix, a second field of view, two dwell times, the full spoke and the centre-out extreme.
CASES = [
    ('reference protocol', {}),
    ('matrix 128', {'matrix': 128}),
    ('fov 200 mm', {'fov_mm': 200.0}),
    ('dwell 10 us', {'dwell_s': 10e-6}),
    ('dwell 40 us', {'dwell_s': 40e-6}),
    ('partial 0.75', {'partial_fourier': 0.75}),
    ('centre-out 0.5', {'partial_fourier': 0.5}),
    ('matrix 96', {'matrix': 96}),
]


def _spokes(overrides: dict) -> tuple[list[dict], list[np.ndarray], dict]:
    """One compiled spoke per angle -- see the adapter for why they are not stacked."""
    rows, trajectories, claims = [], [], {}
    for index, angle in enumerate(ANGLES):
        reference = seqcraft_radial.build(angle_rad=angle, **overrides)
        measured = fingerprint.measure(reference, centre='k-min')
        k_adc = np.asarray(reference.sequence.calculate_kspacePP()[0])
        row = checks.spoke_geometry(measured, k_adc)[0]
        row['index'] = index
        rows.append(row)
        trajectories.append(k_adc)
        claims = reference.semantic
    return rows, trajectories, claims


def main() -> None:
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    results: dict[str, object] = {}

    logging.info('%-20s %-9s %-7s %-10s %-11s %-10s %-9s %s', 'case', 'samples', 'centre',
                 'dk 1/m', 'claimed k0', 'measured', 'rotation', 'invariants')
    for name, overrides in CASES:
        rows, trajectories, claimed = _spokes(overrides)
        wanted = [math.degrees(a) for a in ANGLES]
        verdicts = checks.check_spokes(rows, expected_angles_deg=wanted)

        equivariance = [checks.rotation_equivariance(trajectories[0], trajectories[n],
                                                     math.degrees(ANGLES[n]))
                        for n in range(1, len(ANGLES))]
        verdicts.append(max(equivariance, key=lambda row: row['worst']))

        centre_sample = claimed['center_sample']
        spoke = trajectories[0]
        measured_centre = float(np.linalg.norm(spoke[:, centre_sample]))

        results[name] = {'overrides': overrides, 'claims': claimed,
                         'invariants': {v['invariant']: {'worst': v['worst'], 'pass': v['pass']}
                                        for v in verdicts},
                         'all_pass': all(v['pass'] for v in verdicts)}
        logging.info('%-20s %-9d %-7d %-10.4f %-11.2f %-10.2f %-9.2e %s',
                     name, claimed['num_samples'], centre_sample, claimed['dk_per_m'],
                     claimed['k_first_per_m'], float(spoke[0, 0]),
                     verdicts[-1]['worst'],
                     'all pass' if results[name]['all_pass']
                     else [v['invariant'] for v in verdicts if not v['pass']])
        # The claim the module makes about its own centre, against the trajectory it produced.
        assert measured_centre < 1e-3, (name, measured_centre)

    failed = [name for name, row in results.items() if not row['all_pass']]
    logging.info('\n%d/%d cases pass every radial invariant%s',
                 len(results) - len(failed), len(results),
                 '' if not failed else f'; failed: {failed}')

    out = Path(__file__).with_name('validation_result.json')
    out.write_text(json.dumps(results, indent=1, default=float))
    logging.info('-> %s', out)


if __name__ == '__main__':
    main()
