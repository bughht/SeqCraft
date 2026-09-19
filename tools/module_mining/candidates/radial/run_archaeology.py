"""
The radial fine scan's first measurement pass: what the references actually build.

Three questions this answers with numbers rather than with reading, all of them inputs to the
candidate boundary:

* **where the centre of k-space lands relative to the ADC samples** -- three references derive it
  three different ways, so agreement between them is evidence and disagreement is a finding;
* **whether full-spoke and centre-out are one design or two** -- ``write_ute.py`` has a
  ``readout_asymmetry`` argument, and whether it spans both cleanly decides the question;
* **whether the spokes are rotations of one design**, which is the metamorphic test a radial
  readout has and a Cartesian one does not.

    python tools/module_mining/candidates/radial/run_archaeology.py
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from module_mining import fingerprint  # noqa: E402
from module_mining.candidates.radial import checks  # noqa: E402
from module_mining.candidates.radial.adapters import pypulseq_radial  # noqa: E402

SPOKES = 8
SAMPLES = 64


def _measure(reference):
    measured = fingerprint.measure(reference, centre='k-min')
    k_adc = reference.sequence.calculate_kspacePP()[0]
    return measured, checks.spoke_geometry(measured, np.asarray(k_adc))


def main() -> None:
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    results: dict[str, object] = {}

    # --- the full spoke -------------------------------------------------------------
    full = pypulseq_radial.build('write_radial_gre.py', n_spokes=SPOKES, n_dummy=2, n_x=SAMPLES)
    measured, rows = _measure(full)
    wanted = [180.0 * n / SPOKES for n in range(SPOKES)]
    verdicts = checks.check_spokes(rows, expected_angles_deg=wanted)

    logging.info('write_radial_gre.py -- %d spokes of %d samples, %d blocks',
                 len(rows), rows[0]['centre_sample'] * 2, measured['l0']['n_blocks'])
    logging.info('  %-12s %-10s %-12s %-12s %-10s', 'angle/deg', 'centre', 'k range', 'dk', 'kz')
    for row in rows[:4]:
        logging.info('  %-12.4f %-10s %-12s %-12.4f %-10.1e', row['angle_deg'],
                     f"{row['centre_sample']}/{row['index'] and '' or ''}{SAMPLES}",
                     f"{row['k_min_per_m']:.1f}..{row['k_max_per_m']:.1f}",
                     row['dk_per_m'], row['kz_max_per_m'])
    for verdict in verdicts:
        logging.info('  %-58s %-10.3g %-6s (%s)', verdict['what'], verdict['worst'],
                     'pass' if verdict['pass'] else 'FAIL', verdict['invariant'])
    results['full_spoke'] = {'rows': rows, 'checks': verdicts}

    # --- rotation equivariance, measured against the reference itself ----------------
    k_adc = np.asarray(full.sequence.calculate_kspacePP()[0]).reshape(3, len(rows), -1)
    equivariance = [checks.rotation_equivariance(k_adc[:, 0], k_adc[:, n], wanted[n])
                    for n in range(1, len(rows))]
    worst = max(equivariance, key=lambda row: row['worst'])
    logging.info('  %-58s %-10.3g %-6s (R7)', 'every spoke is spoke 0 rotated',
                 worst['worst'], 'pass' if worst['pass'] else 'FAIL')
    results['rotation'] = equivariance

    # --- full spoke to centre-out, in one script -------------------------------------
    logging.info('\nwrite_ute.py -- one `readout_asymmetry` argument across the family')
    logging.info('  %-12s %-14s %-22s %-10s', 'asymmetry', 'centre sample', 'signed k range', 'dk')
    ute_rows = {}
    for asymmetry in (0.0, 0.25, 0.5, 0.75, 1.0):
        reference = pypulseq_radial.build('write_ute.py', n_spokes=4, n_x=SAMPLES,
                                          readout_asymmetry=asymmetry)
        _, rows_ute = _measure(reference)
        row = rows_ute[0]
        logging.info('  %-12.2f %-14s %-22s %-10.4f', asymmetry,
                     f"{row['centre_sample']}/{2 * SAMPLES}",
                     f"{row['k_min_per_m']:.1f} .. {row['k_max_per_m']:.1f}", row['dk_per_m'])
        ute_rows[asymmetry] = row
    results['ute'] = ute_rows

    out = Path(__file__).with_name('archaeology_result.json')
    out.write_text(json.dumps(results, indent=1, default=float))
    logging.info('\n-> %s', out)


if __name__ == '__main__':
    main()
