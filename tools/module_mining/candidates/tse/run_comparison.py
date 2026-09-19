"""
Experiments E0, E2 and E3 -- the two Python references at one protocol, measured side by side.

**E0, the protocol, is a result and not a setting.**  The three references cannot be run at each
other's defaults, and reconciling them is where the first real facts appear:

* the readout is fixed inside ``write_tse.py`` at 6.4 ms of sampling, so SeqCraft is told
  ``bandwidth_hz_px = 1/6.4e-3``;
* ``write_tse.py``'s 2 ms refocusing pulse peaks at 116 % of its own ``max_b1`` -- it warns and
  continues -- and SeqCraft's ``Refocusing`` **refuses** the same pulse at 130 % of the notebook's
  limit.  Rather than raise the limit to force agreement, the comparison keeps SeqCraft's 4 ms
  pulse and matches the *echo spacing* instead, which is the quantity the physics depends on;
* SeqCraft's shortest train-bound spacing at this geometry is 12.120 ms, so that -- not the
  reference's nominal 12 ms -- is the common ``te``.

What is then compared is L0 inventory, L3 acquisition semantics and L4 TSE physics.  What is
**not** compared is the phase-encode table, the crusher amounts or the block structure: those are
recorded design differences (D1, D2, D6 in ``findings.md``), and a comparator that normalised them
away would be proving that it cannot see them.

    python tools/module_mining/candidates/tse/run_comparison.py
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from module_mining import fingerprint  # noqa: E402
from module_mining.candidates.tse.adapters import pypulseq_tse, seqcraft_fse  # noqa: E402

FOV_MM = 256.0
MATRIX = 64
ECHOES = 16
THICKNESS_MM = 5.0
SAMPLING_S = 6.4e-3
#: SeqCraft's shortest train-bound spacing at this geometry with a 4 ms refocusing pulse.
ECHO_SPACING_S = 12.120e-3


def _summary(rows: list[dict]) -> dict[str, dict]:
    return {row['invariant']: {'worst': row['worst'], 'unit': row['unit'], 'pass': row['pass']}
            for row in rows}


def main() -> None:
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    namespace = seqcraft_fse.notebook_namespace()
    r1 = seqcraft_fse.build(
        echoes=ECHOES, namespace=namespace, dummy_shots=1,
        fov_mm=(FOV_MM, FOV_MM), matrix=(MATRIX, MATRIX), thickness_mm=THICKNESS_MM,
        bandwidth_hz_px=1.0 / SAMPLING_S, echo_spacing_s=ECHO_SPACING_S,
        segments=seqcraft_fse.interleaved(ECHOES, MATRIX))
    r2 = pypulseq_tse.build(fov=FOV_MM / 1e3, n_x=MATRIX, n_y=MATRIX, n_echo=ECHOES,
                            slice_thickness=THICKNESS_MM / 1e3, te=ECHO_SPACING_S, tr=2.0)

    measured = {ref.name: fingerprint.measure(ref) for ref in (r1, r2)}
    dky = 1e3 / FOV_MM
    result = {'protocol': {'fov_mm': FOV_MM, 'matrix': MATRIX, 'echoes': ECHOES,
                           'thickness_mm': THICKNESS_MM, 'sampling_s': SAMPLING_S,
                           'echo_spacing_s': ECHO_SPACING_S},
              'invariants': {}, 'coverage': {}, 'comparison': {}}

    for name, m in measured.items():
        expected = (seqcraft_fse.expected_ky(r1) if name == 'seqcraft_fse' else None)
        result['invariants'][name] = _summary(
            fingerprint.check_invariants(m, expected_ky_per_m=expected))
        result['coverage'][name] = fingerprint.check_ky_coverage(
            m, dky_per_m=dky, n_lines=MATRIX)
        logging.info('%s: esp %.4f ms, TE1 %.4f ms, %d readouts, %d blocks',
                     name, m['shots'][0]['echo_spacing_s'] * 1e3,
                     (m['readouts'][0]['t_echo_s']
                      - m['shots'][m['readouts'][0]['shot']]['t_excitation_s']) * 1e3,
                     m['l0']['n_adc'], m['l0']['n_blocks'])

    # Keyed off the references themselves: the SeqCraft adapter names itself after the
    # implementation it measured -- notebook or package -- so a literal would go stale.
    result['comparison'] = fingerprint.compare(measured[r1.name], measured[r2.name])

    logging.info('\n%-12s %-12s %-12s %s', 'invariant', 'seqcraft', 'pypulseq', 'unit')
    keys = sorted(set(result['invariants'][r1.name]) & set(result['invariants'][r2.name]))
    for key in keys:
        a = result['invariants'][r1.name][key]
        b = result['invariants'][r2.name][key]
        logging.info('%-12s %-12.4g %-12.4g %s', key, a['worst'], b['worst'], a['unit'])

    logging.info('\ncoverage: %s', json.dumps(
        {n: {k: v for k, v in c.items() if k in ('n_acquired', 'n_distinct', 'duplicated',
                                                 'missing', 'lattice_max_error_per_m', 'pass')}
         for n, c in result['coverage'].items()}))
    logging.info('\ncomparison: %s', json.dumps(result['comparison'], default=float))

    out = Path(__file__).with_name('comparison_result.json')
    out.write_text(json.dumps(result, indent=1, default=float))
    logging.info('-> %s', out)


if __name__ == '__main__':
    main()
