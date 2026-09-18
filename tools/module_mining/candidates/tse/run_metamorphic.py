"""
Step 9 -- truths that need no second implementation.

Each of these compares the sequence against *itself* under a transformation, which is worth having
because the two independent references agree on so much that agreement stops being informative.

``M1``  turbo factor 1 is a spin echo.  Not re-run here: SeqCraft's own CI already pins
        ``FSE2D(echoes=1)`` against what ``examples/se_2d/01`` writes, event for event, in
        ``tests/modules/test_se_notebooks.py``.  Cited rather than duplicated.
``M2``  changing the ky ordering must not change the readout waveform.  The table is data; if it
        reaches the readout's design, something that should have been policy became physics.
``M3``  the module's claimed effective TE must equal the measured time from the excitation to the
        sample carrying ``ky = 0``.  This is the one check that tests a *claim* rather than a
        construction, which is the failure mode no k-space comparison can see.
``M4``  increasing the turbo factor must not move the echoes that were already there.

    python tools/module_mining/candidates/tse/run_metamorphic.py
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from module_mining import fingerprint  # noqa: E402
from module_mining.candidates.tse.adapters import seqcraft_fse  # noqa: E402

N_LINES, ECHOES = 128, 16
GRADIENT_TOL_HZ_M = 1e-6
TIME_TOL_S = 1e-9


def _linear(echoes: int, n_lines: int) -> list[list[int]]:
    return [list(range(s * echoes, (s + 1) * echoes)) for s in range(n_lines // echoes)]


def _echo_times(measured: dict) -> np.ndarray:
    """Time from each readout's own excitation to its echo, in play order."""
    return np.asarray([r['t_echo_s'] - measured['shots'][r['shot']]['t_excitation_s']
                       for r in measured['readouts']])


def main() -> None:
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    namespace = seqcraft_fse.notebook_namespace()
    results = {}

    # M2 -- the same instance, two tables.
    a = seqcraft_fse.build(echoes=ECHOES, namespace=namespace)
    b = seqcraft_fse.build(echoes=ECHOES, namespace=namespace,
                           segments=_linear(ECHOES, N_LINES))
    wave_a, wave_b = fingerprint.measure(a)['_waveforms'], fingerprint.measure(b)['_waveforms']
    gx_error = float(np.max(np.abs(np.asarray(wave_a[0]) - np.asarray(wave_b[0]))))
    gy_change = float(np.max(np.abs(np.asarray(wave_a[1]) - np.asarray(wave_b[1]))))
    results['M2'] = {
        'what': 'reordering the table leaves the readout waveform untouched',
        'gx_max_difference_hz_m': gx_error, 'gy_max_difference_hz_m': gy_change,
        'tolerance_hz_m': GRADIENT_TOL_HZ_M,
        'pass': gx_error <= GRADIENT_TOL_HZ_M and gy_change > GRADIENT_TOL_HZ_M,
        'note': 'gy must change -- a table that changed nothing would mean it was never applied',
    }

    # M3 -- the claim against the measurement.
    measured = fingerprint.measure(a)
    dky = a.semantic['dky_per_m']
    centre = min(measured['readouts'], key=lambda r: abs(r['k_at_echo_per_m'][1]))
    measured_te = (centre['t_echo_s']
                   - measured['shots'][centre['shot']]['t_excitation_s'])
    claimed_te = float(a.semantic['te_eff_s'])
    results['M3'] = {
        'what': 'claimed effective TE equals the measured time to the ky = 0 sample',
        'claimed_s': claimed_te, 'measured_s': float(measured_te),
        'error_s': abs(claimed_te - measured_te),
        'ky_at_that_sample_per_m': centre['k_at_echo_per_m'][1],
        'dky_per_m': dky,
        'tolerance_s': float(namespace['opts'].grad_raster_time),
        'pass': abs(claimed_te - measured_te) <= float(namespace['opts'].grad_raster_time),
    }

    # M4 -- a longer train must not move the echoes already there.
    short = fingerprint.measure(seqcraft_fse.build(echoes=8, namespace=namespace))
    long = fingerprint.measure(seqcraft_fse.build(echoes=ECHOES, namespace=namespace))
    short_te, long_te = _echo_times(short)[:8], _echo_times(long)[:8]
    drift = float(np.max(np.abs(short_te - long_te)))
    results['M4'] = {
        'what': 'the first eight echo times are unchanged when the train doubles',
        'max_drift_s': drift, 'tolerance_s': TIME_TOL_S, 'pass': drift <= TIME_TOL_S,
    }

    for name, row in results.items():
        logging.info('%-4s %-62s %s', name, row['what'], 'pass' if row['pass'] else 'FAIL')
    logging.info('  M2 gx %.3g Hz/m, gy %.3g Hz/m', results['M2']['gx_max_difference_hz_m'],
                 results['M2']['gy_max_difference_hz_m'])
    logging.info('  M3 claimed %.6f ms, measured %.6f ms', results['M3']['claimed_s'] * 1e3,
                 results['M3']['measured_s'] * 1e3)
    logging.info('  M4 drift %.3g s', results['M4']['max_drift_s'])

    out = Path(__file__).with_name('metamorphic_result.json')
    out.write_text(json.dumps(results, indent=1, default=float))
    logging.info('-> %s', out)


if __name__ == '__main__':
    main()
