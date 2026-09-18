"""
Experiment E1 -- does the references' construction depend on ``rf_dead_time == rf_ringdown_time``?

Both external TSE references place the RF inside a selection plateau of length
``tEx + ringdown + dead`` and delay the pulse by ``dead``.  The effective centre then sits at
``dead + tEx/2`` and the plateau midpoint at ``(tEx + ringdown + dead)/2``: the two coincide only
when the dead time and the ringdown are equal, which both references set them to be.  When they
are not, the area before the centre and the area after it differ by
``a_sel * (dead - ringdown) / 2``, and because a refocusing pulse conjugates k that residual
**alternates sign echo to echo**.

SeqCraft's ``Refocusing`` removes the assumption: it symmetrises the plateau about the measured
effective centre and then solves the crusher amplitudes by integration.  Its own notebook runs at
``rf_ringdown_time = 30 us``, which is exactly the case in question.

This script runs both implementations at both ringdowns and measures, rather than arguing.  If R2
does **not** degrade at 30 us, the analysis in ``findings.md`` is wrong and the invariant table
needs revising before anything is extracted.

    python tools/module_mining/candidates/tse/run_ringdown.py
"""

from __future__ import annotations

import copy
import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from module_mining import fingerprint  # noqa: E402
from module_mining.candidates.tse.adapters import pypulseq_tse, seqcraft_fse  # noqa: E402

#: 100 us is what both references assume.  40 us is asymmetric *and* keeps the reference's derived
#: windows on the gradient raster, so its construction still runs.  30 us is SeqCraft's own
#: notebook value, at which the reference cannot be built at all -- recorded rather than skipped.
RINGDOWNS_S = (100e-6, 40e-6, 30e-6)
ECHOES = 16

#: The invariants this experiment is about: the slice-axis balance, and the timing half of it.
WATCHED = ('I3', 'I6', 'I4z', 'I1-dc')


def _rows(build) -> dict[str, dict]:
    """Build, measure and check -- or record why the reference could not be built at all."""
    try:
        reference = build()
    except Exception as error:                                       # noqa: BLE001
        return {'_failed': {'error': f'{type(error).__name__}: {error}'}}
    measured = fingerprint.measure(reference)
    rows = {row['invariant']: row for row in fingerprint.check_invariants(measured)}
    rows['_context'] = {'n_readouts': measured['l0']['n_adc'],
                        'echo_spacing_s': measured['shots'][0]['echo_spacing_s']}
    return rows


def main() -> None:
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    namespace = seqcraft_fse.notebook_namespace()
    results: dict[str, dict] = {}

    for ringdown in RINGDOWNS_S:
        # Copy the notebook's own scanner rather than restating it, so only the ringdown differs.
        opts = copy.copy(namespace['opts'])
        opts.rf_ringdown_time = ringdown
        results[f'seqcraft_fse@{ringdown * 1e6:.0f}us'] = _rows(
            lambda opts=opts: seqcraft_fse.build(echoes=ECHOES, namespace=namespace, opts=opts))
        results[f'pypulseq_tse@{ringdown * 1e6:.0f}us'] = _rows(
            lambda r=ringdown: pypulseq_tse.build(
                n_echo=ECHOES, opts_overrides={'rf_ringdown_time': r}))

    logging.info('%-28s %-12s %-14s %-14s %-10s', 'reference', *WATCHED)
    for name, rows in results.items():
        if '_failed' in rows:
            logging.info('%-28s %s', name, rows['_failed']['error'])
            continue
        logging.info('%-28s %-12.3e %-14.3e %-14.3e %-10.3f', name,
                     *(rows[key]['worst'] for key in WATCHED))
    logging.info('\nunits: I3 1/m (kz at the echo), I6 s (echo vs midpoint), '
                 'I4z 1/m (area balance), I1-dc dwells (nearest sample to kx = 0)')

    out = Path(__file__).with_name('ringdown_result.json')
    out.write_text(json.dumps(
        {name: {k: v for k, v in rows.items()} for name, rows in results.items()},
        indent=1, default=float))
    logging.info('-> %s', out)


if __name__ == '__main__':
    main()
