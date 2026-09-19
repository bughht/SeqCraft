"""
`GRE3DTR` against the official 3D reference, and across a sweep.

The non-selective path is compared with `pulseq/pulseq`'s ``matlab/demoSeq/gre3d.seq`` -- shipped
by that repository, so no MATLAB run is needed.

**On lattice and semantic centre, not on sample positions.** The reference prephases with
``-gx.area/2`` and has no sample at the centre of k-space: its two central samples sit at
+-dkx/2. SeqCraft's `CartesianLine` puts one there. Both are correct; comparing raw sample
positions would be comparing that design choice rather than the encoding.

The slab-selective path has no executable reference, so it is validated by the degenerate limits
and the signed-moment arithmetic, which is what `tests/modules/test_gre_3d_tr.py` carries.

    python tools/module_mining/candidates/gre3d/run_validation.py
"""

from __future__ import annotations

import json
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pypulseq as pp

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import seqcraft as sc  # noqa: E402

REFERENCE = Path(os.environ.get(
    'MODULE_MINING_PULSEQ', '/Users/yiyund/Code_mgh/pulseq')) / 'matlab/demoSeq/gre3d.seq'

FOV_MM = (200.0, 200.0, 160.0)
MATRIX = (64, 64, 64)
SAMPLES = 64


def _scanner() -> pp.Opts:
    """The reference's own limits: 20 mT/m, and a fixed 400 us rise time it states as a slew."""
    return pp.Opts(max_grad=20, grad_unit='mT/m', max_slew=120, slew_unit='T/m/s', B0=3.0,
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)


def reference_lattices() -> dict[str, dict]:
    """delta-k, index range and distinct count per axis, measured from the shipped .seq."""
    seq = pp.Sequence()
    seq.read(str(REFERENCE))
    k = np.asarray(seq.calculate_kspacePP()[0]).reshape(3, -1, SAMPLES)
    echo = k[:, :, SAMPLES // 2 - 1]

    out: dict[str, dict] = {}
    for axis_index, (axis, fov_m) in enumerate((('x', 0.2), ('y', 0.2), ('z', 0.16))):
        dk = 1.0 / fov_m
        if axis == 'x':
            out[axis] = {'dk_per_m': dk, 'k_at_central_samples': [float(k[0, 0, SAMPLES // 2 - 1]),
                                                                  float(k[0, 0, SAMPLES // 2])]}
            continue
        indices = np.round(echo[axis_index] / dk).astype(int)
        out[axis] = {'dk_per_m': dk, 'min_index': int(indices.min()),
                     'max_index': int(indices.max()), 'distinct': int(np.unique(indices).size),
                     'off_lattice_per_m': float(np.abs(echo[axis_index] - indices * dk).max())}
    return out


def candidate_lattices(opts) -> dict[str, dict]:
    """The same quantities from `GRE3DTR`, over the centre and both edges of each axis."""
    tr = sc.modules.GRE3DTR(opts=opts, fov_mm=FOV_MM, matrix=MATRIX, bandwidth_hz_px=250.0)
    probes = [(0, 0), (tr.center_line, tr.center_partition), (MATRIX[1] - 1, MATRIX[2] - 1)]

    rows = []
    for line, partition in probes:
        k = sc.kspace(tr(line=line, partition=partition), opts)['k_adc']
        echo = k[:, tr.ro.echo_sample(0)]
        rows.append({'line': line, 'partition': partition,
                     'kx_at_echo': float(echo[0]),
                     'ky': float(echo[1]), 'ky_index': int(round(echo[1] / tr.dk_per_m('y'))),
                     'kz': float(echo[2]), 'kz_index': int(round(echo[2] / tr.dk_per_m('z'))),
                     'kz_held_per_m': float(np.ptp(k[2]))})
    return {'dk_per_m': {a: tr.dk_per_m(a) for a in 'xyz'},
            'center': {'line': tr.center_line, 'partition': tr.center_partition},
            'probes': rows, 'te_s': tr.te_s, 'tr_s': tr.tr_s, 'winder_s': tr.winder_s}


SWEEP = [
    ('reference geometry', {}),
    ('unequal ny/nz', {'matrix': (32, 16, 8)}),
    ('thin slab, many partitions', {'fov_mm': (200.0, 200.0, 40.0), 'matrix': (32, 16, 64)}),
    ('slab = fov_z', {'matrix': (32, 16, 8), 'slab_thickness_mm': 160.0}),
    ('slab < fov_z', {'matrix': (32, 16, 8), 'slab_thickness_mm': 112.0}),
    ('slab > fov_z', {'matrix': (32, 16, 8), 'slab_thickness_mm': 200.0}),
    ('sharper slab', {'matrix': (32, 16, 8), 'slab_thickness_mm': 160.0,
                      'rf_time_bw_product': 8.0}),
]


def main() -> None:
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    opts = _scanner()
    results: dict[str, object] = {}

    if REFERENCE.exists():
        reference = reference_lattices()
        candidate = candidate_lattices(opts)
        results['reference'] = reference
        results['candidate'] = candidate

        logging.info('reference gre3d.seq   dk = %.4f / %.4f / %.4f 1/m (x/y/z)',
                     reference['x']['dk_per_m'], reference['y']['dk_per_m'],
                     reference['z']['dk_per_m'])
        logging.info('candidate GRE3DTR     dk = %.4f / %.4f / %.4f 1/m',
                     *(candidate['dk_per_m'][a] for a in 'xyz'))
        logging.info('reference ky %+d .. %+d (%d distinct) | kz %+d .. %+d (%d distinct)',
                     reference['y']['min_index'], reference['y']['max_index'],
                     reference['y']['distinct'], reference['z']['min_index'],
                     reference['z']['max_index'], reference['z']['distinct'])
        logging.info('reference has no DC sample on x: central samples at %+.3f / %+.3f 1/m',
                     *reference['x']['k_at_central_samples'])
        for row in candidate['probes']:
            logging.info('  candidate line %2d partition %2d -> ky index %+3d, kz index %+3d, '
                         '|kx| at echo %.1e, kz spread %.1e',
                         row['line'], row['partition'], row['ky_index'], row['kz_index'],
                         abs(row['kx_at_echo']), row['kz_held_per_m'])
    else:                                                             # pragma: no cover
        logging.warning('reference not found at %s -- lattice comparison skipped', REFERENCE)

    logging.info('\n%-28s %-9s %-10s %-9s %-9s %-8s', 'sweep case', 'winder', 'limit p',
                 'TE/ms', 'TR/ms', 'A_slab')
    sweep = {}
    for name, overrides in SWEEP:
        settings = {'fov_mm': FOV_MM, 'matrix': MATRIX, 'bandwidth_hz_px': 250.0, **overrides}
        tr = sc.modules.GRE3DTR(opts=opts, **settings)
        sweep[name] = {'winder_s': tr.winder_s, 'limiting_partition': tr.limiting_partition,
                       'te_s': tr.te_s, 'tr_s': tr.tr_s,
                       'slab_rephase_area_per_m': tr.slab_rephase_area_per_m}
        logging.info('%-28s %-9.1f %-10d %-9.3f %-9.3f %-8.2f', name, tr.winder_s * 1e6,
                     tr.limiting_partition, tr.te_s * 1e3, tr.tr_s * 1e3,
                     tr.slab_rephase_area_per_m)
    results['sweep'] = sweep

    out = Path(__file__).with_name('validation_result.json')
    out.write_text(json.dumps(results, indent=1, default=float))
    logging.info('\n-> %s', out)


if __name__ == '__main__':
    main()
