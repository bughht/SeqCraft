"""
Step 8 -- a small, high-information parameter sweep of the SeqCraft implementation.

One nominal protocol proves nothing about a reusable abstraction, and a full Cartesian product
proves little more at much greater cost.  These cases were chosen so that each one can fail on
its own account:

* ``echoes`` 1, 8, 16 and 72 -- a spin echo, two trains, and a single-shot HASTE;
* a second matrix, because the readout's own asymmetry about the echo changes with it;
* ``partial_fourier`` below 1, where the echo stops being near the middle of the readout;
* three orderings on one instance, which must not change the waveform at all;
* a requested echo spacing longer than the minimum, which is the case where the echo would drift
  off the midpoint if the shift that corrects it were missing.

    python tools/module_mining/candidates/tse/run_sweep.py
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from module_mining import fingerprint  # noqa: E402
from module_mining.candidates.tse.adapters import seqcraft_fse  # noqa: E402

#: The invariants reported per case, in the order the table lists them.
WATCHED = ('I1', 'I2', 'I3', 'I4x', 'I4y', 'I4z', 'I5', 'I6', 'I7')


def _centric(echoes: int, n_lines: int) -> list[list[int]]:
    order = sorted(range(n_lines), key=lambda line: (abs(line - n_lines // 2), line))
    shots = n_lines // echoes
    return [[order[n * shots + s] for n in range(echoes)] for s in range(shots)]


def _linear(echoes: int, n_lines: int) -> list[list[int]]:
    return [list(range(s * echoes, (s + 1) * echoes)) for s in range(n_lines // echoes)]


def cases() -> list[tuple[str, dict]]:
    """``(name, kwargs for seqcraft_fse.build)``."""
    return [
        ('echoes=1', {'echoes': 1}),
        ('echoes=8', {'echoes': 8}),
        ('echoes=16', {'echoes': 16}),
        ('matrix=64', {'echoes': 16, 'matrix': (64, 64), 'fov_mm': (256.0, 256.0),
                       'segments': seqcraft_fse.interleaved(16, 64)}),
        ('linear', {'echoes': 16, 'segments': _linear(16, 128)}),
        ('centric', {'echoes': 16, 'segments': _centric(16, 128)}),
        ('esp +2 ms', {'echoes': 16, 'echo_spacing_s': 12.7e-3}),
        ('haste', {'echoes': 72, 'partial_fourier': 0.625,
                   'segments': [[128 // 2 - 8 + n for n in range(72)]]}),
    ]


def main() -> None:
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    namespace = seqcraft_fse.notebook_namespace()
    source = sys.argv[1] if len(sys.argv) > 1 else 'notebook'
    results: dict[str, dict] = {'_source': source}

    header = '%-12s ' + ' '.join(['%-10s'] * len(WATCHED)) + ' %-9s %-8s'
    logging.info(header, 'case', *WATCHED, 'esp/ms', 'bound')
    for name, kwargs in cases():
        reference = seqcraft_fse.build(namespace=namespace, source=source, **kwargs)
        measured = fingerprint.measure(reference)
        rows = {row['invariant']: row for row in fingerprint.check_invariants(
            measured, expected_ky_per_m=seqcraft_fse.expected_ky(reference))}
        results[name] = {
            'parameters': {k: v for k, v in kwargs.items() if k != 'segments'},
            'semantic': {k: reference.semantic[k]
                         for k in ('echo_spacing_s', 'te_first_s', 'te_eff_s', 'bound', 'shift_s')},
            'invariants': {k: {'worst': rows[k]['worst'], 'unit': rows[k]['unit'],
                               'pass': rows[k]['pass']} for k in WATCHED},
            'all_pass': all(rows[k]['pass'] for k in WATCHED),
        }
        logging.info('%-12s ' + ' '.join(['%-10.3g'] * len(WATCHED)) + ' %-9.4f %-8s',
                     name, *(rows[k]['worst'] for k in WATCHED),
                     reference.semantic['echo_spacing_s'] * 1e3, reference.semantic['bound'])

    failed = [name for name, row in results.items()
              if isinstance(row, dict) and not row.get('all_pass', True)]
    logging.info('\n%s: %d/%d cases pass every invariant%s', source,
                 len(results) - 1 - len(failed), len(results) - 1,
                 '' if not failed else f'; failed: {failed}')

    out = Path(__file__).with_name(f'sweep_result_{source}.json')
    out.write_text(json.dumps(results, indent=1, default=float))
    logging.info('-> %s', out)


if __name__ == '__main__':
    main()
