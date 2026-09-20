"""
Read an emitted ``.seq`` and report the physical story it tells, without reading the Python.

The review question this exists to answer:

    If I inspected this file knowing nothing about the module that wrote it, would its RF,
    gradients, ADC and timing tell the same physical story as the API name and the docstring?

That is a different question from the ones the rest of this package asks. A candidate can satisfy
every k-space and timing invariant, compile legally, and still emit the wrong **pulse family** --
which is exactly what happened to the first `GRE3DTR`: its non-selective mode played the slab
path's 3 ms shaped sinc with the selection gradient removed, so the sequence was spatially
non-selective, internally consistent, passing, and not the experiment the mode's name promised.

So this is a review aid, not a gate. It prints; it does not assert.

    python tools/module_mining/inspect_emitted.py path/to/one.seq [more.seq ...]
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pypulseq as pp

#: A magnitude this flat is a hard pulse.  Anything shaped varies by far more.
_FLAT_FRACTION = 0.01


def first_rf_block(seq: pp.Sequence) -> tuple[int, Any] | tuple[None, None]:
    """The first block carrying RF, which is where a pulse family is visible."""
    for index in sorted(seq.block_events):
        block = seq.get_block(index)
        if getattr(block, 'rf', None) is not None:
            return index, block
    return None, None


def describe(path: Path) -> dict[str, Any]:
    """What the file says about its own excitation, gradients and readout."""
    seq = pp.Sequence()
    seq.read(str(path))
    index, block = first_rf_block(seq)
    if block is None:                                                 # pragma: no cover
        return {'path': path.name, 'rf': 'none'}

    rf = block.rf
    magnitude = np.abs(np.asarray(rf.signal))
    on_axis = {axis: getattr(block, f'g{axis}', None) is not None for axis in 'xyz'}
    adc_blocks = sum(1 for i in seq.block_events if int(seq.block_events[i][5]))

    return {
        'path': path.name,
        'rf_samples': int(magnitude.size),
        'rf_shape': ('hard / block' if float(np.ptp(magnitude)) / float(magnitude.max())
                     < _FLAT_FRACTION else 'shaped / soft'),
        'rf_duration_ms': float(rf.t[-1] - rf.t[0]) * 1e3 if magnitude.size > 1 else 0.0,
        'gradients_with_the_rf': sorted(axis for axis, present in on_axis.items() if present),
        'spatially_selective': any(on_axis.values()),
        'rf_block': index,
        # `use` and the carrier offset are the whole physical story for a preparation pulse, and
        # neither is visible in the shape.  A saturation labelled 'excitation', or one offset to
        # the wrong side of water, reads as a perfectly ordinary pulse until these are printed.
        'rf_use': getattr(rf, 'use', None) or 'undefined',
        'rf_freq_offset_hz': float(getattr(rf, 'freq_offset', 0.0) or 0.0),
        'rf_freq_ppm': float(getattr(rf, 'freq_ppm', 0.0) or 0.0),
        'adc_events': adc_blocks,
        'blocks': len(seq.block_events),
        'duration_s': float(seq.duration()[0]),
    }


def main(paths: list[str]) -> None:
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    for name in paths:
        row = describe(Path(name))
        logging.info('%s', row['path'])
        logging.info('   first RF        block %s, %d samples, **%s**, %.3f ms',
                     row['rf_block'], row['rf_samples'], row['rf_shape'], row['rf_duration_ms'])
        logging.info('   with the RF     gradients on %s -> %s',
                     row['gradients_with_the_rf'] or 'no axis',
                     'spatially selective' if row['spatially_selective'] else 'NON-selective')
        offset = row['rf_freq_offset_hz']
        side = 'on resonance' if offset == 0.0 else (
            f'{abs(offset):.1f} Hz {"BELOW" if offset < 0 else "ABOVE"} the carrier')
        logging.info('   declared use    **%s**, %s%s', row['rf_use'], side,
                     f", freq_ppm {row['rf_freq_ppm']:+.3f}" if row['rf_freq_ppm'] else '')
        logging.info('   acquisition     %d ADC events over %d blocks, %.2f s',
                     row['adc_events'], row['blocks'], row['duration_s'])
        logging.info('')


if __name__ == '__main__':
    main(sys.argv[1:] or [])
