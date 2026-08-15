"""
A tree as arrays: ``sample(root, opts)``.

Turning a :class:`~seqcraft.design.logic.LogicBlock` into gradient waveforms on a uniform grid,
plus the RF and ADC spans, with no plotting library in sight.  It was the private front half of
:func:`seqcraft.display.plot_block`, and it is here because the operation is useful without
drawing: checking a moment numerically, feeding a simulator, writing a test that asserts on a
waveform rather than on a picture.

**It samples the tree, not the compiled sequence** -- deliberately.  It shows what the module
meant, before the compiler decided where the block boundaries go, which is what you want when the
question is "did I place this correctly".  What the compiler then made of it is what
:func:`seqcraft.plot_sequence` draws, and the difference between the two pictures is exactly the
block structure.

Sampling is for looking, not for arithmetic
-------------------------------------------
The grid here is uniform, which is convenient and *lossy*: an arbitrary gradient's samples sit at
raster centres and an extended trapezoid's knots are not uniformly spaced at all, so both are
interpolated onto the grid rather than reproduced on it.  For a moment, a split or a sum, use
:func:`~seqcraft.design.events.knots_of` and :func:`~seqcraft.design.events.pwl_moment`, which are
exact; the compiler does, and the difference is a 2.5 % amplitude error on a spiral rather than a
rounding detail.

Examples
--------
>>> import numpy as np, pypulseq as pp, seqcraft as sc
>>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s')
>>> block = sc.LogicBlock('spoiler').add(0.0, pp.make_trapezoid('z', area=500.0, system=opts))
>>> grid, grads, marks = sc.sample(block, opts)
>>> sorted(grads), marks
(['z'], [])
>>> bool(np.allclose(np.diff(grid), opts.grad_raster_time))
True
>>> round(float(np.trapezoid(grads['z'], grid)), 1)          # the area, near enough to see
500.0
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from . import events as ev
from .logic import BARRIER, flatten

if TYPE_CHECKING:
    from pypulseq.opts import Opts

    from .logic import LogicBlock

__all__ = ['sample']


def sample(
    root: LogicBlock,
    opts: Opts,
) -> tuple[np.ndarray, dict[str, np.ndarray], list[tuple]]:
    """
    Sample a tree onto a uniform grid: gradients per axis, plus the RF and ADC spans.

    Parameters
    ----------
    root
        The block to sample.  Nested blocks are flattened, so a component that nests is measured
        by its actual events rather than by its direct children.
    opts
        Supplies the gradient raster, which is the grid spacing.

    Returns
    -------
    grid
        Times in seconds from the start of `root`, spaced by ``opts.grad_raster_time``.
    grads
        ``channel -> amplitudes`` in Hz/m on `grid`, summed across everything on that axis.  Only
        the axes actually used appear.
    marks
        ``(kind, start, end, label)`` for each RF, ADC and barrier, in tree order.  ``kind`` is
        ``'rf'``, ``'adc'`` or ``'barrier'``; a barrier has ``start == end``.
    """
    raster = float(opts.grad_raster_time)
    placed = list(flatten(root))
    if not placed:
        return np.zeros(1), {}, []

    total = root.duration
    n = max(2, int(round(total / raster)) + 1)
    grid = np.arange(n) * raster
    grads: dict[str, np.ndarray] = {}
    marks: list[tuple] = []

    for start, event, path in placed:
        kind = getattr(event, 'type', None)
        delay = float(getattr(event, 'delay', 0.0) or 0.0)
        label = '.'.join(path) or kind or '?'
        if kind in ('trap', 'grad'):
            tt, wf = ev.waveform_of(event, raster)
            values = np.interp(grid - start, tt, wf, left=0.0, right=0.0)
            grads[event.channel] = grads.get(event.channel, np.zeros(n)) + values
        elif kind == 'rf':
            marks.append(('rf', start + delay, start + delay + float(event.shape_dur), label))
        elif kind == 'adc':
            span = float(event.num_samples) * float(event.dwell)
            marks.append(('adc', start + delay, start + delay + span, label))
        elif kind == BARRIER:
            marks.append(('barrier', start, start, label))
    return grid, grads, marks
