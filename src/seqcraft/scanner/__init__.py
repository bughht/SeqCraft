"""
Describing the machine: what it may be driven to, and how the body responds.

Two files, kept together because they answer two halves of one question and are wanted at
different moments:

:mod:`seqcraft.scanner.opts`
    The **limits and timing** a sequence is designed and compiled against.  seqcraft describes
    those with the official :class:`pypulseq.Opts` and nothing else -- there is no wrapper here,
    only the two operations ``Opts`` makes awkward: :func:`~seqcraft.scanner.opts.derate` and
    :func:`~seqcraft.scanner.opts.from_scanner`.

:mod:`seqcraft.scanner.hardware`
    The **response model** for peripheral-nerve-stimulation prediction.  Not a limit, and nothing
    on the compile path reads it: :func:`seqcraft.pns` takes it as a third argument, beside the
    tree and the ``Opts``.

Both are re-exported at the package root, so the shorter spelling is the one to use::

    sc.opts.derate(opts, grad=0.85)
    sc.hardware.synthetic_hardware()
"""

from __future__ import annotations

from . import hardware, opts

__all__ = ['hardware', 'opts']
