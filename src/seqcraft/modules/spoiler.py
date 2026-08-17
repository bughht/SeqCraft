"""
:func:`spoiler` -- a gradient that winds a chosen number of phase turns across one voxel.

A function rather than a :class:`~seqcraft.Module` subclass, and that is the rule the top of
``modules/`` exists for: **the role folders hold module subclasses; the top level holds what is
not one.**  This is one line of arithmetic over ``pp.make_trapezoid``.  It designs nothing per
call, answers no timing question, and holds no state, so a class around it would add an
``__init__``, a ``build`` and a tag for no behaviour at all.

Gradient spoiling and RF spoiling are different mechanisms and are in different places: the
quadratic phase increment is a scalar schedule with no waveform, so it belongs to
:class:`~seqcraft.modules.GRE2D`, which is the layer that knows how many repetitions there are.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pypulseq as pp

from ..design.logic import LogicBlock
from ._support import require_axis, require_positive

if TYPE_CHECKING:
    from pypulseq.opts import Opts

__all__ = ['spoiler']


def spoiler(
    opts: Opts,
    *,
    cycles_per_voxel: float = 4.0,
    voxel_mm: float,
    axis: str = 'z',
) -> LogicBlock:
    r"""
    Return a spoiler gradient that dephases `cycles_per_voxel` turns across one voxel.

    A gradient of area :math:`A` in 1/m imposes :math:`\varphi = 2\pi A \Delta x` across a
    distance :math:`\Delta x`, so the area that buys a whole number of turns is

    .. math:: A = \frac{\text{cycles per voxel}}{\text{voxel}_\mathrm{m}}

    Four to eight cycles is the usual range.  The same choice is quoted in the literature as
    "4\ :math:`\pi` dephasing" and as "2 cycles per voxel" under different halving conventions,
    which is why this takes the formula's input rather than shipping a bare default area.

    Parameters
    ----------
    opts
        The scanner.  Positional, because this is a function and `opts` is what it is *about*.
    cycles_per_voxel
        Full 2\ :math:`\pi` turns of phase wound across one voxel along `axis`.
    voxel_mm
        The voxel dimension **along** `axis`, in millimetres.  **No default.**  For the usual
        ``axis='z'`` that is the *slice thickness*, not the in-plane voxel: passing the in-plane
        size under-spoils by the ratio of the two, and the symptom is faint residual banding
        that looks like anything except a spoiler bug.
    axis
        Logical gradient channel.

    Returns
    -------
    LogicBlock
        One trapezoid, tagged ``'spoiler'``.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s')
    >>> block = spoiler(o, cycles_per_voxel=4.0, voxel_mm=5.0)
    >>> block
    LogicBlock(spoiler, 1 node, 0.74 ms)
    >>> round(float(block.nodes[0].item.area), 6)               # 4 turns across 5 mm
    800.0
    """
    area_per_m = require_positive(cycles_per_voxel, 'cycles_per_voxel') / (
        require_positive(voxel_mm, 'voxel_mm') / 1e3
    )
    gradient = pp.make_trapezoid(channel=require_axis(axis), area=area_per_m, system=opts)
    return LogicBlock('spoiler').add(0.0, gradient)
