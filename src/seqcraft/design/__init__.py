"""
What you build: the tree, the arithmetic under it, and the component convention.

Everything a sequence is expressed *in*, before anything has been compiled.  The dividing line
from :mod:`seqcraft.compiler` is sharp and worth stating: nothing here knows what a pulseq block
is.  A logic block lets anything overlap anything, a raster quantises a time, a unit converts a
number -- and none of them has an opinion about where the boundaries will fall.

:mod:`~seqcraft.design.logic`
    :class:`~seqcraft.design.logic.LogicBlock`, a tree of pulseq events with relative start times.
    Two attributes and one method; overlap is free everywhere.
:mod:`~seqcraft.design.module`
    :class:`~seqcraft.design.module.Module`, the standard shape for a reusable component.  A
    convention, not a gate -- the compiler never asks what produced a tree.
:mod:`~seqcraft.design.events`
    :func:`~seqcraft.design.events.derive`, the one sanctioned way to copy a pypulseq event, plus
    the exact-gradient primitives (``knots_of``, ``pwl_moment``) and the kind vocabulary.
:mod:`~seqcraft.design.timing`
    :class:`~seqcraft.design.timing.Raster` and integer-tick arithmetic, so a time that has been
    added up for minutes still lands where it was asked to.
:mod:`~seqcraft.design.units`
    :func:`~seqcraft.design.units.convert`, one function for every pair of units a sequence deals
    in, and the plausibility bands that catch a millimetre passed as a metre.
``sample`` used to be here, as a sixth module.  It is in :mod:`seqcraft.analysis` now, beside the
other three ways of measuring a tree -- one place to look for "give it a tree, get numbers back",
and one place to state which of them are exact.

Import from :mod:`seqcraft` rather than from here for everyday use.
"""

from __future__ import annotations

from . import events, logic, module, timing, units

__all__ = ['events', 'logic', 'module', 'timing', 'units']
