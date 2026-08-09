"""
The module library: reusable building blocks, grouped by what they do.

Every one of these is designed in ``__init__`` and assembled by ``build()`` into a
:class:`~seqcraft.core.logic.LogicBlock`, on top of the shared :class:`~seqcraft.modules.base.Module`
base.  Adding a new one means writing a class in the right subpackage and exporting it here;
nothing in the core changes.

============  ==========================================================================
Subpackage    Contents
============  ==========================================================================
``base``      :class:`Module`, the optional convenience base these share
``rf``        excitation, refocusing, inversion, saturation
``readout``   Cartesian lines, variable-density spirals, noise scans
``encoding``  phase and partition encodes, prephasers, spoilers, crushers, diffusion
``prep``      fat saturation, inversion recovery
``control``   delays, triggers, barriers, raw events
============  ==========================================================================

**This is a library, not an interface.**  Nothing obliges your own components to live here, to
inherit :class:`Module`, or to be classes at all -- seqcraft's requirement is that whatever you
hand the compiler is a ``LogicBlock``.  What the base buys you is the bookkeeping these modules
all needed anyway: the scanner, the resolved limit regime, a unit check and parameters for
provenance.

Examples
--------
>>> import seqcraft as sc
>>> sc.modules.SpiralVDS.__name__
'SpiralVDS'

There is no registry and no string-keyed lookup: the class *is* the name.  What used to need one --
running the contract suite over the module library -- is :func:`seqcraft.testing.module_subclasses`,
which walks ``Module.__subclasses__()`` and therefore cannot be forgotten:

>>> 'SpiralVDS' in sc.testing.module_subclasses()
True
"""

from __future__ import annotations

from . import base, control, encoding, prep, readout, rf
from .base import Module
from .control import Barrier, Delay, RawEvents, Trigger
from .encoding import (
    ArbitraryDiffusion,
    BipolarDiffusion,
    Crusher,
    MonopolarDiffusion,
    PartitionEncode,
    PhaseEncode,
    Prephaser,
    Spoiler,
    direction_condition_number,
    dti_directions,
)
from .prep import FatSat, InversionRecovery
from .readout import CartesianLine, EPIReadout, NoiseAcquisition, SpiralVDS, vds_trajectory
from .rf import (
    AdiabaticInversion,
    GaussSaturation,
    HardExcitation,
    HardRefocusing,
    RefocusingPulse,
    RFPulse,
    SincExcitation,
    SincRefocusing,
    SlabExcitation,
    SLRExcitation,
    SLRRefocusing,
)

__all__ = [
    'AdiabaticInversion',
    'ArbitraryDiffusion',
    'Barrier',
    'BipolarDiffusion',
    'CartesianLine',
    'Crusher',
    'Delay',
    'EPIReadout',
    'FatSat',
    'GaussSaturation',
    'HardExcitation',
    'HardRefocusing',
    'InversionRecovery',
    'Module',
    'MonopolarDiffusion',
    'NoiseAcquisition',
    'PartitionEncode',
    'PhaseEncode',
    'Prephaser',
    'RFPulse',
    'RawEvents',
    'RefocusingPulse',
    'SLRExcitation',
    'SLRRefocusing',
    'SincExcitation',
    'SincRefocusing',
    'SlabExcitation',
    'Spoiler',
    'SpiralVDS',
    'Trigger',
    'base',
    'control',
    'direction_condition_number',
    'dti_directions',
    'encoding',
    'prep',
    'readout',
    'rf',
    'vds_trajectory',
]
