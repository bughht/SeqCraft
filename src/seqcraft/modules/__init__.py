"""
The module library: reusable building blocks, grouped by what they do.

Every one of these is a plain :class:`~seqcraft.core.module.Module` -- designed in ``__init__``,
assembled by ``build()`` into a :class:`~seqcraft.core.logic.LogicBlock`.  Adding a new one means
writing a class in the right subpackage and exporting it here; nothing in the core changes.

============  ==========================================================================
Subpackage    Contents
============  ==========================================================================
``rf``        excitation, refocusing, inversion, saturation
``readout``   Cartesian lines, variable-density spirals, noise scans
``encoding``  phase and partition encodes, prephasers, spoilers, crushers, diffusion
``prep``      fat saturation, inversion recovery
``control``   delays, triggers, barriers, raw events
============  ==========================================================================

Examples
--------
>>> import seqcraft as sc
>>> sc.modules.SpiralVDS.__name__
'SpiralVDS'

There is no registry and no string-keyed lookup: the class *is* the name.  What used to need one --
running the contract suite over every module -- is :func:`seqcraft.testing.all_modules`, which walks
``Module.__subclasses__()`` and therefore cannot be forgotten:

>>> 'SpiralVDS' in sc.testing.all_modules()
True
"""

from __future__ import annotations

from . import control, encoding, prep, readout, rf
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
from .readout import CartesianLine, NoiseAcquisition, SpiralVDS, vds_trajectory
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
    'FatSat',
    'GaussSaturation',
    'HardExcitation',
    'HardRefocusing',
    'InversionRecovery',
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
    'control',
    'direction_condition_number',
    'dti_directions',
    'encoding',
    'prep',
    'readout',
    'rf',
    'vds_trajectory',
]
