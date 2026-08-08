"""
Core layer: the data model, the compiler, and the arithmetic they rest on.

Everything here is required to get from a :class:`~seqcraft.core.logic.LogicBlock` to a legal,
validated ``.seq`` -- and nothing else is.  The module library (:mod:`seqcraft.modules`),
sequence-programming vocabulary (:mod:`seqcraft.ordering`) and output tooling
(:mod:`seqcraft.provenance`, :mod:`seqcraft.display`) sit one level up, outside this package,
because the compiler never references them.

In particular there is no module abstraction here.  The compiler's input is a ``LogicBlock`` and
it never asks what produced one, so :class:`seqcraft.modules.base.Module` -- an optional
convenience for writing reusable components -- lives with the components rather than with the
compiler.

Import from :mod:`seqcraft` rather than from here for everyday use; this package is the
implementation surface.
"""

from __future__ import annotations

from . import events, timing, units, validate
from .compiler import CompiledSequence, WriteResult, compile_sequence
from .errors import (
    CompileError,
    ConfigurationError,
    DefinitionConflict,
    HardwareLimitError,
    MissingExtraError,
    RasterError,
    SeqCraftError,
    UnitSanityError,
    UnknownFieldError,
)
from .geometry import Geometry
from .logic import Item, LogicBlock, Node, barrier, flatten, span
from .report import Issue, Report, ReportFailed
from .system import Limits, System, load_hardware, synthetic_hardware
from .timing import Raster
from .units import convert

__all__ = [
    'CompileError',
    'CompiledSequence',
    'ConfigurationError',
    'DefinitionConflict',
    'Geometry',
    'HardwareLimitError',
    'Issue',
    'Item',
    'Limits',
    'LogicBlock',
    'MissingExtraError',
    'Node',
    'Raster',
    'RasterError',
    'Report',
    'ReportFailed',
    'SeqCraftError',
    'System',
    'UnitSanityError',
    'UnknownFieldError',
    'WriteResult',
    'barrier',
    'compile_sequence',
    'convert',
    'events',
    'flatten',
    'load_hardware',
    'span',
    'synthetic_hardware',
    'timing',
    'units',
    'validate',
]
