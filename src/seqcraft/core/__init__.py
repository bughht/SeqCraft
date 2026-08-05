"""
Core layer: the data model, the compiler, and the arithmetic they rest on.

Import from :mod:`seqcraft` rather than from here for everyday use; this package is the
implementation surface.
"""

from __future__ import annotations

from . import events, ordering, raster, units, validate
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
from .module import Module
from .registry import lookup, register, registered
from .report import Issue, Report, ReportFailed
from .system import Limits, System, load_hardware, synthetic_hardware

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
    'Module',
    'Node',
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
    'events',
    'flatten',
    'load_hardware',
    'lookup',
    'ordering',
    'raster',
    'register',
    'registered',
    'span',
    'synthetic_hardware',
    'units',
    'validate',
]
