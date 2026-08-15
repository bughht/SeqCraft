"""
Core layer: the data model, the compiler, and the arithmetic they rest on.

Everything here is required to get from a :class:`~seqcraft.core.logic.LogicBlock` to a legal,
validated ``.seq`` -- and nothing else is.  The component contract (:class:`seqcraft.Module`),
the scanner description (:mod:`seqcraft.scanner`), sequence-programming vocabulary
(:mod:`seqcraft.ordering`) and output tooling (:mod:`seqcraft.provenance`,
:mod:`seqcraft.display`) sit one level up, outside this package, because the compiler references
none of them.

Two consequences of that rule are worth stating, because both were once otherwise:

**There is no module abstraction here.**  The compiler's input is a ``LogicBlock`` and it never
asks what produced one, so :class:`seqcraft.Module` lives beside the components rather than with
the compiler.  A test asserts this package does not import it.

**There is no scanner class here.**  The compiler takes a :class:`pypulseq.Opts` -- eight fields of
it -- rather than a seqcraft wrapper, so the same object that configures ``pp.make_trapezoid``
configures the compile.

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
    'LogicBlock',
    'MissingExtraError',
    'Node',
    'Raster',
    'RasterError',
    'Report',
    'ReportFailed',
    'SeqCraftError',
    'UnitSanityError',
    'UnknownFieldError',
    'WriteResult',
    'barrier',
    'compile_sequence',
    'convert',
    'events',
    'flatten',
    'span',
    'timing',
    'units',
    'validate',
]
