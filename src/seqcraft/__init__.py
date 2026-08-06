"""
seqcraft -- composable, verifiable MRI pulse sequence programming on top of pypulseq.

Three concepts, and no more.

:class:`~seqcraft.core.logic.LogicBlock`
    A tree of pulseq events and nested blocks, each with a start time.  Two attributes and one
    method.  Anything may overlap anything.
:class:`~seqcraft.core.module.Module`
    A reusable sequence task.  ``__init__`` designs, ``build()`` returns a logic block, and
    timing a caller needs in order to place it is a plain property.
:func:`~seqcraft.core.compiler.compile_sequence`
    Turns a tree into legal pulseq blocks: finds block boundaries, sums gradients that share an
    axis, and validates the result against the amplifier.

Getting started
---------------
>>> import pypulseq as pp
>>> import seqcraft as sc
>>> system = sc.System.preset('generic_3t')
>>> exc = sc.modules.SincExcitation(system, flip_deg=15, duration_us=1000,
...                                 slice_thickness_mm=5)
>>> ro = sc.modules.CartesianLine(system, fov_ro_mm=250, matrix_ro=64,
...                               readout_duration_us=3200)
>>> seq = sc.LogicBlock('demo')
>>> _ = seq.add(0.0, exc.build())
>>> _ = seq.add(8e-3 - ro.time_to_echo, ro.build())     # this defines TE
>>> out = sc.compile(seq, system)
>>> out.check().ok
True

When to reach for this
----------------------
For one sequence, once, write raw pypulseq -- it is a fine tool for that.  seqcraft pays for
itself when you have a *family* of sequences, a loop over more than two axes, gradients that must
overlap without you hand-splitting blocks, or files that have to be reproducible six months
later.

There are no recipes here on purpose.  A recipe is somebody else's sequence choices baked into
library code, and changing your own sequence should never mean editing a package.  The notebooks in
``examples/`` build their sequences from modules, start to finish -- copy one and edit it.  And when
it does not fit, :class:`~seqcraft.modules.control.basic.RawEvents` and the ``CompiledSequence.seq``
attribute are always there.

Notes
-----
Importing this package is cheap and side-effect free: it does not import matplotlib, does not
touch ``pypulseq.Opts.default``, and prints nothing.
"""

from __future__ import annotations

import importlib

from . import _compat, modules, ordering
from ._version import __version__
from .core import events, timing, units, validate
from .core.compiler import CompiledSequence, WriteResult
from .core.compiler import compile_sequence as compile  # noqa: A001, A004
from .core.errors import (
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
from .core.geometry import Geometry
from .core.logic import Item, LogicBlock, Node, barrier, flatten, span
from .core.module import Module
from .core.report import Issue, Report, ReportFailed
from .core.system import Limits, System, load_hardware, synthetic_hardware
from .core.timing import Raster
from .core.units import convert
from .ordering import (
    bit_reversed_order,
    centric_order,
    golden_angle,
    interleaved_slice_order,
    linear_order,
    rf_spoil_phase,
)

# Fail once with a complete list, rather than letting the first module that needs a missing
# pypulseq function fail with an opaque AttributeError halfway through a build.
_compat.require()

#: Also available as ``sc.compile_sequence``, for callers who would rather not shadow the
#: builtin name in their own namespace.
compile_sequence = compile

#: Attributes resolved on first access rather than at import time.
#:
#: ``display`` is the only module allowed to import matplotlib, and ``import seqcraft`` must not
#: pull matplotlib in -- so ``sc.plot_block(exc.build())`` works without paying that cost until
#: it is used.  ``testing`` is deferred for tidiness: the contract assertions are documented as
#: available downstream, but they are not part of the sequence-building API.
_LAZY: dict[str, tuple[str, str | None]] = {
    'display': ('.display', None),
    'plot_block': ('.display', 'plot_block'),
    'plot_kspace': ('.display', 'plot_kspace'),
    'plot_sequence': ('.display', 'plot_sequence'),
    'plot_trajectory': ('.display', 'plot_trajectory'),
    'provenance': ('.provenance', None),
    'testing': ('.testing', None),
}


def __getattr__(name: str) -> object:
    """Resolve a deferred attribute on first use."""
    entry = _LAZY.get(name)
    if entry is None:
        msg = f'module {__name__!r} has no attribute {name!r}'
        raise AttributeError(msg)
    module_path, attribute = entry
    module = importlib.import_module(module_path, __name__)
    value = module if attribute is None else getattr(module, attribute)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    """Include the lazily-resolved helpers in ``dir(seqcraft)``."""
    return sorted({*globals(), *_LAZY})


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
    'Raster',
    'RasterError',
    'Report',
    'ReportFailed',
    'SeqCraftError',
    'System',
    'UnitSanityError',
    'UnknownFieldError',
    'WriteResult',
    '__version__',
    'barrier',
    'bit_reversed_order',
    'centric_order',
    'compile',
    'compile_sequence',
    'convert',
    'display',
    'events',
    'flatten',
    'golden_angle',
    'interleaved_slice_order',
    'linear_order',
    'load_hardware',
    'modules',
    'ordering',
    'plot_block',
    'plot_kspace',
    'plot_sequence',
    'plot_trajectory',
    'provenance',
    'rf_spoil_phase',
    'span',
    'synthetic_hardware',
    'testing',
    'timing',
    'units',
    'validate',
]
