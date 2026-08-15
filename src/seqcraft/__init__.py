"""
seqcraft -- composable, verifiable MRI pulse sequence programming on top of pypulseq.

Three things, and one of them is pypulseq's.

:class:`~seqcraft.core.logic.LogicBlock`
    A tree of pulseq events and nested blocks, each with a start time.  Two attributes and one
    method.  Anything may overlap anything.
:func:`~seqcraft.core.compiler.compile_sequence`
    Turns a tree into legal pulseq blocks: finds block boundaries, sums gradients that share an
    axis, and validates the result against the amplifier.
:class:`pypulseq.Opts`
    The scanner.  Not wrapped, not subclassed, not hidden behind a seqcraft class -- the same
    object you pass to ``pp.make_trapezoid`` is the one you pass to ``sc.compile``.

:class:`~seqcraft.module.Module` is the standard shape for a *reusable* component that produces
blocks.  It is a convention, not a gate: the compiler never asks what produced a tree, so a plain
function is a component too.

Getting started
---------------
>>> import pypulseq as pp
>>> import seqcraft as sc
>>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
>>> rf, gz, gzr = pp.make_sinc_pulse(flip_angle=0.26, duration=1e-3, slice_thickness=5e-3,
...                                  delay=opts.rf_dead_time, use='excitation',
...                                  system=opts, return_gz=True)
>>> gx = pp.make_trapezoid('x', flat_area=256.0, flat_time=3.2e-3, system=opts)
>>> adc = pp.make_adc(num_samples=64, duration=3.2e-3, delay=gx.rise_time, system=opts)
>>>
>>> tr = sc.LogicBlock('demo')
>>> _ = tr.add(0.0, rf).add(0.0, gz)
>>> _ = tr.add(pp.calc_duration(gz), gzr)
>>> _ = tr.add(5e-3, gx).add(5e-3, adc)             # 5 ms after excitation
>>> out = sc.compile(tr, opts)
>>> out.check().ok
True

**Set the dead times.**  pypulseq defaults ``rf_dead_time``, ``rf_ringdown_time`` and
``adc_dead_time`` to zero, which is wrong on every real scanner: the sequence compiles cleanly,
validates cleanly, and is refused or silently mangled at the console.  They are properties of your
installation, not of the scanner model, so no preset and no vendor database can supply them.
:func:`seqcraft.scanner.opts.from_scanner` looks the amplitudes up in PulseqSystems and requires
the rest.

When to reach for this
----------------------
For one sequence, once, write raw pypulseq -- it is a fine tool for that.  seqcraft pays for
itself when you have a *family* of sequences, a loop over more than two axes, gradients that must
overlap without you hand-splitting blocks, or files that have to be reproducible six months
later.

There are no recipes here on purpose.  A recipe is somebody else's sequence choices baked into
library code, and changing your own sequence should never mean editing a package.  seqcraft ships
**no concrete modules at all**: what it ships is the tree, the compiler, and the contract.

Layout
------
``core/`` is the compile path -- the tree, the compiler, and the arithmetic they rest on --
and holds nothing else.  Everything beside it is deliberately *not* on that path:
:mod:`~seqcraft.scanner` describes the machine, :mod:`~seqcraft.module` is the component
contract, and :mod:`~seqcraft.display`, :mod:`~seqcraft.provenance` and :mod:`~seqcraft.testing`
are tools around it.

Notes
-----
Importing this package is cheap and side-effect free: it does not import matplotlib, does not
touch ``pypulseq.Opts.default``, and prints nothing.
"""

from __future__ import annotations

import importlib

from . import _compat
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
from .core.report import Issue, Report, ReportFailed
from .core.timing import Raster
from .core.units import convert
from .module import Module
from .scanner import hardware, opts

# Fail once with a complete list, rather than letting the first caller that needs a missing
# pypulseq function fail with an opaque AttributeError halfway through a build.
_compat.require()

#: Also available as ``sc.compile_sequence``, for callers who would rather not shadow the
#: builtin name in their own namespace.
compile_sequence = compile

#: Attributes resolved on first access rather than at import time.
#:
#: ``display`` is the only module allowed to import matplotlib, and ``import seqcraft`` must not
#: pull matplotlib in -- so ``sc.plot_block(block, opts)`` works without paying that cost until it
#: is used.  ``testing`` is deferred for tidiness: the contract assertions are documented as
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
    'LogicBlock',
    'MissingExtraError',
    'Module',
    'Node',
    'Raster',
    'RasterError',
    'Report',
    'ReportFailed',
    'SeqCraftError',
    'UnitSanityError',
    'UnknownFieldError',
    'WriteResult',
    '__version__',
    'barrier',
    'compile',
    'compile_sequence',
    'convert',
    'display',
    'events',
    'flatten',
    'hardware',
    'opts',
    'plot_block',
    'plot_kspace',
    'plot_sequence',
    'plot_trajectory',
    'provenance',
    'scanner',
    'span',
    'testing',
    'timing',
    'units',
    'validate',
]
