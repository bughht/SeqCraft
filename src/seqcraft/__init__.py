"""
seqcraft -- composable, verifiable MRI pulse sequence programming on top of pypulseq.

Three things, and one of them is pypulseq's.

:class:`~seqcraft.design.logic.LogicBlock`
    A tree of pulseq events and nested blocks, each with a start time.  Two attributes and one
    method.  Anything may overlap anything.
:func:`~seqcraft.compiler.compile_sequence`
    Turns a tree into legal pulseq blocks: finds block boundaries, sums gradients that share an
    axis, and validates the result against the amplifier.
:class:`pypulseq.Opts`
    The scanner.  Not wrapped, not subclassed, not hidden behind a seqcraft class -- the same
    object you pass to ``pp.make_trapezoid`` is the one you pass to ``sc.compile``.

:class:`~seqcraft.design.module.Module` is the standard shape for a *reusable* component that produces
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
>>> seq = sc.compile(tr, opts)                      # a pypulseq.Sequence, or an exception
>>> len(seq.block_events), round(seq.duration()[0] * 1e3, 3)
(4, 8.24)

**The central contract.**  ``sc.compile(tree, opts)`` returns a :class:`pypulseq.Sequence` and
nothing else.  If the tree cannot become a legal sequence it **raises**; if the compiler had to
change a waveform to make it legal it **warns**.  There is no report object to inspect and no
result wrapper to unpack.

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
library code, and changing your own sequence should never mean editing a package.  What
:mod:`seqcraft.modules` ships instead is a small set of MR building blocks that each leave the
sequence choices with the caller -- ``GRE2D`` takes the list of phase-encode lines to acquire
rather than an acceleration factor -- and every one of them was extracted from a working sequence
rather than designed.

Layout
------
Named for what each layer answers, in the order a sequence passes through them::

    scanner/     what you build against    Opts, and the PNS response model
    design/      what you build            the tree, events, timing, units
    modules/     what you build with       the MR building blocks, each one extracted
    compiler/    the transform             boundaries, legalization, emission, verification
    analysis     measuring a tree          sample, moments, kspace, pns
    display      looking at a tree         plot_block

and :mod:`~seqcraft.errors` beside them, which everything may raise.

The dependencies run one way -- ``errors -> design -> compiler -> analysis -> display`` -- and two
tests assert it, so nothing on the compile path can come to import the display helpers or the
scanner package.

This module is the only global re-export layer.  Import from here.

Notes
-----
Importing this package is side-effect free: it does not touch ``pypulseq.Opts.default`` and prints
nothing.  ``display`` and ``testing`` are resolved on first access rather than at import, so
neither the plotting stack nor the assertions are paid for unless used -- though note that
``import pypulseq`` itself imports matplotlib (in ``Sequence/calc_grad_spectrum.py``), so the
saving is seqcraft's own weight rather than the whole plotting stack's.  ``analysis`` is eager: it
imports numpy and the compiler, both of which are already loaded.
"""

from __future__ import annotations

import importlib

from . import _compat, modules
from ._version import __version__
from .analysis import kspace, moments, pns, sample
from .compiler import compile_sequence as compile  # noqa: A001, A004
from .compiler.errors import CompileError, DefinitionConflict, HardwareLimitError
from .compiler.verification import CompilerContractError
from .design import events, timing, units
from .design.logic import Item, LogicBlock, Node, barrier, flatten, span
from .design.module import Module
from .design.timing import Raster, RasterError
from .design.units import convert
from .errors import (
    ConfigurationError,
    MissingExtraError,
    SeqCraftError,
    SeqCraftWarning,
)
from .scanner import hardware, opts
from .scanner.opts import UnknownFieldError

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
#: is used.  It is the only entry left: ``testing`` used to be here too, and the assertions it
#: held are in the test suite now.
_LAZY: dict[str, tuple[str, str | None]] = {
    'display': ('.display', None),
    'plot_block': ('.display', 'plot_block'),
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
    'CompilerContractError',
    'ConfigurationError',
    'DefinitionConflict',
    'HardwareLimitError',
    'Item',
    'LogicBlock',
    'MissingExtraError',
    'Module',
    'Node',
    'Raster',
    'RasterError',
    'SeqCraftError',
    'SeqCraftWarning',
    'UnknownFieldError',
    '__version__',
    'analysis',
    'barrier',
    'compile',
    'compile_sequence',
    'convert',
    'display',
    'events',
    'flatten',
    'hardware',
    'kspace',
    'moments',
    'modules',
    'opts',
    'plot_block',
    'pns',
    'sample',
    'scanner',
    'span',
    'timing',
    'units',
]
