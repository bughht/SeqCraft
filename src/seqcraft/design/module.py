"""
:class:`Module` -- the standard shape for a reusable sequence component.

A module is **an MR building block that takes parameters and returns one**
:class:`~seqcraft.design.logic.LogicBlock`.  That is the entire idea, and everything here follows
from keeping it that small.

Three layers, three kinds of knowledge, no overlap::

   pypulseq Opts + parameters          Module                LogicBlock              pulseq
   ────────────────────────────  ─►  designs once   ─►   tree of events   ─►   sc.compile()
   what the scanner can do            answers when       relative times        block boundaries
   what this block should be          the echo is        overlap is free       gradient legality

+--------------+------------------------------------------+---------------------------+
| Layer        | Owns                                     | Knows nothing about       |
+==============+==========================================+===========================+
| ``Module``   | **MR physics** -- a b-value is these two | where it will be placed   |
|              | lobes; k = 0 arrives 1.92 ms in          |                           |
+--------------+------------------------------------------+---------------------------+
| ``LogicBlock``| **time** -- a child starts 8.71 ms into | what an echo is           |
|              | its parent                               |                           |
+--------------+------------------------------------------+---------------------------+
| compiler     | **pulseq legality** -- boundaries,       | what produced the tree    |
|              | splits, sums, limits                     |                           |
+--------------+------------------------------------------+---------------------------+

The conventions, and why each is one
------------------------------------
**One call returns one block** -- including when the parts are far apart in time.  A diffusion
encoding is one block with a hole in the middle, and the caller drops the refocusing pulse into
the hole.  Overlap in the tree costs nothing, so nothing has to be reserved and no artificial
duration is added.

**``__init__`` designs, ``build`` assembles.**  Waveforms are created once and stored on ``self``;
a call reads them and derives variants.  Thirty diffusion directions cost one design and thirty
cheap calls.

**Calls are pure.**  A call must not mutate the module or the events on it.  Derive a modified
event with :func:`seqcraft.events.derive`; never assign to ``self.g.amplitude``.  The classic bug
-- ``self.gx.amplitude = -self.gx.amplitude`` in a readout loop -- still compiles, and makes TR 500
differ from TR 1.

**The module declares no duration; the block measures itself.**  :attr:`LogicBlock.duration` is
measured from the nodes, so a module-level ``duration`` would be a second source of truth that can
disagree with the first -- and a build argument the property cannot see is exactly how it comes to
disagree.  Build the block, then read it::

    exc_block = exc()
    tr.add(0.0, exc_block)
    tr.add(exc_block.duration + gap, readout(line=k))

Building is cheap -- the design happened in ``__init__`` -- so nothing is lost by building before
placing.  What is gained is that **a call argument may now change the duration freely**, which a
declared duration had to forbid.

**Semantic offsets are a different question, and they are the module's.**  Where k = 0 falls
inside a readout, or where an RF's effective centre is, cannot be measured from a tree of events:
the block knows times, not meanings.  Those belong to whatever module knows the physics, expressed
however that module's domain wants -- ``ro.time_to_echo(shot=1)``, taking the same arguments as
``build`` if a build argument moves the echo.  The base standardises no names for them.

**Composition is an attribute.**  A module that holds modules just holds them, passes ``opts``
down, and calls them in ``build``.  There is no registration API, and nesting produces the
provenance path for free.

**Nothing else belongs in the base.**  No provenance walker, no unit checker, no registry, no
submodule bookkeeping, no scanner wrapper.  Each is a real feature; none is required by *every*
module, and a base class that carries what only some subclasses need is how a small idea becomes a
framework.

**And the base is not a requirement of the system.**  The compiler takes a ``LogicBlock`` and never
asks what produced it, so a plain function is still a valid component::

    def blip(opts, *, area_per_m, axis='y'):
        g = pp.make_trapezoid(channel=axis, area=area_per_m, system=opts)
        return sc.LogicBlock('blip').add(0.0, g)

``Module`` is the standard shape for a *reusable* component, not a gate.

Examples
--------
>>> import pypulseq as pp
>>> import seqcraft as sc
>>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
>>>
>>> class PhaseEncode(sc.Module):
...     '''A phase-encode blip, designed once at its largest and scaled per line.'''
...
...     def __init__(self, *, opts, fov_mm, matrix, axis='y', tag=None):
...         super().__init__(opts=opts, tag=tag)
...         self.dk = 1e3 / float(fov_mm)
...         self.g = pp.make_trapezoid(channel=axis, area=self.dk * int(matrix) / 2, system=opts)
...
...     def build(self, *, line: int = 0) -> sc.LogicBlock:
...         scale = line * self.dk / float(self.g.area)
...         g = sc.events.derive(self.g,
...                              amplitude=float(self.g.amplitude) * scale,
...                              area=float(self.g.area) * scale,
...                              flat_area=float(self.g.flat_area) * scale)
...         return sc.LogicBlock().add(0.0, g)      # untagged: _finalize names it
>>>
>>> pe = PhaseEncode(opts=opts, fov_mm=250, matrix=64)
>>> pe(line=17)
LogicBlock(PhaseEncode, 1 node, 0.30 ms)
>>> pe(line=-32)                       # a different line, deliberately the same duration
LogicBlock(PhaseEncode, 1 node, 0.30 ms)
>>> PhaseEncode(opts=opts, fov_mm=250, matrix=64, tag='ky')().tag
'ky'

Designing at the **largest** area and scaling down is what keeps every line the same length, so the
caller's placement arithmetic does not depend on which line is being acquired.

Composition needs no API, and auto-tagging pays for itself:

>>> class Pair(sc.Module):
...     def __init__(self, *, opts, gap_s=500e-6, tag=None):
...         super().__init__(opts=opts, tag=tag)
...         self.pe = PhaseEncode(opts=opts, fov_mm=250, matrix=64)   # opts passed down
...         self.gap_s = gap_s
...
...     def build(self, *, line=0) -> sc.LogicBlock:
...         first = self.pe(line=line)                   # place by what was just built
...         return (sc.LogicBlock()
...                 .add(0.0, first)
...                 .add(first.duration + self.gap_s, self.pe(line=-line)))
>>>
>>> out = sc.compile(sc.LogicBlock('tr').add(0.0, Pair(opts=opts)(line=10)), opts)
>>> out.origin(0)
('tr', 'Pair', 'PhaseEncode')

Not one tag string was written anywhere, and every compiled block still traces back to the module
that produced it.

The two failure modes the base catches:

>>> class Broken(sc.Module): pass
>>> Broken(opts=opts)
Traceback (most recent call last):
    ...
TypeError: Can't instantiate abstract class Broken with...abstract method...build...

>>> class Wrong(sc.Module):
...     def build(self): return pp.make_delay(1e-3)
>>> Wrong(opts=opts)()
Traceback (most recent call last):
    ...
TypeError: Wrong.build() returned SimpleNamespace, not a LogicBlock
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any

from ..design.logic import LogicBlock

if TYPE_CHECKING:
    from pypulseq.opts import Opts

__all__ = ['Module']


class Module(ABC):
    """
    A reusable MR building block: parameters in, one ``LogicBlock`` out.

    Parameters
    ----------
    opts
        The scanner, as a :class:`pypulseq.Opts`.  The only scanner input a module gets: rasters,
        dead times, ringdown, gradient and B1 limits, gamma, B0, sample limits.

        **Required rather than defaulted**, because pypulseq's fallback is the *process-global*
        ``Opts.default``, which makes a sequence depend on import order.  Pass it explicitly,
        always, including down to submodules.
    tag
        Optional identity.  Tags become the provenance path :func:`~seqcraft.flatten` builds and
        the ``from ...`` clause in every compiler warning, so two instances of one class doing
        different jobs -- two readouts, two spoilers -- are told apart by it.  Defaults to the
        class name.

    Attributes
    ----------
    opts : pypulseq.Opts
    tag : str or None

    Notes
    -----
    ``build`` is deliberately declared as ``build(*args, **kwargs)`` and narrowed by subclasses to
    named keywords (``line=``, ``direction=``, ``shot=``).  A type checker reports that as an
    override violation; the looseness is the point, since a module's arguments are its own domain,
    so the rule is silenced for module code rather than the signatures weakened.

    A component with more than one output -- a diffusion encoding's ``pre`` and ``post`` -- gives
    ``build`` the output that always exists and returns the rest from ordinary named methods.
    """

    def __init__(self, *, opts: Opts, tag: str | None = None) -> None:
        self.opts = opts
        self.tag = tag

    def __call__(self, *args: Any, **kwargs: Any) -> LogicBlock:
        """
        Build the block, check it, and name it.

        ``module(...)`` is the interface; ``build`` is what a subclass writes.  Two reasons for
        the split: it reads correctly beside ``add`` -- ``tr.add(t, readout(line=17))`` -- and it
        is the single place the framework gets to check and name what came back.
        """
        return self._finalize(self.build(*args, **kwargs))

    @abstractmethod
    def build(self, *args: Any, **kwargs: Any) -> LogicBlock:
        """
        Assemble and return this module's block.  **Write this; call** ``module(...)``.

        Abstract because a subclass that does not produce a block is not a module.
        """

    def _finalize(self, block: LogicBlock) -> LogicBlock:
        """
        Reject a wrong return type, and name an unnamed block.

        Exactly two things.  The type check puts *your* class in the traceback rather than
        producing a failure hundreds of lines later inside ``add``; the naming never overrides a
        tag ``build`` set deliberately.
        """
        if not isinstance(block, LogicBlock):
            msg = (
                f'{type(self).__name__}.build() returned {type(block).__name__}, '
                f'not a LogicBlock'
            )
            raise TypeError(msg)
        if not block.tag:
            block.tag = self.tag or type(self).__name__
        return block
