"""
The module base class: where a sequence's logic lives.

A module encapsulates one sequence task -- an excitation, a readout, a phase encode, a spoiler,
one lobe of a diffusion pair.  It works like a ``torch.nn.Module``: :meth:`Module.__init__` does
the expensive design once, and :meth:`Module.build` is cheap and returns a value.

Three conventions, and no machinery to enforce them
---------------------------------------------------
**``__init__`` designs, ``build`` assembles.**  Waveforms are created in ``__init__`` and stored
on ``self``, so a 30-direction diffusion encoding costs one design, not thirty.

**``build``'s arguments select the variant.**  A diffusion module returns its first or second
lobe according to ``part=``; an excitation applies a slice offset and an RF phase.  Nothing is
declared in advance -- they are ordinary keyword arguments with ordinary defaults, and Python
reports a typo as a ``TypeError``.

**Timing a caller needs in order to place the block is a property.**  ``exc.isodelay``,
``readout.time_to_echo``, ``diff.lobe_duration``.  The module has the domain knowledge, so the
module answers.  The block does not, and cannot: it does not know where it sits, and pinning the
answer into it would stop the same block being reused elsewhere.

Writing one
-----------
>>> import pypulseq as pp
>>> import seqcraft as sc
>>>
>>> class Blip(sc.Module):
...     '''A single phase-encode blip on one axis.'''
...     def __init__(self, system, *, area_per_m, axis='y'):
...         super().__init__(system)
...         self.axis = axis
...         self.area_per_m = area_per_m
...         self.g = pp.make_trapezoid(channel=axis, area=area_per_m, system=self.opts)
...
...     @property
...     def duration(self):
...         return float(pp.calc_duration(self.g))
...
...     def build(self, *, scale=1.0):
...         g = sc.events.derive(self.g, amplitude=self.g.amplitude * scale)
...         return sc.LogicBlock('blip').add(0.0, g)
>>>
>>> blip = Blip(sc.System.preset('generic_3t'), area_per_m=200.0)
>>> blip.build(scale=-1.0)
LogicBlock(blip, 1 node, 0.36 ms)
>>> round(blip.duration * 1e6)
360
"""

from __future__ import annotations

import abc
import functools
from typing import TYPE_CHECKING, Any

from .validate import check_units

if TYPE_CHECKING:
    from pypulseq.opts import Opts

    from .logic import LogicBlock
    from .system import System

__all__ = ['Module']


class Module(abc.ABC):
    """
    Base class for a sequence module.

    Parameters
    ----------
    system
        The scanner.  Held so :attr:`opts` can resolve the limits to design against, and so
        nested modules inherit it.
    regime
        Which of `system`'s named limit regimes this module designs against.  A diffusion
        encoding may run at full amplitude while an EPI readout is derated for peripheral nerve
        stimulation; those are two regimes of one system, and this is how a module picks one.
        The compiler independently validates the *combined* waveform.

    Attributes
    ----------
    system : System
    regime : str

    Notes
    -----
    Subclasses implement :meth:`build`, and nothing else is required.  There are no class
    variables to declare, no categories to register, no hooks to override: what a module *is*,
    is a thing that can build a logic block.

    The one thing that happens behind your back is unit validation.  After a subclass's
    ``__init__`` returns, :func:`seqcraft.validate.check_units` walks the attributes it set and
    rejects any numeric value implausible for the unit its name declares, so ``fov_mm=0.22``
    fails at construction with *"0.22 looks like metres, did you mean fov_mm=220?"* instead of
    producing a sequence with a 22 cm error in it.  Suffix a float attribute with its unit
    (``_mm``, ``_us``, ``_deg``, ``_per_m``) and the check applies; name it without one and it
    is skipped.
    """

    def __init_subclass__(cls, **kwargs: Any) -> None:
        """Wrap the subclass's own ``__init__`` so unit checking runs once it returns."""
        super().__init_subclass__(**kwargs)
        init = cls.__dict__.get('__init__')
        if init is None or getattr(init, '_seqcraft_checked', False):
            return

        @functools.wraps(init)
        def checked(self: Module, *args: Any, **kw: Any) -> None:
            init(self, *args, **kw)
            check_units(self)

        checked._seqcraft_checked = True  # type: ignore[attr-defined]
        cls.__init__ = checked  # type: ignore[method-assign]

    def __init__(self, system: System, *, regime: str = 'default') -> None:
        self.system = system
        self.regime = regime
        # Resolved now rather than lazily, so an unknown regime fails at construction with the
        # list of valid names, not at build time inside a loop over 1700 TRs.
        self._opts: Opts = system.limits(regime)

    @property
    def opts(self) -> Opts:
        """
        The pypulseq ``Opts`` this module designs against.

        Always pass this as ``system=`` to pypulseq's ``make_*`` functions.  Passing ``None``
        would silently fall back to the process-global ``Opts.default``, which makes a sequence
        depend on import order.
        """
        return self._opts

    @abc.abstractmethod
    def build(self, **args: Any) -> LogicBlock:
        """
        Return a :class:`~seqcraft.core.logic.LogicBlock` holding this module's events.

        Parameters
        ----------
        **args
            Variant selection: which lobe, which line, which slice, what phase.  Declare the
            ones the module supports as keyword arguments with defaults.

        Returns
        -------
        LogicBlock
            A fresh block each call.  Must not mutate the module or the events stored on it;
            derive modified events with :func:`seqcraft.events.derive`.
        """

    # -------------------------------------------------------------------------- reporting
    def params(self) -> dict[str, Any]:
        """
        Return this module's scalar parameters, for the provenance sidecar and ``repr``.

        Walks ``__dict__`` for JSON-safe values, so a subclass gets sensible provenance without
        declaring anything.  Pypulseq events and numpy arrays are skipped; nested modules are
        summarised.

        Examples
        --------
        >>> import seqcraft as sc
        >>> exc = sc.modules.SincExcitation(sc.System.preset('generic_3t'), flip_deg=15,
        ...                                 duration_us=1000, slice_thickness_mm=5)
        >>> exc.params()['flip_deg']
        15.0
        """
        out: dict[str, Any] = {}
        for key, value in vars(self).items():
            if key.startswith('_') or key in ('system', 'regime'):
                continue
            if isinstance(value, Module):
                out[key] = {'module': type(value).__name__, **value.params()}
            elif isinstance(value, (bool, int, float, str, type(None))):
                out[key] = value
            elif isinstance(value, (list, tuple)) and all(
                isinstance(v, (bool, int, float, str)) for v in value
            ):
                out[key] = list(value)
        return out

    def submodules(self) -> dict[str, Module]:
        """Return the modules held as attributes, keyed by attribute name."""
        return {k: v for k, v in vars(self).items() if isinstance(v, Module)}

    def __repr__(self) -> str:
        """``ClassName(key=value, ...)`` over the scalar parameters."""
        inner = ', '.join(f'{k}={v!r}' for k, v in self.params().items() if not isinstance(v, dict))
        return f'{type(self).__name__}({inner})'
