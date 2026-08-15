"""
Assertions for testing your own sequence components.

What seqcraft can check about a component is what it *produces*: a
:class:`~seqcraft.core.logic.LogicBlock` that lands on the raster, respects the amplifier, and
compiles to a legal ``.seq``.  Nothing here asks what produced it, so a plain function and a class
of your own design get exactly the same checks a :class:`seqcraft.Module` subclass gets::

    import seqcraft as sc

    def test_my_component():
        sc.testing.assert_output(lambda: my_encoding.pre(direction=(1, 0, 0)), opts)

:func:`assert_all` is the one to point at a module.  It adds the checks that only mean something
for the ``Module`` convention -- a call returns the block, waveforms live on ``self``, ``opts`` is
the scanner -- and reads them off the module itself::

    def test_my_module():
        sc.testing.assert_all(MyModule(opts=opts, ...), line=17)

They are ordinary functions raising ``AssertionError``, so they work with pytest, unittest, or a
bare script.

Examples
--------
>>> import pypulseq as pp
>>> import seqcraft as sc
>>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
>>>
>>> class Spoiler(sc.Module):
...     def __init__(self, *, opts, area_per_m, tag=None):
...         super().__init__(opts=opts, tag=tag)
...         self.g = pp.make_trapezoid('z', area=area_per_m, system=opts)
...
...     def build(self, *, scale: float = 1.0) -> sc.LogicBlock:
...         return sc.LogicBlock().add(0.0, sc.events.derive(
...             self.g, amplitude=float(self.g.amplitude) * scale,
...             area=float(self.g.area) * scale, flat_area=float(self.g.flat_area) * scale))
>>>
>>> sc.testing.assert_all(Spoiler(opts=opts, area_per_m=500.0), scale=-1.0)
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any

from .core import events as ev
from .core.compiler import compile_sequence
from .core.logic import LogicBlock, flatten
from .core.timing import Raster
from .module import Module

if TYPE_CHECKING:
    from collections.abc import Callable

    from pypulseq.opts import Opts

__all__ = [
    'assert_all',
    'assert_block',
    'assert_compiles',
    'assert_deterministic',
    'assert_output',
    'assert_pure',
    'assert_raster',
    'assert_within_limits',
    'module_subclasses',
]


def _name(block: LogicBlock) -> str:
    """Return a block's tag, or a stand-in, for error messages."""
    return block.tag or 'untagged'


# ------------------------------------------------------------------ what any component produces
def assert_block(block: object) -> None:
    """
    Assert that `block` is a well-formed :class:`~seqcraft.core.logic.LogicBlock`.

    This is the whole of seqcraft's contract with a component: whatever it hands the compiler is
    a logic block, and one whose children are events or nested blocks at finite times.
    ``nodes`` is a plain list you are invited to mutate, so the second half is worth checking.

    Raises
    ------
    AssertionError
        Naming what was returned, or which node is malformed.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s')
    >>> g = pp.make_trapezoid('x', area=100.0, system=opts)
    >>> sc.testing.assert_block(sc.LogicBlock('mine').add(0.0, g))
    """
    assert isinstance(block, LogicBlock), (
        f'expected a seqcraft.LogicBlock, got {type(block).__name__} -- a component takes part '
        f'in a sequence by returning one'
    )
    for index, node in enumerate(block.nodes):
        assert math.isfinite(node.start), f'{_name(block)} node {index} starts at {node.start}'
        assert isinstance(node.item, LogicBlock) or getattr(node.item, 'type', None) is not None, (
            f'{_name(block)} node {index} holds a {type(node.item).__name__}, which is neither a '
            f'pulseq event nor a LogicBlock'
        )
    assert math.isfinite(block.duration), f'{_name(block)} has a duration of {block.duration}'


def assert_deterministic(make: Callable[[], LogicBlock]) -> None:
    """
    Assert that two calls of `make` produce the same block.

    Parameters
    ----------
    make
        A callable of no arguments returning a block -- ``lambda: module(line=17)``,
        ``lambda: readout.prephaser()``, or a bare function reference.

    Notes
    -----
    Compares the whole tree, not the direct children: a component that nests would otherwise have
    its actual events skipped and pass vacuously.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s')
    >>> def spoiler():
    ...     return sc.LogicBlock('spoil').add(0.0, pp.make_trapezoid('z', area=500.0,
    ...                                                              system=opts))
    >>> sc.testing.assert_deterministic(spoiler)
    """
    first, second = list(flatten(make())), list(flatten(make()))
    assert len(first) == len(second), (
        f'two calls produced {len(first)} then {len(second)} events'
    )
    for index, ((t_a, a, path), (t_b, b, _)) in enumerate(zip(first, second, strict=True)):
        where = '.'.join(path) or '-'
        assert abs(t_a - t_b) < 1e-12, f'event {index} ({where}) moved from {t_a} to {t_b}'
        assert ev.content_hash(a) == ev.content_hash(b), f'event {index} ({where}) changed'


def assert_raster(block: LogicBlock, opts: Opts) -> None:
    """
    Assert that every gradient in `block` starts on the gradient raster of `opts`.

    Walks the whole tree rather than its direct children, so a nested component's own placement is
    checked too.  Only gradients: an RF or ADC event carries its own dead time in its ``delay`` and
    answers to the RF raster, which pypulseq's own timing check covers.
    """
    raster = Raster(float(opts.grad_raster_time), 'gradient')
    for index, (start, event, path) in enumerate(flatten(block)):
        if getattr(event, 'type', None) not in ('trap', 'grad'):
            continue
        assert raster.holds(start), (
            f'{_name(block)} gradient {index} ({".".join(path) or "-"}) starts at '
            f'{start * 1e6:.4f} us, off the {raster.dt * 1e6:.0f} us raster'
        )


def assert_within_limits(block: LogicBlock, opts: Opts) -> None:
    """
    Assert that `block` respects the amplitude and slew limits of `opts`.

    Per-axis only: the vector norm across simultaneous axes routinely exceeds the per-axis limit and
    is legal on real amplifiers, which is why the compiler reports it as a warning.

    Pass the ``Opts`` the component was *designed* against.  A part designed against derated limits
    is checked against those (:func:`seqcraft.scanner.opts.derate`); the un-derated ceiling is what
    the finished sequence is compiled and validated against.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s')
    >>> block = sc.LogicBlock('spoil').add(0.0, pp.make_trapezoid('z', area=500.0, system=opts))
    >>> sc.testing.assert_within_limits(block, opts)
    """
    # `flatten`, not iteration: a block's direct children may be nested blocks, and a component that
    # nests would otherwise have its actual gradients skipped entirely and pass this vacuously.
    placed = [
        (start, event) for start, event, _ in flatten(block)
        if getattr(event, 'type', None) in ('trap', 'grad')
    ]
    # Node times matter: two lobes of a bipolar pair are both on one axis, and without their
    # starts they would be taken to play simultaneously and sum to zero.
    violations = [
        entry for entry in ev.check_limits(
            [event for _, event in placed],
            opts,
            starts=[start for start, _ in placed],
        )
        if not entry[0].endswith('_norm')
    ]
    assert not violations, f'{_name(block)}: {violations}'


def assert_compiles(block: LogicBlock, opts: Opts) -> None:
    """
    Assert that `block` produces a legal pulseq sequence on its own.

    A component that only works when something else happens to be beside it is not reusable, so this
    compiles the block alone and requires a clean report.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
    ...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> block = sc.LogicBlock('spoil').add(0.0, pp.make_trapezoid('z', area=500.0, system=opts))
    >>> sc.testing.assert_compiles(block, opts)
    """
    if not block.nodes or block.duration == 0.0:
        return
    out = compile_sequence(LogicBlock(f'assert_{_name(block)}').add(0.0, block), opts)
    errors = [
        issue for issue in out.check().errors
        if issue.kind != 'timing' or 'TotalDuration' not in issue.message
    ]
    assert not errors, f'{_name(block)} does not compile cleanly on its own: {errors}'


def assert_output(make: Callable[[], LogicBlock], opts: Opts) -> None:
    """
    Run every block-level assertion against whatever `make` returns.

    The universal one: it takes a callable, so it applies to a module call, to one of several
    methods on a class of your own, or to a bare function.  Nothing here inspects the caller.

    Parameters
    ----------
    make
        A callable of no arguments returning a :class:`~seqcraft.core.logic.LogicBlock`.
    opts
        The scanner to check the raster and limits against, and to compile with.

    Examples
    --------
    An arbitrary function is a sequence component, and gets the whole suite:

    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
    ...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> def crusher(area_per_m=400.0):
    ...     g = pp.make_trapezoid('z', area=area_per_m, system=opts)
    ...     return sc.LogicBlock('crush').add(0.0, g)
    >>> sc.testing.assert_output(crusher, opts)
    """
    block = make()
    assert_block(block)
    assert_deterministic(make)
    assert_raster(block, opts)
    assert_within_limits(block, opts)
    assert_compiles(block, opts)


# ---------------------------------------------------------------- the conventions Module follows
def assert_pure(component: object, make: Callable[[], LogicBlock]) -> None:
    """
    Assert that calling `make` mutates neither `component` nor the events stored on it.

    For a component that designs once and assembles per TR, `make` is called once per TR -- so one
    that mutates itself makes TR 500 differ from TR 1, and the difference is usually a sign flip
    that produces a plausible but wrong image.

    Parameters
    ----------
    component
        The object whose stored pypulseq events must not change.
    make
        The call under test, e.g. ``lambda: pe(line=17)``.

    Raises
    ------
    AssertionError
        Naming the attribute that changed.

    Notes
    -----
    Checked after *each* of two calls, not once at the end.  The canonical bug this exists to
    catch -- the reference implementation's ``self.gx.amplitude = -self.gx.amplitude`` inside its
    readout loop -- is an involution, so comparing only before and against after the second call
    finds it back where it started and reports nothing.
    """
    before = _event_hashes(component)
    for call in (1, 2):
        make()
        after = _event_hashes(component)
        changed = [key for key in before if before[key] != after.get(key)]
        assert not changed, (
            f'{type(component).__name__} mutated {changed} on call {call}; derive modified events '
            f'with seqcraft.events.derive() instead of assigning to them'
        )


def assert_all(module: Any, **build_args: Any) -> None:
    """
    Run every assertion here against `module`, through calling it.

    The convenience wrapper for the :class:`seqcraft.Module` convention: ``module(**args)`` returns
    the block and ``module.opts`` says what to check it against.  Inheriting ``Module`` is not
    required -- any object with those two behaviours passes -- and a component shaped otherwise
    wants :func:`assert_output` instead.

    Parameters
    ----------
    module
        The module to check.
    **build_args
        Passed to every call, so a module whose default variant is degenerate can be checked in a
        meaningful configuration.

    Notes
    -----
    There is no duration check here, and deliberately so.  A module declares no duration -- the
    block it returns measures itself -- so there is no second number that can disagree with the
    first, and nothing left to assert.
    """
    def make() -> LogicBlock:
        block: LogicBlock = module(**build_args)
        return block

    assert_pure(module, make)
    assert_output(make, module.opts)


def _event_hashes(component: object) -> dict[str, str]:
    """Content-hash the pypulseq events a component stores as attributes."""
    return {
        key: ev.content_hash(value)
        for key, value in vars(component).items()
        if getattr(value, 'type', None) is not None
    }


# ------------------------------------------------------------------------------------- discovery
def module_subclasses() -> dict[str, type[Module]]:
    """
    Return every concrete subclass of :class:`seqcraft.Module`, keyed by class name.

    seqcraft ships no concrete modules, so this finds **your** library -- import your package
    first.  It is not the set of valid seqcraft components either: a function or a class that
    never heard of ``Module`` produces logic blocks just as well, and anything generic belongs in
    :func:`assert_output`, which asks nothing about ancestry.

    What this *is* for is parametrising a contract suite over the classes that do follow the
    convention, so one of them gains the whole suite the moment it is written.  Subclassing is the
    registration: there is no decorator to forget, which is what a registry could not guarantee --
    a module that omitted ``@register()`` silently lost its coverage, the failure the registry
    existed to prevent.

    Classes whose name begins with an underscore are private shared bases and are skipped, as are
    those left abstract by an ``abc.abstractmethod``.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> class Crusher(sc.Module):
    ...     def build(self):
    ...         return sc.LogicBlock().add(0.0, pp.make_trapezoid('z', area=400.0,
    ...                                                           system=self.opts))
    >>> 'Crusher' in sc.testing.module_subclasses()
    True
    """
    out: dict[str, type[Module]] = {}

    def walk(cls: type[Module]) -> None:
        for sub in cls.__subclasses__():
            if not sub.__name__.startswith('_') and not getattr(sub, '__abstractmethods__', None):
                out[sub.__name__] = sub
            walk(sub)

    walk(Module)
    return dict(sorted(out.items()))
