"""
Assertions for testing your own sequence components.

What seqcraft can check about a component is what it *produces*: a
:class:`~seqcraft.design.logic.LogicBlock` that is reproducible, respects the amplifier, and compiles
to a legal ``.seq`` on its own.  Nothing here asks what produced it, so a plain function and a class
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

What is deliberately *not* here is any check the compiler already makes with a better message.
A malformed node cannot be constructed -- :meth:`~seqcraft.design.logic.LogicBlock.add` rejects
anything that is neither an event nor a block -- and an off-raster gradient start makes
:func:`~seqcraft.compiler.compile_sequence` raise, naming the nearest raster point above and
below and two ways to fix it.  A second, thinner assertion for either would only be a worse error
message competing with the good one.

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

from typing import TYPE_CHECKING, Any

from .compiler import compile_sequence
from .design import events as ev
from .design.logic import LogicBlock, flatten

if TYPE_CHECKING:
    from collections.abc import Callable

    from pypulseq.opts import Opts

__all__ = [
    'assert_all',
    'assert_deterministic',
    'assert_output',
    'assert_pure',
    'assert_within_limits',
]


def _name(block: LogicBlock) -> str:
    """Return a block's tag, or a stand-in, for error messages."""
    return block.tag or 'untagged'


# ------------------------------------------------------------------ what any component produces
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


def assert_output(make: Callable[[], LogicBlock], opts: Opts) -> None:
    """
    Run every block-level assertion against whatever `make` returns.

    The universal one: it takes a callable, so it applies to a module call, to one of several
    methods on a class of your own, or to a bare function.  Nothing here inspects the caller.

    Three checks, and each is here because nothing else performs it:

    **Determinism** -- two calls agree across the whole tree, by content hash.

    **Per-axis limits** -- against the ``Opts`` the component was *designed* against, with node
    times, which the compiler's own limit check cannot supply because it sees one block at a time.

    **It compiles alone** -- a component that only works when something else happens to be beside
    it is not reusable.  The compile is what checks the raster: an off-raster gradient start makes
    :func:`~seqcraft.compiler.compile_sequence` raise ``CompileError`` naming the nearest
    raster point above and below, which is strictly more than a separate assertion could say.

    Parameters
    ----------
    make
        A callable of no arguments returning a :class:`~seqcraft.design.logic.LogicBlock`.
    opts
        The scanner to check the limits against, and to compile with.

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
    assert_deterministic(make)
    assert_within_limits(block, opts)
    if not block.nodes or block.duration == 0.0:
        return
    out = compile_sequence(LogicBlock(f'assert_{_name(block)}').add(0.0, block), opts)
    errors = [
        issue for issue in out.check().errors
        if issue.kind != 'timing' or 'TotalDuration' not in issue.message
    ]
    assert not errors, f'{_name(block)} does not compile cleanly on its own: {errors}'


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
