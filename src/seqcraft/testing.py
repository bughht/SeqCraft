"""
Assertions for testing your own sequence components.

What seqcraft can check about a component is what it *produces*: a
:class:`~seqcraft.core.logic.LogicBlock` that lands on the raster, respects the amplifier, and
compiles to a legal ``.seq``.  Nothing here asks what produced it, so a plain function and a class
of your own design get exactly the checks a built-in module gets::

    import seqcraft as sc

    def test_my_component():
        sc.testing.assert_output(lambda: my_encoding.pre(direction=(1, 0, 0)), system)

:func:`assert_all` adds the checks that only mean something for the convention the built-in
library follows -- a ``build()`` method, waveforms stored on ``self``, a ``duration`` property --
and is the one to point at a module::

    def test_my_module():
        sc.testing.assert_all(MyModule(system, ...))

Inheriting :class:`~seqcraft.modules.base.Module` is not required by either: ``assert_all`` calls
the attributes it names and skips the ones that are absent.

They are ordinary functions raising ``AssertionError``, so they work with pytest, unittest, or a
bare script.

Examples
--------
>>> import seqcraft as sc
>>> exc = sc.modules.SincExcitation(sc.System.preset('generic_3t'), flip_deg=15,
...                                 duration_us=1000, slice_thickness_mm=5)
>>> sc.testing.assert_all(exc)
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any

from .core import events as ev
from .core.compiler import compile_sequence
from .core.logic import LogicBlock, flatten
from .modules.base import Module

if TYPE_CHECKING:
    from collections.abc import Callable

    from .core.system import System

__all__ = [
    'assert_all',
    'assert_block',
    'assert_compiles',
    'assert_deterministic',
    'assert_duration_is_honest',
    'assert_output',
    'assert_pure',
    'assert_raster',
    'assert_timing_properties_in_range',
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
    >>> system = sc.System.preset('generic_3t')
    >>> g = pp.make_trapezoid('x', area=100.0, system=system.default)
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
        A callable of no arguments returning a block -- ``lambda: module.build(line=17)``,
        ``lambda: readout.prephaser()``, or a bare function reference.

    Notes
    -----
    Compares the whole tree, not the direct children: a component that nests would otherwise have
    its actual events skipped and pass vacuously.

    Examples
    --------
    >>> import seqcraft as sc
    >>> spoil = sc.modules.Spoiler(sc.System.preset('generic_3t'), twists=4, voxel_mm=5)
    >>> sc.testing.assert_deterministic(spoil.build)
    """
    first, second = list(flatten(make())), list(flatten(make()))
    assert len(first) == len(second), (
        f'two calls produced {len(first)} then {len(second)} events'
    )
    for index, ((t_a, a, path), (t_b, b, _)) in enumerate(zip(first, second, strict=True)):
        where = '.'.join(path) or '-'
        assert abs(t_a - t_b) < 1e-12, f'event {index} ({where}) moved from {t_a} to {t_b}'
        assert ev.content_hash(a) == ev.content_hash(b), f'event {index} ({where}) changed'


def assert_raster(block: LogicBlock, system: System) -> None:
    """
    Assert that every gradient in `block` starts on `system`'s gradient raster.

    Walks the whole tree rather than its direct children, so a nested component's own placement is
    checked too.  Only gradients: an RF or ADC event carries its own dead time in its ``delay`` and
    answers to the RF raster, which pypulseq's own timing check covers.
    """
    raster = system.grad_raster
    for index, (start, event, path) in enumerate(flatten(block)):
        if getattr(event, 'type', None) not in ('trap', 'grad'):
            continue
        assert raster.holds(start), (
            f'{_name(block)} gradient {index} ({".".join(path) or "-"}) starts at '
            f'{start * 1e6:.4f} us, off the {raster.dt * 1e6:.0f} us raster'
        )


def assert_within_limits(block: LogicBlock, system: System, *, regime: str = 'default') -> None:
    """
    Assert that `block` respects the amplitude and slew limits of `regime`.

    Per-axis only: the vector norm across simultaneous axes routinely exceeds the per-axis limit and
    is legal on real amplifiers, which is why the compiler reports it as a warning.

    Examples
    --------
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> spoil = sc.modules.Spoiler(system, twists=4, voxel_mm=5)
    >>> sc.testing.assert_within_limits(spoil.build(), system)
    """
    # `flatten`, not iteration: a block's direct children may be nested blocks, and a component that
    # nests -- FatSat, EPIReadout -- would otherwise have its actual gradients skipped entirely and
    # pass this vacuously.
    placed = [
        (start, event) for start, event, _ in flatten(block)
        if getattr(event, 'type', None) in ('trap', 'grad')
    ]
    # Node times matter: two lobes of a bipolar pair are both on one axis, and without their
    # starts they would be taken to play simultaneously and sum to zero.
    violations = [
        entry for entry in ev.check_limits(
            [event for _, event in placed],
            system.limits(regime),
            starts=[start for start, _ in placed],
        )
        if not entry[0].endswith('_norm')
    ]
    assert not violations, f'{_name(block)}: {violations}'


def assert_compiles(block: LogicBlock, system: System, *, regime: str = 'default') -> None:
    """
    Assert that `block` produces a legal pulseq sequence on its own.

    A component that only works when something else happens to be beside it is not reusable, so this
    compiles the block alone and requires a clean report.

    Examples
    --------
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> exc = sc.modules.SincExcitation(system, flip_deg=15, duration_us=1000,
    ...                                 slice_thickness_mm=5)
    >>> sc.testing.assert_compiles(exc.build(), system)
    """
    if not block.nodes or block.duration == 0.0:
        return
    out = compile_sequence(
        LogicBlock(f'assert_{_name(block)}').add(0.0, block), system, regime=regime
    )
    errors = [
        issue for issue in out.check().errors
        if issue.kind != 'timing' or 'TotalDuration' not in issue.message
    ]
    assert not errors, f'{_name(block)} does not compile cleanly on its own: {errors}'


def assert_output(
    make: Callable[[], LogicBlock], system: System, *, regime: str = 'default'
) -> None:
    """
    Run every block-level assertion against whatever `make` returns.

    The universal one: it takes a callable, so it applies to a module's ``build``, to one of
    several methods on a class of your own, or to a bare function.  Nothing here inspects the
    caller.

    Parameters
    ----------
    make
        A callable of no arguments returning a :class:`~seqcraft.core.logic.LogicBlock`.
    system
        The scanner to check the raster and limits against, and to compile with.
    regime
        Which of `system`'s named limit regimes the block was designed against.

    Examples
    --------
    An arbitrary function is a sequence component, and gets the whole suite:

    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> def crusher(area_per_m=400.0):
    ...     g = pp.make_trapezoid('z', area=area_per_m, system=system.default)
    ...     return sc.LogicBlock('crush').add(0.0, g)
    >>> sc.testing.assert_output(crusher, system)
    """
    block = make()
    assert_block(block)
    assert_deterministic(make)
    assert_raster(block, system)
    assert_within_limits(block, system, regime=regime)
    assert_compiles(block, system, regime=regime)


# --------------------------------------------------- the conventions the built-in library follows
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
        The call under test, e.g. ``lambda: pe.build(line=17)``.

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

    Examples
    --------
    >>> import seqcraft as sc
    >>> pe = sc.modules.PhaseEncode(sc.System.preset('generic_3t'), fov_pe_mm=250, matrix_pe=64)
    >>> sc.testing.assert_pure(pe, lambda: pe.build(line=17))
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


def assert_duration_is_honest(component: object, make: Callable[[], LogicBlock]) -> None:
    """
    Assert that ``component.duration`` is at least as long as the block `make` returns.

    Callers place the *next* thing by this number, so a duration that falls short silently overlaps
    whatever comes after.  Exceeding it is fine -- that is padding to the block raster.  A component
    that declares no ``duration`` is skipped: the property is a convenience, not a requirement.

    Examples
    --------
    >>> import seqcraft as sc
    >>> ro = sc.modules.CartesianLine(sc.System.preset('generic_3t'), fov_ro_mm=250,
    ...                               matrix_ro=64, readout_duration_us=3200)
    >>> sc.testing.assert_duration_is_honest(ro, ro.build)
    """
    declared = getattr(component, 'duration', None)
    if declared is None:
        return
    measured = make().duration
    assert declared >= measured - 1e-12, (
        f'{type(component).__name__}.duration is {declared * 1e6:.1f} us but its block is '
        f'{measured * 1e6:.1f} us -- whatever is placed next would overlap it'
    )


def assert_timing_properties_in_range(component: object) -> None:
    """
    Assert that ``isodelay`` and ``time_to_echo``, where present, lie inside ``duration``.

    They are offsets into the block a caller places by, so one outside it is always a bug -- and a
    quiet one, since the sequence still compiles and just puts the echo somewhere else.  Properties
    the component does not declare are skipped.

    Examples
    --------
    >>> import seqcraft as sc
    >>> exc = sc.modules.SincExcitation(sc.System.preset('generic_3t'), flip_deg=15,
    ...                                 duration_us=1000, slice_thickness_mm=5)
    >>> sc.testing.assert_timing_properties_in_range(exc)
    """
    duration = getattr(component, 'duration', None)
    if not duration:
        return
    for name in ('isodelay', 'time_to_echo'):
        value = getattr(component, name, None)
        if value is not None:
            assert 0.0 <= value <= duration + 1e-12, (
                f'{type(component).__name__}.{name} = {value * 1e6:.1f} us is outside the block '
                f'(0 .. {duration * 1e6:.1f} us)'
            )


def assert_all(module: Any, **build_args: Any) -> None:
    """
    Run every assertion here against `module`, through its ``build`` method.

    The convenience wrapper for the convention the built-in library follows: ``build(**args)``
    returns the block, ``system`` and ``regime`` say what to check it against, and ``duration``,
    ``isodelay`` and ``time_to_echo`` are properties where they apply.  Inheriting
    :class:`~seqcraft.modules.base.Module` is not required -- an object with those attributes
    passes -- and a component shaped otherwise wants :func:`assert_output` instead.

    Parameters
    ----------
    module
        The module to check.
    **build_args
        Passed to every ``build`` call, so a module whose default variant is degenerate can be
        checked in a meaningful configuration.

    Examples
    --------
    >>> import seqcraft as sc
    >>> diff = sc.modules.MonopolarDiffusion(sc.System.preset('generic_3t'),
    ...                                      b_value_s_per_mm2=1000, refocus_duration_us=4200)
    >>> sc.testing.assert_all(diff, part='pre', direction=(1, 0, 0))
    """
    def make() -> LogicBlock:
        block: LogicBlock = module.build(**build_args)
        return block

    assert_pure(module, make)
    assert_duration_is_honest(module, make)
    assert_timing_properties_in_range(module)
    assert_output(make, module.system, regime=getattr(module, 'regime', 'default'))


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
    Return every concrete subclass of :class:`~seqcraft.modules.base.Module`, keyed by class name.

    This is the built-in library, plus whatever of your own inherits the optional base.  It is
    **not** the set of valid seqcraft components: that base is a convenience, and a function or a
    class that never heard of it produces logic blocks just as well.  Anything generic belongs in
    :func:`assert_output`, which asks nothing about ancestry.

    What this *is* for is parametrising a contract suite over classes that do follow the
    convention, so one of them gains the whole suite the moment it is written.  Subclassing is the
    registration: there is no decorator to forget, which is what a registry could not guarantee --
    a module that omitted ``@register()`` silently lost its coverage, the failure the registry
    existed to prevent.

    Classes whose name begins with an underscore are private shared bases and are skipped, as are
    those left abstract by an ``abc.abstractmethod``.  Only classes reachable by import are found,
    so import your package first; ``seqcraft.modules`` is imported by ``import seqcraft``, so the
    built-ins are always present.

    Examples
    --------
    >>> import seqcraft as sc
    >>> found = sc.testing.module_subclasses()
    >>> 'SpiralVDS' in found and 'MonopolarDiffusion' in found
    True
    >>> issubclass(found['SincExcitation'], sc.modules.Module)
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
