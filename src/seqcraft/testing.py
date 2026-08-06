"""
Assertions for testing your own modules.

Import these into your own test suite and a module you write outside seqcraft gets the same
contract checks the built-in ones get::

    import seqcraft as sc

    def test_my_module():
        sc.testing.assert_all(MyModule(system, ...))

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

from typing import Any

from .core import events as ev
from .core.compiler import compile_sequence
from .core.logic import LogicBlock
from .core.module import Module

__all__ = [
    'all_modules',
    'assert_all',
    'assert_compiles',
    'assert_deterministic',
    'assert_duration_is_honest',
    'assert_pure',
    'assert_raster',
    'assert_timing_properties_in_range',
    'assert_within_limits',
]


def _event_hashes(module: Module) -> dict[str, str]:
    return {
        key: ev.content_hash(value)
        for key, value in vars(module).items()
        if getattr(value, 'type', None) is not None
    }


def assert_pure(module: Module, **build_args: Any) -> None:
    """
    Assert that ``build`` mutates neither the module nor the events stored on it.

    ``build`` is called once per TR, so a module that mutates itself makes TR 500 differ from TR 1
    -- and the difference is usually a sign flip that produces a plausible but wrong image.

    Raises
    ------
    AssertionError
        Naming the attribute that changed.

    Examples
    --------
    >>> import seqcraft as sc
    >>> pe = sc.modules.PhaseEncode(sc.System.preset('generic_3t'), fov_pe_mm=250, matrix_pe=64)
    >>> sc.testing.assert_pure(pe, line=17)
    """
    before = _event_hashes(module)
    module.build(**build_args)
    module.build(**build_args)
    after = _event_hashes(module)
    changed = [key for key in before if before[key] != after.get(key)]
    assert not changed, (
        f'{type(module).__name__}.build() mutated {changed}; derive modified events with '
        f'seqcraft.events.derive() instead of assigning to them'
    )


def assert_deterministic(module: Module, **build_args: Any) -> None:
    """
    Assert that two builds with the same arguments produce the same block.

    Examples
    --------
    >>> import seqcraft as sc
    >>> spoil = sc.modules.Spoiler(sc.System.preset('generic_3t'), twists=4, voxel_mm=5)
    >>> sc.testing.assert_deterministic(spoil)
    """
    first, second = module.build(**build_args), module.build(**build_args)
    assert len(first) == len(second), (
        f'{type(module).__name__}.build() produced {len(first)} then {len(second)} nodes'
    )
    for index, (a, b) in enumerate(zip(first.nodes, second.nodes)):
        assert abs(a.start - b.start) < 1e-12, f'node {index} moved from {a.start} to {b.start}'
        if getattr(a.item, 'type', None) is not None:
            assert ev.content_hash(a.item) == ev.content_hash(b.item), f'node {index} changed'


def assert_duration_is_honest(module: Module, **build_args: Any) -> None:
    """
    Assert that ``module.duration`` is at least as long as the block it builds.

    Callers place the *next* thing by this number, so a duration that falls short silently overlaps
    whatever comes after.  Exceeding it is fine -- that is padding to the block raster.

    Examples
    --------
    >>> import seqcraft as sc
    >>> ro = sc.modules.CartesianLine(sc.System.preset('generic_3t'), fov_ro_mm=250,
    ...                               matrix_ro=64, readout_duration_us=3200)
    >>> sc.testing.assert_duration_is_honest(ro)
    """
    declared = getattr(module, 'duration', None)
    if declared is None:
        return
    measured = module.build(**build_args).duration
    assert declared >= measured - 1e-12, (
        f'{type(module).__name__}.duration is {declared * 1e6:.1f} us but its block is '
        f'{measured * 1e6:.1f} us -- whatever is placed next would overlap it'
    )


def assert_timing_properties_in_range(module: Module, **build_args: Any) -> None:
    """
    Assert that ``isodelay`` and ``time_to_echo``, where present, lie inside the block.

    They are offsets into the block a caller places by, so one outside it is always a bug -- and a
    quiet one, since the sequence still compiles and just puts the echo somewhere else.
    """
    duration = getattr(module, 'duration', None)
    if not duration:
        return
    for name in ('isodelay', 'time_to_echo'):
        value = getattr(module, name, None)
        if value is not None:
            assert 0.0 <= value <= duration + 1e-12, (
                f'{type(module).__name__}.{name} = {value * 1e6:.1f} us is outside the block '
                f'(0 .. {duration * 1e6:.1f} us)'
            )


def assert_within_limits(module: Module, **build_args: Any) -> None:
    """
    Assert that the built block respects the amplitude and slew limits of its own regime.

    Per-axis only: the vector norm across simultaneous axes routinely exceeds the per-axis limit and
    is legal on real amplifiers, which is why the compiler reports it as a warning.

    Examples
    --------
    >>> import seqcraft as sc
    >>> spoil = sc.modules.Spoiler(sc.System.preset('generic_3t'), twists=4, voxel_mm=5)
    >>> sc.testing.assert_within_limits(spoil)
    """
    block = module.build(**build_args)
    events = [n.item for n in block if getattr(n.item, 'type', None) in ('trap', 'grad')]
    violations = [
        entry for entry in ev.check_limits(events, module.opts, module.system.grad_raster.dt)
        if not entry[0].endswith('_norm')
    ]
    assert not violations, f'{type(module).__name__}: {violations}'


def assert_raster(module: Module, **build_args: Any) -> None:
    """Assert that every node start lands on the gradient raster."""
    raster = module.system.grad_raster
    for index, node in enumerate(module.build(**build_args)):
        assert raster.holds(node.start), (
            f'{type(module).__name__} node {index} starts at {node.start * 1e6:.4f} us, '
            f'off the {raster.dt * 1e6:.0f} us raster'
        )


def assert_compiles(module: Module, **build_args: Any) -> None:
    """
    Assert that the module produces a legal pulseq sequence on its own.

    A module that only works when something else happens to be beside it is not reusable, so this
    compiles it alone and requires a clean report.

    Examples
    --------
    >>> import seqcraft as sc
    >>> exc = sc.modules.SincExcitation(sc.System.preset('generic_3t'), flip_deg=15,
    ...                                 duration_us=1000, slice_thickness_mm=5)
    >>> sc.testing.assert_compiles(exc)
    """
    block = module.build(**build_args)
    if not block.nodes or block.duration == 0.0:
        return
    name = type(module).__name__
    out = compile_sequence(
        LogicBlock(f'assert_{name}').add(0.0, block), module.system, regime=module.regime
    )
    errors = [
        issue for issue in out.check().errors
        if issue.kind != 'timing' or 'TotalDuration' not in issue.message
    ]
    assert not errors, f'{name} does not compile cleanly on its own: {errors}'


def assert_all(module: Module, **build_args: Any) -> None:
    """
    Run every assertion in this module against `module`.

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
    assert_pure(module, **build_args)
    assert_deterministic(module, **build_args)
    assert_duration_is_honest(module, **build_args)
    assert_timing_properties_in_range(module, **build_args)
    assert_within_limits(module, **build_args)
    assert_raster(module, **build_args)
    assert_compiles(module, **build_args)


def all_modules() -> dict[str, type[Module]]:
    """
    Return every concrete :class:`~seqcraft.Module` subclass, keyed by class name.

    Parametrise your own contract suite over this and a module gains the whole suite the moment it
    is written.  Subclassing *is* the registration: there is no decorator to forget, which is what
    a registry could not guarantee -- a module that forgot ``@register()`` silently lost its
    coverage, the failure the registry existed to prevent.

    Only classes reachable by import are found, so import your package first.  ``seqcraft.modules``
    is imported by ``import seqcraft``, so the built-ins are always present.

    Examples
    --------
    >>> import seqcraft as sc
    >>> found = sc.testing.all_modules()
    >>> 'SpiralVDS' in found and 'MonopolarDiffusion' in found
    True
    >>> issubclass(found['SincExcitation'], sc.Module)
    True
    """
    out: dict[str, type[Module]] = {}

    def walk(cls: type[Module]) -> None:
        for sub in cls.__subclasses__():
            if not getattr(sub, '__abstractmethods__', None):
                out[sub.__name__] = sub
            walk(sub)

    walk(Module)
    return dict(sorted(out.items()))
