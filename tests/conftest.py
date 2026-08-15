"""
Shared fixtures, and the two component assertions the package no longer ships.

One rule the whole suite follows: **nothing here mocks pypulseq.**  Every test builds real events
and compiles a real sequence, because the failures worth catching are physics failures, and a mock
cannot have wrong physics.

A second rule, newer: **no fixture builds a sequence out of library modules.**  seqcraft ships no
concrete modules, and the coupling that used to exist -- compiler tests whose realistic trees came
from the module library -- is exactly what made the compiler impossible to test independently of
whatever the library happened to contain.  Trees here are raw pypulseq events.

Why the assertions are here rather than in ``seqcraft.testing``
---------------------------------------------------------------
Because a package that ships assertions has to keep them working, and these are forty lines that
only this repository's own module tests use.  Everything else the module once asserted -- the
raster, the limits, that a block is well formed -- the compiler now checks, with a better message
and on the *summed* waveform rather than one component at a time.

What is left is exactly what the compiler **cannot** see, and the reason is the same for both:
it validates a tree, and it never sees the second call.  A component that designs once and
assembles per TR is called once per TR, so self-mutation, accumulation and nondeterminism are
invisible to ``sc.compile`` by construction -- each individual call hands it a tree that is
perfectly legal.  ``docs/writing_a_module.md`` teaches the three-line version, which is better
documentation than naming a function because it also explains what the check is *for*.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import pytest
from pypulseq.opts import Opts

import seqcraft as sc

if TYPE_CHECKING:
    from collections.abc import Callable


@pytest.fixture(scope='session')
def opts() -> Opts:
    """
    A generic 3 T scanner: 40 mT/m, 150 T/m/s, with realistic dead times.

    Built with the ordinary ``pp.Opts`` constructor, because that is what seqcraft asks a user to
    do -- there is no scanner class to build it from.  The dead times are stated rather than
    defaulted for the reason the docs give: pypulseq defaults all three to zero, and a sequence
    built on those compiles and validates cleanly before the console refuses it.
    """
    return Opts(
        max_grad=40, grad_unit='mT/m',
        max_slew=150, slew_unit='T/m/s',
        B0=3.0,
        rf_dead_time=100e-6,
        rf_ringdown_time=30e-6,
        adc_dead_time=10e-6,
    )


@pytest.fixture(scope='session')
def derated_opts(opts: Opts) -> Opts:
    """The same scanner derated, for tests that need a second, weaker limit set."""
    return sc.opts.derate(opts, grad=0.8, slew=0.5)


# ------------------------------------------------------ the two checks the compiler cannot make
def _event_hashes(component: object) -> dict[str, str]:
    """Content-hash the pypulseq events a component stores as attributes."""
    return {
        key: sc.events.content_hash(value)
        for key, value in vars(component).items()
        if getattr(value, 'type', None) is not None
    }


def assert_deterministic(make: Callable[[], sc.LogicBlock]) -> None:
    """
    Assert that two calls of `make` produce the same block.

    Compares the whole tree, not the direct children: a component that nests would otherwise have
    its actual events skipped and pass vacuously.
    """
    first, second = list(sc.flatten(make())), list(sc.flatten(make()))
    assert len(first) == len(second), (
        f'two calls produced {len(first)} then {len(second)} events'
    )
    for index, ((t_a, a, path), (t_b, b, _)) in enumerate(zip(first, second, strict=True)):
        where = '.'.join(path) or '-'
        assert abs(t_a - t_b) < 1e-12, f'event {index} ({where}) moved from {t_a} to {t_b}'
        assert sc.events.content_hash(a) == sc.events.content_hash(b), (
            f'event {index} ({where}) changed between calls'
        )


def assert_pure(component: object, make: Callable[[], sc.LogicBlock]) -> None:
    """
    Assert that calling `make` mutates neither `component` nor the events stored on it.

    Checked after *each* of two calls, not once at the end.  The canonical bug this exists to
    catch -- ``self.gx.amplitude = -self.gx.amplitude`` inside a readout loop -- is an involution,
    so comparing only before and after the second call finds it back where it started and reports
    nothing.  It compiles cleanly every TR and produces a plausible but wrong image.
    """
    before = _event_hashes(component)
    for call in (1, 2):
        make()
        after = _event_hashes(component)
        changed = [key for key in before if before[key] != after.get(key)]
        assert not changed, (
            f'{type(component).__name__} mutated {changed} on call {call}; derive modified '
            f'events with seqcraft.events.derive() instead of assigning to them'
        )


def assert_output(make: Callable[[], sc.LogicBlock], opts: Opts) -> None:
    """
    Determinism, and that whatever `make` returns compiles on its own.

    A component that only works when something else happens to be beside it is not reusable.  The
    compile is what checks everything else -- raster, limits, block legality -- and it raises with
    a better message than a separate assertion could give.
    """
    block = make()
    assert_deterministic(make)
    if not block.nodes or block.duration == 0.0:
        return
    sc.compile(sc.LogicBlock(f'assert_{block.tag or "untagged"}').add(0.0, block), opts)


def assert_all(module: Any, **build_args: Any) -> None:
    """
    Both checks against `module`, through calling it.

    The convenience wrapper for the :class:`seqcraft.Module` convention: ``module(**args)``
    returns the block and ``module.opts`` says what to compile it against.  Inheriting ``Module``
    is not required -- any object with those two behaviours passes.
    """
    def make() -> sc.LogicBlock:
        block: sc.LogicBlock = module(**build_args)
        return block

    assert_pure(module, make)
    assert_output(make, module.opts)


@pytest.fixture(scope='session')
def component_checks() -> Any:
    """
    The four assertions above, as a fixture, so a test does not import from ``conftest``.

    ``tests/integration`` already imports ``conftest`` directly for the sequence builders, but
    that file sits beside the tests that use it; this one is two directories up.
    """
    from types import SimpleNamespace

    return SimpleNamespace(
        deterministic=assert_deterministic,
        pure=assert_pure,
        output=assert_output,
        all=assert_all,
    )

