"""
The ``Module`` contract: four members, and what each of them guarantees.

The base is small on purpose, so these tests are mostly about what it *does not* do -- it does not
inspect a returned block beyond its type, does not overwrite a tag its subclass set, does not
declare a duration, and does not reach for a global scanner when one is missing.

The previous version of this file asserted the opposite philosophy: that the base declares no
abstract method and that ``core`` holds the module abstraction.  Both are now false, which is why
it was rewritten rather than adjusted.
"""

from __future__ import annotations

import pypulseq as pp
import pytest

import seqcraft as sc


class Blip(sc.Module):
    """A phase-encode blip: designed once at its largest, scaled per line."""

    def __init__(self, *, opts, fov_mm=250.0, matrix=64, axis='y', tag=None):
        super().__init__(opts=opts, tag=tag)
        self.dk = 1e3 / float(fov_mm)
        self.g = pp.make_trapezoid(channel=axis, area=self.dk * matrix / 2, system=opts)

    def build(self, *, line: int = 0) -> sc.LogicBlock:
        scale = line * self.dk / float(self.g.area)
        return sc.LogicBlock().add(0.0, sc.events.derive(
            self.g,
            amplitude=float(self.g.amplitude) * scale,
            area=float(self.g.area) * scale,
            flat_area=float(self.g.flat_area) * scale,
        ))


class Pair(sc.Module):
    """A module that holds a module -- composition with no registration API."""

    def __init__(self, *, opts, gap_s=500e-6, tag=None):
        super().__init__(opts=opts, tag=tag)
        self.inner = Blip(opts=opts)
        self.gap_s = gap_s

    def build(self, *, line=0) -> sc.LogicBlock:
        first = self.inner(line=line)
        return (sc.LogicBlock()
                .add(0.0, first)
                .add(first.duration + self.gap_s, self.inner(line=-line)))


# ------------------------------------------------------------------------- call and build
def test_calling_returns_exactly_what_build_returned(opts) -> None:
    """``__call__`` is a checkpoint, not a transformation: same object, same nodes."""
    module = Blip(opts=opts)
    made = module.build(line=7)
    nodes_before = list(made.nodes)

    returned = module._finalize(made)

    assert returned is made
    assert list(returned.nodes) == nodes_before


def test_call_and_build_agree_on_the_block(opts) -> None:
    module = Blip(opts=opts)

    assert module(line=7).duration == module.build(line=7).duration


def test_build_arguments_reach_build(opts) -> None:
    """``*args, **kwargs`` on the abstract signature is what lets a subclass name its own."""
    module = Blip(opts=opts)

    assert module(line=0).nodes[0].item.amplitude == 0.0
    assert module(line=8).nodes[0].item.amplitude != 0.0


# -------------------------------------------------------------------------------- tagging
def test_an_untagged_block_is_named_after_the_class(opts) -> None:
    assert Blip(opts=opts)().tag == 'Blip'


def test_an_instance_tag_wins_over_the_class_name(opts) -> None:
    """Two instances of one class doing different jobs are told apart by this."""
    assert Blip(opts=opts, tag='ky_pos')().tag == 'ky_pos'


def test_a_tag_set_inside_build_is_never_overridden(opts) -> None:
    """``_finalize`` names an unnamed block; it does not rename a named one."""

    class Deliberate(sc.Module):
        def build(self) -> sc.LogicBlock:
            return sc.LogicBlock('chosen_by_build')

    assert Deliberate(opts=opts, tag='ignored')().tag == 'chosen_by_build'


# ----------------------------------------------------------------------- the two failures
def test_a_subclass_without_build_cannot_be_instantiated(opts) -> None:
    """
    A subclass that does not produce a block is not a module.

    ``abc`` catches it at construction rather than at the call, so the traceback points at the
    class rather than at whatever loop was building a sequence out of it.
    """

    class Broken(sc.Module):
        pass

    with pytest.raises(TypeError, match='abstract'):
        Broken(opts=opts)


def test_build_returning_a_non_block_names_the_class(opts) -> None:
    """
    Without this the failure surfaces hundreds of lines later, inside ``add``.

    Putting the subclass's own name in the message is the whole value: ``add`` can only say that
    it was handed a ``SimpleNamespace``, not which module handed it over.
    """

    class Wrong(sc.Module):
        def build(self):
            return pp.make_delay(1e-3)

    with pytest.raises(TypeError, match=r'Wrong\.build\(\) returned SimpleNamespace'):
        Wrong(opts=opts)()


# -------------------------------------------------------------------------------- opts
def test_opts_is_required(opts) -> None:
    """
    Not defaulted, because pypulseq's fallback is the *process-global* ``Opts.default``.

    A module that silently designed against that would make the sequence depend on import order,
    which is the failure the whole scanner design is arranged to prevent.
    """
    with pytest.raises(TypeError, match='opts'):
        Blip()


def test_opts_is_kept_as_given(opts) -> None:
    """No wrapping and no copying: the object handed to pypulseq is the one that was passed."""
    assert Blip(opts=opts).opts is opts


def test_a_submodule_gets_the_same_opts(opts) -> None:
    assert Pair(opts=opts).inner.opts is opts


# ------------------------------------------------------------------------------ duration
def test_a_build_argument_may_change_the_duration(opts) -> None:
    """
    The rule that was inverted: a declared duration had to forbid this, and now nothing does.

    The block measures itself, so a call that produces a longer block simply reports a longer
    duration -- there is no second number to disagree with it.
    """

    class Variable(sc.Module):
        def build(self, *, duration_s: float) -> sc.LogicBlock:
            return sc.LogicBlock().add(0.0, pp.make_delay(duration_s))

    module = Variable(opts=opts)

    assert module(duration_s=1e-3).duration == pytest.approx(1e-3)
    assert module(duration_s=5e-3).duration == pytest.approx(5e-3)


def test_the_module_declares_no_duration(opts) -> None:
    """A module-level ``duration`` would be a second source of truth.  There is not one."""
    assert not hasattr(Blip(opts=opts), 'duration')


def test_a_constant_duration_across_lines_is_the_modules_own_doing(opts) -> None:
    """
    Designing at the largest area and scaling down is what keeps every line the same length.

    The base does not arrange this and could not; it is stated here because the caller's
    placement arithmetic depends on it and nothing else checks it.
    """
    module = Blip(opts=opts)

    assert module(line=1).duration == module(line=-32).duration


# ---------------------------------------------------------------------------- provenance
def test_nesting_produces_the_provenance_path(opts) -> None:
    """
    The payoff of auto-tagging: not one tag string is written by hand, and every compiled block
    still traces back to the module that produced it.
    """
    tree = sc.LogicBlock('tr').add(0.0, Pair(opts=opts)(line=10))

    out = sc.compile(tree, opts)

    assert out.origin(0) == ('tr', 'Pair', 'Blip')


# ------------------------------------------------------------------------ the layer boundary
def test_core_does_not_import_the_module_contract() -> None:
    """
    The compile path must stay independent of the component contract, in that direction.

    ``core`` takes a ``LogicBlock`` and never asks what produced it.  Importing ``seqcraft.module``
    from anywhere inside ``core`` would make that false by construction, and the dependency would
    be invisible until something tried to reuse the compiler on its own.
    """
    import importlib
    import pkgutil
    import sys

    import seqcraft.core

    for info in pkgutil.walk_packages(seqcraft.core.__path__, 'seqcraft.core.'):
        importlib.import_module(info.name)

    offenders = sorted(
        name for name, mod in sys.modules.items()
        if name.startswith('seqcraft.core')
        and mod is not None
        and any(
            getattr(value, '__module__', '') == 'seqcraft.module'
            or value is sys.modules.get('seqcraft.module')
            for value in vars(mod).values()
        )
    )
    assert not offenders, f'seqcraft.core reaches seqcraft.module from: {offenders}'


def test_a_plain_function_is_still_a_component(opts) -> None:
    """``Module`` is the standard shape for a reusable component, not a gate."""

    def crusher(opts, *, area_per_m=400.0) -> sc.LogicBlock:
        return sc.LogicBlock('crush').add(
            0.0, pp.make_trapezoid('z', area=area_per_m, system=opts))

    out = sc.compile(sc.LogicBlock('tr').add(0.0, crusher(opts)), opts)

    assert out.check().ok
    assert out.origin(0) == ('tr', 'crush')
