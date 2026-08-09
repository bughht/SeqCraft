"""
The architectural boundary: ``LogicBlock`` is the contract, ``Module`` is a convenience.

Every test here exists to hold one line open in seqcraft's design -- that a sequence component owes
the compiler a logic block and nothing else.  Three groups:

* ``core`` knows nothing about modules, and the compiler cannot be made to care;
* a function, and a class that inherits nothing, are first-class components;
* the optional base still earns its keep for the components that do want it.

The physics of the built-in library is ``test_module_physics.py``.  This file is about what
seqcraft does *not* require.
"""

from __future__ import annotations

import importlib
import inspect
import json
import re
from pathlib import Path

import pypulseq as pp
import pytest

import seqcraft as sc

CORE = Path(sc.core.__file__).parent


# --------------------------------------------------------------------------- the core boundary
def test_core_holds_no_module_abstraction() -> None:
    """
    ``core`` is what gets a logic block to a legal ``.seq``, and a module is not on that path.

    The old ``seqcraft.core.module`` is gone rather than shimmed: the package is 0.3.0 alpha and
    a compatibility import would keep asserting exactly the membership this branch removes.
    """
    assert not hasattr(sc.core, 'Module')
    assert not hasattr(sc, 'Module')
    assert not (CORE / 'module.py').exists()
    with pytest.raises(ImportError):
        importlib.import_module('seqcraft.core.module')


def test_core_never_imports_the_module_library() -> None:
    """
    The dependency runs one way: modules import core, and core imports nothing back.

    Asserted on the source rather than on ``sys.modules``, because importing any ``seqcraft.core``
    submodule pulls in the parent package and with it the library, so a runtime probe would prove
    nothing.
    """
    pattern = re.compile(
        r'^\s*(?:from\s+(?:\.\.|seqcraft\.)modules|import\s+seqcraft\.modules)', re.MULTILINE
    )
    offenders = [
        path.name for path in sorted(CORE.glob('*.py'))
        if pattern.search(path.read_text(encoding='utf-8'))
    ]
    assert not offenders, f'core imports the module library from: {offenders}'


def test_the_compiler_signature_mentions_only_blocks_and_systems() -> None:
    """``compile`` takes a tree and a scanner; what produced the tree is not its business."""
    parameters = inspect.signature(sc.compile).parameters
    assert list(parameters)[:2] == ['root', 'system']
    assert parameters['root'].annotation in (sc.LogicBlock, 'LogicBlock')


# -------------------------------------------------------------- components that inherit nothing
def _spoiler(system: sc.System, *, area_per_m: float = 800.0, axis: str = 'z') -> sc.LogicBlock:
    """A plain function is a sequence component: it returns a logic block."""
    g = pp.make_trapezoid(channel=axis, area=area_per_m, system=system.default)
    return sc.LogicBlock('spoil').add(0.0, g)


class _TwoSidedEncoding:
    """
    A gradient pair straddling a refocusing pulse, inheriting nothing at all.

    The shape the old ``build(**args) -> LogicBlock`` contract could not express well: two outputs
    with different meanings, named for what they are rather than distinguished by a keyword.
    """

    def __init__(self, system: sc.System, *, area_per_m: float, axis: str = 'z') -> None:
        self.system = system
        self.area_per_m = float(area_per_m)
        self.lobe = pp.make_trapezoid(channel=axis, area=self.area_per_m, system=system.default)

    @property
    def lobe_duration(self) -> float:
        return float(pp.calc_duration(self.lobe))

    def pre(self) -> sc.LogicBlock:
        return sc.LogicBlock('encode_pre').add(0.0, self.lobe)

    def post(self) -> sc.LogicBlock:
        return sc.LogicBlock('encode_post').add(0.0, self.lobe)


def test_a_plain_function_produces_a_block_that_compiles(system) -> None:
    block = _spoiler(system)
    assert isinstance(block, sc.LogicBlock)
    out = sc.compile(sc.LogicBlock('from_a_function').add(0.0, block), system)
    assert out.check().ok


def test_a_class_that_inherits_nothing_can_expose_several_block_methods(system) -> None:
    """``pre()`` and ``post()``, not ``build(part=...)`` -- and the compiler cannot tell."""
    encoding = _TwoSidedEncoding(system, area_per_m=600.0)
    refoc = sc.modules.HardRefocusing(system, flip_deg=180, duration_us=1000)

    seq = sc.LogicBlock('two_sided')
    seq.add(0.0, encoding.pre())
    seq.add(encoding.lobe_duration, refoc.build())
    seq.add(encoding.lobe_duration + refoc.duration, encoding.post())

    out = sc.compile(seq, system)
    assert out.check().ok
    assert isinstance(encoding.pre(), sc.LogicBlock)
    assert isinstance(encoding.post(), sc.LogicBlock)


def test_a_sequence_may_mix_a_module_a_class_and_a_function(system) -> None:
    """The three coexist because none of them is a category the compiler knows about."""
    exc = sc.modules.SincExcitation(
        system, flip_deg=15, duration_us=1000, slice_thickness_mm=5)
    encoding = _TwoSidedEncoding(system, area_per_m=400.0)

    seq = sc.LogicBlock('mixed')
    seq.add(0.0, exc.build())
    seq.add(exc.duration, encoding.pre())
    seq.add(exc.duration + encoding.lobe_duration, _spoiler(system))

    out = sc.compile(seq, system)
    assert out.check().ok


# --------------------------------------------------------------- the testing helpers are generic
def test_assert_output_asks_nothing_about_ancestry(system) -> None:
    """A bare function and each method of a base-less class get the whole block-level suite."""
    encoding = _TwoSidedEncoding(system, area_per_m=600.0)
    sc.testing.assert_output(lambda: _spoiler(system), system)
    sc.testing.assert_output(encoding.pre, system)
    sc.testing.assert_output(encoding.post, system)


def test_assert_output_rejects_something_that_is_not_a_block(system) -> None:
    with pytest.raises(AssertionError, match='LogicBlock'):
        sc.testing.assert_output(lambda: 'not a block', system)


def test_assert_all_needs_the_attributes_it_names_not_the_base(system) -> None:
    """``build``, ``system``, ``regime``, ``duration`` -- duck-typed, never isinstance-checked."""

    class Duck:
        def __init__(self) -> None:
            self.system = system
            self.regime = 'default'
            self.g = pp.make_trapezoid('z', area=400.0, system=system.default)

        @property
        def duration(self) -> float:
            return float(pp.calc_duration(self.g))

        def build(self) -> sc.LogicBlock:
            return sc.LogicBlock('duck').add(0.0, self.g)

    duck = Duck()
    assert not isinstance(duck, sc.modules.Module)
    sc.testing.assert_all(duck)


def test_assert_pure_still_catches_a_component_mutating_itself(system) -> None:
    """The check that matters most, and it never needed a base class to work."""

    class Leaky:
        def __init__(self) -> None:
            self.g = pp.make_trapezoid('x', area=200.0, system=system.default)

        def readout(self) -> sc.LogicBlock:
            self.g.amplitude = -self.g.amplitude          # the reference implementation's bug
            return sc.LogicBlock('leak').add(0.0, self.g)

    leaky = Leaky()
    with pytest.raises(AssertionError, match='mutated'):
        sc.testing.assert_pure(leaky, leaky.readout)


def test_module_subclasses_is_the_library_not_the_universe(system) -> None:
    """
    Discovery enumerates what inherits the base.  That is a smaller set than "components".

    It is the right set for parametrising the built-in contract suite and the wrong one for
    deciding what seqcraft accepts, which is why nothing outside the tests consults it.
    """
    found = sc.testing.module_subclasses()
    assert 'SpiralVDS' in found
    assert '_TwoSidedEncoding' not in found

    # ... and the component it cannot see compiles exactly like the ones it can.
    encoding = _TwoSidedEncoding(system, area_per_m=600.0)
    sc.testing.assert_output(encoding.pre, system)


def test_discovery_finds_a_subclass_of_your_own() -> None:
    """Subclassing is still the registration, for anyone who wants the contract suite."""

    class Mine(sc.modules.Module):
        def build(self) -> sc.LogicBlock:
            return _spoiler(self.system)

    found = sc.testing.module_subclasses()
    assert found.get('Mine') is Mine
    # And seqcraft's own coverage assertion filters by package, so it cannot be contaminated.
    assert Mine.__module__ != 'seqcraft.modules'


# ------------------------------------------------------------- what the optional base still does
class _Blip(sc.modules.Module):
    """A phase-encode blip, written on top of the base for what the base provides."""

    def __init__(
        self, system: sc.System, *, area_per_m: float, axis: str = 'y', regime: str = 'default'
    ) -> None:
        super().__init__(system, regime=regime)
        self.axis = axis
        self.area_per_m = float(area_per_m)
        self.g = pp.make_trapezoid(channel=axis, area=self.area_per_m, system=self.opts)

    @property
    def duration(self) -> float:
        return float(pp.calc_duration(self.g))

    def build(self, *, scale: float = 1.0) -> sc.LogicBlock:
        g = sc.events.derive(self.g, amplitude=float(self.g.amplitude) * scale)
        return sc.LogicBlock('blip').add(0.0, g)


def test_the_base_declares_no_abstract_method(system) -> None:
    """No ``build`` requirement, so a component may name its outputs for its own domain."""
    assert not getattr(sc.modules.Module, '__abstractmethods__', ())

    class NoBuild(sc.modules.Module):
        def __init__(self, sys_: sc.System) -> None:
            super().__init__(sys_)
            self.g = pp.make_trapezoid('x', area=300.0, system=self.opts)

        def prephaser(self) -> sc.LogicBlock:
            return sc.LogicBlock('pre').add(0.0, sc.events.derive(self.g, amplitude=-float(self.g.amplitude)))

        def readout(self) -> sc.LogicBlock:
            return sc.LogicBlock('ro').add(0.0, self.g)

    component = NoBuild(system)          # a TypeError under the old abstract contract
    assert not hasattr(component, 'build')
    sc.testing.assert_output(component.prephaser, system)
    sc.testing.assert_output(component.readout, system)


def test_the_base_resolves_the_limit_regime_at_construction(derated) -> None:
    quiet = _Blip(derated, area_per_m=200.0, regime='quiet')
    assert quiet.regime == 'quiet'
    assert float(quiet.opts.max_grad) < float(derated.default.max_grad)
    with pytest.raises(sc.ConfigurationError, match='unknown regime'):
        _Blip(derated, area_per_m=200.0, regime='typo')


def test_the_base_still_checks_units_when_a_subclass_init_returns(system) -> None:
    """The one thing that happens behind your back, and it survived the move out of core."""

    class Slab(sc.modules.Module):
        def __init__(self, sys_: sc.System) -> None:
            super().__init__(sys_)
            self.slice_thickness_mm = 0.005          # metres mistaken for millimetres

    with pytest.raises(sc.UnitSanityError, match='slice_thickness_mm'):
        Slab(system)


def test_check_units_takes_any_object_not_just_a_module() -> None:
    """A component that inherits nothing is not shut out of the unit check."""

    class Bare:
        def __init__(self) -> None:
            self.duration_us = 0.0042          # seconds mistaken for microseconds

    with pytest.raises(sc.UnitSanityError, match='duration_us'):
        sc.validate.check_units(Bare())


def test_the_base_reports_json_safe_params_and_a_readable_repr(system) -> None:
    blip = _Blip(system, area_per_m=200.0)
    params = blip.params()
    json.dumps(params)
    assert params['area_per_m'] == 200.0
    assert params['axis'] == 'y'
    assert not {'g', 'system', 'regime'} & set(params)
    assert repr(blip).startswith('_Blip(')


def test_the_base_reports_its_submodules(system) -> None:
    fat = sc.modules.FatSat(system, voxel_mm=5)
    assert set(fat.submodules()) == {'pulse', 'spoiler'}
    assert fat.params()['pulse']['module'] == 'GaussSaturation'


def test_provenance_takes_a_mapping_from_anywhere(system) -> None:
    """Params for the sidecar may come from the base, or from a dict you wrote by hand."""
    encoding = _TwoSidedEncoding(system, area_per_m=600.0)
    side = sc.provenance.build_sidecar({'area_per_m': encoding.area_per_m, 'axis': 'z'})
    assert side['resolved']['area_per_m'] == 600.0
