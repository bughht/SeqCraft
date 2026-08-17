"""
``spoiler`` -- four lines of arithmetic, and the one argument that has no safe default.

Testable in four lines, which is part of the argument for it being a function.
"""

from __future__ import annotations

import pytest

import seqcraft as sc


@pytest.mark.parametrize(('cycles', 'voxel_mm'), [(4.0, 5.0), (8.0, 3.0), (2.5, 1.25)])
def test_the_area_is_cycles_over_the_voxel(opts, cycles, voxel_mm) -> None:
    """``A = cycles_per_voxel / voxel_m``, in 1/m.  The whole module."""
    block = sc.modules.spoiler(opts, cycles_per_voxel=cycles, voxel_mm=voxel_mm)

    assert float(block.nodes[0].item.area) == pytest.approx(cycles / (voxel_mm / 1e3))


def test_the_axis_is_where_it_was_asked_for(opts) -> None:
    assert sc.modules.spoiler(opts, voxel_mm=5.0).nodes[0].item.channel == 'z'
    assert sc.modules.spoiler(opts, voxel_mm=5.0, axis='x').nodes[0].item.channel == 'x'


def test_voxel_mm_has_no_default(opts) -> None:
    """
    Keyword-only and required, because the wrong value is silent.

    For ``axis='z'`` the voxel dimension is the *slice thickness*; a default would invite the
    in-plane size, which under-spoils by the ratio of the two and shows up as faint residual
    banding that reads as anything except a spoiler bug.
    """
    with pytest.raises(TypeError, match='voxel_mm'):
        sc.modules.spoiler(opts)


@pytest.mark.parametrize(('kwargs', 'match'), [
    ({'voxel_mm': 0.0}, 'voxel_mm'),
    ({'voxel_mm': -5.0}, 'voxel_mm'),
    ({'voxel_mm': 5.0, 'cycles_per_voxel': 0.0}, 'cycles_per_voxel'),
    ({'voxel_mm': 5.0, 'axis': 'ky'}, 'axis'),
])
def test_the_refusals_name_the_argument(opts, kwargs, match) -> None:
    with pytest.raises(sc.ConfigurationError, match=match):
        sc.modules.spoiler(opts, **kwargs)


def test_it_compiles_on_its_own(opts, component_checks) -> None:
    """A component that only works beside something else is not reusable."""
    component_checks.output(lambda: sc.modules.spoiler(opts, voxel_mm=5.0), opts)


def test_it_is_not_a_module_subclass() -> None:
    """
    Stated here as well as in the layout test, because this is where a reader will look.

    The top of ``modules/`` is for what is not an ``sc.Module``, and this scores none of the
    four things a class buys: it designs nothing per call, holds no state, answers no timing
    question, and has no variant to derive.
    """
    assert not isinstance(sc.modules.spoiler, type)
