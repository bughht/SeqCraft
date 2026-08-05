"""
Shared fixtures.

One rule the whole suite follows: **nothing here mocks pypulseq.**  Every test builds real events
and compiles a real sequence, because the failures worth catching are physics failures, and a mock
cannot have wrong physics.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

import seqcraft as sc

if TYPE_CHECKING:
    from pypulseq.opts import Opts


@pytest.fixture(scope='session')
def system() -> sc.System:
    """A generic 3 T scanner: 40 mT/m, 150 T/m/s."""
    return sc.System.preset('generic_3t')


@pytest.fixture(scope='session')
def derated() -> sc.System:
    """A 3 T scanner with a derated ``'quiet'`` regime, for testing named regimes."""
    return sc.System.preset('generic_3t').derate('quiet', grad=0.8, slew=0.5)


@pytest.fixture(scope='session')
def opts(system: sc.System) -> Opts:
    """The default regime's pypulseq ``Opts``."""
    return system.default


@pytest.fixture(scope='session')
def geometry() -> sc.Geometry:
    """A small 2D geometry: 64 x 64 over 250 mm, one 5 mm slice."""
    return sc.Geometry(fov_mm=(250.0, 250.0, 5.0), matrix=(64, 64, 1), slice_thickness_mm=5.0)
