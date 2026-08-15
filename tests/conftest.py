"""
Shared fixtures.

One rule the whole suite follows: **nothing here mocks pypulseq.**  Every test builds real events
and compiles a real sequence, because the failures worth catching are physics failures, and a mock
cannot have wrong physics.

A second rule, newer: **no fixture builds a sequence out of library modules.**  seqcraft ships no
concrete modules, and the coupling that used to exist -- compiler tests whose realistic trees came
from the module library -- is exactly what made the compiler impossible to test independently of
whatever the library happened to contain.  Trees here are raw pypulseq events.
"""

from __future__ import annotations

import pytest
from pypulseq.opts import Opts


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
    import seqcraft as sc

    return sc.opts.derate(opts, grad=0.8, slew=0.5)


@pytest.fixture(scope='session')
def geometry():
    """A small 2D geometry: 64 x 64 over 250 mm, one 5 mm slice."""
    import seqcraft as sc

    return sc.Geometry(fov_mm=(250.0, 250.0, 5.0), matrix=(64, 64, 1), slice_thickness_mm=5.0)
