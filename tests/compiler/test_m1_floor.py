"""
The m1 tolerance has a floor, because the relative term alone sat inside the noise.

An arbitrary gradient the compiler has to split and rebuild does not round-trip exactly.  The
residual is proportional to **peak |g|** and to nothing else -- not to knot count, not to duration,
not to ``max_grad``.  The relative term scales as traversed area times the horizon, which grows
faster, so the two cross: short or weak waveforms were refused while long strong ones passed, for
the same physics and with nothing wrong in either.

Before the floor, six of ten legal sign-changing gradients were refused, and which six could not
be predicted from the sequence -- ``n=2000`` at 90 % passed with a residual of 7e-15 while
``n=4000`` at 90 % failed with 1.2e-7.

Nothing here builds a module.  The property under test is *one arbitrary gradient that changes
sign inside itself*, which is a compiler concern that happens to be met first by spirals.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
import seqcraft.compiler as compiler_package
from seqcraft.compiler.verification import _sequence_moments
from seqcraft.design.events import knots_of, pwl_moment

#: Measured residual per unit peak gradient, constant to within 5 % over 2000-8000 knots,
#: 20 % - 90 % of two different ``max_grad`` values.  The floor coefficient clears it by ~120.
MEASURED_FLOOR_PER_PEAK = 8.4e-14


def _opts(max_grad: float = 40.0) -> pp.Opts:
    return pp.Opts(max_grad=max_grad, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)


def _sign_changing(n: int, fraction: float, opts: pp.Opts):
    """One arbitrary gradient that runs out, crosses zero several times, and comes home."""
    r = np.linspace(0.0, 1.0, n)
    amp = fraction * opts.max_grad * r * np.sin(16.0 * n / 4000 * np.pi * r)
    return pp.make_arbitrary_grad(channel='x', waveform=amp, system=opts, delay=0.0,
                                  first=0.0, last=0.0)


def _split_tree(grad, n: int, opts: pp.Opts) -> sc.LogicBlock:
    """The same gradient, cut across a block boundary so the compiler has to rebuild it."""
    duration = n * opts.grad_raster_time
    return sc.LogicBlock('sign_changing').add(0.0, grad).add(duration / 2, sc.barrier())


@pytest.mark.parametrize('n', [1000, 2000, 4000, 8000, 16000])
@pytest.mark.parametrize('fraction', [0.3, 0.9])
def test_a_legal_sign_changing_gradient_compiles(n: int, fraction: float) -> None:
    """
    The regression, stated as the thing that was broken: these are all legal, and all were not.

    Six of these ten combinations were refused before the floor existed.
    """
    opts = _opts()
    sc.compile(_split_tree(_sign_changing(n, fraction, opts), n, opts), opts, name='legal')


def test_the_residual_tracks_peak_amplitude_and_nothing_else() -> None:
    """
    Why the floor is proportional to peak |g| rather than a constant or a bigger coefficient.

    This is the measurement the tolerance is derived from, so it is asserted rather than
    described: if the residual ever stops tracking peak, the floor's form is wrong and this test
    is the one that should say so.
    """
    original = compiler_package.verify_against_tree
    compiler_package.verify_against_tree = lambda *a, **k: None       # the check is under test
    try:
        ratios = []
        for max_grad in (20.0, 40.0):
            opts = _opts(max_grad)
            for n in (2000, 4000, 8000):
                for fraction in (0.2, 0.5):
                    grad = _sign_changing(n, fraction, opts)
                    knots = knots_of(grad, 0.0)
                    seq = sc.compile(_split_tree(grad, n, opts), opts, name='m')
                    residual = abs(_sequence_moments(seq, 1)['x'] - pwl_moment(*knots, 1))
                    ratios.append(residual / float(np.max(np.abs(knots[1]))))
    finally:
        compiler_package.verify_against_tree = original

    assert np.allclose(ratios, MEASURED_FLOOR_PER_PEAK, rtol=0.08), (
        f'residual/peak spread {min(ratios):.3g}..{max(ratios):.3g}, expected '
        f'~{MEASURED_FLOOR_PER_PEAK:.3g}'
    )


def test_the_floor_clears_the_residual_by_two_orders() -> None:
    """The lower margin: the tolerance must sit well above what the representation costs."""
    opts = _opts()
    grad = _sign_changing(4000, 0.9, opts)
    peak = float(np.max(np.abs(knots_of(grad, 0.0)[1])))
    assert 1e-11 * peak / (MEASURED_FLOOR_PER_PEAK * peak) > 100.0


def test_the_floor_stays_far_below_what_the_check_exists_to_catch() -> None:
    """
    The upper margin, and the one that matters: a floor that swallowed the signal would be worse
    than the bug.

    m1 exists to catch a gradient played at the wrong time.  One raster of displacement moves m1
    by the area traversed times the raster, and that must stay orders of magnitude above the
    tolerance.
    """
    opts = _opts()
    for n, fraction in ((1000, 0.3), (16000, 0.9)):
        grad = _sign_changing(n, fraction, opts)
        times, amps = knots_of(grad, 0.0)
        peak = float(np.max(np.abs(amps)))
        traversed = float(np.trapezoid(np.abs(amps), times))
        signal = traversed * opts.grad_raster_time
        assert signal / (1e-11 * peak) > 500.0, f'n={n} fraction={fraction}'


def test_a_displaced_gradient_is_still_caught() -> None:
    """
    The floor did not blind the check.

    A trapezoid moved by one raster is exactly the failure m1 is for, and it must still be
    refused -- with the floor in place, on an axis whose peak is large.
    """
    opts = _opts()
    trap = pp.make_trapezoid(channel='x', area=5000.0, duration=4e-3, system=opts)
    knots = knots_of(trap, 0.0)
    peak = float(np.max(np.abs(knots[1])))
    shifted = abs(pwl_moment(*knots_of(trap, opts.grad_raster_time), 1) - pwl_moment(*knots, 1))
    assert shifted > 1e-11 * peak * 100.0, (
        f'one raster moves m1 by {shifted:.4g}, floor is {1e-11 * peak:.4g}'
    )
