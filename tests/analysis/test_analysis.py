"""
The analysis toolbox: four functions, one entry shape, and which of them are exact.

The distinction between :func:`sc.moments` and :func:`sc.sample` is the whole reason this file
exists.  They now sit in one module and it will look tempting to build the first on the second;
:func:`test_moments_and_sample_diverge_on_a_spiral` is what makes that a failing test rather than
a docstring nobody reads.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.design.events import trapz


def _sampled_area(tree, opts, axis: str) -> float:
    """m0 the lossy way: trapezoidal integration of the uniform-grid samples."""
    grid, grads, _ = sc.sample(tree, opts)
    return float(trapz(grads[axis], grid))


def _raster_centre_lobe(opts, n: int = 13, scale: float = 0.1):
    """
    An arbitrary gradient: samples at raster *centres*, curvature over a few of them.

    The shape that separates the two measurements.  ``make_arbitrary_grad`` puts its samples at
    ``(k + 0.5) * raster``, so not one of them lands on the uniform grid ``sample`` interpolates
    onto -- and with only thirteen samples across the lobe, halving the step across the peak
    visibly rounds it.  Thirteen is close to the sharpest such lobe the slew limit permits; making
    it smoother is what makes the two measurements converge again.
    """
    w = np.sin(np.pi * np.linspace(0.0, 1.0, n)) * scale * opts.max_grad
    g = pp.make_arbitrary_grad('x', waveform=w, first=0.0, last=0.0, system=opts)
    return sc.LogicBlock('lobe').add(0.0, g), g


# ------------------------------------------------------------------------------------ moments
def test_moments_returns_the_analytic_area_of_a_trapezoid(opts) -> None:
    """A trapezoid's area is the number that was asked for, to the last bit pypulseq kept."""
    tree = sc.LogicBlock('spoil').add(0.0, pp.make_trapezoid('z', area=500.0, system=opts))
    assert sc.moments(tree)['z'] == pytest.approx(500.0, rel=1e-9)


def test_moments_sums_over_the_whole_tree(opts) -> None:
    """Nested blocks are flattened, or a component that nests would be measured vacuously."""
    inner = sc.LogicBlock('inner').add(0.0, pp.make_trapezoid('x', area=100.0, system=opts))
    tree = (
        sc.LogicBlock('tr')
        .add(0.0, inner)
        .add(2e-3, pp.make_trapezoid('x', area=-40.0, system=opts))
        .add(2e-3, pp.make_trapezoid('y', area=250.0, system=opts))
    )
    assert sc.moments(tree) == pytest.approx({'x': 60.0, 'y': 250.0}, abs=1e-9)


def test_moments_takes_no_opts_and_does_no_compile(opts) -> None:
    """
    A moment is a property of the waveform, not of the scanner that plays it.

    Stated as a test because the alternative signature -- ``moments(tree, opts)`` -- is the one
    every other function here has, so the difference is worth pinning.
    """
    tree = sc.LogicBlock('t').add(5e-6, pp.make_trapezoid('x', area=100.0, system=opts))
    # 5 us is off the gradient raster, so this tree cannot be compiled at all.
    with pytest.raises(sc.CompileError):
        sc.compile(tree, opts)
    assert sc.moments(tree)['x'] == pytest.approx(100.0, rel=1e-9)


def test_m1_sees_a_displacement_that_m0_cannot(opts) -> None:
    """``m1`` shifts by ``area * dt``; ``m0`` is exactly the quantity a time shift preserves."""
    g = pp.make_trapezoid('x', area=100.0, duration=1e-3, system=opts)
    here = sc.LogicBlock('t').add(0.0, g)
    later = sc.LogicBlock('t').add(1e-3, g)

    assert sc.moments(here)['x'] == pytest.approx(sc.moments(later)['x'], rel=1e-12)
    assert sc.moments(later, 1)['x'] - sc.moments(here, 1)['x'] == pytest.approx(
        100.0 * 1e-3, rel=1e-9
    )


# ------------------------------------------------------------------- exact against approximate
def test_moments_and_sample_agree_on_a_trapezoid(opts) -> None:
    """
    A trapezoid's knots are on the raster, so interpolation has nothing to lose.

    The first half of the distinction: ``sample`` is not *always* wrong, which is why the failure
    it can cause is easy to miss.
    """
    tree = sc.LogicBlock('t').add(0.0, pp.make_trapezoid('x', area=500.0, system=opts))
    exact = sc.moments(tree)['x']
    sampled = _sampled_area(tree, opts, 'x')
    assert abs(sampled - exact) / abs(exact) < 0.05


def test_sample_rounds_the_peak_off_a_raster_centre_waveform(opts) -> None:
    """
    The reason ``moments`` must never be built on ``sample`` -- and the reason m0 cannot say so.

    An arbitrary gradient's samples sit at raster *centres*, so ``sample`` interpolates them onto
    raster *edges* and loses amplitude across a peak: percent-level, and visible in the picture.

    What it does **not** lose is area, and that is the trap.  Linear interpolation errs
    antisymmetrically about each knot, so the two halves cancel under the integral and m0 comes
    back bit-identical.  That is exactly why the compiler's m0 invariant could not see a resampled
    spiral, and it is why a ``moments`` built on ``sample`` would look correct on every test
    anyone would think to write.

    Asserted rather than documented, because the two now live in one file and the shortcut would
    look like a simplification.
    """
    tree, event = _raster_centre_lobe(opts)
    grid, grads, _ = sc.sample(tree, opts)
    peak_exact = float(np.max(np.abs(sc.events.knots_of(event, 0.0)[1])))
    peak_sampled = float(np.max(np.abs(grads['x'])))

    lost = (peak_exact - peak_sampled) / peak_exact
    assert lost > 0.005, (
        f'sample() must visibly round the peak off a raster-centre waveform; it lost {lost:.3%}. '
        f'If this stops being true the exactness argument needs re-checking, not the threshold '
        f'lowering.'
    )

    # ... and the area survives it exactly, which is the whole trap.
    assert float(trapz(grads['x'], grid)) == pytest.approx(
        sc.moments(tree)['x'], rel=1e-9
    ), 'm0 is precisely the quantity that cannot detect this, so it must agree here'


def test_sample_finds_the_rf_and_adc_spans(opts) -> None:
    """The other half of what ``sample`` returns, and what ``plot_block`` shades."""
    rf = pp.make_sinc_pulse(flip_angle=0.5, duration=1e-3, system=opts, use='excitation',
                            slice_thickness=5e-3, apodization=0.5, time_bw_product=4)
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = (
        sc.LogicBlock('tr')
        .add(0.0, rf)
        .add(3e-3, adc)
        .add(5e-3, sc.barrier('done'))
    )
    kinds = [kind for kind, *_ in sc.sample(tree, opts)[2]]
    assert kinds == ['rf', 'adc', 'barrier']


# ------------------------------------------------------------------------------ kspace and pns
def test_kspace_takes_a_tree_and_names_its_return(opts) -> None:
    """
    ``calculate_kspacePP`` returns its tuple in a different order from ``calculate_kspace``.

    Getting that wrong swaps the trajectory for its timebase -- a wrong answer with no error --
    which is the whole reason this wrapper exists.
    """
    gx = pp.make_trapezoid('x', flat_area=128.0, flat_time=1.28e-3, system=opts)
    adc = pp.make_adc(num_samples=64, duration=1.28e-3, delay=gx.rise_time, system=opts)
    rf = pp.make_sinc_pulse(flip_angle=0.5, duration=1e-3, system=opts, use='excitation',
                            slice_thickness=5e-3, apodization=0.5, time_bw_product=4)
    tree = sc.LogicBlock('t').add(0.0, rf).add(3e-3, gx, adc)

    k = sc.kspace(tree, opts)
    assert set(k) == {'k_adc', 't_adc', 'k', 't_k', 't_excitation', 't_refocusing'}
    assert k['k_adc'].shape[0] == 3, 'k_adc is (3, n_samples), not the timebase'
    assert k['t_adc'].ndim == 1
    assert k['k_adc'].shape[1] == k['t_adc'].size


def test_pns_returns_the_whole_curve_not_just_a_verdict(opts) -> None:
    """
    When ``ok`` is False, ``peak`` says how much but ``norm`` and ``t`` say *where*.

    Which is what you need to fix it, so returning ``(ok, peak)`` was the wrong shape.
    """
    tree = sc.LogicBlock('t').add(
        0.0, pp.make_trapezoid('x', area=500.0, duration=2e-3, system=opts))
    r = sc.pns(tree, opts, sc.hardware.synthetic_hardware())

    assert set(r) == {'ok', 'peak', 'norm', 'components', 't'}
    assert isinstance(r['ok'], bool)
    assert 0.0 <= r['peak']
    assert r['norm'].size == r['t'].size


def test_the_synthetic_hardware_says_it_is_not_a_real_scanner() -> None:
    """
    The one safety note in the package, and it used to live in a docstring that was deleted.

    It travels with the model now -- in ``is_synthetic`` and in the ``repr`` -- rather than in the
    docstring of whatever happens to consume it.
    """
    hw = sc.hardware.synthetic_hardware()
    assert hw.is_synthetic is True
    assert 'NOT a real scanner' in repr(hw)
    assert 'human scan' in repr(hw)
