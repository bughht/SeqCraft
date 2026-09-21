r"""
The adapter between physical k and sigpy's coordinates, checked against physics rather than sigpy.

SigPy is the backend.  It is not evidence that the conversion *into* its convention is right, and
the conversion is one line -- which is why it needs this rather than why it does not.  A wrong
scale, a swapped axis, a flipped sign or a half-pixel origin all produce an image that looks like
an image, and a magnitude image hides every one of them.

So the oracle here is an **explicit dense DFT** written from

.. math:: y_j = \sum_r x(r) \exp(-2 \pi i\, k_j \cdot r)

with the sampling grid written out by hand.  It calls no sigpy.  The *tests* of course exercise
the sigpy-backed adapter and operator -- that is the comparison -- but nothing sigpy computes is
ever used as the reference.

Both real consumers drive the same core: ``RadialReadout`` and ``SpiralReadout``.  Nothing under
test asks which one it has.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'examples'))

sigpy = pytest.importorskip('sigpy', reason='needs seqcraft[recon]')
from noncartesian_recon import (  # noqa: E402
    NonCartesianAcquisition,
    dcf,
    encoding_operator,
    reconstruct,
    sigpy_coord,
)

#: Small on purpose.  A dense DFT is O(N^2 * samples), and every property under test is a property
#: of the convention rather than of the size.
MATRIX = 16
FOV_M = 0.24


@pytest.fixture
def opts() -> pp.Opts:
    return pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
                   adc_samples_limit=8192)


# --------------------------------------------------------------- the independent oracle
def image_grid(matrix: int = MATRIX, fov_m: float = FOV_M) -> tuple[np.ndarray, np.ndarray]:
    """
    Return ``(x, y)`` pixel-centre coordinates in metres, matching sigpy's grid convention.

    Written out rather than borrowed: the pixel-centre convention is one of the things under test,
    so taking it from the library being tested would make the comparison vacuous.  An ``N``-point
    axis spans the field of view with samples at ``(i - N/2) * FOV / N``.
    """
    axis = (np.arange(matrix) - matrix // 2) * (fov_m / matrix)
    return np.meshgrid(axis, axis, indexing='ij')


def dense_dft(k_per_m: np.ndarray, image: np.ndarray, fov_m: float = FOV_M) -> np.ndarray:
    """
    Return ``(shots, samples)``: the signal an image produces, by explicit summation.

    ``sum_r x(r) exp(-2 pi i k . r)``, with no FFT, no gridding and no sigpy.  The reference.
    """
    x, y = image_grid(image.shape[0], fov_m)
    flat = image.reshape(-1)
    rx, ry = x.reshape(-1), y.reshape(-1)
    shots, _, samples = k_per_m.shape
    out = np.empty((shots, samples), dtype=np.complex128)
    for shot in range(shots):
        phase = np.exp(-2j * np.pi * (np.outer(k_per_m[shot, 0], rx)
                                      + np.outer(k_per_m[shot, 1], ry)))
        out[shot] = phase @ flat
    return out


# --------------------------------------------------------------------- the two consumers
def _from_blocks(blocks: list, opts: pp.Opts, *, samples: int,
                 echo_times_s: tuple[float, ...]) -> NonCartesianAcquisition:
    """
    Build an acquisition by **measuring** each shot's compiled trajectory.

    The same path a notebook takes: compose the block, run `sc.kspace` over it, keep `k_adc`.  It
    is also the honest one -- what the reconstruction is handed must be what the scanner plays,
    not what a design pass intended.
    """
    k = np.stack([sc.kspace(block, opts)['k_adc'][:2, :samples] for block in blocks])
    times = sc.kspace(blocks[0], opts)['t_adc'][:samples]
    return NonCartesianAcquisition(k_per_m=k, t_adc_s=times - times[0] + (times[1] - times[0]) / 2,
                                   fov_m=FOV_M, matrix=MATRIX, echo_times_s=echo_times_s)


def radial_acquisition(opts: pp.Opts, spokes: int = 12) -> NonCartesianAcquisition:
    """A real `RadialReadout`: spokes stacked on the shot axis."""
    spoke = sc.modules.RadialReadout(opts=opts, fov_mm=FOV_M * 1e3, matrix=MATRIX, dwell_s=20e-6)
    blocks = [spoke(angle_rad=np.pi * i / spokes) for i in range(spokes)]
    dwell = spoke.dwell_s
    return _from_blocks(blocks, opts, samples=spoke.num_samples,
                        echo_times_s=((spoke.center_sample + 0.5) * dwell,))


def spiral_acquisition(opts: pp.Opts) -> NonCartesianAcquisition:
    """A real `SpiralReadout`: interleaves stacked on the shot axis."""
    spiral = sc.modules.SpiralReadout(opts=opts, fov_mm=FOV_M * 1e3, matrix=MATRIX, shots=4,
                                      dwell_s=1e-5, variant='out')
    blocks = [spiral(angle_rad=2.0 * np.pi * i / 4) for i in range(4)]
    offset = spiral.time_to_echo(0) - spiral.sample_times_s()[0] + spiral.dwell_s / 2
    return _from_blocks(blocks, opts, samples=spiral.num_samples, echo_times_s=(offset,))


CONSUMERS = ('radial', 'spiral')


def acquisition(kind: str, opts: pp.Opts) -> NonCartesianAcquisition:
    """Either real consumer, by name."""
    return radial_acquisition(opts) if kind == 'radial' else spiral_acquisition(opts)


# ------------------------------------------------------------------ coordinate convention
@pytest.mark.parametrize('kind', CONSUMERS)
def test_coord_is_cycles_per_fov(kind: str, opts: pp.Opts) -> None:
    """
    Scale: a sample at ``k_max`` lands at the edge of the grid, ``matrix / 2``.

    This is the whole adapter, and getting it wrong scales the image without otherwise breaking it.
    """
    acq = acquisition(kind, opts)
    coord = sigpy_coord(acq)
    assert coord.shape == (acq.shots, acq.samples, 2)
    reach = float(np.abs(coord).max())
    assert reach == pytest.approx(MATRIX / 2, rel=0.05), f'{reach:.3f} against {MATRIX / 2}'


@pytest.mark.parametrize('kind', CONSUMERS)
def test_coord_preserves_axis_order_and_sign(kind: str, opts: pp.Opts) -> None:
    """A transpose or a sign flip here is a mirrored image, which looks entirely plausible."""
    acq = acquisition(kind, opts)
    coord = sigpy_coord(acq)
    assert np.allclose(coord[..., 0], acq.k_per_m[:, 0] * acq.fov_m)
    assert np.allclose(coord[..., 1], acq.k_per_m[:, 1] * acq.fov_m)


# ------------------------------------------------------- the forward model, against physics
@pytest.mark.parametrize('kind', CONSUMERS)
def test_sigpy_forward_matches_a_dense_dft(kind: str, opts: pp.Opts) -> None:
    """
    The central check: the sigpy operator through our adapter reproduces explicit summation.

    If the adapter is wrong in scale, sign, axis order or origin, this fails -- and it is the only
    test here whose reference is computed without sigpy.
    """
    acq = acquisition(kind, opts)
    rng = np.random.default_rng(0)
    image = (rng.standard_normal((MATRIX, MATRIX))
             + 1j * rng.standard_normal((MATRIX, MATRIX))).astype(np.complex128)

    reference = dense_dft(acq.k_per_m, image)
    through_sigpy = np.asarray(encoding_operator(acq)(image))

    error = np.abs(through_sigpy - reference).max() / np.abs(reference).max()
    # 1 % is sigpy's gridding-interpolation accuracy at its default kernel width, not slack for a
    # convention error: the next test shows a wider kernel shrinks this, which a wrong scale,
    # sign or axis order would not.
    assert error < 1e-2, f'relative error {error:.3g}'


@pytest.mark.parametrize('kind', CONSUMERS)
def test_the_residual_is_interpolation_and_not_convention(kind: str, opts: pp.Opts) -> None:
    """
    A wider gridding kernel shrinks the disagreement; a wrong convention would not move.

    This is what separates "the NUFFT is approximate" from "the adapter is wrong", and without it
    the tolerance above would just be a number chosen to pass.
    """
    acq = acquisition(kind, opts)
    rng = np.random.default_rng(4)
    image = (rng.standard_normal((MATRIX, MATRIX))
             + 1j * rng.standard_normal((MATRIX, MATRIX))).astype(np.complex128)
    reference = dense_dft(acq.k_per_m, image)

    def error_at(width: int) -> float:
        operator = float(MATRIX) * sigpy.linop.NUFFT(
            (MATRIX, MATRIX), sigpy_coord(acq), oversamp=2.0, width=width,
        )
        got = np.asarray(operator(image))
        return float(np.abs(got - reference).max() / np.abs(reference).max())

    assert error_at(8) < 0.2 * error_at(3)


@pytest.mark.parametrize('kind', CONSUMERS)
def test_a_point_carries_the_analytic_phase(kind: str, opts: pp.Opts) -> None:
    """
    One off-centre point, and the **complex phase** of every sample compared with ``-2 pi k . r``.

    A magnitude image hides a sign error and an axis swap; this does not.  The point is placed off
    both axes and asymmetrically, so swapping x and y changes the answer.
    """
    acq = acquisition(kind, opts)
    x, y = image_grid()
    row, col = 11, 6                       # off-centre, and not on the diagonal
    image = np.zeros((MATRIX, MATRIX), dtype=np.complex128)
    image[row, col] = 1.0
    position = np.array([x[row, col], y[row, col]])

    analytic = np.exp(-2j * np.pi * (acq.k_per_m[:, 0] * position[0]
                                     + acq.k_per_m[:, 1] * position[1]))
    measured = np.asarray(encoding_operator(acq)(image))

    phase_error = np.abs(np.angle(measured * np.conj(analytic)))
    assert phase_error.max() < 0.02, f'worst phase error {phase_error.max():.4g} rad'


# ------------------------------------------------------------------------- adjointness
@pytest.mark.parametrize('kind', CONSUMERS)
def test_the_operator_is_adjoint(kind: str, opts: pp.Opts) -> None:
    """
    ``<Ax, y> == <x, A^H y>``, numerically.

    Composed from sigpy primitives so this holds by construction -- which is the argument for
    composing rather than hand-writing, and this test is what makes the argument checkable.  A CG
    on a nearly-adjoint pair converges to the wrong image without raising.
    """
    acq = acquisition(kind, opts)
    operator = encoding_operator(acq)
    rng = np.random.default_rng(1)
    x = (rng.standard_normal(operator.ishape) + 1j * rng.standard_normal(operator.ishape))
    y = (rng.standard_normal(operator.oshape) + 1j * rng.standard_normal(operator.oshape))

    left = np.vdot(y, operator(x))
    right = np.vdot(operator.H(y), x)
    assert abs(left - right) / max(abs(left), 1e-30) < 1e-5


def test_adjointness_survives_off_resonance_segmentation(opts: pp.Opts) -> None:
    """The segmented operator is a sum of linops, so it is adjoint too -- asserted, not assumed."""
    acq = spiral_acquisition(opts)
    field = np.full((MATRIX, MATRIX), 40.0)
    operator = encoding_operator(acq, b0_hz=field, segments=6)
    rng = np.random.default_rng(2)
    x = (rng.standard_normal(operator.ishape) + 1j * rng.standard_normal(operator.ishape))
    y = (rng.standard_normal(operator.oshape) + 1j * rng.standard_normal(operator.oshape))
    assert abs(np.vdot(y, operator(x)) - np.vdot(operator.H(y), x)) \
        / max(abs(np.vdot(y, operator(x))), 1e-30) < 1e-5


# ------------------------------------------------------------- ordering and timing semantics
@pytest.mark.parametrize('kind', CONSUMERS)
def test_samples_are_consumed_in_acquisition_order(kind: str, opts: pp.Opts) -> None:
    """
    Permuting the samples permutes the data the same way.

    Not a tautology: it fails if anything inside sorts, reshapes across the shot axis, or pairs
    sample *n* of the trajectory with a different sample of the data.  A reordering produces a
    scrambled image rather than an error.
    """
    acq = acquisition(kind, opts)
    rng = np.random.default_rng(3)
    image = (rng.standard_normal((MATRIX, MATRIX))
             + 1j * rng.standard_normal((MATRIX, MATRIX))).astype(np.complex128)
    straight = np.asarray(encoding_operator(acq)(image))

    order = rng.permutation(acq.samples)
    shuffled = NonCartesianAcquisition(
        k_per_m=acq.k_per_m[:, :, order], t_adc_s=acq.t_adc_s[order], fov_m=acq.fov_m,
        matrix=acq.matrix, echo_times_s=acq.echo_times_s,
    )
    assert np.allclose(np.asarray(encoding_operator(shuffled)(image)), straight[:, order],
                       atol=1e-8)


def test_echo_time_is_explicit_not_inferred(opts: pp.Opts) -> None:
    """
    The container takes the echo instant; it never derives one.

    ``min(t_adc_s)`` is right for a spiral-out for exactly one reason -- its first sample *is* its
    echo -- and wrong for everything else.  ``'in-out'`` makes that concrete: its echo is the seam
    in the middle of the acquisition, so the inferred answer and the true one are far apart.

    ``'in'`` would be the sharper example and cannot be used, because ``SpiralReadout`` currently
    reports **no** origin crossing for it -- see ``plans/spiral/findings.md`` section 9.  That is a
    sequence-side defect, isolated rather than worked around here.
    """
    inout = sc.modules.SpiralReadout(opts=opts, fov_mm=FOV_M * 1e3, matrix=MATRIX, shots=4,
                                     dwell_s=1e-5, variant='in-out')
    times = inout.sample_times_s()
    offset = inout.time_to_echo(0) - times[0] + inout.dwell_s / 2
    acq = _from_blocks([inout(angle_rad=0.0)], opts, samples=inout.num_samples,
                       echo_times_s=(offset,))
    inferred = float(acq.t_adc_s.min())
    assert abs(acq.echo_times_s[0] - inferred) > 0.25 * np.ptp(acq.t_adc_s)
    assert acq.echo_relative_times(0)[0] < 0.0 < acq.echo_relative_times(0)[-1]


def test_absolute_looking_times_are_refused(opts: pp.Opts) -> None:
    """An echo outside the readout means somebody passed sequence time, and it is caught here."""
    acq = radial_acquisition(opts)
    with pytest.raises(ValueError, match='outside the readout'):
        NonCartesianAcquisition(k_per_m=acq.k_per_m, t_adc_s=acq.t_adc_s, fov_m=acq.fov_m,
                                matrix=acq.matrix, echo_times_s=(0.25,))


# ------------------------------------------------------------------------ nothing is trajectory-aware
def test_the_core_treats_both_consumers_identically(opts: pp.Opts) -> None:
    """
    One core, two consumers.  The shared path must not branch on trajectory family.

    Checked by source inspection rather than by behaviour, because the failure this guards against
    is a helper that became "generic" by being renamed.
    """
    source = (Path(__file__).resolve().parents[2] / 'examples' / 'noncartesian_recon.py').read_text()
    body = '\n'.join(line for line in source.splitlines()
                     if not line.lstrip().startswith('#'))
    for word in ('spiral', 'Spiral', 'radial_readout', 'RadialReadout'):
        assert f'if {word}' not in body and f'elif {word}' not in body
    # `'radial'` survives only as the name of one DCF option, never as a branch on the input.
    assert body.count("== 'radial'") <= 1


# ------------------------------------------- how the density estimate is USED, measured not assumed
def _known_image() -> np.ndarray:
    """A small structured phantom: a disc, a square and a point, so errors have somewhere to show."""
    x, y = image_grid()
    image = np.zeros((MATRIX, MATRIX), dtype=np.complex128)
    image[np.hypot(x, y) < 0.25 * FOV_M] = 1.0
    image[3:6, 3:6] = 0.5
    image[11, 6] = 2.0
    return image


def _error_against(truth: np.ndarray, estimate: np.ndarray) -> float:
    """Relative error after removing the global scale a solve is free to choose."""
    scale = np.vdot(estimate, truth) / max(np.vdot(estimate, estimate).real, 1e-30)
    return float(np.linalg.norm(scale * estimate - truth) / np.linalg.norm(truth))


@pytest.mark.parametrize('kind', CONSUMERS)
def test_pipe_menon_weights_are_sane(kind: str, opts: pp.Opts) -> None:
    """Positive, finite, mean one, and heavier away from a centre every trajectory oversamples."""
    acq = acquisition(kind, opts)
    weights = dcf(acq, method='pipe')
    assert weights.shape == (acq.shots, acq.samples)
    assert np.all(np.isfinite(weights)) and np.all(weights > 0)
    assert float(weights.mean()) == pytest.approx(1.0)
    radius = np.hypot(acq.k_per_m[:, 0], acq.k_per_m[:, 1])
    inner, outer = radius < 0.25 * radius.max(), radius > 0.75 * radius.max()
    assert weights[outer].mean() > weights[inner].mean()


@pytest.mark.parametrize('kind', CONSUMERS)
def test_gridding_needs_the_weights_and_a_solve_does_not(kind: str, opts: pp.Opts) -> None:
    """
    The distinction the API exists to make: an estimate, and three ways of using it.

    Gridding -- one adjoint -- is the use that *depends* on density compensation, because nothing
    else corrects the sampling density.  A conjugate gradient does not: it reaches the
    least-squares solution either way, and the weights only change how fast.
    """
    acq = acquisition(kind, opts)
    truth = _known_image()
    data = dense_dft(acq.k_per_m, truth)          # the independent oracle, again
    weights = dcf(acq, method='pipe')

    plain_grid = _error_against(truth, reconstruct(data, acq, weights=None, use='none',
                                                   iterations=0))
    weighted_grid = _error_against(truth, reconstruct(data, acq, weights=weights, use='gridding'))
    assert weighted_grid < plain_grid, f'{weighted_grid:.3f} against {plain_grid:.3f}'

    solved = _error_against(truth, reconstruct(data, acq, weights=None, use='none',
                                               iterations=40))
    assert solved < weighted_grid


@pytest.mark.parametrize('kind', CONSUMERS)
def test_the_two_weighted_spellings_are_one_estimator(kind: str, opts: pp.Opts) -> None:
    """
    **This is the measurement that decides the policy**, rather than inheriting one.

    ``'preconditioner'`` and ``'data'`` are the **same** estimator, and the names say otherwise.
    Both assemble ``A^H W A x = A^H W y``, because ``(W^(1/2)A)^H (W^(1/2)A) = A^H W A`` and
    ``(W^(1/2)A)^H (W^(1/2)y) = A^H W y``.  Checked here on the operator itself, not on a solve,
    so no iteration count or tolerance is involved.

    **This test replaces one that asserted the opposite.**  The first version claimed the two
    converged to different images and compared the two distances with ``>``; they differ by
    about 4e-4 of each other, so it passed locally and failed in CI on floating-point ordering.
    The assertion was not flaky -- it was false, and the flake is what exposed it.  A true
    preconditioner would leave the fixed point at the unweighted solution; this moves it, which
    makes it a weighted least squares whatever it is called.
    """
    acq = acquisition(kind, opts)
    operator = encoding_operator(acq)
    weights = dcf(acq, method='pipe')
    image = _known_image().astype(complex)
    data = dense_dft(acq.k_per_m, _known_image())
    root = np.sqrt(weights)

    # use='preconditioner' assembles these two ...
    normal_w = operator.H(weights * operator(image))
    rhs_w = operator.H(weights * data)
    # ... and use='data' assembles these, by applying sqrt(w) to both sides first.
    normal_root = operator.H(root * (root * operator(image)))
    rhs_root = operator.H(root * (root * data))

    assert np.abs(normal_w - normal_root).max() < 1e-9 * np.abs(normal_w).max()
    assert np.abs(rhs_w - rhs_root).max() < 1e-9 * np.abs(rhs_w).max()


@pytest.mark.parametrize('kind', CONSUMERS)
def test_weighting_matters_only_when_the_data_are_inconsistent(kind: str,
                                                               opts: pp.Opts) -> None:
    """
    Why the DCF policy question has the answer it has.

    On **noiseless, consistent** data the weighted and unweighted solves reach the same image:
    both can drive the residual to zero, and a weight cannot change where zero is.  So density
    compensation is not a better estimator there -- it is a faster route, and section 3 of
    ``examples/gre_radial_2d/02`` measures the route.

    Add noise and the system is no longer consistent, the weights decide *which* residual is
    minimised, and the two answers separate.  That is the regime in which calling this a choice
    of estimator means something.
    """
    acq = acquisition(kind, opts)
    clean = dense_dft(acq.k_per_m, _known_image())
    weights = dcf(acq, method='pipe')
    rng = np.random.default_rng(7)
    noise = rng.standard_normal(clean.shape) + 1j * rng.standard_normal(clean.shape)
    noisy = clean + 0.05 * float(np.abs(clean).mean()) * noise

    def difference(a: np.ndarray, b: np.ndarray) -> float:
        scale = np.vdot(b, a) / max(np.vdot(b, b).real, 1e-30)
        return float(np.linalg.norm(scale * b - a) / np.linalg.norm(a))

    def solve(data, **kwargs):
        return reconstruct(data, acq, iterations=120, **kwargs)

    consistent = difference(solve(clean, weights=None, use='none'),
                            solve(clean, weights=weights, use='preconditioner'))
    inconsistent = difference(solve(noisy, weights=None, use='none'),
                              solve(noisy, weights=weights, use='preconditioner'))

    assert consistent < 0.05
    assert inconsistent > 4.0 * consistent


@pytest.mark.parametrize('kind', CONSUMERS)
def test_every_documented_use_reconstructs_the_image(kind: str, opts: pp.Opts) -> None:
    """Whatever the policy, none of the four options may be quietly broken."""
    acq = acquisition(kind, opts)
    truth = _known_image()
    data = dense_dft(acq.k_per_m, truth)
    weights = dcf(acq, method='pipe')
    for use in ('preconditioner', 'data', 'gridding', 'none'):
        image = reconstruct(data, acq, weights=weights, use=use, iterations=40)
        assert _error_against(truth, image) < 0.5, use


def test_radial_weighting_is_offered_and_differs_from_pipe(opts: pp.Opts) -> None:
    """``|k|`` is exact for uniform spokes; keeping it makes the comparison available."""
    acq = radial_acquisition(opts)
    assert not np.allclose(dcf(acq, method='radial'), dcf(acq, method='pipe'), rtol=0.2)
    assert np.allclose(dcf(acq, method='none'), 1.0)


def test_an_unknown_dcf_method_is_refused(opts: pp.Opts) -> None:
    """The spiral analytic Jacobian is deliberately not here; asking for it says so."""
    with pytest.raises(ValueError, match="must be 'pipe', 'radial' or 'none'"):
        dcf(radial_acquisition(opts), method='analytic')
