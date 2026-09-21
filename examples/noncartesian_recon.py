r"""
The forward model a non-Cartesian readout is reconstructed with, and the solve that inverts it.

Beside :mod:`phantom` and for the same reason: ``phantom.py`` established that shared example code
sits at the top of ``examples/``.  **Not** in ``src/seqcraft/``.  The division of labour is

.. code-block:: text

    SeqCraft      builds pulse sequences
    MRzeroCore    simulates them
    SigPy         reconstructs the resulting data

and a reconstruction module inside the package would put a reconstruction dependency on the
compile path.  This is example infrastructure, not API.

One equation
------------
.. math::

    y_{j}(t) = \sum_r x(r) \exp\bigl( -2 \pi i (k_j(t) \cdot r + \Delta f(r)\, t) \bigr)

``j`` interleaves, and ``t`` measured **from the echo** -- never absolute, and never inferred.
Both terms carry the same sign, which is the scanner's convention; written with a minus on the
second, an off-resonance correction adds exactly the phase it is meant to remove and the result
reads as *"correction does not help"* rather than as a sign error.  That is how a sign error
survives being tested.

What this owns, and what it borrows
-----------------------------------
SigPy is the **backend**.  It is not evidence that the conversion into its convention is right.
What this file owns is one line --

.. code-block:: text

    physical k [1/m]  x  FOV [m]  ->  cycles per FOV  ->  sigpy `coord`

-- and being one line is why it is dangerous rather than why it is safe.  A wrong scale, a swapped
axis, a flipped sign or an off-by-half-a-pixel origin all produce an image that looks like an
image.  :mod:`tests.examples.test_noncartesian_recon` checks it against an **explicit dense DFT**
written from :math:`\exp(-2\pi i k \cdot r)` that calls no SigPy at all, and against the analytic
phase of a single off-centre point, which is where a sign error shows and a magnitude image does
not.

Two consumers, one core
-----------------------
:class:`~seqcraft.modules.RadialReadout` and :class:`~seqcraft.modules.SpiralReadout` both drive
everything here.  Nothing in this file asks which one it has.  Where the mathematics genuinely is
trajectory-specific -- a variable-density spiral's analytic Jacobian, say -- it belongs in the
notebook that needs it, not behind a general name.

**Needs ``seqcraft[recon]``** -- ``sigpy`` -- and nothing else.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import sigpy as sp
import sigpy.mri as smri

__all__ = [
    'NonCartesianAcquisition',
    'dcf',
    'encoding_operator',
    'reconstruct',
    'sigpy_coord',
]


@dataclass(frozen=True)
class NonCartesianAcquisition:
    """
    Where the samples are, when they were taken, and which instant is the echo.

    Trajectory-neutral on purpose.  A radial acquisition is spokes stacked on the shot axis and a
    spiral is interleaves; nothing below distinguishes them.

    Parameters
    ----------
    k_per_m
        ``(shots, 2, samples)`` in 1/m, **in acquisition order**.  The order matters: the
        reconstruction pairs sample *n* of this array with sample *n* of the data, and a
        reordering here is a scrambled image rather than an error.
    t_adc_s
        ``(samples,)`` from the **readout block's** start -- the same origin `echo_times_s` uses.
    fov_m, matrix
        The image this is reconstructed onto.
    echo_times_s
        When ``k = 0`` is crossed, one time per echo, **explicit**.  A tuple because a
        spiral ``'out-in'`` crosses twice, and because inferring it from ``min(t_adc_s)`` is
        right for a spiral-out for exactly one reason -- its first sample *is* its echo -- and
        wrong for ``'in'``, ``'in-out'``, every multi-echo train and every EPI.
    """

    k_per_m: np.ndarray
    t_adc_s: np.ndarray
    fov_m: float
    matrix: int
    echo_times_s: tuple[float, ...] = (0.0,)

    def __post_init__(self) -> None:
        """Refuse a shape that will not broadcast, and times that look absolute."""
        k = np.atleast_3d(np.asarray(self.k_per_m, dtype=np.float64))
        object.__setattr__(self, 'k_per_m', k)
        object.__setattr__(self, 't_adc_s', np.asarray(self.t_adc_s, dtype=np.float64).ravel())
        object.__setattr__(self, 'echo_times_s',
                           tuple(float(t) for t in np.atleast_1d(self.echo_times_s)))
        if k.ndim != 3 or k.shape[1] != 2 or k.shape[2] != self.t_adc_s.size:
            msg = (f'k_per_m must be (shots, 2, samples) matching t_adc_s; got {k.shape} against '
                   f'{self.t_adc_s.size} sample times.')
            raise ValueError(msg)
        span = float(np.ptp(self.t_adc_s)) if self.t_adc_s.size > 1 else 0.0
        for echo, when in enumerate(self.echo_times_s):
            if not self.t_adc_s.min() - span <= when <= self.t_adc_s.max() + span:
                msg = (f'echo_times_s[{echo}] = {when * 1e3:.3f} ms is outside the readout, whose '
                       f'samples run {self.t_adc_s.min() * 1e3:.3f} to '
                       f'{self.t_adc_s.max() * 1e3:.3f} ms.  Both are measured from the readout '
                       f'block, not from the excitation -- absolute times put an extra '
                       f'2 pi df TE on every sample, which varies with position and so reads as '
                       f'a structured artefact rather than as a global phase.')
                raise ValueError(msg)

    @property
    def shots(self) -> int:
        """Interleaves or spokes, stacked on the sample axis by the operator."""
        return int(self.k_per_m.shape[0])

    @property
    def samples(self) -> int:
        """Samples per shot."""
        return int(self.k_per_m.shape[2])

    def echo_relative_times(self, echo: int = 0) -> np.ndarray:
        """Sample times measured **from** `echo` -- what the off-resonance term multiplies."""
        return self.t_adc_s - self.echo_times_s[echo]


def sigpy_coord(acquisition: NonCartesianAcquisition) -> np.ndarray:
    """
    Return the trajectory in sigpy's convention: ``(shots, samples, 2)`` in cycles per FOV.

    **The whole adapter.**  SigPy's ``coord`` counts cycles across the reconstructed field of
    view, so a sample at :math:`k_{\\max} = N / (2\\,\\mathrm{FOV})` lands at :math:`N/2`, the edge
    of an ``N``-point grid.  Multiplying physical k by the field of view is exactly that change of
    units and nothing else -- no sign, no transpose, no offset.

    Everything that could be wrong here is invisible downstream, so it is checked against a dense
    DFT that shares no code with sigpy.
    """
    return np.moveaxis(acquisition.k_per_m, 1, -1) * float(acquisition.fov_m)


def dcf(acquisition: NonCartesianAcquisition, *, method: str = 'pipe') -> np.ndarray:
    """
    Density compensation weights, ``(shots, samples)``, normalised to a mean of one.

    **This estimates the density.  It does not decide how the estimate is used** -- see
    :func:`reconstruct`, where the same weights can act as a gridding weight, as a least-squares
    data weighting or as a preconditioner, and those are three different estimators rather than
    three spellings of one.

    Parameters
    ----------
    method
        ``'pipe'`` -- :func:`sigpy.mri.pipe_menon_dcf`, the default.  It needs no analytic density
        and is what a real pipeline uses, which is also why it is the one a teaching example
        should show.
        ``'radial'`` -- ``|k|``, exact for uniform radial spokes and a useful thing to compare
        against because its failure on a variable-density trajectory is instructive.
        ``'none'`` -- ones.

    Notes
    -----
    A variable-density spiral's analytic Jacobian is ``|dk/dt| / FOV(|k|)``, and it is **not**
    here: it needs the density the path was designed with, which is a property of one trajectory
    family rather than of non-Cartesian sampling.  A notebook that wants it can compute it from
    its own design parameters in four lines.
    """
    k = acquisition.k_per_m
    if method == 'none':
        return np.ones((acquisition.shots, acquisition.samples), dtype=np.float64)
    if method == 'radial':
        return _normalise(np.hypot(k[:, 0], k[:, 1]))
    if method != 'pipe':
        msg = f"method must be 'pipe', 'radial' or 'none', got {method!r}."
        raise ValueError(msg)
    coord = sigpy_coord(acquisition).reshape(-1, 2)
    weights = np.asarray(smri.pipe_menon_dcf(coord, img_shape=(acquisition.matrix,) * 2,
                                             show_pbar=False))
    return _normalise(weights.reshape(acquisition.shots, acquisition.samples))


def _normalise(weights: np.ndarray) -> np.ndarray:
    """Weights with a mean of one, so a regularisation strength means the same thing throughout."""
    weights = np.abs(np.asarray(weights, dtype=np.float64))
    return weights / max(float(weights.mean()), np.finfo(float).tiny)


def encoding_operator(acquisition: NonCartesianAcquisition, *,
                      b0_hz: np.ndarray | float = 0.0, echo: int = 0,
                      segments: int = 1) -> sp.linop.Linop:
    """
    Return ``A``: image ``(matrix, matrix)`` -> data ``(shots, samples)``.

    **Composed from sigpy primitives rather than hand-written**, so sigpy derives ``A.H`` itself
    and the pair is adjoint *by construction*.  Writing both directions by hand is where a forward
    and an adjoint come to disagree -- and a conjugate gradient on a pair that is only *nearly*
    adjoint does not raise.  It converges to the wrong image and looks under-regularised.

    Parameters
    ----------
    b0_hz
        Off-resonance in hertz, either a constant or a ``(matrix, matrix)`` map.
    echo
        Which origin crossing the sample times are measured from.
    segments
        Time segments for the off-resonance term.  ``1`` is no correction at all and is the right
        default: a reader should see the uncorrected image before the correction.

    Notes
    -----
    Interleaves are **one** inverse problem, not several averaged: solving them separately throws
    away exactly the joint conditioning that makes multi-shot work.  That is why the operator maps
    one image to all shots at once.
    """
    coord = sigpy_coord(acquisition)
    shape = (acquisition.matrix, acquisition.matrix)
    # sigpy's NUFFT is built on its *centred, unitary* FFT, so it divides by sqrt(N_pixels).
    # The physical forward model has no such factor -- it is a plain sum over the object -- so
    # the adapter restores it.  Left out, every image comes back scaled by the matrix size, which
    # is invisible in a windowed magnitude image and wrong the moment two reconstructions are
    # compared or a regularisation strength means anything.  Found by the dense-DFT check, which
    # disagreed by exactly 15/16.
    nufft = float(np.sqrt(np.prod(shape))) * sp.linop.NUFFT(shape, coord)

    field = np.asarray(b0_hz, dtype=np.float64)
    if segments <= 1 or not np.any(field):
        return nufft

    times = acquisition.echo_relative_times(echo)
    # Time segmentation: A = sum_l Diag(w_l) . NUFFT . Diag(exp(-2i pi b0 t_l)).  Every factor is
    # an sp.linop, so the sum is one too, and its adjoint comes for free.
    edges = np.linspace(times.min(), times.max(), int(segments))
    width = max(float(edges[1] - edges[0]), np.finfo(float).tiny) if len(edges) > 1 else 1.0
    operator = None
    for centre in edges:
        weight = np.clip(1.0 - np.abs(times - centre) / width, 0.0, 1.0)
        if not np.any(weight):
            continue
        phase = np.exp(-2j * np.pi * field * centre)
        term = sp.linop.Multiply(nufft.oshape, np.broadcast_to(weight, nufft.oshape)) * nufft \
            * sp.linop.Multiply(shape, np.broadcast_to(phase, shape).astype(np.complex128))
        operator = term if operator is None else operator + term
    return operator if operator is not None else nufft


def reconstruct(data: np.ndarray, acquisition: NonCartesianAcquisition, *,
                weights: np.ndarray | None = None, use: str = 'preconditioner',
                b0_hz: np.ndarray | float = 0.0, echo: int = 0, segments: int = 1,
                iterations: int = 30, lamda: float = 0.0) -> np.ndarray:
    """
    Solve ``A x = y`` for the image, by conjugate gradient on the normal equations.

    Parameters
    ----------
    data
        ``(shots, samples)``, in acquisition order.
    weights
        Density compensation from :func:`dcf`, or ``None`` for unweighted.
    use
        **How the weights are used, which is a separate question from what they are.**

        ``'preconditioner'`` -- ``w`` inside the conjugate gradient, changing the path and not the
        fixed point.  The solution is the unweighted least-squares one, reached sooner.
        ``'data'`` -- ``sqrt(w)`` applied to both the operator and the data, which is a *weighted*
        least squares: a different estimator, not a faster route to the same one.
        ``'gridding'`` -- no solve at all, just ``A^H (w y)``.  One adjoint, and what "density
        compensation" classically means.
        ``'none'`` -- ignore the weights.

        Which of these a notebook should use is an evidence question, not a default to inherit;
        ``tests/examples/test_noncartesian_recon.py`` measures them against a known image.
    lamda
        Tikhonov strength, in units where the weights' mean is one.
    """
    operator = encoding_operator(acquisition, b0_hz=b0_hz, echo=echo, segments=segments)
    measured = np.asarray(data, dtype=np.complex128)
    if weights is None or use == 'none':
        weighting = None
    else:
        weighting = np.asarray(weights, dtype=np.float64)

    if weighting is not None and use == 'gridding':
        return np.asarray(operator.H(measured * weighting))

    if weighting is not None and use == 'data':
        root = sp.linop.Multiply(operator.oshape, np.sqrt(weighting))
        operator, measured = root * operator, measured * np.sqrt(weighting)
        weighting = None

    precondition = None
    if weighting is not None and use == 'preconditioner':
        precondition = sp.linop.Multiply(operator.oshape, weighting)
    elif weighting is not None and use != 'preconditioner':
        msg = f"use must be 'preconditioner', 'data', 'gridding' or 'none', got {use!r}."
        raise ValueError(msg)

    normal = operator.H * precondition * operator if precondition is not None \
        else operator.H * operator
    rhs = operator.H(precondition(measured)) if precondition is not None else operator.H(measured)
    if lamda:
        normal = normal + lamda * sp.linop.Identity(operator.ishape)
    image = np.zeros(operator.ishape, dtype=np.complex128)
    if iterations <= 0:
        return np.asarray(rhs)
    # Driven explicitly: in this sigpy, ConjugateGradient is an `Alg` rather than an `App`, so it
    # has `update()` and `done()` and no `run()`.  Writing the loop also makes the iteration count
    # a visible parameter of the comparison rather than something buried in a solver object.
    solver = sp.app.ConjugateGradient(normal, rhs, image, max_iter=int(iterations))
    while not solver.done():
        solver.update()
    return np.asarray(image)
