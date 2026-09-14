"""
The forward model a non-Cartesian readout is reconstructed with, and the solve that inverts it.

Beside :mod:`phantom` and for the same reason: ``phantom.py`` established that shared example code
sits at the top of ``examples/``.  **Not** in ``src/seqcraft/``.  The package builds sequences, its
only optional imports are matplotlib and sigpy-for-SLR, and a reconstruction module in ``src/``
would put an ESPIRiT dependency on the compile path.

One equation, and correcting off-resonance means putting ``df`` *into* it rather than adding a step
after it::

    y_{c,j,e}(t) = sum_r  s_c(r) x_e(r) exp( -2i pi ( k_j(t) . r + df(r) t ) )

``c`` coils, ``j`` interleaves, ``e`` echoes, and ``t`` measured **from that echo**.  Both terms
carry the same sign, which is the scanner's convention; written with a minus on the second, the
correction adds exactly the phase it is meant to remove and the result reads as *"off-resonance
correction does not help"* rather than as a sign error -- which is how a sign error survives being
tested.

Four things this does differently from the usual hand-rolled spiral reconstruction, each of which
was a real failure rather than a preference:

**1. The operator is composed from sigpy primitives, not hand-written.**  Time segmentation is
``A = sum_l Diag(w_l) . NUFFT(coord) . Diag(s_c exp(-2i pi b0 t_l))`` and every factor is an
``sp.linop``, so sigpy derives ``A.H`` itself and it is adjoint **by construction**.  Hand-writing
both ``_apply`` methods is where a forward and an adjoint come to disagree -- and a conjugate
gradient on a pair that is *nearly* adjoint does not raise: it converges to the wrong image and
looks under-regularised.

**2. Interleaves are solved together.**  Four interleaves of one image are one inverse problem, and
solving them separately and averaging throws away exactly the joint conditioning that makes
multi-shot work.  Echoes are the opposite -- separate images, separate solves, one shared coil map.

**3. Density compensation is a preconditioner, not a data weighting.**  ``sqrt(w)`` applied to both
the operator and the data is a *weighted* least squares: a different estimator, not a faster path
to the same one.  On a variable-density spiral the weights vary by two orders of magnitude, so the
difference is not academic.  ``w`` goes inside the CG, where it changes the path and not the fixed
point.  ``weighting='data'`` reproduces the old behaviour and stays reachable, because that is what
makes the comparison a measurement.

**4. The density compensation knows about the density.**  ``'radial'`` weights by ``|k|`` and
``'speed'`` by ``|dk/dt|``, and neither is right on a tapered spiral: the area one sample
represents is the arc spacing times the **radial turn spacing**, and the turn spacing is exactly
what ``density`` changed.  So the analytic Jacobian is ``|dk/dt| / FOV(|k|)``, and the default is
``sigpy.mri.pipe_menon_dcf``, which needs no analytic density and is what a real pipeline uses.

And one thing kept **verbatim in intent** from the code this replaces: :class:`SpiralReadout`
refuses sample times that look absolute.  The off-resonance term is ``exp(-2i pi df(r) t)`` with
``t`` measured *from the echo*; hand it absolute times and every sample carries an extra
``2 pi df TE``, and since ``df`` varies with position that is not a global phase but a structured,
plausible-looking artefact -- 2.2 cycles at TE = 40 ms against a 56 Hz field range.  The old class
inferred the echo from ``min(t_adc_s)``, which is right for a spiral-out for exactly one reason --
its first sample *is* its echo -- and wrong for ``'in'``, ``'in-out'``, every multi-echo train and
every EPI.  Here the echo time is taken explicitly and never inferred.

**Needs ``seqcraft[recon]``** -- ``sigpy`` -- and nothing else.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import sigpy as sp
import sigpy.mri as smri

__all__ = [
    'SpiralReadout', 'dcf', 'encoding_operator', 'field_map', 'reconstruct', 'segments_for',
]


@dataclass(frozen=True)
class SpiralReadout:
    """
    Where the samples are, when they were taken, and which of them is the echo.

    Parameters
    ----------
    k_per_m
        ``(shots, 2, samples)`` in 1/m, in acquisition order.
    t_adc_s
        ``(samples,)`` from the **readout block's** start -- the same origin :attr:`te_s` uses.
    fov_m, matrix
        The image the operator reconstructs onto.
    te_s
        Where ``k = 0`` falls in `t_adc_s`, one time per echo.  **Explicit**, and never inferred
        from ``min(t_adc_s)``: see the module docstring.
    """

    k_per_m: np.ndarray
    t_adc_s: np.ndarray
    fov_m: float
    matrix: int
    te_s: tuple[float, ...] = (0.0,)

    def __post_init__(self) -> None:
        """Refuse a shape that will not broadcast, and times that look absolute."""
        k = np.atleast_3d(np.asarray(self.k_per_m, dtype=np.float64))
        object.__setattr__(self, 'k_per_m', k)
        object.__setattr__(self, 't_adc_s', np.asarray(self.t_adc_s, dtype=np.float64).ravel())
        object.__setattr__(self, 'te_s', tuple(float(t) for t in np.atleast_1d(self.te_s)))
        if k.ndim != 3 or k.shape[1] != 2 or k.shape[2] != self.t_adc_s.size:
            msg = (f'k_per_m must be (shots, 2, samples) matching t_adc_s; got {k.shape} against '
                   f'{self.t_adc_s.size} sample times.')
            raise ValueError(msg)
        span = float(np.ptp(self.t_adc_s)) if self.t_adc_s.size > 1 else 0.0
        for echo, when in enumerate(self.te_s):
            if not self.t_adc_s.min() - span <= when <= self.t_adc_s.max() + span:
                msg = (f'te_s[{echo}] = {when * 1e3:.3f} ms is outside the readout, whose samples '
                       f'run {self.t_adc_s.min() * 1e3:.3f} to {self.t_adc_s.max() * 1e3:.3f} ms. '
                       f'Both are measured from the readout block, not from the excitation.')
                raise ValueError(msg)

    @property
    def shots(self) -> int:
        """Interleaves stacked on the sample axis by the operator."""
        return int(self.k_per_m.shape[0])

    @classmethod
    def from_sidecar(cls, path: str | Path, *, echo: int | None = None) -> SpiralReadout:
        """Load the ``*_nominal.npz`` a ``01_build`` notebook wrote beside its ``.seq`` files."""
        data = np.load(Path(path))
        te = tuple(float(t) for t in np.atleast_1d(data['te_s']))
        return cls(k_per_m=data['k_per_m'], t_adc_s=data['t_adc_s'],
                   fov_m=float(data['fov_mm']) / 1e3, matrix=int(data['matrix']),
                   te_s=te if echo is None else (te[echo],))

    def echo_relative_times(self, echo: int = 0) -> np.ndarray:
        """Sample times measured **from** `echo` -- what the off-resonance term multiplies."""
        return self.t_adc_s - self.te_s[echo]

    def delayed(self, delay_s: float) -> SpiralReadout:
        """
        The same readout with the trajectory shifted `delay_s` later relative to the samples.

        The one calibration neither a sidecar nor a field map supplies: a gradient chain delay
        displaces ``k`` along the arm, which on a spiral is a rotation *and* a radial shift and
        therefore a blur that no phase correction removes.  One parameter, so it is swept rather
        than reasoned about.
        """
        shifted = np.stack([
            np.stack([np.interp(self.t_adc_s + float(delay_s), self.t_adc_s, axis)
                      for axis in shot])
            for shot in self.k_per_m])
        return replace(self, k_per_m=shifted)


#: Segments per cycle of field swing per unit of relative error, **measured** against the exact
#: operator rather than derived: see :func:`segments_for`.
_SEGMENT_CONSTANT = 1.35


def segments_for(b0_hz: np.ndarray, sample_times_s: np.ndarray, *,
                 target_error: float = 0.1) -> int:
    """
    How many time segments the off-resonance operator needs for a given approximation error.

    Derived from an accuracy target rather than a round number of cycles.  Within one segment the
    phase is held at the segment's centre and neighbouring segments are blended linearly, so the
    residual falls with the number of segments and the field swing across the readout,
    ``ptp(df) * T`` in cycles, is the only other thing it depends on::

        relative error  ~  1.35 * ptp(df) * T / segments

    **The constant is measured, not derived.**  Against the exact dense operator on an 18.6 ms
    single-shot arm and a 40 Hz field -- 0.74 cycles of swing -- the relative RMSE runs 1.3e-1 at
    8 segments, 5.9e-2 at 16 and 1.4e-2 at 48, which is that expression to well within the
    accuracy anyone needs from a segment count.  ``examples/gre_spiral_2d/02`` §4 is where it is
    checked, and it is checked because an approximation nobody checks is an assumption.

    `target_error` is the **operator's** relative RMSE, not the image's, and 0.1 is the default
    because that is already far below anything an image shows: the difference between 12 segments
    and 48 on a brain is invisible and costs four times the reconstruction.  Eight is the floor,
    because fewer than that is a piecewise-constant approximation with visible seams whatever the
    field is.
    """
    span_hz = float(np.ptp(np.asarray(b0_hz))) if np.size(b0_hz) else 0.0
    duration = float(np.ptp(np.asarray(sample_times_s)))
    if span_hz <= 0.0 or duration <= 0.0:
        return 1
    return max(8, int(math.ceil(_SEGMENT_CONSTANT * span_hz * duration
                                / float(target_error))))


def dcf(readout: SpiralReadout, *, method: str = 'pipe',
        fov_mm_of_kr: np.ndarray | None = None) -> np.ndarray:
    """
    Density compensation weights, ``(shots, samples)``, normalised to a mean of one.

    Parameters
    ----------
    readout
        The trajectory to weight.
    method
        ``'pipe'``     -- :func:`sigpy.mri.pipe_menon_dcf`, the default and what a real pipeline
        uses; it needs no analytic density at all.
        ``'analytic'`` -- ``|dk/dt| / FOV(|k|)``, the exact Jacobian of a variable-density spiral
        and the thing ``'pipe'`` should be checked against.
        ``'speed'``    -- ``|dk/dt|`` alone, which is ``'analytic'`` with the taper left out.
        ``'radial'``   -- ``|k|``, which is ``'speed'`` with the traversal left out too.
        ``'none'``     -- ones.
    fov_mm_of_kr
        ``(2, n)`` of ``[kr / kmax, FOV in mm]`` -- what ``Spiral2D.fov_mm_of`` produces and what a
        ``01`` notebook writes into its sidecar.  Required by ``'analytic'`` and ignored otherwise.
    """
    k = readout.k_per_m
    if method == 'none':
        return np.ones(k.shape[::2], dtype=np.float64)
    if method == 'pipe':
        coord = _coord(readout).reshape(-1, 2)
        weights = np.asarray(smri.pipe_menon_dcf(coord, img_shape=(readout.matrix,) * 2,
                                                 show_pbar=False))
        return _normalise(weights.reshape(k.shape[0], k.shape[2]))

    radius = np.hypot(k[:, 0], k[:, 1])
    if method == 'radial':
        return _normalise(radius)
    speed = np.abs(np.gradient(k[:, 0] + 1j * k[:, 1], readout.t_adc_s, axis=-1))
    if method == 'speed':
        return _normalise(speed)
    if method != 'analytic':
        msg = f"method must be 'pipe', 'analytic', 'speed', 'radial' or 'none', got {method!r}."
        raise ValueError(msg)
    if fov_mm_of_kr is None:
        msg = "method='analytic' needs fov_mm_of_kr -- the density, sampled; 'pipe' needs nothing."
        raise ValueError(msg)
    kmax = readout.matrix / (2 * readout.fov_m)
    fov = np.interp(radius / kmax, fov_mm_of_kr[0], fov_mm_of_kr[1])
    return _normalise(speed / fov)


def _normalise(weights: np.ndarray) -> np.ndarray:
    """Weights with a mean of one, so `lamda` means the same thing whichever method produced them."""
    weights = np.abs(np.asarray(weights, dtype=np.float64))
    return weights / max(float(weights.mean()), np.finfo(float).tiny)


def _coord(readout: SpiralReadout) -> np.ndarray:
    """The trajectory in sigpy's convention: ``(shots, samples, 2)`` in cycles per FOV."""
    return np.moveaxis(readout.k_per_m, 1, -1) * readout.fov_m


def encoding_operator(readout: SpiralReadout, *, sens: np.ndarray | None = None,
                      b0_hz: np.ndarray | float = 0.0, echo: int = 0,
                      segments: int | None = 8, device: int = -1) -> sp.linop.Linop:
    """
    Return ``A``: image -> ``(coils, shots, samples)``, with off-resonance inside it.

    Parameters
    ----------
    readout
        Where and when the samples are.
    sens
        ``(coils, matrix, matrix)`` receive sensitivities, or ``None`` for a single uniform coil.
    b0_hz
        ``(matrix, matrix)`` field map, or a scalar.  Zero skips the segmentation entirely.
    echo
        Which echo's sample times to measure from.  This is the argument whose absence is a
        structured artefact rather than a global phase.
    segments
        Time segments.  ``None`` is the **exact** dense operator -- ``samples x pixels`` complex
        exponentials per application, the right thing to *check* a reconstruction against and
        useless to perform one with.  It stays because an approximation nobody checks is an
        assumption.
    device
        ``-1`` for CPU, or a sigpy device index.

    Notes
    -----
    Every factor is an ``sp.linop``, so ``A.H`` is derived rather than written and the pair is
    adjoint by construction.  Composing also removes the probe normalisation a hand-written
    operator needs, because the same constant then appears on both sides.
    """
    shape = (readout.matrix, readout.matrix)
    coord = _coord(readout)
    times = readout.echo_relative_times(echo)
    field = np.zeros(shape) + np.asarray(b0_hz, dtype=np.float64)

    coils = sp.linop.Multiply(shape, np.asarray(sens, dtype=np.complex64)) if sens is not None \
        else sp.linop.Reshape((1, *shape), shape)
    nufft = sp.linop.NUFFT(coils.oshape, coord)

    if not np.any(field):
        return nufft * coils
    if segments is None:
        return _ExactOffResonance(coils.oshape, coord, field, times, readout.fov_m) * coils

    # The phase is held at each segment's centre and the segments are blended with the linear
    # interpolation between them.  `position` is **clipped** into the centres' own range first,
    # so the weights sum to one at the two ends of the readout as well as in the middle -- an
    # unclipped triangular kernel decays to a half there, and the resulting amplitude notch at
    # both ends of every arm is larger than the off-resonance error it was meant to remove.
    count = max(int(segments), 1)
    edges = np.linspace(times.min(), times.max(), count + 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    step = centres[1] - centres[0] if count > 1 else 1.0
    position = np.clip((times - centres[0]) / step, 0.0, count - 1)
    operator = None
    for index, centre in enumerate(centres):
        weight = np.clip(1.0 - np.abs(position - index), 0.0, 1.0)
        term = (sp.linop.Multiply(nufft.oshape, weight.astype(np.complex64)) * nufft
                * sp.linop.Multiply(coils.oshape,
                                    np.exp(-2j * np.pi * field * centre).astype(np.complex64)))
        operator = term if operator is None else operator + term
    return operator * coils


class _ExactOffResonance(sp.linop.Linop):
    """The dense ``exp(-2i pi (k.r + df t))`` operator, for checking the segmented one against."""

    def __init__(self, ishape, coord, field_hz, times_s, fov_m):
        matrix = ishape[-1]
        axis = (np.arange(matrix) - matrix // 2) * fov_m / matrix
        grid = np.stack(np.meshgrid(axis, axis, indexing='ij'), axis=-1).reshape(-1, 2)
        # sigpy's NUFFT carries a 1/sqrt(pixels) normalisation, so this one does too -- otherwise
        # the exact operator and the segmented one differ by a constant and the comparison between
        # them measures the constant.
        self._phase = (np.exp(-2j * np.pi * (
            coord.reshape(-1, 2) / fov_m @ grid.T
            + times_s[None, :].repeat(coord.shape[0], 0).reshape(-1, 1)
            * field_hz.reshape(1, -1))) / matrix).astype(np.complex64)
        self._grid_shape = ishape
        self._samples = coord.shape[:2]
        super().__init__((ishape[0], *self._samples), ishape)

    def _apply(self, x):
        flat = np.asarray(x).reshape(self._grid_shape[0], -1)
        return (flat @ self._phase.T).reshape(self.oshape)

    def _adjoint_linop(self):
        return _ExactAdjoint(self)


class _ExactAdjoint(sp.linop.Linop):
    """``_ExactOffResonance`` transposed-conjugate, so the pair is exact rather than nearly so."""

    def __init__(self, forward: _ExactOffResonance):
        self._forward = forward
        super().__init__(forward.ishape, forward.oshape)

    def _apply(self, y):
        flat = np.asarray(y).reshape(self._forward.oshape[0], -1)
        return (flat @ self._forward._phase.conj()).reshape(self.oshape)

    def _adjoint_linop(self):
        return self._forward


def reconstruct(data: np.ndarray, readout: SpiralReadout, *, sens: np.ndarray | None = None,
                b0_hz: np.ndarray | float = 0.0, echo: int = 0, segments: int | None = 0,
                iterations: int = 60, lamda: float = 1e-3, weighting: str = 'precondition',
                dcf_method: str = 'pipe', fov_mm_of_kr: np.ndarray | None = None,
                device: int = -1) -> np.ndarray:
    """
    Solve ``min_x ||A x - y||^2 + lamda ||x||^2`` for one echo, over all interleaves at once.

    Parameters
    ----------
    data
        ``(coils, shots, samples)``, or ``(shots, samples)`` for a single coil.
    readout, sens, b0_hz, echo, device
        As :func:`encoding_operator`.
    segments
        ``0`` sizes it with :func:`segments_for`; ``None`` is the exact operator; an int is that
        many.
    iterations
        60 rather than the 12 a Cartesian-flavoured default would suggest: correcting off-resonance
        *worsens* the operator's conditioning, and at 20 the corrected image carries more
        high-frequency texture than the uncorrected one -- which reads as "correction makes it
        worse" and is non-convergence.
    lamda
        A **fraction of the largest eigenvalue** of ``A^H A``, found by a seeded power iteration.
        Absolute Tikhonov weights are meaningless here: that spectrum scales with the sample count,
        so ``1e-3`` absolute is four orders of magnitude below anything it could affect.
    weighting
        ``'precondition'`` puts the density compensation **inside** the conjugate gradient, as
        the image-domain preconditioner ``1 / (A^H W A 1)`` -- the operator's own row sums
        under the weights, one extra application to compute.  A preconditioner changes the
        path and cannot change the fixed point, so this stays the *unweighted* least-squares
        estimator.  ``'data'`` applies ``sqrt(w)`` to both the operator and the data, which is
        a **weighted** least squares -- a different estimator, and on a variable-density
        spiral the weights vary by two orders of magnitude, so it is a different picture.  It
        stays reachable so the difference is measured rather than asserted.
    dcf_method, fov_mm_of_kr
        Passed to :func:`dcf`.
    """
    y = np.asarray(data, dtype=np.complex64)
    if y.ndim == 2:
        y = y[None]
    times = readout.echo_relative_times(echo)
    count = (segments_for(b0_hz, times) if segments == 0 else segments)
    operator = encoding_operator(readout, sens=sens, b0_hz=b0_hz, echo=echo, segments=count,
                                 device=device)
    weights = dcf(readout, method=dcf_method, fov_mm_of_kr=fov_mm_of_kr).astype(np.float32)

    if weighting == 'data':
        root = sp.linop.Multiply(operator.oshape, np.sqrt(weights).astype(np.complex64))
        operator, y, preconditioner = root * operator, np.sqrt(weights) * y, None
    elif weighting == 'precondition':
        preconditioner = _row_sum_preconditioner(operator, weights)
    else:
        msg = f"weighting must be 'precondition' or 'data', got {weighting!r}."
        raise ValueError(msg)

    app = sp.app.LinearLeastSquares(
        operator, y, lamda=float(lamda) * _spectral_norm(operator),
        solver='ConjugateGradient', P=preconditioner, max_iter=int(iterations),
        show_pbar=False)
    return np.asarray(app.run())


def _row_sum_preconditioner(operator: sp.linop.Linop, weights: np.ndarray) -> sp.linop.Linop:
    """``1 / (A^H W A 1)`` in the image domain: the density compensation, inside the CG.

    One application of the operator and its adjoint, and it is the classic density/intensity
    preconditioner rather than anything invented here.  Floored at a millionth of its own peak
    so a pixel the trajectory never reaches cannot divide by nothing.
    """
    W = sp.linop.Multiply(operator.oshape, weights.astype(np.complex64))
    probe = np.ones(operator.ishape, dtype=np.complex64)
    rows = np.abs(np.asarray(operator.H * (W * (operator * probe))))
    floor = max(float(rows.max()) * 1e-6, np.finfo(np.float32).tiny)
    return sp.linop.Multiply(operator.ishape,
                             (1.0 / np.maximum(rows, floor)).astype(np.complex64))


def _spectral_norm(operator: sp.linop.Linop, *, iterations: int = 15, seed: int = 0) -> float:
    """The largest eigenvalue of ``A^H A``, by a **seeded** power iteration so `lamda` is stable."""
    rng = np.random.default_rng(seed)
    x = (rng.standard_normal(operator.ishape) + 1j * rng.standard_normal(operator.ishape))
    x = x.astype(np.complex64)
    value = 1.0
    for _ in range(iterations):
        x = operator.H * (operator * x)
        value = float(np.linalg.norm(x))
        if value == 0.0:
            return 1.0
        x /= value
    return value


def field_map(echo1: np.ndarray, echo2: np.ndarray, *, delta_te_s: float,
              mask: np.ndarray | None = None) -> np.ndarray:
    """
    ``angle(S1 conj(S2)) / (2 pi dTE)`` in Hz, with the wrap check as an assertion.

    **The sign is the operator's**, not a convention chosen here.  This module's forward model
    accumulates ``exp(-2i pi df t)``, so a spin at ``+df`` arrives at the second echo with
    ``S2 = S1 exp(-2i pi df dTE)`` and the difference has to be read the other way round from the
    one the expression is usually written with.  Getting it backwards produces a map that is the
    exact negation of the field, which is the input that makes an off-resonance correction *worse*
    than doing nothing -- and it is invisible in the map itself, which looks like a plausible
    shim either way.

    A two-point map is unambiguous only while ``|df| * dTE < 0.5``; past that it folds, and a
    folded field map is a *correction* that makes the image worse in exactly the region that
    needed it most.  So the bound is checked inside the mask rather than mentioned in a comment.
    """
    difference = np.angle(np.asarray(echo1) * np.conj(np.asarray(echo2)))
    field = difference / (2 * np.pi * float(delta_te_s))
    inside = np.ones(field.shape, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    worst = float(np.abs(difference[inside]).max() / (2 * np.pi)) if inside.any() else 0.0
    if worst > 0.45:
        limit = 0.5 / float(delta_te_s)
        msg = (f'the phase difference reaches {worst:.3f} cycles inside the mask, so this map is '
               f'at or past its +/-{limit:.0f} Hz wrap. Shorten delta_te_s, or unwrap first.')
        raise ValueError(msg)
    return field
