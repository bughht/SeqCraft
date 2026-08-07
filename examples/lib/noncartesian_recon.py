"""
Off-resonance-corrected reconstruction for readouts whose samples are not on a uniform grid.

The problem is not the trajectory, it is **time**.  Every sample is acquired at a different moment
relative to the echo, so a voxel off resonance by ``df`` accumulates a phase ``2*pi*df*t`` that
varies along the readout.  Where that phase lands decides what the artefact looks like, and the two
readouts here land in different places from the *same* term:

* a **spiral** sweeps k continuously, so the phase spreads over the whole plane and the image
  **blurs** -- worse the longer the readout, which is the direction a diffusion sequence is pushed in
  for SNR;
* an **EPI** train advances one ``ky`` step per echo spacing, so the phase is linear in ``ky`` and the
  image **shifts** along the phase-encode direction -- ``df`` divided by the phase-encode bandwidth
  per pixel, which is five pixels per 100 Hz on the sequence the EPI example builds.

One operator therefore serves both, and no deblurring or unwarping step is needed anywhere: the
forward model is exact, with no segmentation or time-binning,

.. code-block:: text

    y_c(t) = sum_r  s_c(r) x(r) exp(-i phi(r, t))
    phi(r, t) = 2 pi ( k(t) . r  +  df(r) t )  [+ phi_extra(t)]

so correcting off-resonance is a matter of putting ``df`` into the operator and solving.
``phi_extra`` is a spatially uniform, time-varying phase, which is where a concomitant-field term
goes for an off-isocentre slice.

The price is that this is a dense operator: ``n_samples x n_pixels`` complex exponentials per
application.  Being exact makes it the right thing to *check* a sequence against, because any
residual blur is then the sequence's and not the reconstruction's -- and it makes it useless for a
real one.  The spiral DTI sequence has 59 100 samples at 128 x 128, which is 968 million
exponentials per matrix-vector product: about 25 s, so 2.8 hours for one conjugate-gradient solve.

:class:`SegmentedOffresonance` is the version that scales, approximating the off-resonance term
between a few instants in time and leaving the rest to an ordinary NUFFT -- 26 s for the same solve.
Approximations have to be checked rather than trusted, so it is validated against the exact operator
at a matrix where both are affordable; :func:`reconstruct` reaches it through ``segments=``.

Lives in ``examples/lib`` rather than inside seqcraft: the package builds sequences, and
reconstruction is a downstream concern with its own dependency (sigpy).

Usage
-----
::

    import sys; sys.path.insert(0, '../lib')
    from noncartesian_recon import Readout, reconstruct_shot

    readout = Readout.from_sidecar('seq/epi_dwi.traj.npz')
    image = reconstruct_shot(data, readout, b0_hz=field_map, sens=coil_maps, dcf='speed')

`data` is ``(n_coils, n_samples)`` however it was obtained -- simulated, or read from a vendor raw
file -- which is the whole point of the split: the reconstruction never learns where it came from.
:func:`reconstruct` underneath takes the trajectory as plain arrays if the sidecar is not involved.

Two things differ per readout, and both are explicit rather than inferred, because inferring them
from a spiral is what made them wrong for anything else:

* **where the echo is** -- :attr:`Readout.t_echo_s`, since a spiral's echo is its first sample and an
  EPI train's is a third of the way in;
* **which density compensation** -- ``dcf='radial'`` for a spiral, ``'speed'`` for a ramp-sampled
  EPI, whose samples are *sparsest* where k = 0.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from seqcraft.core.errors import MissingExtraError, format_error

if TYPE_CHECKING:
    from collections.abc import Sequence

__all__ = [
    'NoncartesianOffresonance',
    'Readout',
    'SegmentedOffresonance',
    'density_compensation',
    'reconstruct',
    'reconstruct_shot',
    'reconstruct_volume',
    'segments_for',
    'speed_compensation',
]


def _sigpy() -> Any:
    """Import sigpy, with a message naming the extra to install."""
    try:
        import sigpy  # noqa: PLC0415
    except ImportError as err:  # pragma: no cover - depends on the environment
        msg = format_error(
            'reconstruction needs sigpy.',
            {'missing': 'sigpy'},
            ['pip install "seqcraft[recon]"'],
        )
        raise MissingExtraError(msg) from err
    return sigpy


def segments_for(
    b0_hz: np.ndarray | float,
    sample_times_s: np.ndarray,
    *,
    target_error: float = 2e-2,
) -> int:
    """
    Return the number of time segments needed for a given operator accuracy.

    Parameters
    ----------
    b0_hz
        The off-resonance map, in hertz.  Only its **range** matters, not its absolute values:
        :class:`SegmentedOffresonance` removes the mean exactly, so what is left to interpolate spans
        plus or minus half the range.
    sample_times_s
        Sample times; only the span is used.
    target_error
        Worst-case relative error in the encoding operator, over the field's whole range.  The default
        of 2 % sounds loose and is not: it is a bound at the single worst pixel and instant, and the
        error it produces in an *image* is several times smaller -- measured below.

    Notes
    -----
    Time segmentation interpolates ``exp(-2i pi df t)`` linearly between segment centres, and linear
    interpolation of a unit-modulus exponential over a phase step ``dtheta`` is wrong by about
    ``dtheta^2 / 8``.  Inverting that for ``L = 2 pi df_halfrange T / dtheta`` gives this count.

Measured on the spiral DTI sequence -- a 59 ms readout and a field spanning -60 to +112 Hz,
    reconstructed at 128 x 128 with the B0 map properly interpolated:

    .. code-block:: text

        L      dtheta   bound     seconds   RMSE against the reference
        16       1.99   0.50         27       0.111
        34       0.94   0.11         41       0.024
        67       0.48   0.028        68       0.0155
       134       0.24   0.0071      122       0.0151
       268       0.12   0.0018      229       0.0150

    The knee is near ``L = 67``: past it the RMSE moves by 3 % while the cost triples, because
    something other than segmentation sets the floor.  The default target lands just beyond that knee.

    Do not pick ``L`` by eye from an image, and do not pick it by a round number of residual cycles.
    The old rule here aimed at 0.1 cycles per segment, which is ``dtheta = 0.63`` and a **5 % bound**;
    that was genuinely too coarse, and it was also not what caused the artefact it was blamed for.

    A visible ripple in a corrected image has at least three candidate sources, and **which one
    dominates depends on the other two** -- so diagnosing them one at a time finds whichever was looked
    at first:

    - *this* approximation, bounded above and cheap to rule out by making ``L`` finer;
    - the **B0 map's own resampling**, which at a non-integer grid ratio is a moire generator rather
      than a smooth error -- 160 onto 128 cost 0.215 cycles of patterned phase on this sequence,
      against segmentation's 0.028 bound;
    - **solver convergence**, because correcting worsens the operator's conditioning: at 20 conjugate-
      gradient iterations the corrected image carried *more* high-frequency texture than the
      uncorrected one and the map error was worth 9 %, while at 60 the solver contributed nothing and
      the same map error was worth a factor of two and a half.

    Vary all three together before concluding anything, and prefer RMSE against a reference
    reconstruction over looking at the image: structured residuals of very different origin look alike.

    Examples
    --------
    >>> import numpy as np
    >>> times = np.linspace(0.0, 59.1e-3, 1000)
    >>> field = np.array([[-60.0, 112.0]])              # a 172 Hz range about a 26 Hz mean
    >>> segments_for(field, times)
    80
    >>> segments_for(field, times, target_error=5e-2)   # the old rule, in effect
    51
    >>> segments_for(np.zeros((4, 4)), times)           # nothing to correct, nothing to segment
    8
    """
    times = np.asarray(sample_times_s, dtype=np.float64).reshape(-1)
    span = float(times.max() - times.min()) if times.size else 0.0
    field = np.asarray(b0_hz, dtype=np.float64)
    half_range = 0.5 * float(field.max() - field.min()) if field.size else 0.0
    if span <= 0.0 or half_range <= 0.0 or target_error <= 0.0:
        return 8
    step = math.sqrt(8.0 * float(target_error))
    return max(8, int(math.ceil(2.0 * math.pi * half_range * span / step)))


@dataclass(frozen=True)
class Readout:
    """
    What a non-Cartesian reconstruction needs to know about the readout, and nothing else.

    Attributes
    ----------
    k_per_m
        ``(n_interleaves, 2, n_samples)`` trajectory as ``(kx, ky)`` in 1/m, at **ADC sample times**.
    t_adc_s
        ``(n_samples,)`` sample times **relative to the echo**, so they are negative before it.  For
        a spiral that means they start at (near) zero, since a spiral begins at k = 0; for an EPI
        train the echo is a third of the way in and the first third of the array is negative.
    fov_m, matrix
        Geometry, so the reconstruction needs nothing else passed alongside.
    t_echo_s
        Where the echo sits in `t_adc_s`, for a sidecar that recorded the times from the readout
        block's start rather than from the echo.  Subtracted on the way through
        :meth:`interleaf`; leave it at zero when the times are already echo-relative.

    Notes
    -----
    This is the whole interface between a scan and a reconstruction, and it is deliberately four
    arrays rather than a sequence object: what a reconstruction needs is the *trajectory*, and where
    it came from -- seqcraft's sidecar, a vendor trajectory file, a field-camera measurement -- is not
    its business.

    **The echo reference is explicit, and it used to be inferred.**  This class rebased the times on
    ``min(t_adc_s)``, which is right for a spiral for one reason only: a spiral's first sample *is*
    its echo.  An EPI train reaches ``ky = 0`` 17 ms into a 50 ms readout, and rebasing on its first
    sample puts ``2 pi df * 17 ms`` of phase on the operator -- **spatially varying**, since ``df``
    varies with position, so not a constant that comes out in the wash.  On a 172 Hz field range
    that is up to 2.9 cycles of unmodelled phase.

    Orientation, stated once so it cannot drift.  Trajectory in as ``(kx, ky)``; image out as
    ``[y, x]``, rows increasing with y, which is what ``imshow`` and every downstream tool assume.
    :func:`reconstruct` pairs them accordingly and there is nothing to transpose at the call site.

    Off-resonance sign, likewise: ``exp(-2i pi (k.r + df t))``, both terms negative, which is the
    standard convention and what a scanner delivers. A map from a vendor field-map sequence goes in
    as it comes out.
    """

    k_per_m: np.ndarray
    t_adc_s: np.ndarray
    fov_m: float
    matrix: int
    t_echo_s: float = 0.0

    def __post_init__(self) -> None:
        """
        Reject times that look absolute rather than readout-relative.

        The failure this catches is a sidecar carrying sample times measured from the start of the
        *sequence*, which for a spin echo includes TE.  The operator would then apply
        ``2 pi df * (TE + t)`` instead of ``2 pi df * t``, and since ``df`` varies with position the
        surplus is an image artefact rather than a global phase -- structured, plausible-looking, and
        not obviously a units problem.

        The discriminator is that ``t_echo_s`` claims the echo is at the very start while the times
        do not begin anywhere near it.
        """
        times = np.asarray(self.t_adc_s, dtype=np.float64).reshape(-1)
        if times.size < 2:
            return
        span = float(times.max() - times.min())
        lead = float(times.min() - self.t_echo_s)
        if span > 0.0 and lead > 0.05 * span:
            msg = format_error(
                'the sample times start well after the echo, so they look absolute rather than '
                'measured from the readout.',
                {
                    'first_sample_ms': times.min() * 1e3,
                    't_echo_s_given_ms': self.t_echo_s * 1e3,
                    'readout_span_ms': span * 1e3,
                },
                [
                    'pass t_echo_s= the time of k=0 within the readout '
                    '(readout.time_to_echo for a seqcraft module)',
                    'or subtract the readout block start from t_adc_s before constructing this',
                ],
            )
            raise ValueError(msg)

    @classmethod
    def from_sidecar(cls, path: Any) -> Readout:
        """
        Load from the ``.traj.npz`` written beside a ``.seq``.

        Reads ``t_echo_s`` when the sidecar records one, so a readout whose echo is not its first
        sample describes itself rather than relying on a convention at the call site.

        Examples
        --------
        ::

            readout = Readout.from_sidecar('seq/epi_dwi.traj.npz')
            image = reconstruct_shot(twix_data, readout, b0_hz=field_map, sens=coil_maps)
        """
        with np.load(str(path)) as sidecar:
            trajectory = np.asarray(sidecar['k_per_m'], dtype=np.float64)
            if 't_adc_s' in sidecar:
                times = np.asarray(sidecar['t_adc_s'], dtype=np.float64)
            else:
                # Older sidecars stored only the dwell; the samples are uniform, so this is exact.
                dwell = float(sidecar['dwell_s'])
                delay = float(sidecar.get('adc_delay_s', 0.0))
                times = delay + (np.arange(trajectory.shape[-1]) + 0.5) * dwell
            return cls(
                k_per_m=trajectory,
                t_adc_s=times,
                fov_m=float(sidecar['fov_m']),
                matrix=int(sidecar['matrix']),
                t_echo_s=float(sidecar['t_echo_s']) if 't_echo_s' in sidecar else 0.0,
            )

    @property
    def n_interleaves(self) -> int:
        """Number of rotated shots covering k-space together."""
        return int(self.k_per_m.shape[0])

    @property
    def n_samples(self) -> int:
        """ADC samples per interleaf."""
        return int(self.k_per_m.shape[-1])

    def interleaf(self, index: int = 0) -> tuple[np.ndarray, np.ndarray]:
        """
        Return ``(k_per_m, t_relative_to_echo_s)`` for one interleaf.

        Times come back measured from the **echo**, which is what the operator wants: an absolute
        time would carry TE into the off-resonance term and put ``2*pi*df*TE`` of phase on
        everything -- and since ``df`` varies with position, that phase is an image artefact rather
        than a global constant.

        The reference is :attr:`t_echo_s`, taken from the sidecar.  It is *not* inferred from the
        first sample: that inference is correct for a spiral and wrong for anything whose echo falls
        mid-readout.
        """
        return (
            np.asarray(self.k_per_m[index % self.n_interleaves], dtype=np.float64),
            np.asarray(self.t_adc_s, dtype=np.float64) - float(self.t_echo_s),
        )


def reconstruct_shot(
    data: np.ndarray,
    readout: Readout,
    *,
    interleaf: int = 0,
    b0_hz: np.ndarray | float = 0.0,
    sens: np.ndarray | None = None,
    matrix: int | None = None,
    dcf: str = 'radial',
    iterations: int = 80,
    lamda: float = 1e-3,
    segments: int | None = 0,
    device: int = -1,
    progress: bool = False,
) -> np.ndarray:
    """
    Reconstruct one shot from its raw channels and a :class:`Readout`.

    Parameters
    ----------
    data
        ``(n_coils, n_samples)`` complex samples, in acquisition order.  This is what a vendor raw
        reader hands over: with ``twixtools``, one ``mdb.data`` per acquisition; with ``mapVBVD``, one
        column of the image data array.
    readout
        From :meth:`Readout.from_sidecar`.
    interleaf
        Which rotated shot this data belongs to.  Single-shot sequences only have interleaf 0; for a
        multi-shot spiral it comes from the ``SEG`` label in the ``.seq``.
    b0_hz
        ``(matrix, matrix)`` off-resonance map in hertz on the reconstruction grid, or ``0.0`` for no
        correction.  Sign convention as :class:`Readout` documents.
    sens
        ``(n_coils, matrix, matrix)`` coil sensitivities -- ESPIRiT maps from a reference scan.
        Defaults to a single uniform channel, which is correct only if the data is fully sampled.
    matrix
        Overrides the sidecar's matrix, to reconstruct coarser or finer than nominal.
    dcf, iterations, lamda, segments, device, progress
        As :func:`reconstruct`.  ``dcf='speed'`` for a ramp-sampled EPI.

    Returns
    -------
    numpy.ndarray
        ``(matrix, matrix)`` complex image, indexed ``[y, x]``.

    Notes
    -----
    Deliberately thin: it selects the interleaf, rebases the sample times on the echo, and calls
    :func:`reconstruct`.  Everything specific to *where the data came from* stays outside, which is
    what makes the same call work for simulated and scanner data.

    What still has to be got right for real raw data, none of which this function can check:

    - **Channel order** must match `sens`, and a vendor reader may reorder or compress channels.
    - **Samples per acquisition** must match the trajectory.  Vendors often prepend dummy or
      oversampled points; the count disagreeing with ``readout.n_samples`` is the first thing to look
      at when an image comes out as noise.
    - **Gradient delay** shifts the trajectory against the data by a microsecond or two, which a
      spiral shows as a ring or a shading and which no B0 map explains.  It needs measuring, and
      then ``t_adc_s`` shifting to match.
    """
    channels = np.atleast_2d(np.asarray(data, dtype=np.complex64))
    if channels.shape[0] > channels.shape[1]:          # (n_samples, n_coils) given
        channels = channels.T
    k, times = readout.interleaf(interleaf)
    if channels.shape[1] != times.size:
        msg = format_error(
            'the data and the trajectory disagree about how many samples there are.',
            {'data': channels.shape[1], 'trajectory': times.size},
            [
                'check for vendor dummy or oversampled points at the start of the acquisition',
                'check that this is one shot rather than a concatenation of several',
            ],
        )
        raise ValueError(msg)
    return reconstruct(
        channels, k, times,
        fov_m=readout.fov_m,
        matrix=int(matrix if matrix is not None else readout.matrix),
        b0_hz=b0_hz,
        sens=sens,
        n_per_arm=readout.n_samples,
        dcf=dcf,
        iterations=iterations,
        lamda=lamda,
        segments=segments,
        device=device,
        progress=progress,
    )


def reconstruct_volume(
    shots: Any,
    readout: Readout,
    *,
    b0_hz: np.ndarray | float = 0.0,
    sens: np.ndarray | None = None,
    matrix: int | None = None,
    dcf: str = 'radial',
    iterations: int = 80,
    lamda: float = 1e-3,
    segments: int | None = 0,
    device: int = -1,
    progress: bool = False,
) -> np.ndarray:
    """
    Reconstruct one image from every interleaf of a multi-shot spiral.

    Parameters
    ----------
    shots
        The interleaves' data, either a sequence in interleaf order or a mapping from interleaf index
        to ``(n_coils, n_samples)``.  On real data the interleaf comes from the ``SEG`` label.
    readout, b0_hz, sens, matrix, dcf, iterations, lamda, segments, device, progress
        As :func:`reconstruct_shot`.

    Returns
    -------
    numpy.ndarray
        ``(matrix, matrix)`` complex image, indexed ``[y, x]``.

    Notes
    -----
    All interleaves are solved **together**, as one inverse problem over their pooled samples.  That
    is not the same as reconstructing each and averaging: individually each interleaf is severely
    under-determined -- a quarter of k-space for a four-shot spiral -- and averaging four aliased
    images leaves the aliasing.  Pooled, they are exactly the fully sampled problem the sequence was
    designed to be.

    The samples are **sorted by time**, which matters for a reason that is easy to miss.  Every
    interleaf is sampled at the same times relative to its own echo, so concatenating them gives a
    time axis that resets at each interleaf boundary; :class:`SegmentedOffresonance` needs it
    non-decreasing, since it slices contiguous ranges out of it.  Sorting is safe because the
    operator only ever sees ``(k, t)`` pairs and never depends on acquisition order -- as long as the
    data, the trajectory and the times are permuted together, which is what happens here.
    """
    ordered = (
        [np.asarray(shots[index]) for index in sorted(shots)] if hasattr(shots, 'keys')
        else [np.asarray(shot) for shot in shots]
    )
    if not ordered:
        msg = format_error('reconstruct_volume needs at least one interleaf.', {})
        raise ValueError(msg)
    if len(ordered) > readout.n_interleaves:
        msg = format_error(
            'more interleaves given than the trajectory describes.',
            {'given': len(ordered), 'trajectory': readout.n_interleaves},
        )
        raise ValueError(msg)

    channels: list[np.ndarray] = []
    coordinates: list[np.ndarray] = []
    times: list[np.ndarray] = []
    for index, shot in enumerate(ordered):
        block = np.atleast_2d(np.asarray(shot, dtype=np.complex64))
        if block.shape[0] > block.shape[1]:            # (n_samples, n_coils) given
            block = block.T
        k, t = readout.interleaf(index)
        if block.shape[1] != t.size:
            msg = format_error(
                f'interleaf {index}: the data and the trajectory disagree on sample count.',
                {'data': block.shape[1], 'trajectory': t.size},
            )
            raise ValueError(msg)
        channels.append(block)
        coordinates.append(k)
        times.append(t)

    pooled_t = np.concatenate(times)
    order = np.argsort(pooled_t, kind='stable')
    return reconstruct(
        np.concatenate(channels, axis=1)[:, order],
        np.concatenate(coordinates, axis=1)[:, order],
        pooled_t[order],
        fov_m=readout.fov_m,
        matrix=int(matrix if matrix is not None else readout.matrix),
        b0_hz=b0_hz,
        sens=sens,
        n_per_arm=readout.n_samples,
        dcf=dcf,
        iterations=iterations,
        lamda=lamda,
        segments=segments,
        device=device,
        progress=progress,
    )


def density_compensation(k_per_m: np.ndarray, *, n_per_arm: int | None = None) -> np.ndarray:
    """
    Return a ``|k|``-based density-compensation weight per sample, normalised to mean one.

    Parameters
    ----------
    k_per_m
        ``(2, n_samples)`` or ``(3, n_samples)`` k-space coordinates in 1/m.
    n_per_arm
        Samples per interleaf, used only to set the floor near the origin.

    Returns
    -------
    numpy.ndarray
        ``(n_samples,)`` weights.

    Notes
    -----
    A spiral samples the centre of k-space far more densely than the edge, so an unweighted adjoint
    is dominated by the centre and looks heavily low-pass filtered.  Weighting by ``|k|`` is the
    cheap approximation to the true Jacobian; it is not exact for a variable-density spiral, which
    is one reason to solve iteratively rather than trust a single adjoint.

    **For a spiral only.**  It is a radial argument, and it is backwards for any trajectory that does
    not sample its centre densely: a ramp-sampled EPI line passes ``k = 0`` at full gradient, which
    is where its samples are *sparsest*, so ``|k|`` down-weights exactly the samples carrying the
    signal.  Use :func:`speed_compensation` there -- and for a spiral too, where it agrees with this
    to within the variable-density part it cannot see.
    """
    radius = np.sqrt(np.sum(np.asarray(k_per_m[:2], dtype=float) ** 2, axis=0))
    floor = 0.5 * float(radius.max()) / max(int(n_per_arm or radius.size), 1)
    weights = radius + floor
    return (weights / float(np.mean(weights))).astype(np.float32)


def speed_compensation(
    k_per_m: np.ndarray,
    sample_times_s: np.ndarray,
    *,
    floor: float = 1e-3,
) -> np.ndarray:
    """
    Return ``|dk/dt|``-based weights per sample, normalised to mean one.

    Parameters
    ----------
    k_per_m
        ``(2, n_samples)`` or ``(3, n_samples)`` k-space coordinates in 1/m, in acquisition order.
    sample_times_s
        ``(n_samples,)`` sample times.  Only the differences are used.
    floor
        Fraction of the median speed to clamp at, so a stationary instant -- an EPI blip, where the
        readout gradient is at zero -- does not get a zero weight and drop its sample entirely.

    Returns
    -------
    numpy.ndarray
        ``(n_samples,)`` weights.

    Notes
    -----
    The correct Jacobian for a trajectory traversed at varying speed, and correct for **both** the
    readouts here rather than one.  Sampling uniformly in time means the local sample density is
    inversely proportional to how fast k is moving, so the compensating weight is the speed itself.

    For a ramp-sampled EPI line that is the whole story: the gradient is a triangle, so samples pile
    up at the two ends of every line where k crawls and thin out at the apex where k = 0 is crossed
    at full speed.  ``|k|``-weighting gets this exactly the wrong way round -- it down-weights the
    centre of k-space, which is where the signal is.

    For a spiral it reproduces the ``|k|`` argument in the amplitude-limited part, where speed is
    constant and the density falls as ``1/|k|``, and does better near the origin and in the
    variable-density taper, where ``|k|`` is only an approximation.

    Speeds are taken from central differences, so a sample's weight reflects the interval it
    actually represents rather than the one after it.  Across a shot boundary in a pooled multi-shot
    problem the difference is meaningless, which the clamp absorbs.
    """
    k = np.asarray(k_per_m, dtype=np.float64)
    t = np.asarray(sample_times_s, dtype=np.float64).reshape(-1)
    if k.shape[-1] != t.size:
        msg = format_error(
            'speed_compensation needs one time per k-space sample.',
            {'k_samples': k.shape[-1], 'times': t.size},
        )
        raise ValueError(msg)
    if t.size < 2:
        return np.ones(t.size, dtype=np.float32)

    distance = np.linalg.norm(np.diff(k, axis=-1), axis=0)
    interval = np.abs(np.diff(t))
    segment_speed = np.divide(
        distance, interval, out=np.zeros_like(distance), where=interval > 0
    )
    # Central difference: average the two segments a sample sits between, and take the single
    # neighbouring segment at each end.
    speed = np.empty(t.size)
    speed[0] = segment_speed[0]
    speed[-1] = segment_speed[-1]
    speed[1:-1] = 0.5 * (segment_speed[:-1] + segment_speed[1:])

    median = float(np.median(speed[speed > 0])) if np.any(speed > 0) else 1.0
    weights = np.maximum(speed, floor * median)
    return (weights / float(np.mean(weights))).astype(np.float32)


class NoncartesianOffresonance:
    """
    Exact phase-accumulation encoding operator, as a sigpy ``Linop``.

    Parameters
    ----------
    matrix
        ``(ny, nx)`` image shape.
    sens
        ``(n_coils, ny, nx)`` complex coil sensitivities.  Pass a single ones-array for a
        single-channel simulation.
    coord_cycles_per_fov
        ``(n_samples, 2)`` k-space coordinates as ``(ky, kx)`` in **cycles per FOV** -- that is,
        ``k[1/m] * FOV[m]``.  Cycles per FOV rather than 1/m because the spatial grid is then the
        unit square and the operator needs no further knowledge of the geometry.
    sample_times_s
        ``(n_samples,)`` time of each sample **relative to the echo**.  Relative, because the phase
        that matters is the one accumulated since the echo, and an absolute time would put a large
        constant phase on everything.
    b0_hz
        ``(ny, nx)`` off-resonance map in hertz.  Zeros gives the uncorrected reconstruction, which
        is the useful comparison.
    weights_sqrt
        Optional ``(n_samples,)`` square root of the density compensation, applied to both the
        forward and adjoint so the normal equations stay symmetric.
    extra_phase_rad
        Optional ``(n_samples,)`` spatially uniform phase -- a concomitant-field term, for instance.
    device
        sigpy device index; ``-1`` for CPU.
    time_batch
        Samples per batch.  Trades memory for speed; the operator materialises
        ``time_batch x n_pixels`` complex numbers at a time.

    Notes
    -----
    Not a subclass of ``sigpy.linop.Linop`` at class-definition time, because sigpy is an optional
    dependency and importing it at module scope would make ``import seqcraft`` require it.
    :meth:`linop` builds the real thing on demand.
    """

    def __init__(
        self,
        *,
        matrix: tuple[int, int],
        sens: np.ndarray,
        coord_cycles_per_fov: np.ndarray,
        sample_times_s: np.ndarray,
        b0_hz: np.ndarray,
        weights_sqrt: np.ndarray | None = None,
        extra_phase_rad: np.ndarray | None = None,
        device: int = -1,
        time_batch: int = 4096,
    ) -> None:
        self.ny, self.nx = int(matrix[0]), int(matrix[1])
        self.n_pixels = self.ny * self.nx
        self.sens = np.asarray(sens, dtype=np.complex64).reshape(-1, self.n_pixels)
        self.n_coils = self.sens.shape[0]
        self.coord = np.asarray(coord_cycles_per_fov, dtype=np.float32)
        self.times = np.asarray(sample_times_s, dtype=np.float32).reshape(-1)
        self.b0 = np.asarray(b0_hz, dtype=np.float32).reshape(self.n_pixels)
        self.weights_sqrt = (
            None if weights_sqrt is None
            else np.asarray(weights_sqrt, dtype=np.float32).reshape(-1)
        )
        self.extra_phase = (
            None if extra_phase_rad is None
            else np.asarray(extra_phase_rad, dtype=np.float32).reshape(-1)
        )
        self.device = int(device)
        self.time_batch = max(64, int(time_batch))

        if self.coord.ndim != 2 or self.coord.shape[1] != 2:
            msg = format_error(
                'coord must be (n_samples, 2) as (ky, kx).', {'given': str(self.coord.shape)}
            )
            raise ValueError(msg)
        self.n_samples = self.coord.shape[0]
        for name, array in (('sample_times_s', self.times), ('weights_sqrt', self.weights_sqrt),
                            ('extra_phase_rad', self.extra_phase)):
            if array is not None and array.shape[0] != self.n_samples:
                msg = format_error(
                    f'{name} must have one entry per sample.',
                    {'given': array.shape[0], 'samples': self.n_samples},
                )
                raise ValueError(msg)

    def linop(self) -> Any:
        """Return the sigpy ``Linop`` implementing this operator."""
        sp = _sigpy()
        outer = self

        # Spatial grid on the unit square, matching cycles-per-FOV coordinates.
        ry = (np.arange(outer.ny, dtype=np.float32) - outer.ny / 2.0) / float(outer.ny)
        rx = (np.arange(outer.nx, dtype=np.float32) - outer.nx / 2.0) / float(outer.nx)
        grid_y, grid_x = np.meshgrid(ry, rx, indexing='ij')

        device = sp.Device(outer.device)
        state = {
            'sens': sp.to_device(outer.sens, device),
            'coord': sp.to_device(outer.coord, device),
            'times': sp.to_device(outer.times, device),
            'b0': sp.to_device(outer.b0, device),
            'ry': sp.to_device(grid_y.reshape(outer.n_pixels), device),
            'rx': sp.to_device(grid_x.reshape(outer.n_pixels), device),
            'weights': None if outer.weights_sqrt is None
            else sp.to_device(outer.weights_sqrt, device),
            'forward_phase': None if outer.extra_phase is None
            else sp.to_device(np.exp(1j * outer.extra_phase).astype(np.complex64), device),
            'adjoint_phase': None if outer.extra_phase is None
            else sp.to_device(np.exp(-1j * outer.extra_phase).astype(np.complex64), device),
        }

        def phase_block(lo: int, hi: int) -> Any:
            xp = device.xp
            coord = state['coord'][lo:hi]
            times = state['times'][lo:hi]
            gradient = 2.0 * xp.pi * (
                xp.outer(coord[:, 0], state['ry']) + xp.outer(coord[:, 1], state['rx'])
            )
            offresonance = 2.0 * xp.pi * xp.outer(times, state['b0'])
            # Off-resonance carries the *same* sign as k.r, which is the standard MRI convention and
            # what a scanner delivers: s(t) = sum_r m(r) exp(-2i pi (k.r + df t)).  Written with a
            # minus here instead, the correction adds exactly the phase it is meant to remove, and
            # the result looks like "off-resonance correction does not help" rather than like a sign
            # error -- so it survives being tested.
            return (gradient + offresonance).astype(xp.float32, copy=False)

        class Forward(sp.linop.Linop):
            """Image to k-space."""

            def __init__(self) -> None:
                super().__init__((outer.n_coils, outer.n_samples), (outer.ny, outer.nx))

            def _apply(self, x: Any) -> Any:
                xp = device.xp
                x = sp.to_device(x, device).reshape(outer.n_pixels).astype(xp.complex64, copy=False)
                y = xp.zeros((outer.n_coils, outer.n_samples), dtype=xp.complex64)
                weighted = state['sens'] * x[None, :]
                for lo in range(0, outer.n_samples, outer.time_batch):
                    hi = min(lo + outer.time_batch, outer.n_samples)
                    kernel = xp.exp(-1j * phase_block(lo, hi)).astype(xp.complex64, copy=False)
                    y[:, lo:hi] = (kernel @ weighted.T).T
                if state['weights'] is not None:
                    y = y * state['weights'][None, :]
                if state['forward_phase'] is not None:
                    y = y * state['forward_phase'][None, :]
                return y

            def _adjoint_linop(self) -> Any:
                return Adjoint()

        class Adjoint(sp.linop.Linop):
            """k-space to image."""

            def __init__(self) -> None:
                super().__init__((outer.ny, outer.nx), (outer.n_coils, outer.n_samples))

            def _apply(self, y: Any) -> Any:
                xp = device.xp
                y = sp.to_device(y, device).astype(xp.complex64, copy=False)
                if state['weights'] is not None:
                    y = y * state['weights'][None, :]
                if state['adjoint_phase'] is not None:
                    y = y * state['adjoint_phase'][None, :]
                image = xp.zeros(outer.n_pixels, dtype=xp.complex64)
                conjugate = xp.conj(state['sens'])
                for lo in range(0, outer.n_samples, outer.time_batch):
                    hi = min(lo + outer.time_batch, outer.n_samples)
                    kernel = xp.exp(1j * phase_block(lo, hi)).astype(xp.complex64, copy=False)
                    coils = (kernel.T @ y[:, lo:hi].T).T
                    image += xp.sum(conjugate * coils, axis=0)
                return image.reshape(outer.ny, outer.nx)

            def _adjoint_linop(self) -> Any:
                return Forward()

        return Forward()


class SegmentedOffresonance:
    """
    Time-segmented off-resonance operator: a few NUFFTs instead of a dense phase matrix.

    Parameters
    ----------
    matrix
        ``(ny, nx)`` image shape.
    sens
        ``(n_coils, ny, nx)`` complex coil sensitivities.
    coord_cycles_per_fov
        ``(n_samples, 2)`` k-space coordinates as ``(ky, kx)`` in cycles per FOV.  This is exactly
        sigpy's NUFFT convention, and sigpy's kernel is ``exp(-2i pi k.r)`` -- the same sign as
        :class:`NoncartesianOffresonance`, so the two operators are directly comparable.
    sample_times_s
        ``(n_samples,)`` time of each sample **relative to the echo**, non-decreasing.
    b0_hz
        ``(ny, nx)`` off-resonance map in hertz.
    segments
        Number of time segments.  Needs to be large enough that the residual phase *within* a
        segment is small: roughly ``4 * max|df| * T`` is comfortable, and
        :func:`reconstruct` defaults to that.
    weights_sqrt, extra_phase_rad, device
        As :class:`NoncartesianOffresonance`.

    Notes
    -----
    The exact operator is ``n_samples x n_pixels`` complex exponentials per application, which is a
    fine way to *verify* a reconstruction and hopeless as a way to perform one: the real sequence in
    the spiral build notebook has 59 100 samples at 128 x 128, which is 968 million exponentials per matrix-vector
    product -- about 25 s, so 2.8 hours for one conjugate-gradient solve.

    Time segmentation replaces it with the standard approximation.  The off-resonance term is the
    only part of the phase that is not a Fourier transform, and it varies *slowly* in time, so
    interpolate it between a few instants:

    .. code-block:: text

        exp(-2i pi df(r) t)  ~  sum_l  w_l(t) exp(-2i pi df(r) t_l)

        y_c(t) = sum_l w_l(t) . NUFFT[ s_c x exp(-2i pi df t_l) ](t)

    with ``w_l`` the linear hat function on the segment centres.  Each segment is then an ordinary
    NUFFT, which is ``O(N log N)`` in the image and linear in the samples.

    Cost grows only weakly with the number of segments.  With linear interpolation each sample has
    non-zero weight in exactly two segments, and for a spiral the samples are time-ordered, so each
    segment's support is a contiguous slice of about ``2/L`` of the readout -- the *interpolation*
    work is therefore independent of ``L``.  What does scale with ``L`` is one oversampled FFT of the
    image per segment, which is why going from 8 segments to 96 costs about 2.7x rather than 12x.

    Accuracy is set by the phase left over inside a segment and the error falls as roughly its
    square, so :func:`segments_for` derives ``L`` from an explicit accuracy target rather than from a
    round number of residual cycles.  Measured against the exact operator on the spiral DTI sequence
    (59 ms readout, a field spanning -60 to +112 Hz):

    .. code-block:: text

        L      dtheta      bound      measured RMSE
         8       5.20       3.4           0.85
        16       2.60       0.84          0.25
        32       1.30       0.21          0.046
        64       0.65       0.053         0.012
        80       0.39       0.019         -- the default target

    The bound is worst case over the field's whole range while the RMSE averages over an image whose
    pixels mostly sit well inside it, which is why RMSE runs several times better -- and why RMSE is
    the wrong thing to choose ``L`` by.  At ``L = 64`` the operator is 5 % wrong somewhere, and being
    structured rather than random that shows up as a visible ripple even though the RMSE reads 0.012.

    Being an approximation, it has to be checked rather than trusted, which is what the exact
    operator is for: the reconstruction notebook reconstructs the same data both ways at a matrix where both are
    affordable and compares.
    """

    def __init__(
        self,
        *,
        matrix: tuple[int, int],
        sens: np.ndarray,
        coord_cycles_per_fov: np.ndarray,
        sample_times_s: np.ndarray,
        b0_hz: np.ndarray,
        segments: int = 16,
        weights_sqrt: np.ndarray | None = None,
        extra_phase_rad: np.ndarray | None = None,
        device: int = -1,
    ) -> None:
        self.ny, self.nx = int(matrix[0]), int(matrix[1])
        self.sens = np.asarray(sens, dtype=np.complex64).reshape(-1, self.ny, self.nx)
        self.n_coils = self.sens.shape[0]
        self.coord = np.asarray(coord_cycles_per_fov, dtype=np.float32)
        self.times = np.asarray(sample_times_s, dtype=np.float64).reshape(-1)
        field = np.asarray(b0_hz, dtype=np.float32).reshape(self.ny, self.nx)
        # Centre the field before segmenting, and carry its mean exactly.
        #
        # Accuracy is set by the phase left over *inside* a segment, so it depends on how far the
        # field ranges from the demodulation reference -- and the choice of reference is arbitrary.
        # Left uncentred, a map spanning -60 to +112 Hz is segmented as though it spanned +-112,
        # needing 30 % more segments for the same error and getting none of it back.
        #
        # The mean costs nothing to remove: a spatially uniform off-resonance is a phase that depends
        # on time alone, so it multiplies the *data* rather than the image and is applied exactly,
        # with no interpolation anywhere.  That is what `extra_phase_rad` already does, so it folds
        # straight in.
        self.offset_hz = float(np.mean(field))
        self.b0 = (field - self.offset_hz).astype(np.float32)
        self.segments = max(2, int(segments))
        self.weights_sqrt = (
            None if weights_sqrt is None
            else np.asarray(weights_sqrt, dtype=np.float32).reshape(-1)
        )
        uniform = -2.0 * np.pi * self.offset_hz * self.times
        self.extra_phase = (
            uniform.astype(np.float32) if extra_phase_rad is None
            else (np.asarray(extra_phase_rad, dtype=np.float64).reshape(-1)
                  + uniform).astype(np.float32)
        )
        self.device = int(device)
        self.n_samples = self.coord.shape[0]

        if np.any(np.diff(self.times) < -1e-12):
            msg = format_error(
                'sample_times_s must be non-decreasing for time segmentation.',
                {'samples': self.n_samples},
                ['sort the samples by time, or use the exact operator (segments=None)'],
            )
            raise ValueError(msg)

        # Segment centres, and for each the contiguous slice of samples it touches.
        first, last = float(self.times[0]), float(self.times[-1])
        span = max(last - first, 1e-12)
        self.centres = first + span * np.arange(self.segments) / (self.segments - 1)
        step = span / (self.segments - 1)
        self.ranges: list[tuple[int, int, np.ndarray]] = []
        for centre in self.centres:
            lo = int(np.searchsorted(self.times, centre - step, side='left'))
            hi = int(np.searchsorted(self.times, centre + step, side='right'))
            if hi <= lo:
                self.ranges.append((0, 0, np.zeros(0, dtype=np.float32)))
                continue
            weight = 1.0 - np.abs(self.times[lo:hi] - centre) / step
            self.ranges.append((lo, hi, np.clip(weight, 0.0, 1.0).astype(np.float32)))

    def linop(self) -> Any:
        """Return the sigpy ``Linop`` implementing this operator."""
        sp = _sigpy()
        outer = self
        device = sp.Device(outer.device)
        xp = device.xp

        state = {
            'sens': sp.to_device(outer.sens, device),
            'weights': None if outer.weights_sqrt is None
            else sp.to_device(outer.weights_sqrt, device),
            'forward_phase': None if outer.extra_phase is None
            else sp.to_device(np.exp(1j * outer.extra_phase).astype(np.complex64), device),
            'adjoint_phase': None if outer.extra_phase is None
            else sp.to_device(np.exp(-1j * outer.extra_phase).astype(np.complex64), device),
        }
        # exp(-2i pi df t_l), one image per segment: the NUFFT supplies exp(-2i pi k.r), so this
        # carries off-resonance with the same sign, matching the exact operator and the scanner.
        # Precomputed because each is reused every iteration and there are only `segments` of them.
        state['modulation'] = [
            sp.to_device(
                np.exp(-2j * np.pi * outer.b0 * float(centre)).astype(np.complex64), device
            )
            for centre in outer.centres
        ]
        state['coord'] = [
            sp.to_device(outer.coord[lo:hi], device) if hi > lo else None
            for lo, hi, _ in outer.ranges
        ]
        state['hat'] = [
            sp.to_device(weight, device) if weight.size else None
            for _, _, weight in outer.ranges
        ]

        # sigpy's NUFFT is normalised -- it divides by sqrt(n_pixels) and by the interpolation
        # width to the power of the dimension -- while :class:`NoncartesianOffresonance` computes the
        # plain sum ``sum_r x(r) exp(-2i pi k.r)``.  The two therefore differ by a constant, and
        # only in *magnitude*: comparing the phase of a single sample, which is the obvious check,
        # shows perfect agreement and hides it completely.  Left in, the reconstructed image comes
        # out scaled by several hundred, which also shifts what the Tikhonov weight means.
        #
        # Measured rather than copied from sigpy's internals, so it cannot drift with a release: a
        # delta at the grid centre sits at r = 0, so every sample of its transform is exactly 1.
        reference = next((c for c in state['coord'] if c is not None), None)
        scale = 1.0
        if reference is not None:
            probe = np.zeros((outer.ny, outer.nx), dtype=np.complex64)
            probe[outer.ny // 2, outer.nx // 2] = 1.0
            measured = float(
                xp.mean(xp.abs(sp.nufft(sp.to_device(probe, device), reference)))
            )
            scale = measured if measured > 1e-30 else 1.0
        state['scale'] = scale

        class Forward(sp.linop.Linop):
            """Image to k-space."""

            def __init__(self) -> None:
                super().__init__((outer.n_coils, outer.n_samples), (outer.ny, outer.nx))

            def _apply(self, x: Any) -> Any:
                image = sp.to_device(x, device).astype(xp.complex64, copy=False)
                y = xp.zeros((outer.n_coils, outer.n_samples), dtype=xp.complex64)
                weighted = state['sens'] * image[None, ...]
                for index, (lo, hi, _) in enumerate(outer.ranges):
                    if hi <= lo:
                        continue
                    segment = weighted * state['modulation'][index][None, ...]
                    y[:, lo:hi] += (
                        sp.nufft(segment, state['coord'][index]) / state['scale']
                        * state['hat'][index][None, :]
                    )
                if state['weights'] is not None:
                    y = y * state['weights'][None, :]
                if state['forward_phase'] is not None:
                    y = y * state['forward_phase'][None, :]
                return y

            def _adjoint_linop(self) -> Any:
                return Adjoint()

        class Adjoint(sp.linop.Linop):
            """k-space to image."""

            def __init__(self) -> None:
                super().__init__((outer.ny, outer.nx), (outer.n_coils, outer.n_samples))

            def _apply(self, y: Any) -> Any:
                data = sp.to_device(y, device).astype(xp.complex64, copy=False)
                if state['weights'] is not None:
                    data = data * state['weights'][None, :]
                if state['adjoint_phase'] is not None:
                    data = data * state['adjoint_phase'][None, :]
                image = xp.zeros((outer.ny, outer.nx), dtype=xp.complex64)
                conjugate = xp.conj(state['sens'])
                for index, (lo, hi, _) in enumerate(outer.ranges):
                    if hi <= lo:
                        continue
                    piece = data[:, lo:hi] * state['hat'][index][None, :]
                    # Divided by the same real factor as the forward, so the pair stays adjoint.
                    coils = sp.nufft_adjoint(
                        piece, state['coord'][index],
                        oshape=(outer.n_coils, outer.ny, outer.nx),
                    ) / state['scale']
                    image += xp.sum(
                        conjugate * coils * xp.conj(state['modulation'][index])[None, ...], axis=0
                    )
                return image

            def _adjoint_linop(self) -> Any:
                return Forward()

        return Forward()


def reconstruct(
    data: np.ndarray,
    k_per_m: np.ndarray,
    sample_times_s: np.ndarray,
    *,
    fov_m: float,
    matrix: int,
    b0_hz: np.ndarray | float = 0.0,
    sens: np.ndarray | None = None,
    extra_phase_rad: np.ndarray | None = None,
    n_per_arm: int | None = None,
    dcf: str = 'radial',
    iterations: int = 12,
    lamda: float = 1e-3,
    segments: int | None = None,
    device: int = -1,
    progress: bool = False,
) -> np.ndarray:
    """
    Reconstruct one image by conjugate gradient on the exact off-resonance model.

    Parameters
    ----------
    data
        ``(n_coils, n_samples)`` or ``(n_samples,)`` complex k-space data.
    k_per_m
        ``(2, n_samples)`` or ``(3, n_samples)`` k-space coordinates in 1/m.
    sample_times_s
        ``(n_samples,)`` sample times **relative to the echo**.
    fov_m
        Field of view, metres.  Converts the coordinates to cycles per FOV.
    matrix
        Reconstructed matrix, assumed square.
    b0_hz
        Off-resonance map in hertz, ``(matrix, matrix)``, or a scalar.  Pass ``0.0`` for the
        uncorrected reconstruction.
    sens
        ``(n_coils, matrix, matrix)`` coil sensitivities.  Defaults to ones.
    extra_phase_rad
        Optional spatially uniform, time-varying phase -- a concomitant-field term.
    n_per_arm
        Samples per interleaf, for the density-compensation floor.
    dcf
        Which density compensation to weight the samples by.  ``'radial'`` is :func:`density
        compensation <density_compensation>`, the ``|k|`` argument a spiral wants.  ``'speed'`` is
        :func:`speed_compensation`, which is the correct Jacobian for any trajectory sampled
        uniformly in time and the **only** correct one for a ramp-sampled EPI, whose samples are
        sparsest at k = 0.  ``'none'`` weights every sample equally.
    iterations
        Conjugate-gradient iterations.
    lamda
        Tikhonov weight.  Non-zero because the operator is ill-conditioned wherever k-space is
        undersampled, which for a variable-density spiral is the outer ring by design.
    segments
        ``None`` uses the **exact** dense operator: every sample against every pixel, no
        approximation anywhere, and ``n_samples * matrix^2`` complex exponentials per application.
        That is the right thing for verifying a sequence at a small matrix and unusable at a real
        one -- 59 100 samples at 128 x 128 is 2.8 hours per solve.

        A positive integer switches to :class:`SegmentedOffresonance`, which needs a few NUFFTs
        instead -- seconds in place of hours on that same problem.  Pass ``0``, the default, to let
        :func:`segments_for` derive the count from the field's range and the readout length.
    device
        sigpy device; ``-1`` for CPU.
    progress
        Show sigpy's progress bar.

    Returns
    -------
    numpy.ndarray
        ``(matrix, matrix)`` complex image.

    Notes
    -----
    ``sample_times_s`` must be **relative to the echo**, not absolute.  Absolute times put a large
    common phase on every sample, which the solver then has to undo -- and for a spin echo the
    absolute time includes TE, so the off-resonance term would be wrong by ``2*pi*df*TE`` uniformly.

    Iteration count is not a detail.  The operator is ill-conditioned exactly when the trajectory is
    undersampled, so conjugate gradient converges slowly there: on a 4x undersampled shot with a
    12-channel array the error falls from 0.156 at 30 iterations to 0.052 at 200, and stopping early
    leaves a residual that looks like aliasing rather than like non-convergence.  Twelve is enough
    only for a fully sampled, well-conditioned problem.
    """
    sp = _sigpy()

    data = np.atleast_2d(np.asarray(data, dtype=np.complex64))
    if data.shape[0] > data.shape[1]:  # (n_samples, n_coils) given
        data = data.T
    n_coils, n_samples = data.shape

    k = np.asarray(k_per_m, dtype=np.float64)
    coord = np.stack([k[1] * fov_m, k[0] * fov_m], axis=-1).astype(np.float32)  # (ky, kx)
    if dcf == 'radial':
        weights = density_compensation(k, n_per_arm=n_per_arm)
    elif dcf == 'speed':
        weights = speed_compensation(k, sample_times_s)
    elif dcf == 'none':
        weights = np.ones(n_samples, dtype=np.float32)
    else:
        msg = format_error(
            f'unknown density compensation {dcf!r}.',
            {'given': dcf},
            ["use 'radial' for a spiral, 'speed' for a ramp-sampled EPI, or 'none'"],
        )
        raise ValueError(msg)
    field = (
        np.full((matrix, matrix), float(b0_hz), dtype=np.float32)
        if np.isscalar(b0_hz) else np.asarray(b0_hz, dtype=np.float32)
    )
    if sens is None:
        sens = np.ones((1, matrix, matrix), dtype=np.complex64)

    shared = {
        'matrix': (matrix, matrix),
        'sens': sens,
        'coord_cycles_per_fov': coord,
        'sample_times_s': sample_times_s,
        'b0_hz': field,
        'weights_sqrt': np.sqrt(weights),
        'extra_phase_rad': extra_phase_rad,
        'device': device,
    }
    if segments is None:
        operator = NoncartesianOffresonance(**shared).linop()
    else:
        count = int(segments)
        if count <= 0:
            count = segments_for(field, sample_times_s)
        operator = SegmentedOffresonance(**shared, segments=count).linop()

    weighted = (data * np.sqrt(weights)[None, :]).astype(np.complex64)
    solver = sp.app.LinearLeastSquares(
        operator,
        sp.to_device(weighted, sp.Device(device)),
        lamda=lamda,
        solver='ConjugateGradient',
        max_iter=int(iterations),
        show_pbar=progress,
    )
    return np.asarray(sp.to_device(solver.run(), sp.Device(-1)), dtype=np.complex64)
