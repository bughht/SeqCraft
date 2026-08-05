"""
Bridge to MRzeroCore: simulate what the compiled sequence will actually play.

MRzero's own ``.seq`` reader handles **trapezoid gradients only**, so a spiral cannot be imported
through it at all.  :func:`to_mr0` therefore takes the events directly -- either from a compiled
:class:`~seqcraft.core.compiler.CompiledSequence`, or from a ``pypulseq.Sequence`` that has read the
written file back, since pypulseq's reader does handle arbitrary waveforms.

Reading the file back is the stronger of the two.  What gets simulated is then the bytes that will be
handed to the scanner rather than an object that ought to correspond to them, and "ought to" is where
a stale sidecar or a lossy round trip hides.

What is modelled
----------------
MRzero's phase-distribution-graph simulation covers T1, T2, T2', B0 off-resonance, B1 and
**diffusion**, which is what makes it usable for checking a DWI sequence rather than only its
geometry:

.. code-block:: text

    b        = (1/3) (2 pi)^2 dt (k1^2 + k1 k2 + k2^2)     accumulated per event
    signal  *= exp(-1e-9 * D * cumsum(b))                  D in 1e-3 mm^2/s
    phase   *= exp(2i pi (t * B0 + k . r))                 B0 in Hz

So a b-value error, a wrong diffusion direction, or a missing refocusing sign flip all show up as
the wrong attenuation -- and off-resonance shows up as the spiral blurring it really causes.

How the conversion works
------------------------
A repetition begins at every RF pulse, which is MRzero's model.  Within a repetition each pulseq
block is subdivided into events fine enough to resolve its gradients, because ``b`` is accumulated
per event and a single event spanning a 15 ms diffusion lobe would integrate it as one straight
line.  Every ADC sample becomes its own event, since that is where signal is recorded.

This lives in ``examples/lib`` rather than inside seqcraft on purpose.  seqcraft's job is to build
sequences; simulating and reconstructing them are downstream concerns with heavy dependencies of
their own (MRzeroCore, torch), and folding them into the package would make every user pay for them.

Usage
-----
::

    import sys; sys.path.insert(0, 'examples/lib')
    from mr0_bridge import to_mr0, simulate, phantom_with_diffusion

    seq, meta = to_mr0(compiled)                       # straight from the compiler
    played = pypulseq.Sequence(system=system.limits('dwi'))
    played.read('scan.seq')
    seq, meta = to_mr0(played, system=system)          # or from the written file
    signal = simulate(seq, phantom_with_diffusion(size_m=..., matrix=64, d_1e_3_mm2_s=0.8))
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

from seqcraft.core import events as ev
from seqcraft.core.errors import MissingExtraError, format_error

if TYPE_CHECKING:
    from seqcraft.core.compiler import CompiledSequence

__all__ = [
    'Mr0Meta',
    'coil_sensitivities',
    'phantom_with_diffusion',
    'simulate',
    'to_mr0',
]

_AXES = ('x', 'y', 'z')


def _mr0() -> Any:
    """Import MRzeroCore, with a message naming the extra to install."""
    try:
        import MRzeroCore  # noqa: PLC0415
    except ImportError as err:  # pragma: no cover - depends on the environment
        msg = format_error(
            'simulation needs MRzeroCore.',
            {'missing': 'MRzeroCore'},
            ['pip install "seqcraft[sim]"'],
        )
        raise MissingExtraError(msg) from err
    return MRzeroCore


def _torch() -> Any:
    """Import torch, with a message naming the extra to install."""
    try:
        import torch  # noqa: PLC0415
    except ImportError as err:  # pragma: no cover
        msg = format_error(
            'simulation needs torch.', {'missing': 'torch'}, ['pip install "seqcraft[sim]"']
        )
        raise MissingExtraError(msg) from err
    return torch


@dataclass
class Mr0Meta:
    """
    What the conversion learned, and what the reconstruction needs.

    Attributes
    ----------
    n_reps
        Number of MRzero repetitions, one per RF pulse.
    adc_times_s
        Absolute time of every ADC sample, from the start of the sequence.
    adc_rep
        Which repetition each ADC sample belongs to.
    kspace_per_m
        ``(3, n_samples)`` k-space coordinate of every ADC sample, in 1/m, taken from the gradient
        moments actually handed to the simulator -- so it is the trajectory that was simulated, not
        a separate calculation that could disagree.

        Measured **from each shot's first ADC sample**, not from the start of the sequence.  A
        running total would carry the diffusion lobes' and crushers' enormous excursions into the
        readout: for a b=1000 encoding the accumulated k is many times ``k_max``, and the
        refocusing pulse's sign flip is not in it either.  What the reconstruction needs is the k
        that encodes position *during the readout*, which is exactly the increment since sampling
        began.
    echo_times_s
        Absolute time of each repetition's exciting pulse.
    flip_angles_deg
        Flip angle of each repetition's pulse, integrated from the RF waveform.
    """

    n_reps: int
    adc_times_s: np.ndarray
    adc_rep: np.ndarray
    kspace_per_m: np.ndarray
    echo_times_s: np.ndarray
    flip_angles_deg: np.ndarray
    events_per_rep: list[int] = field(default_factory=list)

    @property
    def n_samples(self) -> int:
        """Total ADC samples."""
        return int(self.adc_times_s.size)

    def samples_of_rep(self, index: int) -> np.ndarray:
        """Return the sample indices belonging to repetition `index`."""
        return np.flatnonzero(self.adc_rep == index)


def _flip_and_phase(rf: Any, raster: float) -> tuple[float, float]:
    """
    Return the flip angle in radians and the carrier phase of an RF event.

    The flip is obtained by **Bloch-simulating the pulse on resonance**: decompose it into hard-pulse
    rotations, apply them to ``M = z``, and read ``flip = arccos(Mz)`` off the result.

    The obvious alternative, ``2*pi * integral B1 dt``, is the small-tip approximation and is wrong
    exactly where it matters.  An SLR ``'se'`` refocusing pulse is designed so its *refocusing
    profile* is flat, and its envelope has almost **zero net area** -- integrate it and a 180 degree
    pulse reports as 0 degrees, so the sequence simulates as producing no signal at all with nothing
    resembling a pulse problem to point at.  A sinc's negative side lobes cause a milder version of
    the same error in the other direction.

    Two details the integration has to respect:

    The integration runs over **intervals** of ``rf.t``, not over samples.  ``make_block_pulse``
    stores a rectangular pulse as *two* samples spanning one interval: assume raster spacing and a
    90 degree hard pulse reports as 0.36 degrees; step once per sample instead and it reports as
    180.

    The phase comes from ``rf.phase_offset`` rather than from the simulated axis.  For a good 180 the
    transverse component of the result vanishes, so the rotation axis is not recoverable from it --
    the sign of a numerically tiny residual decides the answer, and it came out a half turn wrong.
    pypulseq puts the carrier phase in that field by construction, so there is nothing to infer.
    """
    signal = np.asarray(rf.signal, dtype=complex)
    times = np.asarray(getattr(rf, 't', None), dtype=float) if hasattr(rf, 't') else None
    if times is None or times.size != signal.size:
        times = np.arange(signal.size, dtype=float) * raster

    # Rotate interval by interval, not sample by sample.  A rectangular pulse is two samples across
    # one interval; stepping per sample applies that interval twice and doubles the flip.
    if signal.size < 2:
        amplitudes = signal
        dwells = np.array([float(getattr(rf, 'shape_dur', raster))])
    else:
        amplitudes = 0.5 * (signal[:-1] + signal[1:])       # trapezoidal, exact for linear segments
        dwells = np.diff(times)

    magnetisation = np.array([0.0, 0.0, 1.0])
    for sample, dwell in zip(amplitudes, dwells):
        omega = 2.0 * math.pi * np.array([float(sample.real), float(sample.imag), 0.0])
        norm = float(np.linalg.norm(omega))
        if norm < 1e-12 or dwell <= 0.0:
            continue
        axis = omega / norm
        angle = norm * float(dwell)
        cos, sin = math.cos(angle), math.sin(angle)
        magnetisation = (
            magnetisation * cos
            + np.cross(axis, magnetisation) * sin
            + axis * float(np.dot(axis, magnetisation)) * (1.0 - cos)
        )

    flip = math.acos(max(-1.0, min(1.0, float(magnetisation[2]))))
    return flip, float(rf.phase_offset)


def _grad_samples(block: Any, duration: float, raster: float) -> np.ndarray:
    """Return ``(n, 3)`` gradient amplitudes in Hz/m on a uniform grid covering the block."""
    n = max(1, int(round(duration / raster)))
    grid = (np.arange(n) + 0.5) * raster
    out = np.zeros((n, 3))
    for axis_index, axis in enumerate(_AXES):
        grad = getattr(block, f'g{axis}', None)
        if grad is None:
            continue
        tt, wf = ev.waveform_of(grad, raster)
        out[:, axis_index] = np.interp(grid, tt, wf, left=0.0, right=0.0)
    return out


def _areas(block: Any, duration: float, raster: float, edges: np.ndarray) -> np.ndarray:
    """
    Return the exact gradient area in each interval of `edges`, as ``(len(edges) - 1, 3)``.

    Integrates the sampled waveform once, cumulatively, and differences the result at the interval
    edges.  Total area is then preserved no matter how the intervals are chosen, which is the whole
    point: the obvious alternative -- resample the gradient *at* the edges and apply the trapezoidal
    rule -- silently loses area whenever an interval is long compared with the gradient's features.

    In the worst case it loses all of it.  A trapezoid gradient starts and ends at zero, so one
    interval spanning a whole prephaser samples zero at both ends and integrates to **exactly zero**:
    the phase encode and readout prephaser vanish, k-space starts at the origin instead of its corner,
    and a Cartesian image comes out as a bright smear with no hint of where it went.  Blocks shorter
    than about twice the event length are the ones at risk, which is most winders.
    """
    fine = _grad_samples(block, duration, raster)
    # Rectangle rule per raster, which is what the hardware plays: each sample is held for one raster.
    cumulative = np.concatenate([np.zeros((1, 3)), np.cumsum(fine, axis=0) * raster], axis=0)
    sample_edges = np.arange(len(fine) + 1) * raster
    at_edges = np.stack([
        np.interp(edges, sample_edges, cumulative[:, axis]) for axis in range(3)
    ], axis=1)
    return np.diff(at_edges, axis=0)


def _drop(moments: np.ndarray, axes: tuple[int, ...]) -> np.ndarray:
    """Zero the given axis columns of a gradient-moment array."""
    if not axes:
        return moments
    out = np.asarray(moments, dtype=float).copy()
    out[..., list(axes)] = 0.0
    return out


def to_mr0(
    compiled: CompiledSequence | Any,
    *,
    system: Any = None,
    event_dt_s: float = 200e-6,
    first_rep: int = 0,
    max_reps: int | None = None,
    ignore_axes: tuple[str, ...] = ('z',),
) -> tuple[Any, Mr0Meta]:
    """
    Convert a compiled sequence, or a ``.seq`` read back from disk, into an MRzero ``Sequence``.

    Parameters
    ----------
    compiled
        Either the result of :func:`~seqcraft.core.compiler.compile_sequence`, or a bare
        ``pypulseq.Sequence`` -- in which case pass `system` as well, since a ``.seq`` file records
        the rasters but not the :class:`~seqcraft.core.system.System` object they came from.

        Reading the written file back is the stronger check of the two: it simulates the bytes that
        will be handed to the scanner rather than an object that ought to correspond to them.  It is
        also the only way to simulate a spiral from a file at all, because MRzero's own ``.seq``
        reader handles trapezoid gradients only, while pypulseq's reads arbitrary waveforms.
    system
        Required when `compiled` is a plain ``pypulseq.Sequence``; ignored otherwise.
    event_dt_s
        Target event length for blocks with no ADC.  Diffusion lobes are integrated per event, so
        this has to be short enough to resolve them: a single event spanning a 15 ms lobe would
        integrate ``b`` as one straight line and under-report it badly.  200 us gives about 75
        events per lobe, which reproduces the analytic b-value to well under a percent.
    first_rep, max_reps
        Simulate repetitions ``[first_rep, first_rep + max_reps)`` and drop the rest.  Cost *and
        memory* scale with ADC samples times phantom voxels -- MRzero materialises that product --
        so a full 19-volume acquisition is far beyond what is needed to check a sequence, and even
        two shots of a 59 100-sample readout against a 160 x 160 phantom asks for 6.7 GB.

        When the shots are independent, which a 6 s TR and a spoiler make them, simulating one at a
        time is exact and halves the peak allocation: ``first_rep=0, max_reps=3`` for the first shot
        here and ``first_rep=3, max_reps=3`` for the second, since a sequence beginning with fat
        saturation spends three repetitions per shot.
    ignore_axes
        Axes whose gradients are dropped, defaulting to the slice axis.

        **A 2D phantom is one voxel thick, so it cannot represent through-slice structure**, and its
        Nyquist limit along that axis is correspondingly tiny -- 125 1/m for a 4 mm slice.  A
        slice-select gradient reaches 1000 1/m and a crusher more, so MRzero's voxel dephasing
        function reads the state as *completely* dephased and the simulated signal is exactly zero.
        That is the right answer to the question asked; it is not the question meant.  A 2D
        simulation assumes ideal slice selection, which means leaving the slice gradients out.

        Pass ``()`` to keep every axis -- correct only if the phantom actually resolves the slice
        direction, which needs several voxels through it.

    Returns
    -------
    sequence, meta
        The MRzero sequence and an :class:`Mr0Meta` describing what went into it.

    Notes
    -----
    Built with ``normalized_grads=False``, so gradient moments are in SI units (1/m) and MRzero
    applies no phantom-size scaling.  That keeps the k-space handed to the simulator identical to
    the one the reconstruction will use.
    """
    mr0 = _mr0()
    torch = _torch()

    seq = getattr(compiled, 'seq', compiled)
    resolved = system if system is not None else getattr(compiled, 'system', None)
    if resolved is None:
        msg = format_error(
            'a bare pypulseq.Sequence needs the System it was built for.',
            {'given': type(compiled).__name__},
            ['to_mr0(seq, system=my_system)'],
        )
        raise ValueError(msg)
    raster = resolved.grad_raster_s
    rf_raster = resolved.rf_raster_s
    dropped = tuple(_AXES.index(a) for a in ignore_axes if a in _AXES)

    # ---- pass one: group blocks into repetitions, one per RF pulse -------------------------
    reps: list[dict[str, Any]] = []
    t = 0.0
    seen = 0                    # RF pulses encountered, including those skipped
    for index in sorted(seq.block_events):
        block = seq.get_block(index)
        duration = float(seq.block_durations[index])
        if getattr(block, 'rf', None) is not None:
            seen += 1
            if seen > first_rep:
                if max_reps is not None and len(reps) >= max_reps:
                    break
                flip, phase = _flip_and_phase(block.rf, rf_raster)
                reps.append({
                    'flip': flip,
                    'phase': phase,
                    'use': str(getattr(block.rf, 'use', 'undefined')),
                    't_pulse': t + float(block.rf.delay) + float(block.rf.center),
                    'blocks': [],
                })
        if reps:
            reps[-1]['blocks'].append((t, index, block, duration))
        t += duration

    if not reps:
        msg = format_error(
            'the sequence contains no RF pulse, so there is nothing to simulate.',
            {'blocks': len(seq.block_events)},
        )
        raise ValueError(msg)

    # ---- pass two: build each repetition's events ------------------------------------------
    usages = {
        'excitation': mr0.PulseUsage.EXCIT,
        'refocusing': mr0.PulseUsage.REFOC,
        'saturation': mr0.PulseUsage.FATSAT,
        'inversion': mr0.PulseUsage.REFOC,
    }
    out = mr0.Sequence(normalized_grads=False)
    adc_times: list[float] = []
    adc_rep: list[int] = []
    kspace: list[np.ndarray] = []
    running_k = np.zeros(3)
    events_per_rep: list[int] = []

    for rep_index, entry in enumerate(reps):
        times: list[float] = []
        moments: list[np.ndarray] = []
        usage: list[int] = []

        for block_start, _index, block, duration in entry['blocks']:
            adc = getattr(block, 'adc', None)
            if adc is not None:
                # One event per ADC sample: that is where signal is recorded, and the gradient
                # moment accumulated over each dwell is what places it in k-space.
                n = int(adc.num_samples)
                dwell = float(adc.dwell)
                delay = float(adc.delay)
                tail = duration - delay - n * dwell

                # One event per ADC sample, plus the lead before sampling starts and any tail after,
                # so every part of the block's area is carried.
                edges = np.concatenate((
                    [0.0] if delay > 0 else [],
                    delay + np.arange(n + 1) * dwell,
                    [duration] if tail > 1e-12 else [],
                ))
                areas = _areas(block, duration, raster, edges)
                offset = 1 if delay > 0 else 0

                if delay > 0:
                    times.append(delay)
                    moments.append(_drop(areas[0], dropped))
                    usage.append(0)
                for sample in range(n):
                    increment = areas[offset + sample]
                    times.append(dwell)
                    moments.append(_drop(increment, dropped))
                    usage.append(1)
                    adc_times.append(block_start + delay + (sample + 0.5) * dwell)
                    adc_rep.append(rep_index)
                    running_k = running_k + increment
                    kspace.append(running_k.copy())
                if tail > 1e-12:
                    times.append(tail)
                    moments.append(_drop(areas[offset + n], dropped))
                    usage.append(0)
            else:
                # No ADC: subdivide finely enough for the diffusion integral to be right, since b
                # accumulates per event and one event spanning a 15 ms lobe would integrate it as a
                # straight line.  The area itself is exact at any subdivision -- see `_areas`.
                n = max(1, int(round(duration / event_dt_s)))
                dt = duration / n
                areas = _areas(block, duration, raster, np.arange(n + 1) * dt)
                for step in range(n):
                    times.append(dt)
                    moments.append(_drop(areas[step], dropped))
                    usage.append(0)
                    running_k = running_k + areas[step]

        if not times:
            continue
        rep = out.new_rep(len(times))
        rep.pulse.usage = usages.get(entry['use'], mr0.PulseUsage.UNDEF)
        rep.pulse.angle = torch.tensor(float(entry['flip']))
        rep.pulse.phase = torch.tensor(float(entry['phase']))
        rep.event_time[:] = torch.tensor(np.asarray(times, dtype=np.float32))
        rep.gradm[:] = torch.tensor(np.asarray(moments, dtype=np.float32))
        rep.adc_usage[:] = torch.tensor(np.asarray(usage, dtype=np.int64))
        rep.adc_phase[:] = math.pi / 2.0 - float(entry['phase'])
        events_per_rep.append(len(times))

    k_absolute = np.stack(kspace, axis=1) if kspace else np.zeros((3, 0))
    adc_rep_array = np.asarray(adc_rep, dtype=int)
    # Re-base each shot's trajectory on its own first sample, which for a spiral is k = 0.
    k_relative = k_absolute.copy()
    for rep in np.unique(adc_rep_array):
        inside = adc_rep_array == rep
        k_relative[:, inside] -= k_absolute[:, inside][:, :1]

    meta = Mr0Meta(
        n_reps=len(out),
        adc_times_s=np.asarray(adc_times),
        adc_rep=adc_rep_array,
        kspace_per_m=k_relative,
        echo_times_s=np.asarray([e['t_pulse'] for e in reps]),
        flip_angles_deg=np.asarray([math.degrees(e['flip']) for e in reps]),
        events_per_rep=events_per_rep,
    )
    return out, meta


def coil_sensitivities(matrix: int, n_coils: int, *, ring_radius: float = 0.72) -> np.ndarray:
    """
    Return ``(n_coils, matrix, matrix)`` complex sensitivities for a ring of loop coils.

    Parameters
    ----------
    matrix
        Grid size, matching the phantom.
    n_coils
        Number of elements, evenly spaced on a circle.
    ring_radius
        Radius of the circle, as a fraction of the FOV.  Slightly outside the object, as a head
        array is.

    Returns
    -------
    numpy.ndarray
        Normalised so the root-sum-of-squares over channels is 1 at the centre, which keeps the
        simulated signal comparable with the single-coil case.

    Notes
    -----
    Without coil sensitivities, undersampling cannot be undone.  A single uniform coil measures one
    number per k-space sample and an undersampled trajectory then has genuinely fewer measurements
    than unknowns -- no reconstruction recovers what was never encoded, and the aliasing that
    results is a property of the experiment, not of the algorithm.  Simulating a 4x undersampled
    spiral with one coil therefore says nothing about how it will behave on a scanner with a
    32-channel array, where the spatial variation between channels supplies exactly the missing
    information.

    The model is the usual one for a surface loop: magnitude falling as ``1/(distance + offset)``,
    with a smoothly varying phase.  It is not a Biot--Savart calculation and is not meant to
    reproduce a specific array; what matters for testing a reconstruction is that the channels are
    spatially distinct and smooth, which is what makes them invertible.
    """
    axis = (np.arange(matrix) - matrix / 2) / matrix
    grid_a, grid_b = np.meshgrid(axis, axis, indexing='ij')
    out = np.empty((n_coils, matrix, matrix), dtype=np.complex64)
    for index in range(n_coils):
        angle = 2.0 * math.pi * index / n_coils
        distance = np.hypot(grid_a - ring_radius * math.cos(angle),
                            grid_b - ring_radius * math.sin(angle))
        magnitude = 1.0 / (distance + 0.35)
        # A phase that varies across the field of view, as a real element's does; without it the
        # channels differ only in magnitude and are far easier to unfold than the real thing.
        phase = 2.0 * math.pi * 0.35 * distance + angle
        out[index] = (magnitude * np.exp(1j * phase)).astype(np.complex64)
    centre = matrix // 2
    scale = float(np.sqrt(np.sum(np.abs(out[:, centre, centre]) ** 2)))
    return (out / max(scale, 1e-12)).astype(np.complex64)


def phantom_with_diffusion(
    *,
    size_m: tuple[float, float, float],
    matrix: int,
    d_1e_3_mm2_s: float = 0.8,
    b0_hz: np.ndarray | float = 0.0,
    t1_s: float = 1.0,
    t2_s: float = 0.08,
    t2dash_s: float = 0.05,
    shape: str = 'blobs',
    n_coils: int = 1,
) -> tuple[Any, dict[str, np.ndarray]]:
    """
    Build a 2D ``SimData`` phantom with a diffusion coefficient and an off-resonance map.

    Parameters
    ----------
    size_m
        Physical extent ``(x, y, z)`` in metres.
    matrix
        In-plane matrix of the phantom grid.  Make it at least twice the reconstructed matrix so
        the simulation is not itself band-limited to the thing being measured.
    d_1e_3_mm2_s
        Diffusion coefficient in **10^-3 mm^2/s**, which is MRzero's unit.  Free water at body
        temperature is about 3.0; white matter is 0.7 to 0.9 along a fibre and lower across it.
    b0_hz
        Off-resonance in hertz: a scalar, or a map matching the phantom grid.  This is what makes a
        spiral blur, so it is the interesting knob.
    t1_s, t2_s, t2dash_s
        Relaxation times.
    shape
        ``'blobs'`` for three asymmetric discs of differing diffusivity -- asymmetric on purpose,
        because a symmetric phantom hides orientation and mirroring errors.  ``'uniform'`` for a
        single disc, which is what you want when measuring attenuation rather than looking at an
        image.
    n_coils
        Number of receive channels.  ``1`` gives a uniform coil, which is right for measuring
        attenuation and **wrong for anything undersampled**: with one channel there is no spatial
        information to unfold aliasing with, so a 4x undersampled trajectory cannot be reconstructed
        no matter what the algorithm does.  Set this to the array you will actually scan with -- 20
        to 32 for a head coil -- and hand the returned ``'coil_sens'`` to the reconstruction.

    Returns
    -------
    sim_data, maps
        The MRzero ``SimData``, and a dict of the ``(matrix, matrix)`` grid maps it was built from:
        ``'PD'``, ``'D'``, ``'B0'``, plus ``'coil_sens'`` with shape ``(n_coils, matrix, matrix)``.

    The maps are returned because ``SimData`` keeps only the voxels above ``PD_threshold``, as a
    flat array -- so there is no way to reshape it back to a picture afterwards.  Anything that
    wants to *display* the phantom needs these.

    Notes
    -----
    Diffusion is given per voxel, so ``'blobs'`` produces regions with genuinely different
    attenuation -- which is the thing a DWI sequence is supposed to measure, and therefore the thing
    worth checking end to end.
    """
    mr0 = _mr0()
    torch = _torch()

    n = int(matrix)
    y, x = np.meshgrid(
        (np.arange(n) - n / 2) / n, (np.arange(n) - n / 2) / n, indexing='ij'
    )
    radius = np.hypot(x, y)

    pd = np.zeros((n, n), dtype=np.float32)
    diffusion = np.zeros((n, n), dtype=np.float32)
    if shape == 'uniform':
        inside = radius < 0.4
        pd[inside] = 1.0
        diffusion[inside] = d_1e_3_mm2_s
    elif shape == 'point':
        # One small disc, off centre along both axes and at different distances, so no reflection or
        # transpose can map it onto itself.  This is the phantom for checking orientation: a disc,
        # or anything symmetric, correlates about as well with its own transpose as with itself and
        # settles nothing.
        inside = np.hypot(x - 0.25, y + 0.15) < 0.07
        pd[inside] = 1.0
        diffusion[inside] = d_1e_3_mm2_s
    else:
        pd[radius < 0.42] = 1.0
        diffusion[radius < 0.42] = d_1e_3_mm2_s
        # Three discs, deliberately not symmetric about either axis.
        for cx, cy, r, scale, weight in (
            (-0.15, 0.12, 0.12, 0.25, 1.15),      # restricted, bright
            (0.18, 0.10, 0.09, 2.5, 0.85),        # nearly free water, dim
            (0.02, -0.20, 0.07, 1.0, 1.0),        # same D as background
        ):
            blob = np.hypot(x - cx, y - cy) < r
            pd[blob] = weight
            diffusion[blob] = d_1e_3_mm2_s * scale

    field = (
        np.full((n, n), float(b0_hz), dtype=np.float32)
        if np.isscalar(b0_hz) else np.asarray(b0_hz, dtype=np.float32)
    )
    if field.shape != (n, n):
        msg = format_error(
            'b0_hz must be a scalar or match the phantom grid.',
            {'given': str(field.shape), 'expected': f'({n}, {n})'},
        )
        raise ValueError(msg)

    def volume(values: np.ndarray) -> Any:
        """
        Hand a ``[y, x]`` map to MRzero, which indexes its grid ``[x, y]``.

        The transpose is the whole of it, and leaving it out is not a display quirk -- it silently
        reflects the object about its own diagonal.  MRzero's affine maps array index 0 to physical
        axis 0, so whatever varies along rows *is* x as far as the simulator is concerned, while
        every convention on the reconstruction side and in ``imshow`` treats rows as y.

        Owning it here rather than at the call site is what keeps it fixed.  Corrected further
        downstream -- by transposing the B0 and sensitivity maps on their way into the operator, say
        -- the maps line up and the *image* still comes out mirrored, which looks like a
        reconstruction bug and gets debugged as one.
        """
        return torch.tensor(np.ascontiguousarray(values.T)[:, :, None].astype(np.float32))

    mask = pd > 0
    shape = (n, n, 1)
    # Voxel-index to world affine, in millimetres, matching MRzero's own convention: a diagonal of
    # size/shape and a translation putting the grid centre at the isocentre.  Required since
    # MRzeroCore 1.0, and it is what makes voxel_pos come out in the right place -- get it wrong and
    # the image is scaled or shifted without anything complaining.
    affine = torch.eye(3, 4)
    for axis in range(3):
        affine[axis, axis] = size_m[axis] / shape[axis] * 1000.0
    affine[:, 3] = -torch.as_tensor(list(size_m)) / 2.0 * 1000.0

    sens = (
        np.ones((1, n, n), dtype=np.complex64) if n_coils <= 1
        else coil_sensitivities(n, int(n_coils))
    )
    phantom = mr0.VoxelGridPhantom(
        PD=volume(pd),
        T1=volume(np.where(mask, t1_s, 1e-6)),
        T2=volume(np.where(mask, t2_s, 1e-6)),
        T2dash=volume(np.where(mask, t2dash_s, 1e-6)),
        D=volume(diffusion),
        B0=volume(field * mask),
        # B1 and coil_sens both carry a leading transmit/receive-channel dimension.
        B1=volume(np.ones_like(pd))[None, ...],
        # Transposed per channel, for the same reason as `volume`.
        coil_sens=torch.tensor(
            np.ascontiguousarray(np.swapaxes(sens, -1, -2))[..., None].astype(np.complex64)
        ),
        size=torch.tensor(list(size_m)),
        affine=affine,
    )
    return phantom.build(), {
        'PD': pd,
        'D': diffusion,
        'B0': field * mask,
        'coil_sens': sens,
    }


def simulate(
    sequence: Any,
    phantom: Any,
    *,
    min_emitted_signal: float = 1e-2,
    min_latent_signal: float = 1e-2,
    progress: bool = False,
) -> np.ndarray:
    """
    Run the MRzero phase-distribution-graph simulation and return the complex signal.

    Parameters
    ----------
    sequence
        From :func:`to_mr0`.
    phantom
        From :func:`phantom_with_diffusion`, or any MRzero ``SimData``.
    min_emitted_signal, min_latent_signal
        Pruning thresholds for the phase-distribution graph.  Lower is more accurate and slower;
        the defaults are MRzero's own and are ample for checking a sequence.
    progress
        Print MRzero's progress bar.

    Returns
    -------
    numpy.ndarray
        ``(n_samples, n_coils)`` complex signal, in acquisition order.
    """
    mr0 = _mr0()
    graph = mr0.compute_graph(sequence, phantom, 200, 1e-3)
    signal = mr0.execute_graph(
        graph, sequence, phantom,
        min_emitted_signal=min_emitted_signal,
        min_latent_signal=min_latent_signal,
        print_progress=progress,
    )
    return signal.detach().cpu().numpy()
