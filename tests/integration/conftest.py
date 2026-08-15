"""
Sequence builders for the integration tests.

These live **in the tests**, not in the package.  seqcraft ships a tree, a compiler and a contract;
a sequence is something you assemble, and baking one into library code would mean changing a
package to change your own scan.

They are also built out of **raw pypulseq events and nothing else**.  That is the point of this
file after the module reform: it is the standing proof that the compile path stands alone, needing
no module layer at all.  Two sequences, a gradient echo and a spin echo, chosen because between
them they exercise everything the compiler has to get right -- three winders overlapping on three
axes, an RF that forces a boundary, an ADC that forbids one, labels that must reach the right
readout, and a refocusing pulse that inverts the sign of everything before it.

Each is deliberately the shortest honest version: enough to check against physics, with none of a
notebook's configurability.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import pypulseq as pp
import pytest
from pypulseq.opts import Opts

import seqcraft as sc

if TYPE_CHECKING:
    from seqcraft.result import CompiledSequence

#: Small enough that the whole tier runs in a few seconds.
MATRIX = 32
FOV_MM = 250.0
SLICE_MM = 5.0

#: The scanner every recipe here is built against.  Stated in full, including the dead times
#: pypulseq would otherwise default to zero.
OPTS = Opts(
    max_grad=40, grad_unit='mT/m',
    max_slew=150, slew_unit='T/m/s',
    B0=3.0,
    rf_dead_time=100e-6,
    rf_ringdown_time=30e-6,
    adc_dead_time=10e-6,
)


def interleaved_slice_order(n: int) -> list[int]:
    """
    Even slices then odd ones, so no two consecutive excitations are neighbours.

    Adjacent slices share their profile's tails, and exciting one right after the other saturates
    that overlap -- which shows up as alternating bright and dark slices in the image.
    """
    return [*range(0, n, 2), *range(1, n, 2)]


def rf_spoil_phase(index: int, increment_rad: float) -> float:
    """
    Quadratic RF-spoiling phase ``phi_n = increment * n(n+1)/2``, wrapped into ``[0, 2 pi)``.

    Evaluated in closed form rather than accumulated: an accumulator makes shot 900 depend on
    shots 1 to 899 having been built first, which is a reproducibility hazard rather than a physics
    one -- but a real one.
    """
    return (increment_rad * index * (index + 1) / 2.0) % (2.0 * math.pi)


def slice_positions_m(n_slices: int, thickness_mm: float, gap_mm: float = 0.0) -> list[float]:
    """
    Slice centre offsets along the slice axis, in metres, symmetric about zero.

    Scaled to metres last: ``(i * pitch_mm) / 1e3`` is exact for the common integer-millimetre
    pitches, whereas ``i * (pitch_mm / 1e3)`` accumulates float noise that then shows up in the
    ``SlicePositions`` definition.
    """
    pitch_mm = thickness_mm + gap_mm
    first = -(n_slices - 1) / 2.0
    return [(first + i) * pitch_mm / 1e3 for i in range(n_slices)]


def geometry_definitions(
    *, fov_mm: float, matrix: int, slice_mm: float, n_slices: int,
) -> dict[str, object]:
    """
    The ``[DEFINITIONS]`` a fully sampled 2D acquisition is responsible for.

    Written out here rather than derived from a ``Geometry`` object, because `compile` takes
    pulseq's own definition keys and nothing else: FOV, matrix and slice order are decisions about
    the scan, and the compiler is indifferent to why the tree looks the way it does.  A geometry
    dataclass that produces exactly this mapping is in ``salvage/geometry.py`` for when a module
    library wants one.

    The one rule worth stating: ``kSpaceCenterLine`` here and the ``LIN`` label values below come
    from the same expression, ``matrix // 2``.  The reference implementation wrote
    ``kSpaceCenterLine = Ny/2 = 73.0`` while its navigator computed 36.5, and the two disagreeing
    is precisely the failure that costs an image.
    """
    return {
        'FOV': [fov_mm / 1e3, fov_mm / 1e3, slice_mm * n_slices / 1e3],
        'SliceThickness': slice_mm / 1e3,
        'SliceGap': 0.0,
        'SlicePositions': slice_positions_m(n_slices, slice_mm),
        'BaseResolution': matrix,
        'PhaseResolution': 1.0,
        'kSpaceCenterLine': matrix // 2,
        'ReadoutOversamplingFactor': 1.0,
    }


def _readout(opts: Opts, *, fov_mm: float, matrix: int, duration_s: float = 3.2e-3):
    """
    A frequency-encoding lobe, its ADC, and where k = 0 falls inside it.

    ``time_to_echo`` is the module-layer question the block cannot answer -- an offset into the
    readout, not a duration -- so it is computed here beside the gradient it describes.  Sampling
    starts when the ramp ends, so the echo is half a flat top later.
    """
    dk = 1e3 / fov_mm                                       # 1/FOV, in 1/m
    gx = pp.make_trapezoid('x', flat_area=matrix * dk, flat_time=duration_s, system=opts)
    adc = pp.make_adc(num_samples=matrix, duration=duration_s, delay=gx.rise_time, system=opts)
    time_to_echo = float(gx.rise_time) + 0.5 * float(gx.flat_time)
    return gx, adc, time_to_echo, dk


def build_gre(
    opts: Opts = OPTS,
    *,
    fov_mm: float = FOV_MM,
    matrix: int = MATRIX,
    slice_mm: float = SLICE_MM,
    te_s: float = 8e-3,
    tr_s: float = 20e-3,
    flip_deg: float = 15.0,
    n_slices: int = 2,
) -> CompiledSequence:
    """
    A spoiled 2D gradient echo.

    The three winders -- slice rephaser on z, phase blip on y, readout prephaser on x -- are placed
    at one time on three axes, which is the overlap the whole design exists to make free.
    """
    defs = geometry_definitions(
        fov_mm=fov_mm, matrix=matrix, slice_mm=slice_mm, n_slices=n_slices)
    positions = slice_positions_m(n_slices, slice_mm)
    centre_line = matrix // 2
    raster = sc.Raster(opts.block_duration_raster, 'block')

    rf, gz, gz_reph = pp.make_sinc_pulse(
        flip_angle=math.radians(flip_deg), duration=1e-3, slice_thickness=slice_mm * 1e-3,
        apodization=0.5, time_bw_product=4, delay=opts.rf_dead_time, use='excitation',
        system=opts, return_gz=True,
    )
    gx, adc, time_to_echo, dk = _readout(opts, fov_mm=fov_mm, matrix=matrix)
    gx_pre = pp.make_trapezoid('x', area=-float(gx.area) / 2.0, duration=1e-3, system=opts)
    # Four twists of phase across a voxel is enough to spoil the transverse magnetisation.
    spoil = pp.make_trapezoid('z', area=4.0 / (slice_mm * 1e-3), system=opts)

    # The RF's effective centre: a symmetric sinc refocuses half a pulse after its own start.
    isodelay = float(rf.delay) + 0.5 * float(rf.shape_dur)
    t_winders = float(pp.calc_duration(gz))
    t_readout = raster.ceil(isodelay + te_s - time_to_echo)
    t_spoil = raster.ceil(t_readout + float(pp.calc_duration(gx)))

    seq = sc.LogicBlock('gre_2d')
    index = 0
    for slice_index in interleaved_slice_order(n_slices):
        z = positions[slice_index]
        # Fully sampled, so the phase-encode table is every line on the recon grid.  Which lines a
        # partial-Fourier or accelerated acquisition takes is a sequence-programming choice;
        # salvage/geometry_pe.py has that arithmetic, including the residue nudge that keeps k = 0
        # on the sampled lattice.
        for line in range(matrix):
            t0 = index * tr_s
            phase = rf_spoil_phase(index, math.radians(117.0))
            # A slice is selected by offsetting the RF's frequency, not by moving the gradient.
            shifted = sc.events.derive(
                rf,
                freq_offset=float(gz.amplitude) * z,
                phase_offset=phase - 2.0 * math.pi * float(gz.amplitude) * z * float(rf.delay),
            )
            pe = pp.make_trapezoid(
                'y', area=(line - centre_line) * dk, duration=1e-3, system=opts)
            seq.add(t0, shifted, gz)
            seq.add(t0 + t_winders, gz_reph, pe, gx_pre)
            seq.add(t0 + t_readout, gx, adc)
            seq.add(t0 + t_spoil, spoil)
            seq.add(t0 + t_readout,
                    pp.make_label('LIN', 'SET', line),
                    pp.make_label('SLC', 'SET', slice_index))
            index += 1

    return sc.compile(
        seq, opts, name='gre_2d',
        definitions={**defs, 'TE': te_s, 'TR': tr_s, 'FlipAngle': flip_deg,
                     'RfSpoilIncrementDeg': 117.0},
    )


def build_se(
    opts: Opts = OPTS,
    *,
    fov_mm: float = FOV_MM,
    matrix: int = MATRIX,
    slice_mm: float = SLICE_MM,
    te_s: float = 20e-3,
    tr_s: float = 500e-3,
    n_slices: int = 2,
) -> CompiledSequence:
    """
    A 2D spin echo.

    Both winders sit before the refocusing pulse and therefore carry the **opposite** sign: a 180
    inverts accumulated phase, so a prephaser ahead of it must prephase the other way.  With the
    conventional sign the readout ends at 3*k_max and the phase encode is mirrored -- the second
    invisible in a k-space extent check, because ``|k|`` is symmetric.
    """
    defs = geometry_definitions(
        fov_mm=fov_mm, matrix=matrix, slice_mm=slice_mm, n_slices=n_slices)
    positions = slice_positions_m(n_slices, slice_mm)
    centre_line = matrix // 2
    raster = sc.Raster(opts.block_duration_raster, 'block')

    rf, gz, _ = pp.make_sinc_pulse(
        flip_angle=math.pi / 2, duration=2e-3, slice_thickness=slice_mm * 1e-3,
        apodization=0.5, time_bw_product=4, delay=opts.rf_dead_time, use='excitation',
        system=opts, return_gz=True,
    )
    rf180, gz180, _ = pp.make_sinc_pulse(
        flip_angle=math.pi, duration=4e-3, slice_thickness=slice_mm * 1e-3,
        apodization=0.5, time_bw_product=4, delay=opts.rf_dead_time, use='refocusing',
        system=opts, return_gz=True,
    )
    # Crushers straddle the 180 so anything it fails to invert is dephased rather than imaged.
    crusher = pp.make_trapezoid('z', area=4.0 / (slice_mm * 1e-3), system=opts)
    gx, adc, time_to_echo, dk = _readout(opts, fov_mm=fov_mm, matrix=matrix)
    gx_pre = pp.make_trapezoid('x', area=float(gx.area) / 2.0, duration=1e-3, system=opts)

    isodelay = float(rf.delay) + 0.5 * float(rf.shape_dur)
    isodelay180 = float(rf180.delay) + 0.5 * float(rf180.shape_dur)
    t_winders = float(pp.calc_duration(gz))
    refoc_duration = float(pp.calc_duration(gz180)) + 2.0 * float(pp.calc_duration(crusher))

    # TE is bounded from below at both halves; take whichever binds.
    first = raster.ceil(t_winders - isodelay + 1e-3 + float(pp.calc_duration(crusher))
                        + isodelay180)
    second = raster.ceil(refoc_duration - isodelay180 + time_to_echo)
    te = max(te_s, 2 * max(first, second))

    t_refoc = raster.ceil(isodelay + te / 2 - isodelay180 - float(pp.calc_duration(crusher)))
    t_readout = raster.ceil(isodelay + te - time_to_echo)

    seq = sc.LogicBlock('se_2d')
    index = 0
    for slice_index in interleaved_slice_order(n_slices):
        z = positions[slice_index]
        for line in range(matrix):                      # fully sampled; see build_gre
            t0 = index * tr_s
            offsets = {'freq_offset': float(gz.amplitude) * z}
            pe = pp.make_trapezoid(
                'y', area=-(line - centre_line) * dk, duration=1e-3, system=opts)
            seq.add(t0, sc.events.derive(rf, **offsets), gz)
            seq.add(t0 + t_winders, pe, gx_pre)
            seq.add(t0 + t_refoc, crusher)
            seq.add(t0 + t_refoc + float(pp.calc_duration(crusher)),
                    sc.events.derive(rf180, freq_offset=float(gz180.amplitude) * z), gz180)
            seq.add(t0 + t_refoc + float(pp.calc_duration(crusher))
                    + float(pp.calc_duration(gz180)), crusher)
            seq.add(t0 + t_readout, gx, adc)
            seq.add(t0 + t_readout,
                    pp.make_label('LIN', 'SET', line),
                    pp.make_label('SLC', 'SET', slice_index))
            index += 1

    return sc.compile(
        seq, opts, name='se_2d',
        definitions={**defs, 'TE': te, 'TR': tr_s},
    )


@pytest.fixture(scope='module')
def gre() -> CompiledSequence:
    return build_gre()


@pytest.fixture(scope='module')
def se() -> CompiledSequence:
    return build_se()


@pytest.fixture(params=['gre', 'se'])
def compiled(request) -> CompiledSequence:
    """Every recipe, for the checks that must hold of all of them."""
    return request.getfixturevalue(request.param)
