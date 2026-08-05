"""
Sequence builders for the integration tests.

These live **in the tests**, not in the package. seqcraft ships modules and a compiler; a sequence is
something you assemble, and baking one into library code would mean changing a package to change
your own scan. The notebooks in ``examples/`` do the same thing at greater length and with the
reasoning spelled out.

Each builder here is deliberately the shortest honest version -- enough to exercise the compiler and
the physics checks, with none of the notebook's configurability.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import pypulseq as pp
import pytest

import seqcraft as sc

if TYPE_CHECKING:
    from seqcraft.core.compiler import CompiledSequence

#: Small enough that the whole tier runs in a few seconds.
MATRIX = 32
FOV_MM = 250.0
SLICE_MM = 5.0


def build_gre(
    system: sc.System,
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
    geometry = sc.Geometry(
        fov_mm=(fov_mm, fov_mm, slice_mm), matrix=(matrix, matrix, 1),
        slice_thickness_mm=slice_mm, n_slices=n_slices,
    )
    exc = sc.modules.SincExcitation(
        system, flip_deg=flip_deg, duration_us=1000, slice_thickness_mm=slice_mm, rephase=False)
    ro = sc.modules.CartesianLine(
        system, fov_ro_mm=fov_mm, matrix_ro=matrix, readout_duration_us=3200, prephase=False)
    pe = sc.modules.PhaseEncode(system, fov_pe_mm=fov_mm, matrix_pe=matrix)
    spoil = sc.modules.Spoiler(system, twists=4, voxel_mm=slice_mm)

    raster = system.block_raster_s
    t_winders = exc.duration
    t_readout = sc.raster.ceil_to(exc.isodelay + te_s - ro.time_to_echo, raster)
    t_spoil = sc.raster.ceil_to(t_readout + ro.duration, raster)

    seq = sc.LogicBlock('gre_2d')
    index = 0
    for slice_index in sc.interleaved_slice_order(n_slices):
        z = geometry.slice_positions_m[slice_index]
        for line in geometry.pe_lines:
            t0 = index * tr_s
            phase = sc.rf_spoil_phase(index, math.radians(117.0))
            seq.add(t0, exc.build(slice_offset_m=z, rf_phase_rad=phase))
            seq.add(t0 + t_winders, exc.rephaser())
            seq.add(t0 + t_winders, pe.build(line=line - geometry.kspace_center_line))
            seq.add(t0 + t_winders, ro.prephaser_block())
            seq.add(t0 + t_readout, ro.build())
            seq.add(t0 + t_spoil, spoil.build())
            seq.add(t0 + t_readout,
                    pp.make_label('LIN', 'SET', line),
                    pp.make_label('SLC', 'SET', slice_index))
            index += 1

    return sc.compile(
        seq, system, geometry=geometry, name='gre_2d',
        definitions={'TE': te_s, 'TR': tr_s, 'FlipAngle': flip_deg,
                     'RfSpoilIncrementDeg': 117.0, **ro.definitions()},
    )


def build_se(
    system: sc.System,
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
    invisible in a k-space extent check, because |k| is symmetric.
    """
    geometry = sc.Geometry(
        fov_mm=(fov_mm, fov_mm, slice_mm), matrix=(matrix, matrix, 1),
        slice_thickness_mm=slice_mm, n_slices=n_slices,
    )
    exc = sc.modules.SincExcitation(
        system, flip_deg=90, duration_us=2000, slice_thickness_mm=slice_mm, rephase=False)
    refoc = sc.modules.SincRefocusing(
        system, flip_deg=180, duration_us=4000, slice_thickness_mm=slice_mm,
        crusher_twists=4, crusher_voxel_mm=slice_mm)
    ro = sc.modules.CartesianLine(
        system, fov_ro_mm=fov_mm, matrix_ro=matrix, readout_duration_us=3200, prephase=False)
    pe = sc.modules.PhaseEncode(system, fov_pe_mm=fov_mm, matrix_pe=matrix)

    raster = system.block_raster_s
    winders = sc.raster.ceil_to(max(pe.duration, ro.prephase_duration), raster)
    first = sc.raster.ceil_to(
        (exc.duration - exc.isodelay) + winders + refoc.isodelay, raster)
    second = sc.raster.ceil_to((refoc.duration - refoc.isodelay) + ro.time_to_echo, raster)
    te_min = 2 * max(first, second)
    te = max(te_s, te_min)

    t_refoc = sc.raster.ceil_to(exc.isodelay + te / 2 - refoc.isodelay, raster)
    t_readout = sc.raster.ceil_to(exc.isodelay + te - ro.time_to_echo, raster)

    seq = sc.LogicBlock('se_2d')
    index = 0
    for slice_index in sc.interleaved_slice_order(n_slices):
        z = geometry.slice_positions_m[slice_index]
        for line in geometry.pe_lines:
            t0 = index * tr_s
            seq.add(t0, exc.build(slice_offset_m=z))
            seq.add(t0 + exc.duration,
                    pe.build(line=line - geometry.kspace_center_line, scale=-1.0))
            seq.add(t0 + exc.duration, ro.prephaser_block(polarity=-1))
            seq.add(t0 + t_refoc, refoc.build(slice_offset_m=z))
            seq.add(t0 + t_readout, ro.build())
            seq.add(t0 + t_readout,
                    pp.make_label('LIN', 'SET', line),
                    pp.make_label('SLC', 'SET', slice_index))
            index += 1

    return sc.compile(
        seq, system, geometry=geometry, name='se_2d',
        definitions={'TE': te, 'TR': tr_s, **ro.definitions()},
    )


def build_dti(
    system: sc.System,
    *,
    fov_mm: float = 240.0,
    matrix: int = MATRIX,
    slice_mm: float = 4.0,
    b_values: tuple[float, ...] = (0.0, 1000.0),
    directions: tuple[tuple[float, float, float], ...] | None = None,
    n_directions: int = 6,
    n_interleaves: int = 1,
    density: float = 0.5,
    n_slices: int = 2,
    slices_per_tr: int = 1,
    tr_s: float = 6.0,
    te_s: float | None = None,
    fat_sat: bool = False,
    spoil_axes: tuple[str, ...] = ('z',),
    rephase: bool = True,
    regime: str = 'default',
) -> CompiledSequence:
    """
    A spin-echo spiral DTI shot set.

    ``rephase`` is exposed so a test can show what happens without the excitation's slice rewinder:
    the refocusing pulse does **not** undo the excitation gradient's tail, because that pulse's own
    slice-select gradient is symmetric about its centre and its lead and tail cancel each other.
    """
    geometry = sc.Geometry(
        fov_mm=(fov_mm, fov_mm, 0.0), matrix=(matrix, matrix, 1),
        slice_thickness_mm=slice_mm, n_slices=n_slices,
    )
    exc = sc.modules.SincExcitation(
        system, flip_deg=90, duration_us=2000, slice_thickness_mm=slice_mm,
        rephase=rephase, regime=regime)
    refoc = sc.modules.SincRefocusing(
        system, flip_deg=180, duration_us=4000, slice_thickness_mm=slice_mm,
        crusher_twists=4, crusher_voxel_mm=slice_mm, regime=regime)
    readout = sc.modules.SpiralVDS(
        system, fov_mm=fov_mm, matrix=matrix, n_interleaves=n_interleaves,
        density=density, regime=regime)
    spoiler = sc.modules.Spoiler(
        system, twists=4, voxel_mm=slice_mm, axes=spoil_axes, regime=regime)
    fat = sc.modules.FatSat(system, voxel_mm=slice_mm, regime=regime) if fat_sat else None

    reference = sc.modules.MonopolarDiffusion(
        system, b_value_s_per_mm2=max(b_values),
        refocus_duration_us=refoc.duration * 1e6, regime=regime)
    encoders = {
        b: sc.modules.MonopolarDiffusion(
            system, b_value_s_per_mm2=b, refocus_duration_us=refoc.duration * 1e6,
            lobe_duration_us=reference.lobe_duration * 1e6, regime=regime)
        for b in b_values
    }

    raster = system.block_raster_s
    delta = reference.lobe_duration
    first = sc.raster.ceil_to((exc.duration - exc.isodelay) + delta + refoc.isodelay, raster)
    second = sc.raster.ceil_to(
        (refoc.duration - refoc.isodelay) + delta + readout.time_to_echo, raster)
    te_min = 2 * max(first, second)
    te = te_min if te_s is None else sc.raster.ceil_to(te_s, raster)
    if te < te_min - 1e-12:
        msg = (f'TE {te * 1e3:.2f} ms is below the minimum {te_min * 1e3:.2f} ms; the diffusion '
               f'lobe contributes {delta * 1e3:.2f} ms to each half')
        raise sc.ConfigurationError(msg)

    t_refoc = sc.raster.ceil_to(exc.isodelay + te / 2 - refoc.isodelay, raster)
    t_lobe1 = sc.raster.ceil_to(t_refoc - delta, raster)
    t_lobe2 = sc.raster.ceil_to(t_refoc + refoc.duration, raster)
    t_readout = sc.raster.ceil_to(exc.isodelay + te - readout.time_to_echo, raster)
    shot_s = sc.raster.ceil_to(t_readout + readout.duration + spoiler.duration, raster)
    prep_s = fat.duration if fat is not None else 0.0
    period_s = sc.raster.ceil_to(prep_s + shot_s, raster)
    if tr_s < period_s * slices_per_tr:
        msg = (f'TR {tr_s * 1e3:.0f} ms cannot hold {slices_per_tr} shot(s) of '
               f'{period_s * 1e3:.2f} ms')
        raise sc.ConfigurationError(msg)

    table = directions if directions is not None else sc.modules.dti_directions(n_directions)
    volumes = [(b, (0.0, 0.0, 1.0)) for b in b_values if b == 0.0]
    volumes += [(b, d) for b in b_values if b != 0.0 for d in table]

    plan = [
        (volume, b, direction, interleaf, slice_index)
        for volume, (b, direction) in enumerate(volumes)
        for interleaf in range(n_interleaves)
        for slice_index in sc.interleaved_slice_order(n_slices)
    ]

    seq = sc.LogicBlock('spiral_dti')
    for shot, (volume, b, direction, interleaf, slice_index) in enumerate(plan):
        period, position = divmod(shot, slices_per_tr)
        t0 = period * tr_s + position * period_s
        z = geometry.slice_positions_m[slice_index]
        diff = encoders[b]
        if fat is not None:
            seq.add(t0, fat.build())
        t = t0 + prep_s
        seq.add(t, exc.build(slice_offset_m=z))
        seq.add(t + t_lobe1, diff.build(part='pre', direction=direction))
        seq.add(t + t_refoc, refoc.build(slice_offset_m=z))
        seq.add(t + t_lobe2, diff.build(part='post', direction=direction))
        seq.add(t + t_readout, readout.build(interleaf=interleaf))
        seq.add(t + t_readout + readout.duration, spoiler.build())
        seq.add(t + t_readout,
                pp.make_label('SLC', 'SET', slice_index),
                pp.make_label('SET', 'SET', volume),
                pp.make_label('SEG', 'SET', interleaf))

    return sc.compile(
        seq, system, geometry=geometry, name='spiral_dti', regime=regime,
        definitions={
            'TE': te, 'TR': tr_s,
            'bValues': [float(b) for b in b_values],
            'AchievedbValues': [round(encoders[b].achieved_b_s_per_mm2(), 2) for b in b_values],
            'DiffusionScheme': 'monopolar',
            'DiffusionDirections': len(table),
            'DiffusionLobeDuration': reference.lobe_duration,
            'DiffusionBigDelta': reference.big_delta,
            'SpiralInterleaves': n_interleaves,
            'SpiralDensity': density,
            'SlicesPerTR': slices_per_tr,
        },
    )


@pytest.fixture(scope='module')
def gre(system: sc.System) -> CompiledSequence:
    return build_gre(system)


@pytest.fixture(scope='module')
def se(system: sc.System) -> CompiledSequence:
    return build_se(system)


@pytest.fixture(scope='module')
def dti() -> CompiledSequence:
    scanner = sc.System.preset('generic_3t').derate('dwi', grad=0.85, slew=0.65)
    return build_dti(scanner, regime='dwi', n_directions=6, n_slices=2)


@pytest.fixture(params=['gre', 'se', 'dti'])
def compiled(request) -> CompiledSequence:
    return request.getfixturevalue(request.param)
