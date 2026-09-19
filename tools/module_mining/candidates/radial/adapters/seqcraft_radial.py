"""
The extracted ``RadialReadout``, measured by the same code that measured the references.

A spoke is not a sequence, so the adapter wraps it in the smallest tree that compiles: **one
spoke**, on its own.  No excitation, no spoiler, no TR -- anything else would be measuring a
composition rather than the leaf.

One spoke per sequence rather than a stack of them, and that is not a style choice.  ``k`` is
integrated from the start of the sequence and reset only by an excitation pulse, so a stack of
bare spokes has each one starting wherever the last ended: the first measurement of this
candidate reported a rotation error of 5.8e+02 1/m for that reason, and the module was fine.
A harness that quietly sums its own trajectories is worse than no harness.
"""

from __future__ import annotations

from typing import Any

import seqcraft as sc

from ....reference import ReferenceSequence

PROVENANCE = {
    'repo': 'bughht/SeqCraft',
    'path': 'src/seqcraft/modules/readout/radial_readout.py',
    'license': 'MIT',
    'role': 'candidate',
}


def build(*, angle_rad: float, opts: Any = None, fov_mm: float = 260.0,
          matrix: int = 64, dwell_s: float = 20e-6, partial_fourier: float = 1.0,
          **overrides: Any) -> ReferenceSequence:
    """Build and compile one spoke at `angle_rad`, alone."""
    import pypulseq as pp  # noqa: PLC0415

    scanner = opts if opts is not None else pp.Opts(
        max_grad=28, grad_unit='mT/m', max_slew=120, slew_unit='T/m/s',
        rf_dead_time=100e-6, rf_ringdown_time=20e-6, adc_dead_time=10e-6,
    )
    spoke = sc.modules.RadialReadout(
        opts=scanner, fov_mm=fov_mm, matrix=matrix, dwell_s=dwell_s,
        partial_fourier=partial_fourier, **overrides)

    tree = sc.LogicBlock('spoke').add(0.0, spoke(angle_rad=float(angle_rad)))

    return ReferenceSequence(
        name='seqcraft_radial',
        sequence=sc.compile(tree, scanner, name='module_mining_radial'),
        tree=tree,
        parameters={'fov_mm': fov_mm, 'matrix': matrix, 'dwell_s': dwell_s,
                    'partial_fourier': partial_fourier, 'angle_rad': float(angle_rad)},
        provenance=PROVENANCE,
        # Unlike every external reference, this one *states* its geometry.  The comparator checks
        # the claims against the trajectory rather than taking them.
        semantic={
            'num_samples': spoke.num_samples,
            'center_sample': spoke.center_sample,
            'dk_per_m': spoke.dk_per_m,
            'k_first_per_m': spoke.k_first_per_m,
            'k_last_per_m': spoke.k_last_per_m,
            'k_max_per_m': spoke.k_max_per_m,
            'time_to_center_s': spoke.time_to_center(),
        },
    )
