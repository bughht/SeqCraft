"""
One small, complete acquisition, written the way the authoring skill says to write one.

This is a **test fixture, not an example**: it exists so the deterministic checks in
``test_emitted_acquisition.py`` have a real emitted ``.seq`` to inspect.  It is deliberately
tiny -- an 8x8 matrix is not a protocol anyone would scan -- because what is under test is the
authoring and validation contract, not the image.

What makes it *complete* rather than a kernel demonstration is the part above the repetition: an
explicit dummy-repetition policy, a stated view order, and a recorded table that the test checks
the emitted file actually walks.  Those are the things the compiler does not check.

Every number here is a deliberate, stated choice.  Nothing in this file is read back by the test
as if it were evidence -- the test derives its expectations from the physics and measures the
emitted events.
"""

from __future__ import annotations

from pathlib import Path

import pypulseq as pp

import seqcraft as sc

#: Conservative illustrative design limits.  Not a scanner specification.
OPTS = pp.Opts(
    max_grad=28, grad_unit='mT/m',
    max_slew=120, slew_unit='T/m/s',
    rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
    B0=3.0,
)

FOV_MM = (200.0, 200.0)
MATRIX = (8, 8)
THICKNESS_MM = 5.0
FLIP_DEG = 15.0
BANDWIDTH_HZ_PX = 500.0
DUMMIES = 2                 # unacquired repetitions, so the first acquired line is not the first


def view_order(ny: int) -> list[int]:
    """Linear, low-to-high.  Stated rather than implied, because the test checks it."""
    return list(range(ny))


def build(path: Path) -> dict[str, object]:
    """Write the acquisition to `path` and return what it claims, for the test to check against."""
    kernel = sc.modules.GRE2DTR(
        opts=OPTS, fov_mm=FOV_MM, matrix=MATRIX, thickness_mm=THICKNESS_MM,
        flip_deg=FLIP_DEG, bandwidth_hz_px=BANDWIDTH_HZ_PX,
    )
    lines = view_order(MATRIX[1])

    tree, at = sc.LogicBlock('reference_acquisition'), 0.0
    for _ in range(DUMMIES):
        tree.add(at, kernel(line=kernel.center_line, acquire=False))
        at += kernel.tr_s
    for line in lines:
        tree.add(at, kernel(line=line))
        at += kernel.tr_s

    seq = sc.compile(tree, OPTS, name='reference_acquisition', definitions={
        'FOV': [FOV_MM[0] * 1e-3, FOV_MM[1] * 1e-3, THICKNESS_MM * 1e-3],
        'Matrix': list(MATRIX),
        'TE': kernel.te_s,
        'TR': kernel.tr_s,
        'FlipAngle': FLIP_DEG,
    })
    seq.write(str(path))

    return {
        'kernel': kernel,
        'tree': tree,
        'lines': lines,
        'dummies': DUMMIES,
        'matrix': MATRIX,
        'fov_mm': FOV_MM,
    }
