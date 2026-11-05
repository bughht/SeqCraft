"""Direct tests of the PlacedEvent-to-ready-block legalization boundary."""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

from seqcraft.compiler.legalization import legalize_blocks, superpose
from seqcraft.compiler.model import PlacedEvent
from seqcraft.compiler.verification import verify_ready_blocks
from seqcraft.design.timing import Raster


def _placed(event, path: tuple[str, ...]) -> PlacedEvent:
    """Build a synthetic placed gradient without involving tree traversal."""
    start = float(event.delay)
    end = float(pp.calc_duration(event))
    return PlacedEvent(0.0, start, end, start, end, event, path)


def test_legalization_returns_complete_ready_blocks_and_explicit_notes(opts) -> None:
    """Same-axis transformation is visible in the result and needs no Sequence."""
    shape = {'duration': 2e-3, 'rise_time': 400e-6, 'system': opts}
    first = _placed(pp.make_trapezoid('x', area=10.0, **shape), ('tr', 'rewinder'))
    second = _placed(pp.make_trapezoid('x', area=20.0, **shape), ('tr', 'readout'))
    raster = Raster(float(opts.block_duration_raster), 'block')

    result = legalize_blocks((0.0, 2e-3), (first, second), {}, opts, raster)

    assert len(result.blocks) == 1
    assert result.blocks[0].kinds == ('grad',)
    assert result.blocks[0].origin == ('tr',)
    assert verify_ready_blocks(result.blocks, expected_start=0.0) == ()
    assert result.notes == (
        ('merge', ('tr.readout+tr.rewinder (axis x)',)),
    )


def _placed_at(event, node_t: float, path: tuple[str, ...]) -> PlacedEvent:
    """A synthetic placed gradient at `node_t`."""
    start = node_t + float(event.delay)
    end = node_t + float(pp.calc_duration(event))
    return PlacedEvent(node_t, start, end, start, end, event, path)


@pytest.mark.parametrize('t0', [1.0, 10.0, 100.0, 1000.0])
def test_superposition_is_the_same_wherever_the_block_sits(opts, t0: float) -> None:
    """
    Values interpolated at a cut, or at another piece's knot, do not depend on the block's position.

    A trapezoid, whose corners sit on raster edges, overlaps an arbitrary gradient, whose samples
    sit at raster centres, and both are cut by the interval.  In ticks from the block's start the
    interpolation is a difference of exact integers, so the knots come out the same bit for bit.
    """
    s = np.linspace(0.0, 1.0, 60)
    curve = np.sin(np.pi * s) * np.cos(2.0 * np.pi * s) * 0.3 * float(opts.max_grad)
    spiral = pp.make_arbitrary_grad('x', waveform=curve, first=0.0, last=0.0, system=opts)
    trap = pp.make_trapezoid('x', area=20.0, duration=300e-6, system=opts)

    def knots(at: float) -> tuple[np.ndarray, np.ndarray]:
        pieces = [_placed_at(spiral, at, ('spiral',)), _placed_at(trap, at + 120e-6, ('trap',))]
        return superpose(pieces, at + 100e-6, at + 500e-6)

    near, far = knots(0.0), knots(t0)
    assert np.array_equal(near[0], far[0])
    assert np.array_equal(near[1], far[1])
