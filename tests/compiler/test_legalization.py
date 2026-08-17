"""Direct tests of the PlacedEvent-to-ready-block legalization boundary."""

from __future__ import annotations

import pypulseq as pp

from seqcraft.compiler.legalization import legalize_blocks
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
