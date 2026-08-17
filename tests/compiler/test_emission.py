"""Direct tests of ready-block-to-PyPulseq mechanical emission."""

from __future__ import annotations

from types import SimpleNamespace

import pypulseq as pp
import pytest

from seqcraft.compiler.emission import emit_blocks
from seqcraft.compiler.errors import CompileError
from seqcraft.compiler.model import PulseqReadyBlock


def test_synthetic_ready_blocks_emit_without_placement_or_legalization(opts) -> None:
    """Emission consumes only ready blocks and preserves their explicit durations."""
    gradient = pp.make_trapezoid('x', area=10.0, duration=1e-3, system=opts)
    blocks = (
        PulseqReadyBlock(0, 0.0, 1e-3, 1e-3, (), (), ()),
        PulseqReadyBlock(
            1,
            1e-3,
            2e-3,
            1e-3,
            (gradient,),
            (('tr', 'readout'),),
            ('tr', 'readout'),
        ),
    )
    sequence = pp.Sequence(system=opts)

    emit_blocks(sequence, blocks)

    assert len(sequence.block_events) == 2
    assert sequence.duration()[0] == pytest.approx(2e-3)
    assert sequence.get_block(2).gx is not None


def test_pypulseq_rejection_is_mapped_with_ready_block_context() -> None:
    """The mechanical stage adds block and provenance context to a backend failure."""
    ready = PulseqReadyBlock(
        3,
        4e-3,
        5e-3,
        1e-3,
        (SimpleNamespace(type='grad'),),
        (('tr', 'readout'),),
        ('tr', 'readout'),
    )

    class RejectingSequence:
        def add_block(self, *events) -> None:
            raise ValueError('backend refused it')  # noqa: TRY003 - simulates third-party text

    with pytest.raises(CompileError) as error:
        emit_blocks(RejectingSequence(), (ready,))

    message = str(error.value)
    assert 'block 3 at 4000.0 us' in message
    assert 'tr.readout' in message
    assert 'backend refused it' in message
