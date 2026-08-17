"""Mechanical emission of ready blocks into a :class:`pypulseq.Sequence`."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import pypulseq as pp

from ..errors import format_error
from .errors import CompileError
from .model import PulseqReadyBlock

if TYPE_CHECKING:
    from collections.abc import Iterable

__all__ = ['emit_blocks']


def emit_blocks(seq: Any, blocks: Iterable[PulseqReadyBlock]) -> None:
    """
    Add each validated ready block to `seq` without making scheduling decisions.

    Legalization owns event assignment, gradient transformation, limits and block duration. This
    stage only states each explicit duration, calls PyPulseq and adds block context to a rejection.
    """
    for ready in blocks:
        if not ready.events:
            seq.add_block(pp.make_delay(ready.duration))
            continue

        # add_block() takes no duration argument. A delay event is how pulseq states an explicit
        # block length, and set_block() takes the max of it and the event extents.
        try:
            seq.add_block(*ready.events, pp.make_delay(ready.duration))
        except (ValueError, RuntimeError) as err:
            origin = ', '.join(
                sorted({'.'.join(path) for path in ready.source_paths if path})
            ) or '?'
            msg = format_error(
                f'pypulseq rejected block {ready.index} at {ready.start * 1e6:.1f} us: {err}',
                {
                    'duration_us': f'{ready.duration * 1e6:.1f}',
                    'events': ', '.join(getattr(event, 'type', '?') for event in ready.events),
                    'from': origin,
                },
                ['this is a compiler bug unless the tree contains raw events built against '
                 'a different Opts -- please report it with the tree that produced it'],
            )
            raise CompileError(msg) from err
