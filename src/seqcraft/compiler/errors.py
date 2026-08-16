"""
The three exceptions only the compiler raises.

An exception lives with the code that raises it; only those raised from more than one package stay
at the root in :mod:`seqcraft.errors`.  That rule is what keeps the root module from becoming a
registry of everything that can go wrong anywhere -- a file you have to edit to add an error to a
package that does not otherwise touch it.

All three are re-exported as ``sc.CompileError``, ``sc.HardwareLimitError`` and
``sc.DefinitionConflict``, so the spelling a caller writes never depends on the layout.
"""

from __future__ import annotations

from ..errors import SeqCraftError

__all__ = ['CompileError', 'DefinitionConflict', 'HardwareLimitError']


class CompileError(SeqCraftError):
    """
    A logic-block tree cannot be expressed as legal pulseq blocks.

    Raised for everything the compiler cannot resolve by scheduling: two RF or ADC events
    overlapping in time, a negative absolute start, a gradient starting off the gradient raster,
    a mandatory gap with nothing in it that may be cut, a block boundary that would fall inside a
    gradient an ADC is sampling, two ADCs writing the same k-space address, and a sequence
    pypulseq's own timing check refuses.

    The message names the event, its provenance path, and usually two concrete remedies.
    """


class HardwareLimitError(SeqCraftError):
    """
    The machine cannot play it.

    The *summed* waveform exceeds ``max_grad`` or ``max_slew`` on an axis, or one ADC or RF event
    carries more samples than the vendor interpreter accepts.

    Summed, because that is the only place the truth is visible: two individually legal gradients
    on one axis can reach 189 % of the slew limit together, and no module can see that in
    isolation.  The vector norm across simultaneous axes is *not* this -- it routinely exceeds the
    per-axis limit and is legal on real amplifiers, so it warns.
    """


class DefinitionConflict(SeqCraftError):
    """
    Two sources claimed the same ``[DEFINITIONS]`` key with different values.

    Last-writer-wins is how a ``.seq`` came to say ``kSpaceCenterLine = 73.0`` while its own
    navigator used 36.5: neither source was wrong about itself, and nothing compared them.
    """
