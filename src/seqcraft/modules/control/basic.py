"""
Control modules: delays, triggers, barriers and raw events.

These carry no encoding.  They exist so that everything a sequence needs to say can be said as a
module, and so the escape hatch out of seqcraft is a first-class object rather than a hole.

Examples
--------
>>> import seqcraft as sc
>>> system = sc.System.preset('generic_3t')
>>> Delay(system, duration_ms=20).build()
LogicBlock(delay, 1 node, 20.00 ms)
>>> Trigger(system, channel='ext1').build()
LogicBlock(trigger, 1 node, 0.01 ms)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pypulseq as pp

from ...core import events as ev
from ...core.errors import ConfigurationError, format_error
from ...core.logic import LogicBlock, barrier
from ...core.module import Module
from ...core.raster import ceil_to
from ...core.registry import register
from ...core.validate import require_in, require_positive

if TYPE_CHECKING:
    from types import SimpleNamespace

    from ...core.system import System

__all__ = ['Barrier', 'Delay', 'RawEvents', 'Trigger']


@register()
class Delay(Module):
    """
    A block that occupies time and does nothing.

    Parameters
    ----------
    system
        The scanner.
    duration_ms
        How long, milliseconds.  Rounded up to the block raster.

    Notes
    -----
    Rarely needed for TR padding -- placing the next TR at ``(i + 1) * tr`` already leaves the
    gap, and the compiler fills it.  It is useful when the *gap itself* is the point: an
    inversion time, a saturation recovery delay, a magnetisation-preparation interval.

    Examples
    --------
    >>> import seqcraft as sc
    >>> ti = Delay(sc.System.preset('generic_3t'), duration_ms=250)
    >>> round(ti.duration * 1e3)
    250
    """

    def __init__(self, system: System, *, duration_ms: float, regime: str = 'default') -> None:
        super().__init__(system, regime=regime)
        self.duration_ms = float(duration_ms)
        require_positive(self, 'duration_ms')

    @property
    def duration(self) -> float:
        """Seconds occupied, rounded up to the block raster."""
        return ceil_to(self.duration_ms / 1e3, self.system.block_raster_s)

    def build(self) -> LogicBlock:
        """Return a block holding a single delay event."""
        return LogicBlock('delay').add(0.0, pp.make_delay(self.duration))


@register()
class Trigger(Module):
    """
    An external trigger or output pulse.

    Parameters
    ----------
    system
        The scanner.
    channel
        Pulseq trigger channel: ``'physio1'``, ``'physio2'``, ``'ext1'`` for inputs, or
        ``'osc0'``, ``'osc1'`` for outputs.
    duration_us
        Pulse length, microseconds.
    delay_us
        Delay from the start of its block.

    Notes
    -----
    A trigger attaches to whichever compiled block contains its start, so it does **not** force a
    block of its own.  When it must land on its own block -- to make an oscilloscope trace
    unambiguous, for instance -- put a :class:`Barrier` beside it.

    pypulseq builds inputs and outputs with different functions: ``make_trigger`` for the physio
    channels the scanner *waits on*, and ``make_digital_output_pulse`` for the oscilloscope
    channels it *drives*.  Both are triggers as far as a sequence is concerned, so this module
    dispatches on the channel name rather than making the caller know which is which.

    Examples
    --------
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> round(Trigger(system, channel='osc0', duration_us=100).duration * 1e6)
    100
    >>> Trigger(system, channel='physio1').is_input
    True
    """

    #: Channels the scanner waits on, built with ``make_trigger``.
    INPUTS = ('physio1', 'physio2')
    #: Channels the scanner drives, built with ``make_digital_output_pulse``.
    OUTPUTS = ('osc0', 'osc1', 'ext1')

    def __init__(
        self,
        system: System,
        *,
        channel: str = 'osc0',
        duration_us: float = 10.0,
        delay_us: float = 0.0,
        regime: str = 'default',
    ) -> None:
        super().__init__(system, regime=regime)
        self.channel = str(channel)
        self.duration_us = float(duration_us)
        self.delay_us = float(delay_us)
        require_in(self, 'channel', (*self.INPUTS, *self.OUTPUTS))
        require_positive(self, 'duration_us')
        make = pp.make_trigger if self.is_input else pp.make_digital_output_pulse
        self.event = make(
            channel=self.channel,
            duration=self.duration_us / 1e6,
            delay=self.delay_us / 1e6,
            system=self.opts,
        )

    @property
    def is_input(self) -> bool:
        """``True`` for a channel the scanner waits on rather than drives."""
        return self.channel in self.INPUTS

    @property
    def duration(self) -> float:
        """Seconds occupied, delay included."""
        return ceil_to((self.delay_us + self.duration_us) / 1e6, self.system.block_raster_s)

    def build(self) -> LogicBlock:
        """Return the trigger event."""
        return LogicBlock('trigger').add(0.0, self.event)


@register()
class Barrier(Module):
    """
    Forces a pulseq block boundary where it is placed.

    Parameters
    ----------
    system
        The scanner.  Unused, but kept so barriers construct like every other module.
    tag
        Recorded in the compile report, so a surprising extra block can be traced back.

    Notes
    -----
    The compiler chooses boundaries from where RF and ADC events fall, which is right almost
    always.  This is the escape hatch for when it is not: a trigger that must sit alone, or a
    gradient you want split at a known instant so a later reconstruction step can find the seam.

    It occupies no time, so placing one never changes the sequence's duration -- only its block
    structure.

    Examples
    --------
    >>> import seqcraft as sc
    >>> b = Barrier(sc.System.preset('generic_3t'), tag='midpoint')
    >>> b.duration
    0.0
    >>> b.build()
    LogicBlock(barrier, 1 node, 0.00 ms)
    """

    def __init__(self, system: System, *, tag: str = 'barrier', regime: str = 'default') -> None:
        super().__init__(system, regime=regime)
        self.tag = str(tag)

    @property
    def duration(self) -> float:
        """Always zero: a barrier changes block structure, not timing."""
        return 0.0

    def build(self) -> LogicBlock:
        """Return the barrier marker."""
        return LogicBlock('barrier').add(0.0, barrier(self.tag))


@register()
class RawEvents(Module):
    """
    Arbitrary pypulseq events wrapped as a module.

    The pressure-relief valve.  When seqcraft has no module for what you need, build the events
    with pypulseq directly and wrap them here: they take part in compilation, provenance and
    validation exactly like anything else, because a node *is* an event.

    Parameters
    ----------
    system
        The scanner.
    events
        Pulseq events, each with an optional ``at_s`` offset supplied through `starts`.
    starts
        Start time of each event within the block, seconds.  Defaults to all zero.
    tag
        Tag for the built block.

    Notes
    -----
    Its existence in the documentation is what makes the abstraction safe to adopt: nothing
    forces you to wait for a module to be written.  The compile report notes how many blocks came
    from raw events, so drift towards hand-built sequences stays visible rather than silent.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> g = pp.make_trapezoid('x', area=100.0, system=system.default)
    >>> raw = RawEvents(system, events=(g,), tag='my_gradient')
    >>> raw.build()
    LogicBlock(my_gradient, 1 node, 0.26 ms)
    """

    def __init__(
        self,
        system: System,
        *,
        events: tuple[SimpleNamespace, ...],
        starts: tuple[float, ...] | None = None,
        tag: str = 'raw',
        regime: str = 'default',
    ) -> None:
        super().__init__(system, regime=regime)
        self.tag = str(tag)
        self.events = tuple(events)
        self.starts = tuple(starts) if starts is not None else (0.0,) * len(self.events)
        if len(self.starts) != len(self.events):
            msg = format_error(
                'starts must have one entry per event.',
                {'events': len(self.events), 'starts': len(self.starts)},
            )
            raise ConfigurationError(msg)
        for e in self.events:
            if getattr(e, 'type', None) is None:
                msg = format_error(
                    f'RawEvents was given a {type(e).__name__}, not a pulseq event.',
                    {},
                    ['build events with pp.make_trapezoid / make_adc / make_sinc_pulse / ...'],
                )
                raise ConfigurationError(msg)

    @property
    def duration(self) -> float:
        """Seconds occupied, measured from the events."""
        return max(
            (start + float(pp.calc_duration(e)) for start, e in zip(self.starts, self.events)),
            default=0.0,
        )

    def build(self) -> LogicBlock:
        """Return the events at their given starts, sanitised of any registration state."""
        out = LogicBlock(self.tag)
        for start, event in zip(self.starts, self.events):
            out.add(start, ev.sanitise(event))
        return out
