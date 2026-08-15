"""
The data model: :class:`LogicBlock`, a tree of pulseq events with relative start times.

A logic block is a list of children, each with a start time measured from the block's own
start.  A child is either a pulseq event or another logic block.  That is the whole model.

Overlap is legal everywhere in the tree -- two gradients on one axis, a gradient across an
RF, anything.  Making it legal *for pulseq* is :mod:`seqcraft.compiler`'s job.  Keeping
those two concerns apart is what lets a module say what it means (a diffusion lobe here, a
readout there) without also having to know where pulseq's block boundaries will land.

Three things this deliberately does not have:

**No marks or anchors.**  A block cannot know where it sits, so it cannot know when its echo
occurs.  The *module* knows, because the module has the domain knowledge, and it says so as a
plain property: ``readout.time_to_echo``.  Pinning that into the block would also stop the
same block being reused elsewhere.

**No labels field.**  A pulseq label is an event, so it is a node like any other.

**No declared duration.**  :attr:`LogicBlock.duration` is measured from the nodes, so a block
cannot advertise a length that disagrees with what it plays.  To make one longer, add a delay
event -- which is how you lengthen a block in pulseq anyway.

Examples
--------
>>> import pypulseq as pp
>>> from pypulseq.opts import Opts
>>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=170, slew_unit='T/m/s')
>>> gx = pp.make_trapezoid(channel='x', area=100.0, system=o)
>>> lb = LogicBlock('readout')
>>> lb.add(0.0, gx)
LogicBlock(readout, 1 node, 0.24 ms)
>>> round(lb.duration * 1e6)                           # measured, not declared
240
>>> outer = LogicBlock('tr').add(0.0, lb).add(1e-3, lb)
>>> round(outer.duration * 1e6)                        # 1000 + 240
1240
>>> [round(t * 1e6) for t, _, _ in flatten(outer)]     # absolute times
[0, 1000]
>>> [path for _, _, path in flatten(outer)]            # provenance is the tree
[('tr', 'readout'), ('tr', 'readout')]
"""

from __future__ import annotations

import numbers
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Union, overload

import pypulseq as pp

from ..errors import ConfigurationError, format_error

if TYPE_CHECKING:
    from collections.abc import Iterator
    from types import SimpleNamespace

    from typing_extensions import Self

__all__ = ['BARRIER', 'Item', 'LogicBlock', 'Node', 'barrier', 'flatten', 'span']

#: A child of a logic block: a pulseq event or a nested block.
Item = Union['SimpleNamespace', 'LogicBlock']

#: ``type`` value of the pseudo-event produced by :func:`barrier`.
BARRIER = 'seqcraft_barrier'


@dataclass
class Node:
    """
    One child of a :class:`LogicBlock`.

    Parameters
    ----------
    start
        Seconds from the start of the enclosing block.  Called `start` rather than "offset"
        because in MR an offset is a frequency or a phase, and reusing the word invites
        exactly the wrong reading.
    item
        A pulseq event or a nested :class:`LogicBlock`.
    """

    start: float
    item: Item


class LogicBlock:
    """
    A tree of pulseq events and nested blocks, each with a start time.

    Parameters
    ----------
    tag
        Optional label, used in plots, error messages and the compiled block's provenance
        path.  Nothing depends on it and it is never required.

    Attributes
    ----------
    nodes : list[Node]
        The children, in insertion order.  A plain list on purpose: create with
        :meth:`add`, read with ``lb.nodes[i]`` or iteration, update with
        ``lb.nodes[i].start += dt``, delete with ``del lb.nodes[i]``.  CRUD is Python's,
        not seqcraft's.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=170, slew_unit='T/m/s')
    >>> g = pp.make_trapezoid(channel='y', area=50.0, system=o)
    >>> lb = LogicBlock('blip')
    >>> _ = lb.add(0.0, g)
    >>> len(lb)
    1
    >>> lb.nodes[0].start += 100e-6            # update
    >>> round(lb.duration * 1e6)
    280
    >>> del lb.nodes[0]                        # delete
    >>> lb.duration
    0.0
    """

    __slots__ = ('nodes', 'tag')

    def __init__(self, tag: str = '') -> None:
        self.tag = tag
        self.nodes: list[Node] = []

    # ------------------------------------------------------------------------ create
    @overload
    def add(self, start: float, *items: Item) -> Self: ...

    @overload
    def add(self, start: list[Any] | tuple[Any, ...], /) -> Self: ...

    def add(self, start: float | list[Any] | tuple[Any, ...], *items: Item) -> Self:
        """
        Add events or nested blocks, starting `start` seconds into this block.

        Several items at one time is the common case -- an RF and its slice-select gradient,
        a readout gradient and its ADC -- so `add` takes any number, the way pulseq's
        ``add_block`` does.

        There is a second call shape for when the schedule is *computed* rather than written
        out: pass a list (or tuple) of ``[time, *items]`` rows and the whole table goes in
        one call.  ``nodes`` ends up exactly as the chained form would leave it -- rows are
        appended in the order given, items within a row keep theirs, and nothing is sorted
        by time.  It is input-side sugar over an unchanged data model.

        Parameters
        ----------
        start
            Seconds from this block's start.  May be zero; negative is allowed here and
            rejected by the compiler if it makes an absolute time negative.

            A ``list`` or ``tuple`` selects the batch form instead, and is then the only
            argument.  It is either a table of rows or -- when its first element is a
            number -- a single bare row, so the outer brackets are optional for one line.
            ``list``/``tuple`` and not "any iterable" on purpose: a :class:`LogicBlock` is
            itself iterable and must never be mistaken for a table.
        *items
            Pulseq events and/or nested :class:`LogicBlock` objects.  Not accepted
            alongside a table.

        Returns
        -------
        LogicBlock
            ``self``, so calls may chain.  It is the same block, not a copy.

        Raises
        ------
        ConfigurationError
            If an item is neither a pulseq event nor a :class:`LogicBlock`, if the first
            argument is an item rather than a time, or if a batch row is not
            ``[time, *items]``.  Catching this here rather than at compile time means the
            traceback points at the line that added the wrong thing, and a batch message
            names the row so a thirty-row table points at the offending one.

        Examples
        --------
        >>> import pypulseq as pp
        >>> from pypulseq.opts import Opts
        >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=170, slew_unit='T/m/s')
        >>> rf = pp.make_block_pulse(flip_angle=1.57, duration=500e-6, system=o,
        ...                          use='excitation')
        >>> gz = pp.make_trapezoid(channel='z', area=100.0, system=o)
        >>> lb = LogicBlock('exc').add(0.0, rf)
        >>> len(lb)
        1

        The batch form, and the chained form it is equal to:

        >>> table = LogicBlock('exc').add([[0.0, rf, gz],
        ...                                [1.2e-3, gz]])
        >>> chained = LogicBlock('exc').add(0.0, rf, gz).add(1.2e-3, gz)
        >>> [(round(n.start * 1e6), n.item.type) for n in table]
        [(0, 'rf'), (0, 'trap'), (1200, 'trap')]
        >>> [n.start for n in table] == [n.start for n in chained]
        True

        A single row needs no outer brackets:

        >>> len(LogicBlock('exc').add([0.0, rf, gz]))
        2
        """
        plan: list[tuple[int | None, Any, tuple[Item, ...]]]

        if isinstance(start, (list, tuple)):
            if items:
                msg = format_error(
                    f'LogicBlock.add() in its batch form takes exactly one argument, got '
                    f'{1 + len(items)}.',
                    {'tag': self.tag or '(untagged)', 'extra_args': len(items)},
                    [
                        'put every row inside the table: add([[t0, rf, gz], [t1, gzr]])',
                        'or add one instant at a time: add(t0, rf, gz)',
                    ],
                )
                raise ConfigurationError(msg)
            # A first element that is a number means this is one bare row, not a table.
            rows = [start] if start and isinstance(start[0], numbers.Real) else start
            plan = []
            for index, row in enumerate(rows):
                if not isinstance(row, (list, tuple)) or not row:
                    msg = format_error(
                        f'LogicBlock.add() batch row {index} is not a [time, *items] list, '
                        f'got {type(row).__name__}.',
                        {'tag': self.tag or '(untagged)', 'row': index},
                        [
                            'every row is [start_seconds, event, ...] -- e.g. [1.2e-3, rf, gz]',
                            'a row may hold a time and nothing else, but never nothing at all',
                        ],
                    )
                    raise ConfigurationError(msg)
                if not isinstance(row[0], numbers.Real):
                    msg = format_error(
                        f'LogicBlock.add() batch row {index} starts with '
                        f'{type(row[0]).__name__}, not a time.',
                        {'tag': self.tag or '(untagged)', 'row': index},
                        ['a row is [start_seconds, *items], so the time comes first'],
                    )
                    raise ConfigurationError(msg)
                plan.append((index, row[0], tuple(row[1:])))
        else:
            if isinstance(start, LogicBlock) or getattr(start, 'type', None) is not None:
                msg = format_error(
                    f'LogicBlock.add() takes the start time first, got '
                    f'{type(start).__name__}.',
                    {'tag': self.tag or '(untagged)'},
                    [
                        'add(0.0, item) places it at the start of this block',
                        'a computed schedule is a table of rows: add([[t0, rf], [t1, gz]])',
                    ],
                )
                raise ConfigurationError(msg)
            plan = [(None, start, items)]

        for row_index, when, row_items in plan:
            for item in row_items:
                if not isinstance(item, LogicBlock) and getattr(item, 'type', None) is None:
                    fields: dict[str, object] = {'tag': self.tag or '(untagged)'}
                    if row_index is not None:
                        fields['row'] = row_index
                    fields['start_us'] = when * 1e6
                    msg = format_error(
                        f'LogicBlock.add() takes pulseq events or LogicBlocks, got '
                        f'{type(item).__name__}.',
                        fields,
                        [
                            'a pulseq event comes from pp.make_trapezoid / make_sinc_pulse / '
                            'make_adc / make_label / make_delay',
                            'to nest, pass another LogicBlock -- usually module.build()',
                        ],
                    )
                    raise ConfigurationError(msg)
                self.nodes.append(Node(float(when), item))
        return self

    def copy(self) -> LogicBlock:
        """
        Return a new block with a new node list, sharing the same items.

        Adding one block object at many times shares it, which is normally what you want --
        1700 references to one spoiler rather than 1700 copies.  Use this when you want a
        variant instead: mutating the copy's ``nodes`` leaves the original alone.

        Examples
        --------
        >>> a = LogicBlock('a')
        >>> b = a.copy()
        >>> b.nodes is a.nodes
        False
        """
        out = LogicBlock(self.tag)
        out.nodes = [Node(n.start, n.item) for n in self.nodes]
        return out

    # -------------------------------------------------------------------------- read
    @property
    def duration(self) -> float:
        """
        Seconds from this block's start to the end of its last node.

        Measured, never declared, so a block cannot claim a length it does not play.  To
        make a block longer than its contents -- a b=0 diffusion volume that must occupy
        the same slot as an encoded one -- add a delay event:

        >>> import pypulseq as pp
        >>> lb = LogicBlock('b0').add(0.0, pp.make_delay(4.2e-3))
        >>> round(lb.duration * 1e3, 1)
        4.2
        """
        return max((n.start + span(n.item) for n in self.nodes), default=0.0)

    def __len__(self) -> int:
        """The number of direct children (not the total number of leaf events)."""
        return len(self.nodes)

    def __iter__(self) -> Iterator[Node]:
        """Iterate over the direct children."""
        return iter(self.nodes)

    def __repr__(self) -> str:
        """One line: tag, child count, duration."""
        label = self.tag or 'untagged'
        plural = '' if len(self.nodes) == 1 else 's'
        return f'LogicBlock({label}, {len(self.nodes)} node{plural}, {self.duration * 1e3:.2f} ms)'

    def describe(self, indent: int = 0) -> str:
        """
        Return a multi-line tree rendering, for debugging and notebooks.

        Examples
        --------
        >>> import pypulseq as pp
        >>> from pypulseq.opts import Opts
        >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=170, slew_unit='T/m/s')
        >>> g = pp.make_trapezoid(channel='x', area=100.0, system=o)
        >>> inner = LogicBlock('pre').add(0.0, g)
        >>> print(LogicBlock('ro').add(0.0, inner).describe())
        ro  0.24 ms
          +0.0 us  pre  0.24 ms
            +0.0 us  trap x
        """
        pad = '  ' * indent
        if indent == 0:
            lines = [f'{self.tag or "untagged"}  {self.duration * 1e3:.2f} ms']
        else:
            lines = []
        for n in sorted(self.nodes, key=lambda n: n.start):
            head = f'{pad}  +{n.start * 1e6:.1f} us  '
            if isinstance(n.item, LogicBlock):
                lines.append(f'{head}{n.item.tag or "untagged"}  {n.item.duration * 1e3:.2f} ms')
                lines.append(n.item.describe(indent + 1))
            else:
                kind = getattr(n.item, 'type', '?')
                ch = f' {n.item.channel}' if hasattr(n.item, 'channel') else ''
                lines.append(f'{head}{kind}{ch}')
        return '\n'.join(line for line in lines if line)


def span(item: Item) -> float:
    """
    Return how long `item` occupies, in seconds.

    For a nested block that is its :attr:`~LogicBlock.duration`; for an event it is
    ``pp.calc_duration``, which already includes the event's own ``delay`` and any ringdown
    or dead time pypulseq accounts for.

    Examples
    --------
    >>> import pypulseq as pp
    >>> round(span(pp.make_delay(2e-3)) * 1e3, 1)
    2.0
    """
    if isinstance(item, LogicBlock):
        return item.duration
    if getattr(item, 'type', None) == BARRIER:
        return 0.0
    return float(pp.calc_duration(item))


def flatten(
    root: LogicBlock,
    t0: float = 0.0,
    path: tuple[str, ...] = (),
) -> Iterator[tuple[float, SimpleNamespace, tuple[str, ...]]]:
    """
    Walk a tree, yielding every leaf event with its absolute start time and tag path.

    This is the compiler's only entry point into the model, and the reason the tree needs no
    methods of its own: reading a tree is not part of being one.

    Parameters
    ----------
    root
        The block to walk.
    t0
        Absolute time of `root`'s start.  Recursion uses it; callers normally do not.
    path
        Enclosing tags accumulated so far.  Recursion uses it; callers normally do not.

    Yields
    ------
    start, event, path
        `start` is seconds from `t0`, summed down the tree.  It does **not** include the
        event's own ``delay``; the compiler adds that, because an RF's delay carries the
        transmit dead time and an ADC's carries the receive dead time, and those are
        physical constraints belonging to the event rather than positioning decisions.
        `path` is the chain of non-empty tags enclosing the event, so provenance is the
        tree and needs no bookkeeping.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=170, slew_unit='T/m/s')
    >>> g = pp.make_trapezoid(channel='x', area=10.0, system=o)
    >>> deep = LogicBlock('a').add(0.0, LogicBlock('b').add(1e-3, LogicBlock('c').add(2e-3, g)))
    >>> [(round(t * 1e3, 1), p) for t, _, p in flatten(deep)]
    [(3.0, ('a', 'b', 'c'))]
    """
    here = (*path, root.tag) if root.tag else path
    for node in root.nodes:
        if isinstance(node.item, LogicBlock):
            yield from flatten(node.item, t0 + node.start, here)
        else:
            yield (t0 + node.start, node.item, here)


def barrier(tag: str = 'barrier') -> SimpleNamespace:
    """
    Return a zero-duration pseudo-event that forces a pulseq block boundary where it sits.

    The compiler chooses boundaries from where RF and ADC events fall, which is right almost
    always.  This is the escape hatch for the times it is not -- a hardware trigger that must
    land on its own block, or a gradient you want split at a known instant so a later
    reconstruction step can find the seam.

    Parameters
    ----------
    tag
        Recorded in the compile report so a surprising extra block can be traced back.

    Examples
    --------
    >>> b = barrier('midpoint')
    >>> b.type, b.tag
    ('seqcraft_barrier', 'midpoint')
    >>> span(b)                       # occupies no time
    0.0
    """
    from types import SimpleNamespace as _NS  # noqa: PLC0415

    return _NS(type=BARRIER, tag=str(tag), delay=0.0)
