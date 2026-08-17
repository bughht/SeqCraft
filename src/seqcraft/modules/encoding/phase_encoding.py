"""
:class:`PhaseEncode` -- one Cartesian phase-encode blip, designed once and scaled per line.

The module exists for one invariant: **every line takes the same time.**  Designing at the
largest area and scaling down with ``pp.scale_grad`` is what buys it, and it is what lets the
caller place the readout without knowing which line is being acquired.  Designing each blip at
its own minimum duration would be shorter on average and would move the echo line by line.

``encoding/`` because it is a gradient with no ADC, imposing a phase you intend to sample later.
The axis is an argument rather than a fact about the class: the same module encodes ``y`` in a 2D
scan and ``z`` as the partition axis when 3D arrives.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pypulseq as pp

from ...design.events import derive
from ...design.logic import LogicBlock
from ...design.module import Module
from ...design.timing import EPS
from ...errors import ConfigurationError, format_error
from .._support import ceil_raster, require_axis, require_positive

if TYPE_CHECKING:
    from pypulseq.opts import Opts

__all__ = ['PhaseEncode']

#: A blip needs two lines to have a direction; one is a scan with no phase encoding in it.
_MIN_MATRIX = 2


class PhaseEncode(Module):
    """
    A phase-encode blip for line `line` of `matrix`, on one axis.

    Parameters
    ----------
    opts
        The scanner.
    fov_mm
        Field of view along `axis`.  Sets the k-space step, ``dk = 1000 / fov_mm`` in 1/m.
    matrix
        Number of encoded lines.
    axis
        Logical gradient channel.
    duration_s
        Stretch every blip to this length.  ``None`` is the shortest legal one, which is
        :attr:`min_duration_s`.  A composite passes the winder maximum here; see
        :class:`~seqcraft.modules.GRE2DTR`.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s')
    >>> pe = PhaseEncode(opts=o, fov_mm=250.0, matrix=64)
    >>> pe.center_line
    32
    >>> pe(line=32).nodes[0].item.area                      # the centre line encodes nothing
    0.0
    >>> round(float(pe(line=40).nodes[0].item.area), 6)     # (40 - 32) * 1000 / 250
    32.0

    Every line takes the same time, which is the whole point:

    >>> {round(pe(line=n).duration * 1e6) for n in (0, 32, 63)}
    {300}

    A rewinder is the same blip with the opposite area and the same duration:

    >>> back = pe(line=40, rewind=True)
    >>> round(float(back.nodes[0].item.area), 6), round(back.duration * 1e6)
    (-32.0, 300)

    Notes
    -----
    **`line` runs ``0 ... matrix - 1``**, with no default, and out-of-range values raise::

        k_per_m = (line - center_line) * (1000 / fov_mm)        # center_line == matrix // 2

    which spans −32 … +32 for a matrix of 65 and −32 … +31 for 64, with no special case for
    parity.  Zero-based is what lines up with array indices, with pulseq's ``LIN`` label, and
    with ``range(matrix)``, so the same integer means the same line in all four places.

    `rewind` is a *build* argument rather than a ``rewinder()`` method because **only ``build``
    may return a block**: ``__call__`` is the single place :meth:`Module._finalize` type-checks
    and tags what came back, and a public method returning a :class:`LogicBlock` would be an
    unchecked second exit.
    """

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: float,
        matrix: int,
        axis: str = 'y',
        duration_s: float | None = None,
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.fov_mm = require_positive(fov_mm, 'fov_mm')
        self.matrix = int(matrix)
        self.axis = require_axis(axis)
        if self.matrix < _MIN_MATRIX:
            msg = format_error(
                f'matrix must be at least {_MIN_MATRIX}, got {self.matrix}.',
                {'matrix': self.matrix},
                ['a single-line scan needs no phase encoding at all'],
            )
            raise ConfigurationError(msg)

        #: k-space step between adjacent lines, 1/m.
        self.dk_per_m = 1e3 / self.fov_mm
        largest = pp.make_trapezoid(
            channel=self.axis, area=self.center_line * self.dk_per_m, system=opts,
        )
        self._min_duration_s = float(pp.calc_duration(largest))

        if duration_s is None:
            self._blip = largest
        else:
            wanted = ceil_raster(
                require_positive(duration_s, 'duration_s'), opts.grad_raster_time,
            )
            # EPS, not exact: a caller passing back the maximum of this minimum and another
            # module's can land one ulp below it after snapping onto the raster.
            if wanted < self._min_duration_s - EPS:
                msg = format_error(
                    f'duration_s = {duration_s * 1e6:.1f} us is shorter than the largest blip '
                    f'this matrix needs.',
                    {
                        'duration_s': duration_s,
                        'min_duration_s': self._min_duration_s,
                        'matrix': self.matrix,
                        'fov_mm': self.fov_mm,
                    },
                    [
                        f'pass duration_s >= {self._min_duration_s:.6g}',
                        'or pass duration_s=None for the shortest legal blip',
                    ],
                )
                raise ConfigurationError(msg)
            self._blip = pp.make_trapezoid(
                channel=self.axis,
                area=self.center_line * self.dk_per_m,
                duration=wanted,
                system=opts,
            )

    # ------------------------------------------------------------------ what it knows
    @property
    def center_line(self) -> int:
        """
        The line that encodes ``k = 0``, as ``matrix // 2``.

        Defined once and read, never recomputed: three consumers must agree on it -- the ``k``
        mapping here, the ``LIN`` label a composite emits, and the ``kSpaceCenterLine`` written
        into the ``.seq`` -- and two of them are in other files.
        """
        return self.matrix // 2

    @property
    def min_duration_s(self) -> float:
        """Shortest legal blip at the largest area, seconds.  Unchanged by `duration_s`."""
        return self._min_duration_s

    def k_per_m(self, line: int) -> float:
        """Return the k-space position line `line` encodes, in 1/m."""
        return (self._check(line) - self.center_line) * self.dk_per_m

    # ----------------------------------------------------------------------- assembly
    def build(self, *, line: int, rewind: bool = False) -> LogicBlock:
        """
        Return the blip for `line`, or its rewinder.

        Parameters
        ----------
        line
            Zero-based phase-encode index, ``0 ... matrix - 1``.
        rewind
            Negate the area: same line, same duration, opposite sign, for the end of the TR.
        """
        area = self.k_per_m(line)
        scale = (-area if rewind else area) / float(self._blip.area)
        # derive() strips pypulseq's registration state from the scaled copy, so the same design
        # can be scaled again after a compile has registered one of its outputs.
        return LogicBlock().add(0.0, derive(pp.scale_grad(self._blip, scale)))

    def _check(self, line: int) -> int:
        """Return `line` as an int, having checked it addresses a line of this matrix."""
        index = int(line)
        if not 0 <= index < self.matrix:
            msg = format_error(
                f'line = {index} is outside 0 ... {self.matrix - 1}.',
                {'line': index, 'matrix': self.matrix, 'center_line': self.center_line},
                [
                    f'lines are zero-based, so the centre of k-space is {self.center_line}',
                    f'the full table is range({self.matrix})',
                ],
            )
            raise ConfigurationError(msg)
        return index
