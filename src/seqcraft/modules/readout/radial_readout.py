r"""
:class:`RadialReadout` -- one radial spoke: prephaser, readout gradient and ADC, at an angle.

``readout/`` because it contains an ADC.  A **leaf**, not a kernel, because it determines every
event it emits from its own parameters: given the field of view, the matrix, the dwell and the
angle, nothing above it is consulted.

Why it is built out of :class:`~seqcraft.modules.CartesianLine`
---------------------------------------------------------------
A 2D radial spoke *is* a Cartesian line pointing somewhere other than along an axis.  The
prephaser that puts the first sample at the right place, the dwell and bandwidth arithmetic, the
ADC that samples the flat top, and -- the hard one -- which sample carries ``k = 0`` are the same
problem, already solved and already tested.  Designing them again here would be a second source
of truth for one piece of arithmetic.

So this module designs a canonical spoke **along x** with ``CartesianLine`` and rotates the
gradients it emits.  What it adds is the rotation, the semantic trajectory geometry a caller
needs, and a contract stated in the vocabulary of a spoke rather than of a line.

``partial_fourier`` spans full spoke to centre-out
---------------------------------------------------
Reusing the existing parameter rather than adding ``readout_asymmetry`` beside it, because it is
the same quantity and it already lands on the right numbers::

    partial_fourier = 1.0    matrix 64 -> 64 samples, centre at 32, k = -32*dk ... +31*dk
    partial_fourier = 0.75   matrix 64 -> 48 samples, centre at 16
    partial_fourier = 0.5    matrix 64 -> 32 samples, centre at  0   -- centre-out

The official UTE example spans the same continuum with an ``readout_asymmetry`` argument running
the other way (``0`` full, ``1`` centre-out) and holds the sample count fixed while shrinking
``dk``; this holds ``dk = 1/FOV`` fixed and drops samples.  Both keep the resolution the protocol
asked for.  The range is deliberately limited to ``0.5 ... 1.0``: below a centre-out spoke there
is no evidence for what the geometry should mean.

**The centre sample is not the middle sample.**  With an even matrix a full spoke is asymmetric
by one -- 64 samples run ``-32*dk ... +31*dk`` -- which is the same convention
:class:`~seqcraft.modules.PhaseEncode` uses for ``center_line = matrix // 2``.  It is not tidied
up here, because consistency between modules is worth more than a symmetric-looking spoke, and
because the official reference does exactly this.

What it does not own
--------------------
Which angles are acquired and in what order -- equal increments, a golden angle, a randomised
schedule -- is acquisition policy, and it stays with the caller.  OpenMRF's radial readout makes
the split visible: its ``phi_mode`` builds the list of angles before any waveform exists.  This
module answers one question, *give me a spoke at this angle*, and asks none.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ...design.events import GRADIENT_KINDS, derive
from ...design.logic import LogicBlock, flatten
from ...design.module import Module
from ...errors import ConfigurationError, format_error
from .cartesian_line import CartesianLine

if TYPE_CHECKING:
    from pypulseq.opts import Opts

    from ...design.events import Event

__all__ = ['RadialReadout']

#: Below this, a direction cosine is treated as exactly zero and its gradient is not emitted, so
#: that a spoke along an axis is one gradient rather than one gradient and a numerical ghost of
#: 1e-11 Hz/m beside it.  Far below any amplitude a scanner can resolve.
_ON_AXIS = 1e-12


class RadialReadout(Module):
    """
    One radial spoke: a prephaser, a readout gradient and an ADC, oriented in the x-y plane.

    Parameters
    ----------
    opts
        The scanner.
    fov_mm
        Field of view across the spoke, millimetres.  Sets ``dk = 1/FOV``.
    matrix
        Samples across the **full** diameter.  With `partial_fourier` below 1 the spoke carries
        fewer than this, and the resolution is unchanged; see :attr:`num_samples`.
    bandwidth_hz_px, dwell_s
        The receiver, given either way.  Exactly one, as for
        :class:`~seqcraft.modules.CartesianLine`.
    partial_fourier
        ``1.0`` is a full spoke through the centre of k-space and ``0.5`` is centre-out, with
        everything between an asymmetric spoke.  Outside that range raises.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    line : CartesianLine
        The canonical spoke, designed along x.  Exposed because a caller reading
        ``radial.line.prephaser_area_per_m`` should not have to rebuild one.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=28, grad_unit='mT/m', max_slew=120, slew_unit='T/m/s',
    ...          rf_dead_time=100e-6, rf_ringdown_time=20e-6, adc_dead_time=10e-6)
    >>> spoke = RadialReadout(opts=o, fov_mm=260.0, matrix=64, dwell_s=20e-6)
    >>> spoke.num_samples, spoke.center_sample
    (64, 32)
    >>> round(spoke.dk_per_m, 4), round(spoke.k_first_per_m, 2), round(spoke.k_last_per_m, 2)
    (3.8462, -123.08, 119.23)
    >>> centre_out = RadialReadout(opts=o, fov_mm=260.0, matrix=64, dwell_s=20e-6,
    ...                            partial_fourier=0.5)
    >>> centre_out.num_samples, centre_out.center_sample
    (32, 0)
    """

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: float,
        matrix: int,
        bandwidth_hz_px: float | None = None,
        dwell_s: float | None = None,
        partial_fourier: float = 1.0,
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.partial_fourier = self._check_partial_fourier(partial_fourier)
        self.line = CartesianLine(
            opts=opts, fov_mm=fov_mm, matrix=matrix, bandwidth_hz_px=bandwidth_hz_px,
            dwell_s=dwell_s, partial_fourier=self.partial_fourier, prephase=True, axis='x',
        )

    # ------------------------------------------------------------------ what it knows
    @property
    def fov_mm(self) -> float:
        """Field of view across the spoke, millimetres."""
        return self.line.fov_mm

    @property
    def matrix(self) -> int:
        """Samples across the full diameter, before `partial_fourier`."""
        return self.line.matrix

    @property
    def num_samples(self) -> int:
        """ADC samples this spoke actually acquires."""
        return self.line.num_samples

    @property
    def dwell_s(self) -> float:
        """Seconds between samples."""
        return self.line.dwell_s

    @property
    def bandwidth_hz_px(self) -> float:
        """Receiver bandwidth per pixel, Hz."""
        return self.line.bandwidth_hz_px

    @property
    def dk_per_m(self) -> float:
        """Sample spacing along the spoke, 1/m.  ``1/FOV``, unchanged by `partial_fourier`."""
        return self.line.dk_per_m

    @property
    def center_sample(self) -> int:
        """
        Which ADC sample carries ``k = 0``.

        **Not the middle of the ADC window**, and the difference is the point.  A full spoke with
        an even matrix has its centre one sample past the middle, and a centre-out spoke has it
        at sample zero.  A caller that assumes ``num_samples // 2`` is right only for the full
        spoke, and wrong silently: a gridding reconstruction still produces an image.

        OpenMRF's radial readout has to compile a probe sequence and read its trajectory back to
        discover this number.  That is the reverse-engineering this property exists to make
        unnecessary.
        """
        return self.line.echo_sample(0)

    def time_to_center(self) -> float:
        """
        Seconds from the start of this module's block to the ``k = 0`` sample.

        Not measurable from the tree: a block knows when its events play, not which instant
        among them is the centre of k-space.
        """
        return self.line.time_to_echo()

    @property
    def k_first_per_m(self) -> float:
        """Signed position of the first sample along the spoke, 1/m.  Negative for a full spoke."""
        return -self.center_sample * self.dk_per_m

    @property
    def k_last_per_m(self) -> float:
        """Signed position of the last sample along the spoke, 1/m."""
        return (self.num_samples - 1 - self.center_sample) * self.dk_per_m

    @property
    def k_max_per_m(self) -> float:
        """The largest radius the spoke reaches, 1/m."""
        return max(abs(self.k_first_per_m), abs(self.k_last_per_m))

    @property
    def prephaser_duration_s(self) -> float:
        """
        Seconds the prephaser occupies.

        Exposed for the reason :class:`~seqcraft.modules.CartesianLine`'s is: a caller overlapping
        it with a slice rephaser has to know how long it is, and measuring it off the block is
        the thing a module should spare them.
        """
        return self.line.prephaser_duration_s

    def direction(self, angle_rad: float) -> tuple[float, float]:
        """The unit vector this spoke points along, as ``(x, y)`` direction cosines."""
        return (math.cos(float(angle_rad)), math.sin(float(angle_rad)))

    # ------------------------------------------------------------------------- assembly
    def build(
        self, *, angle_rad: float, acquire: bool = True, phase_deg: float = 0.0,
    ) -> LogicBlock:
        """
        Return the spoke at `angle_rad`: prephaser, readout gradient and, unless a dummy, the ADC.

        The block this returns is **already oriented**.  A caller does not rotate it afterwards,
        because the module's semantic properties describe the oriented spoke and a trajectory
        whose orientation lived somewhere else would have two owners.

        Parameters
        ----------
        angle_rad
            In-plane angle from the x axis, radians.  Which angles a scan acquires, and in what
            order, is the caller's -- this module places exactly the one it is given.
        acquire
            ``False`` plays identical gradients with no ADC, for a dummy repetition.
        phase_deg
            Receiver phase, degrees.  **Give it the same value the excitation got**: the receiver
            is phase-locked to the transmitter, so an RF-spoiling schedule has to move both.
        """
        cos_a, sin_a = self.direction(angle_rad)
        out = LogicBlock()
        for start, event, _ in flatten(self.line(acquire=acquire, phase_deg=phase_deg)):
            if getattr(event, 'type', '') not in GRADIENT_KINDS:
                out.add(start, event)
                continue
            for channel, scale in (('x', cos_a), ('y', sin_a)):
                if abs(scale) > _ON_AXIS:
                    out.add(start, self._scaled(event, channel, scale))
        return out

    @staticmethod
    def _check_partial_fourier(partial_fourier: float) -> float:
        """
        Return `partial_fourier`, refusing a spoke this module has no evidence for.

        **Closed at both ends**, unlike :func:`~seqcraft.modules._support.require_range`, which is
        half-open because the fractions it bounds are degenerate at their lower end.  0.5 is not
        degenerate here: it is a centre-out spoke, which is a sequence family rather than an edge
        case.  Below it the readout would start on the far side of the centre of k-space and run
        away from it, which is not a spoke, and no reference this was extracted from does it.
        """
        value = float(partial_fourier)
        if not math.isfinite(value) or not 0.5 <= value <= 1.0:
            msg = format_error(
                f'partial_fourier must be in [0.5, 1.0], got {value:g}.',
                {'partial_fourier': value},
                ['1.0 is a full spoke through the centre of k-space',
                 '0.5 is centre-out, the shortest spoke that still reaches it',
                 'between them the spoke is asymmetric about the centre'],
            )
            raise ConfigurationError(msg)
        return value

    def _scaled(self, event: Event, channel: str, scale: float) -> Event:
        """
        One component of a rotated gradient: the same waveform in time, scaled in amplitude.

        Every derived field is scaled with it.  Both radial references that extend a readout's
        flat time warn in a comment that the event's stored area is then wrong; leaving a stale
        ``area`` here would be the same defect, and it would be read by anything asking this
        module what it plays.
        """
        kind = getattr(event, 'type', '')
        if kind == 'trap':
            return derive(
                event, channel=channel,
                amplitude=float(event.amplitude) * scale,
                area=float(event.area) * scale,
                flat_area=float(event.flat_area) * scale,
            )
        if kind == 'grad':
            return derive(
                event, channel=channel,
                waveform=np.asarray(event.waveform, dtype=float) * scale,
                first=float(event.first) * scale,
                last=float(event.last) * scale,
            )
        msg = format_error(                                              # pragma: no cover
            f'cannot rotate a {kind!r} gradient event.',
            {'type': kind, 'channel': getattr(event, 'channel', None)},
            ['a radial spoke is built from trapezoids and arbitrary gradients'],
        )
        raise ConfigurationError(msg)                                    # pragma: no cover
