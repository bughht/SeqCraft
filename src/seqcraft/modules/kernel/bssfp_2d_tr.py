"""
:class:`bSSFP2DTR` -- one balanced repetition of a 2D Cartesian bSSFP acquisition.

``kernel/`` for the same reason :class:`~seqcraft.modules.GRE2DTR` is: it composes leaves from
more than one folder, and it is the repeating unit.

What makes it balanced
----------------------
**The net gradient area on each axis is zero over the RF-to-RF interval.**  That is the
definition.  ``TE = TR/2`` is the standard symmetric realisation of it and the default here,
not the definition -- an asymmetric echo is still balanced if the net area is still zero, and
this module will build one when asked.

The timing origin is the **RF effective centre**
----------------------------------------------
Every time this module declares is measured from there::

    TE               = echo - RF effective centre of repetition n
    TR               = RF effective centre of n+1 - RF effective centre of n
    balance interval = [RF effective centre of n, RF effective centre of n+1]

The effective centre is read from the RF waveform through ``pp.calc_rf_center``, exactly as
:meth:`Excitation.time_to_center` does.  **It is not the waveform's midpoint**, not the event's
start, and not the block's start or end: for a minimum-phase or otherwise asymmetric pulse the
effective centre is nowhere near halfway through the event, and a caller that assumed otherwise
would be declaring a TE it does not achieve.  :meth:`time_to_rf_center` reports where it lands
inside the returned block, so nothing downstream has to reconstruct it from a pulse duration.

The interval is **RF centre to RF centre**, which is not the interval the echo lives in and not
this block's edges either.  For repetition ``n`` it starts at RF centre ``n``, runs through
everything this repetition plays, and ends at RF centre ``n+1`` -- so it contains the
*post*-centre half of this selection gradient and the *pre*-centre half of the **next** one.

Which is why the z axis is settled by four terms contributed by two repetitions::

    A_post(n) + R_tail(n)  +  W_lead(n+1) + A_pre(n+1)  =  0
    \\____________________/     \\____________________/
       played by n                 played by n+1

``A_post`` and ``A_pre`` are the selection gradient's areas after and before the RF centre;
``R_tail`` and ``W_lead`` are the two balancing lobes this module emits.  For a symmetric pulse
``A_post == A_pre == A`` and the realisation is ``R_tail = W_lead = -A``.

**Each balancing lobe cancels one selection half, not two.**  A trailing lobe of ``-(A_post +
A_pre)`` beside a non-zero leading winder counts the neighbour's half twice, and the interval
then carries a net ``-A``.  The extra gradient a balanced repetition carries over a spoiled one
is a second lobe, not a larger one.

The contract is compositional
-----------------------------
No single emitted repetition is self-contained, because the interval ends inside its successor.
What this module guarantees is:

    two consecutive repetitions, each discharging **its own** selection halves, have zero net
    gradient area on each axis over the RF-centre-to-RF-centre interval between them.

That is a weaker condition on the neighbour than it first looks, and the split above is why.
Because ``R_tail`` cancels only ``A_post`` and ``W_lead`` only ``A_pre``, the interval reads::

    A_post(n) - A_post(n)  +  -A_pre(n+1) + A_pre(n+1)   =   0

-- two independent cancellations rather than one shared sum.  Neither repetition is built
against an assumption about the other, so **the geometry does not have to match**: a neighbour
with a different slice thickness, a different pulse duration or a different readout bandwidth
still composes to a balanced interval, as do a different `flip_deg` and a different `line`.
Each of those is measured in the tests rather than argued here.

What the neighbour does have to do is honour the same contract -- discharge its own pre-centre
half before its own RF centre.  A neighbour that plays no leading winder leaves ``A_pre`` in the
interval uncancelled, and a spoiled gradient-echo repetition is exactly such a neighbour.  So
the condition is on the *contract* the neighbour keeps, not on the protocol it was configured
with.

Balance-compatible is not timing-compatible
-------------------------------------------
Those are two properties and they do not come together::

    balance-compatible   neighbouring repetitions satisfy the RF-to-RF M0 contract
    timing-compatible    they can additionally be stacked by block duration while preserving
                         the declared RF-centre-to-RF-centre TR

Stacking by block duration puts repetition ``n+1``'s block start at ``start_n + tr_s``, so the
interval the magnetisation actually sees is::

    block_duration(n) + time_to_rf_center(n+1) - time_to_rf_center(n)

which equals ``tr_s`` only when the two RF-centre offsets agree.  Measured on compiled trains, a
neighbour differing in `flip_deg`, `line` or -- on the protocols tested -- `bandwidth_hz_px`
leaves the offset alone and is both; a neighbour differing in `thickness_mm` or `rf_duration_s`
is balance-compatible and **not** timing-compatible, missing the declared TR by 120 and 500
microseconds respectively.  `bandwidth_hz_px` reaches the shared window and so could move the
offset on another protocol; the reliable test is :meth:`time_to_rf_center` itself, not the
parameter list.

A heterogeneous train is placed from the offsets rather than from the durations::

    start(n+1) = start(n) + tr_s + time_to_rf_center(n) - time_to_rf_center(n+1)

This module does not do that placement for you, and nothing here is an acquisition framework.
The shipped examples stack by duration, because a continuous acquisition and a start-up
flip-angle ramp are both timing-compatible.

**The lumped alternative is where matching would have mattered.**  A trailing lobe of
``-(A_post + A_pre)`` with no leading winder is also balanced, and for a train of identical
repetitions it is indistinguishable -- but it is built against an assumed successor, so it
breaks the moment the successor's pulse differs.  The split realised here is the one that does
not need the assumption.

The symmetric realisation, and why the lobes sit where they do
--------------------------------------------------------------
Zero net area is the balance condition, and it does not by itself say where in the interval the
area was played.  The standard symmetric 2D scheme says more than that: the readout and slice
structures are arranged symmetrically about the echo, which makes them **first-order** balanced
over the interval as well -- the structure the bSSFP flow literature assumes when it says that
only the phase-encode axis varies from repetition to repetition.

That symmetry has to be built, not hoped for.  On ``z`` the four terms of the interval are::

    [+select post][gap_post][-rephase]  ...  [-winder][gap_pre][+select pre]

and they mirror about the midpoint only when ``gap_post + post_half == gap_pre + pre_half``.  The
pulse's halves are **not** equal -- the transmit dead time sits inside the selection lobe but
before the RF -- so this module pads the shorter side.  With that pad the first moment over the
interval is zero on ``z`` to machine precision, and without it the trailing lobe drifts wherever
the echo time puts it.

What each axis then carries over the RF-to-RF interval::

    z   M1 = 0            the two balancing lobes mirror exactly
    x   M1 small, and the SAME for every line -- the prephaser cancels the area before the echo
        and the balancing lobe the area after it, and the half-dwell offset of the echo sample
        makes those differ a little, so the two lobes mirror in time but not quite in area
    y   M1 varies with the line, because the blip and its rewind have opposite signs

That last row is the one a steady state is sensitive to.  A constant per-TR phase is harmless --
the condition is that the phase be the *same* from one repetition to the next -- so what
perturbs a balanced steady state is the phase-encode term, and how much depends on how far apart
in k-space two consecutively acquired lines are.  That is view ordering, and it lives in the
acquisition.  ``examples/bssfp_2d/03_flow_and_motion.ipynb`` measures all of this against Bieri
and Scheffler's criterion.

What this layer does not own
----------------------------
Steady-state establishment of any kind, and every decision attached to it: the start-up or
catalyzation method, how many start-up repetitions there are, segmentation and the segment
count, ``ky`` and acquisition ordering, interruption and restart policy, and reconstruction.

**A balanced repetition is not a repetition in steady state.**  The waveform condition above is
true of the first interval of a two-repetition sequence; the magnetisation reaches steady state
only after a train, and only for a given tissue.  This module can claim the first and never the
second.  Segmentation is a loop around the repetitions and gets no class -- see the shipped
example, where it is ordinary authoring code.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pypulseq as pp

from ...augmentation import FlowCompensation
from ...design import _augment
from ...design import joint as _joint
from ...design import scope as _scope
from ...design.events import Event, derive
from ...design.logic import LogicBlock
from ...design.module import Module
from ...errors import ConfigurationError, format_error
from .._support import area_until, ceil_raster, require_axis, require_pair, require_positive
from ..encoding.phase_encoding import PhaseEncode
from ..readout.cartesian_line import CartesianLine
from ..rf.excitation import Excitation

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from pypulseq.opts import Opts

__all__ = ['bSSFP2DTR']


class bSSFP2DTR(Module):  # noqa: N801 -- `bSSFP` is the domain spelling; see ruff.toml
    """
    One balanced repetition of a 2D Cartesian bSSFP acquisition.

    Parameters
    ----------
    opts
        The scanner.
    fov_mm
        Field of view, millimetres, as ``(readout, phase)`` or one value for both.
    matrix
        ``(nx, ny)``.
    thickness_mm
        Slice thickness, millimetres.  Selective only: a non-selective excitation is a 3D case,
        and 3D is a different public class.
    flip_deg
        Flip angle, degrees.  May differ between adjacent repetitions -- see **Compatible** in
        the module documentation.
    rf_duration_s
        Forwarded to :class:`~seqcraft.modules.Excitation`.  ``None`` leaves that module's own
        default, rather than copying a second one here that could disagree with it.
    te_s
        Echo time, seconds, from the RF effective centre.

        ``None`` does **not** mean "the shortest echo this protocol can reach", which is what it
        means on :class:`~seqcraft.modules.GRE2DTR`.  Here it means ``TR/2``: the symmetric timing
        the standard description of the sequence assumes.  Passing a value asks for an asymmetric-TE realisation, which is
        still balanced and still legal, but is a choice with a contrast consequence rather than
        a default anybody should arrive at by passing nothing.
    tr_s
        Repetition time, seconds, measured RF centre to RF centre.  With ``te_s=None`` this is
        the whole specification: TE follows as ``TR/2``.  ``None`` for both asks for the
        shortest legal **symmetric** repetition.
    bandwidth_hz_px
        Readout bandwidth per pixel.
    partial_fourier
        Readout-direction partial echo, forwarded to
        :class:`~seqcraft.modules.CartesianLine`.  This is one of the cases where an asymmetric
        `te_s` is the point; phase-encode-direction partial Fourier is a different thing and
        does not call for one.
    flow_comp
        :class:`~seqcraft.FlowCompensation`, asking that the first gradient moment be zero
        **at the echo** on the named axes as well as the zeroth.

        **This is echo-time first-moment compensation, and it is not the flow compensation the
        bSSFP literature describes.**  Bieri and Scheffler (Magn Reson Med 2005;54:901) null the
        first moment *between excitations*, and say explicitly that their design "does not null
        the first moment between the excitation and echo", because what perturbs a steady state
        is a phase that varies from one repetition to the next rather than a phase on one
        acquired line.  Asking for this here nulls the echo-time moment and leaves the
        RF-to-RF moment alone -- on the shipped example protocol it makes the
        repetition-to-repetition increment larger, which
        ``examples/bssfp_2d/03_flow_and_motion.ipynb`` measures.

        What it is good for is the echo-time phase itself: a spin moving at constant velocity
        arrives at the echo carrying ``2*pi*v*M1``, and this removes that term from each acquired
        line.  What it is not is a drop-in answer to flow artefacts in a bSSFP acquisition.

        Balance is untouched either way.  The designer reshapes what plays between the RF centre
        and the echo and is handed the same **total** area on each axis, so the RF-to-RF sum is
        the one the conventional repetition had -- which the tests measure rather than assume.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Examples
    --------
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
    ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> tr = bSSFP2DTR(opts=o, fov_mm=250.0, matrix=(128, 128), thickness_mm=5.0,
    ...                flip_deg=35.0, bandwidth_hz_px=800.0)

    The standard symmetric realisation, as closely as the raster permits, and by how much it
    misses -- see :attr:`symmetry_residual_s`:

    >>> abs(tr.te_s - tr.tr_s / 2) <= o.grad_raster_time / 2
    True
    >>> round(tr.symmetry_residual_s * 1e6, 3)
    2.1

    The block lasts exactly one repetition time, which is what lets a caller stack them:

    >>> tr(line=64).duration == tr.tr_s
    True

    Notes
    -----
    **The receiver is phase-locked to the transmitter.**  `phase_deg` reaches the RF and the ADC
    from the same argument in the same place, because a phase progression applied to one and not
    the other writes itself into ``ky``.  The progression itself -- the alternation bSSFP is
    usually run with, or any other schedule -- belongs to whatever knows how many repetitions
    there are, exactly as RF spoiling does on :class:`~seqcraft.modules.GRE2DTR`.

    **There are no spoilers.**  A spoiler is a deliberate non-zero net area, which is the
    condition this module exists to not have.
    """

    _PE_AXIS = 'y'

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: float | tuple[float, float],
        matrix: tuple[int, int],
        thickness_mm: float,
        flip_deg: float = 35.0,
        rf_duration_s: float | None = None,
        te_s: float | None = None,
        tr_s: float | None = None,
        bandwidth_hz_px: float = 400.0,
        partial_fourier: float = 1.0,
        flow_comp: FlowCompensation | None = None,
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        _augment.require_intent_type(flow_comp, FlowCompensation, 'flow_comp')
        fov_x, fov_y = require_pair(fov_mm, 'fov_mm')
        nx, ny = require_pair(matrix, 'matrix')
        self.fov_mm = (fov_x, fov_y)
        self.matrix = (int(nx), int(ny))
        self.thickness_mm = require_positive(thickness_mm, 'thickness_mm')

        pulse = {} if rf_duration_s is None else {'duration_s': rf_duration_s}
        self.exc = Excitation(opts=opts, flip_deg=flip_deg, thickness_mm=self.thickness_mm,
                              **pulse)

        # The two selection halves, measured on the gradient rather than assumed equal.  A
        # minimum-phase pulse's effective centre is nowhere near its midpoint, so `A_pre == A_post`
        # is a property of the pulse and not something to build in.
        self._a_pre = area_until(self.exc.gz, self.exc.time_to_center())
        self._a_post = -self.exc.rephaser_area_per_m

        # The earliest anything may follow the pulse, in the excitation block's own frame: the
        # transmit chain has to be clear of it and the selection gradient has to be done.
        raster = float(opts.grad_raster_time)
        self._grad_start_s = ceil_raster(
            max(float(pp.calc_duration(self.exc.gz)), float(pp.calc_duration(self.exc.rf))),
            raster,
        )

        # The two slice balancing lobes, placed so that the z structure is **symmetric about the
        # echo**.  That symmetry is what makes the slice axis first-order compensated over the
        # RF-to-RF interval, which is the structure the balanced-SSFP literature describes for a
        # symmetric 2D scheme -- and zero net area alone does not give it.
        #
        # Each lobe has to sit the same distance from its own RF centre:
        #
        #     ... [+select post][gap_post][-rephase] ... [-winder][gap_pre][+select pre] ...
        #
        # mirrors about the midpoint only when `gap_post + post_half == gap_pre + pre_half`.  The
        # pulse's own halves are not equal -- the transmit dead time sits inside the selection
        # lobe but before the RF -- so the shorter side is padded to match the longer one.
        self._z_window_s = ceil_raster(
            max(self._lobe_duration_s(-self._a_pre), self._lobe_duration_s(-self._a_post)), raster)
        post_half = self._grad_start_s - self.exc.time_to_center()
        pre_half = self.exc.time_to_center()
        #: Gap between the selection gradient and the trailing lobe, and between the leading lobe
        #: and the next selection gradient.  One of the two is always zero.
        self._z_gap_s = ceil_raster(max(pre_half - post_half, 0.0), raster)
        self._lead_gap_s = ceil_raster(max(post_half - pre_half, 0.0), raster)
        self._lead_s = self._z_window_s + self._lead_gap_s
        self._z_lead = self._lobe('z', -self._a_pre, self._z_window_s)

        # Who owns a moment requirement on each axis.  All three are this repetition's: `x` is
        # the readout's own prephaser, `y` the encode blip, `z` the trailing balancing lobe.
        self._moment_owners: dict[str, str] = {'x': 'readout', 'y': 'joint', 'z': 'joint'}
        routed = self._route(flow_comp)
        readout_moment = {'_null_moment_order': 1} if routed.get('x') else {}

        probe_ro = CartesianLine(opts=opts, fov_mm=fov_x, matrix=self.matrix[0], axis='x',
                                 bandwidth_hz_px=bandwidth_hz_px,
                                 partial_fourier=partial_fourier, **readout_moment)
        probe_pe = PhaseEncode(opts=opts, fov_mm=fov_y, matrix=self.matrix[1], axis='y')
        # One window serves every axis between the pulse and the readout, and one more after it:
        # axes overlap for free, so each window is its longest participant rather than their sum.
        self.winder_s = ceil_raster(
            max(probe_ro.prephaser_duration_s, probe_pe.min_duration_s,
                self._lobe_duration_s(-probe_ro.area_after_echo_per_m)),
            raster,
        )
        # The jointly designed winder, when something claimed `y` or `z`.  It reshapes what plays
        # between the RF centre and the echo and leaves the **total** area on each axis alone,
        # which is why the balance condition is unaffected by it.
        self._joint: dict[str, object] = {}
        #: The delay the designed schedule puts between the winder and the readout, when an
        #: explicit echo time needs one.  Zero unless the scope is in use.
        self._joint_fill_s = 0.0
        jointly = {a for a, owner in routed.items() if owner == 'joint'}
        claims, self._encoding_states = _augment.claims_and_states(
            _augment.flow_comp_for(flow_comp, jointly), None)
        if claims:
            self._solve_with_scope(claims, probe_ro, fov_y, bandwidth_hz_px, partial_fourier,
                                   readout_moment, opts)

        self.ro = CartesianLine(opts=opts, fov_mm=fov_x, matrix=self.matrix[0], axis='x',
                                bandwidth_hz_px=bandwidth_hz_px, partial_fourier=partial_fourier,
                                prephaser_duration_s=self.winder_s, **readout_moment)
        self.pe = PhaseEncode(opts=opts, fov_mm=fov_y, matrix=self.matrix[1], axis='y',
                              duration_s=self.winder_s)

        #: ``None`` when the scope designed `z`: the waveform it produced carries the same total
        #: area, and is emitted in its place -- and in the scope's own window rather than in the
        #: symmetric slot, because a jointly designed waveform is solved against that schedule.
        self._z_tail = (None if 'z' in self._joint
                        else self._lobe('z', -self._a_post, self._z_window_s))
        self._x_balance = self._lobe('x', -self.ro.area_after_echo_per_m, self.winder_s)

        self._ro_duration_s = self.ro().duration
        self._tail_s = self.winder_s
        self._te_s, self._tr_s = self._resolve_timing(te_s, tr_s)
        if self._joint and self._head_fill_s > 1e-12:          # pragma: no cover - guarded above
            raise ConfigurationError(format_error(
                'the designed winder would be emitted away from the echo it was designed for.',
                {'head_fill_s': self._head_fill_s, 'joint_fill_s': self._joint_fill_s},
                ['this is an internal invariant; please report it with the protocol'],
            ))
    # ------------------------------------------------------------------ what it knows
    @property
    def center_line(self) -> int:
        """The phase-encode line that encodes ``k = 0``.  One convention, defined once."""
        return self.pe.center_line

    def voxel_mm(self, axis: str) -> float:
        """The voxel dimension along `axis`, millimetres."""
        return {
            'x': self.fov_mm[0] / self.matrix[0],
            'y': self.fov_mm[1] / self.matrix[1],
            'z': self.thickness_mm,
        }[require_axis(axis)]

    @property
    def min_te_s(self) -> float:
        """
        Shortest echo time this repetition can reach, seconds, **as an asymmetric realisation**.

        This is the earliest the echo can physically arrive.  It is not the default, and on a
        symmetric protocol it is shorter than :attr:`te_s`: see :attr:`min_symmetric_te_s`.
        """
        return self._head_min_s

    @property
    def min_symmetric_te_s(self) -> float:
        """
        Shortest echo time reachable **with** ``TE = TR/2``, seconds.

        The larger of the two halves: an echo earlier than this cannot be the midpoint of any
        repetition, because the second half would have to be shorter than the readout tail, the
        rewinders and the next excitation's leading winder together.  This is the TE a default
        ``bSSFP2DTR`` achieves, and half of :attr:`min_tr_s`, both to within
        :attr:`symmetry_residual_s`.
        """
        head, _ = self._symmetric_fills()
        return self._head_min_s + head

    @property
    def symmetry_residual_s(self) -> float:
        """
        How far this repetition's echo sits from the exact midpoint, seconds.  Signed.

        ``TE = TR/2`` is the symmetric-timing target, and the default aims at it.  **For a given
        readout geometry the exact midpoint may not lie on the available timing lattice**, and
        when it does not, the default takes the nearest symmetric realisation and reports what
        it achieved here rather than rounding the declared TE to ``TR/2``.

        Where the offset comes from, on the protocols measured: the echo is the ADC sample where
        ``k = 0`` and that sample sits half a dwell off the readout window's centre, while both
        fills are whole numbers of the gradient raster -- so the difference between the two
        halves moves in raster steps and the quantity it has to cancel need not be one of them.
        That makes this at most half a gradient raster with the readout placement
        :class:`~seqcraft.modules.CartesianLine` currently offers: 5 microseconds at the usual
        10, against a TE of some thousands.  A readout geometry whose echo does land on the
        raster gives exactly zero.

        It is reported rather than hidden, because a quantity that is zero except when it is not
        is exactly the kind a caller should be able to assert on.

        Examples
        --------
        >>> from pypulseq.opts import Opts
        >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
        ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
        >>> tr = bSSFP2DTR(opts=o, fov_mm=250.0, matrix=(128, 128), thickness_mm=5.0,
        ...                flip_deg=35.0, bandwidth_hz_px=800.0)
        >>> abs(tr.symmetry_residual_s) <= o.grad_raster_time / 2
        True
        """
        return self._te_s - self._tr_s / 2.0

    @property
    def te_s(self) -> float:
        """The echo time this repetition achieves, seconds, from the RF effective centre."""
        return self._te_s

    @property
    def min_tr_s(self) -> float:
        """
        Shortest **symmetric** repetition time, seconds.

        Deliberately not "the shortest repetition that can hold this readout": the symmetric
        realisation is the one this module defaults to, so the minimum it reports is the minimum
        of *that* realisation.  An explicit `te_s` asks a different question and relaxes it,
        which is what :attr:`min_te_s` reports against.
        """
        head, tail = self._symmetric_fills()
        return self._span_s(head, tail)

    @property
    def tr_s(self) -> float:
        """
        The repetition time, seconds: **RF effective centre to RF effective centre**.

        That is the definition.  A block of this module also happens to *last* ``tr_s``, which
        is what lets a homogeneous train be stacked at ``n * tr_s`` -- but that is a convenience
        of this composition rather than what TR means, and it stops delivering the declared TR
        the moment two neighbours disagree on :meth:`time_to_rf_center`.
        """
        return self._tr_s

    def time_to_rf_center(self) -> float:
        """
        Seconds from the start of this module's block to the **RF effective centre**.

        The timing origin of everything this module declares: :attr:`te_s` is measured from
        here, :attr:`tr_s` is the distance from here to the next repetition's, and the balance
        interval runs between the same two instants.

        Read from the waveform through :meth:`Excitation.time_to_center`, which is
        ``pp.calc_rf_center`` plus the pulse's own delay.  **Not the midpoint of the RF event**
        -- for a minimum-phase pulse those are far apart, and the difference is the error a
        caller makes by deriving the origin from a pulse duration instead of asking.

        It exists so that nothing downstream has to reconstruct it.  Two repetitions with
        different offsets cannot be stacked by block duration without moving the RF-to-RF
        interval; see **Balance-compatible is not timing-compatible** above for the placement a
        heterogeneous train needs, and the identity below for the one a homogeneous one gets.

        Examples
        --------
        >>> from pypulseq.opts import Opts
        >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
        ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
        >>> tr = bSSFP2DTR(opts=o, fov_mm=250.0, matrix=(128, 128), thickness_mm=5.0,
        ...                flip_deg=35.0, bandwidth_hz_px=800.0)
        >>> abs(tr.time_to_echo() - (tr.time_to_rf_center() + tr.te_s)) < 1e-12
        True
        """
        return self._lead_s + self.exc.time_to_center()

    def time_to_echo(self) -> float:
        """
        Seconds from the start of this module's block to ``k = 0``.

        ``time_to_rf_center() + te_s``, and asserted to be, because TE is defined from the
        effective centre and this is the only other place the same origin is used.
        """
        return self.time_to_rf_center() + self._te_s

    # ----------------------------------------------------------------------- assembly
    def build(self, *, line: int, phase_deg: float = 0.0, acquire: bool = True,
              center_mm: tuple[float, float, float] = (0.0, 0.0, 0.0)) -> LogicBlock:
        """
        Return one balanced repetition.

        Parameters
        ----------
        line
            Zero-based phase-encode index.
        phase_deg
            RF carrier phase for this repetition, degrees -- and the receiver phase, which is
            set from the same value here so the two cannot disagree.
        acquire
            ``False`` drops the ADC and the ``LIN`` label and changes nothing else, so a
            start-up repetition loads the gradients exactly as an acquired one does.  **It does
            not make the repetition a preparation**: what start-up means is the caller's, and
            this argument only says whether to sample.
        center_mm
            ``(x, y, z)`` centre of the imaging volume, millimetres.  ``y`` must be ``0.0``.
        """
        x, y, z = (float(v) for v in center_mm)
        if y != 0.0:
            msg = format_error(
                f'center_mm y = {y:g} mm is out of scope: a phase-encode shift is a per-line '
                f'ADC phase, not a per-TR frequency.',
                {'center_mm': center_mm},
                ['shift the slice with z, or the readout FOV with x'],
            )
            raise ConfigurationError(msg)

        # The winder starts after whatever delay an explicit echo time needs, because that is
        # where the designed schedule put it -- emitting it earlier would translate a solved
        # waveform out from under its own moment.  Zero unless the scope is in use.
        winder = self._lead_s + self._grad_start_s + self._head_fill_s + self._joint_fill_s
        tail = winder + self._ro_duration_s
        out = (
            LogicBlock()
            # The leading winder first, then the excitation: this lobe belongs to the interval
            # that ENDS at the RF centre just after it.
            .add(0.0, self._z_lead)
            .add(self._lead_s, self.exc(phase_deg=phase_deg, position_mm=z, rephase=False))
            .add(winder, self._encode(line))
            .add(winder, self.ro(acquire=acquire, phase_deg=phase_deg, offset_mm=x))
            .add(tail, self.pe(line=line, rewind=True))
            .add(tail, self._x_balance)
        )
        if self._z_tail is not None:
            # Immediately after the selection gradient plus whatever pad the symmetry needs --
            # not in the shared winder window, which sits wherever the echo time puts it.
            out.add(self._lead_s + self._grad_start_s + self._z_gap_s, self._z_tail)
        designed = self._joint_block('z', (0, _augment.ONE_STATE))
        if designed is not None:
            out.add(winder, designed)
        if acquire:
            out.add(winder, pp.make_label(type='SET', label='LIN', value=int(line)))
        fill = self._tr_s - (tail + self._tail_s)
        if fill > 1e-9:
            out.add(tail + self._tail_s, pp.make_delay(fill))
        return out

    # ----------------------------------------------------------------------- design
    def _lobe(self, channel: str, area_per_m: float, duration_s: float) -> Event:
        """One balancing lobe of `area_per_m`, stretched to the shared window."""
        return derive(pp.make_trapezoid(channel=channel, area=area_per_m,
                                        duration=duration_s, system=self.opts))

    def _lobe_duration_s(self, area_per_m: float) -> float:
        """
        The shortest this lobe can be, seconds -- a participant in the shared window.

        The channel does not enter: `opts` holds one amplitude and one slew, so the minimum is a
        property of the area.  Measured rather than derived, so a change in how pypulseq rounds
        a ramp onto the raster reaches the window instead of being approximated here.
        """
        if abs(area_per_m) < 1e-12:
            return 0.0
        return float(pp.calc_duration(
            pp.make_trapezoid(channel='z', area=area_per_m, system=self.opts)))

    @property
    def _z_tail_end_s(self) -> float:
        """Seconds from the RF effective centre to the end of the trailing slice lobe."""
        return (self._grad_start_s - self.exc.time_to_center()
                + self._z_gap_s + self._z_window_s)

    @property
    def _head_floor_s(self) -> float:
        """
        The smallest head fill the geometry allows, seconds.

        The winder window may start while the trailing slice lobe is still playing -- different
        axes overlap for free -- but the lobe has to be **finished by the echo**, or slice
        gradient would still be running when ``k_z`` is supposed to be zero.  Usually zero: the
        readout is longer than the lobe on every protocol shipped here.
        """
        head_without_fill = (self._grad_start_s - self.exc.time_to_center()
                             + self._joint_fill_s + self.ro.time_to_echo())
        return max(self._z_tail_end_s - head_without_fill, 0.0)

    @property
    def _head_min_s(self) -> float:
        """TE at the shortest: RF centre to echo with only the fill the geometry forces."""
        return (self._grad_start_s - self.exc.time_to_center()
                + self._joint_fill_s + self.ro.time_to_echo() + self._head_floor_s)

    @property
    def _tail_min_s(self) -> float:
        """Echo to the NEXT RF centre, at the shortest -- the other half of the interval."""
        return ((self._ro_duration_s - self.ro.time_to_echo()) + self._tail_s
                + self._lead_s + self.exc.time_to_center())

    @property
    def _raster_s(self) -> float:
        """The granularity every fill moves in: both rasters have to be satisfied at once."""
        return max(float(self.opts.grad_raster_time), float(self.opts.block_duration_raster))

    def _span_s(self, head_fill_s: float, tail_fill_s: float) -> float:
        """The RF-to-RF duration these two fills produce, seconds -- also the block duration."""
        return (self._lead_s + self._grad_start_s + head_fill_s + self._joint_fill_s
                + self._ro_duration_s + self._tail_s + tail_fill_s)

    def _symmetric_fills(self) -> tuple[float, float]:  # noqa: D401
        """
        The shortest pair of fills that puts the echo at the midpoint, as closely as it can go.

        The two halves differ by ``(head - tail) - delta``, where ``delta`` is what the fixed
        parts of the repetition already owe.  Both fills are raster multiples, so that difference
        moves in whole rasters and `delta` generally is not one -- hence `round` rather than
        `ceil`, which leaves the echo within half a raster of the midpoint instead of always
        past it.  See :attr:`symmetry_residual_s`.
        """
        delta = self._tail_min_s - self._head_min_s
        steps = round(delta / self._raster_s) * self._raster_s
        head, tail = (steps, 0.0) if steps >= 0 else (0.0, -steps)
        return self._head_floor_s + head, tail

    def _resolve_timing(self, te_s: float | None, tr_s: float | None) -> tuple[float, float]:
        """Return ``(te, tr)``, having set the two fills that realise them."""
        if te_s is None and tr_s is not None:
            head, tail = self._symmetric_fills_within(self._require_tr(tr_s))
        elif te_s is None:
            head, tail = self._symmetric_fills()
        else:
            head = self._require_te(te_s)
            tail = (0.0 if tr_s is None
                    else self._require_asymmetric_tr(tr_s, self._head_min_s + head))
        self._head_fill_s = head
        self._tail_fill_s = tail
        return (self._grad_start_s - self.exc.time_to_center() + self._joint_fill_s
                + self.ro.time_to_echo() + head), self._span_s(head, tail)

    def _symmetric_fills_within(self, tr_s: float) -> tuple[float, float]:
        """
        Split a requested TR into two fills that put the echo as near the midpoint as it goes.

        The caller fixed the total, so only the split is free: `head` is rounded to the raster
        and `tail` takes the remainder, which keeps the requested TR exact rather than trading
        it for a symmetry the raster cannot deliver anyway.
        """
        floor = self._head_floor_s
        total = tr_s - self._span_s(floor, 0.0)
        delta = self._tail_min_s - self._head_min_s
        head = min(max(round((total + delta) / 2.0 / self._raster_s) * self._raster_s, 0.0),
                   total)
        return floor + head, total - head

    def _require_te(self, te_s: float) -> float:
        """Return the head fill that reaches a requested TE, or refuse and say what can."""
        wanted = require_positive(te_s, 'te_s')
        if wanted < self.min_te_s - 1e-12:
            raise ConfigurationError(format_error(
                f'te_s = {wanted * 1e3:.3f} ms is shorter than this repetition can achieve.',
                {'te_s': wanted, 'min_te_s': self.min_te_s,
                 'min_symmetric_te_s': self.min_symmetric_te_s,
                 'bandwidth_hz_px': self.ro.bandwidth_hz_px},
                [f'pass te_s >= {self.min_te_s:.6g}',
                 'or te_s=None for the symmetric TE = TR/2',
                 'a higher bandwidth_hz_px shortens the readout, and with it min_te_s'],
            ))
        return self._head_floor_s + ceil_raster(wanted - self.min_te_s, self._raster_s)

    def _require_tr(self, tr_s: float) -> float:
        """Return a requested symmetric TR on the raster, or refuse and say what can."""
        wanted = require_positive(tr_s, 'tr_s')
        if wanted < self.min_tr_s - 1e-12:
            raise ConfigurationError(format_error(
                f'tr_s = {wanted * 1e3:.3f} ms is shorter than a symmetric repetition can '
                f'achieve.',
                {'tr_s': wanted, 'min_tr_s': self.min_tr_s,
                 'min_symmetric_te_s': self.min_symmetric_te_s},
                [f'pass tr_s >= {self.min_tr_s:.6g}',
                 'or tr_s=None for the shortest symmetric repetition',
                 'an explicit te_s asks for an asymmetric realisation, which relaxes this'],
            ))
        return ceil_raster(wanted, self._raster_s)

    def _require_asymmetric_tr(self, tr_s: float, te: float) -> float:
        """Return the tail fill that reaches a requested TR beside an explicit TE."""
        wanted = require_positive(tr_s, 'tr_s')
        floor = self._span_s(te - self._head_min_s, 0.0)
        if wanted < floor - 1e-12:
            raise ConfigurationError(format_error(
                f'tr_s = {wanted * 1e3:.3f} ms is shorter than te_s = {te * 1e3:.3f} ms and the '
                f'rest of the repetition need together.',
                {'tr_s': wanted, 'te_s': te, 'min_tr_for_this_te_s': floor},
                [f'pass tr_s >= {floor:.6g}',
                 'or a shorter te_s, which shortens the first half with it',
                 'or te_s=None, which lets this module choose the symmetric pair'],
            ))
        return ceil_raster(wanted - floor, self._raster_s)

    # ------------------------------------------------------- the jointly designed winder
    def _solve_with_scope(self, claims, probe_ro, fov_y, bandwidth_hz_px, partial_fourier,
                          readout_moment, opts) -> None:
        """
        Size the designed winder, then re-design it at the echo time it will actually be emitted
        at.

        Two passes, because the two are coupled: the window sets the shortest symmetric TE, and a
        waveform whose first moment is nulled at one echo time is not nulled at another.  The
        designer takes the echo as a request and reports the delay it inserted **between the
        winder and the readout** to reach it; putting the delay anywhere earlier would translate
        the solved waveform out from under its own moment.

        Iterated rather than done exactly twice, because a design made for a later echo may
        choose a different window and move the echo again.  The loop runs until the first half is
        at least as long as the second, which is the condition under which the symmetric solve
        needs no padding of its own -- and that is what keeps every microsecond of TE inside the
        designer's schedule instead of in front of it.
        """
        raster = float(opts.grad_raster_time)
        wanted: float | None = None
        for _ in range(8):
            window, fill, designs = self._design_joint(
                claims, probe_ro, fov_y, self.winder_s, wanted, opts)
            self.winder_s, self._joint_fill_s, self._joint = (
                ceil_raster(window, raster), fill, designs)
            probe = CartesianLine(opts=opts, fov_mm=self.fov_mm[0], matrix=self.matrix[0],
                                  axis='x', bandwidth_hz_px=bandwidth_hz_px,
                                  partial_fourier=partial_fourier,
                                  prephaser_duration_s=self.winder_s, **readout_moment)
            head, tail = self._halves_with(probe)
            if head >= tail - 1e-12:
                return
            wanted = ceil_raster(tail, raster)

    def _halves_with(self, probe: CartesianLine) -> tuple[float, float]:
        """``(first half, second half)`` of the repetition at the current window, seconds."""
        head = (self._grad_start_s - self.exc.time_to_center()
                + self._joint_fill_s + probe.time_to_echo())
        tail = (probe().duration - probe.time_to_echo() + self.winder_s + self._lead_s
                + self.exc.time_to_center())
        return head, tail

    def _route(self, flow_comp: FlowCompensation | None) -> dict[str, str]:
        """Which owner handles each requested axis.  Reads only an intent's `axes`."""
        routed: dict[str, str] = {}
        component = type(self).__name__
        for intent, axis in _augment.axes_claimed((flow_comp,)):
            owner = self._moment_owners.get(axis)
            if owner is None:
                _augment.refuse_unowned_axis(
                    component, intent, axis, self._moment_owners,
                    f'{axis!r} is not an axis this repetition plays an adjustable gradient on')
            else:
                routed[axis] = owner
        _augment.require_owners_can_serve(
            component, (flow_comp,), routed,
            lambda: tuple(a for a, o in self._moment_owners.items() if o == 'joint'))
        return routed

    def _encode(self, line: int) -> LogicBlock:
        """The phase-encode waveform: the ordinary blip, or the jointly designed one."""
        return self._joint_block('y', (line, _augment.ONE_STATE)) or self.pe(line=line)

    def _joint_block(self, axis: str, state: object) -> LogicBlock | None:
        """One axis' jointly designed waveform, or ``None`` when nothing claimed that axis."""
        designed = self._joint.get(axis)
        if designed is None:
            return None
        out = LogicBlock('joint')
        for at, event in _joint.events_for(designed, state):
            out.add(at - designed.schedule.window_start_s, event)
        return out

    def _design_joint(self, claims: Sequence[object], probe_ro: CartesianLine, fov_y: float,
                      local_min_s: float, te_request: float | None,
                      opts: Opts) -> tuple[float, float, dict[str, object]]:
        """
        Adapt this repetition to the shared physical designer, and return the window it chose.

        The designed region is the winder between the RF centre and the readout, exactly as it is
        for a spoiled gradient echo -- the balance condition lives outside it and is untouched by
        what happens inside, because the designer is given the same **total** area on each axis
        and only redistributes it in time.
        """
        geometry = _scope.ScopeGeometry(
            origin_s=self._lead_s + self.exc.time_to_center(),
            window_start_s=self._lead_s + self._grad_start_s,
            tail_s=probe_ro.time_to_echo() - probe_ro.prephaser_duration_s,
        )
        requirements = [
            self._axis_requirement(axis, group, probe_ro, fov_y, opts)
            for axis, group in _joint.group_by_axis(claims).items()
        ]
        designed = _scope.design_scope(geometry, requirements, opts, min_window_s=local_min_s,
                                       te_request_s=te_request)
        return designed.window_s, designed.fill_s, dict(designed.designs)

    def _axis_requirement(self, axis: str, group: Sequence[object], probe_ro: CartesianLine,
                          fov_y: float, opts: Opts) -> _scope.AxisRequirement:
        """What this axis wants at the echo, and what already plays between the two instants."""
        claimed = {getattr(claim, 'order', None) for claim in group}
        encode = PhaseEncode(opts=opts, fov_mm=fov_y, matrix=self.matrix[1], axis='y')
        base = ((lambda index: float(encode.k_per_m(index))) if axis == self._PE_AXIS
                else (lambda index: 0.0))
        indices = tuple(range(self.matrix[1])) if axis == self._PE_AXIS else (0,)

        targets: dict[object, tuple[float, float | None]] = {}
        for index in indices:
            resolved = {
                order: _joint.resolve_claims(
                    group, self._encoding_states, axis=axis, order=order,
                    base={key: (base(index) if order == 0 else 0.0)
                          for key in self._encoding_states})
                for order in _joint.ORDERS
            }
            for key in self._encoding_states:
                targets[(index, key)] = (
                    resolved[0][key], resolved[1][key] if 1 in claimed else None)

        readout = probe_ro if axis == probe_ro.axis else None
        # The leading winder and the excitation, as this repetition emits them.  `rephase=False`
        # always: this module plays both z lobes itself, so the excitation never carries one.
        before = (LogicBlock()
                  .add(0.0, self._z_lead)
                  .add(self._lead_s, self.exc(rephase=False)))

        def fixed(state: object, schedule: _joint.Schedule) -> tuple[float, float]:
            """
            What already plays on this axis between the semantic origin and the echo.

            Integrated from the RF **effective centre**, so the leading winder and the pre-centre
            half of the selection gradient fall outside it -- they act on magnetisation that does
            not exist yet.  They are in `before` so that the block is the one this repetition
            emits, and the integration bounds are what excludes them.
            """
            moments = [
                _joint.measure_moment(before, order, axis, origin_s=schedule.origin_s,
                                      start_s=schedule.origin_s, end_s=schedule.endpoint_s)
                for order in _joint.ORDERS
            ]
            if readout is not None:
                lobe = _joint.placed(readout.gx, schedule.window_start_s + schedule.window_s)
                for order in _joint.ORDERS:
                    moments[order] += _joint.measure_moment(
                        lobe, order, axis, origin_s=schedule.origin_s,
                        start_s=schedule.origin_s, end_s=schedule.endpoint_s)
            return (moments[0], moments[1])

        return _scope.AxisRequirement(
            axis=axis, states=tuple(targets), design_states=self._design_states(axis, targets),
            target=lambda state: targets[state], fixed=fixed,
        )

    def _design_states(self, axis: str, targets: Mapping[object, object]) -> tuple[object, ...]:
        """Which states to size the schedule from; the designer verifies the rest."""
        states = tuple(targets)
        if axis != self._PE_AXIS:
            return states
        lines = sorted({index for index, _key in states})
        ends = {lines[0], lines[-1]}
        return tuple(state for state in states if state[0] in ends)
