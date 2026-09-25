"""
:class:`GRE2DTR` -- one repetition of a spoiled 2D gradient echo.

``kernel/`` because it composes modules from more than one leaf folder, and because it is **the
repeating unit**.  ``readout/``'s membership rule -- *contains an ADC* -- would admit a whole
scan, so the leaf rules cannot tell a leaf from a composite.  Two more rules can: a kernel
composes leaves, and ``imaging/`` composes kernels.

This is the level worth reusing.  MPRAGE is an inversion followed by a train of GRE repetitions
and MP2RAGE is two such trains at different inversion times; neither reuses a whole 2D GRE, and
both reuse one TR.

What this layer owns, and what it does not
------------------------------------------
It owns the arithmetic no caller should repeat: the winder coupling, where each block goes
relative to ``time_to_center`` and ``time_to_echo``, TE and TR fill, the rewinder, the spoiler,
and the ``LIN`` label.  It does not own the loop, the line ordering, the sampling pattern or the
RF-spoiling schedule -- every one of those is a property of a whole acquisition, and belongs to
:class:`~seqcraft.modules.GRE2D`.

The winder is three axes wide
-----------------------------
The slice rephaser is on ``z``, the phase-encode blip on ``y`` and the readout prephaser on ``x``,
and **all three play at once**.  Overlap on different axes costs nothing in the tree and nothing
on the hardware, so waiting for the rephaser before starting the other two adds its whole duration
to TE and buys nothing.  On the reference protocol here that is 320 microseconds of TE, most of a
tenth of the echo time, given away by a placement decision rather than by any physics.

It is also why :meth:`~seqcraft.modules.Excitation.time_to_rephaser` and
:attr:`~seqcraft.modules.Excitation.rephaser_duration_s` exist: the excitation reports where its
own tail begins and how long it is, this layer starts the other two axes there, and all three are
stretched to the longest of them.  The rephaser is a *participant* in the winder coupling, not a
phase before it.

Spoiling is three axes wide too, and for a different reason
-----------------------------------------------------------
``spoil_axis`` defaults to ``('x', 'z')`` -- the readout and the slice -- because what a spoiled
GRE has to destroy is *transverse magnetisation carried from one repetition into the next*, and a
gradient only destroys it along its own direction.  Winding four cycles on ``z`` leaves a
distribution that is perfectly coherent in ``x``, so the residual FID and the stimulated echoes
that follow it survive; they arrive at the next readout with a ky that does not belong to it, and
the symptom is a ghost along the phase-encode direction rather than anything that looks like a
spoiling problem.  It is worst exactly where an inversion-prepared train is worst -- a long train
whose first repetitions carry far more magnetisation than the steady state.

The **readout axis** is the cheap one and the one most often missed.  After the echo, ``kx`` has
only half the readout's area left on it -- half a cycle across one voxel -- so the readout does
almost none of its own spoiling despite looking like the largest gradient in the repetition.

The **phase-encode axis** already has a gradient in the tail: the rewinder, which returns ``ky``
to zero so that every repetition starts from the same place.  Naming ``'y'`` in `spoil_axis` adds
a spoiler *beside* the rewinder rather than instead of it, and the compiler sums the two on ``y``
-- so the net is a fixed dephasing identical every repetition, rather than one that depends on
which line was just acquired.  It is off by default because the rewinder is what the phase-encode
axis is for, and because ``x`` and ``z`` between them already cover two directions.

Cycles are counted **per axis, across that axis's own voxel**: the slice thickness on ``z``, and
the in-plane voxel on ``x`` and ``y``.  Those differ by a factor of three on the reference
protocol here, so one shared area would be a different amount of spoiling on each axis while
reading as the same number.

Geometry is per axis here and scalar in the leaves
--------------------------------------------------
``matrix`` is ``(nx, ny)`` -- readout then phase, matching image-array order -- and ``fov_mm`` is
``(x, y)`` or a scalar meaning square.  The pair **stops here**.  ``PhaseEncode`` and
``CartesianLine`` keep a scalar ``matrix``, because each acts on exactly one axis and does not
know the other exists, which is what lets ``PhaseEncode`` serve ``y`` today and ``z`` when 3D
arrives.  This is the layer that knows the scan is two-dimensional, so this is the layer that
splits the pair.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pypulseq as pp

from ...augmentation import FlowCompensation, VelocityEncoding
from ...design.events import AXES
from ...design.logic import LogicBlock
from ...design.module import Module
from ...design.timing import EPS
from ...errors import ConfigurationError, format_error
from .. import _augment, _joint
from .._support import ceil_raster, require_axis, require_pair, require_positive
from ..encoding.phase_encoding import PhaseEncode
from ..readout.cartesian_line import CartesianLine
from ..rf.excitation import Excitation
from ..spoiler import spoiler

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from pypulseq.opts import Opts

__all__ = ['GRE2DTR']


class GRE2DTR(Module):
    """
    One repetition: excitation, phase encode, readout, rewinder, spoiler, TR fill.

    Parameters
    ----------
    opts
        The scanner.  Passed down to every part.
    fov_mm
        ``(x, y)`` in millimetres, or a scalar meaning square.
    matrix
        ``(nx, ny)`` -- readout then phase.
    thickness_mm
        Slice thickness.  Also the voxel dimension a ``z`` spoiler winds across; see
        :meth:`voxel_mm`.
    flip_deg
        Flip angle, degrees.
    rf_duration_s
        Excitation pulse length, seconds.  ``None`` is
        :class:`~seqcraft.modules.Excitation`'s own default, and the argument exists because this
        is usually the largest term in :attr:`min_te_s`: a 3 ms pulse puts 1.5 ms into TE before
        any gradient is played, which a short readout cannot shorten.  Forwarded unchanged, and it
        adds no arithmetic here -- the winder and TE measure the excitation rather than predict it.
    te_s, tr_s
        Requested echo and repetition times.  ``None`` is "as short as possible".  A request
        below :attr:`min_te_s` or :attr:`min_tr_s` raises rather than being silently lengthened;
        a request above is rounded **up** onto the gradient raster, and the achieved values are
        readable back as :attr:`te_s` and :attr:`tr_s`.
    bandwidth_hz_px
        Receive bandwidth per pixel.
    partial_fourier
        Readout-direction partial echo.  The phase-encode direction is undersampled by acquiring
        fewer lines, which is a *build* argument of :class:`~seqcraft.modules.GRE2D`.
    echoes, polarity, echo_spacing_s
        Forwarded unchanged to :class:`~seqcraft.modules.CartesianLine`, and **that is the entire
        change to this class**: a multi-echo GRE is this repetition with a readout that is read
        more than once.  `te_s` here stays the *first* echo's, which is what TE means; the array
        is ``tr.ro.te_s``.  Nothing else moves -- the winder coupling, :attr:`min_te_s`,
        :attr:`min_tr_s`, the rewinder and the spoiler all take a longer readout without a line
        of new arithmetic, because each of them measures the block rather than predicting it.
    spoil_cycles_per_voxel
        Gradient spoiling, in turns of phase wound across one voxel **along each spoiled axis** --
        the slice thickness on ``z``, the in-plane voxel on ``x`` and ``y``.
    spoil_axis
        One axis or several.  Default ``('x', 'z')``: a gradient spoils only along its own
        direction, so a single axis leaves the residual coherent in the other two.  Naming ``'y'``
        adds a spoiler beside the rewinder rather than instead of it.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    exc : Excitation
    pe : PhaseEncode
    ro : CartesianLine
        The three leaf modules, held as plain attributes.  Composition needs no API.
    spoilers : dict[str, LogicBlock]
        The spoiler on each axis of `spoil_axis`, keyed by axis.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
    ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> tr = GRE2DTR(opts=o, fov_mm=250.0, matrix=(128, 64), thickness_mm=5.0)
    >>> tr.center_line
    32
    >>> block = tr(line=40)
    >>> abs(block.duration - tr.tr_s) < 1e-12          # the TR fill makes the block exact
    True

    The label is emitted here, and only for a repetition that samples:

    >>> [n.item.value for n in block if getattr(n.item, 'type', '') == 'labelset']
    [40]
    >>> [n.item.value for n in tr(line=40, acquire=False)
    ...  if getattr(n.item, 'type', '') == 'labelset']
    []

    Notes
    -----
    **The ``LIN`` label is this layer's job, not ``CartesianLine``'s.**  A Cartesian line is
    identical every TR and does not know *which* line it is; the repetition does.  It is emitted
    into the block holding the ADC, so every readout carries the k-space line it belongs to --
    the independent record a reconstruction's gridding can be checked against.  ``acquire=False``
    suppresses it as well as the ADC: a dummy that emitted ``LIN`` would put a line index in the
    file for a readout that was never sampled.

    ``SLC`` joins it when multi-slice arrives.  Nothing else: ``AVG``, ``REP`` and ``SET`` have
    no meaning in this sequence, and a label nothing reads is noise in the file.

    **`center_mm` is two unrelated mechanisms on two different events**, and this is the layer
    that holds both because it is the layer that holds both events -- ``z`` becomes a frequency
    on the excitation's RF, ``x`` a frequency on the readout's ADC.  A slice offset moves *which
    spins are excited*, before any signal exists; a readout offset moves *where the received band
    sits*, after.  **`y` must be zero and a non-zero value raises**: a phase-encode shift is a
    linear phase across ``ky``, so it is a per-line ADC phase rather than a per-TR frequency.
    That is implementable and standard, and deliberately out of scope here; the argument is a
    three-tuple rather than a pair so that adding it later changes no signature.
    """

    #: The axis `PhaseEncode` acts on here.  A class constant because the routing table is built
    #: from it before any leaf exists: the coupled design runs first and may take that axis'
    #: winder over.
    _PE_AXIS = 'y'

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: float | tuple[float, float],
        matrix: tuple[int, int],
        thickness_mm: float,
        flip_deg: float = 15.0,
        rf_duration_s: float | None = None,
        te_s: float | None = None,
        tr_s: float | None = None,
        bandwidth_hz_px: float = 200.0,
        partial_fourier: float = 1.0,
        echoes: int = 1,
        polarity: str | None = None,
        echo_spacing_s: float | None = None,
        spoil_cycles_per_voxel: float = 4.0,
        spoil_axis: str | Iterable[str] = ('x', 'z'),
        flow_comp: FlowCompensation | None = None,
        velocity_encode: VelocityEncoding | None = None,
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        _augment.require_intent_type(flow_comp, FlowCompensation, 'flow_comp')
        _augment.require_intent_type(velocity_encode, VelocityEncoding, 'velocity_encode')
        fov_x, fov_y = require_pair(fov_mm, 'fov_mm')
        nx, ny = require_pair(matrix, 'matrix')
        self.fov_mm = (fov_x, fov_y)
        self.matrix = (int(nx), int(ny))
        self.thickness_mm = require_positive(thickness_mm, 'thickness_mm')
        self.spoil_cycles_per_voxel = require_positive(
            spoil_cycles_per_voxel, 'spoil_cycles_per_voxel',
        )
        self.spoil_axis = _spoil_axes(spoil_axis)

        # ``None`` rather than a repeated ``3e-3``: the pulse length is Excitation's decision and a
        # default copied to here is a second one that can disagree with it.  This layer only
        # forwards, and the pulse is worth forwarding because it is usually the largest term in
        # min_te_s -- at 800 Hz/px the readout is 1.3 ms and the default pulse is 3.
        pulse = {} if rf_duration_s is None else {'duration_s': rf_duration_s}
        self.exc = Excitation(opts=opts, flip_deg=flip_deg, thickness_mm=self.thickness_mm, **pulse)

        # The winder coupling, resolved once, over **three** participants.  The slice rephaser is
        # on z, the phase-encode blip on y and the readout prephaser on x, and all three play at
        # the same time -- overlap on different axes costs nothing.  Waiting for the rephaser
        # before starting the other two would add its whole duration to TE for nothing.
        #
        # Each leaf reports its own minimum and accepts an override, so the composite takes the
        # maximum and passes it down.  Stretching the short ones keeps TE at its minimum, which
        # inserting a delay would not -- and no leaf knows the other two exist.
        # The train's three arguments go to the probe as well as to the real one.  The prephaser
        # is unchanged by `echoes` -- there is still exactly one and it still cancels the *first*
        # lobe's pre-echo area -- but a bipolar train moves the dwell, and the dwell moves the
        # prephaser's area and therefore its minimum duration.  A probe built without them would
        # be measuring a different readout.
        train = dict(echoes=echoes, polarity=polarity, echo_spacing_s=echo_spacing_s)

        # Who owns a moment requirement on each axis, and therefore where an intent is routed.
        # `y` is this module's own winder; `x` is the readout's, which solves its first moment
        # itself; `z` is the slice rephaser, which `Excitation` realises and this module never
        # reshapes -- so it is absent, and a claim there is refused rather than silently dropped.
        self._moment_owners: dict[str, str] = {self._PE_AXIS: 'joint', 'x': 'readout'}
        routed = self._route(flow_comp, velocity_encode)

        # The readout route has to reach the probe as well as the real readout.  A compensated
        # prephaser is longer than an ordinary one, and the probe is what sets `winder_s` -- a
        # probe built uncompensated would size the winder for a waveform nobody plays.
        readout_moment = {'_null_moment_order': 1} if routed.get('x') else {}
        probe_ro = CartesianLine(opts=opts, fov_mm=fov_x, matrix=self.matrix[0], axis='x',
                                 bandwidth_hz_px=bandwidth_hz_px, partial_fourier=partial_fourier,
                                 **readout_moment, **train)
        probe_pe = PhaseEncode(opts=opts, fov_mm=fov_y, matrix=self.matrix[1], axis='y')
        self.winder_s = ceil_raster(
            max(probe_ro.prephaser_duration_s, probe_pe.min_duration_s,
                self.exc.rephaser_duration_s),
            opts.grad_raster_time,
        )
        # **Provisional, internal.**  When an augmentation claims a moment on the phase-encode
        # axis, that axis' winder stops being a blip this module can design alone: the encode
        # area and the claimed first moment constrain the same window, and the window's length
        # moves the echo every fixed moment is measured to.  So the coupled part is handed to
        # `_joint`, and what comes back is a window length folded into this module's own -- the
        # same way three local minima are already folded together above.
        jointly = {a for a, owner in routed.items() if owner == 'joint'}
        joint_claims, self.encoding_states = _augment.claims_and_states(
            _augment.flow_comp_for(flow_comp, jointly), velocity_encode)
        self._joint: dict[str, _joint.JointDesign] = {}
        if joint_claims:
            # AUTO first, so the minimum is known before an explicit TE is judged
            # against it -- and so the refusal for a too-short TE stays this module's,
            # naming `te_s`, rather than the designer's naming an axis.
            window, self._joint = self._design_joint(
                joint_claims, probe_ro, fov_y, self.winder_s, None, opts)
            self.winder_s = ceil_raster(window, opts.grad_raster_time)
            if te_s is not None and te_s >= self._reachable_te_s(probe_ro) - EPS:
                window, self._joint = self._design_joint(
                joint_claims, probe_ro, fov_y, self.winder_s, te_s, opts)
                self.winder_s = ceil_raster(window, opts.grad_raster_time)

        self.ro = CartesianLine(opts=opts, fov_mm=fov_x, matrix=self.matrix[0], axis='x',
                                bandwidth_hz_px=bandwidth_hz_px, partial_fourier=partial_fourier,
                                prephaser_duration_s=self.winder_s,
                                **readout_moment, **train)
        self.pe = PhaseEncode(opts=opts, fov_mm=fov_y, matrix=self.matrix[1], axis='y',
                              duration_s=self.winder_s)
        # One per axis, each winding its cycles across *that axis's* voxel.  The in-plane voxel is
        # 1.7 mm where the slice is 5 mm on the reference protocol, so a shared area would be
        # three times the spoiling on x that it is on z while reading as the same number.
        self.spoilers = {
            axis: spoiler(opts, cycles_per_voxel=self.spoil_cycles_per_voxel,
                          voxel_mm=self.voxel_mm(axis), axis=axis)
            for axis in self.spoil_axis
        }

        # Building is cheap -- the design happened above -- so measure rather than declare.
        # The winder starts where the slice rephaser does, not where the excitation block ends:
        # the transmit chain has to be clear of the pulse (calc_duration includes the ringdown),
        # and after that the three axes work together.
        self._winder_start_s = ceil_raster(
            max(self.exc.time_to_rephaser(), float(pp.calc_duration(self.exc.rf))),
            opts.grad_raster_time,
        )
        self._ro_duration_s = self.ro().duration
        self._tail_s = max(self.pe(line=0, rewind=True).duration,
                           *(block.duration for block in self.spoilers.values()))
        self._te_fill_s = self._resolve_te(te_s)
        self._tr_s = self._resolve_tr(tr_s)

    # ------------------------------------------------------------------ what it knows
    @property
    def center_line(self) -> int:
        """The phase-encode line that encodes ``k = 0``.  One convention, defined once."""
        return self.pe.center_line

    def voxel_mm(self, axis: str) -> float:
        """
        The voxel dimension along `axis`, millimetres.

        This is the layer that knows the scan is two-dimensional, so this is the layer that can
        answer it: ``x`` and ``y`` are ``fov / matrix`` on their own axis and ``z`` is the slice
        thickness.  It exists because a spoiler takes the voxel dimension **along its own axis**,
        and passing the in-plane size for a ``z`` spoiler under-spoils by the ratio of the two --
        which is faint residual banding rather than anything that looks like an error.
        """
        return {
            'x': self.fov_mm[0] / self.matrix[0],
            'y': self.fov_mm[1] / self.matrix[1],
            'z': self.thickness_mm,
        }[require_axis(axis)]

    @property
    def min_te_s(self) -> float:
        """
        Shortest achievable echo time, seconds.  A feasibility fact known at design time.

        **TE means the first echo**, in a multi-echo protocol as in any other, which is what
        ``ro.time_to_echo()`` already returns.  The rest of the train follows at
        ``ro.echo_spacing_s`` and a caller reads the whole array off ``ro.te_s`` -- so `echoes`
        adds no arithmetic here, and it adds none to :attr:`min_tr_s` either, because that is
        built from a block that measures itself.
        """
        return self._winder_start_s - self.exc.time_to_center() + self.ro.time_to_echo()

    @property
    def te_s(self) -> float:
        """The echo time this repetition achieves, seconds."""
        return self.min_te_s + self._te_fill_s

    @property
    def min_tr_s(self) -> float:
        """Shortest repetition time **for the configured TE**, seconds, on the block raster."""
        return ceil_raster(self._tail_start_s + self._tail_s, self.opts.block_duration_raster)

    @property
    def tr_s(self) -> float:
        """The repetition time this block occupies, seconds."""
        return self._tr_s

    def time_to_echo(self) -> float:
        """
        Seconds from the start of this module's block to k = 0.

        The same question :meth:`CartesianLine.time_to_echo` answers, asked of the whole
        repetition.  TE is this minus :meth:`Excitation.time_to_center`, which is exactly what
        :attr:`te_s` reports.
        """
        return self._readout_start_s + self.ro.time_to_echo()

    # ----------------------------------------------------------------------- assembly
    def build(
        self,
        *,
        line: int,
        phase_deg: float = 0.0,
        acquire: bool = True,
        center_mm: tuple[float, float, float] = (0.0, 0.0, 0.0),
        encoding_state: object = None,
    ) -> LogicBlock:
        """
        Return one repetition.

        Parameters
        ----------
        line
            Zero-based phase-encode index.  The two things that vary between repetitions are
            this and `phase_deg`.
        phase_deg
            RF carrier phase for this repetition, degrees.
        acquire
            ``False`` drops the ADC and the ``LIN`` label and changes nothing else, so a dummy
            loads the gradients exactly as a real repetition does.
        center_mm
            ``(x, y, z)`` centre of the imaging volume, millimetres.  ``y`` must be ``0.0``.
        """
        x, y, z = (float(v) for v in center_mm)
        if y != 0.0:
            msg = format_error(
                f'center_mm y = {y:g} mm is out of scope: a phase-encode shift is a per-line '
                f'ADC phase, not a per-TR frequency.',
                {'center_mm': center_mm},
                [
                    'shift the slice with z, or the readout FOV with x',
                    'a phase-encode offset needs adc.phase_offset += 2*pi*k_y(line)*y, which '
                    'this module deliberately does not do yet',
                ],
            )
            raise ConfigurationError(msg)

        state = self._resolve_state(encoding_state)
        start = self._readout_start_s
        tail = self._tail_start_s
        out = (
            LogicBlock()
            .add(0.0, self.exc(phase_deg=phase_deg, position_mm=z))
            .add(start, self._encode(line, state))
            # The same phase to both: the receiver is phase-locked to the transmitter, so an
            # RF-spoiling schedule that moves one and not the other writes its quadratic phase
            # into ky.  This is the layer that holds both events, so this is where they agree.
            .add(start, self.ro(acquire=acquire, phase_deg=phase_deg, offset_mm=x))
            # The rewinder stays whatever `spoil_axis` says: it returns ky to zero so that every
            # repetition starts from the same place, and a y spoiler beside it sums to a *fixed*
            # dephasing rather than to one that depends on the line just acquired.
            .add(tail, self.pe(line=line, rewind=True))
        )
        for block in self.spoilers.values():
            out.add(tail, block)
        if acquire:
            out.add(start, pp.make_label(type='SET', label='LIN', value=int(line)))
        # The block measures itself, so the fill is what makes its duration exactly TR -- which
        # is what lets the layer above stack repetitions at n * tr_s and nothing else.  Snapped
        # onto the block raster because every term already is one, and the last bit of a float
        # subtraction is otherwise enough to make two identical sequences hash differently.
        fill = ceil_raster(self.tr_s - tail - self._tail_s, self.opts.block_duration_raster)
        if fill > 1e-9:
            out.add(tail + self._tail_s, pp.make_delay(fill))
        return out

    def _resolve_state(self, encoding_state: object) -> object:
        """
        Turn what the caller passed into the state key this repetition designed against.

        Refuses in **both** directions, because silence is the worse failure either way.  A
        missing state would realise an arbitrary one of a pair; an ignored state hands the caller
        two identical repetitions, which they find out about when the subtraction comes back
        zero.
        """
        states = self.encoding_states
        if states == (_augment.ONE_STATE,):
            if encoding_state is not None:
                _augment.refuse_state_mismatch(type(self).__name__, states, encoding_state)
            return _augment.ONE_STATE
        if encoding_state not in states:
            _augment.refuse_state_mismatch(type(self).__name__, states, encoding_state)
        return encoding_state

    # --------------------------------------------------------------------- capability
    @property
    def _joint_axes(self) -> tuple[str, ...]:
        """The axes this repetition designs jointly, derived from who owns what."""
        return tuple(a for a, owner in self._moment_owners.items() if owner == 'joint')

    def _route(self, flow_comp: FlowCompensation | None,
               velocity_encode: VelocityEncoding | None) -> dict[str, str]:
        """
        Which owner handles each intent, or a refusal naming the axis.

        Reads nothing about an intent but its `axis` and, for the one case that needs it, whether
        it is a difference between states.  There is no branch on augmentation class here, and a
        third augmentation would need no change to this method.
        """
        routed: dict[str, str] = {}
        component = type(self).__name__
        for intent, axis in [(i, a) for i in (flow_comp, velocity_encode) if i is not None
                             for a in i.axes]:
            owner = self._moment_owners.get(axis)
            if owner is None:
                _augment.refuse_unowned_axis(
                    component, axis, self._moment_owners, self._why_unowned(axis))
            elif owner != 'joint' and isinstance(intent, VelocityEncoding):
                # Not an impossibility: the readout's local solve nulls its first moment rather
                # than aiming it at a value, and a difference between states needs a target.
                _augment.refuse_unsupported_route(component, axis, self._joint_axes)
            else:
                routed[axis] = owner
        return routed

    def _why_unowned(self, axis: str) -> str:
        """The sequence-specific reason an axis has no owner here."""
        if axis == 'z':
            return ('this repetition\'s z gradient is the slice rephaser, which Excitation '
                    'realises for itself; GRE3DTR does own z, because its z winder carries the '
                    'partition encode')
        return f'{axis!r} is not an axis this repetition plays an adjustable gradient on'

    def _encode(self, line: int, encoding_state: object) -> LogicBlock:
        """The phase-encode waveform: the ordinary blip, or the jointly designed one."""
        return self._joint_block(self._PE_AXIS, (line, encoding_state)) or self.pe(line=line)

    def _joint_block(self, axis: str, state: object) -> LogicBlock | None:
        """One axis' jointly designed waveform, or ``None`` when nothing claimed that axis."""
        designed = self._joint.get(axis)
        if designed is None:
            return None
        out = LogicBlock('joint')
        for at, event in _joint.events_for(designed, state):
            out.add(at - designed.schedule.window_start_s, event)
        return out

    def _reachable_te_s(self, probe_ro: CartesianLine) -> float:
        """The echo time the current winder reaches, before any explicit-TE fill."""
        raster = float(self.opts.grad_raster_time)
        start = ceil_raster(
            max(self.exc.time_to_rephaser(), float(pp.calc_duration(self.exc.rf))), raster)
        return (start + self.winder_s + probe_ro.time_to_echo()
                - probe_ro.prephaser_duration_s - self.exc.time_to_center())

    # ------------------------------------------------------------- the coupled design
    def _design_joint(self, claims: Sequence[object], probe_ro: CartesianLine, fov_y: float,
                      local_min_s: float, te_request: float | None,
                      opts: Opts) -> tuple[float, dict[str, _joint.JointDesign]]:
        """
        Design every claimed axis against **one common candidate schedule**, and return the
        window that serves them all.

        Two things this must not do, both of which change the first moment whenever ``m0`` is
        non-zero, because ``m1' = m1 - dt m0``:

        * solve one axis at *its* minimum and let another axis' minimum widen the winder
          afterwards -- the echo moves and the solved moment goes stale;
        * design at the minimum TE and then let an explicit ``te_s`` insert fill in front of the
          winder -- the whole waveform translates, which is the same error by another route.

        So a candidate is a window **and** the fill an explicit TE would need at that window, and
        every axis is solved against the schedule that results.  A candidate is accepted only if
        every claimed axis is feasible at it.
        """
        raster = float(opts.grad_raster_time)
        origin = self.exc.time_to_center()
        start = ceil_raster(
            max(self.exc.time_to_rephaser(), float(pp.calc_duration(self.exc.rf))), raster)
        echo_in_lobe = probe_ro.time_to_echo() - probe_ro.prephaser_duration_s
        excitation = self.exc()
        by_axis = _joint.group_by_axis(claims)
        problems = {axis: self._axis_problem(axis, group, probe_ro, fov_y, excitation, opts)
                    for axis, group in by_axis.items()}

        exhausted = True
        for steps in range(int(round(local_min_s / raster)), _joint.SEARCH_LIMIT_WINDOWS):
            window = steps * raster
            reachable_te = start + window + echo_in_lobe - origin
            fill = 0.0 if te_request is None else te_request - reachable_te
            if fill < -EPS:
                # This window already overshoots the requested TE; a longer one only overshoots
                # further, so the request is infeasible and the kernel's own refusal will say so.
                exhausted = False
                break
            schedule = _joint.Schedule(
                origin_s=origin,
                endpoint_s=start + max(fill, 0.0) + window + echo_in_lobe,
                window_start_s=start + max(fill, 0.0),
                window_s=window,
            )
            designs = {}
            for axis, problem in problems.items():
                found = _joint.attempt(problem, schedule, opts)
                if found is None:
                    break
                designs[axis] = found
            if len(designs) == len(problems):
                return window, designs
        # Two different failures, and only one of them is about physics.  Running out of
        # candidate windows means the ceiling was reached without trying what lies beyond it;
        # saying "infeasible" there would claim something the search never established.
        first = next(iter(problems.values()))
        if exhausted:
            return _joint.refuse_search_exhausted(
                first, steps=_joint.SEARCH_LIMIT_WINDOWS, raster_s=raster)
        return _joint.refuse_infeasible(first)

    def _axis_problem(self, axis: str, group: Sequence[object], probe_ro: CartesianLine,
                      fov_y: float, excitation: LogicBlock, opts: Opts) -> _joint.JointProblem:
        """
        One axis' physical problem: what is wanted per state, and what already plays on it.

        **No augmentation is named here.**  The claims arrive resolved into absolute targets by
        :func:`~seqcraft.modules._joint.resolve_claims`, so a third augmentation needs no change
        to this method -- which is the whole point of the boundary.
        """
        claimed = {getattr(claim, 'order', None) for claim in group}
        base_m0 = self._base_moment(axis, fov_y, probe_ro, opts)
        indices = tuple(range(self.matrix[1])) if axis == self._PE_AXIS else (0,)
        targets: dict[object, tuple[float, float | None]] = {}
        for index in indices:
            resolved = {
                order: _joint.resolve_claims(
                    group, self.encoding_states, axis=axis, order=order,
                    base={key: (base_m0(index) if order == 0 else 0.0)
                          for key in self.encoding_states})
                for order in _joint.ORDERS
            }
            for key in self.encoding_states:
                targets[(index, key)] = (
                    resolved[0][key], resolved[1][key] if 1 in claimed else None)

        readout = probe_ro if axis == probe_ro.axis else None

        def fixed(state: object, schedule: _joint.Schedule) -> tuple[float, float]:
            """What already plays on this axis between the two instants, as emitted."""
            moments = [
                _joint.measure_moment(excitation, order, axis, origin_s=schedule.origin_s,
                                      start_s=0.0, end_s=schedule.endpoint_s)
                for order in _joint.ORDERS
            ]
            if readout is not None:
                # The readout lobe's own contribution up to the echo, placed where this schedule
                # puts it.  Measured from its knots rather than derived, so partial Fourier and
                # ramp sampling come out of it without a second derivation.
                lobe = _joint.placed(
                    readout.gx, schedule.window_start_s + schedule.window_s)
                for order in _joint.ORDERS:
                    moments[order] += _joint.measure_moment(
                        lobe, order, axis, origin_s=schedule.origin_s,
                        start_s=0.0, end_s=schedule.endpoint_s)
            return (moments[0], moments[1])

        return _joint.JointProblem(axis=axis, targets=targets, fixed=fixed)

    def _base_moment(self, axis: str, fov_y: float, probe_ro: CartesianLine, opts: Opts):
        """What the sequence itself wants at order 0 on `axis`, per index."""
        if axis == self._PE_AXIS:
            encode = PhaseEncode(opts=opts, fov_mm=fov_y, matrix=self.matrix[1], axis=axis)
            return lambda index: float(encode.k_per_m(index))
        # Readout and slice axes want k = 0 at the echo; what gets them there is the fixed part.
        return lambda index: 0.0

    # ------------------------------------------------------------------------ timing
    @property
    def _readout_start_s(self) -> float:
        """When the winder -- slice rephaser, phase-encode blip, readout prephaser -- starts."""
        return self._winder_start_s + self._te_fill_s

    @property
    def _tail_start_s(self) -> float:
        """When the rewinder and the spoiler start."""
        return self._readout_start_s + self._ro_duration_s

    def _resolve_te(self, te_s: float | None) -> float:
        """Return the delay that lengthens TE to the request, or zero for the shortest."""
        if te_s is None:
            return 0.0
        wanted = require_positive(te_s, 'te_s')
        if wanted < self.min_te_s - 1e-12:
            msg = format_error(
                f'te_s = {wanted * 1e3:.3f} ms is shorter than this repetition can achieve.',
                {'te_s': wanted, 'min_te_s': self.min_te_s,
                 'bandwidth_hz_px': self.ro.bandwidth_hz_px},
                [
                    f'pass te_s >= {self.min_te_s:.6g}',
                    'or te_s=None for the shortest achievable echo time',
                    'a higher bandwidth_hz_px shortens the readout, and with it min_te_s',
                ],
            )
            raise ConfigurationError(msg)
        return ceil_raster(wanted - self.min_te_s, self.opts.grad_raster_time)

    def _resolve_tr(self, tr_s: float | None) -> float:
        """
        Return the repetition time this block will occupy, seconds.

        The achieved value rather than a fill: ``min_tr_s + fill`` and ``ceil(requested)`` are
        the same number in exact arithmetic and differ in the last bit in float, which is enough
        to make two identically-specified sequences hash differently.
        """
        if tr_s is None:
            return self.min_tr_s
        wanted = require_positive(tr_s, 'tr_s')
        if wanted < self.min_tr_s - 1e-12:
            msg = format_error(
                f'tr_s = {wanted * 1e3:.3f} ms is shorter than this repetition can achieve.',
                {'tr_s': wanted, 'min_tr_s': self.min_tr_s, 'te_s': self.te_s},
                [
                    f'pass tr_s >= {self.min_tr_s:.6g}',
                    'or tr_s=None for the shortest achievable repetition time',
                    'a shorter te_s shortens min_tr_s with it',
                ],
            )
            raise ConfigurationError(msg)
        return ceil_raster(wanted, self.opts.block_duration_raster)


def _spoil_axes(value: str | Iterable[str]) -> tuple[str, ...]:
    """
    Return `spoil_axis` as an ordered, de-duplicated tuple of logical axes.

    A bare string is one axis rather than an iterable of characters, which is the trap in accepting
    both: ``spoil_axis='xz'`` would otherwise silently mean ``('x', 'z')`` on one implementation
    and raise on another, so it raises here and the tuple is written out.

    Examples
    --------
    >>> _spoil_axes('z')
    ('z',)
    >>> _spoil_axes(('z', 'x', 'z'))
    ('x', 'z')
    >>> _spoil_axes(())
    Traceback (most recent call last):
        ...
    seqcraft.errors.ConfigurationError: spoil_axis names no axis, so nothing would be spoiled.
    >>> _spoil_axes('xz')
    Traceback (most recent call last):
        ...
    seqcraft.errors.ConfigurationError: spoil_axis must be one of 'x', 'y', 'z', got 'xz'.
    """
    names = (value,) if isinstance(value, str) else tuple(value)
    if not names:
        msg = format_error(
            'spoil_axis names no axis, so nothing would be spoiled.',
            {'spoil_axis': value},
            [
                "spoil_axis=('x', 'z') is the default: a gradient spoils only along its own "
                'direction',
                'RF spoiling is a separate mechanism, and rf_spoil=False turns that one off',
            ],
        )
        raise ConfigurationError(msg)
    ordered = sorted({require_axis(name, 'spoil_axis') for name in names}, key=AXES.index)
    return tuple(ordered)
