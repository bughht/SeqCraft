r"""
:class:`Refocusing` -- a 180 and its crusher pair, as one waveform symmetric about the pulse.

``rf/`` because the pulse's ``use`` is ``'refocusing'``.  That is the folder's second membership
value and until now it had no member; MRzero maps the same field onto its ``PulseUsage``, and
pypulseq's ``calculate_kspacePP`` reads it to decide where to **conjugate k**.  Which is the whole
reason this module is more than a wrapper around ``make_sinc_pulse``.

The one rule
------------
A refocusing pulse conjugates k, so between consecutive refocusing centres each axis's gradient
area *before* the echo must equal its area *after* it.  Writing :math:`t^{\text{ref}}_n` for the
effective centre of pulse *n*,

.. math:: k_n = -\big(k_{n-1} + A(\text{echo}_{n-1} \to t^{\text{ref}}_n)\big)
                 + A(t^{\text{ref}}_n \to \text{echo}_n)

so :math:`k_n = k_{n-1} = 0` for every *n* exactly when the two areas are equal.  Nothing else is
required, and nothing less will do.

On the selection axis both crushers straddle **this module's own pulse**, so this module can
balance them itself -- and it is the only axis where that is true.  On the readout axis the lobe
*after* readout *n* and the lobe *before* readout *n+1* are the pair that has to balance, and they
live in two different readout blocks, so the caller places those.

This module therefore crushes on the **selection axis only**, which is also where a crusher is
actually wanted: it is what dephases the FID an imperfect 180 leaves behind.  A readout-axis crusher
is a free choice rather than a requirement -- the conjugation cancels whatever the two lobes share,
so its area has no k-space consequence at all and it is paid for in echo spacing.  ``writeTSE.m``
adds one (``rSpoilFac``); ``examples/se_2d/01`` measures what it costs and declines it.

What "equal-area crushers" does not mean
----------------------------------------
Not "the same trapezoid twice".  The balance point is the RF's **effective centre**, and the
selection plateau's own halves are unequal whenever

* the transmit dead time and the ringdown differ -- the plateau covers both, the pulse sits between
  them, and ``(dead - ringdown)/2`` of plateau lands on one side only; or
* the pulse is asymmetric -- a minimum-phase SLR 180's centre is nowhere near its midpoint.

On a scanner with ``rf_dead_time = 100 us`` and ``rf_ringdown_time = 30 us``, a 4 ms TBW-4 sinc
over a 6.25 mm slab leaves a residual of ``a_sel * (dead - ringdown) / 2`` = **5.6 1/m**, 0.035
cycles across the slab.  Small -- and it **alternates sign echo to echo**, because
:math:`k_n = -k_{n-1} + \delta`, which is exactly the odd/even modulation an FSE is famous for and
which reads as a hardware fault.  The pulseq and pypulseq TSE demos get it right *by accident*:
they set the dead time and the ringdown to the same 100 us, so the plateau is symmetric.  Change
either and the sequence is quietly wrong.

Two things fix it and this module does both, because they fix different halves of the invariant:

1. **The plateau is symmetrised about the effective centre**, rather than merely made long enough
   to cover the pulse: ``t_plateau = 2 * max(dead + centre, duration - centre + ringdown)``.  That
   is what makes ``time_to_center() == duration / 2``, the *time* half, which the area solve cannot
   buy and which is what puts the echo at the midpoint between two refocusing pulses.  It costs
   ``|dead - ringdown|`` for a symmetric pulse and twice the pulse's own asymmetry for a
   minimum-phase one; :attr:`plateau_padding_s` reports it rather than hiding it.
2. **The crusher areas are then solved by integration**, not by formula: build the candidate,
   integrate to :meth:`time_to_center` with :func:`~seqcraft.modules._support.area_until`, and set
   the two amplitudes so that ``area_until(g, centre) == total - area_until(g, centre)``.  With a
   symmetrised plateau that returns two *equal* amplitudes, so step 2 looks redundant -- and it is
   exactly the check that says step 1 worked, on a pulse shape nobody has tried yet.  Deriving the
   amplitudes and asserting the areas from one formula would compare a number with itself.

Why the gradient is one event
-----------------------------
``writeTSE.m`` hand-splits the selection axis into ``GS1 ... GS7`` because a pulseq block holds one
gradient per axis, so a waveform spanning an RF has to be cut at every block boundary by hand.
:func:`seqcraft.compile` does that cutting and carries ``first``/``last`` across each seam, so what
this module emits is **one continuous waveform** that starts and ends at zero --
``sc.compile(sc.LogicBlock('probe').add(0.0, refoc()), opts)`` is therefore a complete contract
check.  A module that only works when something happens to be beside it is not reusable.

The B1 refusal, which lives here rather than in the compiler
------------------------------------------------------------
``sc.compile`` checks gradient amplitude and slew and **does not check RF amplitude at all**.  Peak
B1 scales with flip angle at a fixed shape, so a 180 needs exactly twice a 90's: measured, a 2 ms
sinc 180 at TBW 4 asks **130 % of a 20 uT ``max_b1``** and 2.7 ms is the floor.  A 180 is where
that first bites, so it is refused here, with the duration that fixes it.
:class:`~seqcraft.modules.Excitation` has the same hole and it is deliberately not patched from
here: widening a new module's change to fix an old one is how a one-module change becomes five.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

import numpy as np
import pypulseq as pp

from ...design.events import derive
from ...design.logic import LogicBlock
from ...design.module import Module
from ...design.timing import EPS
from ...errors import ConfigurationError, SeqCraftWarning, format_error
from .._support import area_until, ceil_raster, require_axis, require_positive, shift_slice

if TYPE_CHECKING:
    from pypulseq.opts import Opts

    from ...design.events import Event

__all__ = ['Refocusing']

#: The four pulse shapes, and the pypulseq factory each names -- the same table
#: :class:`~seqcraft.modules.Excitation` resolves, for the same reason: a pypulseq missing one
#: fails through :func:`seqcraft._compat.require`'s single complete list rather than through an
#: import-time ``AttributeError`` here.
_FACTORIES: dict[str, str] = {
    'sinc': 'make_sinc_pulse',
    'slr': 'make_slr_pulse',
    'gauss': 'make_gauss_pulse',
    'block': 'make_block_pulse',
}

#: Arguments this module owns; `pulse_opts` may not reach past it to set one of them.
_RESERVED = frozenset({
    'delay', 'duration', 'flip_angle', 'freq_offset', 'phase_offset', 'return_delay',
    'return_gz', 'slice_thickness', 'system', 'time_bw_product', 'use',
})

#: Iterations the crusher-amplitude fixed point is allowed.  It is monotone -- a longer ramp lowers
#: the amplitude, which shortens the ramp it needs -- so this bounds a loop that converges in three
#: or four, rather than trading accuracy for time.
_MAX_SOLVER_STEPS = 12

#: Amplitude agreement, Hz/m, that ends the fixed point.  Far below one part in 10^6 of any crusher
#: this designs, and the areas are measured by integration afterwards regardless.
_SOLVER_TOL_HZ_M = 1e-6


class Refocusing(Module):
    """
    A refocusing pulse and its crusher pair, balanced in time and in area about the pulse.

    Parameters
    ----------
    opts
        The scanner.
    thickness_mm
        Slice thickness, millimetres, or ``None`` for a non-selective pulse.  **No default**, as
        for :class:`~seqcraft.modules.Excitation` and :class:`~seqcraft.modules.IRPrep`.  A caller
        normally wants this 20-25 % thicker than the excitation's, because a refocusing profile's
        edges are worse than an excitation's and a slice refocused only in its middle loses the
        edges twice per echo.  That factor belongs to whatever holds both pulses.
    flip_deg
        Flip angle, degrees.  An argument rather than a constant because 120-150 degree trains are
        standard for SAR and are the entry point for variable-flip-angle work.  Contrast
        :class:`~seqcraft.modules.IRPrep`, where 180 is *not* an argument, because an inversion by
        another angle is a different preparation.
    duration_s
        Pulse duration.  ``4e-3`` rather than ``Excitation``'s ``3e-3``, and that is the point: a
        180 needs exactly twice a 90's peak B1 at the same shape.  Over ``opts.max_b1`` raises,
        naming the duration that fixes it.
    time_bw_product
        ``None`` defers to the factory's own default, which differs by shape.  Passing it with
        ``pulse='block'`` raises.
    pulse
        One of ``'sinc'``, ``'slr'``, ``'gauss'``, ``'block'``.  An unknown value raises, listing
        them.
    pulse_opts
        Forwarded to the chosen factory -- ``{'filter_type': 'pm'}`` for an SLR pulse.  An
        unrecognised key raises, naming what this `pulse` accepts.
    crush_cycles_per_voxel
        Turns of phase wound across one voxel along `crush_axis`, **counted from the start of the
        block to the RF's effective centre** -- so on the selection axis it includes the selection
        plateau's own half.  That is what makes it the number the balance is solved against rather
        than one lobe's area.  A value below the plateau's own contribution -- about
        ``time_bw_product / 2`` cycles across the slice -- makes the leading crusher change sign:
        legal, no longer anything a reader would recognise as a crusher, and warned about with the
        number.
    crush_voxel_mm
        The voxel dimension along `crush_axis` that the cycles are counted across.  Defaults to
        `thickness_mm`, and is **required** when the pulse is non-selective, exactly as
        :class:`~seqcraft.modules.IRPrep`'s ``spoil_voxel_mm`` is and for the same reason.
    crush_duration_s
        The crusher window, one side.  ``None`` is :attr:`min_crush_duration_s`.  A composite
        passes the **maximum over three axes** here, in the same shape
        :class:`~seqcraft.modules.CartesianLine`'s ``prephaser_duration_s`` and
        :class:`~seqcraft.modules.PhaseEncode`'s ``duration_s`` take.
    crush_axis
        Defaults to the selection axis, which is the only axis whose crusher this module can fuse
        with anything.  Crushing on ``'x'`` is the readout's business and on ``'y'`` is
        :class:`~seqcraft.modules.PhaseEncode`'s; naming either here would put two gradients on one
        axis in one interval for no gain.  With no selection axis it defaults to ``'z'``, matching
        ``IRPrep.spoil_axis``.
    rise_time_s
        ``None`` derives each transition's ramp from ``opts.max_slew`` by a short fixed point.
        **It buys no echo spacing at all** -- measured, exactly zero, because the spacing depends on
        the crusher *window* and not on its ramps -- so this argument exists for slew headroom.  The
        references' fixed ``250e-6`` is that number frozen at one system.
    axis
        Selection axis.  ``None`` means ``'z'`` when selective, and is the only legal value when
        not.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    rf : Event
        The pulse, with its ``delay`` already carrying the lead-in that puts its effective centre
        at the middle of the block.
    gz : Event or None
        The selection gradient **as the factory designed it**, kept for its amplitude, which is
        what a slice offset is biased against.  It is not what gets emitted when the crusher shares
        its axis: that is one fused waveform, in :attr:`gradients`.
    gradients : tuple[tuple[float, Event], ...]
        ``(start, event)`` for every gradient this module emits, in the order ``build`` adds them.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s',
    ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> refoc = Refocusing(opts=o, thickness_mm=6.25, crush_voxel_mm=5.0)
    >>> refoc()                                   # the pulse, and one continuous waveform
    LogicBlock(Refocusing, 2 nodes, 5.60 ms)

    The invariant, both halves of it:

    >>> abs(refoc.time_to_center() - refoc().duration / 2) < 1e-12
    True
    >>> round(refoc.area_to_center_per_m, 6), round(refoc.area_from_center_per_m, 6)
    (600.0, 600.0)

    The crusher window this asks for is the number a composite maximises over three axes:

    >>> round(refoc.min_crush_duration_s * 1e6)
    700

    Where the trailing crusher starts is the earliest another axis may begin work, and putting the
    readout there rather than after the block is worth 1.4 ms per echo:

    >>> round(refoc.time_to_crusher() * 1e3, 3), round(refoc().duration * 1e3, 3)
    (4.9, 5.6)

    Symmetrising the plateau costs exactly ``|dead - ringdown|`` for a symmetric pulse:

    >>> round(refoc.plateau_padding_s * 1e6)
    70

    Notes
    -----
    **`time_to_center` is the conjugation instant**, and every echo time in a spin echo is measured
    from it.  Same definition as ``Excitation``'s and ``IRPrep``'s -- ``rf.delay +
    pp.calc_rf_center(rf)[0]``, from the start of the block ``build`` returns -- because these
    compose by addition only if they share an origin.

    **The crusher size has no k-space consequence.**  The conjugation cancels whatever the two
    lobes share, so what a larger crusher costs is the *window*, and through it the echo spacing,
    rather than accuracy.  What it buys is dephasing of the FID an imperfect 180 leaves behind.

    **`axis` refuses any value when non-selective; `position_mm` refuses only a non-zero one.**
    The rule ``Excitation`` and ``IRPrep`` both follow: ``None`` for an argument with no identity
    value, and the identity value itself when there is one.
    """

    def __init__(
        self,
        *,
        opts: Opts,
        thickness_mm: float | None,
        flip_deg: float = 180.0,
        duration_s: float = 4e-3,
        time_bw_product: float | None = None,
        pulse: str = 'sinc',
        pulse_opts: dict[str, Any] | None = None,
        crush_cycles_per_voxel: float = 3.0,
        crush_voxel_mm: float | None = None,
        crush_duration_s: float | None = None,
        crush_axis: str | None = None,
        rise_time_s: float | None = None,
        axis: str | None = None,
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.pulse = self._check_pulse(pulse)
        self.flip_deg = require_positive(flip_deg, 'flip_deg')
        self.duration_s = require_positive(duration_s, 'duration_s')
        self.selective = thickness_mm is not None
        self.thickness_mm = (
            require_positive(thickness_mm, 'thickness_mm') if self.selective else None
        )
        self.axis = self._check_axis(axis)
        self.crush_axis = require_axis(
            (self.axis or 'z') if crush_axis is None else crush_axis, 'crush_axis',
        )
        self.crush_cycles_per_voxel = require_positive(
            crush_cycles_per_voxel, 'crush_cycles_per_voxel',
        )
        self.crush_voxel_mm = self._check_crush_voxel(crush_voxel_mm)
        self.rise_time_s = (
            None if rise_time_s is None else require_positive(rise_time_s, 'rise_time_s')
        )
        self._check_combination(time_bw_product)

        kwargs: dict[str, Any] = {
            'flip_angle': np.deg2rad(self.flip_deg),
            'duration': self.duration_s,
            'system': opts,
            'use': 'refocusing',
            'delay': float(opts.rf_dead_time),
        }
        if time_bw_product is not None:
            kwargs['time_bw_product'] = require_positive(time_bw_product, 'time_bw_product')
        if self.selective:
            kwargs['slice_thickness'] = self.thickness_mm / 1e3
            kwargs['return_gz'] = True
        kwargs.update(self._check_pulse_opts(pulse_opts))

        factory = getattr(pp, _FACTORIES[self.pulse])
        with warnings.catch_warnings():
            # pypulseq warns about max_b1 and hands the pulse back anyway.  This module raises
            # instead, four lines below, with the duration that fixes it -- so the warning would
            # only be noise ahead of a better message.
            warnings.simplefilter('ignore', UserWarning)
            if self.selective:
                # The third return is a rephaser, dropped on purpose: the crusher pair takes its
                # place, and the balance those two are solved for is what rephases the slice.
                rf, gz, _rephaser = factory(**kwargs)
                self.gz = derive(gz, channel=self.axis) if self.axis != 'z' else gz
            else:
                rf, self.gz = factory(**kwargs), None
        self._check_b1(rf)

        #: Area on `crush_axis` between the start of the block and the RF's effective centre, 1/m.
        self.area_to_center_per_m = self.crush_cycles_per_voxel / (self.crush_voxel_mm / 1e3)
        self._min_crush_duration_s = float(pp.calc_duration(pp.make_trapezoid(
            channel=self.crush_axis, area=self.area_to_center_per_m, system=opts,
        )))
        self._crush_duration_s = self._resolve_crush_duration(crush_duration_s)

        raster = float(opts.grad_raster_time)
        center = float(pp.calc_rf_center(rf)[0])
        dead, ringdown = float(opts.rf_dead_time), float(opts.rf_ringdown_time)
        # The window the pulse occupies, **symmetrised about its effective centre** rather than
        # merely made long enough to cover it.  This is the time half of the invariant, and the
        # area solve below cannot buy it.
        self._plateau_s = 2 * ceil_raster(
            max(dead + center, self.duration_s - center + ringdown), raster,
        )
        self._bare_plateau_s = ceil_raster(dead + self.duration_s + ringdown, raster)
        self.gradients, self._plateau_start_s, self._duration_s = self._design_gradients()
        self.rf = derive(rf, delay=self._plateau_start_s + self._plateau_s / 2 - center)

    # ------------------------------------------------------------------ what it knows
    def time_to_center(self) -> float:
        """
        Seconds from the start of this module's block to the RF's effective centre.

        **The conjugation instant**, which is what makes it load-bearing rather than bookkeeping:
        pypulseq's ``calculate_kspacePP`` flips the sign of k here, so every echo time in a spin
        echo and every area balance in a train is measured from this number.

        From ``pp.calc_rf_center``, plus the pulse's own delay, which carries both the transmit
        dead time and the lead-in that symmetrises the plateau.  Identical in form to
        :meth:`Excitation.time_to_center` and :meth:`IRPrep.time_to_center`; a method with a
        different origin would break every timeline built by addition.
        """
        return float(self.rf.delay) + float(pp.calc_rf_center(self.rf)[0])

    def time_to_crusher(self) -> float:
        """
        Seconds to where the trailing crusher window begins.

        **The earliest another axis may start work, and where a caller drops the readout block.**
        It exists for exactly the reason :meth:`Excitation.time_to_rephaser` does, and it is the
        one number here worth real time: laying the readout after this module's block instead of at
        this instant costs a measured **1.4 ms per echo -- 13 %** of the minimum echo spacing, and
        30 extra pulseq blocks per eight-echo shot.

        By symmetry the *leading* crusher window is three axes wide too, and the axis that wants it
        is the readout in the first interval -- which is where the one-off dephaser goes.
        """
        return self._duration_s - self._crush_duration_s

    @property
    def area_from_center_per_m(self) -> float:
        """
        Area on `crush_axis` from the RF's effective centre to the end of the block, 1/m.

        Equal to :attr:`area_to_center_per_m` by construction, and **measured** rather than
        returned from the same expression: it integrates the waveform this module actually emits,
        so a test comparing the two compares a design against a result rather than a number with
        itself.
        """
        total = sum(
            area_until(event, float(pp.calc_duration(event)))
            for _, event in self.gradients if event.channel == self.crush_axis
        )
        return float(total) - self._area_on_crush_axis_until(self.time_to_center())

    @property
    def crush_duration_s(self) -> float:
        """The crusher window this module achieves, one side, seconds."""
        return self._crush_duration_s

    @property
    def min_crush_duration_s(self) -> float:
        """
        The shortest crusher window, seconds.  Unchanged by `crush_duration_s`.

        The shortest plain trapezoid carrying :attr:`area_to_center_per_m` on `crush_axis` -- exact
        when the crusher has an axis to itself, and conservative when it is fused with the
        selection plateau, since the plateau then carries part of that area.  One number that is
        legal either way, because `crush_axis` is an argument.
        """
        return self._min_crush_duration_s

    @property
    def plateau_padding_s(self) -> float:
        """
        What symmetrising the plateau about the effective centre cost, seconds.

        ``|rf_dead_time - rf_ringdown_time|`` for a symmetric pulse, and twice the pulse's own
        asymmetry for a minimum-phase one.  Reported rather than hidden, because it is the price of
        the time half of the invariant and a reader should be able to see what it was.
        """
        return self._plateau_s - self._bare_plateau_s

    # ----------------------------------------------------------------------- assembly
    def build(self, *, phase_deg: float = 0.0, position_mm: float = 0.0) -> LogicBlock:
        """
        Return the pulse and its gradients, as one block that starts and ends at zero.

        Parameters
        ----------
        phase_deg
            Carrier phase, degrees.  **The CPMG relation belongs to whatever holds both pulses**:
            a refocusing pulse does not know it is in a train.  Excite about ``+y`` and refocus
            about ``+x`` -- a refocusing pulse maps transverse phase to ``2*phase_ref - phase``, and
            90 degrees between the two axes is the fixed point of that map.
        position_mm
            Slice offset from isocentre along `axis`.  Requires a selective pulse; ``0.0`` is valid
            either way, and is what a non-selective pulse delivers.
        """
        rf = self.rf
        phase_rad = float(np.deg2rad(phase_deg))
        if phase_rad:
            rf = derive(rf, phase_offset=float(rf.phase_offset) + phase_rad)
        if position_mm:
            rf = shift_slice(rf, self._selection_gradient(), position_m=float(position_mm) / 1e3)

        out = LogicBlock().add(0.0, rf)
        for start, event in self.gradients:
            out.add(start, event)
        return out

    # -------------------------------------------------------------------------- design
    def _design_gradients(self) -> tuple[tuple[tuple[float, Event], ...], float, float]:
        """
        Return ``((start, event), ...)``, where the plateau starts, and how long the block is.

        Three shapes, one invariant.  When the crusher is on the selection axis the whole thing is
        **one** continuous waveform and the plateau's own half-area is part of what the crushers
        are solved against.  Otherwise the crushers are a pair of plain trapezoids on an axis of
        their own -- equal, so balanced by construction -- and the selection gradient, if there is
        one, is symmetric about the centre on its own.
        """
        t_c = self._crush_duration_s
        if self.selective and self.crush_axis == self.axis:
            fused = self._design_fused(t_c)
            return ((0.0, fused),), t_c, float(pp.calc_duration(fused))

        crusher = pp.make_trapezoid(
            channel=self.crush_axis, area=self.area_to_center_per_m, duration=t_c,
            system=self.opts,
        )
        if not self.selective:
            return (
                ((0.0, crusher), (t_c + self._plateau_s, crusher)),
                t_c, 2 * t_c + self._plateau_s,
            )
        # The selection gradient carries its own ramps here rather than the fused waveform's, so
        # the plateau starts one rise time in.  Equal ramps make its two halves equal without a
        # solve; the balance this module has to work for is the fused case's.
        gz = pp.make_trapezoid(
            channel=self.axis, amplitude=float(self.gz.amplitude), flat_time=self._plateau_s,
            system=self.opts,
        )
        gz_s = float(pp.calc_duration(gz))
        return (
            ((0.0, crusher), (t_c, gz), (t_c + gz_s, crusher)),
            t_c + float(gz.rise_time), 2 * t_c + gz_s,
        )

    def _design_fused(self, t_c: float) -> Event:
        """
        Return the one waveform: crusher, selection plateau, crusher, starting and ending at zero.

        The two crusher amplitudes come out of a fixed point rather than a formula.  Each depends on
        the ramps into and out of it and each ramp depends on the amplitude, and it converges
        because it is monotone: a longer ramp lowers the amplitude, which shortens the ramp it
        needs.  What the loop solves for is *the area to the effective centre*, plateau half
        included, which is the only definition of "balanced" the conjugation respects.
        """
        a_sel = float(self.gz.amplitude)
        t_p, target = self._plateau_s, self.area_to_center_per_m
        self._warn_if_the_plateau_already_carries_it(a_sel * t_p / 2)

        a_pre = a_post = target / t_c
        for _ in range(_MAX_SOLVER_STEPS):
            r_in, r_end = self._ramp(a_pre), self._ramp(a_post)
            r_mid, r_out = self._ramp(a_pre - a_sel), self._ramp(a_post - a_sel)
            # Of the ramp between the two plateaus, the part proportional to a_sel is fixed and the
            # rest scales with the amplitude being solved for -- hence the halves either side.
            new_pre = (target - a_sel * (t_p + r_mid) / 2) / (t_c - r_in / 2 - r_mid / 2)
            new_post = (target - a_sel * (t_p + r_out) / 2) / (t_c - r_end / 2 - r_out / 2)
            settled = (abs(new_pre - a_pre) < _SOLVER_TOL_HZ_M
                       and abs(new_post - a_post) < _SOLVER_TOL_HZ_M)
            a_pre, a_post = new_pre, new_post
            if settled:
                break
        else:                                                          # pragma: no cover
            msg = format_error(
                'the crusher amplitude did not converge.',
                {'crush_duration_s': t_c, 'crush_cycles_per_voxel': self.crush_cycles_per_voxel,
                 'rise_time_s': self.rise_time_s},
                ['pass a longer crush_duration_s', 'or rise_time_s, which fixes every ramp'],
            )
            raise RuntimeError(msg)

        r_in, r_end = self._ramp(a_pre), self._ramp(a_post)
        r_mid, r_out = self._ramp(a_pre - a_sel), self._ramp(a_post - a_sel)
        times = np.array([0.0, r_in, t_c - r_mid, t_c, t_c + t_p, t_c + t_p + r_out,
                          2 * t_c + t_p - r_end, 2 * t_c + t_p])
        amps = np.array([0.0, a_pre, a_pre, a_sel, a_sel, a_post, a_post, 0.0])
        # Limits are checked by sc.compile on the *summed* waveform, which is the only place the
        # truth is.  The ramps here are built from opts.max_slew, so a converged one sits exactly on
        # it and a second check in float would refuse this module's own arithmetic.
        return pp.make_extended_trapezoid(
            channel=self.crush_axis, times=times, amplitudes=amps, system=self.opts,
            max_grad=np.inf, max_slew=np.inf, skip_check=True,
        )

    def _ramp(self, delta: float) -> float:
        """The shortest legal ramp for an amplitude step, or the fixed one if one was asked for."""
        raster = float(self.opts.grad_raster_time)
        if self.rise_time_s is not None:
            return ceil_raster(self.rise_time_s, raster)
        return ceil_raster(max(abs(delta) / float(self.opts.max_slew), raster), raster)

    def _area_on_crush_axis_until(self, t_end: float) -> float:
        """The area on `crush_axis` from the start of the block to `t_end`, 1/m."""
        return float(sum(
            area_until(event, t_end - start)
            for start, event in self.gradients if event.channel == self.crush_axis
        ))

    def _resolve_crush_duration(self, requested_s: float | None) -> float:
        """Return the crusher window, refusing one shorter than the area it has to carry."""
        if requested_s is None:
            return self._min_crush_duration_s
        wanted = ceil_raster(
            require_positive(requested_s, 'crush_duration_s'), self.opts.grad_raster_time,
        )
        # EPS, not exact: a composite passing back `max(refoc.min_crush_duration_s, ...)` hands us
        # our own minimum, and snapping that onto the raster can land one ulp below where it came
        # from.  Times a nanosecond apart are the same time.
        if wanted < self._min_crush_duration_s - EPS:
            msg = format_error(
                f'crush_duration_s = {requested_s * 1e6:.1f} us is shorter than the minimum this '
                f'crusher needs.',
                {'crush_duration_s': requested_s,
                 'min_crush_duration_s': self._min_crush_duration_s,
                 'area_to_center_per_m': self.area_to_center_per_m},
                [
                    f'pass crush_duration_s >= {self._min_crush_duration_s:.6g}',
                    'or None for the shortest legal crusher',
                    'a composite passes the maximum over its three axes here',
                ],
            )
            raise ConfigurationError(msg)
        return wanted

    # -------------------------------------------------------------------- the refusals
    def _check_b1(self, rf: Event) -> None:
        """
        Refuse a pulse over ``opts.max_b1``, naming the duration that fixes it.

        ``sc.compile`` checks gradient amplitude and slew and never looks at RF amplitude, and
        pypulseq warns and hands the pulse back.  A 180 is where this first bites, because peak B1
        scales with flip angle at a fixed shape.  It scales as ``1 / duration`` too, which is what
        makes the fix quotable rather than a search.
        """
        max_b1 = float(getattr(self.opts, 'max_b1', 0.0) or 0.0)
        peak = float(np.abs(rf.signal).max())
        if max_b1 <= 0.0 or peak <= max_b1:
            return
        # Rounded up onto a tenth of a millisecond: the relation is exact for a fixed shape, and a
        # floor quoted to the nanosecond would be refused again by its own last digit.
        floor_ms = float(np.ceil(self.duration_s * peak / max_b1 * 1e4) / 10.0)
        msg = format_error(
            f'a {self.duration_s * 1e3:g} ms {self.flip_deg:g} degree {self.pulse} pulse peaks at '
            f'{peak / max_b1 * 100:.0f} % of max_b1.',
            {'duration_s': self.duration_s, 'flip_deg': self.flip_deg,
             'peak_b1_hz': peak, 'max_b1': max_b1},
            [
                f'pass duration_s >= {floor_ms:.1f} ms, which is where this shape fits',
                'or lower flip_deg: peak B1 scales with the flip angle at a fixed shape',
                'sc.compile checks gradients and not RF amplitude, so this is refused here',
            ],
        )
        raise ConfigurationError(msg)

    def _warn_if_the_plateau_already_carries_it(self, plateau_half: float) -> None:
        """
        Warn when the leading crusher has to change sign to reach the requested area.

        Legal, and occasionally what somebody means -- but the selection plateau's own half-area is
        about ``time_bw_product / 2`` cycles across the slice, so a smaller
        `crush_cycles_per_voxel` produces a *negative* leading lobe, which is no longer anything a
        reader would recognise as a crusher.
        """
        if plateau_half <= self.area_to_center_per_m:
            return
        warnings.warn(
            format_error(
                f'crush_cycles_per_voxel = {self.crush_cycles_per_voxel:g} is below the '
                f'{plateau_half * self.crush_voxel_mm / 1e3:.2f} cycles the selection plateau '
                f'already carries, so the leading crusher changes sign.',
                {'crush_cycles_per_voxel': self.crush_cycles_per_voxel,
                 'area_to_center_per_m': self.area_to_center_per_m,
                 'plateau_half_area_per_m': plateau_half},
                ['the balance still holds; the lobe is simply no longer a crusher',
                 'raise crush_cycles_per_voxel above time_bw_product / 2 to keep it one'],
            ),
            SeqCraftWarning,
            stacklevel=4,
        )

    def _selection_gradient(self) -> Event:
        """Return the selection gradient, or refuse an offset there is nothing to impose."""
        if not self.selective:
            msg = format_error(
                'position_mm needs a selection gradient, and this pulse is non-selective.',
                {'thickness_mm': None, 'pulse': self.pulse},
                [
                    'pass thickness_mm to make the pulse slice-selective',
                    'or leave position_mm at 0.0, which is what a non-selective pulse delivers',
                ],
            )
            raise ConfigurationError(msg)
        return self.gz

    def _check_pulse(self, pulse: str) -> str:
        """Return `pulse` having checked it names one of the four factories."""
        if pulse not in _FACTORIES:
            listed = ', '.join(repr(name) for name in _FACTORIES)
            msg = format_error(
                f'pulse must be one of {listed}, got {pulse!r}.', {'pulse': pulse},
            )
            raise ConfigurationError(msg)
        return pulse

    def _check_axis(self, axis: str | None) -> str | None:
        """Return the selection axis: ``'z'`` by default when selective, ``None`` when not."""
        if not self.selective:
            if axis is not None:
                msg = format_error(
                    f'axis={axis!r} was passed with thickness_mm=None, so it cannot take effect.',
                    {'axis': axis, 'thickness_mm': None},
                    [
                        'pass thickness_mm to select a slice along axis',
                        'or drop axis: a non-selective pulse has no selection axis',
                        'crush_axis is a different argument and applies either way',
                    ],
                )
                raise ConfigurationError(msg)
            return None
        return require_axis('z' if axis is None else axis)

    def _check_crush_voxel(self, crush_voxel_mm: float | None) -> float:
        """Return the length the crusher's cycles are counted across, refusing to invent one."""
        if crush_voxel_mm is not None:
            return require_positive(crush_voxel_mm, 'crush_voxel_mm')
        if self.thickness_mm is not None:
            return self.thickness_mm
        msg = format_error(
            'crush_voxel_mm is needed: a non-selective refocusing pulse has no thickness to count '
            'crusher cycles across.',
            {'thickness_mm': None, 'crush_axis': self.crush_axis,
             'crush_cycles_per_voxel': self.crush_cycles_per_voxel},
            [
                'pass crush_voxel_mm = the imaging voxel dimension along crush_axis',
                'or pass thickness_mm for a slice-selective pulse, which then supplies it',
            ],
        )
        raise ConfigurationError(msg)

    def _check_combination(self, time_bw_product: float | None) -> None:
        """Refuse the two `pulse`-dependent combinations that cannot mean anything."""
        if self.pulse != 'block':
            return
        if self.selective:
            msg = format_error(
                "pulse='block' cannot be slice-selective.",
                {'pulse': 'block', 'thickness_mm': self.thickness_mm},
                [
                    'pass thickness_mm=None for a hard refocusing pulse',
                    "or pulse='sinc' / 'slr' / 'gauss' to select a slice",
                    "a rectangular envelope's slice profile is a sinc with sidelobes, so "
                    'pp.make_block_pulse offers neither slice_thickness nor return_gz',
                ],
            )
            raise ConfigurationError(msg)
        if time_bw_product is not None:
            msg = format_error(
                "time_bw_product has no meaning for pulse='block' at a fixed duration.",
                {'pulse': 'block', 'time_bw_product': time_bw_product},
                ['drop time_bw_product', "or pulse='sinc', whose bandwidth it sets"],
            )
            raise ConfigurationError(msg)

    def _check_pulse_opts(self, pulse_opts: dict[str, Any] | None) -> dict[str, Any]:
        """Return `pulse_opts` having checked every key against the chosen factory."""
        if not pulse_opts:
            return {}
        import inspect  # noqa: PLC0415  (only reached when an escape hatch is actually used)

        factory = getattr(pp, _FACTORIES[self.pulse])
        accepted = set(inspect.signature(factory).parameters) - _RESERVED
        unknown = sorted(set(pulse_opts) - accepted)
        if unknown:
            msg = format_error(
                f'pulse_opts key(s) {", ".join(repr(k) for k in unknown)} are not accepted by '
                f'pulse={self.pulse!r}.',
                {'pulse': self.pulse, 'accepted': ', '.join(sorted(accepted))},
                ['pulse_opts forwards pulse-design parameters only; this module owns the rest'],
            )
            raise ConfigurationError(msg)
        return dict(pulse_opts)
