"""
:class:`CartesianLine` -- prephaser, readout gradient and ADC designed as one thing.

One line of k-space: a prephaser that moves `k` to the start of the line, a readout gradient that
traverses it, and an ADC that samples while it does.

The **line index is not an argument**.  A Cartesian line is identical every TR, and which line of
k-space it lands on is set by the phase encoding on another axis.  What this module does own is
where the echo falls inside the readout, which nothing downstream can recover -- a block of events
knows times, not meanings.

The one piece of arithmetic worth reading twice
-----------------------------------------------
The readout gradient accumulates moment from the instant it starts *ramping*, so k at the echo is
the prephaser plus the ramp-up area plus whatever flat top precedes the echo::

    A_pre = -( 1/2 * G * t_rise  +  G * t_flat_before_echo )

For the symmetric case -- ADC over the flat top, equal rise and fall -- that collapses to a
formula which differs from the classic error by a single word::

    prephaser_area = -gx.flat_area / 2      # WRONG: the ramp-up is missing
    prephaser_area = -gx.area / 2           # nearly right: total area, ramps included

``pp.make_trapezoid`` supplies both fields, one letter apart in spelling and a whole ramp apart in
meaning.  Getting it wrong offsets the whole readout in k by the ramp area: the echo stops
falling at the ADC centre, and the image picks up a linear phase ramp -- or, if the shift is
large, ringing and signal loss that read as a hardware fault.

**The second one is still not exact**, and the residue is worth naming rather than absorbing:
k = 0 belongs at sample index ``matrix // 2``, whose centre time is *half a dwell* after the
midpoint of the ADC window.  That is the convention a centred inverse FFT assumes, and it is what
puts the sampled positions on ``-N/2 ... N/2-1`` times ``dk`` instead of on a half-integer grid.
So the exact answer is ``-gx.area / 2 - dk / 2``, and a formula derived from the trapezoid's
geometry alone cannot know that.

This module computes no closed form at all.  It integrates the gradient's own knots up to the
echo (:func:`~seqcraft.modules._support.area_until`), which produces the half-dwell term for
free and generalises without a second derivation to the two cases where even the corrected
closed form stops holding: **ramp sampling**, where the ADC opens during the ramp, and **partial
echo**, where the echo is not at the centre of the readout.  Hand-deriving each of those
separately is how the ramp gets dropped from one of them.

``prephase=False``, for a spin-echo readout
-------------------------------------------
Dropping the prephaser event and changing nothing else is what a **spin-echo** readout needs.  In a CPMG train the readout's dephasing is not cancelled by an event of its own: it is half
of a pair that straddles a *refocusing pulse*, and the refocusing pulse's conjugation does the
cancelling.  So there is no prephaser to design -- but everything else here is identical, down to
the last microsecond: the dwell resolution and its raster snap, ``num_samples`` and
``pre_echo_samples`` from `partial_fourier`, the amplitude ``dk / dwell``, the flat top rounded onto
the gradient raster, the ADC-dead-time overrun fix, and the exact pre-echo area.

:attr:`~CartesianLine.area_to_echo_per_m` is what makes that work.  It is the *physics* number --
what the readout accumulates by the echo -- rather than the *event's* number, so it is what a
prephaser cancels when there is one, and what a symmetric crusher pair is balanced around when
there is not.

``echoes`` and ``polarity``, for a multi-echo train
---------------------------------------------------
Reading the same line more than once after one excitation is a multi-echo gradient echo.  It is
still one line of k-space, read repeatedly, which is why ``ECO`` varies across the block and
``LIN`` does not.

Three things change with a second echo, and each of them is a place where the obvious
implementation is wrong in a way that compiles, passes every k-space extent check, and produces a
plausible image.

**1.  A reverse lobe needs a different lobe.**  The flat top this module rounds *up* onto the
gradient raster leaves its slack after the last sample, which changes no k position -- correct for
one line, and correct for a monopolar train, where every lobe is that same lobe.  It is wrong the
moment a reverse lobe has to sample the same grid, because the window is then off-centre by half
the rounding.  Under ``'bipolar'`` the lobe is redesigned so that::

    T = 2*guard + num_samples*dwell         # exactly, not "at least"

Measured at 220 mm, 128, 250 Hz/px: today's lobe has a 20.00 us lead guard and a 26.40 us trailing
one, and a two-echo bipolar train built out of it samples two grids **0.9324 1/m = 0.2051 dk**
apart.  With the centred lobe the same measurement is **1.7e-13 1/m**.  A fifth of a dk is not a
rounding error; it is a linear phase ramp between odd and even echoes, which is the exact thing a
field map measures -- so it never appears as an artefact, only as a wrong number in the map.

**2.  The fly-back cancels the whole lobe.**  ``'monopolar'`` returns k to where the prephaser left
it between echoes, so lobe *n+1* traverses the trajectory lobe *n* did.  That area is the lobe's
**total** area, and the three plausible wrong answers all compile::

    flyback_area = -gx.flat_area                # WRONG -- both ramps missing
    flyback_area = -gx.area / 2                 # WRONG -- that is the prephaser's job
    flyback_area = -self.area_to_echo_per_m     # WRONG -- the PRE-echo part only

The third is the dangerous one, because it is a number this module already computes and already
exposes.  Measured at the reference protocol: the lobe's total area is **585.664336 1/m** and
``area_to_echo_per_m`` is **294.638695 1/m** -- the same order of magnitude, roughly a factor of
two apart, and each echo would start half a k-space further along than the last.

**3.  Each seam carries a barrier.**  A Pulseq block holds at most one gradient per axis and one
ADC, so a block boundary has to fall somewhere between consecutive echoes.  ``'bipolar'`` emits
one :func:`~seqcraft.barrier` at each seam and ``'monopolar'`` two, at the fly-back's start and at
its end -- one per gradient edge -- so the train's block structure is part of what this module
emits rather than something a caller has to arrange.

And the fourth thing, which is a fact about the caller rather than about the waveform:
**the echo times are not the echo spacing**.  ``k = 0`` is a *sample*, not an instant, and it sits
at ``pre_echo_samples`` on a forward lobe and at ``num_samples - 1 - pre_echo_samples`` on a
reverse one -- two different distances into two identical lobes.  So a bipolar train's echo times
alternate about the seam-to-seam period by one dwell either way, measured **32.00 us** at 500 Hz/px
and **1953.00 us** at ``partial_fourier=0.75``.  :attr:`~CartesianLine.te_s` is the sanctioned
source of echo times and :attr:`~CartesianLine.echo_spacing_s` is documented as the period and
nothing else, because a two-point field map divided by the period instead of by
``te_s[1] - te_s[0]`` is a **0.75 % scale error in every voxel** -- and a scale error in a field
map is exactly the kind of wrong no image inspection finds.
"""

from __future__ import annotations

from math import pi
from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ...design.events import derive, knots_of, pwl_moment
from ...design.logic import LogicBlock, barrier
from ...design.module import Module
from ...design.timing import EPS, from_ticks, to_ticks
from ...errors import ConfigurationError, format_error
from .._support import (
    area_until,
    ceil_raster,
    dwell_quantum,
    halve_onto,
    require_axis,
    require_count,
    require_positive,
    require_range,
)

if TYPE_CHECKING:
    from pypulseq.opts import Opts

    from ...design.events import Event

__all__ = ['CartesianLine']

#: What each polarity costs, in one line each.  Quoted in the refusal, because the two produce
#: files that differ in the dwell, the echo times, the k ordering of every second echo and the
#: block count -- and look identical in a protocol printout.
_POLARITIES = {
    'monopolar': 'every lobe the same sign, with a fly-back between them: one k grid and '
                 'uniformly spaced echo times, at ~25% of the period spent flying back',
    'bipolar': 'alternating sign and no fly-back: ~12% shorter period and half the blocks, at '
               'a reversed grid on every second echo and echo times that alternate by one dwell',
}


class CartesianLine(Module):
    """
    One Cartesian readout: prephaser, readout gradient, and the ADC that samples it.

    Parameters
    ----------
    opts
        The scanner.
    fov_mm
        Field of view along `axis`.  Sets the k-space step, ``dk = 1000 / fov_mm`` in 1/m.
    matrix
        Samples across the **full** k-space extent.  With `partial_fourier` below 1 the ADC takes
        fewer than this; the extent one sample spans is unchanged, which is what keeps the image
        resolution attached to `matrix` rather than to how much of k-space was sampled.
    bandwidth_hz_px, dwell_s
        The sampling rate, either way round.  **Exactly one is required**; both or neither raises.
        ``dwell = 1 / (bandwidth_hz_px * matrix)``, then snapped onto the ADC raster -- or, under
        ``polarity='bipolar'``, up onto the coarser quantum that keeps the guard on the RF raster
        -- so the achieved values are readable back as :attr:`dwell_s` and :attr:`bandwidth_hz_px`.
    partial_fourier
        Fraction of the full k-space extent that is sampled, with the missing part taken off the
        **pre-echo** side.  ``1.0`` is a symmetric full echo.  ``0 < partial_fourier <= 1``.
    echoes
        How many times the line is read after one excitation.  ``1`` is a single-echo readout and
        is byte-identical to this module before the train existed.  ``>= 1``.
    polarity
        ``'monopolar'`` -- every lobe the same sign, with a fly-back between them.  ``'bipolar'``
        -- alternating sign, no fly-back.  **No default, and required when `echoes` is above 1;**
        passing it with ``echoes=1`` raises.  There is no default for the same reason
        :class:`~seqcraft.modules.Excitation` has none for ``thickness_mm``: it is the design's
        most consequential property, a reconstruction has to know which one it got, and neither
        answer is safe to assume.
    echo_spacing_s
        Requested **seam-to-seam** period.  ``None`` is :attr:`min_echo_spacing_s`; below it
        raises, naming the minimum.  Raises with ``echoes=1``.  Under ``'monopolar'`` the extra
        time lengthens the *fly-back*, never the lobe; under ``'bipolar'`` it lengthens the lobe
        symmetrically, which keeps the sampling window centred in it.
    prephase
        ``False`` drops the prephaser event, leaving the readout gradient and its ADC starting at
        zero.  That is the **spin-echo** readout: its dephasing is half of a pair straddling a
        refocusing pulse, so there is no event of its own to cancel it and the caller places the
        two lobes.  Nothing else moves -- see :attr:`area_to_echo_per_m`.  ``prephase=False``
        works with a **monopolar** train; the bipolar combination raises, because balancing
        reversed lobes across a refocusing pulse is a question nobody has measured yet.
    prephaser_duration_s
        Lengthen the prephaser beyond its own minimum.  ``None`` is the minimum.  A composite
        passes the winder maximum here; see :class:`~seqcraft.modules.GRE2DTR`.  Passing it with
        ``prephase=False`` raises, in the shape this library already uses for an argument that
        cannot take effect.  At ``null_moment_order=1`` it is the **total** across the two
        winder lobes, and the areas are re-solved against it rather than scaled.
    null_moment_order
        The highest gradient moment nulled at the echo on this axis.

        ``0``
            The ordinary prephaser: one lobe, ``m0 = 0`` at the echo, which is what puts
            ``k = 0`` there.  **The default**, and unchanged behaviour.
        ``1``
            Velocity compensation: two winder lobes of opposite sign instead of one, nulling
            ``m0`` **and** ``m1`` at the echo, so a spin moving at constant velocity along `axis`
            arrives with the phase it would have had standing still.

        Needs ``prephase=True``, because the winder is the waveform being reshaped.  It costs
        time -- three lobes before the echo rather than two -- and :attr:`prephaser_duration_s`
        and :meth:`time_to_echo` report what it came to.

        The moment is nulled at the **first** echo.  Later echoes of a train accumulate their own
        first moment from the lobes between them, which this does not compensate.  Nulling the
        second moment as well needs a fourth lobe and is not implemented.
    axis
        Logical gradient channel.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
    ...          adc_dead_time=10e-6)
    >>> ro = CartesianLine(opts=o, fov_mm=250.0, matrix=128, bandwidth_hz_px=250.0)
    >>> ro.num_samples, round(ro.dwell_s * 1e9)
    (128, 31200)

    The prephaser cancels everything the readout has accumulated by the echo, ramp included:

    >>> gx = ro.gx
    >>> from seqcraft.modules._support import area_until
    >>> accumulated = area_until(gx, ro.time_to_echo() - ro.prephaser_duration_s)
    >>> abs(ro.prephaser_area_per_m + accumulated) < 1e-9
    True

    Partial echo moves the echo forward and shortens the ADC, in one factor with no new branch:

    >>> pf = CartesianLine(opts=o, fov_mm=250.0, matrix=128, bandwidth_hz_px=250.0,
    ...                    partial_fourier=0.75)
    >>> pf.num_samples
    96
    >>> round(pf.pre_echo_samples / pf.num_samples, 3)          # (2*pf - 1) / (2*pf) = 1/3
    0.333

    ``prephase=False`` changes no designed number -- only which events come out:

    >>> se = CartesianLine(opts=o, fov_mm=250.0, matrix=128, bandwidth_hz_px=250.0,
    ...                    prephase=False)
    >>> (se.dwell_s, se.num_samples, se.pre_echo_samples) == (
    ...     ro.dwell_s, ro.num_samples, ro.pre_echo_samples)
    True
    >>> se()                                    # the readout gradient and its ADC, and that is all
    LogicBlock(CartesianLine, 2 nodes, 4.06 ms)
    >>> abs(se.area_to_echo_per_m + ro.prephaser_area_per_m) < 1e-9
    True

    Eight echoes off one excitation, both ways round.  Monopolar keeps the single-echo lobe
    exactly and pays a fly-back; bipolar buys that time back and pays a reversed grid:

    >>> mono = CartesianLine(opts=o, fov_mm=250.0, matrix=128, bandwidth_hz_px=250.0,
    ...                      echoes=8, polarity='monopolar')
    >>> from seqcraft.design.events import content_hash
    >>> content_hash(mono.gx) == content_hash(ro.gx)            # the same lobe, byte for byte
    True
    >>> bi = CartesianLine(opts=o, fov_mm=250.0, matrix=128, bandwidth_hz_px=250.0,
    ...                    echoes=8, polarity='bipolar')
    >>> bi.echo_spacing_s < mono.echo_spacing_s
    True
    >>> [bi.polarity_of(n) for n in range(4)]
    [1, -1, 1, -1]

    And the echo times are **not** the spacing, which is what :attr:`te_s` exists to say:

    >>> import numpy as np
    >>> round(float(np.ptp(np.diff(mono.te_s))) * 1e9)          # uniform, to the nanosecond
    0
    >>> round(float(np.ptp(np.diff(bi.te_s))) / bi.dwell_s, 6)  # alternating, by two dwells
    2.0

    Notes
    -----
    **`acquire=False` plays the same gradients with no ADC**, for dummy repetitions.  The
    duration is unchanged, which is the point: a dummy has to load the gradients exactly as a real
    TR does, or the steady state it is establishing is not the one that will be acquired.  It
    suppresses the ``ECO`` and ``REV`` labels with it, for the reason
    :class:`~seqcraft.modules.GRE2DTR` suppresses ``LIN``: a dummy that emitted an echo index
    would put one in the file for a readout that was never sampled.

    **`ECO` and `REV` are emitted on every acquired echo of a train, unconditionally**, including
    echo 0 and including ``'monopolar'`` -- and *not at all* at ``echoes=1``, which is what keeps
    that case byte-identical.  Pulseq labels are stateful: ``SET`` writes a value that persists
    until something changes it, so skipping ``ECO`` on echo 0 as "the default" leaves it wearing
    the previous repetition's last value, and omitting ``REV`` under ``'monopolar'`` inherits
    whatever an EPI train, a bipolar train or a navigator left standing.  A monopolar train that
    inherits ``REV = 1`` is mirrored by the reconstruction, and against a symmetric phantom that
    looks entirely correct.  ``SLC``, ``SEG``, ``REP`` and ``AVG`` are still not this module's
    business.

    **There is deliberately no `k_read_per_m` array.**  :class:`~seqcraft.modules.EPI2D` exposes
    one because ramp sampling makes the trajectory underivable; here the ADC never opens on a
    ramp, so sample *i* of echo *e* is at ``(i - echo_sample(e)) * dk_per_m * polarity_of(e)`` and
    an array would be a second copy of two methods that already exist.

    **`phase_deg` is the receiver phase, and it is not optional in a spoiled sequence.**  The
    receiver is phase-locked to the transmitter: whatever carrier phase the excitation was given,
    the ADC has to demodulate at.  Leaving it at zero while the RF runs a spoiling schedule
    writes the schedule's quadratic phase straight into ``ky`` -- measured on a single off-centre
    voxel, that scattered a point source over the whole phase-encode direction and moved its peak
    by thirteen pixels, with the readout direction perfectly correct beside it.  The composite
    that owns the schedule passes the same number to both; this module only applies it.  One
    receiver phase covers a whole train, because one train is one excitation.

    **`offset_mm` shifts the FOV along the readout axis** by demodulating at a shifted frequency,
    ``adc.freq_offset = gx.amplitude * offset_m``, with the phase referenced to the echo so the
    shift contributes no phase at k = 0.  It is the receive-side twin of the slice offset an
    :class:`~seqcraft.modules.Excitation` applies to its RF, and the sign is a trap in both:
    against a symmetric phantom a mirrored FOV looks entirely correct.  Note that MRzero's
    ``Sequence.import_file`` **parses ``freq`` and then ignores it**, so a simulation through a
    written ``.seq`` cannot see this offset at all -- checking it needs a scanner or a simulator
    that demodulates.
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
        echoes: int = 1,
        polarity: str | None = None,
        echo_spacing_s: float | None = None,
        prephase: bool = True,
        prephaser_duration_s: float | None = None,
        null_moment_order: int = 0,
        axis: str = 'x',
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.fov_mm = require_positive(fov_mm, 'fov_mm')
        self.matrix = int(matrix)
        self.axis = require_axis(axis)
        self.prephase = bool(prephase)
        self.partial_fourier = require_range(partial_fourier, 'partial_fourier', low=0.0, high=1.0)
        if self.matrix < 1:
            msg = format_error(f'matrix must be at least 1, got {self.matrix}.',
                               {'matrix': self.matrix})
            raise ConfigurationError(msg)

        #: How many times the line is read after one excitation.
        self.echoes = require_count(echoes, 'echoes', hint='1 is a single-echo readout')
        #: ``'monopolar'``, ``'bipolar'``, or ``None`` at :attr:`echoes` of 1.
        self.polarity = self._resolve_polarity(polarity)
        self._refuse_bipolar_spin_echo()

        #: k-space step between adjacent samples, 1/m.
        self.dk_per_m = 1e3 / self.fov_mm
        # The sample count comes *before* the dwell, because under 'bipolar' it is what the dwell
        # is quantised against -- see `_resolve_dwell`.  It does not depend on the dwell either
        # way round, so the order costs the single-echo path nothing.
        self.num_samples, self.pre_echo_samples = self._resolve_samples()
        self.dwell_s = self._resolve_dwell(bandwidth_hz_px, dwell_s)
        #: Achieved receive bandwidth per pixel, Hz.  Reads back the snapped dwell.
        self.bandwidth_hz_px = 1.0 / (self.dwell_s * self.matrix)

        self._min_period_s: float | None = None
        self._period_s: float | None = None
        self._flyback: Event | None = None
        self.gx, self.adc = (
            self._design_lobe_centred(echo_spacing_s) if self.polarity == 'bipolar'
            else self._design_readout()
        )
        #: Seconds from the readout gradient's own start to k = 0.
        self._echo_in_gx = float(self.adc.delay) + (self.pre_echo_samples + 0.5) * self.dwell_s
        self._lobe_s = float(pp.calc_duration(self.gx))
        self._gx_reverse = (
            derive(self.gx, amplitude=-float(self.gx.amplitude))
            if self.polarity == 'bipolar' else None
        )
        if self.polarity == 'monopolar':
            self._flyback = self._design_flyback(echo_spacing_s)
        self._refuse_echo_spacing(echo_spacing_s)
        #: The highest gradient moment nulled at the echo on this axis.
        self.null_moment_order = self._check_null_moment_order(null_moment_order)
        self._prephasers: tuple[Event, ...] = (
            self._design_prephaser(prephaser_duration_s) if self.prephase
            else self._refuse_prephaser_duration(prephaser_duration_s)
        )

    # ------------------------------------------------------------------ what it knows
    @property
    def area_to_echo_per_m(self) -> float:
        """
        What the **readout gradient** has accumulated by the echo, 1/m.  Always available.

        The physics number rather than an event's: it is what the prephaser cancels when there is
        one (``prephaser_area_per_m == -area_to_echo_per_m``), and it is the number a symmetric
        crusher pair straddling a refocusing pulse is balanced around when there is not.  Exposing
        it is what lets one module serve a gradient echo and a spin echo.

        Integrated from the readout's own knots, so partial Fourier, ramp sampling and the
        half-dwell sample offset all come out of it without a second derivation.

        In a train this describes the **first** lobe, which is the one the prephaser cancels and
        the one a crusher pair is balanced around.  It is emphatically *not* the fly-back's area;
        see :attr:`flyback_area_per_m`.
        """
        return area_until(self.gx, self._echo_in_gx)

    @property
    def area_after_echo_per_m(self) -> float:
        """
        What the **readout gradient** still has to play after the echo, 1/m.  Always available.

        The other half of :attr:`area_to_echo_per_m`, and exposed for the same reason: it is an
        intrinsic fact about this readout's own geometry, and a caller that needs it should not
        have to re-derive it from ``gx.area`` and a sign convention.

        The lobe's own area less what precedes the echo, so that the two sum to the lobe exactly
        rather than to within a rounding: partial Fourier, the half-dwell sample offset and the
        slack left after the last sample all come out of :attr:`area_to_echo_per_m`, and this is
        the remainder.

        In a train this describes the **first** lobe, exactly as :attr:`area_to_echo_per_m` does.
        """
        return float(self.gx.area) - self.area_to_echo_per_m

    @property
    def echo_moment_imbalance_per_m(self) -> float:
        """
        How asymmetric this readout is about its echo, 1/m: the area after minus the area before.

        Zero would mean the echo sits exactly halfway through the gradient's area.  It never quite
        does -- the flat top is rounded up onto the gradient raster with the slack left after the
        last sample, the echo sample sits half a dwell off the window's centre, and under partial
        Fourier it is nowhere near the middle -- so this is a small number that is not zero, and
        its sign matters.

        **What this module does not do is compensate for it.**  A spin-echo caller straddles the
        readout with a lobe either side of a refocusing pulse and needs the area before the echo
        to equal the area after it; the two lobes therefore differ by this imbalance.  But *how*
        to split it between them depends on the crusher window, the refocusing geometry and the
        echo spacing, none of which this module can see.  The geometry is owned here; the
        compensation belongs to whatever holds both components --
        :class:`~seqcraft.modules.TSEShot` takes half of this as the difference between its two
        readout-axis lobes.
        """
        return self.area_after_echo_per_m - self.area_to_echo_per_m

    @property
    def prephaser_duration_s(self) -> float:
        """
        Seconds the prephaser occupies -- its own minimum unless one was requested.

        At ``null_moment_order=1`` the prephaser is two lobes played back to back and this is the
        pair's total, which is what the readout is placed after.
        """
        return float(sum(pp.calc_duration(lobe) for lobe in self._require_prephasers()))

    @property
    def prephaser_area_per_m(self) -> float:
        """
        The prephaser's area, 1/m -- exactly minus :attr:`area_to_echo_per_m`.

        The **total** across its lobes, which is what puts ``k = 0`` at the echo.  At
        ``null_moment_order=1`` the two lobes have opposite signs and individually are neither of
        these; :attr:`prephaser_lobes` is where to read them.

        Exposed so a test can assert that identity without reaching into the block.
        """
        return float(sum(float(lobe.area) for lobe in self._require_prephasers()))

    @property
    def prephaser_lobes(self) -> tuple[Event, ...]:
        """
        The prephaser's own gradient events, in the order they play.

        One at ``null_moment_order=0``.  Two at ``1``, of opposite sign, whose areas sum to
        :attr:`prephaser_area_per_m` and whose first moments about the echo cancel the readout's.

        A leaf exposing its events, which is what lets a composing module measure what it needs
        over the interval it cares about rather than trusting a cached number.
        """
        return self._require_prephasers()

    @property
    def echo_spacing_s(self) -> float:
        """
        The **seam-to-seam period** of the train, seconds -- identical for every echo.

        Bipolar: the lobe.  Monopolar: the lobe plus the fly-back.  **This is not the spacing
        between echo times**, and the difference is the whole reason :attr:`te_s` exists: under
        ``'bipolar'`` the echo times alternate about this number by one dwell either way, and
        under ``partial_fourier`` below 1 they alternate by nearly the period itself.  A
        reconstruction that divides a phase difference by this rather than by ``te_s[1] -
        te_s[0]`` returns a **scaled** field map, in every voxel, with nothing to look wrong.
        """
        return self._require_train('echo_spacing_s')

    @property
    def min_echo_spacing_s(self) -> float:
        """The shortest legal period, seconds.  A feasibility fact known at design time."""
        self._require_train('min_echo_spacing_s')
        # Set by whichever design method ran, in the same breath as `_period_s` -- so the refusal
        # above is what makes this reachable only when there is a number to return.
        return float(self._min_period_s or 0.0)

    @property
    def te_s(self) -> tuple[float, ...]:
        """
        Seconds from the block's start to ``k = 0``, **per echo**.  ``len == echoes``.

        The array a T2\\* or B0 fit is run against, and the **only sanctioned source of echo
        times**.  ``k = 0`` is a *sample* and not an instant -- index :attr:`pre_echo_samples` on
        a forward lobe, ``num_samples - 1 - pre_echo_samples`` on a reverse one -- so a bipolar
        train's echo times are not uniformly spaced even though its lobes are.  Measured over 8
        echoes at 500 Hz/px: monopolar spreads 0.00 us and bipolar 32.00 us, which is two dwells;
        at ``partial_fourier=0.75`` bipolar spreads **1953.00 us**, comparable to the spacing
        itself.

        A ``uniform_te`` flag was considered and rejected: it invites a branch at the call site
        whose false arm is the one that is wrong, where reading this is correct in both.
        """
        start = self.prephaser_duration_s if self.prephase else 0.0
        period = self._period_s or 0.0
        return tuple(
            start + n * period + self._echo_in_lobe_s(n) for n in range(self.echoes)
        )

    def time_to_echo(self, echo: int = 0) -> float:
        """
        Seconds from the start of this module's block to ``k = 0`` of `echo`.

        The no-argument call is unchanged and still means the **first** echo, which is what TE
        means in a multi-echo protocol; the rest of the train is :attr:`te_s`.

        Measured from the readout gradient's own start when ``prephase=False``, because that is
        then where the block begins.  Not measurable from the tree either way: the block knows when
        its events play, not which instant among them is the echo.
        """
        return self.te_s[self._require_echo(echo)]

    def polarity_of(self, echo: int) -> int:
        """
        ``+1`` or ``-1``: the sign of `echo`'s readout lobe.  Always ``+1`` under ``'monopolar'``.

        Examples
        --------
        >>> import pypulseq as pp
        >>> from pypulseq.opts import Opts
        >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s',
        ...          adc_dead_time=10e-6)
        >>> line = CartesianLine(opts=o, fov_mm=220.0, matrix=128, bandwidth_hz_px=500.0,
        ...                      echoes=4, polarity='bipolar')
        >>> [line.polarity_of(n) for n in range(4)]
        [1, -1, 1, -1]
        """
        return -1 if self.polarity == 'bipolar' and self._require_echo(echo) % 2 else 1

    def echo_sample(self, echo: int) -> int:
        """
        Which ADC sample of `echo` carries ``k = 0``.

        ``pre_echo_samples`` forward, ``num_samples - 1 - pre_echo_samples`` reverse.  It exists
        because every caller otherwise writes that expression by hand, and the wrong branch is not
        an error -- it is a **mirrored echo**, which against a symmetric phantom looks correct.

        Examples
        --------
        >>> import pypulseq as pp
        >>> from pypulseq.opts import Opts
        >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s',
        ...          adc_dead_time=10e-6)
        >>> line = CartesianLine(opts=o, fov_mm=220.0, matrix=128, bandwidth_hz_px=500.0,
        ...                      partial_fourier=0.75, echoes=4, polarity='bipolar')
        >>> [line.echo_sample(n) for n in range(4)]
        [32, 63, 32, 63]
        """
        if self.polarity_of(echo) > 0:
            return self.pre_echo_samples
        return self.num_samples - 1 - self.pre_echo_samples

    @property
    def flyback_area_per_m(self) -> float:
        """
        The fly-back's area, 1/m -- exactly minus the **whole** lobe's area.  Monopolar only.

        Not ``-area_to_echo_per_m``, which is the pre-echo part and is the wrong answer this
        attribute is spelled out to exclude: measured at 220 mm, 128, 250 Hz/px the two are
        -585.664336 and -294.638695 1/m.
        """
        return float(self._require_flyback().area)

    @property
    def flyback_duration_s(self) -> float:
        """Seconds the fly-back occupies.  Monopolar only; ``echo_spacing_s`` lengthens it."""
        return float(pp.calc_duration(self._require_flyback()))

    # ----------------------------------------------------------------------- assembly
    def build(
        self, *, acquire: bool = True, phase_deg: float = 0.0, offset_mm: float = 0.0,
    ) -> LogicBlock:
        """
        Return prephaser, readout gradient and -- unless `acquire` is false -- the ADC.

        With `echoes` above 1 the readout gradient becomes the whole train: one prephaser, the
        lobes, the fly-backs a monopolar train needs, one ADC per echo, the barriers that keep
        every lobe a ``trap``, and the ``ECO`` and ``REV`` labels a reconstruction reads the train
        back with.

        Parameters
        ----------
        acquire
            ``False`` plays identical gradients with no ADC and no labels, for a dummy repetition.
        phase_deg
            Receiver phase for this repetition, degrees.  **Give it the same value the
            excitation got**: the receiver is phase-locked to the transmitter, so an RF-spoiling
            schedule that moves one has to move the other.  One value covers the whole train.
        offset_mm
            Shift the FOV along `axis`, by demodulating at a shifted frequency.
        """
        out = LogicBlock()
        start = 0.0
        if self.prephase:
            at = 0.0
            for lobe in self._prephasers:
                out.add(at, lobe)
                at += float(pp.calc_duration(lobe))
            start = at
        adc = (
            self._adc_for(float(offset_mm) / 1e3, float(np.deg2rad(phase_deg))) if acquire
            else None
        )
        period = self._period_s or 0.0
        for echo in range(self.echoes):
            t0 = start + echo * period
            out.add(t0, self.gx if self.polarity_of(echo) > 0 else self._gx_reverse)
            if adc is not None:
                out.add(t0, adc)
                self._label(out, t0, echo)
            if echo + 1 < self.echoes:
                self._seam(out, t0 + self._lobe_s, start + (echo + 1) * period)
        return out

    def _label(self, out: LogicBlock, t0: float, echo: int) -> None:
        """
        Emit ``ECO`` and ``REV`` for one acquired echo, and nothing at ``echoes=1``.

        Nothing at ``echoes=1`` because that case is byte-identical to this module before the
        train existed, and two extra events are not byte-identical.  Both labels on every echo of
        a train, though -- including echo 0 and including monopolar -- because pulseq labels are
        **stateful** and an omitted ``SET`` inherits rather than defaults.
        """
        if self.echoes == 1:
            return
        out.add(t0, pp.make_label(type='SET', label='ECO', value=int(echo)))
        out.add(t0, pp.make_label(type='SET', label='REV',
                                  value=int(self.polarity_of(echo) < 0)))

    def _seam(self, out: LogicBlock, seam_s: float, next_lobe_s: float) -> None:
        """
        State the block boundary between two echoes, one barrier per gradient edge.

        Monopolar gets **two** -- at the fly-back's start and at its end -- because a block may
        hold only one gradient per axis, so a boundary is needed on each side of the fly-back
        rather than only at the seam.  Bipolar gets one, because there is only the seam.

        **These are not load-bearing against this compiler, and pretending otherwise would be
        worse than leaving them out.**  Measured with two, one and no barriers, at three
        bandwidths and both polarities: identical block counts, every gradient still a ``trap``,
        no ``merge`` warning.  :func:`~seqcraft.compiler.boundaries.find_boundaries` already
        offers every gradient edge as a boundary candidate and accepts it when it falls strictly
        inside no gradient, and every seam in this design is exactly that -- unlike
        :class:`~seqcraft.modules.EPI2D`'s, which a centred blip covers.

        What they buy is that the structure is *stated*.  The rule that makes them redundant is
        documented as a preference rather than a guarantee, and a barrier is this library's
        documented way of saying where a block ends.  Zero events in the file, zero difference in
        the output, and the design no longer depends on a preference holding.
        """
        out.add(seam_s, barrier('seam'))
        if self._flyback is not None:
            out.add(seam_s, self._flyback)
            out.add(next_lobe_s, barrier('lobe'))

    def _adc_for(self, offset_m: float, phase_rad: float) -> Event:
        """Return the ADC with the receiver phase and the FOV shift this call asked for."""
        if offset_m == 0.0 and phase_rad == 0.0:
            return self.adc
        freq_offset = float(self.gx.amplitude) * offset_m
        # Referenced to the echo rather than to the first sample, so the shift adds no phase at
        # k = 0.  Omitting this term is a linear phase across the image that reads as a gradient
        # delay -- the same correction shift_slice applies at an RF pulse's effective centre.
        phase_offset = phase_rad - 2 * pi * freq_offset * (
            self._echo_in_gx - float(self.adc.delay)
        )
        return derive(self.adc, freq_offset=freq_offset, phase_offset=phase_offset)

    # -------------------------------------------------------------------------- design
    def _resolve_samples(self) -> tuple[int, int]:
        """
        Return ``(num_samples, pre_echo_samples)``, **even** under ``'bipolar'``.

        Even because the sampling window has to be exactly centred in its lobe and the lobe has to
        land on the gradient raster; an odd count leaves the guard half a raster short, and the
        two polarities then sample grids that are not the same grid.  At ``matrix=128`` and
        ``partial_fourier`` of 1.0 or 0.75 the count already is even -- the forcing is for the odd
        cases, and it is one ``max(2, 2*round(n/2))``.
        """
        count = max(1, int(round(self.partial_fourier * self.matrix)))
        if self.polarity == 'bipolar':
            count = max(2, 2 * int(round(count / 2)))
        pre_echo = count - (self.matrix - self.matrix // 2)
        if pre_echo < 0:
            msg = format_error(
                f'partial_fourier = {self.partial_fourier:g} leaves no pre-echo samples at '
                f'matrix = {self.matrix}.',
                {'partial_fourier': self.partial_fourier, 'matrix': self.matrix,
                 'num_samples': count},
                ['partial_fourier below ~0.55 has too little pre-echo data for a phase estimate'],
            )
            raise ConfigurationError(msg)
        return self._check_sample_count(count), pre_echo

    def _check_sample_count(self, count: int) -> int:
        """
        Return `count`, refusing a sampling geometry the receiver cannot digitise.

        A scanner accepts ADC sample counts in multiples of ``opts.adc_samples_divisor``, and the
        requirement is on the **count**, not on the parity of the matrix: ``matrix=65`` at
        ``partial_fourier=1.0`` gives 65 samples, and an even matrix gives an illegal count just
        as easily -- 0.6 of 128 is 77, which is why ``examples/fse_2d`` runs HASTE at 0.625.

        Checked here because this is the layer that **computes** the count.  The compiler still
        refuses an illegal sequence, but it does so after the caller has built one, naming a block
        index rather than the two arguments that produced it.

        **Nothing is rounded.**  Moving 65 samples to 64 or 68 would quietly change the partial
        Fourier fraction, the k-space extent, the semantic centre sample and the sampling density
        -- four things this module exists to be exact about.
        """
        divisor = int(getattr(self.opts, 'adc_samples_divisor', 1) or 1)
        if divisor <= 1 or count % divisor == 0:
            return count

        below, above = (count // divisor) * divisor, (count // divisor + 1) * divisor
        # A suggestion has to be one this module would itself accept: at or below a full readout,
        # and leaving the post-echo half intact.  Offering `partial_fourier=1.05` would send the
        # caller into the next refusal.
        floor = self.matrix - self.matrix // 2
        fixes = [
            f'partial_fourier={n / self.matrix:.6g} gives {n} samples at matrix={self.matrix}'
            for n in (below, above) if floor <= n <= self.matrix
        ]
        fixes.append(f'matrix={above} at partial_fourier=1.0 gives {above} samples')
        msg = format_error(
            f'the requested geometry produces {count} ADC samples, but this scanner requires a '
            f'count divisible by {divisor}.',
            {'matrix': self.matrix, 'partial_fourier': self.partial_fourier,
             'num_samples': count, 'adc_samples_divisor': divisor},
            fixes,
        )
        raise ConfigurationError(msg)

    def _resolve_dwell(self, bandwidth_hz_px: float | None, dwell_s: float | None) -> float:
        """
        Return the dwell time, from whichever of the two was given.

        On the ADC raster, except under ``'bipolar'``, where it is snapped **up** onto the coarser
        quantum :func:`~seqcraft.modules._support.dwell_quantum` solves for -- because a centred
        window makes the guard an ADC ``delay``, and pypulseq's timing check divides those by
        ``rf_raster_time`` and not by the ten-times-finer ADC raster.  At 128 samples that is a
        multiple of 500 ns, so **the same request produces a different bandwidth in the two
        modes**: 500 Hz/px asked for is 500.8 achieved monopolar and 488.3 bipolar.  That is a
        fact about the file, it belongs in the protocol table, and :attr:`bandwidth_hz_px` reads
        it back -- which is why that attribute has always been the achieved value.

        Up rather than to nearest, for the reason :class:`~seqcraft.modules.EPI2D` gives: rounding
        down raises the bandwidth, which shortens the lobe, which is the direction that turns a
        feasible protocol into a refusal.
        """
        if (bandwidth_hz_px is None) == (dwell_s is None):
            msg = format_error(
                'exactly one of bandwidth_hz_px and dwell_s is required.',
                {'bandwidth_hz_px': bandwidth_hz_px, 'dwell_s': dwell_s},
                [
                    'bandwidth_hz_px=200 is how a protocol quotes it',
                    'dwell_s=5e-6 is how the receiver is actually programmed',
                    'they are one number: dwell = 1 / (bandwidth_hz_px * matrix)',
                ],
            )
            raise ConfigurationError(msg)
        wanted = (
            float(dwell_s) if dwell_s is not None
            else 1.0 / (require_positive(bandwidth_hz_px, 'bandwidth_hz_px') * self.matrix)
        )
        if self.polarity == 'bipolar':
            quantum = to_ticks(dwell_quantum(self.num_samples, opts=self.opts))
            asked = to_ticks(require_positive(wanted, 'dwell_s'))
            return from_ticks(max(quantum, -(-asked // quantum) * quantum))
        raster = float(self.opts.adc_raster_time)
        snapped = round(require_positive(wanted, 'dwell_s') / raster) * raster
        if snapped <= 0.0:
            msg = format_error(
                f'a dwell of {wanted * 1e9:.3f} ns rounds to zero on the '
                f'{raster * 1e9:.0f} ns ADC raster.',
                {'dwell_s': wanted, 'adc_raster_time': raster},
                ['lower bandwidth_hz_px, or reduce matrix'],
            )
            raise ConfigurationError(msg)
        return snapped

    def _readout_amplitude(self) -> float:
        """
        Return ``dk / dwell``, having checked the amplifier can reach it.

        **No amplitude re-solve, in either lobe design.**  :class:`~seqcraft.modules.EPI2D`
        measures its amplitude into place because ramp sampling makes the window's traversal not
        ``G*N*tau``; here the window is entirely on the flat top, so ``dk / dwell`` is exact.
        Saying so is the point of this docstring, or the next reader adds the solve back.
        """
        amplitude = self.dk_per_m / self.dwell_s
        if abs(amplitude) > self.opts.max_grad:
            need = self.dk_per_m / self.opts.max_grad
            msg = format_error(
                f'the readout needs {amplitude:.0f} Hz/m, above the '
                f'{self.opts.max_grad:.0f} Hz/m limit.',
                {'fov_mm': self.fov_mm, 'dwell_s': self.dwell_s,
                 'bandwidth_hz_px': self.bandwidth_hz_px},
                [
                    f'lower the bandwidth to at most {1.0 / (need * self.matrix):.0f} Hz/px',
                    'or widen fov_mm, which lowers dk and so lowers the amplitude',
                ],
            )
            raise ConfigurationError(msg)
        return amplitude

    def _design_readout(self) -> tuple[Event, Event]:
        """
        Return the readout gradient and its ADC, sampling the flat top.

        The single-echo lobe, and -- unchanged, deliberately, not re-derived -- the monopolar
        train's lobe too.  Its flat top is rounded *up* onto the gradient raster and the slack
        left after the last sample, which changes no k position when every lobe in the train is
        this same lobe.  ``'bipolar'`` is the case where that stops being true; see
        :meth:`_design_lobe_centred`.
        """
        amplitude = self._readout_amplitude()
        # The flat top is the ADC window rounded up onto the gradient raster; the few spare
        # microseconds sit after the last sample, where they change no k-space position.
        flat_time = ceil_raster(self.num_samples * self.dwell_s, self.opts.grad_raster_time)
        gx = pp.make_trapezoid(
            channel=self.axis, amplitude=amplitude, flat_time=flat_time, system=self.opts,
        )
        # The sampling offset goes in the ADC's own delay, and the node goes where the gradient
        # starts: pp.make_adc raises a delay below adc_dead_time up to it, and seqcraft preserves
        # an event's delay rather than folding it away, so a node offset would add to it twice.
        adc = pp.make_adc(
            num_samples=self.num_samples, dwell=self.dwell_s, delay=gx.rise_time, system=self.opts,
        )
        # make_adc raises a delay below adc_dead_time up to it, so on a scanner whose receive
        # dead time exceeds the ramp the window starts inside the flat top and would otherwise
        # run off the end of it -- sampling the fall ramp, which is not what the k arithmetic
        # here assumes.  Lengthen the flat top instead; the ramp time is unaffected.
        overrun = (
            float(adc.delay) - float(gx.rise_time) + self.num_samples * self.dwell_s - flat_time
        )
        if overrun > 0.0:
            gx = pp.make_trapezoid(
                channel=self.axis,
                amplitude=amplitude,
                flat_time=ceil_raster(flat_time + overrun, self.opts.grad_raster_time),
                system=self.opts,
            )
        return gx, adc

    def _design_lobe_centred(self, requested_s: float | None) -> tuple[Event, Event]:
        """
        Return the **bipolar** lobe and its ADC, with the sampling window exactly centred.

        Rule 1, and it holds with equality::

            T == 2*guard + num_samples*dwell

        A trapezoid is symmetric about its own midpoint, so sample *i* of a forward lobe and
        sample ``N-1-i`` of a reverse one land on the same k position **if and only if** that
        equation is exact.  Rounding the flat top up and leaving the slack at one end -- which is
        what :meth:`_design_readout` correctly does for a line -- puts the two polarities 0.2051
        ``dk`` apart at the reference protocol, and no k-space extent check notices.

        The guard is whatever is left after the window, split equally between the two ends.  It
        has to clear the receiver dead time *and* the rise time, so that the window is entirely on
        the flat top; the flat top absorbs the raster rounding instead of the tail of the window.
        A requested `echo_spacing_s` lengthens ``T``, which grows both guards symmetrically and
        keeps Rule 1 rather than breaking it.
        """
        amplitude = self._readout_amplitude()
        sample_s = from_ticks(self.num_samples * to_ticks(self.dwell_s))
        rise_s = ceil_raster(amplitude / float(self.opts.max_slew), self.opts.grad_raster_time)
        guard_min = max(float(self.opts.adc_dead_time), rise_s)
        total_s = self._resolve_period(
            requested_s, ceil_raster(sample_s + 2 * guard_min, self.opts.grad_raster_time),
        )
        guard_s = halve_onto(total_s, sample_s, self.opts.rf_raster_time, opts=self.opts)
        gx = pp.make_trapezoid(
            channel=self.axis, amplitude=amplitude, rise_time=rise_s,
            flat_time=total_s - 2 * rise_s, system=self.opts,
        )
        adc = pp.make_adc(
            num_samples=self.num_samples, dwell=self.dwell_s, delay=guard_s, system=self.opts,
        )
        if abs(float(adc.delay) - guard_s) > EPS:
            msg = format_error(
                f'the receiver needs {float(adc.delay) * 1e6:.1f} us before the first sample and '
                f'the lobe leaves {guard_s * 1e6:.1f} us.',
                {'adc_dead_time': float(self.opts.adc_dead_time), 'guard_s': guard_s,
                 'bandwidth_hz_px': self.bandwidth_hz_px},
                [
                    'lower bandwidth_hz_px, which lengthens the lobe and the guard with it',
                    'or pass a longer echo_spacing_s, which does the same without moving the '
                    'sampling rate',
                ],
            )
            raise ConfigurationError(msg)
        # Rule 1, stated as an invariant rather than as an intention.  Measured 0.000 ns; the only
        # way it can fail is a scanner whose gradient raster is not an even number of RF rasters,
        # which is worth naming here rather than discovering in the field map.
        if to_ticks(total_s) - 2 * to_ticks(guard_s) - to_ticks(sample_s):
            msg = format_error(
                f'the sampling window is not centred in its lobe: '
                f'T - 2*guard - N*dwell = '
                f'{from_ticks(to_ticks(total_s) - 2 * to_ticks(guard_s) - to_ticks(sample_s)) * 1e9:.3f} ns.',
                {'echo_spacing_s': total_s, 'guard_s': guard_s, 'dwell_s': self.dwell_s,
                 'num_samples': self.num_samples},
                ['bipolar needs grad_raster_time to be an even multiple of rf_raster_time',
                 "polarity='monopolar' has no reverse lobe and does not need Rule 1"],
            )
            raise ConfigurationError(msg)
        self._period_s = total_s
        return gx, adc

    def _design_flyback(self, requested_s: float | None) -> Event:
        """
        Return the **monopolar** fly-back: minus the lobe's *total* area, ramps included.

        It returns k to where the prephaser left it, so that lobe *n+1* traverses the same
        trajectory lobe *n* did.  Three wrong answers compile, and the module docstring measures
        all three; the dangerous one is ``-area_to_echo_per_m``, because it is a number this
        module already exposes and it is only a factor of two out.

        The shortest trapezoid at that area unless `echo_spacing_s` asks for a longer period, in
        which case it is **stretched** rather than followed by a gap: a longer fly-back at the same
        area is a lower amplitude and a lower slew rate, which is the eddy-current-cheaper of the
        two ways to spend the time.
        """
        area = -area_until(self.gx, float(pp.calc_duration(self.gx)))
        shortest = pp.make_trapezoid(channel=self.axis, area=area, system=self.opts)
        floor_s = float(pp.calc_duration(shortest))
        self._period_s = self._resolve_period(
            requested_s, ceil_raster(self._lobe_s + floor_s, self.opts.grad_raster_time),
        )
        wanted = self._period_s - self._lobe_s
        if wanted <= floor_s + EPS:
            return shortest
        return pp.make_trapezoid(
            channel=self.axis, area=area, duration=wanted, system=self.opts,
        )

    def _resolve_period(self, requested_s: float | None, minimum_s: float) -> float:
        """Return the seam-to-seam period, refusing one below `minimum_s` and naming it.

        The minimum is an argument rather than read off ``self`` because the two lobe designs
        compute it differently -- bipolar from the window and the guards, monopolar from the lobe
        and the shortest fly-back -- and each knows it before it has anywhere to store it.
        """
        self._min_period_s = minimum_s
        if requested_s is None:
            return minimum_s
        wanted = ceil_raster(
            require_positive(requested_s, 'echo_spacing_s'), self.opts.grad_raster_time,
        )
        if wanted < minimum_s - EPS:
            msg = format_error(
                f'echo_spacing_s = {float(requested_s) * 1e6:.1f} us is shorter than the '
                f'minimum this train needs.',
                {'echo_spacing_s': requested_s, 'min_echo_spacing_s': minimum_s,
                 'polarity': self.polarity, 'bandwidth_hz_px': self.bandwidth_hz_px},
                [
                    f'pass echo_spacing_s >= {minimum_s:.6g}',
                    'or echo_spacing_s=None for the shortest legal period',
                    'a higher bandwidth_hz_px shortens the lobe, and with it the minimum',
                ],
            )
            raise ConfigurationError(msg)
        return wanted

    def _design_prephaser(self, requested_s: float | None) -> tuple[Event, ...]:
        """Return the prephaser lobes, in play order."""
        if self.null_moment_order >= 1:
            return self._design_velocity_compensated(requested_s)
        return (self._design_prephaser_m0(requested_s),)

    def _design_prephaser_m0(self, requested_s: float | None) -> Event:
        """Return the prephaser that puts k = 0 at the echo, stretched if one was requested."""
        area = -self.area_to_echo_per_m
        shortest = pp.make_trapezoid(channel=self.axis, area=area, system=self.opts)
        self._min_prephaser_duration_s = float(pp.calc_duration(shortest))
        if requested_s is None:
            return shortest
        wanted = ceil_raster(
            require_positive(requested_s, 'prephaser_duration_s'), self.opts.grad_raster_time,
        )
        # EPS, not exact: a caller passing `max(ro.prephaser_duration_s, pe.min_duration_s)`
        # back in hands us its own minimum, and snapping that onto the raster can land one ulp
        # below where it came from.  Times a nanosecond apart are the same time.
        if wanted < self._min_prephaser_duration_s - EPS:
            msg = format_error(
                f'prephaser_duration_s = {requested_s * 1e6:.1f} us is shorter than the '
                f'minimum this readout needs.',
                {'prephaser_duration_s': requested_s,
                 'min_prephaser_duration_s': self._min_prephaser_duration_s,
                 'prephaser_area_per_m': area},
                [
                    f'pass prephaser_duration_s >= {self._min_prephaser_duration_s:.6g}',
                    'or pass None for the shortest legal prephaser',
                ],
            )
            raise ConfigurationError(msg)
        return pp.make_trapezoid(
            channel=self.axis, area=area, duration=wanted, system=self.opts,
        )

    def _design_velocity_compensated(self, requested_s: float | None) -> tuple[Event, Event]:
        """
        Return the two prephaser lobes that null ``m0`` **and** ``m1`` at the echo.

        The construction is D1's velocity-compensated readout: three lobes before the echo rather
        than two, here the two winders plus the readout's own ramp-up and half flat top.  The
        readout lobe is fixed, so what is solved for is the pair of winder areas.

        Taking the echo as the time origin and writing the readout's own contribution up to it as
        ``m0_ro`` and ``m1_ro``, two lobes of equal duration `D` played back to back and ending
        where the readout begins have centroids at ``t_rs - 3D/2`` and ``t_rs - D/2``.  A
        symmetric trapezoid's first moment about an external origin is its area times its
        centroid, exactly, so the two conditions are **linear** in the areas::

            A1 + A2                                  = -m0_ro
            A1 (t_rs - 3D/2) + A2 (t_rs - D/2)       = -m1_ro

        which solves to ``A1 = (S t_rs - M - D S / 2) / D`` with ``S = -m0_ro`` and
        ``M = -m1_ro``, and ``A2 = S - A1``.  Both ``m0_ro``, ``m1_ro`` and ``t_rs`` are
        properties of the readout lobe alone, so they are computed once.

        `D` itself is not closed form: the areas depend on it and the largest area a trapezoid of
        duration `D` can carry depends on it too, so the shortest legal pair is found by
        bisection on the gradient raster and then checked to be minimal.  The **physics** is the
        linear solve above; the search only picks the shortest duration that fits it.
        """
        raster = float(self.opts.grad_raster_time)
        shortest = self._shortest_compensated_half_s(raster)
        self._min_prephaser_duration_s = 2.0 * shortest
        half = shortest
        if requested_s is not None:
            wanted = ceil_raster(
                require_positive(requested_s, 'prephaser_duration_s'), 2.0 * raster,
            )
            if wanted < self._min_prephaser_duration_s - EPS:
                self._refuse_short_compensated_prephaser(requested_s)
            half = wanted / 2.0
        return tuple(  # type: ignore[return-value]
            pp.make_trapezoid(channel=self.axis, area=area, duration=half, system=self.opts)
            for area in self._compensated_areas(half)
        )

    def _compensated_areas(self, half_s: float) -> tuple[float, float]:
        """The two winder areas, in play order, for a pair of lobes `half_s` long each."""
        m0_ro, m1_ro = (self._readout_moment_to_echo(order) for order in (0, 1))
        target_area, target_moment = -m0_ro, -m1_ro
        first = (target_area * -self._echo_in_gx - target_moment
                 - 0.5 * half_s * target_area) / half_s
        return first, target_area - first

    def _readout_moment_to_echo(self, order: int) -> float:
        """
        The readout lobe's own `order`-th moment up to the echo, about the echo.

        Integrated from the lobe's knots, so ramp sampling, partial Fourier and the half-dwell
        sample offset are all already in it -- the same reason :attr:`area_to_echo_per_m` is
        measured rather than derived.  It does not depend on where the lobe is placed, because
        the origin travels with it.
        """
        times, amps = knots_of(self.gx)
        cut = int(np.searchsorted(times, self._echo_in_gx))
        edge = float(np.interp(self._echo_in_gx, times, amps))
        times = np.append(times[:cut], self._echo_in_gx) - self._echo_in_gx
        return float(pwl_moment(times, np.append(amps[:cut], edge), order))

    def _limit_areas(self) -> tuple[float, float]:
        """
        What the two winder areas tend to as the lobes are made longer, 1/m.

        Writing ``S = -m0_ro`` and ``C = S t_rs - M``, the solve in
        :meth:`_design_velocity_compensated` is ``A1(D) = C/D - S/2`` and ``A2(D) = 3S/2 - C/D``.
        So each area is a constant plus a term in ``1/D``, and these are those constants.
        """
        target_area = -self._readout_moment_to_echo(0)
        return -0.5 * target_area, 1.5 * target_area

    def _largest_trapezoid_area(self, duration_s: float) -> float:
        """
        The largest area a symmetric trapezoid of `duration_s` can carry, 1/m.

        Triangular while the ramps have not reached ``max_grad``, trapezoidal after.  Strictly
        increasing in the duration, which is half of why :meth:`_compensated_pair_fits` is
        monotone.
        """
        amplitude = min(float(self.opts.max_grad), float(self.opts.max_slew) * duration_s / 2.0)
        return amplitude * (duration_s - amplitude / float(self.opts.max_slew))

    def _compensated_pair_fits(self, half_s: float) -> bool:
        """
        Whether a pair of lobes `half_s` long can carry the areas the solve asks of them.

        **This predicate is monotone in `half_s`, and that is what licenses the bisection in
        :meth:`_shortest_compensated_half_s`.**  Each area is ``A_i(D) = L_i + K_i / D``, so as
        `D` grows ``A_i`` moves monotonically from where it is towards its limit ``L_i`` and never
        past it -- which puts every later value inside the closed interval between the two, and
        therefore bounds ``|A_i(D')| <= max(|A_i(D)|, |L_i|)`` for every ``D' >= D``.  Requiring
        the **limit** areas to fit as well as the current ones makes that bound a statement about
        this predicate: if it holds at `D` it holds at every longer duration, because
        :meth:`_largest_trapezoid_area` only grows.

        Without the limit condition the predicate would not be monotone in general, because
        ``A_i`` may cross zero on its way to ``L_i`` and come back out larger.  On every readout
        geometry measured it does **not** bind -- the solved areas dominate their limits
        throughout the feasible range -- so it changes no answer here and is carried to make the
        up-set property hold without assuming a sign relationship between ``A_i`` and ``L_i``
        that a geometry nobody has tried might not have.
        """
        largest = self._largest_trapezoid_area(half_s)
        if largest <= 0.0:
            return False
        areas = self._compensated_areas(half_s) + self._limit_areas()
        return all(abs(area) <= largest for area in areas)

    def _shortest_compensated_half_s(self, raster: float) -> float:
        """
        The shortest legal lobe duration for the pair, on the raster.

        Bracketed by doubling and then bisected, which is valid because
        :meth:`_compensated_pair_fits` is monotone -- see its docstring for why.  The answer is
        then checked against the raster step below it, so it is the minimum on the lattice rather
        than merely a feasible point.

        A final loop widens it if the trapezoid pypulseq actually builds at that duration lands
        outside the limits: the predicate works from the smooth largest-area curve, and a rise
        time rounded onto the raster can cost a little of it at the boundary.
        """
        high = raster
        for _ in range(40):
            if self._compensated_pair_fits(high):
                break
            high *= 2.0
        else:  # pragma: no cover - 2^40 rasters is not a reachable gradient system
            self._refuse_impossible_compensation()
        low = high / 2.0
        while high - low > raster:
            middle = ceil_raster((low + high) / 2.0, raster)
            if middle >= high:
                break
            if self._compensated_pair_fits(middle):
                high = middle
            else:
                low = middle
        while not self._built_pair_is_legal(high):  # pragma: no cover - raster boundary only
            high += raster
        return float(high)

    def _built_pair_is_legal(self, half_s: float) -> bool:
        """Whether the lobes pypulseq builds at `half_s` are inside the gradient limits."""
        for area in self._compensated_areas(half_s):
            try:
                lobe = pp.make_trapezoid(
                    channel=self.axis, area=area, duration=half_s, system=self.opts,
                )
            except (ValueError, ZeroDivisionError):
                return False
            amplitude = abs(float(lobe.amplitude))
            if amplitude > float(self.opts.max_grad) * (1.0 + 1e-9):
                return False
            if amplitude > float(self.opts.max_slew) * float(lobe.rise_time) * (1.0 + 1e-9):
                return False
        return True

    # -------------------------------------------------------------------- the refusals
    def _echo_in_lobe_s(self, echo: int) -> float:
        """Seconds from a lobe's own start to ``k = 0``, which differs with the polarity."""
        if self.polarity_of(echo) > 0:
            return self._echo_in_gx
        return float(self.adc.delay) + (self.echo_sample(echo) + 0.5) * self.dwell_s

    def _require_echo(self, echo: int) -> int:
        """Return `echo` as an index into this train, refusing one it does not have."""
        index = int(echo)
        if not 0 <= index < self.echoes:
            msg = format_error(
                f'echo {index} is outside this train, which has {self.echoes}.',
                {'echo': index, 'echoes': self.echoes},
                ['echoes are zero-based: 0 is the first, echoes - 1 the last',
                 'te_s is the whole array, and len(te_s) == echoes'],
            )
            raise ConfigurationError(msg)
        return index

    def _resolve_polarity(self, polarity: str | None) -> str | None:
        """
        Return the validated `polarity`, refusing one that cannot take effect and one that is not.

        Required above one echo and refused at one, both named in the message.  ``None`` has no
        default above one echo on purpose: the two modes produce files that differ in the dwell,
        the echo times, the ordering of every second echo and the block count, and look identical
        in a protocol printout.  Writing it out at every call site costs eleven characters.
        """
        if self.echoes == 1:
            if polarity is not None:
                msg = format_error(
                    f'polarity = {polarity!r} was passed with echoes=1, so it cannot take '
                    f'effect: a single echo has no train to alternate.',
                    {'echoes': 1, 'polarity': polarity},
                    ['drop polarity for a single-echo readout',
                     'or pass echoes >= 2, which is what polarity describes'],
                )
                raise ConfigurationError(msg)
            return None
        if polarity is None:
            msg = format_error(
                f'polarity is required with echoes={self.echoes}, and has no default.',
                {'echoes': self.echoes, 'polarity': None},
                [f"polarity='{name}': {why}" for name, why in _POLARITIES.items()]
                + ['a reconstruction has to know which one it got, and neither is safe to assume'],
            )
            raise ConfigurationError(msg)
        if polarity not in _POLARITIES:
            msg = format_error(
                f'polarity must be one of {", ".join(repr(k) for k in _POLARITIES)}, '
                f'got {polarity!r}.',
                {'polarity': polarity, 'echoes': self.echoes},
                [f"'{name}': {why}" for name, why in _POLARITIES.items()],
            )
            raise ConfigurationError(msg)
        return polarity

    def _refuse_bipolar_spin_echo(self) -> None:
        """
        Refuse ``polarity='bipolar'`` with ``prephase=False``, which nobody has measured.

        Legal in principle and out of scope here: a multi-echo *spin* echo with alternating
        readout polarity has to balance the reversed lobes across the refocusing pulse as well,
        which is a different question from the one ``prephase=False`` answers today.  Monopolar
        with ``prephase=False`` is built and tested.
        """
        if self.polarity == 'bipolar' and not self.prephase:
            msg = format_error(
                "polarity='bipolar' with prephase=False is not supported.",
                {'polarity': 'bipolar', 'prephase': False},
                [
                    "polarity='monopolar' is the multi-echo spin-echo readout that is measured",
                    'balancing reversed lobes across a refocusing pulse needs the crusher pair '
                    'to be a pair for both polarities, which nobody has measured here',
                ],
            )
            raise ConfigurationError(msg)

    def _refuse_echo_spacing(self, requested_s: float | None) -> None:
        """Refuse `echo_spacing_s` at one echo, where there is no seam for it to describe."""
        if self.echoes == 1 and requested_s is not None:
            msg = format_error(
                f'echo_spacing_s = {float(requested_s) * 1e6:.1f} us was passed with echoes=1, '
                f'so it cannot take effect: there is no second echo to space.',
                {'echoes': 1, 'echo_spacing_s': requested_s},
                ['drop echo_spacing_s for a single-echo readout',
                 'or pass echoes >= 2 and a polarity, whose period it then sets'],
            )
            raise ConfigurationError(msg)

    def _require_train(self, name: str) -> float:
        """Return the period, or refuse a question about a train that was not designed."""
        if self._period_s is None:
            msg = format_error(
                f'this readout has one echo, so {name} describes nothing.',
                {'echoes': self.echoes, 'time_to_echo': self.time_to_echo()},
                [
                    'read time_to_echo() instead: with one echo the block has one echo time',
                    'or pass echoes >= 2 and a polarity, which is what makes a period exist',
                ],
            )
            raise ConfigurationError(msg)
        return self._period_s

    def _require_flyback(self) -> Event:
        """Return the fly-back, or refuse a question about an event that was not designed."""
        if self._flyback is None:
            why = (
                "polarity='bipolar' alternates the lobes instead, so there is no fly-back"
                if self.polarity == 'bipolar'
                else 'this readout has one echo, so there is nothing to fly back between'
            )
            msg = format_error(
                f'this readout has no fly-back: {why}.',
                {'echoes': self.echoes, 'polarity': self.polarity},
                [
                    'echo_spacing_s is what lengthens the period under either polarity',
                    "polarity='monopolar' with echoes >= 2 is what designs a fly-back",
                ],
            )
            raise ConfigurationError(msg)
        return self._flyback

    def _require_prephasers(self) -> tuple[Event, ...]:
        """Return the prephaser lobes, or refuse a question about events not designed."""
        if not self._prephasers:
            msg = format_error(
                'this readout has no prephaser, because prephase=False.',
                {'prephase': False, 'area_to_echo_per_m': self.area_to_echo_per_m},
                [
                    'read area_to_echo_per_m instead: it is what the readout accumulates by the '
                    'echo, which is the number a crusher pair is balanced around',
                    'or pass prephase=True for a readout that cancels its own dephasing',
                ],
            )
            raise ConfigurationError(msg)
        return self._prephasers

    def _check_null_moment_order(self, order: int) -> int:
        """Return `order` having checked this readout can null moments up to it."""
        order = int(order)
        if order == 0:
            return order
        if order == 1 and self.prephase:
            return order
        if order == 1:
            msg = format_error(
                'null_moment_order=1 needs a prephaser to reshape, and prephase=False.',
                {'null_moment_order': order, 'prephase': self.prephase},
                ['pass prephase=True, whose winder is what carries the compensation',
                 'a readout with no prephaser has no waveform here to null a moment with'],
            )
            raise ConfigurationError(msg)
        msg = format_error(
            f'null_moment_order={order} is out of range: 0 or 1.',
            {'null_moment_order': order},
            ['0 nulls the zeroth moment at the echo, which is the ordinary prephaser',
             '1 also nulls the first, which is velocity compensation',
             'nulling the second needs a fourth lobe and is not implemented'],
        )
        raise ConfigurationError(msg)

    def _refuse_short_compensated_prephaser(self, requested_s: float) -> None:
        """Refuse a compensated prephaser too short for the areas it has to carry."""
        msg = format_error(
            f'prephaser_duration_s = {requested_s * 1e6:.1f} us is shorter than the minimum a '
            f'velocity-compensated prephaser needs for this readout.',
            {'prephaser_duration_s': requested_s,
             'min_prephaser_duration_s': self._min_prephaser_duration_s,
             'null_moment_order': self.null_moment_order,
             'prephaser_area_per_m': -self.area_to_echo_per_m},
            [
                f'pass prephaser_duration_s >= {self._min_prephaser_duration_s:.6g}',
                'or pass None for the shortest legal pair',
                'compensation costs time: two winders carry more area than the one that only '
                'nulls the zeroth moment',
            ],
        )
        raise ConfigurationError(msg)

    def _refuse_impossible_compensation(self) -> None:  # pragma: no cover - unreachable
        msg = format_error(
            'no velocity-compensated prephaser fits this gradient system.',
            {'null_moment_order': self.null_moment_order,
             'area_to_echo_per_m': self.area_to_echo_per_m},
            ['pass null_moment_order=0', 'or widen the readout, which lowers the area to cancel'],
        )
        raise ConfigurationError(msg)

    def _refuse_prephaser_duration(self, requested_s: float | None) -> tuple[Event, ...]:
        """Refuse a prephaser duration for a prephaser that will not exist, and return no event."""
        if requested_s is not None:
            msg = format_error(
                f'prephaser_duration_s = {requested_s * 1e6:.1f} us was passed with '
                f'prephase=False, so it cannot take effect.',
                {'prephase': False, 'prephaser_duration_s': requested_s},
                [
                    'drop prephaser_duration_s: with prephase=False there is no prephaser to '
                    'stretch',
                    'or pass prephase=True, whose winder duration it then sets',
                ],
            )
            raise ConfigurationError(msg)
        return ()
