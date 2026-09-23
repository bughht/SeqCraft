r"""
:class:`T2Prep` -- a T2-weighting preparation: tip down, refocus, tip up, spoil.

What it is
----------
A preparation period that converts :math:`T_2` into longitudinal contrast before an imaging train
that has none of its own.  Magnetisation is tipped into the transverse plane, held there for a
declared time while a train of refocusing pulses reverses static dephasing, tipped back to
:math:`+z`, and whatever is left transverse is crushed:

.. code-block:: text

    90x  ->  composite 180 x4  ->  composite tip-up  ->  spoiler
         |<------ prep_time_s ------>|

What survives is weighted by :math:`e^{-\tau/T_2}` rather than by :math:`e^{-\tau/T_2^*}`, and
that distinction is the whole module: without the refocusing train the same interval would decay
at :math:`T_2^*`, which is a property of the magnet rather than of the tissue.

The weighting is a **factor**.  The signal is :math:`M_z\,e^{-\tau/T_2}`, where :math:`M_z` is the
longitudinal magnetisation that arrives -- which under a steady state is not the equilibrium
value, and is not something this module can know.  The claim here is the exponential, not the
absolute signal.

``prep_time_s``, and where it is measured from
----------------------------------------------
The public quantity is the **transverse period**: from the effective tip-down rotation instant to
the effective tip-up rotation instant.  For this mode's hard and composite pulses that reduces to
RF centre to RF centre, and a composite is booked at the midpoint of the group -- symmetric for
the refocusing composites by construction, and stated for the asymmetric tip-up rather than left
to be inferred.  :meth:`T2Prep.time_to_tip_down` and :meth:`T2Prep.time_to_tip_up` report both
ends, so the interval can be measured rather than trusted.

Implementations differ on this bookkeeping -- some book edge to edge, tip-down end to tip-up start
-- and the two differ by one tip-pulse duration, which is milliseconds against a preparation of
tens.  This module takes the physical definition: the period the magnetisation is transverse.

Every RF group starts on the gradient raster, because the seam between two consecutive pulses is
where a block boundary has to be cut and only raster instants can be cut.  So the placement is
quantised first and ``prep_time_s`` is **measured off it** -- it is what the emitted sequence
delivers, within half a raster of what was asked for.  :attr:`T2Prep.requested_prep_time_s` keeps
the request.

The one mode
------------
``mlev4-hard``, and deliberately only that:

============================  ==================================================================
tip-down                      hard 90 about +x
refocusing train              four composites, each 90x-180y-90x, at 1/8, 3/8, 5/8 and 7/8 of
                              ``prep_time_s``
phase pattern                 MLEV-4, ``+ + - -`` across the four composites
tip-up                        composite: 270 about +x, then 360 about -x.  Net -90
selection                     none.  The family is non-selective
spoiler                       after the tip-up, entirely outside ``prep_time_s``
============================  ==================================================================

Every pulse is a hard block pulse at **one** :math:`B_1` amplitude, set by the 180's duration, so
a flip angle is a duration: a 90 is half a 180 and a 360 is twice one.  That is what makes the
composites composites, and it makes the peak-:math:`B_1` question a single number for the whole
block.

Out of scope, each a different physical claim rather than a parameter of this one: adiabatic
BIR-4/AHP tip pairs and the :math:`B_1`-robustness claim that lives with them; a two-composite
``mlev2-hard`` train; an arbitrary refocusing count; the :math:`B_1`-robust variant whose
rephasing lobe lands inside the *imaging* sequence, which no single module can hold; and a
diffusion preparation, which is the same silhouette making a different claim.

This mode makes **no** quantitative :math:`B_0` or :math:`B_1` robustness claim.

References
----------
Bernstein, King & Zhou, *Handbook of MRI Pulse Sequences*, Elsevier 2004, §17.4.
Levitt & Freeman, *J. Magn. Reson.* **43**, 65 (1981) -- the composite pulse.
Brittain et al., *Magn. Reson. Med.* **33**, 689 (1995) -- the preparation in use.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ...design.events import derive
from ...design.logic import LogicBlock
from ...design.module import Module
from ...design.timing import Raster
from ...errors import ConfigurationError, format_error
from .._support import check_peak_b1, require_axis, require_positive, require_usable_max_b1
from ..spoiler import spoiler

if TYPE_CHECKING:
    from pypulseq.opts import Opts

    from ...design.events import Event

__all__ = ['T2Prep']

#: Where the four composite refocusing pulses sit, as fractions of ``prep_time_s``.  Equal
#: spacing with half-intervals at each end is what makes every gap either side of every pulse the
#: same, which is the refocusing condition.
_REFOCUS_FRACTIONS: tuple[float, ...] = (0.125, 0.375, 0.625, 0.875)

#: MLEV-4.  The pattern is over the composites, not over their sub-pulses: ``True`` adds 180
#: degrees to every pulse of that composite.
_MLEV4_INVERTED: tuple[bool, ...] = (False, False, True, True)

#: One composite refocusing pulse, as (flip in degrees, phase in degrees) relative to the
#: composite's own reference.  Levitt & Freeman's 90x-180y-90x.
_COMPOSITE_REFOCUS: tuple[tuple[float, float], ...] = ((90.0, 0.0), (180.0, 90.0), (90.0, 0.0))

#: The tip-up: 270 about +x, then 360 about -x.  Net -90, which returns what the tip-down sent.
_COMPOSITE_TIP_UP: tuple[tuple[float, float], ...] = ((270.0, 0.0), (360.0, 180.0))

#: The tip-down.  A single hard 90 about +x.
_TIP_DOWN: tuple[tuple[float, float], ...] = ((90.0, 0.0),)


class T2Prep(Module):
    r"""
    A T2-weighting preparation: 90, a refocusing train, a tip-up, and a spoiler.

    Parameters
    ----------
    opts
        The scanner.
    prep_time_s
        **The physical target.**  The interval the magnetisation spends transverse, measured
        between the effective tip-down and tip-up rotation instants.  Rounded up onto eight
        gradient rasters so that the four refocusing centres land exactly;
        :attr:`requested_prep_time_s` keeps the request.
    refocus_duration_s
        Duration of the 180 in each composite.  It sets the :math:`B_1` amplitude every pulse in
        the block shares, so a shorter one is a stronger pulse: halve it and peak :math:`B_1`
        doubles.
    spoil_cycles_per_voxel, spoil_axis, spoil_voxel_mm
        The crusher that destroys what the tip-up left transverse.  It lies entirely **after**
        ``prep_time_s``.  `spoil_voxel_mm` is required: the preparation is non-selective, so
        there is no length in this module to count cycles across.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    prep_time_s : float
        The realised transverse interval, in seconds -- **measured** off the placement rather than
        declared, so it is exactly what the emitted sequence delivers.
    requested_prep_time_s : float
        What the caller asked for, before rounding onto the raster.
    min_prep_time_s : float
        The shortest preparation the pulses themselves leave room for.
    b1_hz : float
        The single amplitude every pulse in the block is played at.
    pulses : tuple
        Every RF event, in time order: the tip-down, then four composites of three, then the
        two-pulse tip-up.  Fifteen in all.
    spoiler : LogicBlock
        The crusher, placed after the tip-up.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
    ...                rf_dead_time=100e-6, rf_ringdown_time=30e-6)
    >>> prep = sc.modules.T2Prep(opts=opts, prep_time_s=50e-3, spoil_voxel_mm=4.0)
    >>> len(prep.pulses)
    15
    >>> prep.time_to_tip_up() - prep.time_to_tip_down() == prep.prep_time_s
    True
    >>> round(prep.prep_time_s * 1e3, 3)
    50.004
    """

    def __init__(
        self,
        *,
        opts: Opts,
        prep_time_s: float,
        refocus_duration_s: float = 1e-3,
        spoil_cycles_per_voxel: float = 4.0,
        spoil_axis: str = 'z',
        spoil_voxel_mm: float,
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.requested_prep_time_s = require_positive(prep_time_s, 'prep_time_s')
        self.refocus_duration_s = require_positive(refocus_duration_s, 'refocus_duration_s')

        # One B1 for the whole block, set by the 180.  A hard pulse's flip angle is
        # 2 pi B1 tau, so at a fixed B1 the duration is proportional to the angle -- which is what
        # lets a 90, a 270 and a 360 belong to the same composite.
        self.b1_hz = 0.5 / self.refocus_duration_s
        self._rf_raster = Raster(float(opts.rf_raster_time), 'rf_raster_time')
        self._grad_raster = Raster(float(opts.grad_raster_time), 'grad_raster_time')

        require_usable_max_b1(self.opts, described=self._described())
        self._tip_down = self._group(_TIP_DOWN)
        self._composites = tuple(self._group(_COMPOSITE_REFOCUS, inverted=flip)
                                 for flip in _MLEV4_INVERTED)
        self._tip_up = self._group(_COMPOSITE_TIP_UP)
        self.pulses = tuple(rf for group in (self._tip_down, *self._composites, self._tip_up)
                            for rf in group)
        self._check_b1()

        self.min_prep_time_s = self._shortest_preparation()
        if self.requested_prep_time_s < self.min_prep_time_s:
            self._refuse_short_preparation()
        # Every RF group starts on the gradient raster, because the seam between two consecutive
        # pulses is where the compiler must cut a block boundary and only raster instants are
        # cuttable.  So the placement is quantised first and the realised interval is *measured*
        # off it, rather than declared and then approximated.
        self._tip_up_start = float(self._grad_raster.nearest(
            self.time_to_tip_down() + self.requested_prep_time_s
            - self._effective_centre(self._tip_up, start_s=0.0)))
        self.prep_time_s = self.time_to_tip_up() - self.time_to_tip_down()
        self._composite_starts = tuple(
            float(self._grad_raster.nearest(
                self.time_to_tip_down() + fraction * self.prep_time_s
                - self._effective_centre(group, start_s=0.0)))
            for fraction, group in zip(_REFOCUS_FRACTIONS, self._composites))

        self.spoil_axis = require_axis(spoil_axis, 'spoil_axis')
        self.spoil_cycles_per_voxel = require_positive(spoil_cycles_per_voxel,
                                                       'spoil_cycles_per_voxel')
        self.spoil_voxel_mm = require_positive(spoil_voxel_mm, 'spoil_voxel_mm')
        self.spoiler = spoiler(opts, cycles_per_voxel=self.spoil_cycles_per_voxel,
                               voxel_mm=self.spoil_voxel_mm, axis=self.spoil_axis)

    # ------------------------------------------------------------------ what it knows
    def time_to_tip_down(self) -> float:
        """Seconds from the block's start to the **effective tip-down rotation instant**."""
        return self._effective_centre(self._tip_down, start_s=0.0)

    def time_to_tip_up(self) -> float:
        """
        Seconds from the block's start to the **effective tip-up rotation instant**.

        A composite has no single instant at which it rotates, so the convention has to be stated
        rather than assumed.  This module books the **amplitude-weighted centre** of the group:
        each pulse contributes its own centre, weighted by its flip angle, which at a fixed
        :math:`B_1` is the area under its envelope.

        It reduces correctly at both ends of the family.  For one hard pulse it is that pulse's
        centre.  For the symmetric 90-180-90 refocusing composite it is the centre of the 180, by
        symmetry.  For the asymmetric 270/-360 tip-up it is a number, which is the point.
        """
        return self._effective_centre(self._tip_up, start_s=self._tip_up_start)

    def time_to_refocus(self, index: int) -> float:
        """
        Seconds from the block's start to composite refocusing pulse `index`'s effective centre.

        As realised, not as intended: the group's start is quantised onto the gradient raster, so
        this is within half a raster of ``fraction * prep_time_s`` rather than exactly on it.  The
        refocusing condition is that the gaps either side of each pulse are **equal**, which
        survives that quantisation, and :attr:`interval_gaps_s` is how it is checked.
        """
        return self._effective_centre(self._composites[index],
                                      start_s=self._composite_starts[index])

    @property
    def interval_gaps_s(self) -> tuple[float, ...]:
        """
        The five intervals between consecutive rotation instants, in seconds.

        Tip-down to the first composite, composite to composite, and the last composite to the
        tip-up.  Refocusing works when the dephasing either side of each pulse is equal, so the
        first and last of these should be half of each middle one -- which is what the fractions
        1/8, 3/8, 5/8 and 7/8 buy.
        """
        instants = [self.time_to_tip_down(),
                    *(self.time_to_refocus(k) for k in range(len(_REFOCUS_FRACTIONS))),
                    self.time_to_tip_up()]
        return tuple(b - a for a, b in zip(instants, instants[1:]))

    def time_to_spoiler(self) -> float:
        """
        Seconds from the block's start to the spoiler, which begins after the tip-up.

        Rounded up onto the gradient raster: the RF that precedes it is on the finer RF raster and
        carries dead and ringdown times, so the instant it ends is not generally a gradient start.
        The rounding lands the spoiler **later**, never inside the preparation.
        """
        return float(self._grad_raster.ceil(
            self._tip_up_start + self._duration_of(self._tip_up)))

    @property
    def refocus_fractions(self) -> tuple[float, ...]:
        """Where the composites sit, as fractions of :attr:`prep_time_s`."""
        return _REFOCUS_FRACTIONS

    # ----------------------------------------------------------------------- assembly
    def build(self, *, phase_deg: float = 0.0) -> LogicBlock:
        """
        Return the preparation: fifteen pulses, then the spoiler.

        Parameters
        ----------
        phase_deg
            A carrier phase added to **every** pulse, which rotates the whole preparation about
            `z` and leaves every phase *difference* -- the tip pair's inverse relationship and the
            MLEV pattern -- untouched.  Phase-cycling a preparation across shots is a real thing
            to want.
        """
        out = LogicBlock()
        for start, group in self._placed_groups():
            at = start
            for rf in group:
                out.add(at, derive(rf, phase_offset=float(rf.phase_offset)
                                   + float(np.deg2rad(phase_deg))) if phase_deg else rf)
                at += self._slot_of(rf)
        return out.add(self.time_to_spoiler(), self.spoiler)

    # ------------------------------------------------------------------- the geometry
    def _placed_groups(self) -> tuple[tuple[float, tuple[Event, ...]], ...]:
        """Every RF group with the instant its first pulse starts, in time order."""
        return ((0.0, self._tip_down),
                *zip(self._composite_starts, self._composites),
                (self._tip_up_start, self._tip_up))

    def _group(self, angles: tuple[tuple[float, float], ...], *,
               inverted: bool = False) -> tuple[Event, ...]:
        """Build one hard-pulse group, every pulse at the block's single B1."""
        turn = 180.0 if inverted else 0.0
        return tuple(self._hard(flip_deg, phase_deg + turn) for flip_deg, phase_deg in angles)

    def _hard(self, flip_deg: float, phase_deg: float) -> Event:
        """One hard block pulse: the duration follows the flip angle at a fixed B1."""
        duration_s = self._rf_raster.ceil(flip_deg / 360.0 / self.b1_hz)
        return pp.make_block_pulse(
            flip_angle=float(np.deg2rad(flip_deg)),
            duration=float(duration_s),
            phase_offset=float(np.deg2rad(phase_deg)),
            delay=float(self.opts.rf_dead_time),
            system=self.opts,
            use='preparation',
        )

    def _slot_of(self, rf: Event) -> float:
        """
        The wall time one pulse occupies, rounded up onto the gradient raster.

        Pulses inside a group are played back to back, so the seam between two of them is where
        the compiler must cut a block boundary -- a Pulseq block holds one RF.  An RF event is on
        the finer RF raster and carries dead and ringdown times, so its natural end is not
        generally a legal boundary.  Rounding each pulse's slot up puts every seam on the gradient
        raster, at a cost of under one raster per pulse.
        """
        return float(self._grad_raster.ceil(float(pp.calc_duration(rf))))

    def _duration_of(self, group: tuple[Event, ...]) -> float:
        """The wall time a whole group occupies."""
        return float(sum(self._slot_of(rf) for rf in group))

    def _effective_centre(self, group: tuple[Event, ...], *, start_s: float) -> float:
        """The group's amplitude-weighted rotation instant -- see :meth:`time_to_tip_up`."""
        weighted, total, at = 0.0, 0.0, start_s
        for rf in group:
            centre = at + float(rf.delay) + float(pp.calc_rf_center(rf)[0])
            # The integral of the envelope, which is the flip angle -- NOT the sum of the
            # samples.  `make_block_pulse` emits two samples whatever the duration, so a sum
            # would weight a 360 and a 90 equally and put the centre in the wrong place.
            weight = abs(float(np.trapezoid(np.abs(np.asarray(rf.signal)),
                                            np.asarray(rf.t))))
            weighted += weight * centre
            total += weight
            at += self._slot_of(rf)
        return weighted / total

    def _shortest_preparation(self) -> float:
        """
        Return the shortest ``prep_time_s`` the pulses leave room for, in seconds.

        Three things can bind, and which one does depends on the pulse durations rather than on
        the preparation time, so all three are solved rather than searched:

        * the first composite must start after the tip-down ends;
        * consecutive composites are ``prep_time_s / 4`` apart and must not overlap;
        * the last composite must end before the tip-up begins.
        """
        tip_down = self._duration_of(self._tip_down)
        composite = self._duration_of(self._composites[0])
        first = _REFOCUS_FRACTIONS[0]
        after_tip_down = (tip_down - self.time_to_tip_down() + composite / 2.0) / first
        between = 4.0 * composite
        before_tip_up = (composite / 2.0 + self._effective_centre(self._tip_up, start_s=0.0)
                         ) / (1.0 - _REFOCUS_FRACTIONS[-1])
        return float(self._grad_raster.ceil(max(after_tip_down, between, before_tip_up) / 8.0)
                     * 8.0)

    # -------------------------------------------------------------------- the refusals
    def _described(self) -> str:
        return (f'a {self.refocus_duration_s * 1e3:g} ms hard-pulse T2 preparation at '
                f'{self.b1_hz:.0f} Hz')

    def _check_b1(self) -> None:
        """Refuse the block when any pulse in it exceeds ``opts.max_b1``.

        Every pulse shares one amplitude, so one of them failing means all of them do -- and the
        remedy is the 180's duration, which is what sets that amplitude.
        """
        limit = float(getattr(self.opts, 'max_b1', 0.0) or 0.0)
        remedies = ['lengthen refocus_duration_s'] if not limit > 0.0 else [
            f'pass refocus_duration_s >= {0.5 / limit * 1e3:.2f} ms, which is where this B1 fits',
            'every pulse in the block shares one B1, so the 180 is what sets it',
        ]
        for rf in self.pulses:
            check_peak_b1(rf, self.opts, described=self._described(), remedies=remedies)

    def _refuse_short_preparation(self) -> None:
        msg = format_error(
            f'prep_time_s = {self.requested_prep_time_s * 1e3:g} ms is shorter than these '
            f'pulses fit into.',
            {'prep_time_s_ms': self.requested_prep_time_s * 1e3,
             'min_prep_time_ms': self.min_prep_time_s * 1e3,
             'refocus_duration_ms': self.refocus_duration_s * 1e3},
            [f'pass prep_time_s >= {self.min_prep_time_s * 1e3:g} ms',
             'or shorten refocus_duration_s, which shortens every pulse in the block',
             'a T2 preparation is normally tens of milliseconds -- it is competing with T2'],
        )
        raise ConfigurationError(msg)
