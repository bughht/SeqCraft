r"""
:class:`T2Prep` -- a T2-weighting preparation: tip down, refocus, tip up, spoil.

What it is
----------
A preparation period that converts :math:`T_2` into longitudinal contrast before a gradient-echo
imaging train that would otherwise follow :math:`T_2^*` decay.  Magnetisation is tipped into the
transverse plane, held there for a declared time while a train of refocusing pulses reverses
static dephasing, tipped back to :math:`+z`, and whatever is left transverse is crushed:

.. code-block:: text

    +90  ->  refocusing train  ->  tip-up  ->  spoiler
         |<---- prep_time_s ---->|

What survives is weighted by :math:`e^{-\tau/T_2}` rather than by :math:`e^{-\tau/T_2^*}`.
Without the refocusing train, static dephasing adds to the tissue's intrinsic :math:`T_2` decay
and the transverse signal follows :math:`T_2^*` instead -- which is faster, and magnet-dependent
as well as tissue-dependent.

The weighting is a **factor**.  The signal is :math:`M_z\,e^{-\tau/T_2}`, where :math:`M_z` is the
longitudinal magnetisation present when the preparation begins -- which under a steady state is
not the equilibrium value.  The preparation contributes the exponential; the absolute signal also
depends on that :math:`M_z`, which this module does not set.

The pulses
----------
Hard, non-adiabatic block pulses throughout, non-selective, with an MLEV-4 refocusing train:

==============================  ===============================================================
tip-down                        hard 90 about +x
refocusing train                four composites, each 90x-180y-90x, at 1/8, 3/8, 5/8 and 7/8 of
                                ``prep_time_s``
phase pattern                   MLEV-4, ``+ + - -`` across the four composites
``tip_up='simple'``             a single hard -90.  **The default**
``tip_up='composite_270_360'``  270 about +x then 360 about -x, also a net -90
selection                       none.  The family is non-selective
spoiler                         after the tip-up, entirely outside ``prep_time_s``
==============================  ===============================================================

Composite refocusing pulses and phase cycling are used so that pulse imperfections tend to cancel
across the train rather than accumulate.  **No quantitative** :math:`B_0` **or**
:math:`B_1` **robustness is claimed for either** ``tip_up`` **value**;
``examples/t2prep_gre_2d/02`` compares them under :math:`B_1` and :math:`B_0` offsets and finds no
consistent advantage for the composite form, which costs about 2.6 ms more in RF.

``prep_time_s``
---------------
The public quantity is the **transverse period**, and both of its ends are RF centres:

==================================  ===========================================================
from                                the centre of the ``+90`` tip-down
to, ``tip_up='simple'``             the centre of the ``-90``
to, ``tip_up='composite_270_360'``  the centre of the **270**
==================================  ===========================================================

Every RF group starts on the gradient raster, because the seam between two consecutive pulses is
where a block boundary has to be cut and only raster instants can be cut.  So the placement is
quantised first and ``prep_time_s`` is **measured off it**;
:attr:`T2Prep.requested_prep_time_s` keeps the request.

One :math:`B_1`, and the durations
----------------------------------
Every pulse is a hard block pulse at **one** :math:`B_1`, and every duration is a whole number of
base durations -- 1 : 2 : 3 : 4 for 90 : 180 : 270 : 360.  The base is quantised onto the RF
raster, so the shared amplitude is exact for any requested ``refocus_duration_s`` rather than only
for convenient ones; :attr:`T2Prep.refocus_duration_s` reports what was realised.

Scope
-----
Hard-pulse ``mlev4-hard`` only.  Not covered: adiabatic BIR-4/AHP tip pairs; a two-composite
train or an arbitrary refocusing count; the :math:`B_1`-robust variant whose rephasing lobe lands
inside the *imaging* sequence; and diffusion preparation.

References
----------
Bernstein, King & Zhou, *Handbook of MRI Pulse Sequences*, Elsevier 2004, §17.4 -- the family,
its non-selectivity, and the ``+90`` / refocusing / ``-90`` / spoiler structure.
Levitt & Freeman, *J. Magn. Reson.* **43**, 65 (1981) -- composite pulses and the MLEV cycles.
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

#: The pulse whose centre is the composite's refocusing instant: the 180 it is built around.
_REFOCUS_CENTRE_INDEX = 1

#: The tip-down.  A single hard 90 about +x.
_TIP_DOWN: tuple[tuple[float, float], ...] = ((90.0, 0.0),)

#: The tip-up forms, and for each the index of the pulse whose centre ends ``prep_time_s``.
#:
#: ``'simple'`` is a single hard -90, the tip-down undone.  ``'composite_270_360'`` is 270 about
#: +x then 360 about -x, also a net -90.  MLEV-4 governs the refocusing train's phase cycling and
#: does not reach the tip-up, so both forms carry the same train.
#:
#: The index is where the preparation ends: the only pulse for the simple form, and the **270**
#: for the composite one.
_TIP_UPS: dict[str, tuple[tuple[tuple[float, float], ...], int]] = {
    'simple': (((90.0, 180.0),), 0),
    'composite_270_360': (((270.0, 0.0), (360.0, 180.0)), 0),
}


class T2Prep(Module):
    r"""
    A T2-weighting preparation: ``+90``, a refocusing train, a tip-up, and a spoiler.

    Parameters
    ----------
    opts
        The scanner.
    prep_time_s
        **The physical target.**  The interval the magnetisation spends transverse, between the
        tip-down's RF centre and the tip-up's.  Every pulse group starts on the gradient raster,
        so the realised value is quantised; :attr:`prep_time_s` reports it.
    tip_up
        Which tip-up to build.

        ``'simple'``
            A single hard ``-90``.  **The default**, and the shorter of the two.
        ``'composite_270_360'``
            ``270`` about +x then ``360`` about -x, also a net ``-90``.  About 2.6 ms more RF, and
            no quantitative robustness is claimed for it.

        The choice moves where ``prep_time_s`` ends; see :meth:`time_to_tip_up`.
    refocus_duration_s
        Duration of the 180 in each composite.  It sets the :math:`B_1` amplitude every pulse in
        the block shares, so a shorter one is a stronger pulse: halve it and peak :math:`B_1`
        doubles.  Rounded up so that half of it lands on the RF raster, which is what keeps the
        1 : 2 : 3 : 4 duration ratios -- and therefore the single :math:`B_1` -- exact.
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
    tip_up : str
        The tip-up in use.
    refocus_duration_s, requested_refocus_duration_s : float
        The realised 180 duration and the request it came from.
    min_prep_time_s : float
        The shortest preparation the pulse durations and the raster allow.  A shorter request
        raises.
    b1_hz : float
        The single amplitude every pulse in the block is played at.
    pulses : tuple
        Every RF event, in time order: the tip-down, four composites of three, then the tip-up --
        fourteen with ``'simple'`` and fifteen with ``'composite_270_360'``.
    spoiler : LogicBlock
        The crusher, placed after the tip-up.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
    ...                rf_dead_time=100e-6, rf_ringdown_time=30e-6)
    >>> prep = sc.modules.T2Prep(opts=opts, prep_time_s=50e-3, spoil_voxel_mm=4.0)
    >>> prep.tip_up
    'simple'
    >>> len(prep.pulses)
    14
    >>> prep.time_to_tip_up() - prep.time_to_tip_down() == prep.prep_time_s
    True
    >>> round(prep.prep_time_s * 1e3, 6)
    50.0

    The other tip-up is one argument, and it is longer:

    >>> composite = sc.modules.T2Prep(opts=opts, prep_time_s=50e-3, spoil_voxel_mm=4.0,
    ...                               tip_up='composite_270_360')
    >>> len(composite.pulses)
    15
    >>> composite().duration > prep().duration
    True
    """

    def __init__(
        self,
        *,
        opts: Opts,
        prep_time_s: float,
        tip_up: str = 'simple',
        refocus_duration_s: float = 1e-3,
        spoil_cycles_per_voxel: float = 4.0,
        spoil_axis: str = 'z',
        spoil_voxel_mm: float,
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.requested_prep_time_s = require_positive(prep_time_s, 'prep_time_s')
        self.tip_up = self._check_tip_up(tip_up)
        require_positive(refocus_duration_s, 'refocus_duration_s')
        self._rf_raster = Raster(float(opts.rf_raster_time), 'rf_raster_time')
        self._grad_raster = Raster(float(opts.grad_raster_time), 'grad_raster_time')

        # **Every duration is a whole number of base durations**, and the base is a 90.  A hard
        # pulse's flip angle is 2 pi B1 tau, so at one B1 the durations are in the ratio of the
        # angles: 1 : 2 : 3 : 4 for 90 : 180 : 270 : 360.  Quantising the *base* onto the RF
        # raster keeps every one of them exact; quantising each pulse separately would not, and
        # the shared B1 would then hold only for durations that happen to divide evenly.
        self._base_duration_s = float(self._rf_raster.ceil(refocus_duration_s / 2.0))
        self.refocus_duration_s = 2.0 * self._base_duration_s
        self.requested_refocus_duration_s = float(refocus_duration_s)
        self.b1_hz = 0.25 / self._base_duration_s

        require_usable_max_b1(self.opts, described=self._described())
        self._tip_down = self._group(_TIP_DOWN)
        self._composites = tuple(self._group(_COMPOSITE_REFOCUS, inverted=flip)
                                 for flip in _MLEV4_INVERTED)
        angles, self._tip_up_index = _TIP_UPS[self.tip_up]
        self._tip_up = self._group(angles)
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
            - self._offset_to_centre(self._tip_up, self._tip_up_index)))
        self.prep_time_s = self.time_to_tip_up() - self.time_to_tip_down()
        self._composite_starts = tuple(
            float(self._grad_raster.nearest(
                self.time_to_tip_down() + fraction * self.prep_time_s
                - self._offset_to_centre(group, _REFOCUS_CENTRE_INDEX)))
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
        return self._offset_to_centre(self._tip_down, 0)

    def time_to_tip_up(self) -> float:
        """
        Seconds from the block's start to the instant ``prep_time_s`` ends.

        ==============================  =========================================
        ``tip_up='simple'``             the centre of the single -90
        ``tip_up='composite_270_360'``  the centre of the **270**
        ==============================  =========================================
        """
        return self._tip_up_start + self._offset_to_centre(self._tip_up, self._tip_up_index)

    def time_to_refocus(self, index: int) -> float:
        """
        Seconds from the block's start to composite refocusing pulse `index`'s effective centre.

        As realised, not as intended: the group's start is quantised onto the gradient raster, so
        this is within half a raster of ``fraction * prep_time_s`` rather than exactly on it.  The
        refocusing condition is that the gaps either side of each pulse are **equal**, which
        survives that quantisation, and :attr:`interval_gaps_s` is how it is checked.
        """
        return self._composite_starts[index] + self._offset_to_centre(
            self._composites[index], _REFOCUS_CENTRE_INDEX)

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
        Return the preparation: every pulse of the chosen tip-up, then the spoiler.

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
        """One hard block pulse.  Its duration is a whole number of base durations."""
        duration_s = round(flip_deg / 90.0) * self._base_duration_s
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

    def _offset_to_centre(self, group: tuple[Event, ...], index: int) -> float:
        """Seconds from a group's start to the RF centre of the pulse at `index`."""
        at = sum(self._slot_of(rf) for rf in group[:index])
        rf = group[index]
        return at + float(rf.delay) + float(pp.calc_rf_center(rf)[0])

    def _shortest_preparation(self) -> float:
        """
        Return the shortest ``prep_time_s`` the pulses leave room for, in seconds.

        Three things can bind, and which one does depends on the pulse durations rather than on
        the preparation time, so all three are solved rather than searched:

        * the first composite must start after the tip-down ends;
        * consecutive composites are ``prep_time_s / 4`` apart and must not overlap;
        * the last composite must end before the tip-up begins.

        The two outer constraints are measured from the composite's **refocusing instant** -- the
        centre of its 180 -- to each end of the group, and those two distances are **not equal**.
        A group is a run of raster-ceiled slots and every pulse carries a dead time in front and a
        ringdown behind, so the refocusing instant is the group's midpoint only when those times
        happen to make it one.
        """
        tip_down = self._duration_of(self._tip_down)
        composite = self._duration_of(self._composites[0])
        lead = self._offset_to_centre(self._composites[0], _REFOCUS_CENTRE_INDEX)
        trail = composite - lead
        first, last = _REFOCUS_FRACTIONS[0], _REFOCUS_FRACTIONS[-1]
        after_tip_down = (tip_down - self.time_to_tip_down() + lead) / first
        between = composite / (_REFOCUS_FRACTIONS[1] - first)
        before_tip_up = (trail + self._offset_to_centre(self._tip_up, self._tip_up_index)
                         ) / (1.0 - last)
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

    def _check_tip_up(self, tip_up: str) -> str:
        """Return `tip_up` having checked it names a tip-up this module has."""
        if tip_up in _TIP_UPS:
            return tip_up
        listed = ', '.join(repr(name) for name in _TIP_UPS)
        msg = format_error(
            f'tip_up must be one of {listed}, got {tip_up!r}.',
            {'tip_up': tip_up},
            ["'simple' is a single -90, and the default",
             "'composite_270_360' is 270 about +x then 360 about -x, also a net -90"],
        )
        raise ConfigurationError(msg)

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
