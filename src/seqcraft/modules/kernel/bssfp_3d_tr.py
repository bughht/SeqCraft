"""
:class:`bSSFP3DTR` -- one balanced repetition of a 3D Cartesian bSSFP acquisition.

A **sibling** of :class:`~seqcraft.modules.bSSFP2DTR` rather than a wrapper around it, for the
reason :class:`~seqcraft.modules.GRE3DTR` is a sibling of :class:`~seqcraft.modules.GRE2DTR`: the
``z`` axis does a different job here.  In 2D it carries slice selection and its balancing.  Here
it carries slab selection, its balancing, **and** the partition encoding -- and the way those
three share one axis is the whole content of this module.

What makes it balanced
----------------------
The same definition as in 2D, and it is not weakened by the extra axis::

    TE               = echo - RF effective centre of repetition n
    TR               = RF effective centre of n+1 - RF effective centre of n
    balance interval = [RF effective centre of n, RF effective centre of n+1]

    net area on x, y and z over that interval = 0

``TE = TR/2`` remains the standard symmetric realisation and the default, not the definition.

Each repetition discharges its own partition encoding
-----------------------------------------------------
This is the part 3D adds, and it is a **compositional** requirement rather than an arithmetic
one.  The balance interval ends inside the *next* repetition, so a decomposition can satisfy
``M0 = 0`` for a train of identical repetitions and still be wrong::

    tail(n) = -A_post(n) + K(p_n)        # balance + encode, before the echo
    lead(n) = -A_pre(n)  - K(p_n)        # balance + rewind, before the NEXT pulse

    interval = A_post(n) + tail(n) + lead(n+1) + A_pre(n+1)
             = K(p_n) - K(p_{n+1})

That is zero only while ``p_n == p_{n+1}``.  A 3D acquisition chooses its view order, so a
repetition whose balance depends on **which partition its successor happens to acquire** is not
a repetition at all.  The realisation here instead gives the encoding and its rewind to the same
repetition::

    lead(n)   = -A_pre(n)                before RF n      -- no partition dependence at all
    tail(n)   = -A_post(n) + K(p_n)      after RF n       -- selection balance + encode
    rewind(n) = -K(p_n)                  after the echo   -- this repetition's own rewind

    interval = A_post(n) + tail(n) + rewind(n) + lead(n+1) + A_pre(n+1)
             = A_post - A_post + K(p_n) - K(p_n) - A_pre + A_pre
             = 0                                           -- for any p_{n+1}

``tests/modules/test_bssfp_3d_tr.py`` asserts that on compiled waveforms with every neighbouring
pair changing **both** indices, which is the only probe that can tell the two decompositions
apart.

Why the encoding is lumped into the balancing lobe
--------------------------------------------------
``tail`` carries one gradient of ``-A_post + K(p)``, not a balancing lobe beside an encoding
blip.  The z axis is already occupied between the pulse and the echo, and a second independent z
gradient in that window merges with the first; at the edge of the partition table the merged
waveform exceeds the slew limit.  One lobe carrying the signed sum is the same moment in a legal
shape, and it is what :class:`~seqcraft.modules.GRE3DTR` does with ``A_slab + A_partition``.

The limiting partition is a *result* of the signed combination and never a property of the index:
``tail`` is extreme wherever ``-A_post + K(p)`` is, which moves to the other edge of k-space when
the slab gradient reverses.  Every supported partition is enumerated.

One window, shared by every partition
-------------------------------------
Letting each partition take its shortest lobe would make TE and TR functions of ``kz``: a
contrast gradient across the volume that no k-space check would show.  So the selection-side
window is one duration for the whole table, sized by whichever partition needs the longest --
and the leading lobe takes that same window, which is what keeps the z structure symmetric about
the echo.

**That symmetry is a choice, not part of the balance definition.**  Balance is the zero-area
condition above.  What the symmetric realisation buys is the clean limit at ``kz = 0``, where z
carries nothing but slab-selection balancing and the first moment over the interval vanishes;
away from the centre the partition encoding contributes its own first-moment term, exactly as
the phase encoding already does on ``y``.

Selective and non-selective are two excitation modes
----------------------------------------------------
``slab_thickness_mm=None`` is not "the selective design with the selection gradient removed".  It
is a different physical excitation: no selection gradient **and** a hard pulse.  That rule is
:func:`~seqcraft.modules._support.resolve_excitation_mode`, shared with
:class:`~seqcraft.modules.GRE3DTR` so the two kernels cannot answer it differently.

What this module does not own
-----------------------------
RF phase cycling, the receiver phase progression, start-up preparation, dummy repetitions,
segmentation and the view order are all the caller's, exactly as on
:class:`~seqcraft.modules.bSSFP2DTR`.  **A balanced repetition is not a repetition in steady
state**, and this module builds the first of those.

There are no spoilers, for the reason there are none in 2D: a spoiler is a deliberate non-zero
net area, which is the condition this module exists not to have.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pypulseq as pp

from ...design.events import derive
from ...design.logic import LogicBlock
from ...design.module import Module
from ...design.validation import require_axis, require_positive
from ...errors import ConfigurationError, format_error
from .._support import area_until, ceil_raster, resolve_excitation_mode
from ..encoding.phase_encoding import PhaseEncode
from ..readout.cartesian_line import CartesianLine
from ..rf.excitation import Excitation

if TYPE_CHECKING:
    from pypulseq.opts import Opts

    from ...design.events import Event

__all__ = ['bSSFP3DTR']

#: Moments below this are treated as zero rather than handed to ``make_trapezoid``, which cannot
#: design a gradient of no area.  One part in 10^9 of a millimetre-scale ``dk``.
_NEGLIGIBLE_PER_M = 1e-9


class bSSFP3DTR(Module):  # noqa: N801 -- `bSSFP` is the domain spelling; see ruff.toml
    """
    One balanced repetition of a 3D Cartesian bSSFP acquisition.

    Parameters
    ----------
    opts
        The scanner.
    fov_mm
        ``(x, y, z)`` field of view, millimetres.  ``z`` is the extent the partitions encode,
        which is not the same number as `slab_thickness_mm`: the excited slab is normally the
        thicker of the two.
    matrix
        ``(nx, ny, nz)``.
    slab_thickness_mm
        The excited slab, millimetres.  ``None`` excites **non-selectively**, which means no
        selection gradient *and* a hard pulse -- see `rf_pulse`.
    flip_deg
        Flip angle.  The bSSFP default is larger than a spoiled gradient echo's.
    rf_pulse
        ``None`` lets the excitation **mode** choose: a shaped pulse for a slab, a hard block
        pulse when there is no slab.  Those are two different physical excitations, and the
        second is not the first with its gradient deleted.  An explicit value is always honoured,
        so a caller who wants a shaped but spatially non-selective pulse -- a spectrally
        selective excitation, say -- still gets one.
    rf_duration_s, rf_time_bw_product
        Forwarded to :class:`~seqcraft.modules.Excitation`, except that a block pulse takes
        :data:`~seqcraft.modules._support.HARD_PULSE_S` rather than a shaped pulse's default.  A
        time-bandwidth product with a block pulse is refused: it has no meaning for that family.
    te_s, tr_s
        ``None`` for the symmetric realisation at its shortest.  See :attr:`min_symmetric_te_s`.
    bandwidth_hz_px
        Readout bandwidth per pixel.
    partial_fourier
        Readout partial-Fourier fraction.
    tag
        Optional label for the returned block.

    Examples
    --------
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s',
    ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> tr = bSSFP3DTR(opts=o, fov_mm=(200.0, 200.0, 160.0), matrix=(64, 64, 16),
    ...                slab_thickness_mm=100.0)
    >>> tr.center_line, tr.center_partition
    (32, 8)
    >>> round(tr.partition_area_per_m(tr.center_partition), 9)
    0.0
    >>> abs(tr.symmetry_residual_s) <= o.grad_raster_time / 2
    True
    """

    _PE_AXIS = 'y'
    _PAR_AXIS = 'z'

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: tuple[float, float, float],
        matrix: tuple[int, int, int],
        slab_thickness_mm: float | None = None,
        flip_deg: float = 35.0,
        rf_pulse: str | None = None,
        rf_duration_s: float | None = None,
        rf_time_bw_product: float | None = None,
        te_s: float | None = None,
        tr_s: float | None = None,
        bandwidth_hz_px: float = 400.0,
        partial_fourier: float = 1.0,
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        fov = self._require_triple(fov_mm, 'fov_mm')
        mat = self._require_triple(matrix, 'matrix')
        self.fov_mm = (fov[0], fov[1], fov[2])
        self.matrix = (int(mat[0]), int(mat[1]), int(mat[2]))
        self.selective = slab_thickness_mm is not None
        self.slab_thickness_mm = (
            require_positive(slab_thickness_mm, 'slab_thickness_mm') if self.selective else None
        )

        # The excitation MODE, not merely its geometry: the rule is shared with GRE3DTR so the
        # two 3D kernels cannot answer it differently.
        self.exc = Excitation(
            opts=opts, flip_deg=flip_deg, thickness_mm=self.slab_thickness_mm,
            **resolve_excitation_mode(
                selective=self.selective, rf_pulse=rf_pulse, rf_duration_s=rf_duration_s,
                rf_time_bw_product=rf_time_bw_product,
                selection=('slab_thickness_mm', slab_thickness_mm),
            ),
        )

        # The two selection halves, measured on the gradient rather than assumed equal: a
        # minimum-phase pulse's effective centre is nowhere near its midpoint.  Both are zero
        # when nothing selects, which is what makes the non-selective path the same arithmetic
        # rather than a second code path.
        self._a_pre = (0.0 if self.exc.gz is None
                       else area_until(self.exc.gz, self.exc.time_to_center()))
        self._a_post = -self.exc.rephaser_area_per_m

        raster = float(opts.grad_raster_time)
        self._grad_start_s = ceil_raster(
            max(float(pp.calc_duration(self.exc.gz)) if self.exc.gz is not None else 0.0,
                float(pp.calc_duration(self.exc.rf))),
            raster,
        )

        nx, ny, nz = self.matrix
        self.pe_z = PhaseEncode(opts=opts, fov_mm=self.fov_mm[2], matrix=nz, axis='z')

        # ONE selection-side window for the whole partition table, taken by both lobes.  Letting
        # each partition take its shortest would make TE a function of kz; letting the two lobes
        # differ would break the symmetry about the echo that the default realisation wants.
        #
        # The trailing lobe's demand is ENUMERATED over the signed combination: which partition
        # makes `-A_post + K(p)` longest depends on the sign of the slab gradient, so reversing
        # it moves the limit to the other edge of k-space.
        self._z_window_s = ceil_raster(
            max(self._lobe_duration_s(-self._a_pre),
                *(self._lobe_duration_s(-self._a_post + self.partition_area_per_m(p))
                  for p in range(nz))),
            raster,
        )
        post_half = self._grad_start_s - self.exc.time_to_center()
        pre_half = self.exc.time_to_center()
        self._z_gap_s = ceil_raster(max(pre_half - post_half, 0.0), raster)
        self._lead_gap_s = ceil_raster(max(post_half - pre_half, 0.0), raster)
        self._lead_s = self._z_window_s + self._lead_gap_s
        self._z_lead = self._lobe('z', -self._a_pre, self._z_window_s)

        probe_ro = CartesianLine(opts=opts, fov_mm=self.fov_mm[0], matrix=nx, axis='x',
                                 bandwidth_hz_px=bandwidth_hz_px,
                                 partial_fourier=partial_fourier)
        probe_pe = PhaseEncode(opts=opts, fov_mm=self.fov_mm[1], matrix=ny, axis='y')
        # One window between the pulse and the readout, and one more after it.  Axes overlap for
        # free, so each is its longest participant rather than their sum -- and the post-echo one
        # now also has to hold this repetition's own partition rewind.
        self.winder_s = ceil_raster(
            max(probe_ro.prephaser_duration_s, probe_pe.min_duration_s,
                self.pe_z.min_duration_s,
                self._lobe_duration_s(-probe_ro.area_after_echo_per_m)),
            raster,
        )

        self.ro = CartesianLine(opts=opts, fov_mm=self.fov_mm[0], matrix=nx, axis='x',
                                bandwidth_hz_px=bandwidth_hz_px, partial_fourier=partial_fourier,
                                prephaser_duration_s=self.winder_s)
        self.pe = PhaseEncode(opts=opts, fov_mm=self.fov_mm[1], matrix=ny, axis='y',
                              duration_s=self.winder_s)
        self._x_balance = self._lobe('x', -self.ro.area_after_echo_per_m, self.winder_s)
        self._ro_duration_s = self.ro().duration
        self._tail_s = self.winder_s
        self._te_s, self._tr_s = self._resolve_timing(te_s, tr_s)

    # ------------------------------------------------------------------ what it knows
    @property
    def center_line(self) -> int:
        """The phase-encode line that encodes ``ky = 0``."""
        return self.pe.center_line

    @property
    def center_partition(self) -> int:
        """The partition that encodes ``kz = 0``.  ``matrix[2] // 2``, as for a line."""
        return self.pe_z.center_line

    def partition_area_per_m(self, partition: int) -> float:
        """The signed z moment `partition` encodes, 1/m -- ``K(p)`` in the module docstring."""
        return self.pe_z.k_per_m(partition)

    def dk_per_m(self, axis: str) -> float:
        """The k-space step between neighbouring samples on `axis`, 1/m."""
        return {'x': self.ro.dk_per_m, 'y': self.pe.dk_per_m,
                'z': self.pe_z.dk_per_m}[require_axis(axis)]

    def voxel_mm(self, axis: str) -> float:
        """The voxel dimension along `axis`, millimetres."""
        return {'x': self.fov_mm[0] / self.matrix[0],
                'y': self.fov_mm[1] / self.matrix[1],
                'z': self.fov_mm[2] / self.matrix[2]}[require_axis(axis)]

    @property
    def min_te_s(self) -> float:
        """
        Shortest echo time this repetition can reach, seconds, **as an asymmetric realisation**.

        The earliest the echo can physically arrive.  Not the default, and on a symmetric
        protocol shorter than :attr:`te_s`: see :attr:`min_symmetric_te_s`.
        """
        return self._head_min_s

    @property
    def min_symmetric_te_s(self) -> float:
        """
        Shortest echo time reachable **with** ``TE = TR/2``, seconds.

        The larger of the two halves: an echo earlier than this cannot be the midpoint of any
        repetition.  This is the TE a default ``bSSFP3DTR`` achieves, and half of
        :attr:`min_tr_s`, both to within :attr:`symmetry_residual_s`.
        """
        head, _ = self._symmetric_fills()
        return self._head_min_s + head

    @property
    def symmetry_residual_s(self) -> float:
        """
        How far this repetition's echo sits from the exact midpoint, seconds.  Signed.

        ``TE = TR/2`` is the target and the default aims at it, but for a given readout geometry
        the exact midpoint need not lie on the timing lattice.  The default then takes the
        nearest symmetric realisation and reports what it achieved rather than rounding the
        declared TE.  At most half a gradient raster.
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

        The symmetric realisation is what this module defaults to, so the minimum it reports is
        that realisation's.  An explicit `te_s` asks a different question and relaxes it.
        """
        head, tail = self._symmetric_fills()
        return self._span_s(head, tail)

    @property
    def tr_s(self) -> float:
        """The repetition time, seconds: **RF effective centre to RF effective centre**."""
        return self._tr_s

    def time_to_rf_center(self) -> float:
        """
        Seconds from the start of this block to the **RF effective centre**.

        The timing origin of everything this module declares.  Read from the waveform rather
        than from the event's midpoint, which for an asymmetric pulse is a different instant.
        """
        return self._lead_s + self.exc.time_to_center()

    def time_to_echo(self) -> float:
        """Seconds from the start of this block to ``k = 0``: ``time_to_rf_center() + te_s``."""
        return self.time_to_rf_center() + self._te_s

    # ----------------------------------------------------------------------- assembly
    def build(
        self,
        *,
        line: int,
        partition: int,
        phase_deg: float = 0.0,
        acquire: bool = True,
        center_mm: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> LogicBlock:
        """
        Return one balanced repetition.

        Parameters
        ----------
        line
            Zero-based phase-encode index.
        partition
            Zero-based partition index.  **Neighbouring repetitions may use any partitions**:
            each discharges its own encoding, so balance never depends on what the successor
            acquires.  That is what makes a caller-owned view order possible.
        phase_deg
            RF carrier phase, and the receiver phase, set from the same value so the two cannot
            disagree.  The progression itself is the caller's.
        acquire
            ``False`` drops the ADC and the ``LIN``/``PAR`` labels and changes nothing else, so
            a start-up repetition loads identical gradients.  **It does not make the repetition
            a preparation**: what start-up means is the caller's.
        center_mm
            ``(x, y, z)`` centre of the imaging volume, millimetres.  ``y`` must be ``0.0``.
        """
        self._require_partition(partition)
        x, y, z = (float(v) for v in center_mm)
        if y != 0.0:
            raise ConfigurationError(format_error(
                f'center_mm y = {y:g} mm is out of scope: a phase-encode shift is a per-line '
                f'ADC phase, not a per-TR frequency.',
                {'center_mm': center_mm},
                ['shift the slab with z, or the readout FOV with x'],
            ))

        winder = self._lead_s + self._grad_start_s + self._head_fill_s
        tail = winder + self._ro_duration_s
        out = (
            LogicBlock()
            # The leading lobe first, then the excitation: this lobe belongs to the interval that
            # ENDS at the RF centre just after it, and it carries no partition term because the
            # repetition that owns the encoding has already discharged it.
            .add(0.0, self._z_lead)
            .add(self._lead_s, self.exc(phase_deg=phase_deg, position_mm=z, rephase=False))
            .add(winder, self.pe(line=line))
            .add(winder, self.ro(acquire=acquire, phase_deg=phase_deg, offset_mm=x))
            .add(tail, self.pe(line=line, rewind=True))
            .add(tail, self._x_balance)
        )
        z_tail = self._z_tail_for(partition)
        if z_tail is not None:
            out.add(self._lead_s + self._grad_start_s + self._z_gap_s, z_tail)
        rewind = self._partition_rewind(partition)
        if rewind is not None:
            out.add(tail, rewind)
        if acquire:
            out.add(winder, pp.make_label(type='SET', label='LIN', value=int(line)))
            out.add(winder, pp.make_label(type='SET', label='PAR', value=int(partition)))
        fill = self._tr_s - (tail + self._tail_s)
        if fill > 1e-9:
            out.add(tail + self._tail_s, pp.make_delay(fill))
        return out

    # ----------------------------------------------------------------------- design
    def _z_tail_for(self, partition: int) -> Event | None:
        """The one pre-echo z lobe: slab-selection balance and the partition encode together."""
        area = -self._a_post + self.partition_area_per_m(partition)
        if abs(area) < _NEGLIGIBLE_PER_M:
            return None
        return self._lobe('z', area, self._z_window_s)

    def _partition_rewind(self, partition: int) -> Event | None:
        """The post-echo lobe that returns ``kz`` to zero inside **this** repetition."""
        area = -self.partition_area_per_m(partition)
        if abs(area) < _NEGLIGIBLE_PER_M:
            return None
        return self._lobe('z', area, self.winder_s)

    def _lobe(self, channel: str, area_per_m: float, duration_s: float) -> Event:
        """One balancing lobe of `area_per_m`, stretched to the shared window."""
        return derive(pp.make_trapezoid(channel=channel, area=area_per_m,
                                        duration=duration_s, system=self.opts))

    def _lobe_duration_s(self, area_per_m: float) -> float:
        """The shortest this lobe can be, seconds -- a participant in a shared window."""
        if abs(area_per_m) < _NEGLIGIBLE_PER_M:
            return 0.0
        return float(pp.calc_duration(
            pp.make_trapezoid(channel='z', area=area_per_m, system=self.opts)))

    @property
    def _z_tail_end_s(self) -> float:
        """Seconds from the RF effective centre to the end of the trailing z lobe."""
        return (self._grad_start_s - self.exc.time_to_center()
                + self._z_gap_s + self._z_window_s)

    @property
    def _head_floor_s(self) -> float:
        """
        The smallest head fill the geometry allows, seconds.

        The winder window may start while the trailing z lobe is still playing -- different axes
        overlap for free -- but that lobe has to be **finished by the echo**, or a z gradient
        would still be running when ``kz`` is supposed to have reached its encoded value.
        """
        head_without_fill = (self._grad_start_s - self.exc.time_to_center()
                             + self.ro.time_to_echo())
        return max(self._z_tail_end_s - head_without_fill, 0.0)

    @property
    def _head_min_s(self) -> float:
        """TE at the shortest: RF centre to echo with only the fill the geometry forces."""
        return (self._grad_start_s - self.exc.time_to_center()
                + self.ro.time_to_echo() + self._head_floor_s)

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
        return (self._lead_s + self._grad_start_s + head_fill_s
                + self._ro_duration_s + self._tail_s + tail_fill_s)

    def _symmetric_fills(self) -> tuple[float, float]:
        """
        The shortest pair of fills that puts the echo at the midpoint, as closely as it can go.

        Both fills are raster multiples, so their difference moves in whole rasters while the
        quantity it has to cancel generally is not one -- hence `round` rather than `ceil`, which
        leaves the echo within half a raster of the midpoint instead of always past it.
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
        return (self._grad_start_s - self.exc.time_to_center()
                + self.ro.time_to_echo() + head), self._span_s(head, tail)

    def _symmetric_fills_within(self, tr_s: float) -> tuple[float, float]:
        """Split a requested TR into the two fills that put the echo nearest the midpoint."""
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

    def _require_partition(self, partition: int) -> None:
        """Refuse a partition outside the table, naming the centre as the useful landmark."""
        if not 0 <= int(partition) < self.matrix[2]:
            raise ConfigurationError(format_error(
                f'partition {partition} is outside the table this repetition encodes.',
                {'partition': partition, 'matrix': self.matrix},
                [f'pass 0 <= partition < {self.matrix[2]}',
                 f'the centre of k-space is partition {self.center_partition}'],
            ))

    @staticmethod
    def _require_triple(value: object, name: str) -> tuple[float, ...]:
        """Return `value` as a 3-tuple, refusing the 2-tuple a 2D kernel would take."""
        try:
            out = tuple(float(v) for v in value)  # type: ignore[union-attr]
        except TypeError:
            out = ()
        if len(out) != 3:
            raise ConfigurationError(format_error(
                f'{name} must have three components for a 3D acquisition.',
                {name: value},
                ['pass (x, y, z)',
                 'a two-component value is a 2D protocol, which bSSFP2DTR takes'],
            ))
        return out
