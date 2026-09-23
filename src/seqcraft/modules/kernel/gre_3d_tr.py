r"""
:class:`GRE3DTR` -- one repetition of a 3D Cartesian gradient echo.

A sibling of :class:`~seqcraft.modules.GRE2DTR`: the two compose the same leaves independently,
and a slice-selective 2D acquisition is not the centre partition of a 3D slab.

The z axis
----------
On x and y a 3D repetition is a 2D one.  On z it is not, and there are two cases:

**Non-selective excitation** -- what every official Pulseq 3D reference does
(``writeGradientEcho3D.m``, the MPRAGE variants, ``write_3Dt1_mprage.py``, all block pulses).  The
z axis carries a partition encode and nothing else.

**Slab-selective excitation.**  Now two moments land on one axis inside one window: the rephasing
the slab selection implies, and the partition encoding.  Played one after the other they cost two
windows of echo time; added together and played as a single gradient they cost one, which is what
this module does.

So this module solves

.. math:: A_z(p) = A_{\text{slab}} + A_{\text{partition}}(p)

for every partition and realises it as **one** z winder of fixed duration, which is what keeps TE
constant across the volume.  It owns that jointly because neither leaf can:
:class:`~seqcraft.modules.Excitation` knows the rephasing its own slab implies and nothing about
partitions, and :class:`~seqcraft.modules.PhaseEncode` knows the moment a partition index wants
and nothing about a slab.

Both terms are signed, and that matters
----------------------------------------
The limiting partition is a *result* of the coupling, never a property of the index::

    A_slab = -120   partitions -200 | 0 | +200   ->  A_z  -320 | -120 |  +80   low edge limits
    A_slab = +120                                ->  A_z   -80 | +120 | +320   high edge limits

Reversing the slab gradient's polarity moves the limit to the other edge of k-space.  So nothing
here assumes the last partition is worst, or that the largest ``|A_partition|`` before combining
is worst, or that a larger index means a more positive ``kz`` -- the Pulseq 3D demo reverses its z
lattice relative to y for a vendor reconstruction, so index order is not even a reliable guide to
sign.  The sign comes from :meth:`~seqcraft.modules.PhaseEncode.k_per_m`, which is semantics.

Every supported partition is enumerated.  :math:`A_z(p)` is affine in ``p``, so on symmetric
hardware with a contiguous table the extreme is at one of the two edges -- but the shortcut buys
nothing at real matrix sizes and stops being true for partial ``kz`` or a supplied partition
table.

One winder window, shared by every partition
--------------------------------------------
Every partition uses **one** winder duration.  Letting each take its own shortest would make TE a
function of ``kz``: a contrast gradient across the volume that no k-space check would show and no
reconstruction expects.

So when the limiting partition needs more time, the shared window is lengthened for **every**
partition, and TE stays constant across the volume.

With ``te_s=None`` and ``tr_s=None`` the result is the shortest legal design.  An explicit request
below the achievable minimum raises, naming the partition responsible.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pypulseq as pp

from ...design.events import derive
from ...design.logic import LogicBlock
from ...design.module import Module
from ...errors import ConfigurationError, format_error
from .._support import ceil_raster, require_positive
from ..encoding.phase_encoding import PhaseEncode
from ..readout.cartesian_line import CartesianLine
from ..rf.excitation import Excitation
from ..spoiler import spoiler

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from pypulseq.opts import Opts

__all__ = ['GRE3DTR']

#: The non-selective pulse duration, when the caller does not give one.  A hard pulse is short by
#: definition, and ``Excitation``'s own 3 ms default is a *shaped selective* pulse's default -- it
#: would spend three milliseconds of echo time on a pulse that selects nothing.  0.2 ms is what
#: ``writeGradientEcho3D.m`` uses, and it is short enough that the RF is not the term that sets TE.
_HARD_PULSE_S = 0.2e-3

#: Moments below this are treated as zero rather than handed to ``make_trapezoid``, which cannot
#: design a gradient of no area.  One part in 10^9 of a millimetre-scale ``dk``.
_NEGLIGIBLE_PER_M = 1e-9


def _limiting_index(areas: Sequence[float], opts: Opts, axis: str = 'z') -> tuple[int, float]:
    """
    Return the index of the largest-demand moment in `areas`, and the duration it needs.

    Pure, and separate from the class, because it is the part worth testing on its own: given a
    list of **signed** moments it enumerates all of them and reports which one needs the longest
    legal gradient.  A test can hand it a slab term of either sign and watch the answer move from
    one end of the list to the other, which is the behaviour a "the last partition is worst"
    shortcut would get wrong exactly when the slab term is large enough to matter.

    Zero-area entries are skipped rather than passed to ``make_trapezoid``, which cannot design a
    gradient of no area.
    """
    worst_index, worst_s = 0, 0.0
    for index, area in enumerate(areas):
        if abs(area) <= _NEGLIGIBLE_PER_M:
            continue
        duration = float(pp.calc_duration(
            pp.make_trapezoid(channel=axis, area=area, system=opts)))
        if duration > worst_s:
            worst_index, worst_s = index, duration
    return worst_index, worst_s


class GRE3DTR(Module):
    """
    One repetition of a 3D Cartesian gradient echo: excitation, encoding, readout, rewind, spoil.

    Parameters
    ----------
    opts
        The scanner.
    fov_mm
        ``(x, y, z)`` field of view, millimetres.  The z component is the **encoded** extent, and
        ``dkz = 1/fov_z`` follows from it.
    matrix
        ``(nx, ny, nz)`` -- readout samples, phase-encode lines, partitions.
    slab_thickness_mm
        ``None`` excites non-selectively, which is what every official Pulseq 3D reference does.
        A thickness excites a slab and **selects the selective path**, where the rephasing it
        implies is solved together with the partition encoding.

        A physical quantity rather than a ``slab_selective=True`` flag, mirroring
        :class:`~seqcraft.modules.Excitation`'s own ``thickness_mm``: two arguments that can
        contradict each other are two chances to be wrong.

        **Independent of ``fov_mm[2]``, in both directions.**  The only requirement is that it is
        positive.  A slab smaller than the encoded FOV is a calibration or inner-volume
        acquisition; a slab larger than it buys transition band and alias protection without
        changing the reconstructed geometry, and a conventional full-volume protocol usually wants
        at least the nominal FOV and often somewhat more.  All three are legitimate, so none is
        enforced.
    flip_deg
        Flip angle, degrees.
    rf_pulse
        The pulse shape, forwarded to :class:`~seqcraft.modules.Excitation`.  ``None`` lets the
        **excitation mode choose**: ``'sinc'`` for a slab, and ``'block'`` when there is no slab.

        That second half matters.  A non-selective excitation with a shaped pulse is the worst of
        both -- it spends a soft pulse's duration and selects nothing -- and every official Pulseq
        3D reference uses a hard block pulse for exactly this reason.  An advanced caller who
        wants a shaped but spatially non-selective pulse, for a spectrally selective excitation
        say, can still ask for one here.
    rf_duration_s, rf_time_bw_product
        Forwarded to :class:`~seqcraft.modules.Excitation`.  ``None`` defers to its default --
        **a number here would be a second default that can drift from the first** -- except for a
        block pulse, where this module supplies :data:`_HARD_PULSE_S` because ``Excitation``'s
        default is a shaped pulse's.  A sharper slab profile is a deliberate choice: raise
        `rf_time_bw_product` explicitly when the transition band matters.
    te_s, tr_s
        ``None`` for the shortest achievable.  A request below what the design can reach raises,
        naming what is limiting.
    bandwidth_hz_px, partial_fourier
        Forwarded to :class:`~seqcraft.modules.CartesianLine`.
    spoil_cycles_per_voxel, spoil_axis
        The tail spoiler, as :func:`~seqcraft.modules.spoiler` takes it.  ``('x',)`` by default,
        where :class:`~seqcraft.modules.GRE2DTR` uses ``('x', 'z')`` -- **because in 3D the z axis
        already has a gradient in the tail**: the partition rewinder.  A spoiler beside it sums
        with it, and two gradients each designed at the full slew limit do not add up to a legal
        one.  ``writeGradientEcho3D.m`` spoils on x for the same reason, its tail being
        ``gyReph, gzReph, gxSpoil``.  Asking for ``'z'`` here is allowed and may need a derated
        design; see :func:`seqcraft.opts.derate`.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    exc : Excitation
    ro : CartesianLine
    pe : PhaseEncode
    pe_z : PhaseEncode
        The leaves.  ``pe_z`` is used for its **semantics** -- the centre convention, ``dk``, the
        signed moment a partition index wants and the range check -- while the z gradient this
        module emits is its own, because the moment it has to realise is not the one ``pe_z``
        would.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=20, grad_unit='mT/m', max_slew=120, slew_unit='T/m/s',
    ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> tr = GRE3DTR(opts=o, fov_mm=(200.0, 200.0, 160.0), matrix=(64, 64, 64))
    >>> tr.center_line, tr.center_partition
    (32, 32)
    >>> round(tr.dk_per_m('z'), 4)
    6.25
    >>> slab = GRE3DTR(opts=o, fov_mm=(200.0, 200.0, 160.0), matrix=(64, 64, 64),
    ...                slab_thickness_mm=180.0)
    >>> slab.selective, round(slab.slab_rephase_area_per_m, 3)
    (True, -11.148)
    """

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: tuple[float, float, float],
        matrix: tuple[int, int, int],
        slab_thickness_mm: float | None = None,
        flip_deg: float = 15.0,
        rf_duration_s: float | None = None,
        rf_time_bw_product: float | None = None,
        rf_pulse: str | None = None,
        te_s: float | None = None,
        tr_s: float | None = None,
        bandwidth_hz_px: float = 200.0,
        partial_fourier: float = 1.0,
        spoil_cycles_per_voxel: float = 4.0,
        spoil_axis: str | Iterable[str] = ('x',),
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.fov_mm = self._require_triple(fov_mm, 'fov_mm')
        self.matrix = tuple(int(n) for n in self._require_triple(matrix, 'matrix'))
        self.selective = slab_thickness_mm is not None
        self.slab_thickness_mm = (
            require_positive(slab_thickness_mm, 'slab_thickness_mm') if self.selective else None
        )
        self.spoil_cycles_per_voxel = require_positive(
            spoil_cycles_per_voxel, 'spoil_cycles_per_voxel')
        self.spoil_axis = (
            (spoil_axis,) if isinstance(spoil_axis, str) else tuple(spoil_axis)
        )

        self.exc = Excitation(opts=opts, flip_deg=flip_deg, thickness_mm=self.slab_thickness_mm,
                              **self._resolve_pulse(rf_pulse, rf_duration_s, rf_time_bw_product))

        nx, ny, nz = self.matrix
        self.pe = PhaseEncode(opts=opts, fov_mm=self.fov_mm[1], matrix=ny, axis='y')
        self.pe_z = PhaseEncode(opts=opts, fov_mm=self.fov_mm[2], matrix=nz, axis='z')
        probe_ro = CartesianLine(opts=opts, fov_mm=self.fov_mm[0], matrix=nx, axis='x',
                                 bandwidth_hz_px=bandwidth_hz_px,
                                 partial_fourier=partial_fourier)

        # The one quantity that makes this a kernel: the slab's rephasing requirement, read as a
        # signed moment rather than as somebody's event.  Zero when non-selective, which is what
        # makes the non-selective path the same code rather than a second algorithm.
        self.slab_rephase_area_per_m = self.exc.rephaser_area_per_m

        self._z_worst_partition, z_worst_s = self._solve_z()
        self.winder_s = ceil_raster(
            max(probe_ro.prephaser_duration_s, self.pe.min_duration_s,
                self.exc.rephaser_duration_s if not self.selective else 0.0,
                z_worst_s),
            opts.grad_raster_time,
        )

        self.ro = CartesianLine(opts=opts, fov_mm=self.fov_mm[0], matrix=nx, axis='x',
                                bandwidth_hz_px=bandwidth_hz_px,
                                partial_fourier=partial_fourier,
                                prephaser_duration_s=self.winder_s)
        self.pe = PhaseEncode(opts=opts, fov_mm=self.fov_mm[1], matrix=ny, axis='y',
                              duration_s=self.winder_s)
        # One z design at the winder duration, carrying the largest combined moment, scaled per
        # partition -- the same shape PhaseEncode uses, so every partition's gradient has the same
        # ramps and only its amplitude differs.
        self._gz_reference_area_per_m = self.combined_z_area_per_m(self._z_worst_partition)
        self._gz_combined = (
            pp.make_trapezoid(channel='z', area=self._gz_reference_area_per_m,
                              duration=self.winder_s, system=opts)
            if abs(self._gz_reference_area_per_m) > _NEGLIGIBLE_PER_M else None
        )

        self.spoilers = {
            axis: spoiler(opts, cycles_per_voxel=self.spoil_cycles_per_voxel,
                          voxel_mm=self.voxel_mm(axis), axis=axis)
            for axis in self.spoil_axis
        }

        self._winder_start_s = ceil_raster(
            max(self.exc.time_to_rephaser(), float(pp.calc_duration(self.exc.rf))),
            opts.grad_raster_time,
        )
        self._ro_duration_s = self.ro().duration
        self._tail_s = max(
            self.pe(line=0, rewind=True).duration,
            self.pe_z(line=0, rewind=True).duration,
            *(block.duration for block in self.spoilers.values()),
        )
        self._te_fill_s = self._resolve_te(te_s)
        self._tr_s = self._resolve_tr(tr_s)

    # ------------------------------------------------------------------ what it knows
    @property
    def center_line(self) -> int:
        """The phase-encode line that encodes ``ky = 0``."""
        return self.pe.center_line

    @property
    def center_partition(self) -> int:
        """The partition that encodes ``kz = 0``.  ``matrix[2] // 2``, as for a line."""
        return self.pe_z.center_line

    def dk_per_m(self, axis: str) -> float:
        """Sample spacing along `axis`, 1/m.  On z this is ``1/fov_z``, from the encoded extent."""
        return {'x': self.ro.dk_per_m, 'y': self.pe.dk_per_m, 'z': self.pe_z.dk_per_m}[axis]

    def voxel_mm(self, axis: str) -> float:
        """The voxel dimension along `axis`.  On z the **partition** thickness, not the slab."""
        return {'x': self.fov_mm[0] / self.matrix[0], 'y': self.fov_mm[1] / self.matrix[1],
                'z': self.fov_mm[2] / self.matrix[2]}[axis]

    def combined_z_area_per_m(self, partition: int) -> float:
        """
        The signed z moment this repetition must realise for `partition`, 1/m.

        ``A_slab + A_partition(p)``, and the whole point of the class.  ``A_slab`` is zero for a
        non-selective excitation, so the non-selective case is this same expression rather than a
        second code path.

        Examples
        --------
        >>> import pypulseq as pp
        >>> from pypulseq.opts import Opts
        >>> o = Opts(max_grad=20, grad_unit='mT/m', max_slew=120, slew_unit='T/m/s',
        ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
        >>> tr = GRE3DTR(opts=o, fov_mm=(200.0, 200.0, 160.0), matrix=(64, 64, 64))
        >>> tr.combined_z_area_per_m(tr.center_partition)            # non-selective, centre
        0.0
        >>> round(tr.combined_z_area_per_m(tr.center_partition + 1), 4)
        6.25
        """
        return self.slab_rephase_area_per_m + self.pe_z.k_per_m(partition)

    @property
    def limiting_partition(self) -> int:
        """
        The partition whose combined z moment needs the longest gradient.

        **Found by enumerating every partition after the signed combination**, not by assuming an
        edge: which edge it is depends on the sign of the slab's rephasing, and reversing the slab
        gradient moves it to the other one.
        """
        return self._z_worst_partition

    @property
    def min_te_s(self) -> float:
        """The shortest echo time this design achieves, seconds."""
        # Not ``+ winder_s``: the readout block *starts* with its prephaser, which was designed
        # at exactly the winder duration, so ``ro.time_to_echo()`` already carries it.  The
        # phase encode and the combined z winder play beside it on their own axes.
        return self._winder_start_s + self.ro.time_to_echo() - self.exc.time_to_center()

    @property
    def te_s(self) -> float:
        """The echo time this repetition achieves, seconds."""
        return self.min_te_s + self._te_fill_s

    @property
    def min_tr_s(self) -> float:
        """The shortest repetition time this design achieves, seconds."""
        return ceil_raster(self._tail_start_s + self._tail_s, self.opts.block_duration_raster)

    @property
    def tr_s(self) -> float:
        """The repetition time this block occupies.  The block is exactly this long."""
        return self._tr_s

    def time_to_echo(self) -> float:
        """Seconds from the start of this module's block to ``k = 0``."""
        return self._readout_start_s + self.ro.time_to_echo()

    # ------------------------------------------------------------------------- assembly
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
        Return one repetition, as a block exactly :attr:`tr_s` long.

        Parameters
        ----------
        line, partition
            Zero-based phase-encode and partition indices.  Which ones are acquired, and in what
            order, is the caller's.
        phase_deg
            RF carrier phase, degrees.  The receiver is given the same value.
        acquire
            ``False`` drops the ADC and the labels and changes nothing else.
        center_mm
            ``(x, y, z)`` centre of the imaging volume, millimetres.  ``y`` must be ``0.0``.
        """
        x, y, z = (float(v) for v in center_mm)
        if y != 0.0:
            msg = format_error(
                f'center_mm y = {y:g} mm is out of scope: a phase-encode shift is a per-line ADC '
                f'phase, not a per-TR frequency.',
                {'center_mm': center_mm},
                ['shift the slab with z, or the readout FOV with x'],
            )
            raise ConfigurationError(msg)
        self.pe_z.k_per_m(partition)                      # range check, with pe_z's own message

        start = self._readout_start_s
        tail = self._tail_start_s
        # rephase=False when selective: this module realises that requirement inside the z
        # winder below, and a standalone rephaser beside it would apply the slab term twice.
        out = (
            LogicBlock()
            .add(0.0, self.exc(phase_deg=phase_deg, position_mm=z, rephase=not self.selective))
            .add(start, self.pe(line=line))
            .add(start, self.ro(acquire=acquire, phase_deg=phase_deg, offset_mm=x))
            .add(tail, self.pe(line=line, rewind=True))
            # Only the partition encoding is unwound: the slab rephasing was a one-way
            # compensation and is already spent.
            .add(tail, self.pe_z(line=partition, rewind=True))
        )
        winder = self._z_winder(partition)
        if winder is not None:
            out.add(start, winder)
        for block in self.spoilers.values():
            out.add(tail, block)
        if acquire:
            out.add(start, pp.make_label(type='SET', label='LIN', value=int(line)))
            out.add(start, pp.make_label(type='SET', label='PAR', value=int(partition)))
        fill = ceil_raster(self.tr_s - tail - self._tail_s, self.opts.block_duration_raster)
        if fill > 1e-9:
            out.add(tail + self._tail_s, pp.make_delay(fill))
        return out

    def _z_winder(self, partition: int):
        """The one z gradient carrying ``A_slab + A_partition(p)``, or ``None`` when it is zero."""
        area = self.combined_z_area_per_m(partition)
        if self._gz_combined is None or abs(area) <= _NEGLIGIBLE_PER_M:
            return None
        scale = area / self._gz_reference_area_per_m
        return derive(pp.scale_grad(self._gz_combined, scale))

    # --------------------------------------------------------------------------- design
    def _solve_z(self) -> tuple[int, float]:
        """
        Return the limiting partition and the shortest z winder that serves **every** partition.

        Enumerated, in the only order that is correct: signed partition moment, combined with the
        slab's, *then* the duration that moment needs.  Taking the worst partition first and
        adding the slab term afterwards finds the wrong edge whenever the slab term is large
        enough to move the maximum, which is exactly when it matters.
        """
        combined = [self.combined_z_area_per_m(p) for p in range(self.matrix[2])]
        return _limiting_index(combined, self.opts)

    # ------------------------------------------------------------------------- timing
    @property
    def _readout_start_s(self) -> float:
        """When the winder -- readout prephaser, phase encode, combined z -- starts."""
        return self._winder_start_s + self._te_fill_s

    @property
    def _tail_start_s(self) -> float:
        """When the rewinders and the spoilers start."""
        return self._readout_start_s + self._ro_duration_s

    def _resolve_te(self, te_s: float | None) -> float:
        """Return the delay that lengthens TE to the request, or zero for the shortest."""
        if te_s is None:
            return 0.0
        wanted = require_positive(te_s, 'te_s')
        if wanted < self.min_te_s - 1e-12:
            msg = format_error(
                f'te_s = {wanted * 1e3:.3f} ms is shorter than this repetition can achieve.',
                {
                    'te_s': wanted,
                    'min_te_s': self.min_te_s,
                    'winder_s': self.winder_s,
                    'limiting_partition': self.limiting_partition,
                    'combined_z_area_per_m': self.combined_z_area_per_m(self.limiting_partition),
                    'slab_rephase_area_per_m': self.slab_rephase_area_per_m,
                },
                [
                    f'pass te_s >= {self.min_te_s:.6g}',
                    'or te_s=None for the shortest achievable echo time',
                    'a higher bandwidth_hz_px shortens the readout, and with it min_te_s',
                    'fewer partitions, or a larger fov_mm[2], shorten the combined z winder',
                ],
            )
            raise ConfigurationError(msg)
        return ceil_raster(wanted - self.min_te_s, self.opts.grad_raster_time)

    def _resolve_tr(self, tr_s: float | None) -> float:
        """Return the repetition time this block will occupy, seconds."""
        if tr_s is None:
            return self.min_tr_s
        wanted = require_positive(tr_s, 'tr_s')
        if wanted < self.min_tr_s - 1e-12:
            msg = format_error(
                f'tr_s = {wanted * 1e3:.3f} ms is shorter than this repetition can achieve.',
                {'tr_s': wanted, 'min_tr_s': self.min_tr_s, 'te_s': self.te_s,
                 'limiting_partition': self.limiting_partition},
                [
                    f'pass tr_s >= {self.min_tr_s:.6g}',
                    'or tr_s=None for the shortest achievable repetition time',
                    'a shorter te_s shortens min_tr_s with it',
                ],
            )
            raise ConfigurationError(msg)
        return ceil_raster(wanted, self.opts.block_duration_raster)

    def _resolve_pulse(self, rf_pulse: str | None, rf_duration_s: float | None,
                       rf_time_bw_product: float | None) -> dict[str, float | str]:
        """
        Return the excitation arguments, letting the mode choose the pulse when nothing was asked.

        A slab wants a shaped pulse; no slab wants a hard one.  Defaulting both to ``'sinc'``
        would give the non-selective path a three-millisecond pulse that selects nothing and puts
        every millisecond of it into TE.
        """
        pulse = rf_pulse if rf_pulse is not None else ('sinc' if self.selective else 'block')
        out: dict[str, float | str] = {'pulse': pulse}

        if rf_duration_s is not None:
            out['duration_s'] = rf_duration_s
        elif pulse == 'block':
            out['duration_s'] = _HARD_PULSE_S

        if rf_time_bw_product is not None:
            if pulse == 'block':
                msg = format_error(
                    'a time-bandwidth product needs a shaped pulse, and this excitation is a hard '
                    'block pulse.',
                    {'rf_time_bw_product': rf_time_bw_product, 'rf_pulse': pulse,
                     'slab_thickness_mm': self.slab_thickness_mm},
                    ["pass rf_pulse='sinc' for a shaped pulse",
                     'or drop rf_time_bw_product, which a block pulse has no use for'],
                )
                raise ConfigurationError(msg)
            out['time_bw_product'] = rf_time_bw_product
        return out

    @staticmethod
    def _require_triple(value: tuple[float, float, float], name: str) -> tuple[float, ...]:
        """Return `value` as a 3-tuple, refusing the 2-tuple a 2D kernel would take."""
        try:
            out = tuple(float(v) for v in value)
        except TypeError:
            out = ()
        if len(out) != 3:
            msg = format_error(
                f'{name} must have three components for a 3D acquisition.',
                {name: value},
                ['(x, y, z) -- the z component is the encoded extent, not the excited slab'],
            )
            raise ConfigurationError(msg)
        return out
