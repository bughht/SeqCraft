r"""
:class:`TSEShot` -- one excitation and its train of refocused Cartesian readouts.

One **shot** rather than one echo, because the excitation cannot be separated from the train
without losing arithmetic: the first interval's length, the initial dephaser's area *and* its
placement, and the shift that keeps the echo at the midpoint all couple the excitation's effective
centre to the first refocusing centre.  The layer above decides how many shots to run and what
each one samples.

What this module owns
---------------------
The **coupling between leaves** -- the values no single leaf has enough information to determine:

* :class:`~seqcraft.modules.CartesianLine` knows how asymmetric it is about its own echo
  (:attr:`~seqcraft.modules.CartesianLine.echo_moment_imbalance_per_m`) but not what window it
  will be given, so it cannot design the lobes that compensate it.
* :class:`~seqcraft.modules.Refocusing` balances its crushers about its own effective centre, on
  the selection axis, which is the only axis whose pair straddles *its* pulse.  The readout-axis
  pair straddles the pulse too, but its two halves live in two different readout blocks and their
  size depends on the readout's geometry -- which a refocusing pulse has no business knowing, and
  which would make the module wrong for SE-EPI, a diffusion spin echo or a spectroscopy pulse.
* :class:`~seqcraft.modules.PhaseEncode` knows its own shortest blip and nothing about the rest.

This module is the first layer that holds all three at once, so it designs one crusher window wide
enough for the slice crusher, the readout lobes and the phase blip, and it is where the echo
spacing that window implies is decided.  That is the same coupling ``GRE2DTR`` performs under the
name "winder coupling", and the mathematics is written here as two small pure functions in case a
third kernel ever wants it -- but it is deliberately **not** promoted to a shared abstraction on
two examples.

The three invariants a train has to hold
----------------------------------------
All three are properties of the emitted waveform, so they can be checked on a compiled sequence
rather than taken on trust.

1. **k is exact at every echo.**  A refocusing pulse conjugates k, so the gradient area between
   one echo and the next pulse's effective centre must equal the area between that centre and the
   following echo, on the readout and selection axes.  The selection axis is
   ``Refocusing``'s own business; the readout axis is this module's, and it is what
   :attr:`lobe_in` and :attr:`lobe_out` are for.  On the phase axis the requirement is different
   -- an encode is *meant* to move k -- and what must hold is that ``ky`` returns to zero at each
   pulse, which is why the blip goes after the refocusing pulse and the rewinder before the next
   one.  Putting the blip before the pulse, which a single-echo spin echo may legally do, makes
   every later pulse flip the encode.

2. **The refocusing pulses are uniformly spaced**, and

3. **every echo sits at the midpoint between its two pulses** -- because an FSE signal is a
   mixture of primary and stimulated echoes, and they coincide only there.  When a bound other
   than the train sets the spacing, the echo would otherwise drift; :attr:`shift_s` is the delay
   that prevents it, and it is not always zero.

What it does not own
--------------------
Which ``ky`` each echo acquires, how many shots there are, how many dummies precede them and which
echo carries the centre of k-space are all acquisition policy, and they belong to
:class:`~seqcraft.modules.FSE2D`.  This module takes one shot's worth of lines as data and plays
them.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pypulseq as pp

from ...design.logic import LogicBlock
from ...design.module import Module
from ...errors import ConfigurationError, format_error
from .._support import ceil_raster, require_pair, require_positive
from ..encoding.phase_encoding import PhaseEncode
from ..readout.cartesian_line import CartesianLine
from ..rf.excitation import Excitation
from ..rf.refocusing import Refocusing
from ..spoiler import spoiler

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

    from pypulseq.opts import Opts

__all__ = ['TSEShot']


def _shared_window_s(minima: Mapping[str, float], raster_s: float) -> float:
    """
    The one window several axes share: the largest of their minima, on the raster.

    Pure, and separate from the class, because nothing about it is specific to a spin-echo train
    -- it is the same statement ``GRE2DTR`` makes about its winders.  Not promoted to a shared
    helper on two examples; see the module docstring.
    """
    return ceil_raster(max(minima.values()), raster_s)


def _widened_window_s(window_s: float, *, base_period_s: float, requested_period_s: float,
                      raster_s: float) -> float:
    """
    The window that delivers `requested_period_s`, given what `window_s` already delivers.

    Widening the window by ``d`` lengthens the period by exactly ``2d``: the echo moves by ``d``
    and the distance from the refocusing centre to the trailing crusher does not depend on the
    window at all.  A request below what the window already delivers is not an error and does not
    shorten it -- the minimum stands.
    """
    return ceil_raster(window_s + max(0.0, (requested_period_s - base_period_s) / 2.0), raster_s)


class TSEShot(Module):
    """
    A 2D turbo-spin-echo shot: one excitation and `echoes` refocused Cartesian readouts.

    Parameters
    ----------
    opts
        The scanner.
    fov_mm
        Field of view, millimetres, as ``(x, y)`` or one number for both.
    matrix
        ``(nx, ny)``.
    thickness_mm
        Slice thickness, millimetres.
    echoes
        Echoes per shot -- the turbo factor.  ``1`` is a conventional spin echo, and is the same
        composition rather than a special case.
    echo_spacing_s
        Requested spacing between refocusing centres.  ``None`` takes the shortest this geometry
        allows; a request below that is ignored rather than refused, and
        :attr:`echo_spacing_s` reports what was achieved.
    tr_s
        Shot repetition time.  ``None`` is :attr:`min_tr_s`.  Shorter raises, because the next
        excitation would land inside the train.
    flip_deg, refocus_flip_deg
        Excitation and refocusing flip angles, degrees.  Reduced refocusing flips are standard for
        SAR and are the entry point for variable-flip work; the value is passed through to
        :class:`~seqcraft.modules.Refocusing`, which refuses one over ``max_b1``.
    bandwidth_hz_px
        Receiver bandwidth per pixel.
    partial_fourier
        Readout partial Fourier, as :class:`~seqcraft.modules.CartesianLine` takes it.  A long
        train with `partial_fourier` below 1 and one shot is what a vendor calls HASTE.
    excitation_duration_s, refocus_duration_s
        Pulse durations.  The refocusing default is longer than the excitation's because a 180
        needs twice a 90's peak B1 at the same shape.
    refocus_thickness_factor
        How much thicker the refocusing slab is than the excited slice.  A refocusing profile's
        edges are worse than an excitation's, and a slice refocused only in its middle loses its
        edges once per echo.
    crush_cycles_slice
        Turns of phase wound across the slice by the refocusing crusher, which is what dephases
        the FID an imperfect 180 leaves behind.
    crush_cycles_readout
        Turns wound across one in-plane voxel by the readout-axis lobe pair, **on top of** the
        compensation those lobes must carry regardless.  Default ``0.0``: the conjugation cancels
        whatever the two lobes share, so a readout crusher has no k-space consequence at all and
        is paid for in echo spacing.  The Pulseq and PyPulseq demos set the equivalent of ``1.0``.
    spoil_cycles_per_voxel, spoil_axis
        The tail spoiler, as :func:`~seqcraft.modules.spoiler` takes it.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    exc : Excitation
    refoc : Refocusing
    pe : PhaseEncode
    ro : CartesianLine
        The leaves this shot couples.  Exposed because a caller reading
        ``shot.ro.echo_sample(0)`` should not have to rebuild one to ask.
    lobe_in, lobe_out : Event
        The readout-axis pair straddling each refocusing pulse.  Their mean is
        `crush_cycles_readout` and their difference is fixed by the readout's own asymmetry about
        the echo -- see :meth:`~seqcraft.modules.CartesianLine.echo_moment_imbalance_per_m`.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s',
    ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> shot = TSEShot(opts=o, fov_mm=256.0, matrix=(128, 128), thickness_mm=5.0, echoes=4)
    >>> round(shot.echo_spacing_s * 1e3, 3)
    10.7
    >>> [round(shot.te_s(n) * 1e3, 3) for n in range(4)]
    [10.702, 21.402, 32.102, 42.802]
    """

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: float | tuple[float, float],
        matrix: tuple[int, int],
        thickness_mm: float,
        echoes: int = 1,
        echo_spacing_s: float | None = None,
        tr_s: float | None = None,
        flip_deg: float = 90.0,
        refocus_flip_deg: float = 180.0,
        bandwidth_hz_px: float = 200.0,
        partial_fourier: float = 1.0,
        excitation_duration_s: float = 2.5e-3,
        refocus_duration_s: float = 4e-3,
        refocus_thickness_factor: float = 1.25,
        crush_cycles_slice: float = 3.0,
        crush_cycles_readout: float = 0.0,
        spoil_cycles_per_voxel: float = 4.0,
        spoil_axis: str | Iterable[str] = ('x', 'z'),
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        fov_x, fov_y = require_pair(fov_mm, 'fov_mm')
        nx, ny = (int(matrix[0]), int(matrix[1]))
        self.fov_mm, self.matrix = (fov_x, fov_y), (nx, ny)
        self.thickness_mm = require_positive(thickness_mm, 'thickness_mm')
        self.echoes = self._require_echoes(echoes)
        raster_s = float(opts.grad_raster_time)
        block_raster_s = float(opts.block_duration_raster)

        self.exc = Excitation(opts=opts, flip_deg=flip_deg, thickness_mm=self.thickness_mm,
                              duration_s=excitation_duration_s)
        # prephase=False: in a CPMG train the readout's dephasing is half of a pair that straddles
        # a refocusing pulse, so there is no event of its own to cancel it.
        self.ro = CartesianLine(opts=opts, fov_mm=fov_x, matrix=nx, axis='x',
                                bandwidth_hz_px=bandwidth_hz_px,
                                partial_fourier=partial_fourier, prephase=False)
        # The readout-axis pair is a *balance* first and a crusher second: their mean is
        # `crush_cycles_readout` and their difference is the readout's own asymmetry about the
        # echo, which the readout knows and this module compensates.
        self.lobe_area_per_m = crush_cycles_readout / (self.voxel_mm('x') / 1e3)
        self.skew_per_m = self.ro.echo_moment_imbalance_per_m / 2.0

        def refocusing(window_s: float | None) -> Refocusing:
            return Refocusing(
                opts=opts, thickness_mm=self.thickness_mm * float(refocus_thickness_factor),
                flip_deg=refocus_flip_deg, duration_s=refocus_duration_s,
                crush_cycles_per_voxel=crush_cycles_slice, crush_voxel_mm=self.thickness_mm,
                crush_duration_s=window_s)

        def encode(window_s: float | None) -> PhaseEncode:
            return PhaseEncode(opts=opts, fov_mm=fov_y, matrix=ny, axis='y', duration_s=window_s)

        # One window, three axes.  Each axis states its own minimum; the widest wins.
        self.crush_min_s = {
            'z': refocusing(None).min_crush_duration_s,
            'x': max(float(pp.calc_duration(pp.make_trapezoid(
                'x', area=self.lobe_area_per_m + sign * self.skew_per_m, system=opts)))
                for sign in (1.0, -1.0)),
            'y': encode(None).min_duration_s,
        }
        crush_s = _shared_window_s(self.crush_min_s, raster_s)
        if echo_spacing_s is not None:
            probe = refocusing(crush_s)
            crush_s = _widened_window_s(
                crush_s,
                base_period_s=2 * (probe.time_to_crusher() - probe.time_to_center() + crush_s
                                   + self.ro.time_to_echo()),
                requested_period_s=require_positive(echo_spacing_s, 'echo_spacing_s'),
                raster_s=raster_s)

        self.crush_s = crush_s
        self.refoc, self.pe = refocusing(crush_s), encode(crush_s)
        self.lobe_in = pp.make_trapezoid('x', area=self.lobe_area_per_m + self.skew_per_m,
                                         duration=crush_s, system=opts)
        self.lobe_out = pp.make_trapezoid('x', area=self.lobe_area_per_m - self.skew_per_m,
                                          duration=crush_s, system=opts)
        self.gx_s = float(pp.calc_duration(self.ro.gx))
        self.readout_s = 2 * crush_s + self.gx_s
        self.ro_start_s = self.refoc.time_to_crusher()

        # Three minima, and `bound` reports which one applied.  The last is not a timing
        # preference: below it, z is on -- and the next RF may be on -- while the ADC is open,
        # which is signal loss rather than an illegal block.
        self.dephase_area_per_m = float(self.lobe_in.area) + self.ro.area_to_echo_per_m
        dephaser_min_s = float(pp.calc_duration(pp.make_trapezoid(
            'x', area=self.dephase_area_per_m, system=opts)))
        train_s = 2 * (self.ro_start_s - self.refoc.time_to_center() + crush_s
                       + self.ro.time_to_echo())
        self.bounds_s = {
            'train': train_s,
            'first interval': 2 * ((self.exc().duration - self.exc.time_to_center())
                                   + max(0.0, dephaser_min_s - crush_s)
                                   + self.refoc.time_to_center()),
            'adc clearance': 2 * (self.ro_start_s + self.readout_s - crush_s - train_s / 2),
        }
        self.bound = max(self.bounds_s, key=lambda name: self.bounds_s[name])
        half_s = ceil_raster(max(self.bounds_s.values()) / 2, raster_s)
        self._echo_spacing_s = 2 * half_s
        # Whatever the train does not need, the echo still has to sit at the midpoint between its
        # two pulses -- the primary and stimulated pathways only coincide there.
        self.shift_s = ceil_raster(half_s - train_s / 2, raster_s)
        self.t_refoc_s = ceil_raster(
            self.exc.time_to_center() + half_s - self.refoc.time_to_center(), raster_s)
        # The dephaser ends where the refocusing plateau starts, so it lives inside the *leading*
        # crusher window and only its overhang costs anything.
        dephaser_s = ceil_raster(
            max(dephaser_min_s, self.t_refoc_s + crush_s - self.exc().duration), raster_s)
        self.dephaser = pp.make_trapezoid('x', area=self.dephase_area_per_m,
                                          duration=dephaser_s, system=opts)

        self.spoilers = {
            axis: spoiler(opts, cycles_per_voxel=spoil_cycles_per_voxel,
                          voxel_mm=self.voxel_mm(axis), axis=axis)
            for axis in ((spoil_axis,) if isinstance(spoil_axis, str) else tuple(spoil_axis))
        }
        self._shot_s = (self.t_refoc_s + (self.echoes - 1) * self._echo_spacing_s
                        + self.ro_start_s + self.shift_s + self.readout_s
                        + max(block.duration for block in self.spoilers.values()))
        self._min_tr_s = ceil_raster(self._shot_s, block_raster_s)
        self._tr_s = self._resolve_tr(tr_s, block_raster_s)

    # ------------------------------------------------------------------ what it knows
    @property
    def center_line(self) -> int:
        """The phase-encode line that encodes ``k = 0``.  Read from the encode."""
        return self.pe.center_line

    @property
    def echo_spacing_s(self) -> float:
        """The spacing between refocusing centres this shot achieves, seconds."""
        return self._echo_spacing_s

    @property
    def shot_s(self) -> float:
        """Seconds the shot's own events occupy, spoiler included, before any TR fill."""
        return self._shot_s

    @property
    def min_tr_s(self) -> float:
        """The shortest legal :attr:`tr_s`: :attr:`shot_s` on the block raster."""
        return self._min_tr_s

    @property
    def tr_s(self) -> float:
        """The shot repetition time, seconds.  The block :meth:`build` returns is exactly this."""
        return self._tr_s

    def voxel_mm(self, axis: str) -> float:
        """The voxel dimension along `axis`: fov/matrix in plane, the slice along z."""
        return {'x': self.fov_mm[0] / self.matrix[0], 'y': self.fov_mm[1] / self.matrix[1],
                'z': self.thickness_mm}[axis]

    def time_to_echo(self, echo: int = 0) -> float:
        """
        Seconds from the start of the shot to ``k = 0`` of `echo`.

        Not measurable from the tree: a block knows when its events play, not which instant among
        them is the echo.
        """
        return (self.t_refoc_s + self._require_echo(echo) * self._echo_spacing_s
                + self.ro_start_s + self.shift_s + self.crush_s + self.ro.time_to_echo())

    def te_s(self, echo: int = 0) -> float:
        """Echo time of `echo`: the excitation's effective centre to that echo, seconds."""
        return self.time_to_echo(echo) - self.exc.time_to_center()

    # ------------------------------------------------------------------------- assembly
    def build(
        self,
        *,
        lines: Iterable[int],
        acquire: bool = True,
        segment_index: int | None = None,
        center_mm: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> LogicBlock:
        """
        Return one shot -- an excitation and its echo train -- as a block exactly `tr_s` long.

        Parameters
        ----------
        lines
            One phase-encode index per echo, in echo order.  ``len(lines)`` must be
            :attr:`echoes`, because the shot plays that many refocusing pulses whether or not
            they sample.  **Which** lines these are is the caller's choice; see
            :class:`~seqcraft.modules.FSE2D`.
        acquire
            ``False`` drops the ADC and the labels and changes nothing else, so a dummy shot
            loads the same gradients a real one does.
        segment_index
            The value for this shot's ``SEG`` label, or ``None`` to emit none.  **The value is the
            caller's and the placement is this module's**: a shot cannot know how many segments an
            acquisition has, and an acquisition should not have to compute where the first readout
            starts.  ``ECO`` is deliberately never emitted -- every echo of an FSE shot goes into
            *one* image at a different ``ky``, so an echo counter would tell a reconstruction to
            split what belongs together.
        center_mm
            ``(x, y, z)`` centre of the imaging volume, millimetres.  ``y`` must be ``0.0``.
        """
        table = self._check(lines)
        x, y, z = (float(v) for v in center_mm)
        if y != 0.0:
            msg = format_error(
                f'center_mm y = {y:g} mm is out of scope: a phase-encode shift is a per-line ADC '
                f'phase, not a per-shot frequency.',
                {'center_mm': center_mm},
                ['shift the slice with z, or the readout FOV with x'],
            )
            raise ConfigurationError(msg)

        # Excite about +y, refocus about +x -- 90 degrees apart is the fixed point of the map a
        # refocusing pulse applies to transverse phase, and the receiver stays at 0.
        out = LogicBlock('shot').add(0.0, self.exc(phase_deg=90.0, position_mm=z))
        out.add(self.t_refoc_s + self.crush_s - float(pp.calc_duration(self.dephaser)),
                self.dephaser)
        for n, line in enumerate(table):
            t_ref = self.t_refoc_s + n * self._echo_spacing_s
            t_ro = t_ref + self.ro_start_s + self.shift_s
            out.add(t_ref, self.refoc(phase_deg=0.0, position_mm=z))
            # The blip goes *after* the 180 and the rewinder *before the next one*, so ky is zero
            # at every refocusing centre.  Putting the blip before the pulse -- which a
            # single-echo reference legally does -- makes every later pulse flip the encode.
            out.add(t_ref + self.ro_start_s, self.pe(line=line))
            readout = LogicBlock('readout').add(0.0, self.lobe_in)
            readout.add(self.crush_s, self.ro(acquire=acquire, offset_mm=x))
            readout.add(self.crush_s + self.gx_s, self.lobe_out)
            out.add(t_ro, readout)
            # Not at the echo: that instant is an ADC-raster quantity, gradient starts are not,
            # and sc.compile refuses it by name.
            out.add(t_ro + self.readout_s - self.crush_s, self.pe(line=line, rewind=True))
            if acquire:
                out.add(t_ro, pp.make_label(type='SET', label='LIN', value=int(line)))
                if n == 0 and segment_index is not None:
                    out.add(t_ro, pp.make_label(type='SET', label='SEG',
                                                value=int(segment_index)))
        for block in self.spoilers.values():
            out.add(t_ro + self.readout_s, block)
        fill = ceil_raster(self._tr_s - self._shot_s, float(self.opts.block_duration_raster))
        if fill > 1e-9:
            out.add(self._shot_s, pp.make_delay(fill))
        return out

    # ---------------------------------------------------------------------- the refusals
    def _require_echoes(self, echoes: int) -> int:
        """Return `echoes` having checked there is at least one."""
        count = int(echoes)
        if count < 1:
            msg = format_error(
                f'echoes = {count} is not a train: a shot plays at least one refocusing pulse.',
                {'echoes': echoes}, ['echoes=1 is a conventional spin echo'])
            raise ConfigurationError(msg)
        return count

    def _require_echo(self, echo: int) -> int:
        """Return `echo` having checked it names one of this shot's echoes."""
        index = int(echo)
        if not 0 <= index < self.echoes:
            msg = format_error(
                f'echo {index} is outside a {self.echoes}-echo train.',
                {'echo': echo, 'echoes': self.echoes}, [])
            raise ConfigurationError(msg)
        return index

    def _resolve_tr(self, tr_s: float | None, block_raster_s: float) -> float:
        """Return the TR, refusing one the train does not fit inside."""
        if tr_s is not None and float(tr_s) < self._min_tr_s - 1e-12:
            msg = format_error(
                f'tr_s = {float(tr_s) * 1e3:.3f} ms is shorter than the '
                f'{self._min_tr_s * 1e3:.3f} ms a {self.echoes}-echo shot occupies, so the next '
                f'excitation would land inside the train.',
                {'tr_s': tr_s, 'min_tr_s': self._min_tr_s, 'echoes': self.echoes},
                [f'tr_s={self._min_tr_s:.6f} or longer', 'fewer echoes', 'a shorter readout'])
            raise ConfigurationError(msg)
        return ceil_raster(max(float(tr_s or 0.0), self._min_tr_s), block_raster_s)

    def _check(self, lines: Iterable[int]) -> tuple[int, ...]:
        """Return `lines` as a tuple, refusing a table this shot cannot play."""
        table = tuple(int(line) for line in lines)
        if len(table) != self.echoes:
            msg = format_error(
                f'{len(table)} line(s) for a {self.echoes}-echo shot: every shot plays '
                f'{self.echoes} refocusing pulses whether or not they sample.',
                {'lines': table[:8], 'echoes': self.echoes},
                ['pass one line per echo', 'or build the shot with a matching `echoes`'])
            raise ConfigurationError(msg)
        outside = sorted({line for line in table if not 0 <= line < self.matrix[1]})
        if outside:
            msg = format_error(
                f'{len(outside)} line(s) fall outside 0 ... {self.matrix[1] - 1}: {outside[:8]}.',
                {'lines': table[:8], 'matrix': self.matrix}, [])
            raise ConfigurationError(msg)
        return table
