"""
:class:`ContinuousEPI2D` -- a blipped echo-planar train sampled by one ADC window.

``readout/`` because it contains an ADC and no RF.  It is the sibling of
:class:`~seqcraft.modules.EPI2D`: the same alternating lobes, the same blips centred on the zero
crossings, the same line tables and the same interleaved-shot timing.  What differs is the
receiver.  One ADC event opens at the start of the first lobe and closes at the end of the last,
so every instant of the train is sampled -- flat tops, ramps, turnarounds and blips alike.

What one window changes
-----------------------
**The lobe is the sampling window.**  There is no guard: a lobe lasts a whole number of dwells,
``echo_spacing_s == samples_per_echo * dwell_s`` exactly, on the gradient raster.  The echo
spacing therefore moves in steps of ``lcm(dwell_s, grad_raster_time)``, and the samples per echo
are an output of the design rather than an input.

**Both polarities sample one kx grid.**  The lobe is a symmetric trapezoid and the window opens on
its first edge, so sample ``i`` of a forward lobe and sample ``S - 1 - i`` of the reverse lobe
after it sit at mirror times and have the same kx.

**``k = 0`` is an instant.**  The readout prephaser is exactly half a lobe's area, so each lobe
crosses ``kx = 0`` at its midpoint.  That is on a sample only when :attr:`samples_per_echo` is odd,
and halfway between two when it is even, so :meth:`time_to_echo` and :meth:`time_to_center_line`
report times and nothing here reports an echo *sample*.

**Every sample is placed by the trajectory.**  Ramp samples are not uniformly spaced in kx, and
samples taken while a blip plays lie between ky lines.  A reconstruction grids from
:meth:`k_per_m` at :meth:`sample_times_s`.  Each half-blip lies inside a readout ramp, so the
flat-top samples of every lobe sit on their line; the samples within ``blip_duration_s / 2`` of a
seam are the ones between lines.  Sample ``j`` belongs to lobe ``j // samples_per_echo``, which
acquired ``lines[j // samples_per_echo]`` with :meth:`polarity` of that lobe.

**The window's size is the scanner's limit.**  One train is one ADC event of ``len(lines) *
samples_per_echo`` samples.  It has to fit ``opts.adc_samples_limit``, the per-event limit the
interpreter accepts, and be a multiple of ``opts.adc_samples_divisor``.  A table above the limit
is refused; the window is never split into several events, and an ``Opts`` with no limit set is
refused at construction.

Out of scope: navigator echoes, readout partial echo, a readout field-of-view offset, ``REF`` and
``IMA`` labels, segmentation of the window inside the interpreter, and reconstruction.
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ...design.events import derive, knots_of
from ...design.logic import LogicBlock
from ...design.module import Module
from ...design.timing import EPS, from_ticks, to_ticks
from ...errors import ConfigurationError, format_error
from .._support import (
    area_until,
    ceil_raster,
    check_epi_table,
    design_blip,
    require_axis,
    require_count,
    require_pair,
    require_positive,
    shot_spacings,
)
from ..encoding.phase_encoding import PhaseEncode

if TYPE_CHECKING:
    from collections.abc import Iterable

    from pypulseq.opts import Opts

    from ...design.events import Event

__all__ = ['ContinuousEPI2D']

#: How many echo spacings, one lattice step apart, the lobe search tries before refusing.
_MAX_SPACING_STEPS = 4096

#: How many ADC rasters either side of `dwell_s` the divisor refusal searches for a dwell whose
#: samples per echo the divisor divides.
_DWELL_SEARCH_STEPS = 64


class ContinuousEPI2D(Module):
    """
    A blipped 2D echo-planar readout sampled by one ADC window across the whole train.

    Parameters
    ----------
    opts
        The scanner.  ``opts.adc_samples_limit`` must be positive: a train's one ADC event has to
        fit it, and ``0``, pypulseq's "no limit", is refused.
    fov_mm
        ``(readout, blip)`` in millimetres, or a scalar meaning square.  Sets both k-space steps.
    matrix
        ``(nx, ny)`` -- readout then blip axis.  Each lobe traverses ``nx * dk_read_per_m``.
    dwell_s
        The receiver dwell, a whole number of ``opts.adc_raster_time``.  Never rounded.
    oversampling
        A bound on the readout sampling density, at least ``1``: every kx step between
        consecutive samples is at most ``dk_read_per_m / oversampling``.  The largest step is
        on the flat top, so the lobe is designed to satisfy::

            readout_amplitude_hz_m * dwell_s <= dk_read_per_m / oversampling

        ``1`` samples the flat top at the field of view, ``2`` at least twice as densely.  It is
        not :class:`~seqcraft.modules.EPI2D`'s sample-count multiplier: the samples per echo
        are an output here.  It bounds the spacing along kx only, and says nothing about whether
        the 2D sampling -- the ky lines, and the samples taken during the blips -- is adequate
        for a given reconstruction.
    blip_lines
        The largest ``|delta ky|``, in lines, that any table this instance builds steps between
        consecutive echoes.  It sets the blip, which every step scales and none lengthens.
    shots
        Interleaved shots.  Above ``1``, :meth:`build` delays shot ``s`` by
        ``s * echo_spacing_s / shots`` and pads its tail by the complement, so the acquisition
        time is affine in ``ky`` across shots.  The echo spacing is chosen so that delay is a
        whole number of block rasters.
    blip_duration_s
        Lengthen the blip beyond its own minimum; ``None`` is the minimum.  Each readout ramp
        holds half a blip, so a longer blip can lengthen the ramps and the echo spacing with them.
    prephase
        ``False`` omits both prephasers: the spin-echo placement, where the caller dephases
        before a refocusing pulse.  The block then opens with the receiver's dead time.  What the
        reported trajectory means in that case is under :meth:`k_per_m`.
    axis, blip_axis
        Logical gradient channels for the readout and the blips.  Must differ.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    gx : SimpleNamespace
        The positive readout lobe; odd lobes play its negation.
    pe : PhaseEncode
        The blip-axis prephaser, which also owns the line convention and :attr:`center_line`.
    samples_per_echo : int
        Samples per lobe.  ``echo_spacing_s == samples_per_echo * dwell_s`` exactly.
    echo_spacing_s : float
        The lobe duration, seconds.
    dk_read_per_m, dk_blip_per_m : float
        The k-space steps along the readout and the blip axis, 1/m.
    shot_delay_s : float
        ``echo_spacing_s / shots``, and ``0.0`` at one shot.

    Examples
    --------
    >>> import seqcraft as sc
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s',
    ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
    ...          adc_samples_limit=8192)
    >>> epi = ContinuousEPI2D(opts=o, fov_mm=220.0, matrix=(64, 64), dwell_s=5e-6,
    ...                       oversampling=1.5)
    >>> epi.samples_per_echo, round(epi.echo_spacing_s * 1e6)
    (118, 590)
    >>> epi.samples_per_echo * epi.dwell_s == epi.echo_spacing_s
    True

    The flat-top step is the bound, here exactly:

    >>> round(epi.readout_amplitude_hz_m * epi.dwell_s / epi.dk_read_per_m, 4)
    0.6667

    One train, one ADC event, every lobe inside it:

    >>> shot = epi(lines=range(64))
    >>> [int(e.num_samples) for _, e, _ in sc.flatten(shot) if e.type == 'adc']
    [7552]

    The echo time is measured to the midpoint of the lobe that acquires the centre line:

    >>> round(epi.time_to_center_line(range(64)) * 1e3, 3)
    19.525
    >>> epi.time_to_center_line(range(64)) == epi.time_to_echo(epi.center_line)
    True

    Notes
    -----
    **The lobe design.**  With ``A = nx * dk_read_per_m``, the echo spacing is the shortest ``T``
    on ``lcm(dwell_s, grad_raster_time)`` -- and, with interleaved shots, divisible into
    `shots` block rasters -- for which a trapezoid of duration ``T`` and ramp ``r`` on the
    gradient raster, with amplitude ``G = A / (T - r)``, meets all of::

        G <= max_grad      G / r <= max_slew      G * dwell_s <= dk_read_per_m / oversampling
        r >= blip_duration_s / 2

    The ramp is the shortest that meets them, which gives the lowest amplitude at that ``T``.
    """

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: float | tuple[float, float],
        matrix: tuple[int, int],
        dwell_s: float,
        oversampling: float,
        blip_lines: int = 1,
        shots: int = 1,
        blip_duration_s: float | None = None,
        prephase: bool = True,
        axis: str = 'x',
        blip_axis: str = 'y',
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self._limit = self._require_limit()
        self._divisor = int(getattr(opts, 'adc_samples_divisor', 1) or 1)
        fov_read, fov_blip = require_pair(fov_mm, 'fov_mm')
        nx, ny = require_pair(matrix, 'matrix')
        self.fov_mm = (require_positive(fov_read, 'fov_mm'), require_positive(fov_blip, 'fov_mm'))
        self.matrix = (require_count(int(nx), 'matrix', low=2), int(ny))
        self.axis = require_axis(axis)
        self.blip_axis = require_axis(blip_axis, 'blip_axis')
        self.prephase = bool(prephase)
        self.dwell_s = self._require_dwell(dwell_s)
        self.oversampling = self._require_oversampling(oversampling)
        self.blip_lines = require_count(
            blip_lines, 'blip_lines', hint='1 is single-shot and blocked multi-shot',
        )
        self.shots = require_count(
            shots, 'shots', hint='1 is single-shot, blocked, or accelerated single-shot',
        )
        self._check_axes()

        self.dk_read_per_m = 1e3 / self.fov_mm[0]
        self.dk_blip_per_m = 1e3 / self.fov_mm[1]
        self._blip = design_blip(
            opts=opts, axis=self.blip_axis, area_per_m=self.blip_lines * self.dk_blip_per_m,
            blip_lines=self.blip_lines, requested_s=blip_duration_s,
        )

        spacing, ramp, amplitude = self._solve(to_ticks(self.dwell_s), self.oversampling)
        self.gx = pp.make_trapezoid(
            channel=self.axis, amplitude=amplitude, rise_time=from_ticks(ramp),
            flat_time=from_ticks(spacing - 2 * ramp), fall_time=from_ticks(ramp), system=opts,
        )
        self._gx_reverse = derive(pp.scale_grad(self.gx, -1.0))
        self.echo_spacing_s = from_ticks(spacing)
        self.samples_per_echo = spacing // to_ticks(self.dwell_s)
        self.shot_delay_s = self.echo_spacing_s / self.shots if self.shots > 1 else 0.0

        # The ADC node sits this long before the first lobe, with the same lead as its own delay,
        # so the first sample is half a dwell into the first lobe.  Rounded onto the gradient
        # raster because the lobes start on it.
        self._lead_s = ceil_raster(float(opts.adc_dead_time), opts.grad_raster_time)
        self.pe, self._prephaser, self._prephase_s = self._design_prephasers()

    # ------------------------------------------------------------------ what it knows
    @property
    def center_line(self) -> int:
        """The blip-axis line that encodes ``ky = 0``: ``ny // 2``, as `PhaseEncode` defines it."""
        return self.pe.center_line

    @property
    def blip_duration_s(self) -> float:
        """Seconds one blip occupies, whatever its step.  Every step takes the same time."""
        return float(pp.calc_duration(self._blip))

    @property
    def blip_area_per_m(self) -> float:
        """The designed blip's area, 1/m: ``blip_lines * dk_blip_per_m``, as achieved."""
        return float(self._blip.area)

    @property
    def readout_amplitude_hz_m(self) -> float:
        """Peak amplitude of a readout lobe, Hz/m.  Times `dwell_s`, the largest kx step."""
        return float(self.gx.amplitude)

    @property
    def area_to_echo_per_m(self) -> float:
        """
        What the first lobe accumulates by its midpoint, 1/m: half the lobe's area.

        It is what the readout prephaser cancels, and what a caller dephasing before a
        refocusing pulse supplies with the opposite sign.  Available with ``prephase=False``.
        """
        return area_until(self.gx, self.echo_spacing_s / 2)

    def k_blip_per_m(self, line: int) -> float:
        """Return the blip-axis k the train must start at to acquire `line` first, in 1/m."""
        return self.pe.k_per_m(line)

    def polarity(self, echo: int) -> int:
        """Return ``+1`` for an even lobe, which sweeps kx upwards, and ``-1`` for an odd one."""
        return 1 if self._check_echo(echo) % 2 == 0 else -1

    def shot_offset_s(self, segment: int | None) -> float:
        """
        Seconds shot `segment` waits before its train: ``segment * echo_spacing_s / shots``.

        Zero at ``shots=1`` and for ``None``, which is what a dummy replaying shot zero wants.
        As in :meth:`EPI2D.shot_offset_s <seqcraft.modules.EPI2D.shot_offset_s>`, the delay makes
        the acquisition time affine in ``ky`` across an interleaved plan.
        """
        if segment is None or self.shots == 1:
            return 0.0
        index = int(segment)
        if not 0 <= index < self.shots:
            msg = format_error(
                f'segment = {index} is not one of the {self.shots} interleaved shots this '
                f'readout was designed for.',
                {'segment': index, 'shots': self.shots, 'shot_delay_s': self.shot_delay_s},
                [
                    f'segments are zero-based, so this instance takes 0 ... {self.shots - 1}',
                    'segment sets the train delay as well as the SEG label once shots > 1, so a '
                    'segment outside the plan would place a shot at a time no ky line wants',
                    'pass shots=1 if the index is a label rather than an interleave position -- '
                    'a blocked table numbers its bands without interleaving them',
                ],
            )
            raise ConfigurationError(msg)
        return index * self.shot_delay_s

    def time_to_echo(self, echo: int = 0, *, segment: int | None = None) -> float:
        """
        Seconds from the start of this module's block to lobe `echo`'s ``kx = 0``: its midpoint.

        An instant, which falls on a sample only when :attr:`samples_per_echo` is odd.  Like
        :meth:`k_per_m`, it assumes the train starts at ``kx = -area_to_echo_per_m``.
        """
        index = self._check_echo(echo)
        return (self.shot_offset_s(segment) + self._train_start_s
                + (index + 0.5) * self.echo_spacing_s)

    def time_to_center_line(self, lines: Iterable[int], *, segment: int | None = None) -> float:
        """
        Seconds from the start of this block to the origin crossing: the echo time's endpoint.

        The midpoint of the lobe that acquires :attr:`center_line`, where ``kx = 0`` and
        ``ky = 0`` together, because the blips play only on the ramps.

        Parameters
        ----------
        lines, segment
            The same table and segment :meth:`build` takes; the answer moves with both.

        Raises
        ------
        ConfigurationError
            On any table :meth:`build` refuses, and on one that omits :attr:`center_line`.
        """
        table = self._check_train(lines)
        if self.center_line not in table:
            msg = format_error(
                f'lines omits the centre of k-space (line {self.center_line}), so there is no '
                f'k = 0 to measure to.',
                {'center_line': self.center_line, 'lines': table[:8]},
                [
                    'include the centre line in this shot',
                    'under segmentation only the shot that carries it has a meaningful echo '
                    'time, and this is the query that says which one that is',
                ],
            )
            raise ConfigurationError(msg)
        return self.time_to_echo(table.index(self.center_line), segment=segment)

    def sample_times_s(self, lines: Iterable[int], *, segment: int | None = None) -> np.ndarray:
        """
        Return the time of every ADC sample of a train, seconds from the start of its block.

        ``len(lines) * samples_per_echo`` of them, read off the ADC event :meth:`build` emits:
        its node, its own ``delay``, and the centre of each dwell.
        """
        table = self._check_train(lines)
        adc = self._adc_for(len(table), 0.0)
        start_s = self.shot_offset_s(segment) + self._prephase_s + float(adc.delay)
        return start_s + (np.arange(int(adc.num_samples)) + 0.5) * float(adc.dwell)

    def k_per_m(self, lines: Iterable[int]) -> np.ndarray:
        """
        Return ``(2, len(lines) * samples_per_echo)``: readout and blip-axis k at every sample, 1/m.

        Integrated exactly from the lobes and blips :meth:`build` emits for this table.  It is a
        **nominal** trajectory: it starts the train at ``kx = -area_to_echo_per_m`` and
        ``ky = k_blip_per_m(lines[0])``, which is what this module's own prephasers deliver.
        With ``prephase=False`` it is what the composition must deliver for it to hold, and it
        is not a measurement of any particular composition: that has to be checked on the
        assembled sequence.
        """
        table = self._check_train(lines)
        x_times, x_amps, y_times, y_amps = self._train_knots(table)
        at = (np.arange(len(table) * self.samples_per_echo) + 0.5) * self.dwell_s
        kx = -self.area_to_echo_per_m + _integral_at(x_times, x_amps, at)
        ky = self.k_blip_per_m(table[0]) + _integral_at(y_times, y_amps, at)
        return np.vstack([kx, ky])

    # ----------------------------------------------------------------------- assembly
    def build(
        self,
        *,
        lines: Iterable[int],
        acquire: bool = True,
        segment: int | None = None,
        phase_deg: float = 0.0,
    ) -> LogicBlock:
        """
        Return one shot: the prephasers, then one lobe per entry of `lines`, sampled by one ADC.

        Parameters
        ----------
        lines
            Blip-axis indices, **in acquisition order**, for this shot.  Consecutive entries may
            step by at most `blip_lines`, and ``len(lines) * samples_per_echo`` has to fit
            ``opts.adc_samples_limit`` and be a multiple of ``opts.adc_samples_divisor``.
        acquire
            ``False`` plays identical gradients with no ADC and no label, over the same span,
            for a dummy shot.
        segment
            The shot's index.  Emitted as a ``SEG`` label on the ADC, and above ``shots=1`` it
            also sets :meth:`shot_offset_s`.
        phase_deg
            The receiver phase of the train's ADC: the carrier phase its excitation was given.

        Notes
        -----
        The block holds the shot offset; the two prephasers, or with ``prephase=False`` nothing;
        the receiver's leading dead time; ``len(lines)`` lobes alternating from positive, with
        a blip centred on each seam; the receiver's trailing dead time; and the interleave tail,
        so every shot of a plan spans the same time.  There is no per-echo label: one ADC event
        cannot carry ``LIN`` or ``REV`` per lobe.
        """
        table = self._check_train(lines)
        offset_s = self.shot_offset_s(segment)
        out = LogicBlock()
        if self.prephase:
            out.add(offset_s, self._prephaser)
            out.add(offset_s, self.pe(line=table[0]))

        train_s = offset_s + self._train_start_s
        for n, line in enumerate(table):
            t0 = train_s + n * self.echo_spacing_s
            out.add(t0, self.gx if n % 2 == 0 else self._gx_reverse)
            if n + 1 < len(table):
                seam_s = t0 + self.echo_spacing_s
                out.add(seam_s - self.blip_duration_s / 2, self._blip_for(table[n + 1] - line))

        node_s = offset_s + self._prephase_s
        adc = self._adc_for(len(table), float(phase_deg))
        if acquire:
            out.add(node_s, adc)
            if segment is not None:
                out.add(node_s, pp.make_label(type='SET', label='SEG', value=int(segment)))
        else:
            # The receiver's trailing dead time, which the ADC event would have reserved.
            end_s = train_s + len(table) * self.echo_spacing_s
            out.add(end_s, pp.make_delay(float(adc.dead_time)))

        tail_s = (self.shots - 1) * self.shot_delay_s - offset_s
        if tail_s > EPS:
            out.add(out.duration, pp.make_delay(tail_s))
        return out

    def _adc_for(self, lines: int, phase_deg: float) -> Event:
        """Return the train's one ADC: every lobe of a `lines`-lobe train, at `phase_deg`."""
        return pp.make_adc(
            num_samples=lines * self.samples_per_echo, dwell=self.dwell_s, delay=self._lead_s,
            phase_offset=float(np.deg2rad(phase_deg)), system=self.opts,
        )

    def _blip_for(self, step: int) -> Event:
        """Return the blip that steps `step` lines along the blip axis, at the fixed duration."""
        scale = step * self.dk_blip_per_m / float(self._blip.area)
        return derive(pp.scale_grad(self._blip, scale))

    def _train_knots(
        self, table: tuple[int, ...],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Return the readout and blip knots of `table`'s train, from its first lobe's start."""
        x_times, x_amps, y_times, y_amps = [], [], [], []
        for n, line in enumerate(table):
            t0 = n * self.echo_spacing_s
            times, amps = knots_of(self.gx if n % 2 == 0 else self._gx_reverse, t0)
            x_times.append(times)
            x_amps.append(amps)
            if n + 1 < len(table):
                start = t0 + self.echo_spacing_s - self.blip_duration_s / 2
                times, amps = knots_of(self._blip_for(table[n + 1] - line), start)
                y_times.append(times)
                y_amps.append(amps)
        empty = np.zeros(0)
        return (np.concatenate(x_times), np.concatenate(x_amps),
                np.concatenate(y_times) if y_times else empty,
                np.concatenate(y_amps) if y_amps else empty)

    @property
    def _train_start_s(self) -> float:
        """When the first lobe starts, from the start of an unshifted block."""
        return self._prephase_s + self._lead_s

    # -------------------------------------------------------------------------- design
    def _solve(self, dwell_ticks: int, oversampling: float) -> tuple[int, int, float]:
        """
        Return ``(spacing, ramp, amplitude)`` of the shortest lobe, the two times in ticks.

        A pure function of the dwell and the bound, so a refusal can price an alternative
        without building it.  Feasibility only improves as the spacing grows -- every limit is
        on ``A / (T - r)`` at a fixed ramp -- so the first spacing that works is the answer.
        """
        raster = to_ticks(self.opts.grad_raster_time)
        lattice = math.lcm(dwell_ticks, raster)
        area = self.matrix[0] * self.dk_read_per_m
        block_s = float(self.opts.block_duration_raster)
        shortest_ramp = max(raster, to_ticks(self.blip_duration_s) // 2)
        # The sampling bound, G * dwell <= dk / oversampling, is T - r >= nx * oversampling * dwell.
        sampled = Fraction(self.matrix[0]) * Fraction(oversampling) * dwell_ticks
        first = max(1, math.ceil((sampled + shortest_ramp) / lattice))
        for k in range(first, first + _MAX_SPACING_STEPS):
            spacing = k * lattice
            if shot_spacings(from_ticks(spacing), self.shots, block_s) is not None:
                continue
            ramp = self._shortest_ramp(spacing, shortest_ramp, raster, area)
            if ramp is None or spacing - ramp < sampled:
                continue
            amplitude = area / from_ticks(spacing - ramp)
            if amplitude <= float(self.opts.max_grad):
                return spacing, ramp, amplitude
        longest = (first + _MAX_SPACING_STEPS - 1) * lattice
        msg = format_error(
            f'no echo spacing up to {from_ticks(longest) * 1e3:.1f} ms lets a lobe traverse '
            f'{area:.1f} 1/m within the gradient limits and the oversampling bound.',
            {'dwell_s': from_ticks(dwell_ticks), 'oversampling': oversampling,
             'matrix': self.matrix, 'fov_mm': self.fov_mm, 'shots': self.shots},
            ['widen fov_mm or lower matrix[0], which lowers the area a lobe traverses'],
        )
        raise ConfigurationError(msg)

    def _shortest_ramp(self, spacing: int, shortest: int, raster: int, area: float) -> int | None:
        """
        Return the shortest ramp, in ticks, at which a `spacing` lobe meets the slew limit.

        ``G / r <= max_slew`` with ``G = A / (T - r)`` is ``r * (T - r) >= A / max_slew``, which
        grows with ``r`` up to ``T / 2``.  The closed form seeds the search, and the raster check
        decides.  ``None`` when no ramp up to half the lobe is enough.
        """
        slew = float(self.opts.max_slew)
        total = from_ticks(spacing)

        def fits(ramp: int) -> bool:
            return area / from_ticks(spacing - ramp) <= slew * from_ticks(ramp) * (1 + 1e-9)

        discriminant = total * total - 4 * area / slew
        if discriminant < 0.0:
            return None
        seed = to_ticks((total - math.sqrt(discriminant)) / 2)
        ramp = max(shortest, -(-seed // raster) * raster)
        while ramp - raster >= shortest and fits(ramp - raster):
            ramp -= raster
        while 2 * ramp <= spacing and not fits(ramp):
            ramp += raster
        return ramp if 2 * ramp <= spacing else None

    def _design_prephasers(self) -> tuple[PhaseEncode, Event | None, float]:
        """
        Return the blip-axis prephaser, the readout-axis one, and the window both share.

        Each axis reports its own minimum and the longer wins, so the two play at once.  With
        ``prephase=False`` there is no window, and the blip-axis module is kept only for its
        line convention.
        """
        probe = PhaseEncode(opts=self.opts, fov_mm=self.fov_mm[1], matrix=self.matrix[1],
                            axis=self.blip_axis)
        if not self.prephase:
            return probe, None, 0.0
        area = -self.area_to_echo_per_m
        shortest = pp.make_trapezoid(channel=self.axis, area=area, system=self.opts)
        window_s = ceil_raster(
            max(float(pp.calc_duration(shortest)), probe.min_duration_s),
            self.opts.grad_raster_time,
        )
        prephaser = pp.make_trapezoid(
            channel=self.axis, area=area, duration=window_s, system=self.opts,
        )
        blip_axis_prephaser = PhaseEncode(
            opts=self.opts, fov_mm=self.fov_mm[1], matrix=self.matrix[1], axis=self.blip_axis,
            duration_s=window_s,
        )
        return blip_axis_prephaser, prephaser, window_s

    # -------------------------------------------------------------------- the refusals
    def _check_train(self, lines: Iterable[int]) -> tuple[int, ...]:
        """Return `lines` as a tuple, refusing a table or an ADC event this train cannot play."""
        table = check_epi_table(
            lines, matrix=self.matrix, center_line=self.center_line,
            blip_lines=self.blip_lines, blip_duration_s=self.blip_duration_s,
            step_consequence='it lengthens the blip, and the readout ramps that must each hold '
                             'half of it',
        )
        count = len(table) * self.samples_per_echo
        if count > self._limit:
            self._refuse_limit(len(table), count)
        if count % self._divisor:
            self._refuse_divisor(len(table), count)
        return table

    def _refuse_limit(self, lines: int, count: int) -> None:
        """Refuse a train whose one ADC event is above ``opts.adc_samples_limit``."""
        per_echo, limit = self.samples_per_echo, self._limit
        fits = limit // per_echo
        longest = next((n for n in range(fits, 0, -1) if (n * per_echo) % self._divisor == 0), 0)
        if longest == fits:
            first = f'a table of at most {fits} lines fits one event at this design'
        elif longest:
            first = (f'a table of at most {longest} lines builds at this design: {fits} fit one '
                     f'event, and {longest} is the longest whose count '
                     f'opts.adc_samples_divisor = {self._divisor} divides')
        else:
            first = f'no table of lines of {per_echo} samples builds within this limit'
        hints = [
            first,
            'fewer lines per excitation changes the acquisition: more interleaved shots (shots '
            'with blip_lines) acquire the rest in other excitations, and a partial-ky table '
            'leaves lines unacquired for a partial-Fourier reconstruction to fill',
        ]
        if self.oversampling > 1:
            lowest = self._solve(to_ticks(self.dwell_s), 1.0)[0] // to_ticks(self.dwell_s)
            if lowest < per_echo:
                hints.append(
                    f'oversampling=1 gives {lowest} samples per echo, so {limit // lowest} '
                    f'lines fit -- the flat-top sample spacing then reaches dk_read_per_m',
                )
        hints += [
            'EPI2D samples the same train with one ADC event per echo',
            'this readout never splits its window into several events, which would put '
            'receiver dead time inside the train',
        ]
        msg = format_error(
            f'{count} ADC samples in this train\'s one ADC event ({lines} lines x {per_echo} '
            f'per echo), above the {limit} that opts.adc_samples_limit allows.',
            {'lines': lines, 'samples_per_echo': per_echo, 'adc_samples_limit': limit,
             'dwell_s': self.dwell_s, 'oversampling': self.oversampling},
            hints,
        )
        raise ConfigurationError(msg)

    def _refuse_divisor(self, lines: int, count: int) -> None:
        """Refuse a train whose one ADC event ``opts.adc_samples_divisor`` does not divide."""
        per_echo, divisor = self.samples_per_echo, self._divisor

        def builds(n: int, samples: int) -> bool:
            return (n * samples) % divisor == 0 and n * samples <= self._limit

        below = next((n for n in range(lines - 1, 0, -1) if builds(n, per_echo)), None)
        above = next((n for n in range(lines + 1, self.matrix[1] + 1) if builds(n, per_echo)),
                     None)
        lengths = [str(n) for n in (below, above) if n is not None]
        hints = []
        if lengths:
            hints.append(
                f'a table of {" or ".join(lengths)} lines divides -- the table is the ky '
                f'coverage, so a different length acquires a different set of lines',
            )
        found = self._dwell_with_divisible_samples(lines)
        if found is not None:
            dwell_s, samples, spacing_s = found
            hints.append(
                f'or dwell_s={dwell_s:.6g}: {samples} samples per echo, a multiple of '
                f'{divisor}, so every table length divides -- the echo spacing becomes '
                f'{spacing_s * 1e6:.1f} us against {self.echo_spacing_s * 1e6:.1f} us here, '
                f'and the receiver bandwidth changes with the dwell',
            )
        msg = format_error(
            f'{count} ADC samples ({lines} lines x {per_echo} per echo) is not a multiple of '
            f'opts.adc_samples_divisor = {divisor}, which PyPulseq\'s timing check requires of '
            f'every ADC event.',
            {'lines': lines, 'samples_per_echo': per_echo, 'adc_samples_divisor': divisor},
            hints,
        )
        raise ConfigurationError(msg)

    def _dwell_with_divisible_samples(self, lines: int) -> tuple[float, int, float] | None:
        """
        Return the nearest dwell whose samples per echo the divisor divides, or ``None``.

        As ``(dwell_s, samples_per_echo, echo_spacing_s)``, and only one at which a train of
        `lines` lines also fits the per-event limit.
        """
        raster = to_ticks(self.opts.adc_raster_time)
        here = to_ticks(self.dwell_s)
        for k in range(1, _DWELL_SEARCH_STEPS + 1):
            for dwell in (here + k * raster, here - k * raster):
                if dwell <= 0:
                    continue
                spacing = self._solve(dwell, self.oversampling)[0]
                samples = spacing // dwell
                if samples % self._divisor == 0 and lines * samples <= self._limit:
                    return from_ticks(dwell), samples, from_ticks(spacing)
        return None

    def _check_echo(self, echo: int) -> int:
        """Return `echo` as a lobe index, refusing one before the train."""
        index = int(echo)
        if index < 0:
            msg = format_error(
                f'echo = {index} is before the start of the train.',
                {'echo': index},
                ['echoes are counted from zero, and this train builds no navigator echoes'],
            )
            raise ConfigurationError(msg)
        return index

    def _check_axes(self) -> None:
        """Refuse a readout and a blip on one channel, which would sum rather than encode."""
        if self.axis == self.blip_axis:
            msg = format_error(
                f'axis and blip_axis are both {self.axis!r}, so the blips would be summed into '
                f'the readout rather than encoding anything.',
                {'axis': self.axis, 'blip_axis': self.blip_axis},
                ["the usual pair is axis='x', blip_axis='y'"],
            )
            raise ConfigurationError(msg)

    def _require_limit(self) -> int:
        """Return ``opts.adc_samples_limit``, refusing an ``Opts`` that sets none."""
        limit = int(getattr(self.opts, 'adc_samples_limit', 0) or 0)
        if limit <= 0:
            msg = format_error(
                f'opts.adc_samples_limit is {limit}, so there is no per-event sample limit for '
                f'this readout\'s one ADC event to fit.',
                {'adc_samples_limit': limit},
                [
                    'set it to the per-event limit of the interpreter the sequence is for, e.g. '
                    'pp.Opts(adc_samples_limit=8192) or sc.opts.from_scanner(..., '
                    'adc_samples_limit=8192); the examples here use 8192',
                    'this readout never splits its window into several events, so that limit is '
                    'what bounds a train',
                ],
            )
            raise ConfigurationError(msg)
        return limit

    def _require_dwell(self, dwell_s: float) -> float:
        """Return `dwell_s`, refusing one that is not a whole number of ADC rasters."""
        dwell = require_positive(dwell_s, 'dwell_s')
        raster = to_ticks(self.opts.adc_raster_time)
        ticks = to_ticks(dwell)
        if ticks % raster:
            legal = [from_ticks(n) for n in ((ticks // raster) * raster,
                                             -(-ticks // raster) * raster) if n > 0]
            msg = format_error(
                f'dwell_s = {dwell * 1e6:.4g} us is not a whole number of adc_raster_time '
                f'({from_ticks(raster) * 1e6:g} us).',
                {'dwell_s': dwell, 'adc_raster_time': from_ticks(raster)},
                [f'the nearest legal dwells are {" and ".join(f"{d:.6g}" for d in legal)}'],
            )
            raise ConfigurationError(msg)
        return from_ticks(ticks)

    def _require_oversampling(self, oversampling: float) -> float:
        """Return `oversampling` as a float, refusing a bound below the field of view."""
        value = float(oversampling)
        if not np.isfinite(value) or value < 1.0:
            msg = format_error(
                f'oversampling = {value:g} would step more than dk_read_per_m between samples on '
                f'the flat top, which undersamples the field of view.',
                {'oversampling': value},
                [
                    'oversampling=1 samples the flat top at the field of view',
                    'oversampling=2 samples it at least twice as densely',
                ],
            )
            raise ConfigurationError(msg)
        return value


def _integral_at(times: np.ndarray, amps: np.ndarray, at: np.ndarray) -> np.ndarray:
    """
    Return the exact running integral of piecewise-linear knots at each time in `at`.

    Zero before the first knot and the total after the last.  Knots may repeat a time, which is
    how two gradients that meet at zero are joined.
    """
    if times.size < 2:
        return np.zeros_like(at)
    span = np.diff(times)
    total = np.concatenate([[0.0], np.cumsum(0.5 * (amps[1:] + amps[:-1]) * span)])
    index = np.clip(np.searchsorted(times, at, side='right') - 1, 0, times.size - 2)
    width = span[index]
    into = np.clip(at - times[index], 0.0, width)
    slope = np.divide(amps[index + 1] - amps[index], width,
                      out=np.zeros_like(into), where=width > 0)
    return total[index] + into * (amps[index] + 0.5 * slope * into)
