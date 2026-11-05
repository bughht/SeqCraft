"""
:class:`WaveReadout` -- a Cartesian readout with Wave gradients on the other two axes.

One Cartesian line, as :class:`~seqcraft.modules.CartesianLine` reads it, plus a sinusoid on one
phase-encoding axis and a cosinusoid on the other while the ADC samples.  The sampled trajectory
is a corkscrew about the line, which spreads aliasing along the readout and is what Wave-CAIPI
reconstructs.  This module owns the readout, the ADC timing, the Wave gradients, their realisation
under `opts`, and the moment the corkscrew's centring requires before the readout.  Which line and
partition are encoded, the acceleration and CAIPI pattern, flow compensation, reconstruction and
PSF calibration all belong to the caller.

The timing convention
---------------------
The Wave is anchored to the **ADC dwell window**, ``[t_a, t_a + N * dwell]``, where ``t_a`` is the
ADC start (the readout lobe's start plus ``adc.delay``) and ``N`` is the raw, oversampled sample
count::

    omega = 2 * pi * cycles / (N * dwell)
    G_sine(t)   = G0 * sin(omega * (t - t_a))        on the window, 0 outside
    G_cosine(t) = G0 * cos(omega * (t - t_a))        on the window, with ramps outside

Sample ``n`` is at ``t_a + (n + 1/2) * dwell`` (pypulseq's convention), and the echo is sample
``N // 2``.  `cycles` is therefore an integer number of periods over the acquisition, and the
sample phases repeat exactly every ``N / cycles`` samples.

The ADC and the gradients live on **independent rasters**.  The dwell is on the ADC raster, the ADC
delay on the RF raster, and the Wave's knots on gradient-raster edges.  Nothing requires the
window, a dwell or a Wave period to be a whole number of gradient rasters, or a sample to fall on
a knot.

How the Wave is realised
------------------------
Each Wave axis is one extended trapezoid whose knots are the continuous prescription evaluated at
every gradient-raster edge covering it: the sinusoid, plus on the cosine axis a linear entrance
ramp ``0 -> G0`` ending at ``t_a`` and an exit ramp ``G0 -> 0`` starting at the window's end.
Each ramp is the shortest one on the raster at ``max_slew``, and the exit starts as soon as the
window closes.  The emitted gradient is the piecewise-linear interpolant of the prescription, so
its peak is at most ``G0`` and its slope at most the prescription's.  An amplitude whose
prescription is within `opts` is therefore emitted within `opts`, wherever the window falls on the
raster.

Every trajectory number below is integrated from the emitted knots and evaluated at the actual
sample times, so it describes what is played, not the analytic prescription.

Centring, and what ``prephase`` means
-------------------------------------
``centering`` sets where the corkscrew sits relative to the encoded line, as a condition on the
ADC-sampled trajectory:

* ``'sinzero'`` -- k at the echo sample equals the encoding target on both Wave axes.
* ``'sinctr'`` -- on the sine axis, the mean of k over the samples equals the encoding target, so
  the corkscrew is centred on the line.  The cosine axis keeps its echo sample on target.

Both are met by a **moment in place before the readout**, never by changing the active Wave, which
is identical under either.  That moment, per axis, is always stated
(:meth:`prephaser_area_per_m`).  ``prephase`` says only who plays it:

* ``True`` -- this module.  The readout axis gets ``CartesianLine``'s prephaser, the sine axis one
  lobe, and the cosine axis one event from rest to ``G0`` that is both its centring moment and its
  entrance ramp.  The block then puts k on the contract by itself.
* ``False`` -- the caller, folded into its own lobes: by hand from :meth:`prephaser_area_per_m`, or
  through :class:`~seqcraft.PhysicalDesignScope` with ``k_at_echo = target +``
  :meth:`echo_offset_per_m`.

The cosine's entrance and exit ramps are emitted under both.  They are part of the waveform, not
of the prephasing.
"""

from __future__ import annotations

from math import pi
from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ...design.events import knots_of
from ...design.logic import LogicBlock, flatten
from ...design.module import Module
from ...design.timing import TICKS_PER_SECOND, from_ticks, to_ticks
from ...errors import ConfigurationError, format_error
from .._support import ceil_raster, require_axis, require_count, require_positive
from .cartesian_line import CartesianLine

if TYPE_CHECKING:
    from pypulseq.opts import Opts

    from ...design.events import Event

__all__ = ['WaveReadout']

#: The two centring conventions, each with what it puts on target.
_CENTERINGS = {
    'sinzero': 'k at the echo sample on target on both Wave axes',
    'sinctr': 'the sine corkscrew centred on the line; the cosine echo sample on target',
}

#: How many raster steps the preparation window may grow by before the search gives up.
_PREPARATION_SEARCH_STEPS = 2000

#: Relative headroom for pypulseq's strict ``>`` limit comparison.  A ramp designed at exactly
#: ``max_slew`` can exceed it by one ulp in floating point.
_ULP_HEADROOM = 1e-9


class WaveReadout(Module):
    """
    One Cartesian readout with Wave gradients: a sine on `sine_axis`, a cosine on `cosine_axis`.

    Parameters
    ----------
    opts
        The scanner.  The only hardware constraint: every emitted gradient is within its
        ``max_grad`` and ``max_slew``, and its rasters are the ones used.
    fov_mm, matrix, bandwidth_hz_px, dwell_s
        The readout, exactly as :class:`~seqcraft.modules.CartesianLine` takes them.  `matrix` is
        reconstructed pixels; `dwell_s` is the ADC's own dwell.
    cycles
        Whole Wave periods over the ADC dwell window.  Any positive integer, odd or even.
    amplitude_mT_m
        The Wave amplitude ``G0``, mT/m, on both Wave axes.  A number is a prescription: if
        ``G0 > max_grad`` or ``G0 * omega > max_slew`` it raises, naming the largest that fits,
        and it is never reduced.  ``'max'`` (the default) is that largest value,
        ``min(max_grad, max_slew / omega)``.  :attr:`amplitude_mT_m` reads back what was used.
    oversampling
        ADC samples per reconstructed pixel, as for ``CartesianLine``.  ``4`` by default.
    axis
        The readout axis.
    sine_axis, cosine_axis
        The Wave axes, each different from `axis` and from each other.  Either may be ``None`` for
        a single-axis Wave; not both.
    centering
        ``'sinzero'`` (default) or ``'sinctr'``; see the module docstring.  ``'sinctr'`` needs a
        `sine_axis`.
    prephase
        Whether this module plays the moment its centring requires before the readout, on every
        axis it owns.  ``False`` leaves that to the caller; it never removes the cosine's ramps.
    prephaser_duration_s
        With ``prephase=True``, the time from the block's start to the readout gradient, which is
        when the readout's prephaser ends.  ``None`` is the shortest in which every preparation
        lobe is legal; shorter raises, naming it.  The Wave lobes end at the ADC start.
    tag
        Optional identity.

    Notes
    -----
    **The trajectory is reported on the ADC samples.**  :meth:`sample_times_s` is the sample time
    vector and :meth:`k_per_m` the ``(3, N)`` positions, relative to the encoding target, on that
    vector -- one row per logical axis, the readout's included.  Both assume the stated moment is
    in place, so they are the same for either `prephase`.

    **Single echo, full echo.**  Partial Fourier, readout trains and flow compensation are not
    parameters.  Flow compensation belongs to the repetition, which sees the Wave as a fixed
    contribution through :class:`~seqcraft.PhysicalDesignScope`.

    **One block.**  Each Wave axis is one event spanning the whole readout, so the compiler places
    everything this module emits in one Pulseq block.  Under ``prephase=True`` that includes the
    readout prephaser and lobe, which it merges into one readout-axis gradient, with a
    ``merge`` warning saying so.

    **No PNS or acoustic property is claimed** for any amplitude, including ``'max'``.

    Examples
    --------
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=80, grad_unit='mT/m', max_slew=200, slew_unit='T/m/s',
    ...          adc_dead_time=20e-6)
    >>> wave = WaveReadout(opts=o, fov_mm=250.0, matrix=250, dwell_s=5e-6, cycles=25,
    ...                    sine_axis='y', cosine_axis='x', axis='z', amplitude_mT_m=5.093)
    >>> wave.line.num_samples, round(wave.window_s * 1e3, 6)       # 250 x 4 samples, 5 ms
    (1000, 5.0)

    k at the echo sample is on target on both Wave axes, for odd cycles too:

    >>> k = wave.k_per_m()
    >>> echo = wave.line.echo_sample(0)
    >>> bool(abs(k[0, echo]) < 1e-9 and abs(k[1, echo]) < 1e-9)
    True
    """

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: float,
        matrix: int,
        cycles: int,
        bandwidth_hz_px: float | None = None,
        dwell_s: float | None = None,
        oversampling: int = 4,
        amplitude_mT_m: float | str = 'max',
        axis: str = 'x',
        sine_axis: str | None = 'y',
        cosine_axis: str | None = 'z',
        centering: str = 'sinzero',
        prephase: bool = True,
        prephaser_duration_s: float | None = None,
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.axis = require_axis(axis)
        self.sine_axis = None if sine_axis is None else require_axis(sine_axis)
        self.cosine_axis = None if cosine_axis is None else require_axis(cosine_axis)
        self._check_axes()
        self.cycles = require_count(cycles, 'cycles', hint='whole periods over the ADC window')
        self.centering = self._resolve_centering(centering)
        self.prephase = bool(prephase)

        readout = dict(opts=opts, fov_mm=fov_mm, matrix=matrix, bandwidth_hz_px=bandwidth_hz_px,
                       dwell_s=dwell_s, oversampling=oversampling, axis=self.axis)
        line = CartesianLine(**readout, prephase=self.prephase,
                             prephaser_duration_s=prephaser_duration_s)
        # The sampling does not depend on the preparation window, so it is fixed here.
        self._n, self._dwell, self._echo_index = (line.num_samples, line.dwell_s,
                                                  line.echo_sample(0))
        self._dt = to_ticks(float(opts.grad_raster_time))
        self._window = line.num_samples * to_ticks(line.dwell_s)
        #: The ADC dwell window the Wave is anchored to, seconds: ``num_samples * dwell_s``.
        self.window_s = from_ticks(self._window)
        #: Wave angular frequency, rad/s: ``2 * pi * cycles / window_s``.
        self.angular_frequency_rad_s = 2 * pi * self.cycles / self.window_s
        self._check_raster_density()
        self._g0 = self._resolve_amplitude(amplitude_mT_m)
        #: The Wave amplitude used, mT/m.
        self.amplitude_mT_m = self._g0 / float(opts.gamma) * 1e3
        self._ramp = (
            to_ticks(ceil_raster(self._g0 / float(opts.max_slew) * (1 + _ULP_HEADROOM),
                                 float(opts.grad_raster_time)))
            if self.cosine_axis else 0
        )

        adc_delay = to_ticks(float(line.adc.delay))
        if self.prephase:
            self._line_start = 0
            preparation, self._waves = self._design_prepared(
                to_ticks(line.prephaser_duration_s), adc_delay,
                fixed=prephaser_duration_s is not None,
            )
            if preparation != to_ticks(line.prephaser_duration_s):
                line = CartesianLine(**readout, prephase=True,
                                     prephaser_duration_s=from_ticks(preparation))
            self._t_a = preparation + adc_delay
        else:
            # The readout lobe starts later when the cosine's entrance ramp is longer than the
            # time between the lobe's start and the ADC's.
            self._line_start = self._raster_ceil(max(0, self._ramp - adc_delay))
            self._t_a = self._line_start + adc_delay
            self._waves = {
                wave_axis: self._emit(wave_axis, *self._mandatory(wave_axis, self._t_a))
                for wave_axis in self._wave_axes
            }
        #: The Cartesian readout this Wave readout plays: lobe, ADC and, under ``prephase=True``,
        #: the readout-axis prephaser.  ``num_samples``, ``dwell_s`` and ``echo_sample`` live here.
        self.line = line
        self._samples, self._echo = self._clock(self._t_a)
        self._mandatory_knots = {
            wave_axis: self._mandatory(wave_axis, self._t_a) for wave_axis in self._wave_axes
        }

    # ------------------------------------------------------------------ what it knows
    def time_to_echo(self) -> float:
        """Seconds from the start of this module's block to the echo sample."""
        return float(self._echo)

    def sample_times_s(self) -> np.ndarray:
        """
        Return every ADC sample's time from the start of this module's block, seconds.

        ``t_a + (n + 1/2) * dwell`` for ``n = 0 .. num_samples - 1``: the sample index and time
        vector the raw data are acquired on.
        """
        return self._clock(self._t_a)[0]

    def k_per_m(self) -> np.ndarray:
        """
        Return ``(3, num_samples)`` k at every ADC sample, 1/m, rows x, y, z.

        Relative to the encoding target the caller puts on each axis, with the stated moment
        (:meth:`prephaser_area_per_m`) in place; an axis this module leaves idle is a row of
        zeros.  Integrated exactly from the emitted knots of :meth:`build` and evaluated at
        :meth:`sample_times_s`.  The readout row is ``(n - echo) * dk_per_m / oversampling``.
        """
        out = np.zeros((3, self.line.num_samples))
        for start, event, _ in flatten(self.build()):
            if getattr(event, 'type', None) not in ('trap', 'grad'):
                continue
            out['xyz'.index(event.channel)] += _integral_at(self._samples,
                                                            *knots_of(event, start))
        if not self.prephase:
            for row, name in enumerate('xyz'):
                if name in self._owned_axes:
                    out[row] += self.prephaser_area_per_m(name)
        return out

    def prephaser_area_per_m(self, axis: str) -> float:
        """
        The moment, 1/m, that must be in place before this block on `axis` for its contract.

        **The requirement, stated whatever `prephase` is.**  ``prephase=True`` plays it inside
        the block; ``prephase=False`` leaves it to the caller.  On the readout axis it is
        ``CartesianLine``'s prephaser area; on a Wave axis it cancels the cosine's entrance ramp and
        the Wave up to the echo sample, plus the ``'sinctr'`` offset.  Zero on an axis the module
        leaves idle.
        """
        name = require_axis(axis)
        if name == self.axis:
            return -self.line.area_to_echo_per_m
        if name not in self._wave_axes:
            return 0.0
        return -self._functional(name, *self._mandatory_knots[name], (self._samples, self._echo))

    def echo_offset_per_m(self, axis: str) -> float:
        """
        k at the echo sample minus the encoding target on `axis`, 1/m, under the contract.

        What a :class:`~seqcraft.PhysicalDesignScope` composition adds to its ``k_at_echo``.  Zero
        everywhere under ``'sinzero'``; under ``'sinctr'``, on the sine axis, k at the echo minus
        the mean of k over the samples.
        """
        name = require_axis(axis)
        if self.centering != 'sinctr' or name != self.sine_axis:
            return 0.0
        times, amps = self._mandatory_knots[name]
        k = _integral_at(self._samples, times, amps)
        return float(_integral_at(np.array([self._echo]), times, amps)[0] - np.mean(k))

    def area_after_echo_per_m(self, axis: str) -> float:
        """
        What this block plays on `axis` after the echo sample, 1/m.

        The other half of a rewinder or spoiler's arithmetic: k leaves the block at
        ``target + echo_offset_per_m + area_after_echo_per_m``.
        """
        name = require_axis(axis)
        if name == self.axis:
            return self.line.area_after_echo_per_m
        if name not in self._wave_axes:
            return 0.0
        times, amps = self._mandatory_knots[name]
        ends = np.array([self._echo, times[-1]])
        before, total = _integral_at(ends, times, amps)
        return float(total - before)

    @property
    def prephaser_duration_s(self) -> float:
        """Seconds from the block's start to the readout gradient.  ``prephase=True`` only."""
        return self.line.prephaser_duration_s

    # ----------------------------------------------------------------------- assembly
    def build(
        self, *, acquire: bool = True, phase_deg: float = 0.0, offset_mm: float = 0.0,
    ) -> LogicBlock:
        """
        Return the readout, its ADC and the Wave gradients.

        Parameters
        ----------
        acquire
            ``False`` plays identical gradients with no ADC, for a dummy repetition.
        phase_deg
            Receiver phase, degrees: the excitation's phase, as for ``CartesianLine``.
        offset_mm
            FOV shift along the readout axis, as for ``CartesianLine``.
        """
        out = LogicBlock()
        out.add(from_ticks(self._line_start),
                self.line(acquire=acquire, phase_deg=phase_deg, offset_mm=offset_mm))
        for start, event in self._waves.values():
            out.add(start, event)
        return out

    # -------------------------------------------------------------------------- design
    @property
    def _wave_axes(self) -> tuple[str, ...]:
        return tuple(a for a in (self.sine_axis, self.cosine_axis) if a is not None)

    @property
    def _owned_axes(self) -> tuple[str, ...]:
        return (self.axis, *self._wave_axes)

    def _raster_ceil(self, ticks: int) -> int:
        return -(-ticks // self._dt) * self._dt

    def _mandatory(self, axis: str, t_a: int) -> tuple[np.ndarray, np.ndarray]:
        """Return the emitted knots without preparation: the Wave, and the cosine's ramps."""
        before = (np.array([t_a - self._ramp, t_a]), np.array([0.0, self._g0]))
        return self._realise(axis, t_a, before if axis == self.cosine_axis else None)

    def _clock(self, t_a: int) -> tuple[np.ndarray, float]:
        """Return the sample times and the echo sample's time for an ADC start of `t_a` ticks."""
        samples = from_ticks(t_a) + (np.arange(self._n) + 0.5) * self._dwell
        return samples, float(samples[self._echo_index])

    def _realise(
        self, axis: str, t_a: int, before: tuple[np.ndarray, np.ndarray] | None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Return ``(times_s, amps)``: the prescription on `axis` at every covering raster edge.

        `before` is the piecewise-linear prescription that ends at ``t_a`` (ticks, values), or
        ``None`` for rest.  The cosine's exit ramp is added here; the sine needs none.
        """
        cosine = axis == self.cosine_axis
        end = t_a + self._window
        after = (np.array([end, end + self._ramp]), np.array([self._g0, 0.0])) if cosine else None
        first = int(before[0][0]) if before is not None else t_a
        last = int(after[0][-1]) if after is not None else end
        edges = np.arange(first // self._dt, -(-last // self._dt) + 1, dtype=np.int64) * self._dt

        since = edges - t_a
        active = (since >= 0) & (since < self._window)
        phase = self.angular_frequency_rad_s * since[active] / TICKS_PER_SECOND
        amps = np.zeros(edges.size)
        amps[active] = self._g0 * (np.cos(phase) if cosine else np.sin(phase))
        if before is not None:
            early = since < 0
            amps[early] = np.interp(edges[early], before[0], before[1], left=0.0)
        if after is not None:
            late = since >= self._window
            amps[late] = np.interp(edges[late], after[0], after[1], right=0.0)
        return edges / TICKS_PER_SECOND, amps

    def _functional(self, axis: str, times: np.ndarray, amps: np.ndarray,
                    clock: tuple[np.ndarray, float]) -> float:
        """The k on `axis` that centring sets to zero, from emitted knots starting at rest."""
        samples, echo = clock
        if self.centering == 'sinctr' and axis == self.sine_axis:
            return float(np.mean(_integral_at(samples, times, amps)))
        return float(_integral_at(np.array([echo]), times, amps)[0])

    def _design_prepared(
        self, shortest: int, adc_delay: int, *, fixed: bool,
    ) -> tuple[int, dict[str, tuple[float, Event]]]:
        """
        Return the preparation window (ticks) and the emitted Wave events under ``prephase=True``.

        The window starts at `shortest` -- the readout prephaser's, or the one requested -- and
        grows by one raster until every Wave lobe is legal.  A requested (`fixed`) window that
        is too short raises, naming the shortest that works.
        """
        window = shortest
        for _ in range(_PREPARATION_SEARCH_STEPS):
            designed = self._try_prepared(window + adc_delay)
            if designed is not None:
                break
            window += self._dt
        else:
            msg = format_error(
                'no preparation window within the search limit makes every Wave lobe legal.',
                {'cycles': self.cycles, 'amplitude_mT_m': self.amplitude_mT_m},
                ['prephase=False, with a composition that plays prephaser_area_per_m'],
            )
            raise ConfigurationError(msg)
        if fixed and window != shortest:
            msg = format_error(
                f'the Wave preparation lobes do not fit in a '
                f'{from_ticks(shortest) * 1e6:.6g} us prephaser window.',
                {'prephaser_duration_s': from_ticks(shortest), 'cycles': self.cycles,
                 'amplitude_mT_m': self.amplitude_mT_m},
                [f'prephaser_duration_s={from_ticks(window):.6g} is the shortest that fits',
                 'or None, for the shortest'],
            )
            raise ConfigurationError(msg)
        return window, {axis: self._emit(axis, *knots) for axis, knots in designed.items()}

    def _try_prepared(self, t_a: int) -> dict[str, tuple[np.ndarray, np.ndarray]] | None:
        """
        Return each Wave axis's emitted knots for an ADC start of `t_a`, or ``None`` if illegal.

        Each axis gets one three-knot lobe from rest at the block's start, through a middle knot
        on the raster, to its value at the ADC start: ``0`` on the sine, ``G0`` on the cosine,
        which makes it the cosine's entrance ramp as well.  The emitted waveform is linear in the
        middle amplitude, so the centring condition is affine in it and one division solves it.
        """
        middle = int(round(t_a / 2 / self._dt)) * self._dt
        if middle <= 0 or middle >= t_a:
            return None
        clock = self._clock(t_a)
        shape = np.array([0, middle, t_a])
        out = {}
        for axis in self._wave_axes:
            end = self._g0 if axis == self.cosine_axis else 0.0
            times, base = self._realise(axis, t_a, (shape, np.array([0.0, 0.0, end])))
            unit = self._realise(axis, t_a, (shape, np.array([0.0, 1.0, end])))[1] - base
            middle_amplitude = (-self._functional(axis, times, base, clock)
                                / self._functional(axis, times, unit, clock))
            amps = base + middle_amplitude * unit
            if not self._legal(times, amps):
                return None
            out[axis] = (times, amps)
        return out

    def _legal(self, times: np.ndarray, amps: np.ndarray) -> bool:
        grad = float(np.max(np.abs(amps)))
        slew = float(np.max(np.abs(np.diff(amps) / np.diff(times))))
        return (grad <= self.opts.max_grad * (1 + _ULP_HEADROOM)
                and slew <= self.opts.max_slew * (1 + _ULP_HEADROOM))

    def _emit(self, axis: str, times: np.ndarray, amps: np.ndarray) -> tuple[float, Event]:
        """Return ``(node_time_s, event)``: the knots as one extended trapezoid."""
        times, amps = _drop_collinear(times, amps, tol=1e-12 * self._g0)
        start = float(times[0])
        event = pp.make_extended_trapezoid(
            channel=axis, times=times - start, amplitudes=amps, system=self.opts,
            max_grad=float(self.opts.max_grad) * (1 + _ULP_HEADROOM),
            max_slew=float(self.opts.max_slew) * (1 + _ULP_HEADROOM),
        )
        return start, event

    # ---------------------------------------------------------------------- validation
    def _check_axes(self) -> None:
        axes = [a for a in (self.axis, self.sine_axis, self.cosine_axis) if a is not None]
        if self.sine_axis is None and self.cosine_axis is None:
            msg = format_error('a Wave readout needs at least one Wave axis.',
                               {'sine_axis': None, 'cosine_axis': None},
                               ["sine_axis='y', cosine_axis='z' with axis='x'"])
            raise ConfigurationError(msg)
        if len(set(axes)) != len(axes):
            msg = format_error(
                'the readout axis and the Wave axes must all differ.',
                {'axis': self.axis, 'sine_axis': self.sine_axis, 'cosine_axis': self.cosine_axis},
                ["axis='x', sine_axis='y', cosine_axis='z' (transverse)",
                 "axis='z', sine_axis='y', cosine_axis='x' (sagittal)"],
            )
            raise ConfigurationError(msg)

    def _resolve_centering(self, centering: str) -> str:
        if centering not in _CENTERINGS:
            msg = format_error(
                f'centering={centering!r} is not one of {sorted(_CENTERINGS)}.',
                {'centering': centering},
                [f'{name!r}: {meaning}' for name, meaning in _CENTERINGS.items()],
            )
            raise ConfigurationError(msg)
        if centering == 'sinctr' and self.sine_axis is None:
            msg = format_error("centering='sinctr' centres the sine axis, and there is none.",
                               {'centering': centering, 'sine_axis': None},
                               ["centering='sinzero'", 'or pass a sine_axis'])
            raise ConfigurationError(msg)
        return centering

    def _check_raster_density(self) -> None:
        """Refuse a Wave period the gradient raster cannot represent: two rasters or fewer."""
        period = self._window / self.cycles
        if period <= 2 * self._dt:
            msg = format_error(
                f'{self.cycles} cycles over a {self.window_s * 1e3:.4g} ms window is a '
                f'{period / TICKS_PER_SECOND * 1e6:.4g} us period, which the '
                f'{from_ticks(self._dt) * 1e6:.4g} us gradient raster cannot represent.',
                {'cycles': self.cycles, 'window_s': self.window_s},
                ['fewer cycles', 'a longer window (lower bandwidth_hz_px)'],
            )
            raise ConfigurationError(msg)

    def _resolve_amplitude(self, amplitude_mT_m: float | str) -> float:
        """Return ``G0`` in Hz/m: the largest legal one for ``'max'``, else the prescription."""
        max_grad, max_slew = float(self.opts.max_grad), float(self.opts.max_slew)
        limit = min(max_grad, max_slew / self.angular_frequency_rad_s)
        to_mT = 1e3 / float(self.opts.gamma)
        if isinstance(amplitude_mT_m, str):
            if amplitude_mT_m != 'max':
                msg = format_error(f'amplitude_mT_m={amplitude_mT_m!r} is not a number or "max".',
                                   {'amplitude_mT_m': amplitude_mT_m},
                                   ['amplitude_mT_m="max"', 'or an amplitude in mT/m'])
                raise ConfigurationError(msg)
            return limit
        wanted = require_positive(amplitude_mT_m, 'amplitude_mT_m') / to_mT
        if wanted > limit * (1 + _ULP_HEADROOM):
            bound = 'max_grad' if max_grad <= max_slew / self.angular_frequency_rad_s else 'max_slew'
            msg = format_error(
                f'a {amplitude_mT_m:g} mT/m Wave at {self.cycles} cycles over '
                f'{self.window_s * 1e3:.4g} ms exceeds {bound}; the largest that fits is '
                f'{limit * to_mT:.6g} mT/m.',
                {'amplitude_mT_m': amplitude_mT_m, 'cycles': self.cycles,
                 'window_s': self.window_s, 'max_grad_mT_m': max_grad * to_mT,
                 'max_slew_T_m_s': max_slew * to_mT},
                [f'amplitude_mT_m <= {limit * to_mT:.6g}', 'amplitude_mT_m="max"',
                 'fewer cycles, or a longer window, if the slew limit binds'],
            )
            raise ConfigurationError(msg)
        return wanted


def _integral_at(times: np.ndarray, knot_t: np.ndarray, knot_g: np.ndarray) -> np.ndarray:
    """
    Return the exact integral of a piecewise-linear waveform from its first knot to each time.

    Zero before the first knot and the total after the last: the cumulative trapezoid over the
    knots, plus the partial trapezoid into the interval each time falls in.
    """
    cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (knot_g[1:] + knot_g[:-1])
                                                  * np.diff(knot_t))])
    index = np.clip(np.searchsorted(knot_t, times, side='right') - 1, 0, len(knot_t) - 2)
    span = knot_t[index + 1] - knot_t[index]
    into = np.clip(times - knot_t[index], 0.0, span)
    slope = (knot_g[index + 1] - knot_g[index]) / span
    return cumulative[index] + into * (knot_g[index] + 0.5 * slope * into)


def _drop_collinear(times: np.ndarray, amps: np.ndarray, *, tol: float
                    ) -> tuple[np.ndarray, np.ndarray]:
    """Drop interior knots on a straight line through their neighbours (uniform `times`)."""
    keep = np.ones(times.size, dtype=bool)
    keep[1:-1] = np.abs(amps[2:] - 2 * amps[1:-1] + amps[:-2]) > tol
    return times[keep], amps[keep]
