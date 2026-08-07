"""
Echo-planar readout train: one shot collects many k-space lines from a single echo.

The train is **one gradient event per axis**, not one per echo.  That is the first thing to know
about this module, because it is what the compiler wants: placing 96 abutting trapezoids on x puts
the tail of echo *n* and the head of echo *n+1* into the same compiled block, and the compiler warns
whenever a block holds more than one gradient on an axis -- 96 warnings per shot, 1824 for a
19-volume acquisition, all of them describing correct output.  One
:func:`pypulseq.make_extended_trapezoid` per axis gives every block exactly one piece to split, and
the compiled result is identical.

Why EPI for diffusion
---------------------
A single echo collects the whole plane, as a spiral does, but on a Cartesian grid -- so the
reconstruction is a gridding problem rather than a resampling one, and the acquisition reaches a
useful resolution in one shot.  The cost is the direction the errors take.  Off-resonance does not
blur an EPI image, it **displaces** it: the phase accrued between successive lines is
``2 pi df * echo_spacing``, which is linear in ``ky``, and a linear phase in k-space is a shift in
image space.  :attr:`pe_bandwidth_per_pixel_hz` is how far -- 20 Hz per pixel on the sequence the
examples build, so 100 Hz of off-resonance moves a voxel five pixels.

That is also why the same reconstruction operator serves both.  A spiral's blur and an EPI's
displacement are the *same* ``exp(-2i pi df t)`` term; only the way ``t`` maps onto k-space differs.

What sets the echo spacing
--------------------------
For a lobe of fixed area the minimum-time trapezoid is a **triangle** at ``G = sqrt(A S)``, whenever
that is under the amplitude limit.  At 240 mm and a 128 matrix it is 48.8 mT/m against a 170 mT/m
limit, so **the amplitude limit never binds and the slew limit always does**.  The echo spacing is
``2 sqrt(A / S)``, and the only levers on it are the slew rate, ``k_max`` and
:attr:`partial_echo`.  Reaching for a readout duration is reaching for the wrong knob, which is why
`flat_time_us` is an override rather than the parameter -- and why an oscillating readout is the
worst case for peripheral nerve stimulation.  Measure it with
:meth:`~seqcraft.core.compiler.CompiledSequence.pns` against your own hardware descriptor; on a
200 T/m/s system this train is not automatically runnable.

Ramp sampling is not optional here
----------------------------------
A minimum-time lobe has no flat top at all, so sampling only the flat top would sample nothing.
Every sample is on a ramp, the trajectory is non-uniform in ``kx``, and the reconstruction must grid
it -- which is what makes the "non-Cartesian" reconstruction in ``examples/lib`` the right one for a
Cartesian sequence.  :meth:`EPIReadout.trajectory` gives the sample positions.

Examples
--------
>>> import seqcraft as sc
>>> system = sc.System.preset('cima_x').derate('epi', grad=0.85, slew=1.0)
>>> ro = EPIReadout(system, fov_ro_mm=240, matrix_ro=128, fov_pe_mm=240, matrix_pe=128,
...                 partial_fourier_pe=0.75, regime='epi')
>>> ro.n_echoes, round(ro.echo_spacing * 1e6)
(96, 520)
>>> ro.echo_index                                    # ky = 0 is not the middle of the train
32
>>> round(ro.train_duration * 1e3, 2), round(ro.time_to_echo * 1e3, 2)
(49.92, 17.26)
>>> round(ro.pe_bandwidth_per_pixel_hz, 2)           # the distortion figure
20.03
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ...core import events as ev
from ...core.errors import ConfigurationError, format_error
from ...core.geometry import round_half_up
from ...core.logic import LogicBlock
from ...core.module import Module
from ...core.timing import Raster
from ...core.units import convert
from ...core.validate import Range, require_in, require_in_range, require_int_in, require_positive

if TYPE_CHECKING:
    from types import SimpleNamespace

    from ...core.system import System

__all__ = ['EPIReadout']

_FOV_RANGE = Range(0.5, 2000.0, 'mm', ((1e3, 'm'), (10.0, 'cm')))


def _unit_moment(t: float, ramp: float, flat: float) -> float:
    """
    Area under a **unit-amplitude** trapezoid from its own start to `t`.

    Amplitude factors out of every question this module asks -- where the prephaser has to leave k,
    where k crosses zero, what the sampled extent is -- so the shape is integrated once at unit
    amplitude and scaled.  That is what keeps the prephaser and :attr:`EPIReadout.time_to_echo`
    from disagreeing, the failure `CartesianLine` documents for its own ramp.

    Examples
    --------
    >>> _unit_moment(0.0, 1e-4, 0.0)
    0.0
    >>> _unit_moment(1e-4, 1e-4, 0.0) == 0.5e-4          # half the rising ramp
    True
    >>> _unit_moment(2e-4, 1e-4, 0.0) == 1e-4            # the whole triangle: ramp * 1
    True
    >>> _unit_moment(3e-4, 1e-4, 1e-4) == 2e-4           # with a flat top: ramp + flat
    True
    """
    if t <= 0.0:
        return 0.0
    if t < ramp:
        return 0.5 * t * t / ramp
    if t <= ramp + flat:
        return 0.5 * ramp + (t - ramp)
    past = min(t - ramp - flat, ramp)
    return 0.5 * ramp + flat + past * (1.0 - past / (2.0 * ramp))


def _unit_moment_inverse(value: float, ramp: float, flat: float) -> float:
    """
    Return the time at which :func:`_unit_moment` reaches `value`.

    Monotonic, so the inverse is well defined; solved in closed form on whichever of the three
    segments contains it rather than by search, so ``time_to_echo`` is exact.

    Examples
    --------
    >>> round(_unit_moment_inverse(0.5e-4, 1e-4, 0.0) * 1e6, 9)
    100.0
    >>> t = _unit_moment_inverse(0.8e-4, 1e-4, 0.0)
    >>> abs(_unit_moment(t, 1e-4, 0.0) - 0.8e-4) < 1e-18
    True
    """
    if value <= 0.0:
        return 0.0
    if value < 0.5 * ramp:
        return math.sqrt(2.0 * value * ramp)
    if value <= 0.5 * ramp + flat:
        return ramp + (value - 0.5 * ramp)
    rest = value - 0.5 * ramp - flat
    # On the falling ramp, m = ramp/2 + flat + p - p^2/(2 ramp); the root inside [0, ramp].
    disc = max(0.0, ramp * ramp - 2.0 * ramp * rest)
    return ramp + flat + (ramp - math.sqrt(disc))


class EPIReadout(Module):
    """
    Echo-planar readout train, with its own prephasers, blips and labels.

    Parameters
    ----------
    system
        The scanner.
    fov_ro_mm, matrix_ro
        Readout geometry.  Sets ``k_max = matrix / (2 FOV)``.
    fov_pe_mm, matrix_pe
        Phase-encode geometry.  Sets the blip area, ``n_shots / FOV``.
    n_shots
        Interleaved segments.  Shot *s* takes every ``n_shots``-th line of the table, so the blip
        is ``n_shots`` times larger and the train ``n_shots`` times shorter.  Playing **one** shot
        of an ``n_shots`` design is exactly an ``R = n_shots`` accelerated acquisition, which is
        how in-plane acceleration is expressed here rather than as a separate parameter.
    partial_fourier_pe
        Fraction of phase-encode k-space acquired, in ``(0.5, 1.0]``, truncating the **low** end as
        :class:`~seqcraft.core.geometry.Geometry` does.  This is the parameter that buys TE: at
        0.75 the ``ky = 0`` echo arrives a third of the way into the train instead of halfway.

        Snapped so the line count divides `n_shots`, always by acquiring *more* lines rather than
        fewer, and the attribute is then updated to the fraction achieved.
    partial_echo
        Fraction of the leading half of **readout** k-space acquired, in ``(0.5, 1.0]``.  Shortens
        every echo and so the whole train.  Resolution is unchanged; what is given up is the
        conjugate-symmetric part, on every line.
    adc_dwell_ns
        Sampling interval, nanoseconds.  ``None`` takes the longest dwell that satisfies Nyquist
        along the readout at the lobe's peak -- ``G * dwell <= 1 / FOV`` -- which is the binding
        constraint when sampling the ramps.
    flat_time_us
        Force a flat top on the readout lobe.  ``None`` solves for the shortest echo, which on any
        realistic system is a pure triangle.  Use this to derate deliberately: a longer, lower
        lobe is quieter and easier on peripheral nerve stimulation.
    blip_duration_us
        Force the blip duration.  ``None`` takes the shortest legal blip: long enough for its own
        area at the slew limit, long enough that ``blip / 2`` covers the ADC dead time, and an
        **even** number of gradient rasters so that ``blip / 2`` lands on the raster too.
    axes
        Which logical axes carry readout and phase encode.
    prephase
        Include the prephasers in :meth:`build`.  As :class:`CartesianLine`, this is a constructor
        argument and not a build one, because dropping them changes the block's duration and would
        silently invalidate :attr:`time_to_echo`.
    labels
        Emit ``LIN`` per echo, ``SEG`` per shot, and ``REV`` on the reversed echoes.
    regime
        Which limit regime to design against.

    Properties
    ----------
    echo_spacing, n_echoes, echo_index, train_duration
    samples_per_echo, n_samples, dwell_s, sample_offset
    time_to_echo, duration_after_echo, duration, prephase_duration
    k_max_ro_per_m, dk_pe_per_m, resolution_mm
    bandwidth_per_pixel_hz, pe_bandwidth_per_pixel_hz, total_readout_time_s
    peak_amplitude_Hz_per_m, blip_duration, ramp_time, flat_time

    Build arguments
    ---------------
    shot : int, default 0
        Which interleaved segment.  Changes the prephaser and the labels, never the timing.
    pe_polarity : {1, -1}, default 1
        Reverse the phase-encode direction -- the blip-up/blip-down pair a distortion correction
        needs.  Both the blips and the PE prephaser flip together, so the duration is untouched.
    rf_phase_rad : float, default 0.0
        Receiver phase, which must be the same value passed to the excitation.

    Notes
    -----
    **The blip is outside every ADC window, and that costs k-space.**  The blip is centred on the
    junction where the readout gradient crosses zero, and the ADC skips the first and last
    ``blip / 2`` so that ``ky`` is *constant* while a line is sampled -- measured drift 0.0 1/m.
    What must therefore cover the k-space extent is the **sampled** area, not the lobe's own:

    .. code-block:: text

        sampled = G * (m(t_end) - m(t_start))  =  (1 + partial_echo) * k_max

    On the example sequence the lobe's own area is 540.53 1/m against 533.33 1/m of sampled extent.
    A prephaser of minus ``k_max`` -- the tempting value -- would leave every shot displaced in
    k-space by the difference, 1.7 ``dk``, and a k-space offset is a linear phase ramp across the
    image: the magnitude image looks perfect and every phase-derived quantity is wrong.

    **The sampling offset is the ADC event's own delay**, and the ADC node sits at the lobe start.
    Two traps at once.  An ADC node off the *block* raster is snapped by the compiler with a
    ``raster`` warning, which moves the ADC against the gradient by up to half a raster -- 10 1/m
    of ``kx`` here; placing the node at the lobe start puts it on the raster by construction.  And
    ``pp.make_adc`` silently raises a ``delay`` below ``adc_dead_time`` up to it, while seqcraft
    *preserves* an event's own delay rather than folding it away -- so carrying the offset in the
    node time instead made the two add, and sampling began 40 us into the lobe rather than 30.

    **`echo_index` is shot-independent.**  It is ``p // n_shots``, where ``p`` is the centre line's
    position in the full table, so every shot has identical timing -- which is required, since a
    build argument may not change a timing property, and is also the physics: segmented EPI plays
    every shot at the same TE so they carry the same T2 weighting.  For the one shot that samples
    the centre line it is exactly the ``ky = 0`` echo; for the others it is within
    ``(n_shots - 1) * dk`` of it.

    Examples
    --------
    >>> import seqcraft as sc
    >>> system = sc.System.preset('cima_x').derate('epi', grad=0.85, slew=1.0)
    >>> ro = EPIReadout(system, fov_ro_mm=240, matrix_ro=128, fov_pe_mm=240, matrix_pe=128,
    ...                 partial_fourier_pe=0.75, regime='epi')
    >>> ro.samples_per_echo, round(ro.dwell_s * 1e9)
    (228, 2000)
    >>> tuple(round(v, 3) for v in ro.resolution_mm)
    (1.875, 1.875)

    Partial Fourier buys TE, and nothing else moves:

    >>> full = EPIReadout(system, fov_ro_mm=240, matrix_ro=128, fov_pe_mm=240, matrix_pe=128,
    ...                   regime='epi')
    >>> full.n_echoes, full.echo_index
    (128, 64)
    >>> round((full.time_to_echo - ro.time_to_echo) * 1e3, 2)     # ms of TE saved
    16.64
    >>> full.resolution_mm == ro.resolution_mm and full.echo_spacing == ro.echo_spacing
    True

    Two shots halve the train and double the phase-encode bandwidth:

    >>> two = EPIReadout(system, fov_ro_mm=240, matrix_ro=128, fov_pe_mm=240, matrix_pe=128,
    ...                  partial_fourier_pe=0.75, n_shots=2, regime='epi')
    >>> two.n_echoes, two.echo_index
    (48, 16)
    >>> two.pe_bandwidth_per_pixel_hz > ro.pe_bandwidth_per_pixel_hz
    True
    """

    def __init__(
        self,
        system: System,
        *,
        fov_ro_mm: float,
        matrix_ro: int,
        fov_pe_mm: float,
        matrix_pe: int,
        n_shots: int = 1,
        partial_fourier_pe: float = 1.0,
        partial_echo: float = 1.0,
        adc_dwell_ns: float | None = None,
        flat_time_us: float | None = None,
        blip_duration_us: float | None = None,
        axes: tuple[str, str] = ('x', 'y'),
        prephase: bool = True,
        labels: bool = True,
        regime: str = 'default',
    ) -> None:
        super().__init__(system, regime=regime)
        self.fov_ro_mm = float(fov_ro_mm)
        self.matrix_ro = int(matrix_ro)
        self.fov_pe_mm = float(fov_pe_mm)
        self.matrix_pe = int(matrix_pe)
        self.n_shots = int(n_shots)
        self.partial_fourier_pe = float(partial_fourier_pe)
        self.partial_echo = float(partial_echo)
        self.adc_dwell_ns = None if adc_dwell_ns is None else float(adc_dwell_ns)
        self.flat_time_us = None if flat_time_us is None else float(flat_time_us)
        self.blip_duration_us = None if blip_duration_us is None else float(blip_duration_us)
        self.axes = tuple(axes)
        self.prephase = bool(prephase)
        self.labels = bool(labels)

        require_positive(self, 'fov_ro_mm', 'fov_pe_mm')
        require_in_range(self, 'fov_ro_mm', _FOV_RANGE.lo, _FOV_RANGE.hi, unit='mm')
        require_in_range(self, 'fov_pe_mm', _FOV_RANGE.lo, _FOV_RANGE.hi, unit='mm')
        require_in_range(self, 'partial_fourier_pe', 0.5, 1.0, unit='')
        require_in_range(self, 'partial_echo', 0.5, 1.0, unit='')
        require_int_in(self, 'matrix_ro', lo=4, hi=8192)
        require_int_in(self, 'matrix_pe', lo=2, hi=8192)
        require_int_in(self, 'n_shots', lo=1, hi=512)
        for axis in self.axes:
            require_in(type('_A', (), {'axis': axis})(), 'axis', ('x', 'y', 'z'))
        if self.axes[0] == self.axes[1]:
            msg = format_error(
                'the readout and phase-encode axes must differ.',
                {'axes': list(self.axes)},
            )
            raise ConfigurationError(msg)

        self._plan_lines()
        self._design_blip()
        self._design_train()
        self._design_prephasers()

    # ------------------------------------------------------------------------------ the plan
    def _plan_lines(self) -> None:
        """
        Choose the phase-encode table, and where in it the ``ky = 0`` echo falls.

        Partial Fourier truncates the low end, and the count is then nudged **up** until it
        divides `n_shots` -- always acquiring more lines, never fewer, so the request is a floor
        rather than a promise.  ``partial_fourier_pe`` is rewritten to what that achieved, the
        convention ``CartesianLine.partial_echo`` already follows.
        """
        if self.matrix_pe % self.n_shots:
            msg = format_error(
                f'{self.matrix_pe} phase-encode lines cannot be split across '
                f'{self.n_shots} shots.',
                {'matrix_pe': self.matrix_pe, 'n_shots': self.n_shots},
                [
                    f'use a matrix_pe that is a multiple of {self.n_shots}',
                    f'or n_shots dividing {self.matrix_pe}',
                ],
            )
            raise ConfigurationError(msg)

        first = self.matrix_pe - round_half_up(self.partial_fourier_pe * self.matrix_pe)
        while first > 0 and (self.matrix_pe - first) % self.n_shots:
            first -= 1
        self._first_line = first
        self._n_echoes = (self.matrix_pe - first) // self.n_shots
        self.partial_fourier_pe = (self.matrix_pe - first) / self.matrix_pe
        # Position of the k-space centre in the full table, then within one shot's own train.
        self._echo_index = (self.centre_line - first) // self.n_shots

    def _design_blip(self) -> None:
        """
        Design the phase-encode blip, and say which constraint set its length.

        Three floors, and which one binds is worth knowing rather than assuming.  On a Cima.X at
        240 mm the blip's **own slew limit** gives 60 us while the ADC dead time asks for only 20;
        on a lower-slew system, or with the larger blips a segmented acquisition needs, the
        ordering can swap.

        Quantised to an **even** number of gradient rasters, because the blip is centred on the
        junction and ``make_extended_trapezoid`` rejects a knot off the raster -- the first version
        of this failed with *"The last time point must be on a gradient raster"* on a 50 us blip.
        """
        raster = self.system.grad_raster
        quantum = Raster(raster.at(2), 'blip')
        area = self.dk_pe_per_m * self.n_shots
        slew_floor = 2.0 * math.sqrt(area / float(self.opts.max_slew))
        amplitude_floor = 2.0 * area / float(self.opts.max_grad)
        dead_floor = 2.0 * float(self.opts.adc_dead_time)
        self._blip_floors = {
            'slew': quantum.ceil(slew_floor),
            'amplitude': quantum.ceil(amplitude_floor),
            'adc_dead_time': quantum.ceil(dead_floor),
        }
        shortest = max(self._blip_floors.values())
        if self.blip_duration_us is None:
            self._blip = shortest
        else:
            self._blip = quantum.ceil(convert(self.blip_duration_us, 'us', 's'))
            if self._blip < shortest - 1e-15:
                binding = max(self._blip_floors, key=lambda k: self._blip_floors[k])
                msg = format_error(
                    f'a {convert(self._blip, "s", "us"):.0f} us blip is shorter than the '
                    f'{convert(shortest, "s", "us"):.0f} us its '
                    f'{binding.replace("_", " ")} limit needs.',
                    {
                        f'{name} floor_us': convert(value, 's', 'us')
                        for name, value in self._blip_floors.items()
                    },
                    [
                        f'use blip_duration_us >= {convert(shortest, "s", "us"):.0f}',
                        'or leave blip_duration_us=None to take the shortest legal blip',
                        'a larger n_shots needs a longer blip, since its area scales with it',
                    ],
                )
                raise ConfigurationError(msg)
        self._blip_area = area

    def _design_train(self) -> None:
        """
        Solve the readout lobe, the echo spacing and the sampling, together.

        They are coupled: the lobe's amplitude follows from the *sampled* extent, which depends on
        where sampling starts, which depends on how many samples fit in the echo spacing, which is
        set by the lobe.  Each step below is exact; only the search over `ramp` iterates.

        The search goes **upward from the first feasible ramp** rather than to a fix point.  A plain
        fix point oscillates -- 250 us needs a G that 250 us cannot reach, and rounding up to 260
        lowers the G it needs, so the two alternate forever.  Feasibility is monotone in `ramp`
        (``G`` falls roughly as ``1/ramp`` while ``S ramp`` rises), so the first ramp that fits is
        the shortest.

        Beyond ``G_max / S`` a longer ramp cannot help, because the lobe is then amplitude-limited
        rather than slew-limited: a flat top grows instead, which is the classic trapezoid readout.
        """
        raster = self.system.grad_raster
        adc_raster = self.system.adc_raster
        s_max = float(self.opts.max_slew)
        target = (1.0 + self.partial_echo) * self.k_max_ro_per_m
        fov_ro_m = convert(self.fov_ro_mm, 'mm', 'm')
        forced_flat = (
            None if self.flat_time_us is None
            else raster.ceil(convert(self.flat_time_us, 'us', 's'))
        )
        dwell = (
            None if self.adc_dwell_ns is None
            else adc_raster.nearest(convert(self.adc_dwell_ns, 'ns', 's'))
        )
        if dwell is not None and dwell <= 0.0:
            msg = format_error(
                f'adc_dwell_ns={self.adc_dwell_ns:g} rounds to zero on the '
                f'{adc_raster.dt * 1e9:.0f} ns ADC raster.', {},
            )
            raise ConfigurationError(msg)

        ramp_floor = max(raster.ceil(math.sqrt(target / s_max)), raster.dt)
        if dwell is not None:
            solved = self._solve_lobe(target, dwell, forced_flat, ramp_floor)
        else:
            solved = self._longest_legal_dwell(target, forced_flat, ramp_floor, fov_ro_m)

        self._ramp = float(solved['ramp'])
        self._flat = float(solved['flat'])
        self._amplitude = float(solved['amplitude'])
        self._echo_spacing = float(solved['echo_spacing'])
        self._offset = float(solved['offset'])
        self._dwell = float(solved['dwell'])
        # The solve carries everything as floats so one dict can hold it; the sample count is the one
        # entry that is genuinely an integer, and callers index arrays with it.
        self._samples = int(solved['samples'])

        # Nyquist along the readout, at the apex where k moves fastest.  Violating it aliases along
        # the readout direction, which looks like a fold rather than like a sampling mistake.
        bound = 1.0 / (fov_ro_m * self._amplitude)
        if self._dwell > bound * (1 + 1e-9):
            msg = format_error(
                f'a {self._dwell * 1e9:.0f} ns dwell aliases along the readout: at '
                f'{self.system.convert(self._amplitude, "Hz/m", "mT/m"):.1f} mT/m the sample '
                f'spacing is {self._amplitude * self._dwell:.3f} 1/m against a Nyquist limit of '
                f'{1.0 / fov_ro_m:.3f}.',
                {'bound_ns': adc_raster.floor(bound) * 1e9},
                [
                    f'use adc_dwell_ns <= {adc_raster.floor(bound) * 1e9:.0f}',
                    'or leave adc_dwell_ns=None to take the longest dwell that fits',
                    'or lengthen the lobe with flat_time_us, which lowers its amplitude',
                ],
            )
            raise ConfigurationError(msg)

        limit = int(getattr(self.opts, 'adc_samples_limit', 0) or 0)
        if limit and self._samples > limit:
            msg = format_error(
                f'{self._samples} samples per echo, above the {limit} the interpreter accepts.',
                {'dwell_ns': self._dwell * 1e9},
                ['raise adc_dwell_ns', 'or reduce matrix_ro'],
            )
            raise ConfigurationError(msg)

        self.adc = pp.make_adc(
            num_samples=self._samples, dwell=self._dwell, delay=self._offset, system=self.opts
        )
        if abs(float(self.adc.delay) - self._offset) > 1e-12:  # pragma: no cover - guard
            msg = format_error(
                'pypulseq moved the ADC delay, so sampling would not start where the design says.',
                {'asked_us': self._offset * 1e6, 'got_us': float(self.adc.delay) * 1e6},
            )
            raise ConfigurationError(msg)

        self._gx = pp.make_extended_trapezoid(
            channel=self.axes[0], amplitudes=self._readout_amplitudes(),
            times=self._readout_times(), system=self.opts,
        )
        self._gy = pp.make_extended_trapezoid(
            channel=self.axes[1], amplitudes=self._blip_amplitudes(),
            times=self._blip_times(), system=self.opts,
        )

    def _longest_legal_dwell(
        self,
        target: float,
        forced_flat: float | None,
        ramp_floor: float,
        fov_ro_m: float,
    ) -> dict[str, float]:
        """
        Return the design that samples the **most** of each echo without aliasing.

        Longest, because past Nyquist a finer dwell buys no resolution and costs data: the sample
        count is what the reconstruction and the simulator both scale with, and it is the *k-space
        extent* that sets resolution, which the dwell does not touch.

        The alternative objective -- maximise the sampled span, ``samples * dwell``, so the gradient
        amplitude for a given extent is as low as possible -- is measurably the wrong trade.  On the
        128-matrix single-shot design it picks 1000 ns and 460 samples per echo against 2000 ns and
        228, doubling the data to buy 4 us more span and **0.2 %** less gradient.

        Three constraints interact, and none is monotone in the dwell, so the choice is searched
        rather than derived:

        * Nyquist along the readout caps the dwell at ``1 / (FOV * G)``, and ``G`` itself depends on
          how much of the echo ends up sampled;
        * the sample count must be a multiple of the ADC divisor;
        * the leftover between the window and the samples must split evenly onto the RF raster, so
          that the sampled span stays centred on the lobe -- see :meth:`_sampling` for what an
          asymmetric one costs.
        """
        adc_raster = self.system.adc_raster
        s_max = float(self.opts.max_slew)

        def legal(candidate: float) -> dict[str, float] | None:
            if candidate <= 0.0:
                return None
            found = self._solve_lobe(target, candidate, forced_flat, ramp_floor)
            bound = 1.0 / (fov_ro_m * found['amplitude'])
            return found if candidate <= bound * (1 + 1e-12) else None

        # Seed above the answer: the smallest amplitude any lobe of this extent could need is the
        # corner-free triangle's, so the dwell it permits is an upper bound on the useful range.
        top = max(adc_raster.count(adc_raster.ceil(1.0 / (fov_ro_m * math.sqrt(target * s_max)))), 1)
        best: dict[str, float] | None = None
        for ticks in range(top, max(top - 4096, 0), -1):
            best = legal(adc_raster.at(ticks))
            if best is not None:
                break
        if best is None:
            msg = format_error(
                'no ADC dwell on this raster samples the readout without aliasing.',
                {
                    'fov_ro_mm': self.fov_ro_mm,
                    'matrix_ro': self.matrix_ro,
                    'adc_raster_ns': adc_raster.dt * 1e9,
                },
                ['reduce matrix_ro', 'increase fov_ro_mm', 'lengthen the lobe with flat_time_us'],
            )
            raise ConfigurationError(msg)
        return best

    def _solve_lobe(
        self,
        target: float,
        dwell: float,
        forced_flat: float | None,
        ramp_floor: float,
    ) -> dict[str, float]:
        """Return the shortest lobe delivering `target` of sampled extent at this `dwell`."""
        raster = self.system.grad_raster
        g_max, s_max = float(self.opts.max_grad), float(self.opts.max_slew)
        ramp_cap = raster.ceil(g_max / s_max)
        flat = 0.0 if forced_flat is None else forced_flat

        ramp = ramp_floor
        while True:
            trial = self._sampling(ramp, flat, dwell, target)
            if trial is not None and trial['amplitude'] <= min(g_max, s_max * ramp) * (1 + 1e-12):
                return trial
            if forced_flat is None and ramp >= ramp_cap:
                break
            if ramp > ramp_cap * 64.0:  # pragma: no cover - guard against a runaway search
                break
            ramp = raster.at(raster.count(ramp) + 1)

        # Amplitude-limited: hold the ramp at its full-amplitude value and grow a flat top.
        ramp = ramp_cap
        for extra in range(1_000_000):
            trial = self._sampling(ramp, raster.at(extra), dwell, target)
            if trial is not None and trial['amplitude'] <= g_max * (1 + 1e-12):
                return trial
        msg = format_error(  # pragma: no cover - needs contradictory limits
            'could not fit an EPI echo inside the gradient limits.',
            {'fov_ro_mm': self.fov_ro_mm, 'matrix_ro': self.matrix_ro},
            ['reduce matrix_ro', 'increase fov_ro_mm', 'design against a stronger regime'],
        )
        raise ConfigurationError(msg)

    def _sampling(
        self,
        ramp: float,
        flat: float,
        dwell: float,
        target: float,
    ) -> dict[str, float] | None:
        """
        Fit the ADC into one echo of this shape, and return the amplitude that shape then needs.

        Returns ``None`` when no sample group fits, which is how the ramp search knows to keep
        looking.  The leftover between the window and the samples is split between the two ends, so
        the sampled span stays centred on the lobe -- and the amplitude is solved from the span
        actually achieved rather than from an assumed symmetry, so a 100 ns rounding cannot leave
        the extent short.
        """
        raster = self.system.grad_raster
        divisor = self.system.adc_samples_divisor
        echo_spacing = raster.ceil(2.0 * ramp + flat)
        half = self._blip / 2.0
        window = echo_spacing - 2.0 * half
        if window <= 0.0:
            return None
        # The sampled span has to be **symmetric about the lobe centre**, so the leftover between the
        # window and the samples must split evenly onto the RF raster.  Otherwise the leading and
        # trailing unsampled corners differ, they no longer cancel between a forward echo and the
        # reversed one after it, and the reversed echoes' kx drifts: measured at 800 ns of asymmetry
        # on a two-shot design, which is 0.28 1/m -- 6.8 % of dk, and invisible in a single-shot
        # sequence, where the leftover happened to be even.
        quantum = 2.0 * self.system.rf_raster.dt
        samples = int(window / dwell + 1e-9) // divisor * divisor
        while samples >= divisor:
            leftover = window - samples * dwell
            if abs(leftover / quantum - round(leftover / quantum)) < 1e-9:
                break
            samples -= divisor
        if samples < divisor:
            return None
        # Three rasters conspire on the offset, and each was found by being caught:
        #
        #   * the **RF** raster, 1 us, because pypulseq's own timing check requires an ADC `delay`
        #     to land on it -- not on the 100 ns ADC raster, which only the dwell answers to.  A
        #     33.2 us offset produced one `RASTER` error per echo, 64 of them for 64 echoes.
        #   * the **gradient** raster, via `half`: the blip is centred on the junction, so it is
        #     quantised to an even number of gradient rasters and `half` is a multiple of one.
        #   * the **block** raster, via the ADC's node, which `build` places at the lobe start so
        #     it lands there by construction.
        #
        # Exact on the RF raster by the choice of `samples` above, so no rounding happens here and
        # the span stays centred.  `half` is a multiple of the RF raster too, being an even number of
        # gradient rasters, so the offset can never fall below it and let the blip into the sampling.
        offset = half + 0.5 * (window - samples * dwell)
        span = (
            _unit_moment(offset + samples * dwell, ramp, flat)
            - _unit_moment(offset, ramp, flat)
        )
        if span <= 0.0:
            return None
        return {
            'ramp': ramp,
            'flat': flat,
            'echo_spacing': echo_spacing,
            'offset': offset,
            'dwell': dwell,
            'samples': samples,
            'amplitude': target / span,
        }

    def _design_prephasers(self) -> None:
        """
        Design the readout and phase-encode prephasers, at the largest area either will need.

        Both keep one duration for every shot and polarity, so the train always starts at the same
        instant -- the same reason :class:`~seqcraft.modules.encoding.cartesian.PhaseEncode`
        designs at its largest blip.
        """
        biggest_line = max(abs(self.signed_lines(shot)[0]) for shot in range(self.n_shots))
        area_ro = -(self.partial_echo * self.k_max_ro_per_m + self._leading_moment)
        area_pe = max(biggest_line * self.dk_pe_per_m, self.dk_pe_per_m)
        reference = (
            pp.make_trapezoid(channel=self.axes[0], area=area_ro, system=self.opts),
            pp.make_trapezoid(channel=self.axes[1], area=area_pe, system=self.opts),
        )
        duration = self.system.grad_raster.ceil(float(pp.calc_duration(*reference)))
        self._pre_ro = pp.make_trapezoid(
            channel=self.axes[0], area=area_ro, duration=duration, system=self.opts
        )
        self._pre_pe = pp.make_trapezoid(
            channel=self.axes[1], area=area_pe, duration=duration, system=self.opts
        )

    # ------------------------------------------------------------------------- the waveforms
    def _readout_times(self) -> np.ndarray:
        """Knot times of the whole readout train, seconds from its start."""
        times = [0.0]
        for echo in range(self._n_echoes):
            start = echo * self._echo_spacing
            times.append(start + self._ramp)
            if self._flat > 0.0:
                times.append(start + self._ramp + self._flat)
            times.append(start + 2.0 * self._ramp + self._flat)
            if 2.0 * self._ramp + self._flat < self._echo_spacing:
                times.append(start + self._echo_spacing)
        return np.asarray(times)

    def _readout_amplitudes(self) -> np.ndarray:
        """Knot amplitudes of the readout train: a triangle wave of alternating sign."""
        amps = [0.0]
        for echo in range(self._n_echoes):
            peak = self._amplitude * (1.0 if echo % 2 == 0 else -1.0)
            amps.append(peak)
            if self._flat > 0.0:
                amps.append(peak)
            amps.append(0.0)
            if 2.0 * self._ramp + self._flat < self._echo_spacing:
                amps.append(0.0)
        return np.asarray(amps)

    def _blip_times(self) -> np.ndarray:
        """Knot times of the blip train: one triangular blip centred on each junction."""
        half = self._blip / 2.0
        times = [0.0]
        for junction in range(1, self._n_echoes):
            centre = junction * self._echo_spacing
            times += [centre - half, centre, centre + half]
        times.append(self._n_echoes * self._echo_spacing)
        return np.asarray(times)

    def _blip_amplitudes(self) -> np.ndarray:
        """Knot amplitudes of the blip train, at unit polarity."""
        peak = self._blip_area / (self._blip / 2.0)
        amps = [0.0]
        for _ in range(1, self._n_echoes):
            amps += [0.0, peak, 0.0]
        amps.append(0.0)
        return np.asarray(amps)

    # --------------------------------------------------------------------------- the k table
    @property
    def centre_line(self) -> int:
        """Reconstruction-grid index of the ``ky = 0`` line, ``matrix_pe // 2``."""
        return self.matrix_pe // 2

    def lines(self, shot: int = 0) -> tuple[int, ...]:
        """
        Reconstruction-grid line indices this shot acquires, in acquisition order.

        The index space is the one ``LIN`` and ``kSpaceCenterLine`` use, so a label value and this
        table cannot disagree.

        Examples
        --------
        >>> import seqcraft as sc
        >>> system = sc.System.preset('cima_x').derate('epi', grad=0.85, slew=1.0)
        >>> ro = EPIReadout(system, fov_ro_mm=240, matrix_ro=64, fov_pe_mm=240, matrix_pe=64,
        ...                 n_shots=2, regime='epi')
        >>> ro.lines(0)[:4], ro.lines(1)[:4]
        ((0, 2, 4, 6), (1, 3, 5, 7))
        >>> ro.lines(0)[ro.echo_index]                  # the centre line, for the shot that has it
        32
        """
        if not 0 <= shot < self.n_shots:
            msg = format_error(
                f'shot must be in [0, {self.n_shots - 1}].', {'got': shot},
            )
            raise ConfigurationError(msg)
        return tuple(range(self._first_line + shot, self.matrix_pe, self.n_shots))

    def signed_lines(self, shot: int = 0) -> tuple[int, ...]:
        """Return :meth:`lines` signed about the centre -- which is what the blips encode."""
        return tuple(index - self.centre_line for index in self.lines(shot))

    # --------------------------------------------------------------------------- properties
    @property
    def dk_pe_per_m(self) -> float:
        """Phase-encode step between adjacent **lines**, ``1 / FOV``, in 1/m."""
        return 1e3 / self.fov_pe_mm

    @property
    def k_max_ro_per_m(self) -> float:
        """Outer readout k-space radius, ``matrix / (2 FOV)``, in 1/m."""
        return self.matrix_ro / (2.0 * convert(self.fov_ro_mm, 'mm', 'm'))

    @property
    def resolution_mm(self) -> tuple[float, float]:
        """Nominal in-plane resolution ``(readout, phase)``, millimetres."""
        return (self.fov_ro_mm / self.matrix_ro, self.fov_pe_mm / self.matrix_pe)

    @property
    def n_echoes(self) -> int:
        """Echoes in one shot's train."""
        return self._n_echoes

    @property
    def echo_index(self) -> int:
        """
        Index of the echo :attr:`time_to_echo` points at -- **not** the middle of the train.

        With partial Fourier the ``ky = 0`` echo arrives early, which is the entire reason to use
        it.  Shot-independent; see the class notes.
        """
        return self._echo_index

    @property
    def echo_spacing(self) -> float:
        """Seconds between successive echoes.  Slew-limited on any realistic system."""
        return self._echo_spacing

    @property
    def ramp_time(self) -> float:
        """Seconds of ramp at each end of one readout lobe."""
        return self._ramp

    @property
    def flat_time(self) -> float:
        """Seconds of flat top on one readout lobe; zero for a minimum-time echo."""
        return self._flat

    @property
    def blip_duration(self) -> float:
        """Seconds occupied by one phase-encode blip."""
        return self._blip

    @property
    def blip_floors(self) -> dict[str, float]:
        """
        The three floors on the blip duration, seconds, so the binding one is visible.

        Examples
        --------
        >>> import seqcraft as sc
        >>> system = sc.System.preset('cima_x').derate('epi', grad=0.85, slew=1.0)
        >>> ro = EPIReadout(system, fov_ro_mm=240, matrix_ro=128, fov_pe_mm=240, matrix_pe=128,
        ...                 regime='epi')
        >>> max(ro.blip_floors, key=lambda k: ro.blip_floors[k])       # not the dead time
        'slew'
        """
        return dict(self._blip_floors)

    @property
    def peak_amplitude_Hz_per_m(self) -> float:
        """Peak readout gradient amplitude, Hz/m."""
        return self._amplitude

    @property
    def samples_per_echo(self) -> int:
        """ADC samples in one echo."""
        return self._samples

    @property
    def n_samples(self) -> int:
        """ADC samples in one whole shot."""
        return self._samples * self._n_echoes

    @property
    def dwell_s(self) -> float:
        """ADC dwell time, seconds."""
        return self._dwell

    @property
    def sample_offset(self) -> float:
        """
        Seconds from a lobe's start to its first sample -- the ADC event's own ``delay``.

        At least ``blip / 2``, so the blip falls entirely outside the sampling and ``ky`` is
        constant while a line is read.
        """
        return self._offset

    @property
    def bandwidth_per_pixel_hz(self) -> float:
        """Readout bandwidth per pixel, ``1 / (dwell * samples_per_echo)``."""
        return 1.0 / (self._dwell * self._samples)

    @property
    def pe_bandwidth_per_pixel_hz(self) -> float:
        """
        Phase-encode bandwidth per pixel, ``1 / (n_echoes * echo_spacing)``.

        **The distortion figure.**  Off-resonance displaces a voxel along the phase-encode
        direction by ``df / this``, in pixels -- so 20 Hz per pixel means 100 Hz of off-resonance
        moves it five pixels.  Nothing else in the sequence predicts that number.
        """
        return 1.0 / (self._n_echoes * self._echo_spacing)

    @property
    def total_readout_time_s(self) -> float:
        """
        ``(n_echoes - 1) * echo_spacing`` -- what FSL ``topup`` and ``eddy`` call readout time.

        Their convention is the time between the first and last echo, not the train's duration, and
        it is referenced to the *acquired* line count rather than the reconstructed matrix.
        """
        return (self._n_echoes - 1) * self._echo_spacing

    @property
    def train_duration(self) -> float:
        """Seconds occupied by the echo train alone, prephasers excluded."""
        return self._n_echoes * self._echo_spacing

    @property
    def prephase_duration(self) -> float:
        """Seconds occupied by the prephasers, or zero when this readout does not own them."""
        if not self.prephase:
            return 0.0
        return self.system.block_raster.ceil(float(pp.calc_duration(self._pre_ro, self._pre_pe)))

    @property
    def duration(self) -> float:
        """Seconds occupied by the whole built block."""
        return self.system.block_raster.ceil(self.prephase_duration + self.train_duration)

    @property
    def _leading_moment(self) -> float:
        """Area the readout lobe accumulates before its first sample, in 1/m."""
        return self._amplitude * _unit_moment(self._offset, self._ramp, self._flat)

    def _time_to_kx_zero(self, echo: int) -> float:
        """
        Seconds from a lobe's start to ``kx = 0`` within it.

        Which is not the middle of the sampling window unless `partial_echo` is 1, and not the same
        for even and odd echoes: an even echo enters at ``-partial_echo * k_max`` and an odd one at
        ``+k_max``, so they have different distances to travel.
        """
        entering = self.partial_echo * self.k_max_ro_per_m if echo % 2 == 0 else self.k_max_ro_per_m
        wanted = entering / self._amplitude + _unit_moment(self._offset, self._ramp, self._flat)
        return _unit_moment_inverse(wanted, self._ramp, self._flat)

    @property
    def time_to_echo(self) -> float:
        """
        Seconds from the start of the built block to **k = 0** -- what TE is measured from.

        Found from the moment, not assumed to be halfway through anything.  With partial Fourier
        0.75 on a 128 matrix it is 17.24 ms into a 49.92 ms train; taking the train's midpoint
        would put TE 7.7 ms late, silently.
        """
        return (
            self.prephase_duration
            + self._echo_index * self._echo_spacing
            + self._time_to_kx_zero(self._echo_index)
        )

    @property
    def duration_after_echo(self) -> float:
        """Seconds from k=0 to the end of the built block."""
        return self.duration - self.time_to_echo

    # -------------------------------------------------------------------------- the k-space
    def sample_times(self) -> np.ndarray:
        """
        ``(n_samples,)`` sample times **relative to the echo**, seconds.

        Negative before it, which is the whole difference from a spiral and the reason a
        reconstruction must be told where the echo is rather than inferring it from the first
        sample.  Identical for every shot, because the train is.

        Examples
        --------
        >>> import numpy as np
        >>> import seqcraft as sc
        >>> system = sc.System.preset('cima_x').derate('epi', grad=0.85, slew=1.0)
        >>> ro = EPIReadout(system, fov_ro_mm=240, matrix_ro=64, fov_pe_mm=240, matrix_pe=64,
        ...                 partial_fourier_pe=0.75, regime='epi')
        >>> t = ro.sample_times()
        >>> bool(t.min() < 0 < t.max())                  # the echo is inside the train
        True
        >>> bool(np.all(np.diff(t) > 0))                 # non-decreasing, as the operator needs
        True
        """
        within = self._offset + (np.arange(self._samples) + 0.5) * self._dwell
        starts = np.arange(self._n_echoes) * self._echo_spacing
        absolute = (starts[:, None] + within[None, :]).reshape(-1)
        echo_at = (
            self._echo_index * self._echo_spacing + self._time_to_kx_zero(self._echo_index)
        )
        return absolute - echo_at

    def trajectory(self, shot: int = 0, *, pe_polarity: int = 1) -> np.ndarray:
        """
        ``(2, n_samples)`` k-space coordinates at ADC sample times, in 1/m.

        Computed from the same moment integral the prephaser uses, so the trajectory handed to a
        reconstruction is the one the waveform produces.  Cross-check it against
        ``CompiledSequence.kspace()['k_adc']``, which the tests do -- an independent oracle is the
        only way to know a sidecar is not describing a different sequence.

        Examples
        --------
        >>> import numpy as np
        >>> import seqcraft as sc
        >>> system = sc.System.preset('cima_x').derate('epi', grad=0.85, slew=1.0)
        >>> ro = EPIReadout(system, fov_ro_mm=240, matrix_ro=64, fov_pe_mm=240, matrix_pe=64,
        ...                 regime='epi')
        >>> k = ro.trajectory(0)
        >>> k.shape == (2, ro.n_samples)
        True
        >>> round(float(np.abs(k[0]).max()), 1) <= round(ro.k_max_ro_per_m, 1)
        True

        Every echo is read at a constant ``ky``, because the blip sits outside the ADC window:

        >>> per_echo = k[1].reshape(ro.n_echoes, ro.samples_per_echo)
        >>> float(np.ptp(per_echo, axis=1).max())
        0.0
        """
        if pe_polarity not in (1, -1):
            msg = format_error(
                f'pe_polarity must be 1 or -1, got {pe_polarity!r}.', {},
                ['use -1 for the blip-down half of a distortion-correction pair'],
            )
            raise ConfigurationError(msg)
        within = self._offset + (np.arange(self._samples) + 0.5) * self._dwell
        travelled = np.array(
            [_unit_moment(float(t), self._ramp, self._flat) for t in within]
        ) - _unit_moment(self._offset, self._ramp, self._flat)

        kx = np.empty((self._n_echoes, self._samples))
        for echo in range(self._n_echoes):
            if echo % 2 == 0:
                kx[echo] = -self.partial_echo * self.k_max_ro_per_m + self._amplitude * travelled
            else:
                kx[echo] = self.k_max_ro_per_m - self._amplitude * travelled
        ky = (
            np.asarray(self.signed_lines(shot), dtype=float)
            * self.dk_pe_per_m
            * float(pe_polarity)
        )
        return np.stack([kx.reshape(-1), np.repeat(ky, self._samples)])

    def definitions(self) -> dict[str, float]:
        """Return the ``.seq`` definitions this readout is responsible for."""
        return {
            'BandwidthPerPixelHz': self.bandwidth_per_pixel_hz,
            'PhaseEncodeBandwidthPerPixelHz': self.pe_bandwidth_per_pixel_hz,
            'EchoSpacing': self.echo_spacing,
            'TotalReadoutTime': self.total_readout_time_s,
            'EPIFactor': self._n_echoes,
            'NumberOfShots': self.n_shots,
            'PartialFourierPE': self.partial_fourier_pe,
            'PartialEcho': self.partial_echo,
        }

    # ------------------------------------------------------------------------------- build
    def build(
        self,
        *,
        shot: int = 0,
        pe_polarity: int = 1,
        rf_phase_rad: float = 0.0,
    ) -> LogicBlock:
        """
        Return the prephasers, the two train gradients, one ADC per echo, and the labels.

        Parameters
        ----------
        shot
            Which interleaved segment.  Selects the phase-encode prephaser and the labels; the
            timing is identical for every shot.
        pe_polarity
            ``-1`` reverses the phase-encode direction, flipping both the blips and the PE
            prephaser -- the blip-down half of a distortion-correction pair.
        rf_phase_rad
            Receiver phase, which must be **the same value passed to the excitation**.  Invisible
            in a magnitude image and fatal in anything reading the phase.

        Examples
        --------
        >>> import seqcraft as sc
        >>> system = sc.System.preset('cima_x').derate('epi', grad=0.85, slew=1.0)
        >>> ro = EPIReadout(system, fov_ro_mm=240, matrix_ro=64, fov_pe_mm=240, matrix_pe=64,
        ...                 partial_fourier_pe=0.75, labels=False, regime='epi')
        >>> block = ro.build()
        >>> len(block.nodes)                    # the two prephasers, and the train
        3
        >>> ro.n_echoes
        48
        >>> down = ro.build(pe_polarity=-1)     # the same duration, mirrored in ky
        >>> abs(down.duration - block.duration) < 1e-12
        True
        """
        if pe_polarity not in (1, -1):
            msg = format_error(
                f'pe_polarity must be 1 or -1, got {pe_polarity!r}.', {},
                ['use -1 for the blip-down half of a distortion-correction pair'],
            )
            raise ConfigurationError(msg)
        lines = self.lines(shot)
        signed = self.signed_lines(shot)

        out = LogicBlock('epi')
        lead = 0.0
        if self.prephase:
            lead = self.prephase_duration
            out.add(
                0.0,
                self._pre_ro,
                _retarget(self._pre_pe, signed[0] * self.dk_pe_per_m * pe_polarity),
            )

        train = LogicBlock('epi_train')
        train.add(0.0, self._gx)
        train.add(0.0, self._gy if pe_polarity == 1 else _flip(self._gy))
        adc = (
            self.adc if rf_phase_rad == 0.0
            else ev.derive(self.adc, phase_offset=float(rf_phase_rad))
        )
        for echo in range(self._n_echoes):
            train.add(echo * self._echo_spacing, adc)
        if self.labels:
            # A label is an event, so it is a node.  The compiler attaches each to the block holding
            # the first ADC at or after its own time, so a boundary landing early cannot make a LIN
            # address the previous echo.
            train.add(0.0, pp.make_label('SEG', 'SET', shot))
            for echo in range(self._n_echoes):
                at = echo * self._echo_spacing
                train.add(at, pp.make_label('LIN', 'SET', int(lines[echo])))
                train.add(at, pp.make_label('REV', 'SET', echo % 2 == 1))
        out.add(lead, train)
        return out

    def prephaser_block(self, *, shot: int = 0, pe_polarity: int = 1) -> LogicBlock:
        """
        Return just the prephasers, for placing them separately.

        Which is the point of ``prephase=False``: the readout prephaser on x, the phase-encode
        prephaser on y and a slice rewinder on z coincide, and the compiler puts all three in one
        block because they are on different axes.

        Examples
        --------
        >>> import seqcraft as sc
        >>> system = sc.System.preset('cima_x').derate('epi', grad=0.85, slew=1.0)
        >>> ro = EPIReadout(system, fov_ro_mm=240, matrix_ro=64, fov_pe_mm=240, matrix_pe=64,
        ...                 prephase=False, regime='epi')
        >>> ro.prephaser_block()
        LogicBlock(epi_pre, 2 nodes, ... ms)
        """
        signed = self.signed_lines(shot)
        return LogicBlock('epi_pre').add(
            0.0,
            self._pre_ro,
            _retarget(self._pre_pe, signed[0] * self.dk_pe_per_m * pe_polarity),
        )


def _flip(grad: SimpleNamespace) -> SimpleNamespace:
    """Return an arbitrary or extended gradient with its sign reversed."""
    return ev.derive(
        grad,
        waveform=-np.asarray(grad.waveform),
        first=-float(getattr(grad, 'first', 0.0)),
        last=-float(getattr(grad, 'last', 0.0)),
    )


def _retarget(reference: SimpleNamespace, area: float) -> SimpleNamespace:
    """
    Scale a trapezoid to a new area, keeping its duration.

    Keeping the duration is what makes every shot's train start at the same instant, so
    :attr:`EPIReadout.time_to_echo` holds whichever shot is built.
    """
    factor = 0.0 if reference.area == 0 else area / float(reference.area)
    return ev.derive(
        reference,
        amplitude=float(reference.amplitude) * factor,
        area=float(reference.area) * factor,
        flat_area=float(getattr(reference, 'flat_area', 0.0)) * factor,
    )
