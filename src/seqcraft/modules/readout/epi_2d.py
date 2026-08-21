"""
:class:`EPI2D` -- one excitation's worth of k-space, as a single blipped echo-planar train.

``readout/`` because it contains ADCs -- one per echo -- and because it is a *leaf*: it holds no
RF and knows nothing about what excited the magnetisation it is sampling.  A gradient echo, a spin
echo and an inversion-prepared train differ only in what is played before it, which is exactly the
split that lets ``examples/gre_epi_2d/`` and ``examples/se_epi_2d/`` share one readout.

What this module is, against :class:`~seqcraft.modules.CartesianLine`
--------------------------------------------------------------------
A Cartesian line is *one* line and does not know which; the caller supplies a phase encoding and
repeats it.  An EPI train is the opposite: **the whole of k-space in one shot**, so the ordering
is inside the readout rather than around it, and three things follow that a per-line readout never
has to think about.

**The lobes alternate sign**, so every other line is sampled from ``+k`` to ``-k`` and has to be
reversed before it can be transformed.  That is not bookkeeping -- it is the reason an EPI ghosts:
anything that shifts the two polarities differently in ``k`` (a gradient delay, an eddy current, a
constant off-resonance) becomes a modulation of period two along ``ky``, and a modulation of
period two is a **replica of the object at half the field of view**.  So the *sequence itself* has
to be exactly symmetric, or the ghost being measured is the sequence's rather than the scanner's;
:attr:`k_read_per_m` and the guard below are what buy that.

**There is no prephaser per line.**  ``ky`` is stepped by a **blip** between echoes, so the
encoding is a train of *differences* rather than a train of positions -- which is why
:class:`~seqcraft.modules.PhaseEncode`, whose whole vocabulary is "line `n` of `matrix`", designs
the blip-axis *prephaser* here and not the blips.  A blip is the same trapezoid with a different
meaning, and the difference is which of the two the ordering table indexes.

**Time runs through k-space.**  Every echo is a different point on whatever decay and phase
evolution the magnetisation is undergoing, so the ordering is a filter across ``ky`` in exactly
the sense ``examples/fse_2d/`` measures for a turbo spin echo -- and off-resonance, which a spin
warp acquisition converts into a per-voxel phase and nothing else, here becomes geometric
distortion along ``ky`` and a ghost along it.  That is the artefact, and it is a property of the
*ordering*, not of the hardware.

**Segmentation makes that a statement about shots as well as echoes.**  An interleaved table puts
consecutive ``ky`` lines in *different* shots, so anything that differs between shots is a
modulation of ``ky`` with the shot plan's period -- a replica of the object at ``Ny/shots``, which
reads as an under-sampling artefact and has nothing to do with sampling.  That is why
:meth:`build` takes ``phase_deg``: a spoiled acquisition needs the receiver phase-locked to a
transmitter that is advancing its carrier, and the flip angle and the run-in that put every shot
at the same point of the steady state belong to whatever composite owns the schedule.
``examples/gre_epi_2d/02`` section 8 measures all three.

The one piece of arithmetic worth reading twice
-----------------------------------------------
**The sampling window has to be exactly centred in its lobe.**  Write ``T`` for the lobe duration,
``g`` for the un-sampled guard at each end, ``N`` for the sample count and ``tau`` for the dwell::

    T = 2*g + N*tau                     # exactly, not "at least"

The reason is the reversal.  A trapezoid is symmetric about its own midpoint, so its running area
satisfies ``A(t) + A(T - t) = A_total``; sample `i` of a forward lobe sits at ``t_i = g + (i+1/2)*tau``
and sample ``N-1-i`` of the reverse lobe that follows it sits at ``T - t_i`` **if and only if**
``2*g + N*tau`` is exactly ``T``.  When it is, the two polarities sample the *same* set of
``k`` positions in opposite order, to floating-point; when it is not -- when the flat top was
rounded up onto the gradient raster and the slack left at one end, which is the natural way to
write it -- they sample two grids a fraction of a sample apart, and that fraction is an N/2 ghost
the sequence made itself.

Everything else here follows from that one equation:

- ``g`` is not free.  It must hold **half a blip** at each end, because the blip is centred on the
  seam between two lobes, and it must hold the receiver's own **dead time**, which is also what
  guarantees the compiler can cut a block boundary between consecutive ADCs.  So
  ``g = max(blip_duration/2, adc_dead_time)``; see :meth:`EPI2D._seam` for why the blip goes on
  the zero crossing and why the barrier that follows it is not optional.
- With ``T`` fixed by ``g`` and the dwell, the **amplitude is what is left to solve**.  With
  ``ramp_sampling`` the window reaches into both ramps, so the peak gradient is *higher* than the
  nominal ``dk/tau`` and the average over the window is what equals it.  ``G*T - G^2/S - S*g^2 =
  N*dk`` is that statement in closed form, and the smaller root is the answer.
- The solved amplitude is then **measured rather than trusted**: the lobe is built, the window
  area is integrated from the event's own knots with
  :func:`~seqcraft.modules._support.area_until`, and the amplitude is rescaled by the ratio.  Area
  is exactly linear in amplitude at fixed ramp times, so one rescale is exact -- and it absorbs
  every raster rounding in one place rather than leaving each of them in the trajectory.

``salvage/epi_moment.py`` holds the unit-amplitude version of that integral and a closed-form
inverse, lifted out of a deleted EPI module for this one to start from.  The integral is here as
:func:`~seqcraft.modules._support.area_until`, which is more general -- it works on any pulseq
gradient and carries the event's own delay.  The *inverse* turned out not to be needed at all: the
echo of an EPI lobe is a **sample**, not a time, so where ``k`` crosses zero is decided by the
sample grid and read off with the same integral.  Solving for the crossing and then finding the
nearest sample is how the two come to disagree by half a dwell.

Ramp sampling is not a bandwidth trick
--------------------------------------
``ramp_sampling=False`` puts the window on the flat top, where ``k`` advances uniformly and no
reconstruction step is needed.  ``ramp_sampling=True`` opens it as wide as the guard allows, which
shortens the lobe by ``2*(rise - g)`` -- on the reference protocol here, 630 us becomes 510 us, a
19 % shorter echo spacing and therefore 19 % less distortion and 19 % less T2* decay across the
train.  What it costs is that the samples are **not uniformly spaced in k**: the ones on the ramps
are closer together.  :attr:`k_read_per_m` is where they are, and a reconstruction that ignores it
gets an image stretched at the edges of the readout direction and a point spread with tails.

That array is the module's real output alongside the ADC, which is why it is a public attribute
rather than something the notebook re-derives.

And it is why `oversampling` defaults to 2 rather than 1.  Twice the samples at half the dwell is
the *same* sampling duration, so the lobe, the guard, the peak gradient and the bandwidth per
reconstructed pixel are all unchanged -- oversampling costs no echo spacing at all.  What it buys
is the sample *spacing*: a ramp-sampled lobe at ``oversampling=1`` leaves gaps of 1.12 times the
reconstruction ``dk`` across its flat top, so the readout does not sample the grid it is
reconstructed onto and no interpolator can recover what was never measured.  At 2 the worst gap is
0.56 and the regridding is a solved problem.
"""

from __future__ import annotations

import math
from math import gcd
from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ...design.events import derive
from ...design.logic import LogicBlock, barrier
from ...design.module import Module
from ...design.timing import EPS, from_ticks, to_ticks
from ...errors import ConfigurationError, format_error
from .._support import (
    area_until,
    ceil_raster,
    require_axis,
    require_pair,
    require_positive,
    require_range,
)
from ..encoding.phase_encoding import PhaseEncode

if TYPE_CHECKING:
    from collections.abc import Iterable

    from pypulseq.opts import Opts

    from ...design.events import Event

__all__ = ['EPI2D']

#: How many gradient rasters the ramp search may add before it gives up.  A bound rather than a
#: ``while True``: the search converges when a longer ramp costs less amplitude than it buys, and
#: refusing with a bandwidth to try beats spinning on a protocol that cannot be built.
_MAX_RAMP_STEPS = 256


class EPI2D(Module):
    """
    A blipped 2D echo-planar readout: prephasers, alternating lobes, blips, one ADC per echo.

    Parameters
    ----------
    opts
        The scanner.
    fov_mm
        ``(readout, blip)`` in millimetres, or a scalar meaning square.  Sets both k-space steps.
    matrix
        ``(nx, ny)`` -- readout then blip axis, matching image-array order.
    bandwidth_hz_px, dwell_s
        The sampling rate, either way round.  **Exactly one is required**; both or neither raises.
        `bandwidth_hz_px` is per *reconstructed* pixel, so `oversampling` does not move it.  With
        `ramp_sampling` it is the *nominal* rate -- the one the flat top would need -- and the
        peak gradient comes out higher; see :attr:`readout_amplitude_hz_m`.  The dwell is snapped
        **up** onto a quantum coarser than the ADC raster; :attr:`bandwidth_hz_px` reads back what
        was achieved rather than what was asked for.
    partial_fourier
        Readout-direction partial echo: the fraction of the full ``kx`` extent that is sampled,
        with the missing part taken off the **pre-echo** side.  ``1.0`` is a full echo.  The
        sample count is rounded to an **even** number, because the forward/reverse symmetry above
        needs one.
    oversampling
        Samples per reconstructed pixel along the readout.  It multiplies :attr:`num_samples` and
        divides :attr:`dwell_s`, so the sampling duration -- and therefore the lobe, the guard,
        the peak gradient and :attr:`bandwidth_hz_px` -- does not move.  ``2`` is the default
        because ``1`` leaves the ramp-sampled flat top under-sampled against the grid it is
        reconstructed onto; see the module docstring.
    ramp_sampling
        Open the ADC during the ramps.  ``True`` is the shorter echo spacing and a non-uniform
        :attr:`k_read_per_m`; ``False`` samples the flat top only and leaves ``k`` uniform.
    blip_lines
        The largest ``|delta ky|``, **in lines**, that any shot this instance builds will step
        between consecutive echoes.  It sets the blip's duration and therefore the guard, the lobe
        and the echo spacing -- which is why it is a design-time number rather than something read
        off the table: an echo spacing that changed with the ordering would make every timing
        query below depend on the shot.  ``1`` is single-shot and blocked multi-shot; interleaved
        multi-shot steps by the number of shots.
    blip_duration_s
        Lengthen the blip beyond its own minimum.  ``None`` is the minimum.  Lengthening it
        lengthens the guard, and so the echo spacing, by twice as much -- it is the knob for
        trading echo spacing against blip-axis slew, which is what peripheral nerve stimulation
        responds to first in an EPI.
    navigator_echoes
        Readout lobes played **before** the blip-axis prephaser, with no blips, an ADC each and
        ``NAV = 1``.  Three is the usual number.  They sample the same ``kx`` grid as the imaging
        echoes at ``ky = 0``, which is what makes them a measurement of the odd/even phase
        difference and nothing else.  ``0`` builds no navigator and emits no ``NAV`` label.
    prephase
        ``False`` drops both prephasers, leaving the train starting wherever the caller left
        ``k``.  That is the **spin-echo** placement: the dephasers go *before* the refocusing
        pulse with the opposite sign, where they cost no echo time because the interval is
        already there.  See :attr:`area_to_echo_per_m` and :meth:`k_blip_per_m`, which are what a
        caller then builds them from.  With navigators it also means they sample at the caller's
        ``ky`` rather than at zero.
    prephaser_duration_s
        Stretch both prephaser windows to this length.  ``None`` is the coupled minimum -- the
        longer of the two axes' own minima, since with no navigators they play at the same
        instant.
    axis, blip_axis
        Logical gradient channels for the readout and the blips.  Must differ.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    pe : PhaseEncode
        The blip-axis prephaser, held as a plain attribute.  It is a *position* on ``ky``, which
        is what that module is for; the blips are steps and are designed here.
    gx : SimpleNamespace
        The positive readout lobe.  Every other echo plays its negation.
    adc : SimpleNamespace
        The ADC, identical at every echo -- same dwell, same delay, same sample count.
    k_read_per_m : numpy.ndarray
        Where the samples of a **forward** lobe land along the readout axis, in 1/m, relative to
        ``k = 0`` at :attr:`pre_echo_samples`.  A reverse lobe samples the same positions in the
        opposite order, so ``k_read_per_m[::-1]`` is its trajectory.

    Examples
    --------
    >>> import numpy as np
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=180, slew_unit='T/m/s',
    ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> epi = EPI2D(opts=o, fov_mm=220.0, matrix=(64, 64), dwell_s=3.5e-6)
    >>> epi.num_samples, epi.pre_echo_samples, round(epi.echo_spacing_s * 1e6)
    (128, 64, 510)

    The sampling window is exactly centred in its lobe, which is what makes the reversal exact:

    >>> abs(2 * epi.guard_s + epi.num_samples * epi.dwell_s - epi.echo_spacing_s) < 1e-12
    True
    >>> float(epi.k_read_per_m[epi.echo_sample(0)])           # k = 0 is a sample, not a moment
    0.0
    >>> epi.polarity(0), epi.polarity(1), epi.echo_sample(1)
    (1, -1, 63)

    Ramp sampling shortens the lobe, raises the peak gradient above the nominal, and leaves the
    extent alone -- it is the *spacing* that stops being uniform:

    >>> flat = EPI2D(opts=o, fov_mm=220.0, matrix=(64, 64), dwell_s=3.5e-6, ramp_sampling=False)
    >>> round(flat.echo_spacing_s * 1e6), round(epi.echo_spacing_s * 1e6)
    (630, 510)
    >>> [round(m.readout_amplitude_hz_m / o.gamma * 1e3, 2) for m in (flat, epi)]
    [15.25, 17.07]
    >>> np.round([np.ptp(np.diff(flat.k_read_per_m)), np.ptp(np.diff(epi.k_read_per_m))], 3)
    array([0.   , 1.666])

    Oversampling is free: twice the samples at half the dwell is the same sampling duration, so
    the lobe does not move and only the sample *spacing* does -- which is the whole point, since
    at ``oversampling=1`` the widest gap is above one reconstruction ``dk``:

    >>> one = EPI2D(opts=o, fov_mm=220.0, matrix=(64, 64), bandwidth_hz_px=epi.bandwidth_hz_px,
    ...             oversampling=1)
    >>> (one.echo_spacing_s, one.guard_s) == (epi.echo_spacing_s, epi.guard_s)
    True
    >>> [round(float(np.diff(m.k_read_per_m).max() / m.dk_read_per_m), 3) for m in (one, epi)]
    [1.119, 0.559]

    One instance builds any ordering, because the table is a build argument -- and an ordering
    that steps further needs a longer blip, which is what `blip_lines` prices at design time:

    >>> single = epi(lines=range(64))
    >>> len(single), round(single.duration * 1e3, 2)
    (384, 32.94)
    >>> four = EPI2D(opts=o, fov_mm=220.0, matrix=(64, 64), dwell_s=3.5e-6, blip_lines=4)
    >>> round(four.echo_spacing_s * 1e6), [len(four(lines=range(s, 64, 4))) for s in range(4)]
    (550, [96, 96, 96, 96])

    An echo time is measured to the centre of k-space, so it is a property of the table:

    >>> round(epi.time_to_center_line(range(64)) * 1e3, 3)
    16.877
    >>> round(epi.time_to_center_line(list(range(64))[::-1]) * 1e3, 3)
    16.363

    Notes
    -----
    **`lines` is one shot, in acquisition order**, and nothing here knows how many shots there
    are.  That is the same choice :class:`~seqcraft.modules.GRE2D` makes for its own table and for
    the same reason: which lines to acquire, and in which order, is a sequence-programming
    decision, and a module that generated it would be handing back the thing taking a list was
    meant to avoid.  ``segment=`` is the shot's index if the caller wants it in the file.

    **The blip is designed once at `blip_lines` and scaled per step**, exactly as
    :class:`~seqcraft.modules.PhaseEncode` designs one blip and scales it: every step then takes
    the same time, so the echo spacing is a constant of the instance rather than a property of
    the table.  A step larger than `blip_lines` is refused rather than lengthened, because
    lengthening it silently would move every echo after it.

    **The polarity of the first imaging echo depends on `navigator_echoes`.**  Lobes alternate
    across the *whole* train, navigators included, because ``k`` on the readout axis has to arrive
    at the imaging train from wherever the navigators left it.  Three navigators therefore start
    the imaging train on a reverse lobe, which is what the ``REV`` label records so no
    reconstruction has to work it out.
    """

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: float | tuple[float, float],
        matrix: tuple[int, int],
        bandwidth_hz_px: float | None = None,
        dwell_s: float | None = None,
        partial_fourier: float = 1.0,
        oversampling: int = 2,
        ramp_sampling: bool = True,
        blip_lines: int = 1,
        blip_duration_s: float | None = None,
        navigator_echoes: int = 0,
        prephase: bool = True,
        prephaser_duration_s: float | None = None,
        axis: str = 'x',
        blip_axis: str = 'y',
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        fov_read, fov_blip = require_pair(fov_mm, 'fov_mm')
        nx, ny = require_pair(matrix, 'matrix')
        self.fov_mm = (require_positive(fov_read, 'fov_mm'), require_positive(fov_blip, 'fov_mm'))
        self.matrix = (int(nx), int(ny))
        self.axis = require_axis(axis)
        self.blip_axis = require_axis(blip_axis, 'blip_axis')
        self.ramp_sampling = bool(ramp_sampling)
        self.prephase = bool(prephase)
        self.partial_fourier = require_range(partial_fourier, 'partial_fourier', low=0.0, high=1.0)
        self.oversampling = _require_count(oversampling, 'oversampling', hint='1 turns it off')
        self.blip_lines = _require_count(
            blip_lines, 'blip_lines', hint='1 is single-shot and blocked multi-shot',
        )
        self.navigator_echoes = _require_count(
            navigator_echoes, 'navigator_echoes', low=0, hint='3 is the usual number',
        )
        self._check_axes()

        #: k-space step between adjacent *reconstructed* samples along `axis`, 1/m.  With
        #: `oversampling` the ADC samples are that much closer together than this.
        self.dk_read_per_m = 1e3 / self.fov_mm[0]
        #: k-space step between adjacent lines along `blip_axis`, 1/m.
        self.dk_blip_per_m = 1e3 / self.fov_mm[1]
        # The sample count comes *before* the dwell, because it is what the dwell is quantised
        # against -- see `_resolve_dwell`.
        self.num_samples, self.pre_echo_samples = self._resolve_samples()
        self.dwell_s = self._resolve_dwell(bandwidth_hz_px, dwell_s)
        #: Nominal receive bandwidth per *reconstructed* pixel, Hz.  Reads back the snapped dwell.
        self.bandwidth_hz_px = 1.0 / (self.dwell_s * self.matrix[0] * self.oversampling)

        # The blip next: it sets the guard, the guard sets the lobe, and the lobe sets everything
        # that is measured from it.  Designed at the largest step and scaled per call.
        self._blip = self._design_blip(blip_duration_s)
        self._guard_min_s = max(self.blip_duration_s / 2, float(self.opts.adc_dead_time))

        self.gx, self.adc, self.guard_s = self._design_lobe()
        self._gx_reverse = derive(pp.scale_grad(self.gx, -1.0))
        #: The lobe duration, which is the echo spacing: ``2*guard_s + num_samples*dwell_s``.
        self.echo_spacing_s = float(pp.calc_duration(self.gx))
        self._echo_in_lobe_s = float(self.adc.delay) + (self.pre_echo_samples + 0.5) * self.dwell_s
        self.k_read_per_m = self._sampled_k()

        self.pe, self.winder_s, self._prephaser = self._design_prephasers(prephaser_duration_s)

    # ------------------------------------------------------------------ what it knows
    @property
    def center_line(self) -> int:
        """The blip-axis line that encodes ``k = 0``.  One convention, defined by `PhaseEncode`."""
        return self.pe.center_line

    @property
    def blip_duration_s(self) -> float:
        """Seconds one blip occupies, whatever its step.  Every step takes the same time."""
        return float(pp.calc_duration(self._blip))

    @property
    def blip_area_per_m(self) -> float:
        """The designed blip's area, 1/m -- ``blip_lines * dk_blip_per_m``, as achieved."""
        return float(self._blip.area)

    @property
    def readout_amplitude_hz_m(self) -> float:
        """
        Peak amplitude of a readout lobe, Hz/m.

        Above the nominal ``dk_read_per_m / dwell_s`` whenever `ramp_sampling` is on, because the
        window then averages the ramps in: it is the *mean* over the sampling window that equals
        the nominal amplitude, and the peak is what the amplifier has to deliver.
        """
        return float(self.gx.amplitude)

    @property
    def area_to_echo_per_m(self) -> float:
        """
        What the **first readout lobe** accumulates by its ``k = 0`` sample, 1/m.

        The physics number rather than an event's, in the same sense as
        :attr:`~seqcraft.modules.CartesianLine.area_to_echo_per_m`: it is what the readout-axis
        prephaser cancels when there is one, and it is what a caller placing the dephaser before
        a refocusing pulse has to supply with the opposite sign.  Always available, including
        with ``prephase=False``.
        """
        return area_until(self.gx, self._echo_in_lobe_s)

    @property
    def prephaser_duration_s(self) -> float:
        """
        Seconds each prephaser window occupies.

        One number for both windows even when navigators separate them, so that
        `prephaser_duration_s` means one thing.  Raises with ``prephase=False``, where no window
        exists.
        """
        return self._require_prephase('prephaser_duration_s')

    def k_blip_per_m(self, line: int) -> float:
        """Return the blip-axis k the train must start at to acquire `line` first, in 1/m."""
        return self.pe.k_per_m(line)

    def polarity(self, echo: int) -> int:
        """
        Return ``+1`` or ``-1`` for echo `echo` -- which direction it sweeps ``k``.

        Counted across the whole train, so `navigator_echoes` moves it.  A ``-1`` echo's samples
        run from ``+k`` to ``-k``, so a reconstruction reverses them; the ``REV`` label carries
        the same answer in the written file.

        **`echo` is an offset from the first imaging echo, so a navigator is negative.**  The
        lobes are one alternating sequence and the navigators are the ones before the imaging
        train, so ``-1`` is the last navigator and ``-navigator_echoes`` the first -- which is
        what lets a reconstruction write ``range(-navigator_echoes, len(lines))`` and get the
        whole file's readouts in order.  Anything earlier than that raises.
        """
        return 1 if (self.navigator_echoes + self._check_echo(echo)) % 2 == 0 else -1

    def echo_sample(self, echo: int) -> int:
        """
        Return **which ADC sample** of `echo` carries ``k = 0`` on the readout axis.

        :attr:`pre_echo_samples` on a forward lobe and ``num_samples - 1 - pre_echo_samples`` on a
        reverse one, mirrored about the window's centre because the lobe is.  A method rather than
        an expression each caller writes because it has five call sites -- :meth:`time_to_echo`,
        the tests, and every reconstruction -- and its wrong branch is an N/2 ghost.

        Negative indices address navigators, as :meth:`polarity` describes.
        """
        return (self.pre_echo_samples if self.polarity(echo) > 0
                else self.num_samples - 1 - self.pre_echo_samples)

    def time_to_echo(self, echo: int = 0) -> float:
        """
        Seconds from the start of this module's block to the ``k = 0`` sample of `echo`.

        Not measurable from the tree, which is the usual reason a module answers a time: the block
        knows when its ADCs open, not which sample inside one is the echo.  Which sample that is
        is :meth:`echo_sample`, so this is one expression rather than two cases.

        Negative indices address navigators, as :meth:`polarity` describes -- and they are not
        simply ``echo * echo_spacing_s`` before imaging echo zero, because the blip-axis prephaser
        plays *between* the two.  That offset is the reason this is a method and not a formula.
        """
        index = self._check_echo(echo)
        start_s = (self._train_start_s + index * self.echo_spacing_s if index >= 0 else
                   self._nav_start_s + (index + self.navigator_echoes) * self.echo_spacing_s)
        return start_s + float(self.adc.delay) + (self.echo_sample(index) + 0.5) * self.dwell_s

    def time_to_center_line(self, lines: Iterable[int]) -> float:
        """
        Seconds from the start of this block to the acquisition of ``k = 0`` in both directions.

        Parameters
        ----------
        lines
            **The same table** :meth:`build` takes, because the answer moves with it: the centre
            of k-space is acquired at whichever echo the table puts it at, and that is the whole
            difference between a linear ordering and a centric one.

        Raises
        ------
        ConfigurationError
            On any table :meth:`build` refuses, and on one that omits :attr:`center_line` -- a
            shot that does not acquire the centre has no echo time to report, and returning one
            would be worse than refusing.

        Notes
        -----
        **This is what an echo time means in an EPI**, and it is not the first echo.  A single
        shot with a linear ordering reaches the centre halfway through the train, so TE is
        roughly half the readout; a centric ordering reaches it at the first echo, and the same
        sequence then has a TE shorter by tens of milliseconds.  Placing a refocusing pulse or an
        inversion against the first echo instead is an error of exactly that size.
        """
        table = self._check(lines)
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
        return self.time_to_echo(table.index(self.center_line))

    # ----------------------------------------------------------------------- assembly
    def build(
        self,
        *,
        lines: Iterable[int],
        acquire: bool = True,
        segment: int | None = None,
        reference: bool | None = None,
        phase_deg: float = 0.0,
    ) -> LogicBlock:
        """
        Return one shot: prephasers, any navigators, then one lobe and one ADC per entry of `lines`.

        Parameters
        ----------
        lines
            Blip-axis indices, **in acquisition order**, for this shot.  Any iterable of ints.
            Consecutive entries may step by at most `blip_lines`.
        acquire
            ``False`` plays identical gradients with no ADCs and no labels, for a dummy shot.  The
            duration is unchanged, which is the point: a dummy has to load the gradients exactly
            as a real shot does, or the eddy-current state it is establishing is not the one that
            will be acquired.
        segment
            The shot's index, emitted as a ``SEG`` label on its first acquired readout.  ``None``
            emits none -- a single-shot acquisition has nothing to segment.
        reference
            Whether this shot is autocalibration data rather than imaging data, emitted as the
            complementary ``REF``/``IMA`` pair on its first acquired readout.  ``None`` emits
            neither.  It is the **only** thing parallel imaging adds here: the acceleration is
            `blip_lines` and the table is `lines`, and both already exist for their own reasons.
            A reconstruction has to tell a calibration shot from an imaging one, and counting
            shots is how it comes to guess wrong.
        phase_deg
            The receiver phase, and it is not optional in a spoiled sequence.  The receiver is
            phase-locked to the transmitter: whatever carrier phase the excitation was given,
            every ADC of the train it opens has to demodulate at.  One train is one excitation,
            so one number covers every echo in it.
        """
        table = self._check(lines)
        adc = self._adc_for(float(phase_deg))
        out = LogicBlock()
        if self.prephase:
            out.add(0.0, self._prephaser)
            # With navigators the blip-axis prephaser waits until they are finished, so that they
            # sample ky = 0 rather than the edge of k-space.  A navigator at the edge measures the
            # same odd/even phase difference multiplied by whatever signal is there, which on a
            # real object is close to none.
            out.add(self._train_start_s - self.winder_s, self.pe(line=table[0]))

        for n in range(self.navigator_echoes):
            t0 = self._nav_start_s + n * self.echo_spacing_s
            out.add(t0, self._lobe(n))
            if acquire:
                out.add(t0, adc)
                if n == 0:
                    out.add(t0, pp.make_label(type='SET', label='NAV', value=1))
                out.add(t0, pp.make_label(type='SET', label='REV', value=int(n % 2)))
            self._seam(out, t0 + self.echo_spacing_s, 0)

        for n, line in enumerate(table):
            t0 = self._train_start_s + n * self.echo_spacing_s
            index = self.navigator_echoes + n
            out.add(t0, self._lobe(index))
            if acquire:
                out.add(t0, adc)
                out.add(t0, pp.make_label(type='SET', label='LIN', value=int(line)))
                out.add(t0, pp.make_label(type='SET', label='REV', value=int(index % 2)))
                if n == 0:
                    self._shot_labels(out, t0, segment, reference)
            if n + 1 < len(table):
                self._seam(out, t0 + self.echo_spacing_s, table[n + 1] - line)
        return out

    def _adc_for(self, phase_deg: float) -> Event:
        """Return the train's ADC demodulating at the carrier phase its excitation was given."""
        if phase_deg == 0.0:
            return self.adc
        return derive(self.adc, phase_offset=float(np.deg2rad(phase_deg)))

    def _shot_labels(
        self, out: LogicBlock, t0: float, segment: int | None, reference: bool | None,
    ) -> None:
        """Emit what describes the *shot* rather than the echo, on its first acquired readout."""
        if self.navigator_echoes:
            out.add(t0, pp.make_label(type='SET', label='NAV', value=0))
        if segment is not None:
            out.add(t0, pp.make_label(type='SET', label='SEG', value=int(segment)))
        if reference is not None:
            out.add(t0, pp.make_label(type='SET', label='REF', value=int(bool(reference))))
            out.add(t0, pp.make_label(type='SET', label='IMA', value=int(not reference)))

    def _seam(self, out: LogicBlock, seam_s: float, step: int) -> None:
        """
        Put the blip on the readout's zero crossing, and force the block boundary there.

        Two statements, and the second is not optional -- which is why they are one method rather
        than two lines in :meth:`build` that a later edit could separate.

        The blip is **centred** on the seam, where the readout gradient passes through zero.  That
        is the instant at which the two axes ask least of the amplifier together, and it is the
        only placement for which the guard need be half a blip rather than a whole one: 510 us of
        echo spacing against 570 at the reference protocol.

        It therefore *straddles* the block boundary the compiler has to cut between two ADCs, and
        there is no instant in that gap which cuts nothing -- every instant strictly inside it is
        in the fall ramp, and the seam itself is inside the blip.  Left to choose, the compiler
        takes the midpoint of the gap and **splits every readout lobe in the train** into two
        arbitrary gradients.  That compiles, passes every k-space check and simulates correctly;
        the only report of it is a ``merge`` warning naming the *readout* axis.

        :func:`~seqcraft.barrier` is the documented escape hatch for exactly this -- *a gradient
        you want split at a known instant so a later reconstruction step can find the seam*.  With
        it the boundary **is** the seam, the lobes stay whole, and the blip splits into the two
        halves ``mr.splitGradientAt`` produces by hand in pulseq's own ``writeEpiRS.m``.  The one
        consequence to expect, not to fix, is a ``merge`` warning on the **blip** axis, once per
        compile: that is the two halves being summed back together, reported instead of hidden.
        """
        if step:
            out.add(seam_s - self.blip_duration_s / 2, self._blip_for(step))
        out.add(seam_s, barrier('seam'))

    # -------------------------------------------------------------------------- design
    def _resolve_dwell(self, bandwidth_hz_px: float | None, dwell_s: float | None) -> float:
        """
        Return the dwell, snapped so that the **guard lands on the RF raster**.

        Not the ADC raster, which is the obvious answer and the wrong one.  The guard is the ADC's
        ``delay``, and pypulseq's timing check divides an RF, ADC or output event's delay by
        ``rf_raster_time`` -- 1 us on Siemens, ten times coarser than the ADC raster the *dwell*
        lives on.  So Rule 1 tightens: ``T - N*dwell`` must be an even number of RF rasters, and
        since ``T`` is on the gradient raster the whole burden falls on ``N*dwell``.

        Solved in integer ticks rather than searched.  With ``a`` the ADC raster and
        ``u = 2*rf_raster``, a dwell of ``k*a`` works exactly when ``N*k*a`` is a multiple of
        ``u``, i.e. when ``k`` is a multiple of ``u / gcd(N*a, u)`` -- at 128 samples, a multiple
        of 500 ns.  Coarse, visible in the protocol table, and much better found here than in a
        timing-check failure that names a block and not a femtosecond.

        Rounded **up**, not to nearest: rounding down raises the bandwidth, which shortens the
        lobe, which is the direction that turns a feasible protocol into a refusal.  A longer
        dwell always builds; the achieved rate is read back off :attr:`bandwidth_hz_px`.
        """
        if (bandwidth_hz_px is None) == (dwell_s is None):
            msg = format_error(
                'exactly one of bandwidth_hz_px and dwell_s is required.',
                {'bandwidth_hz_px': bandwidth_hz_px, 'dwell_s': dwell_s},
                [
                    'bandwidth_hz_px=2841 is how a protocol quotes it',
                    'dwell_s=5.5e-6 is how the receiver is actually programmed',
                    'they are one number: dwell = 1 / (bandwidth_hz_px * matrix[0])',
                    'with ramp_sampling it is the *nominal* rate -- the peak gradient is higher',
                ],
            )
            raise ConfigurationError(msg)
        wanted = (
            float(dwell_s) if dwell_s is not None
            else 1.0 / (require_positive(bandwidth_hz_px, 'bandwidth_hz_px')
                        * self.matrix[0] * self.oversampling)
        )
        adc_ticks = to_ticks(float(self.opts.adc_raster_time))
        pair_ticks = 2 * to_ticks(float(self.opts.rf_raster_time))
        quantum = adc_ticks * (pair_ticks // gcd(self.num_samples * adc_ticks, pair_ticks))
        asked = to_ticks(require_positive(wanted, 'dwell_s'))
        return from_ticks(max(quantum, -(-asked // quantum) * quantum))

    def _resolve_samples(self) -> tuple[int, int]:
        """
        Return ``(num_samples, pre_echo_samples)``: an **even** count, times `oversampling`.

        Even because the sampling window has to be exactly centred in its lobe and the lobe has to
        land on the gradient raster: ``T - N*dwell`` is then an even number of ADC rasters, so the
        guard at each end is a legal ADC delay.  An odd count leaves it half a raster short, and
        the two polarities sample grids 50 ns apart -- which is an N/2 ghost of the sequence's own
        making, in the one measurement an EPI exists to get right.

        Both numbers scale with `oversampling`, because both count *ADC samples* while the
        arithmetic above them counts reconstructed pixels.
        """
        wanted = self.partial_fourier * self.matrix[0]
        count = self.oversampling * max(2, 2 * int(round(wanted / 2)))
        pre_echo = count - self.oversampling * (self.matrix[0] - self.matrix[0] // 2)
        if pre_echo < 0:
            msg = format_error(
                f'partial_fourier = {self.partial_fourier:g} leaves no pre-echo samples at '
                f'matrix[0] = {self.matrix[0]}.',
                {'partial_fourier': self.partial_fourier, 'matrix': self.matrix,
                 'num_samples': count},
                ['partial_fourier below ~0.55 has too little pre-echo data for a phase estimate'],
            )
            raise ConfigurationError(msg)
        return count, pre_echo

    def _design_blip(self, requested_s: float | None) -> Event:
        """
        Return the blip at the largest step, on an **even** number of gradient rasters.

        Even because it is centred on the seam: it starts at ``seam - blip/2``, the seam is on the
        gradient raster, and :func:`~seqcraft.compile` refuses a gradient that starts off it.  An
        odd raster count is the one way to make this design fail at compile time rather than in
        the image -- the good kind of failure, but still worth not having.
        """
        area = self.blip_lines * self.dk_blip_per_m
        floor_s = float(pp.calc_duration(
            pp.make_trapezoid(channel=self.blip_axis, area=area, system=self.opts)))
        wanted = floor_s
        if requested_s is not None:
            wanted = require_positive(requested_s, 'blip_duration_s')
            if wanted < floor_s - EPS:
                msg = format_error(
                    f'blip_duration_s = {wanted * 1e6:.1f} us is shorter than the minimum a '
                    f'{self.blip_lines}-line step needs.',
                    {'blip_duration_s': wanted, 'min_blip_duration_s': floor_s,
                     'blip_lines': self.blip_lines, 'blip_area_per_m': area},
                    [
                        f'pass blip_duration_s >= {floor_s:.6g}',
                        'or blip_duration_s=None for the shortest legal blip',
                    ],
                )
                raise ConfigurationError(msg)
        duration = 2 * ceil_raster(wanted / 2, self.opts.grad_raster_time)
        return pp.make_trapezoid(
            channel=self.blip_axis, area=area, duration=duration, system=self.opts,
        )

    def _design_lobe(self) -> tuple[Event, Event, float]:
        """
        Return the positive readout lobe, its ADC, and the guard the two agreed on.

        The lobe duration is fixed first -- ``T = 2*guard + num_samples*dwell``, rounded up onto
        the gradient raster with the rounding split evenly between the two guards -- and the
        amplitude is what is then solved for.  Doing it the other way round, amplitude first and
        duration second, is what makes the window land off-centre.
        """
        sample_time = from_ticks(self.num_samples * to_ticks(self.dwell_s))
        # The window traverses one dk per *reconstructed* sample, so oversampling changes the
        # sample spacing and not the extent -- which is what makes it free in echo spacing.
        k_window = self.num_samples / self.oversampling * self.dk_read_per_m
        slew = float(self.opts.max_slew)

        if self.ramp_sampling:
            total_s = ceil_raster(sample_time + 2 * self._guard_min_s, self.opts.grad_raster_time)
            guard_s = self._halve(total_s, sample_time)
            rise_s = self._solve_ramp(total_s, guard_s, k_window)
        else:
            # Flat-top sampling: the amplitude is the nominal one by definition, and the ramp is
            # whichever is longer -- what the slew rate needs, or what the blip needs to fit
            # beside a closed receiver.
            amplitude = self.dk_read_per_m / (self.dwell_s * self.oversampling)
            rise_s = ceil_raster(
                max(amplitude / slew, self._guard_min_s), self.opts.grad_raster_time,
            )
            total_s = 2 * rise_s + ceil_raster(sample_time, self.opts.grad_raster_time)
            guard_s = self._halve(total_s, sample_time)

        gx = self._lobe_at(total_s, rise_s, guard_s, k_window)
        self._check_hardware(gx, rise_s, guard_s, total_s)
        adc = pp.make_adc(
            num_samples=self.num_samples, dwell=self.dwell_s, delay=guard_s, system=self.opts,
        )
        if abs(float(adc.delay) - guard_s) > EPS:
            msg = format_error(
                f'the receiver needs {float(adc.delay) * 1e6:.1f} us before the first sample and '
                f'the lobe leaves {guard_s * 1e6:.1f} us.',
                {'adc_dead_time': float(self.opts.adc_dead_time), 'guard_s': guard_s,
                 'blip_duration_s': self.blip_duration_s},
                [
                    'lengthen the blip with blip_duration_s, which lengthens the guard with it',
                    'or lower bandwidth_hz_px, which lengthens the lobe',
                ],
            )
            raise ConfigurationError(msg)
        return gx, adc, guard_s

    def _halve(self, total_s: float, sample_time: float) -> float:
        """
        Return ``(T - N*dwell) / 2`` **in integer ticks**, and refuse a half that is not a delay.

        Not a float subtraction.  ``(470 us - 409.6 us) / 2`` evaluates to 30.200000000000003 us,
        which is 302.00000000000006 ADC rasters, and pypulseq's timing check rejects it -- naming
        the block rather than the femtosecond.  This is the drift
        :mod:`seqcraft.design.timing` exists to remove, and it is the one place in this module
        where the arithmetic has to leave floating point.

        The refusal is Rule 1's parity condition made concrete: the guard is half of a gap that
        has to land on the ADC raster, so the sampling duration must be an **even** number of ADC
        rasters.  An even :attr:`num_samples` gives that for any legal dwell, so nothing this
        module builds should reach it -- which is the point of stating it rather than assuming it.
        """
        raster_ticks = to_ticks(float(self.opts.adc_raster_time))
        gap = to_ticks(total_s) - to_ticks(sample_time)
        if gap % (2 * raster_ticks):
            msg = format_error(
                f'the guard would be {from_ticks(gap) / 2 * 1e6:.4f} us, which is not a whole '
                f'ADC raster.',
                {'num_samples': self.num_samples, 'dwell_s': self.dwell_s,
                 'echo_spacing_s': total_s},
                ['num_samples * dwell_s must be an even number of ADC rasters'],
            )
            raise ConfigurationError(msg)
        return from_ticks(gap // 2)

    def _lobe_at(self, total_s: float, rise_s: float, guard_s: float, k_window: float) -> Event:
        """
        Return the lobe of duration `total_s` whose sampling window traverses exactly `k_window`.

        The amplitude is *measured* into place rather than derived: build at any amplitude,
        integrate the window from the event's own knots, and rescale by the ratio.  Area is
        exactly linear in amplitude at fixed ramp times, so the second build is the answer -- and
        it absorbs every raster rounding above it in one step instead of leaving each of them
        somewhere in the trajectory.
        """
        flat_s = total_s - 2 * rise_s
        probe = self._trapezoid(k_window / total_s, rise_s, flat_s)
        got = area_until(probe, guard_s + self.num_samples * self.dwell_s) - area_until(
            probe, guard_s)
        return self._trapezoid(float(probe.amplitude) * k_window / got, rise_s, flat_s)

    def _trapezoid(self, amplitude: float, rise_s: float, flat_s: float) -> Event:
        """
        Build one readout lobe, checking the amplitude **before** pypulseq does.

        ``pp.make_trapezoid`` refuses an over-limit amplitude with ``ValueError: Amplitude
        violation (133%)`` and nothing else -- no bandwidth, no field of view, no note that ramp
        sampling puts the peak above the nominal ``dk/dwell``, and no indication of which of this
        module's several trapezoids it came from.  :meth:`_check_amplitude` says all of that, so
        it has to run first: a refusal that arrives as somebody else's exception is a refusal the
        caller cannot act on.
        """
        self._check_amplitude(abs(float(amplitude)))
        return pp.make_trapezoid(
            channel=self.axis, amplitude=amplitude, rise_time=rise_s, flat_time=flat_s,
            fall_time=rise_s, system=self.opts,
        )

    def _solve_ramp(self, total_s: float, guard_s: float, k_window: float) -> float:
        """
        Return the ramp time of a ramp-sampled lobe, in seconds.

        The closed form is the seed and the search is the answer.  Over the window
        ``[g, T - g]`` a symmetric trapezoid of amplitude ``G`` and ramp ``G/S`` traverses::

            G*T - G^2/S - S*g^2

        -- the ``S*g^2`` term being the two ramp corners the receiver misses, which is independent
        of ``G`` -- so setting that equal to ``N*dk`` is a quadratic whose smaller root is the
        gentler of the two lobes that would work.  Rounding its ramp up onto the gradient raster
        then *raises* the amplitude the window needs, which can raise it past the slew rate, so
        the raster search below is not a formality: it is the difference between a lobe that
        clears the limit and one that misses it by a per cent.
        """
        slew = float(self.opts.max_slew)
        raster = float(self.opts.grad_raster_time)
        discriminant = total_s * total_s - 4 * (k_window + slew * guard_s * guard_s) / slew
        if discriminant < 0.0:
            reach = slew * total_s * total_s / 4 - slew * guard_s * guard_s
            # The largest dwell-independent statement that can be made about it: what the same
            # window would reach at the amplitude that makes the discriminant vanish.
            feasible = (2 * math.sqrt((k_window + slew * guard_s * guard_s) / slew)
                        - 2 * guard_s) / self.num_samples
            msg = format_error(
                f'the readout cannot reach {k_window:.1f} 1/m inside a {total_s * 1e6:.0f} us '
                f'lobe -- the most this slew rate can traverse there is {reach:.1f} 1/m.',
                {'bandwidth_hz_px': self.bandwidth_hz_px, 'dwell_s': self.dwell_s,
                 'guard_s': guard_s, 'fov_mm': self.fov_mm[0], 'matrix': self.matrix},
                [
                    f'lower bandwidth_hz_px to at most '
                    f'{1.0 / (feasible * self.matrix[0] * self.oversampling):.0f}, which '
                    f'lengthens the lobe',
                    'or *lengthen* blip_duration_s: a longer blip is a longer guard and a longer '
                    'lobe, and the lobe grows faster than the ramp corners the receiver misses',
                    'or ramp_sampling=False, which is longer still but asks less of the ramps',
                ],
            )
            raise ConfigurationError(msg)

        amplitude = 0.5 * slew * (total_s - math.sqrt(discriminant))
        rise_s = max(raster, ceil_raster(amplitude / slew, raster))
        for _ in range(_MAX_RAMP_STEPS):
            if 2 * rise_s > total_s + EPS:
                break
            lobe = self._lobe_at(total_s, rise_s, guard_s, k_window)
            if abs(float(lobe.amplitude)) <= slew * rise_s * (1 + 1e-9):
                return rise_s
            rise_s += raster
        msg = format_error(
            f'no ramp between {raster * 1e6:.0f} us and {total_s * 1e6 / 2:.0f} us reaches '
            f'{k_window:.1f} 1/m inside the lobe without exceeding the slew rate.',
            {'bandwidth_hz_px': self.bandwidth_hz_px, 'guard_s': guard_s,
             'max_slew': float(self.opts.max_slew)},
            [
                'lower bandwidth_hz_px, which lengthens the lobe and leaves more room',
                'or ramp_sampling=False',
            ],
        )
        raise ConfigurationError(msg)

    def _check_hardware(
        self, gx: Event, rise_s: float, guard_s: float, total_s: float,
    ) -> None:
        """Refuse a lobe the amplifier cannot play, before the compiler sees the whole train."""
        amplitude = abs(float(gx.amplitude))
        self._check_amplitude(amplitude)
        if amplitude > float(self.opts.max_slew) * rise_s * (1 + 1e-9):
            msg = format_error(
                f'the readout lobe slews at {amplitude / rise_s:.0f} Hz/m/s, above the '
                f'{float(self.opts.max_slew):.0f} Hz/m/s limit.',
                {'rise_time_s': rise_s, 'lobe_duration_s': total_s, 'guard_s': guard_s},
                ['lower bandwidth_hz_px, which lengthens the lobe'],
            )
            raise ConfigurationError(msg)

    def _check_amplitude(self, amplitude: float) -> None:
        """Refuse a readout amplitude above `max_grad`, naming what to change about the protocol."""
        if amplitude <= float(self.opts.max_grad) * (1 + 1e-9):
            return
        nominal = self.dk_read_per_m / self.dwell_s
        msg = format_error(
            f'the readout needs {amplitude:.0f} Hz/m, above the '
            f'{float(self.opts.max_grad):.0f} Hz/m limit.',
            {'bandwidth_hz_px': self.bandwidth_hz_px, 'fov_mm': self.fov_mm[0],
             'nominal_amplitude_hz_m': nominal, 'ramp_sampling': self.ramp_sampling},
            [
                'lower bandwidth_hz_px, or widen fov_mm',
                'with ramp_sampling the peak is above the nominal dk/dwell, because the '
                'window averages the ramps in',
            ],
        )
        raise ConfigurationError(msg)

    def _sampled_k(self) -> np.ndarray:
        """
        Return where the samples of a forward lobe land in ``k``, relative to the echo sample.

        Integrated per sample from the lobe's own knots rather than from ``dk`` times an index,
        because with ramp sampling those are different numbers -- and the whole point of exposing
        this is that a reconstruction can grid from it.  Uniform sampling comes out of the same
        expression as the special case it is.
        """
        times = float(self.adc.delay) + (np.arange(self.num_samples) + 0.5) * self.dwell_s
        origin = area_until(self.gx, self._echo_in_lobe_s)
        return np.array([area_until(self.gx, float(t)) - origin for t in times])

    def _design_prephasers(
        self, requested_s: float | None,
    ) -> tuple[PhaseEncode, float, Event | None]:
        """
        Return the blip-axis prephaser, the common winder duration and the readout-axis one.

        The two are coupled the way :class:`~seqcraft.modules.GRE2DTR`'s three are: each reports
        its own minimum, the longer wins, and the shorter is stretched to match rather than
        followed by a delay.  With no navigators they play at the same instant, so there is one
        window; with navigators there are two, and they are given the same length so that
        `prephaser_duration_s` means one thing.
        """
        area = -self.area_to_echo_per_m
        shortest = pp.make_trapezoid(channel=self.axis, area=area, system=self.opts)
        probe = PhaseEncode(opts=self.opts, fov_mm=self.fov_mm[1], matrix=self.matrix[1],
                            axis=self.blip_axis)
        floor_s = ceil_raster(
            max(float(pp.calc_duration(shortest)), probe.min_duration_s),
            self.opts.grad_raster_time,
        )
        if not self.prephase:
            if requested_s is not None:
                msg = format_error(
                    f'prephaser_duration_s = {requested_s * 1e6:.1f} us was passed with '
                    f'prephase=False, so it cannot take effect.',
                    {'prephase': False, 'area_to_echo_per_m': self.area_to_echo_per_m},
                    [
                        'drop prephaser_duration_s: with prephase=False there is no prephaser '
                        'to stretch',
                        'or pass prephase=True, whose winder duration it then sets',
                    ],
                )
                raise ConfigurationError(msg)
            return probe, 0.0, None

        if requested_s is None:
            wanted = floor_s
        else:
            wanted = ceil_raster(
                require_positive(requested_s, 'prephaser_duration_s'), self.opts.grad_raster_time,
            )
            # EPS, not exact: a caller passing back a maximum of this minimum and another
            # module's can land one ulp below it after snapping onto the raster.
            if wanted < floor_s - EPS:
                msg = format_error(
                    f'prephaser_duration_s = {requested_s * 1e6:.1f} us is shorter than the '
                    f'minimum this readout and its blip axis need.',
                    {'prephaser_duration_s': requested_s, 'min_prephaser_duration_s': floor_s,
                     'prephaser_area_per_m': area},
                    [
                        f'pass prephaser_duration_s >= {floor_s:.6g}',
                        'or pass None for the shortest legal pair',
                    ],
                )
                raise ConfigurationError(msg)
        prephaser = pp.make_trapezoid(
            channel=self.axis, area=area, duration=wanted, system=self.opts,
        )
        blip_axis_prephaser = PhaseEncode(
            opts=self.opts, fov_mm=self.fov_mm[1], matrix=self.matrix[1], axis=self.blip_axis,
            duration_s=wanted,
        )
        return blip_axis_prephaser, wanted, prephaser

    # ------------------------------------------------------------------------ timing
    @property
    def _nav_start_s(self) -> float:
        """When the first lobe of the train -- navigator or imaging -- starts."""
        return self.winder_s if self.prephase else 0.0

    @property
    def _train_start_s(self) -> float:
        """When the first *imaging* lobe starts."""
        blip_winder = self.winder_s if (self.prephase and self.navigator_echoes) else 0.0
        return self._nav_start_s + self.navigator_echoes * self.echo_spacing_s + blip_winder

    # ------------------------------------------------------------------------ per call
    def _lobe(self, index: int) -> Event:
        """Return the readout lobe for train position `index`, sign alternating from zero."""
        return self.gx if index % 2 == 0 else self._gx_reverse

    def _blip_for(self, step: int) -> Event:
        """Return the blip that steps `step` lines along the blip axis, at the fixed duration."""
        scale = step * self.dk_blip_per_m / float(self._blip.area)
        return derive(pp.scale_grad(self._blip, scale))

    # -------------------------------------------------------------------- the refusals
    def _check(self, lines: Iterable[int]) -> tuple[int, ...]:
        """Return `lines` as a tuple, refusing what this train cannot acquire."""
        table = tuple(int(line) for line in lines)
        ny = self.matrix[1]
        if not table:
            msg = format_error(
                'lines is empty, so there is nothing to acquire.',
                {'matrix': self.matrix},
                [f'lines=range({ny}) is a fully sampled single shot'],
            )
            raise ConfigurationError(msg)

        outside = sorted({line for line in table if not 0 <= line < ny})
        if outside:
            msg = format_error(
                f'{len(outside)} line(s) fall outside 0 ... {ny - 1}.',
                {'offending': outside[:8], 'matrix': self.matrix},
                [
                    'lines are zero-based blip-axis indices, not k-space positions',
                    f'the centre of k-space is line {self.center_line}',
                ],
            )
            raise ConfigurationError(msg)

        repeated = sorted({line for line in table if table.count(line) > 1})
        if repeated:
            msg = format_error(
                f'{len(repeated)} line(s) are acquired more than once in this shot.',
                {'repeated': repeated[:8]},
                [
                    'two readouts writing one k-space address is what the compiler rejects '
                    'downstream, with less to say about it than here',
                    'averaging is a separate acquisition and a separate label',
                ],
            )
            raise ConfigurationError(msg)

        steps = [b - a for a, b in zip(table, table[1:])]
        too_far = sorted({step for step in steps if abs(step) > self.blip_lines})
        if too_far:
            msg = format_error(
                f'{len(too_far)} step(s) of this table exceed the {self.blip_lines}-line blip '
                f'this instance was designed for.',
                {'offending_steps': too_far[:8], 'blip_lines': self.blip_lines,
                 'blip_duration_s': self.blip_duration_s},
                [
                    f'pass blip_lines={max(abs(step) for step in too_far)} to design for this '
                    f'table -- it lengthens the blip, the guard and the echo spacing with it',
                    'or reorder the shot so that consecutive echoes are closer in ky',
                    'a blip is scaled per step and never lengthened per step, because an echo '
                    'spacing that varied with the table would move every echo after it',
                ],
            )
            raise ConfigurationError(msg)
        return table

    def _check_echo(self, echo: int) -> int:
        """Return `echo` as an offset from the first imaging echo, refusing one before the train."""
        index = int(echo)
        if index < -self.navigator_echoes:
            msg = format_error(
                f'echo = {index} is before the start of the train.',
                {'echo': index, 'navigator_echoes': self.navigator_echoes},
                [
                    'imaging echoes are counted from zero',
                    f'a navigator is a negative offset from the first imaging echo, so this '
                    f'train addresses them from {-self.navigator_echoes} to -1'
                    if self.navigator_echoes else 'this train builds no navigator echoes',
                ],
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

    def _require_prephase(self, what: str) -> float:
        """Return the winder duration, or refuse a question about a window that was not built."""
        if not self.prephase:
            msg = format_error(
                f'this readout has no prephaser, because prephase=False, so {what} is undefined.',
                {'prephase': False, 'area_to_echo_per_m': self.area_to_echo_per_m},
                [
                    'read area_to_echo_per_m and k_blip_per_m(line) instead: they are what the '
                    'two dephasers have to cancel, and a caller placing them before a refocusing '
                    'pulse needs them with the opposite sign',
                    'or pass prephase=True for a train that prephases itself',
                ],
            )
            raise ConfigurationError(msg)
        return self.winder_s


def _require_count(value: int, name: str, *, low: int = 1, hint: str = '') -> int:
    """
    Return `value` as an int at or above `low`.

    Three arguments here are counts with a floor -- `oversampling`, `blip_lines` and
    `navigator_echoes` -- and one function is better than three that differ only in the noun.
    Private rather than in ``_support`` because those three are its only callers, which is the
    rule that keeps that file the shared things and not the leftovers.

    Examples
    --------
    >>> _require_count(4, 'blip_lines')
    4
    >>> _require_count(-1, 'navigator_echoes', low=0)
    Traceback (most recent call last):
        ...
    seqcraft.errors.ConfigurationError: navigator_echoes must not be negative, got -1.
    """
    count = int(value)
    if count < low:
        floor = 'not be negative' if low == 0 else f'be at least {low}'
        msg = format_error(f'{name} must {floor}, got {count}.', {name: count},
                           [hint] if hint else ())
        raise ConfigurationError(msg)
    return count
