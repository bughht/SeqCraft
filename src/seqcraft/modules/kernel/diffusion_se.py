r"""
A diffusion-weighted spin echo: the Stejskal--Tanner pair, designed against a requested `b`.

The physical family
-------------------
From Bernstein, King & Zhou, *Handbook of MRI Pulse Sequences* (2004), §9.1
"Diffusion-Weighting Gradients", book pp. 274--281:

    "A diffusion-weighting gradient typically consists of **two lobes with equal area**.  In pulse
    sequences based on spin echoes, the two lobes have the **same polarity** and are placed at
    either side of a refocusing RF pulse."

Same polarity, because the refocusing pulse conjugates what the first lobe accumulated.  That is
why this is a kernel and not a leaf: **the encoding is defined around a pulse that belongs to
another module**, and neither lobe can be designed without knowing where that pulse is and how
long it lasts.

The `b`-value of a trapezoid pair, §9.1 Figure 9.3, with :math:`\delta` the lobe width including
one ramp, :math:`\Delta` the separation of the lobe centres and :math:`\varepsilon` the ramp time:

.. math::

    b = \gamma^2 G^2 \left[ \delta^2\!\left(\Delta - \frac{\delta}{3}\right)
        + \frac{\varepsilon^3}{30} - \frac{\delta\,\varepsilon^2}{6} \right]

The ramp terms are the part the reference implementations leave out -- ``pulseq``'s own
``writeEpiDiffusionRS.m`` carries a helper documented "for trapezoid gradients: TODO" -- and they
are not negligible at the durations diffusion uses.

What this owns, and what it does not
------------------------------------
It owns the coupled design: the excitation, the refocusing pulse, and the pair of lobes whose
width and amplitude follow from the requested `b` and from the windows those two pulses leave.
It does **not** own the readout, the direction schedule, the b-value list, or the averaging -- all
of which are acquisition policy, and all of which a caller writes in four lines.

Realization is analytic, and that is a finding rather than a preference
----------------------------------------------------------------------
With the geometry fixed by the timing, `b` is quadratic in the amplitude, so the amplitude that
reaches a requested `b` is a square root and nothing more.  The minimum echo time for a requested
`b` is then a bisection on one monotone scalar -- as TE grows the windows grow, so the achievable
`b` grows.  **No numerical optimizer is involved**, and none is needed for this family.

That is not a claim that diffusion never needs one.  Several moment orders at once, eddy-current
constraints, fixed waveform segments or a twice-refocused design with its own cancellation
conditions are a different problem, and `tools/module_mining/plans/diffusion/` records where the
line falls.

Validation is independent of the design
---------------------------------------
:func:`seqcraft.b_value` integrates the emitted gradients and applies the refocusing conjugation
itself.  It shares no code with this module, and it measures the **whole tree** -- so the slice
lobes, the crushers and anything a caller adds are all counted, which is what `b` actually means.
This module's own arithmetic is never the thing that checks it.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ...design.logic import LogicBlock
from ...design.module import Module
from ...errors import ConfigurationError, format_error
from .._support import ceil_raster, require_positive
from ..rf.excitation import Excitation
from ..rf.refocusing import Refocusing

if TYPE_CHECKING:
    from pypulseq.opts import Opts

__all__ = ['DiffusionSEPrep']

#: Gyromagnetic ratio of the proton, Hz/T -- the same value `pypulseq` uses.
_GAMMA_HZ_T = 42.576e6

#: The slew fraction diffusion lobes are designed at, by default.  The handbook, §9.1 p. 280:
#: "using the maximum gradient slew rate can also reduce the TE value, [but] the improvement is
#: usually much less than employing the maximum gradient amplitude because the plateau duration
#: typically greatly exceeds the ramp width.  In addition, employing the maximum slew rate
#: contributes to the overall dB/dt value, which may induce unwanted peripheral nerve stimulation
#: or excessive eddy currents.  For these reasons, maximum slew rate is often not used in
#: diffusion-weighting gradients, especially for imaging human subjects."
#:
#: So this is a physical default with a source, not a safety margin someone guessed.
_SLEW_FRACTION = 0.7

#: How close the minimum-TE search gets before it stops, in seconds.  One gradient raster would be
#: pointless -- the result is rounded onto the block raster anyway.
_TE_TOLERANCE_S = 5e-6


def b_of_trapezoid_pair(amplitude_hz_m: float, delta_s: float, separation_s: float,
                        ramp_s: float) -> float:
    r"""
    Return the `b`-value of one trapezoid pair, in s/mm\ :sup:`2`.

    Handbook §9.1, Figure 9.3.  `delta_s` is the lobe width **including one ramp**, `separation_s`
    is the distance between the two lobe centres, and `ramp_s` is the ramp time.  Setting
    ``ramp_s=0`` recovers the rectangular-lobe expression, which is the check the figure's own
    caption suggests.

    A free function because it is the physics, not the module: the tests use it as an independent
    expectation, and it takes no ``self`` that could carry a stale assumption.
    """
    gamma_rad = 2.0 * np.pi * _GAMMA_HZ_T
    grad_t_m = float(amplitude_hz_m) / _GAMMA_HZ_T
    bracket = (delta_s ** 2 * (separation_s - delta_s / 3.0)
               + ramp_s ** 3 / 30.0 - delta_s * ramp_s ** 2 / 6.0)
    return float(gamma_rad ** 2 * grad_t_m ** 2 * bracket / 1e6)


class DiffusionSEPrep(Module):
    r"""
    A diffusion-weighted spin-echo preparation: 90, a lobe, 180, a lobe, and the echo.

    The block runs from the excitation to the end of the second diffusion lobe.  The spin echo is
    at :meth:`time_to_echo`, which is usually **after** the block ends -- a caller places a readout
    so that the readout's own echo lands there, and the gap between is what the readout's
    prephaser needs.

    Parameters
    ----------
    opts
        Scanner limits.  The diffusion lobes are designed at ``max_grad`` and at `slew_fraction`
        of ``max_slew``; everything else uses the full limits.
    thickness_mm
        Slice thickness for the excitation.  The refocusing pulse selects
        ``refocus_thickness_factor`` times this, for the reason ``SE2D`` records.
    b_s_per_mm2
        **The physical target.**  What a caller asks for; the amplitude and the echo time follow.
        ``0.0`` is the b = 0 reference image and is built with no diffusion lobes at all rather
        than with zero-amplitude ones.
    axis
        Which logical axis carries the encoding, or several for an oblique direction.  The
        handbook notes each axis is treated independently, so a direction is a per-axis amplitude
        and the `b` of the set is their sum.
    te_s
        ``None`` -- the default -- uses the shortest echo time that reaches `b_s_per_mm2`.  A
        value is used as given and **refused** if it is too short, rather than silently lengthened.

    Attributes
    ----------
    te_s : float
        The realised echo time, excitation centre to spin echo.
    min_te_s : float
        The shortest echo time that reaches the requested `b`.
    delta_s, separation_s, ramp_s : float
        The lobe width including one ramp, the lobe-centre separation, and the ramp time --
        :math:`\delta`, :math:`\Delta` and :math:`\varepsilon` in the handbook's expression.
    amplitude_hz_m : float
        Per-axis lobe amplitude actually used.  At the minimum echo time this is ``max_grad``; at
        a longer echo time it is whatever reaches the requested `b`.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> opts = pp.Opts(max_grad=80, grad_unit='mT/m', max_slew=200, slew_unit='T/m/s',
    ...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> prep = sc.modules.DiffusionSEPrep(opts=opts, thickness_mm=4.0, b_s_per_mm2=1000.0)
    >>> round(prep.te_s * 1e3, 2) == round(prep.min_te_s * 1e3, 2)
    True
    """

    def __init__(
        self,
        *,
        opts: Opts,
        thickness_mm: float,
        b_s_per_mm2: float,
        axis: str | tuple[str, ...] = 'x',
        te_s: float | None = None,
        flip_deg: float = 90.0,
        refocus_flip_deg: float = 180.0,
        excitation_duration_s: float = 2.5e-3,
        refocus_duration_s: float = 4e-3,
        refocus_thickness_factor: float = 1.25,
        crush_cycles_slice: float = 3.0,
        slew_fraction: float = _SLEW_FRACTION,
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.thickness_mm = require_positive(thickness_mm, 'thickness_mm')
        self.b_s_per_mm2 = float(b_s_per_mm2)
        if self.b_s_per_mm2 < 0.0:
            msg = format_error('b_s_per_mm2 must not be negative.',
                               {'b_s_per_mm2': self.b_s_per_mm2},
                               ['use 0.0 for the unweighted reference image'])
            raise ConfigurationError(msg)
        self.axis = self._check_axes(axis)
        self.slew_fraction = require_positive(slew_fraction, 'slew_fraction')

        self.exc = Excitation(opts=opts, flip_deg=flip_deg, thickness_mm=self.thickness_mm,
                              duration_s=excitation_duration_s)
        self.refoc = Refocusing(opts=opts, thickness_mm=self.thickness_mm * refocus_thickness_factor,
                                flip_deg=refocus_flip_deg, duration_s=refocus_duration_s,
                                crush_cycles_per_voxel=crush_cycles_slice,
                                crush_voxel_mm=self.thickness_mm)

        # --- every number below comes from the two leaves' CONSTRUCTION-TIME properties.  Nothing
        # here calls build(), so the coupled design is complete before any event exists.  That is
        # the architectural claim this module was written to test.
        raster = float(opts.grad_raster_time)
        self._t90c = float(self.exc.time_to_center())
        self._exc_end = ceil_raster(
            float(self.exc.time_to_rephaser()) + float(self.exc.rephaser_duration_s), raster)
        self._refoc_to_center = float(self.refoc.time_to_center())
        self._refoc_block_s = ceil_raster(
            float(self.refoc.time_to_crusher()) + float(self.refoc.crush_duration_s), raster)

        self._amplitude_cap = float(opts.max_grad)
        self.ramp_s = 0.0 if not self.b_s_per_mm2 else ceil_raster(
            self._amplitude_cap / (float(opts.max_slew) * self.slew_fraction), raster)

        # The cubic's root, then rounded UP onto the gradient raster so the lobe starts and ends
        # on it -- an unquantised width puts the first lobe's start off the raster and the
        # compiler refuses the whole tree.  Rounding up over-reaches the requested `b` by a
        # fraction of a per cent, so the amplitude is then trimmed to land on it exactly: `b` is
        # quadratic in the amplitude at fixed geometry, so that is a square root and stays under
        # the cap by construction.
        exact = self._delta_for(self.b_s_per_mm2)
        self.delta_s = 0.0 if exact <= 0.0 else ceil_raster(exact, raster)
        self.separation_s = (0.0 if self.delta_s <= 0.0
                             else self._refoc_block_s + self.delta_s + self.ramp_s)
        self.amplitude_hz_m = self._amplitude_for(self.b_s_per_mm2)
        self.min_te_s = self._shortest_te(self.delta_s, raster)
        self.te_s = self.min_te_s if te_s is None else float(te_s)
        if te_s is not None and self.te_s < self.min_te_s - 1e-12:
            msg = format_error(
                f'te_s = {self.te_s * 1e3:.3f} ms is shorter than this b-value allows.',
                {'te_s_ms': self.te_s * 1e3, 'min_te_ms': self.min_te_s * 1e3,
                 'b_s_per_mm2': self.b_s_per_mm2},
                [f'use te_s >= {self.min_te_s * 1e3:.3f} ms',
                 'or lower b_s_per_mm2',
                 'or raise max_grad -- amplitude buys more than slew rate here'],
            )
            raise ConfigurationError(msg)


    # ------------------------------------------------------------------ the coupled design
    #
    # All of it closed form, and all of it before any event exists.
    #
    # With the amplitude pinned at the scanner maximum -- which is what the handbook says is done,
    # §9.1 p. 280 -- the lobe width follows from the requested `b` by solving a cubic, which is
    # exactly the calculation of the handbook's own Example 9.2.  Substituting the separation
    # `Delta = C + delta + eps`, where `C` is the refocusing block, into Figure 9.3's expression:
    #
    #     b / (gamma^2 G^2) = (2/3) delta^3 + (C + eps) delta^2 - (eps^2/6) delta + eps^3/30
    #
    # The echo time then does not enter the lobe design at all -- it only decides whether the lobe
    # FITS.  Both windows grow as TE/2, so the shortest echo time that holds the pair is a max of
    # two linear expressions.  No search, no optimizer, no iteration.
    def _delta_for(self, b: float) -> float:
        """Return the lobe width reaching `b` at full amplitude, in seconds.  The cubic's root."""
        if b <= 0.0:
            return 0.0
        gamma_rad = 2.0 * np.pi * _GAMMA_HZ_T
        grad_t_m = self._amplitude_cap / _GAMMA_HZ_T
        per_axis = b / len(self.axis)
        scale = per_axis * 1e6 / (gamma_rad ** 2 * grad_t_m ** 2)
        eps, block = self.ramp_s, self._refoc_block_s
        roots = np.roots([2.0 / 3.0, block + eps, -eps ** 2 / 6.0, eps ** 3 / 30.0 - scale])
        real = sorted(float(r.real) for r in roots if abs(r.imag) < 1e-12 and r.real > 0.0)
        if not real:
            msg = format_error(
                f'no lobe width reaches b = {b:.0f} s/mm^2 on this scanner.',
                {'b_s_per_mm2': b, 'max_grad_hz_m': self._amplitude_cap},
                ['lower b_s_per_mm2', 'or raise max_grad'])
            raise ConfigurationError(msg)
        return real[0]

    def _amplitude_for(self, b: float) -> float:
        """
        Return the per-axis amplitude that lands on `b` with the quantised geometry, in Hz/m.

        `b` is quadratic in the amplitude once the widths are fixed, so this is one square root.
        It is always at or below the cap, because the width was rounded UP.
        """
        if b <= 0.0 or self.delta_s <= 0.0:
            return 0.0
        at_cap = b_of_trapezoid_pair(self._amplitude_cap, self.delta_s, self.separation_s,
                                     self.ramp_s) * len(self.axis)
        if at_cap <= 0.0:
            return 0.0
        return self._amplitude_cap * float(np.sqrt(min(1.0, b / at_cap)))

    def _shortest_te(self, delta: float, raster: float) -> float:
        """
        Return the shortest echo time whose windows hold a lobe on each side of the refocusing
        pulse.  A lobe occupies ``delta + ramp`` of wall time, and both windows grow as ``TE/2``.
        """
        occupies = 0.0 if delta <= 0.0 else delta + self.ramp_s
        before = occupies + self._exc_end - self._t90c + self._refoc_to_center
        after = occupies - self._refoc_to_center + self._refoc_block_s
        return ceil_raster(2.0 * max(before, after), raster)

    def _windows(self, te_s: float) -> tuple[float, float]:
        """``(before, after)`` -- the room either side of the refocusing block, at an echo time."""
        refoc_start = self._t90c + te_s / 2.0 - self._refoc_to_center
        return (refoc_start - self._exc_end,
                (self._t90c + te_s) - (refoc_start + self._refoc_block_s))

    # ------------------------------------------------------------------ what it knows
    def time_to_center(self) -> float:
        """Seconds from the block's start to the excitation's effective centre."""
        return self._t90c

    def time_to_echo(self) -> float:
        """
        Seconds from the block's start to the **spin echo**.

        Usually later than the block's own duration: the block ends with the second diffusion
        lobe, and the gap that follows is where a readout's prephaser goes.
        """
        return self._t90c + self.te_s

    def time_to_refocus(self) -> float:
        """Seconds from the block's start to the refocusing pulse's effective centre."""
        return self._t90c + self.te_s / 2.0

    @property
    def lobe_duration_s(self) -> float:
        """How long one diffusion lobe lasts, including both its ramps."""
        return 0.0 if self.delta_s <= 0.0 else self.delta_s + self.ramp_s

    def build(self) -> LogicBlock:
        """Return the preparation: excitation, lobe, refocusing, lobe."""
        raster = float(self.opts.grad_raster_time)
        out = LogicBlock('diffusion_se_prep')
        out.add(0.0, self.exc())
        refoc_start = ceil_raster(self.time_to_refocus() - self._refoc_to_center, raster)
        out.add(refoc_start, self.refoc())
        if self.delta_s > 0.0 and self.amplitude_hz_m > 0.0:
            flat = self.delta_s - self.ramp_s
            for axis in self.axis:
                lobe = pp.make_trapezoid(axis, amplitude=self.amplitude_hz_m,
                                         rise_time=self.ramp_s, flat_time=flat, system=self.opts)
                out.add(refoc_start - self.lobe_duration_s, lobe)
                out.add(refoc_start + self._refoc_block_s, lobe)
        return out

    # ------------------------------------------------------------------ the refusals
    def _check_axes(self, axis: str | tuple[str, ...]) -> tuple[str, ...]:
        """Return the encoding axes, having checked they name real ones."""
        axes = (axis,) if isinstance(axis, str) else tuple(axis)
        bad = [a for a in axes if a not in ('x', 'y', 'z')]
        if bad or not axes:
            msg = format_error(
                f'axis must name one or more of x, y, z; got {axis!r}.', {'axis': axis},
                ["a single axis as 'x', or a direction as ('x', 'y')"])
            raise ConfigurationError(msg)
        return axes
