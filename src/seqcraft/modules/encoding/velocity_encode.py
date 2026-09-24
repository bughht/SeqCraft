r"""
:class:`VelocityEncode` -- a bipolar pair that encodes velocity and leaves stationary spins alone.

What it is
----------
Two gradient lobes of equal area and opposite polarity on one logical axis:

.. code-block:: text

    -A          +A                 polarity = +1
    ____        ####
        \      /
         \____/                    zeroth moment zero, first moment +delta_m1 / 2

A spin standing still sees a net area of zero and comes out with the phase it went in with.  A
spin moving along the axis is somewhere else for the second lobe than it was for the first, so the
two lobes do not cancel and it comes out with a phase proportional to its velocity.

The module emits **one** of two toggles per call.  Subtracting the phase of an acquisition made
with ``polarity=+1`` from one made with ``polarity=-1`` cancels everything that is not velocity --
coil phase, off-resonance, susceptibility -- and leaves a map of it.  That subtraction is the
caller's, and so is which toggle is played when.

The public quantity, and the factor of two
------------------------------------------
``venc_m_s`` is the velocity that accumulates exactly :math:`\pi` of phase **difference between
the toggles**.  The quantity that produces it is the *change* in first moment:

.. math::

    \Delta m_1 = \frac{1}{2\,\mathrm{venc}}

Each toggle carries half of that, :math:`\pm \Delta m_1 / 2`, so **a caller who sets one toggle's
first moment to** :math:`\Delta m_1` **is wrong by a factor of two** -- which is the single easiest
mistake to make here, and it produces a perfectly plausible velocity map at half the stated VENC.
The module takes ``venc_m_s`` and never asks for a moment, so that the mistake has nowhere to
happen.

No :math:`\gamma` appears in that relation, because SeqCraft's gradients are in Hz/m rather than
T/m.  The phase a moving spin accumulates is :math:`\phi = 2\pi\,(m_0 x_0 + m_1 v)`, so with
:math:`m_0 = 0` a first moment of :math:`\Delta m_1 = 1/(2\,\mathrm{venc})` s/m gives
:math:`\Delta\phi = \pi` at :math:`v = \mathrm{venc}` for **any** nucleus.  The Handbook writes the
same number as :math:`\pi / (\gamma\,\mathrm{VENC})` because its gradients are in T/m; the
:math:`\gamma` is the unit conversion, not physics this module has to know.

``polarity``, and what its sign means
-------------------------------------
``polarity`` is **the sign of the emitted first moment**, not the sign of a velocity and not the
sign of a reconstructed phase:

.. code-block:: text

    polarity = +1    m1 = +delta_m1 / 2    the negative lobe is played first
    polarity = -1    m1 = -delta_m1 / 2    the positive lobe is played first

Defining it on the emitted waveform is what makes it checkable: :func:`~seqcraft.moments` reads
the sign back off the events, and nothing has to agree a convention with the reconstruction to do
it.  Which sign of measured flow a reader calls *forward* is a reconstruction choice, and it is
downstream of this module.

Why the first moment needs no time origin
-----------------------------------------
A first moment normally depends on where you start the clock.  This one does not: once
:math:`m_0 = 0`, shifting the origin changes :math:`m_1` by :math:`m_0 \Delta t = 0`.  So the pair
can be placed anywhere in a repetition and still deliver the same :math:`\Delta m_1`, and the
module never has to be told where the echo is.  That is the property that makes this a leaf rather
than something that needs to negotiate with its neighbours.

How the pair is designed
------------------------
Shortest first, in closed form rather than by search.  Each lobe is a trapezoid of area `A` with
rise time `r` and flat time `f`, and the two are played back to back, so their centres are
:math:`\Delta T = 2r + f` apart and

.. math::

    |m_1| = A \, \Delta T = \frac{\Delta m_1}{2}

Two regimes, and which one applies is decided rather than searched:

============================  =================================================================
triangular, ``f = 0``         when the target fits below the gradient maximum.  The amplitude is
                              the real root of :math:`2 g^3 / s^2 = A\,\Delta T`
trapezoidal, ``g = max_grad`` otherwise.  `f` is the positive root of
                              :math:`f^2 + 3 r f + 2 r^2 - A\Delta T / g = 0`
============================  =================================================================

Both times are then **rounded up** onto the gradient raster and the amplitude is solved again
against those realised times, which lands it below the maximum rather than above it and makes
:math:`\Delta m_1` exact instead of approximately right.  :attr:`VelocityEncode.amplitude_hz_per_m`
reports what that came to.

Ramping down through zero and back up at the seam costs exactly as much as ramping straight from
:math:`+g` to :math:`-g` would -- two rise times either way -- so the two-lobe form is not paying
for its simplicity.

Out of scope, each a different claim rather than a parameter of this one: the **merged** designs,
where the encoding lobe and an already-compensated imaging lobe are combined into one waveform,
which reshapes gradients this module does not own; acceleration and higher-order encoding;
concomitant-field phase, which the simple bipolar cancels automatically and which is still not
claimed because it was not measured; and any statement about achievable minimum TE.

References
----------
Bernstein, King & Zhou, *Handbook of MRI Pulse Sequences*, Elsevier 2004, §9.2 "Flow-Encoding
Gradients", pp. 281-291 -- the bipolar pair, the VENC relation, and per-axis independence.
Simonetti et al., *J. Comput. Assist. Tomogr.* **15**, 1051 (1991), via §10.4 p. 343 -- the first
moment is origin-independent once the zeroth is zero.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ...design.logic import LogicBlock
from ...design.module import Module
from ...design.timing import Raster
from ...errors import ConfigurationError, format_error
from .._support import require_axis, require_positive

if TYPE_CHECKING:
    from pypulseq.opts import Opts

    from ...design.events import Event

__all__ = ['VelocityEncode']

#: The two toggles this module knows how to emit, as the sign of the first moment they carry.
_POLARITIES = (1, -1)


class VelocityEncode(Module):
    r"""
    A bipolar velocity-encoding pair on one axis, emitted one toggle at a time.

    Parameters
    ----------
    opts
        The scanner.  ``max_grad`` and ``max_slew`` are what the closed form is solved against.
    venc_m_s
        **The physical target.**  The velocity along `axis` that accumulates :math:`\pi` of phase
        *difference* between the two toggles.  A protocol normally sets it somewhat above the
        expected peak velocity, because flow faster than ``venc_m_s`` wraps.
    axis
        Logical gradient channel.  Each axis is encoded independently, so a through-plane study
        uses one of these and a 4D-flow study uses three, one module each.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    venc_m_s : float
        What was asked for, and what is delivered -- the design solves for it exactly.
    delta_m1_s_per_m : float
        The change in first moment between the toggles, ``1 / (2 * venc_m_s)``.  Each toggle
        carries half of it, with the sign of its `polarity`.
    axis : str
        The channel the pair plays on.
    amplitude_hz_per_m : float
        The lobe amplitude the closed form arrived at, at or below ``opts.max_grad``.
    rise_time_s, flat_time_s : float
        One lobe's ramp and plateau, both on the gradient raster.  ``flat_time_s`` is zero in the
        triangular regime.
    lobe_area_1_per_m : float
        One lobe's area.  The two differ only in sign, which is what makes ``m0`` zero.
    duration_s : float
        The whole pair, both lobes.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s')
    >>> ve = sc.modules.VelocityEncode(opts=opts, venc_m_s=1.5, axis='z')
    >>> round(ve.delta_m1_s_per_m, 6)                       # 1 / (2 * 1.5)
    0.333333

    ``polarity`` is the sign of the first moment, read back off the emitted events:

    >>> round(sc.moments(ve(polarity=+1), 1)['z'], 6)
    0.166667
    >>> round(sc.moments(ve(polarity=-1), 1)['z'], 6)
    -0.166667

    Stationary spins see nothing, which is the other half of the contract:

    >>> {round(sc.moments(ve(polarity=p), 0)['z'], 12) for p in (1, -1)}
    {0.0}
    """

    def __init__(
        self,
        *,
        opts: Opts,
        venc_m_s: float,
        axis: str = 'z',
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.venc_m_s = require_positive(venc_m_s, 'venc_m_s')
        self.axis = require_axis(axis)
        self._raster = Raster(float(opts.grad_raster_time), 'grad_raster_time')

        # delta_m1 is the CHANGE across the toggles, so one toggle carries half of it -- and that
        # half is what the geometry has to produce.  See the module docstring on the factor of two.
        self.delta_m1_s_per_m = 1.0 / (2.0 * self.venc_m_s)
        self._target = self.delta_m1_s_per_m / 2.0
        self.rise_time_s, self.flat_time_s = self._solve_times()
        self.amplitude_hz_per_m = self._amplitude_for(self.rise_time_s, self.flat_time_s)
        self._check_hardware()

        self.lobe_area_1_per_m = self.amplitude_hz_per_m * (self.flat_time_s + self.rise_time_s)
        self._lobe_s = 2.0 * self.rise_time_s + self.flat_time_s
        self.duration_s = 2.0 * self._lobe_s

    # ------------------------------------------------------------------ what it knows
    def first_moment_s_per_m(self, polarity: int = 1) -> float:
        """
        The first moment one toggle carries, in s/m -- ``polarity * delta_m1_s_per_m / 2``.

        No time origin is named because none is needed: the pair's zeroth moment is zero, so
        shifting the origin changes the first moment by ``m0 * dt == 0``.  That is why this is a
        number the module can state rather than a question about where it was placed.
        """
        return self._check_polarity(polarity) * self._target

    def lobe_separation_s(self) -> float:
        """Seconds between the two lobes' centres, which is one lobe's duration."""
        return self._lobe_s

    # ----------------------------------------------------------------------- assembly
    def build(self, *, polarity: int = 1) -> LogicBlock:
        """
        Return one toggle of the pair.

        Parameters
        ----------
        polarity
            ``+1`` or ``-1``, **the sign of the emitted first moment**.  ``+1`` plays the negative
            lobe first, because a negative lobe followed by a positive one is what puts a moving
            spin ahead rather than behind.
        """
        sign = self._check_polarity(polarity)
        out = LogicBlock()
        at = 0.0
        for lobe in (-sign, sign):
            out.add(at, self._lobe(lobe))
            at += self._lobe_s
        return out

    def _lobe(self, sign: int) -> Event:
        """One trapezoid, at the solved amplitude and on the solved times."""
        return pp.make_trapezoid(
            channel=self.axis,
            amplitude=sign * self.amplitude_hz_per_m,
            rise_time=self.rise_time_s,
            flat_time=self.flat_time_s,
            fall_time=self.rise_time_s,
            system=self.opts,
        )

    # ------------------------------------------------------------------- the closed form
    def _solve_times(self) -> tuple[float, float]:
        """
        Return the rise and flat times of one lobe, rounded up onto the gradient raster.

        Shortest first, and the regime is decided rather than searched.  With ``A = g (f + r)``
        and a centre-to-centre separation of ``dT = 2r + f``, the requirement ``A dT = target``
        is a cubic in `g` when ``f == 0`` and a quadratic in `f` when `g` is pinned at the
        maximum.  The triangular form is enough exactly when a full-amplitude triangle already
        overshoots the target.
        """
        g_max, s_max = float(self.opts.max_grad), float(self.opts.max_slew)
        if self._target <= 2.0 * g_max ** 3 / s_max ** 2:
            rise = (self._target * s_max ** 2 / 2.0) ** (1.0 / 3.0) / s_max
            flat = 0.0
        else:
            rise = g_max / s_max
            flat = (-3.0 * rise + np.sqrt(rise * rise + 4.0 * self._target / g_max)) / 2.0
        return float(self._raster.ceil(rise)), float(self._raster.ceil(flat))

    def _amplitude_for(self, rise_s: float, flat_s: float) -> float:
        """
        Solve the amplitude against the **realised** times, so that ``delta_m1`` is exact.

        Rounding the times up and then re-solving is what makes this a design rather than an
        approximation: the target is met to floating point, and because both times grew, the
        amplitude this returns is below the one they were solved from -- so a rounding can never
        push the gradient over the limit.
        """
        return float(self._target / ((flat_s + rise_s) * (2.0 * rise_s + flat_s)))

    # -------------------------------------------------------------------- the refusals
    def _check_polarity(self, polarity: int) -> int:
        """Return `polarity` having checked it is one of the two toggles."""
        if polarity in _POLARITIES:
            return int(polarity)
        msg = format_error(
            f'polarity must be +1 or -1, got {polarity!r}.',
            {'polarity': polarity},
            ['polarity is the sign of the emitted first moment, not of a velocity',
             'a phase-contrast pair is one acquisition at each, subtracted'],
        )
        raise ConfigurationError(msg)

    def _check_hardware(self) -> None:
        """Refuse a pair the amplifier cannot play, naming the limit that binds."""
        g_max, s_max = float(self.opts.max_grad), float(self.opts.max_slew)
        slew = self.amplitude_hz_per_m / self.rise_time_s
        if self.amplitude_hz_per_m <= g_max and slew <= s_max:
            return
        msg = format_error(
            f'venc_m_s = {self.venc_m_s:g} m/s needs a bipolar pair this gradient system '
            f'cannot play.',
            {'venc_m_s': self.venc_m_s,
             'amplitude_hz_per_m': self.amplitude_hz_per_m, 'max_grad': g_max,
             'slew_hz_per_m_per_s': slew, 'max_slew': s_max},
            ['raise venc_m_s -- a larger VENC needs a smaller first moment',
             'or relax opts.max_grad / opts.max_slew if the scanner really allows it'],
        )
        raise ConfigurationError(msg)
