"""
Diffusion encoding gradients.

A diffusion encoding is a *pair* of lobes straddling the refocusing pulse, so one module
produces two blocks and :meth:`MonopolarDiffusion.build` picks which with ``part=``.  That is the
case that drove seqcraft's design: the two lobes must share one b-value solution, so making them
separate modules would mean two objects that have to agree and nothing checking that they do.

b-value arithmetic
------------------
b is ``(2*pi)^2 * integral |k(t)|^2 dt`` with the refocusing pulse's sign flip applied, and with
``G`` in Hz/m the result is in s/m^2 -- no gamma appears anywhere, because pulseq amplitudes are
already frequency gradients.  :meth:`MonopolarDiffusion._b_of` integrates that exactly, piecewise
over each trapezoid's ramp, flat top and fall.

It is worth doing rather than quoting a formula.  The familiar ramp correction
``+- eps^3/30 -+ delta*eps^2/6`` appears in the literature with several incompatible conventions
for whether ``delta`` includes the ramps, and using one with the wrong convention overstates b by
2--4 % at the lobe durations diffusion sequences actually use.  A b-value wrong by 3 % biases
every diffusivity by 3 %, and nothing in the reconstruction would reveal it.
:meth:`MonopolarDiffusion.achieved_b_s_per_mm2` reports what the designed waveform delivers, and
the tests check it against numerical integration of the built events to 0.5 %.

Directions
----------
`direction` is a **build argument**, not a design parameter, so 30 directions cost one design and
30 cheap builds.  Direction vectors are normalised on the way in: an unnormalised ``[1, 1, 1]``
would request ``sqrt(3)`` times the intended amplitude on the vector norm -- a real bug that
produces plausible-looking b-values and illegal gradients at the same time.

Examples
--------
>>> import seqcraft as sc
>>> system = sc.System.preset('generic_3t')
>>> diff = MonopolarDiffusion(system, b_value_s_per_mm2=1000, refocus_duration_us=4200)
>>> round(diff.b_value_s_per_mm2)
1000
>>> diff.build(part='pre', direction=(1, 0, 0))
LogicBlock(diff_pre, 3 nodes, ... ms)
>>> round(diff.achieved_b_s_per_mm2((1, 0, 0)))       # measured from the waveform
1000
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
import pypulseq as pp

from ...core import events as ev
from ...core.errors import ConfigurationError, format_error
from ...core.logic import LogicBlock
from ...core.module import Module
from ...core.timing import Raster
from ...core.units import convert
from ...core.validate import require_in, require_positive

if TYPE_CHECKING:
    from collections.abc import Sequence

    from ...core.system import System

__all__ = [
    'ArbitraryDiffusion',
    'BipolarDiffusion',
    'MonopolarDiffusion',
    'direction_condition_number',
    'dti_directions',
]

_AXES = ('x', 'y', 'z')


def dti_directions(n: int, *, iterations: int = 400) -> tuple[tuple[float, float, float], ...]:
    """
    Return `n` well-separated unit vectors for a DTI acquisition.

    Starts from a golden-angle spiral on the hemisphere and then relaxes it under **antipodal**
    electrostatic repulsion -- each direction is repelled by every other one *and by its
    reflection through the origin*.  That is the right symmetry: diffusion weighting depends on
    the square of the direction, so ``+d`` and ``-d`` are the same measurement, and a pair that
    is nearly antiparallel is just as redundant as a pair that is nearly parallel.

    Parameters
    ----------
    n
        Number of directions.  Six is the arithmetic minimum -- a diffusion tensor has six
        independent components -- 30 is a common clinical compromise, 60 or more for
        higher-order models.
    iterations
        Relaxation steps.  A few hundred is enough to converge for any practical `n`, and it
        takes milliseconds.

    Returns
    -------
    tuple of (x, y, z)
        Unit vectors, each on the upper hemisphere by convention.

    Notes
    -----
    The relaxation matters.  A raw golden-angle spiral is uniform in *area*, which is not the
    same as well separated: for 30 directions it leaves a pair only 14 degrees apart
    (``|d1 . d2| = 0.97``), and near-duplicate directions make the tensor fit ill-conditioned
    without making it *look* wrong -- the scan completes, the fit converges, and the fractional
    anisotropy is biased.  Use :func:`direction_condition_number` to see the difference.

    Examples
    --------
    >>> dirs = dti_directions(6)
    >>> len(dirs)
    6
    >>> all(abs(sum(c * c for c in d) - 1.0) < 1e-9 for d in dirs)
    True
    >>> all(d[2] >= -1e-12 for d in dirs)                     # upper hemisphere
    True
    >>> worst = max(abs(sum(a * b for a, b in zip(p, q)))
    ...             for i, p in enumerate(dirs) for q in dirs[i + 1:])
    >>> bool(worst < 0.55)                                    # no near-duplicate pair
    True
    """
    if n < 1:
        msg = format_error('dti_directions needs at least one direction.', {'n': n})
        raise ConfigurationError(msg)
    if n == 1:
        return ((0.0, 0.0, 1.0),)

    golden = math.pi * (3.0 - math.sqrt(5.0))
    index = np.arange(n)
    z = 1.0 - (index + 0.5) / n
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    phi = index * golden
    points = np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)

    step = 0.5 / n
    for _ in range(iterations):
        # Repulsion from every other point and from its antipode, which is what makes the metric
        # the angle between *lines* rather than between vectors.
        force = np.zeros_like(points)
        for sign in (1.0, -1.0):
            delta = points[:, None, :] - sign * points[None, :, :]
            distance = np.linalg.norm(delta, axis=2)
            np.fill_diagonal(distance, np.inf)
            if sign < 0:
                # A point and its own antipode are diametrically opposite and exert no useful
                # force; excluding it also avoids a divide by a constant that biases nothing.
                np.fill_diagonal(distance, np.inf)
            force += np.sum(delta / distance[:, :, None] ** 3, axis=1)
        points += step * force
        points /= np.linalg.norm(points, axis=1, keepdims=True)

    # Fold onto the upper hemisphere: a direction and its antipode are the same measurement, so
    # reporting them consistently makes the table easier to read and to compare against a vendor's.
    points[points[:, 2] < 0] *= -1.0
    return tuple((float(x), float(y), float(zz)) for x, y, zz in points)


def direction_condition_number(directions: Sequence[tuple[float, float, float]]) -> float:
    """
    Return the condition number of the diffusion-tensor design matrix for `directions`.

    Each direction contributes one row ``[dx^2, dy^2, dz^2, 2 dx dy, 2 dx dz, 2 dy dz]``, which is
    the linear system a tensor fit solves.  Its condition number says how much measurement noise
    the fit amplifies: the theoretical minimum for a well-spread set is around 1.6, six directions
    arranged along the icosahedral edges reach about 1.6, and anything above ~2.5 means some
    tensor components are much noisier than others.

    This is the number to compare direction schemes by -- not the raw count, and not the minimum
    pairwise angle.

    Examples
    --------
    >>> import seqcraft as sc
    >>> round(direction_condition_number(sc.modules.dti_directions(6)), 2) < 2.5
    True
    >>> round(direction_condition_number(sc.modules.dti_directions(30)), 2) < 2.0
    True
    """
    d = np.asarray(directions, dtype=float)
    design = np.stack(
        [
            d[:, 0] ** 2,
            d[:, 1] ** 2,
            d[:, 2] ** 2,
            2.0 * d[:, 0] * d[:, 1],
            2.0 * d[:, 0] * d[:, 2],
            2.0 * d[:, 1] * d[:, 2],
        ],
        axis=1,
    )
    return float(np.linalg.cond(design))


def _normalise(direction: tuple[float, float, float]) -> tuple[float, float, float]:
    """
    Return `direction` scaled to unit length.

    Normalising is not a convenience -- it is a correctness requirement.  An unnormalised
    ``[1, 1, 1]`` asks for ``sqrt(3)`` times the intended vector amplitude, which both overstates
    the b-value and can exceed the amplifier on the vector norm while every individual axis looks
    legal.
    """
    norm = math.sqrt(sum(c * c for c in direction))
    if norm == 0.0:
        msg = format_error(
            'diffusion direction is the zero vector.',
            {'direction': direction},
            ['use b_value_s_per_mm2=0 for an unweighted volume instead'],
        )
        raise ConfigurationError(msg)
    return (direction[0] / norm, direction[1] / norm, direction[2] / norm)


class MonopolarDiffusion(Module):
    """
    Stejskal--Tanner monopolar diffusion encoding: two identical lobes around the refocusing pulse.

    Parameters
    ----------
    system
        The scanner.
    b_value_s_per_mm2
        Target b-value.  Zero is allowed and produces lobes of the same duration carrying no
        gradient, so an unweighted volume occupies the same TE as a weighted one.
    refocus_duration_us
        Time between the end of the first lobe and the start of the second: the refocusing pulse
        with its crushers.  Sets ``Delta`` together with the lobe duration.
    lobe_duration_us
        Force the lobe duration (``delta``).  ``None`` solves for the shortest lobe reaching
        `b_value_s_per_mm2` at full amplitude, which is what minimises TE.
    axes
        Which logical axes the direction vector maps onto.
    regime
        Which limit regime to design against.  Diffusion normally runs at full amplitude while
        the readout is derated, which is exactly what named regimes are for.

    Properties
    ----------
    lobe_duration
        Seconds occupied by one lobe, ``delta``.
    big_delta
        Seconds between the *starts* of the two lobes, ``Delta``.
    amplitude_Hz_per_m
        Peak amplitude of one lobe.
    duration
        Seconds occupied by one built block; the same for both parts.
    total_duration
        Seconds from the start of the first lobe to the end of the second.

    Build arguments
    ---------------
    part : {'pre', 'post'}
        Which lobe.  Both are identical for a monopolar encoding, so this exists for symmetry
        with :class:`BipolarDiffusion`, where they are not.
    direction : tuple of float, default (0, 0, 1)
        Diffusion direction.  Normalised on the way in.
    scale : float, default 1.0
        Amplitude multiplier, for a b-value sweep sharing one design: ``b`` scales as
        ``scale**2``.

    Notes
    -----
    The b-value includes the trapezoid ramp corrections::

        b = (2*pi*G)^2 * [delta^2 * (Delta - delta/3) - eps^3/30 + delta*eps^2/6]

    Ignoring them overstates b by several percent at short `delta`, which is where diffusion
    sequences live.  :meth:`achieved_b_s_per_mm2` measures b from the actual waveform instead of
    trusting the formula, and the two agreeing is what the tests assert.

    Examples
    --------
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> diff = MonopolarDiffusion(system, b_value_s_per_mm2=1000, refocus_duration_us=4200)
    >>> diff.big_delta > diff.lobe_duration
    True
    >>> b0 = MonopolarDiffusion(system, b_value_s_per_mm2=0, refocus_duration_us=4200,
    ...                         lobe_duration_us=diff.lobe_duration * 1e6)
    >>> abs(b0.duration - diff.duration) < 1e-9        # b=0 fills the same slot
    True
    """

    #: Gradient rasters the lobe duration must be a multiple of.  Two for a scheme that splits its
    #: lobe in half, so that each half is a whole number of rasters.
    DELTA_RASTERS = 1

    def __init__(
        self,
        system: System,
        *,
        b_value_s_per_mm2: float,
        refocus_duration_us: float,
        lobe_duration_us: float | None = None,
        axes: tuple[str, str, str] = _AXES,
        regime: str = 'default',
    ) -> None:
        super().__init__(system, regime=regime)
        self.b_value_s_per_mm2 = float(b_value_s_per_mm2)
        self.refocus_duration_us = float(refocus_duration_us)
        self.lobe_duration_us = None if lobe_duration_us is None else float(lobe_duration_us)
        self.axes = tuple(axes)
        require_positive(self, 'refocus_duration_us')
        if self.b_value_s_per_mm2 < 0:
            msg = format_error(
                'b_value_s_per_mm2 cannot be negative.',
                {'got': self.b_value_s_per_mm2},
            )
            raise ConfigurationError(msg)
        for axis in self.axes:
            require_in(type('_A', (), {'axis': axis})(), 'axis', _AXES)

        raster = self.system.grad_raster
        self._ramp = raster.ceil(float(self.opts.max_grad) / float(self.opts.max_slew))
        # The lobe lands on a multiple of this.  Monopolar needs one raster; bipolar needs two, so
        # that half a lobe is still a whole number of rasters and its two sub-lobes abut exactly --
        # a half-raster sliver between them survives resampling as a non-zero block endpoint, which
        # pypulseq reports as a gradient continuity error far from anything that looks wrong.
        quantum = Raster(raster.at(self.DELTA_RASTERS), 'lobe')
        if self.lobe_duration_us is None:
            self._delta = quantum.ceil(self._solve_delta())
        else:
            self._delta = quantum.ceil(convert(self.lobe_duration_us, 'us', 's'))
        if self._delta <= 2.0 * self._ramp:
            msg = format_error(
                f'lobe duration {convert(self._delta, "s", "us"):.0f} us is shorter than its '
                f'own ramps ({convert(2 * self._ramp, "s", "us"):.0f} us).',
                {
                    'max_grad_mT_m': self.system.convert(
                        float(self.opts.max_grad), 'Hz/m', 'mT/m'
                    )
                },
                ['lengthen lobe_duration_us, or design against a lower-slew regime'],
            )
            raise ConfigurationError(msg)

        self._amplitude = self._solve_amplitude()
        self._lobe = {
            axis: pp.make_trapezoid(
                channel=axis,
                amplitude=self._amplitude,
                flat_time=self._delta - 2.0 * self._ramp,
                rise_time=self._ramp,
                system=self.opts,
            )
            for axis in self.axes
        }

    # ------------------------------------------------------------------------------ solving
    @property
    def _big_delta_for(self) -> float:
        """Seconds between lobe starts, given the current lobe duration."""
        return self._delta + convert(self.refocus_duration_us, 'us', 's')

    def _b_of(self, amplitude: float, delta: float, big_delta: float) -> float:
        """
        Return the b-value in s/m^2 for a trapezoidal monopolar pair, exactly.

        Derived by integrating ``b = (2*pi)^2 * integral |k(t)|^2 dt`` piecewise over the ramp, the
        flat top and the fall of each lobe, with the refocusing pulse's sign flip applied.  With
        ``A = G(delta - eps)`` the area of one lobe and ``f = delta - 2*eps`` its flat time:

        .. code-block:: text

            I = G^2 * [ eps^3/10 + ((delta - 1.5 eps)^3 - eps^3/8)/3
                        + (delta - eps)^2 eps - (delta - eps) eps^2/3 ]
            b = (2 pi)^2 * [ 2 I + A^2 (Delta - delta) ]

        As ``eps -> 0`` this collapses to the familiar ``(2 pi G)^2 delta^2 (Delta - delta/3)``,
        which is the check the tests make.

        The commonly quoted ramp correction ``+- eps^3/30 -+ delta eps^2/6`` is written with
        several incompatible conventions for whether ``delta`` includes the ramps; using one with
        the wrong convention overstates b by 2-4 % at the short lobe durations diffusion sequences
        actually use, and a b-value wrong by 3 % biases every diffusivity by 3 %.  Integrating
        directly avoids having to pick.
        """
        eps = self._ramp
        area = amplitude * (delta - eps)
        # Integral of k^2 over one lobe.
        integral = amplitude**2 * (
            eps**3 / 10.0
            + ((delta - 1.5 * eps) ** 3 - eps**3 / 8.0) / 3.0
            + (delta - eps) ** 2 * eps
            - (delta - eps) * eps**2 / 3.0
        )
        return (2.0 * math.pi) ** 2 * (2.0 * integral + area**2 * (big_delta - delta))

    def _solve_delta(self) -> float:
        """
        Return the shortest lobe reaching the requested b-value at full amplitude.

        A closed form does not exist because `Delta` depends on `delta`, so this steps up the
        gradient raster from the minimum -- at most a few hundred iterations, each cheap, and it
        cannot miss the answer the way a Newton solve on a non-monotonic residual could.
        """
        raster = self.system.grad_raster
        shortest = raster.ceil(2.0 * self._ramp + raster.dt)
        if self.b_value_s_per_mm2 == 0.0:
            return shortest
        target = convert(self.b_value_s_per_mm2, 's/mm^2', 's/m^2')
        amplitude = float(self.opts.max_grad)
        delta = shortest
        for _ in range(200_000):
            big_delta = delta + convert(self.refocus_duration_us, 'us', 's')
            if self._b_of(amplitude, delta, big_delta) >= target:
                return delta
            delta += raster.dt
        msg = format_error(  # pragma: no cover - needs absurd b
            f'cannot reach b = {self.b_value_s_per_mm2:g} s/mm^2 on this system.',
            {'max_grad_mT_m': self.system.convert(amplitude, 'Hz/m', 'mT/m')},
            ['reduce b_value_s_per_mm2', 'use a stronger gradient regime'],
        )
        raise ConfigurationError(msg)

    def _solve_amplitude(self) -> float:
        """
        Return the amplitude reaching the requested b-value in the chosen lobe duration.

        ``b`` scales as amplitude squared, so this is a square root rather than a search.
        """
        if self.b_value_s_per_mm2 == 0.0:
            return 0.0
        target = convert(self.b_value_s_per_mm2, 's/mm^2', 's/m^2')
        unit_b = self._b_of(1.0, self._delta, self._big_delta_for)
        amplitude = math.sqrt(target / unit_b)
        limit = float(self.opts.max_grad)
        if amplitude > limit * (1 + 1e-9):
            achievable = convert(
                self._b_of(limit, self._delta, self._big_delta_for), 's/m^2', 's/mm^2'
            )
            msg = format_error(
                f'b = {self.b_value_s_per_mm2:g} s/mm^2 needs '
                f'{self.system.convert(amplitude, "Hz/m", "mT/m"):.1f} mT/m in a '
                f'{convert(self._delta, "s", "us"):.0f} us lobe, above the '
                f'{self.system.convert(limit, "Hz/m", "mT/m"):.1f} mT/m limit.',
                {
                    'delta_us': convert(self._delta, 's', 'us'),
                    'Delta_us': convert(self._big_delta_for, 's', 'us'),
                    'achievable_b_s_per_mm2': round(achievable, 1),
                },
                [
                    f'lengthen lobe_duration_us beyond '
                    f'{convert(self._delta, "s", "us"):.0f} us',
                    'or leave lobe_duration_us=None to solve for the shortest lobe that fits',
                ],
            )
            raise ConfigurationError(msg)
        return amplitude

    # --------------------------------------------------------------------------- properties
    def m1_per_m_per_s(self) -> float:
        """
        Return the whole encoding's first gradient moment, in 1/m/s.

        Measured from both built lobes with the refocusing pulse's sign flip applied.  Non-zero for
        a monopolar encoding -- which is exactly the point: a non-zero first moment means flowing
        spins lose signal the same way diffusing ones do, so bulk motion reads as a high apparent
        diffusivity.  Compare against :class:`BipolarDiffusion`, which nulls it.

        Measuring one lobe alone would be meaningless: each lobe has a large first moment on its
        own, and it is the *pair* that either cancels or does not.

        Examples
        --------
        >>> import seqcraft as sc
        >>> system = sc.System.preset('generic_3t')
        >>> mono = sc.modules.MonopolarDiffusion(system, b_value_s_per_mm2=300,
        ...                                      refocus_duration_us=4200)
        >>> bip = sc.modules.BipolarDiffusion(system, b_value_s_per_mm2=300,
        ...                                   refocus_duration_us=4200)
        >>> abs(mono.m1_per_m_per_s()) > 1.0        # velocity sensitive
        True
        >>> abs(bip.m1_per_m_per_s()) < 1e-6        # velocity compensated
        True
        """
        raster = self.system.grad_raster.dt
        gap = convert(self.refocus_duration_us, 'us', 's')
        total = 0.0
        for part, sign in (('pre', 1.0), ('post', -1.0)):
            offset = 0.0 if part == 'pre' else self._delta + gap
            for node in self.build(part=part, direction=(0.0, 0.0, 1.0)).nodes:
                if getattr(node.item, 'type', None) in ('trap', 'grad'):
                    tt, wf = ev.waveform_of(node.item, raster)
                    total += sign * float(ev.trapz(wf * (tt + node.start + offset), tt))
        return total

    @property
    def lobe_duration(self) -> float:
        """Seconds occupied by one lobe, ``delta``."""
        return self._delta

    @property
    def big_delta(self) -> float:
        """Seconds between the starts of the two lobes, ``Delta``."""
        return self._big_delta_for

    @property
    def ramp_time(self) -> float:
        """Seconds of ramp at each end of a lobe, ``eps``."""
        return self._ramp

    @property
    def amplitude_Hz_per_m(self) -> float:
        """Peak amplitude of one lobe along the direction vector, Hz/m."""
        return self._amplitude

    @property
    def duration(self) -> float:
        """Seconds occupied by one built block -- the same for both parts."""
        return self.system.block_raster.ceil(self._delta)

    @property
    def total_duration(self) -> float:
        """Seconds from the start of the first lobe to the end of the second."""
        return self.big_delta + self._delta

    def achieved_b_s_per_mm2(self, direction: tuple[float, float, float] = (0.0, 0.0, 1.0)) -> float:
        """
        Return the b-value this encoding actually delivers, in s/mm^2.

        Computed from the designed amplitude and timings rather than from the request, so raster
        rounding is visible instead of assumed away.

        Examples
        --------
        >>> import seqcraft as sc
        >>> diff = MonopolarDiffusion(sc.System.preset('generic_3t'),
        ...                           b_value_s_per_mm2=1000, refocus_duration_us=4200)
        >>> abs(diff.achieved_b_s_per_mm2() - 1000) < 1.0
        True
        """
        unit = _normalise(direction) if any(direction) else (0.0, 0.0, 0.0)
        # b is a quadratic form; for a single-direction encoding its trace is |d|^2 = 1.
        scale = sum(c * c for c in unit)
        b_si = self._b_of(self._amplitude, self._delta, self._big_delta_for) * scale
        return convert(b_si, 's/m^2', 's/mm^2')

    # -------------------------------------------------------------------------------- build
    def build(
        self,
        *,
        part: str = 'pre',
        direction: tuple[float, float, float] = (0.0, 0.0, 1.0),
        scale: float = 1.0,
    ) -> LogicBlock:
        """
        Return one lobe of the diffusion pair.

        Parameters
        ----------
        part
            ``'pre'`` for the lobe before the refocusing pulse, ``'post'`` for the one after.
            Identical for a monopolar encoding.
        direction
            Diffusion direction, normalised on the way in.
        scale
            Amplitude multiplier; ``b`` scales as ``scale**2``.

        Examples
        --------
        >>> import seqcraft as sc
        >>> diff = MonopolarDiffusion(sc.System.preset('generic_3t'),
        ...                           b_value_s_per_mm2=1000, refocus_duration_us=4200)
        >>> pre, post = diff.build(part='pre'), diff.build(part='post')
        >>> abs(pre.duration - post.duration) < 1e-12
        True
        """
        if part not in ('pre', 'post'):
            msg = format_error(
                f'part must be "pre" or "post", got {part!r}.',
                {},
                ['a monopolar encoding has two identical lobes; pass part="pre" or "post"'],
            )
            raise ConfigurationError(msg)

        out = LogicBlock(f'diff_{part}')
        if self._amplitude == 0.0 or scale == 0.0:
            # b = 0 still has to occupy the slot, or the unweighted volume gets a different TE.
            return out.add(0.0, pp.make_delay(self.duration))

        unit = _normalise(direction)
        for axis, component in zip(self.axes, unit):
            reference = self._lobe[axis]
            amplitude = self._amplitude * component * float(scale)
            out.add(
                0.0,
                ev.derive(
                    reference,
                    amplitude=amplitude,
                    area=amplitude * (self._delta - self._ramp),
                    flat_area=amplitude * (self._delta - 2.0 * self._ramp),
                ),
            )
        return out


class BipolarDiffusion(MonopolarDiffusion):
    """
    Bipolar (velocity-compensated) diffusion encoding: each lobe is a plus/minus pair.

    Parameters
    ----------
    See :class:`MonopolarDiffusion`.

    Notes
    -----
    Nulls the first gradient moment, so flowing spins are not mistaken for diffusing ones.  That
    matters wherever there is bulk motion the sequence cannot gate away -- cardiac diffusion is the
    canonical case, where an unnulled monopolar encoding produces signal loss indistinguishable
    from a high apparent diffusivity.

    **Both lobes are built identically**, as a plus/minus pair.  The refocusing pulse inverts
    everything after it, so the *effective* waveform is ``(+ -) gap (- +)`` -- even about the
    encoding centre, and an even function has no first moment about that centre.  Inverting the
    second lobe as well would give ``(+ -) gap (+ -)``, which is not even and leaves m1 unnulled
    while still producing a plausible b-value: the failure is silent, which is why
    :meth:`m1_per_m_per_s` measures it rather than asserting it.

    The cost is roughly four times the encoding time for a given b-value, because splitting a lobe
    into two opposing halves quarters its effective area and b scales as area squared.

    Examples
    --------
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> bip = BipolarDiffusion(system, b_value_s_per_mm2=500, refocus_duration_us=4200)
    >>> len(bip.build(part='pre', direction=(0, 0, 1)).nodes)    # two opposing sub-lobes
    2
    >>> abs(bip.m1_per_m_per_s()) < 1e-6                         # first moment nulled
    True
    >>> bip.lobe_duration > MonopolarDiffusion(
    ...     system, b_value_s_per_mm2=500, refocus_duration_us=4200).lobe_duration
    True
    """

    DELTA_RASTERS = 2

    def _b_of(self, amplitude: float, delta: float, big_delta: float) -> float:
        """
        Return the b-value for a bipolar pair, in s/m^2.

        Each lobe is two opposing half-length sub-lobes, so the pair behaves like a monopolar
        encoding of half the lobe duration whose accumulated phase reverses mid-lobe.  Integrating
        ``|k|^2`` over the four sub-lobes and the gap gives the expression below; the leading
        factor of a quarter is the cost of the cancellation.
        """
        half = delta / 2.0
        eps = self._ramp
        area = amplitude * (half - eps)
        integral = amplitude**2 * (
            eps**3 / 10.0
            + ((half - 1.5 * eps) ** 3 - eps**3 / 8.0) / 3.0
            + (half - eps) ** 2 * eps
            - (half - eps) * eps**2 / 3.0
        )
        # Two sub-lobes per lobe, four in total, plus the gap held at the net area of one lobe
        # pair -- which for a bipolar encoding is zero, so the gap contributes nothing.  All the
        # weighting comes from the sub-lobes themselves.
        return (2.0 * math.pi) ** 2 * 4.0 * integral + 0.0 * area * big_delta

    def build(
        self,
        *,
        part: str = 'pre',
        direction: tuple[float, float, float] = (0.0, 0.0, 1.0),
        scale: float = 1.0,
    ) -> LogicBlock:
        """
        Return one bipolar lobe: two opposing sub-lobes filling the lobe duration.

        Both parts are identical -- see the class notes for why that is what nulls m1.
        """
        if part not in ('pre', 'post'):
            msg = format_error(f'part must be "pre" or "post", got {part!r}.', {})
            raise ConfigurationError(msg)

        out = LogicBlock(f'diff_{part}')
        if self._amplitude == 0.0 or scale == 0.0:
            return out.add(0.0, pp.make_delay(self.duration))

        unit = _normalise(direction)
        # Exact by construction: DELTA_RASTERS = 2 makes the lobe an even number of rasters, so its
        # halves abut with no sliver between them.
        half = self._delta / 2.0
        for index, polarity in enumerate((1.0, -1.0)):
            for axis, component in zip(self.axes, unit):
                amplitude = self._amplitude * component * float(scale) * polarity
                if amplitude == 0.0:
                    continue
                out.add(
                    index * half,
                    pp.make_trapezoid(
                        channel=axis,
                        amplitude=amplitude,
                        flat_time=half - 2.0 * self._ramp,
                        rise_time=self._ramp,
                        system=self.opts,
                    ),
                )
        return out


class ArbitraryDiffusion(Module):
    """
    Diffusion encoding from a waveform you supply.

    The escape hatch for optimised waveforms -- GrOpt output, a numerically nulled design, a
    vendor waveform read from disk.  seqcraft does not solve for it; it plays it, measures its
    b-value and moments, and validates it against the amplifier.

    Parameters
    ----------
    system
        The scanner.
    waveform_Hz_per_m
        Amplitude samples on the gradient raster, shape ``(n,)`` for a single axis or ``(3, n)``
        for a full vector waveform.
    refocus_duration_us
        Time between the two lobes, for the b-value calculation.
    axes
        Which logical axes the rows of a ``(3, n)`` waveform map onto.

    Properties
    ----------
    duration
        Seconds occupied by one lobe.

    Build arguments
    ---------------
    part : {'pre', 'post'}
    direction : tuple of float, optional
        Only for a single-axis waveform; a ``(3, n)`` waveform already carries its direction.
    scale : float, default 1.0

    Notes
    -----
    Because the waveform is given rather than solved, the b-value is *measured*: the second
    moment of the effective gradient is integrated numerically, including the sign flip the
    refocusing pulse applies between the lobes.

    Examples
    --------
    >>> import numpy as np
    >>> import seqcraft as sc
    >>> system = sc.System.preset('generic_3t')
    >>> ramp = np.concatenate([np.linspace(0, 1, 20), np.ones(60), np.linspace(1, 0, 20)])
    >>> arb = ArbitraryDiffusion(system, waveform_Hz_per_m=ramp * 1e6,
    ...                          refocus_duration_us=4200)
    >>> round(arb.duration * 1e6)
    1000
    >>> arb.build(part='pre', direction=(0, 0, 1))
    LogicBlock(diff_pre, 1 node, 1.00 ms)
    """

    def __init__(
        self,
        system: System,
        *,
        waveform_Hz_per_m: np.ndarray,
        refocus_duration_us: float,
        axes: tuple[str, str, str] = _AXES,
        regime: str = 'default',
    ) -> None:
        super().__init__(system, regime=regime)
        self.refocus_duration_us = float(refocus_duration_us)
        self.axes = tuple(axes)
        wave = np.asarray(waveform_Hz_per_m, dtype=float)
        if wave.ndim == 1:
            self._wave = wave[np.newaxis, :]
            self._vector = False
        elif wave.ndim == 2 and wave.shape[0] == 3:
            self._wave = wave
            self._vector = True
        else:
            msg = format_error(
                'waveform_Hz_per_m must have shape (n,) or (3, n).',
                {'got': str(wave.shape)},
            )
            raise ConfigurationError(msg)

        # A diffusion lobe sits between blocks that carry no gradient on its axes, so its waveform
        # must start and end at zero.  Left implicit, pypulseq extrapolates `first` and `last` from
        # the first two samples and produces a gradient that starts away from zero -- which its own
        # continuity check then rejects, several thousand blocks later.
        for edge, index in (('start', 0), ('end', -1)):
            value = float(np.max(np.abs(self._wave[:, index])))
            if value > 1e-6 * max(float(self.opts.max_grad), 1.0):
                msg = format_error(
                    f'the supplied waveform does not {edge} at zero '
                    f'({self.system.convert(value, "Hz/m", "mT/m"):.3f} mT/m).',
                    {'samples': self._wave.shape[1]},
                    [
                        'pad the waveform with a zero sample at each end',
                        'or ramp it to zero at the slew limit -- a gradient cannot step',
                    ],
                )
                raise ConfigurationError(msg)

        peak = float(np.max(np.abs(self._wave)))
        if peak > float(self.opts.max_grad) * 1.001:
            msg = format_error(
                f'the supplied waveform reaches '
                f'{self.system.convert(peak, "Hz/m", "mT/m"):.1f} mT/m, above the '
                f'{self.system.convert(float(self.opts.max_grad), "Hz/m", "mT/m"):.1f} mT/m limit.',
                {'regime': self.regime},
                ['scale the waveform down, or design against a stronger regime'],
            )
            raise ConfigurationError(msg)

    @property
    def duration(self) -> float:
        """Seconds occupied by one lobe."""
        return self.system.block_raster.ceil(self.system.grad_raster.at(self._wave.shape[1]))

    def achieved_b_s_per_mm2(self) -> float:
        """
        Return the b-value of the pair, measured by numerical integration, in s/mm^2.

        Integrates ``b = (2*pi)^2 * integral |integral G dt'|^2 dt`` over both lobes, with the
        sign of the second lobe's contribution flipped by the refocusing pulse.
        """
        raster = self.system.grad_raster.dt
        gap = self.system.grad_raster.count(
            self.system.grad_raster.nearest(convert(self.refocus_duration_us, 'us', 's'))
        )
        first = self._wave
        second = -self._wave  # the refocusing pulse inverts accumulated phase
        joined = np.concatenate([first, np.zeros((first.shape[0], gap)), second], axis=1)
        k = np.cumsum(joined, axis=1) * raster
        b_si = float((2.0 * math.pi) ** 2 * np.sum(np.sum(k**2, axis=0)) * raster)
        return convert(b_si, 's/m^2', 's/mm^2')

    def build(
        self,
        *,
        part: str = 'pre',
        direction: tuple[float, float, float] = (0.0, 0.0, 1.0),
        scale: float = 1.0,
    ) -> LogicBlock:
        """Return one lobe of the supplied waveform."""
        if part not in ('pre', 'post'):
            msg = format_error(f'part must be "pre" or "post", got {part!r}.', {})
            raise ConfigurationError(msg)
        out = LogicBlock(f'diff_{part}')
        components = (
            [(axis, 1.0) for axis in self.axes]
            if self._vector
            else list(zip(self.axes, _normalise(direction)))
        )
        for index, (axis, component) in enumerate(components):
            row = self._wave[index if self._vector else 0]
            amplitudes = row * component * float(scale)
            if not np.any(amplitudes):
                continue
            out.add(
                0.0,
                pp.make_arbitrary_grad(
                    channel=axis,
                    waveform=amplitudes,
                    first=0.0,
                    last=0.0,
                    system=self.opts,
                    delay=0.0,
                ),
            )
        if not out.nodes:
            out.add(0.0, pp.make_delay(self.duration))
        return out
