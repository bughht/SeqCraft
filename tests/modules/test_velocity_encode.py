"""
``VelocityEncode`` -- every claim integrated off the emitted waveform, never asked of the module.

The organising rule is the candidate record's rule E: the validator integrates the events, it does
not ask the module what it thinks it built.  :func:`seqcraft.moments` is that integrator, and it
takes no ``Opts`` and does no compile -- a moment is a property of a waveform, so a test that
reads one back is reading the same thing a scanner would play.

The second rule is that the **pair** is the subject.  ``m0 = 0`` is a claim about one toggle, but
``delta_m1`` is a claim about the difference between two, and a module that got both toggles wrong
in the same direction would satisfy every per-toggle check.  So the difference is measured
directly wherever it is the thing being claimed.
"""

from __future__ import annotations

import copy

import numpy as np
import pytest

import seqcraft as sc
from seqcraft.design.logic import flatten

#: A spread of VENCs that crosses the triangular/trapezoidal boundary in both directions.
VENCS = (0.05, 0.2, 0.5, 1.5, 5.0, 50.0)

#: Moments are areas, so "zero" is relative to the lobe area rather than to 1.0.
RELATIVE = 1e-9


def encode(opts, **overrides):
    kwargs = {'opts': opts, 'venc_m_s': 1.5, 'axis': 'z'}
    kwargs.update(overrides)
    return sc.modules.VelocityEncode(**kwargs)


def _moment(module, order, polarity=1, at=0.0):
    """One toggle's moment on `module`'s axis, integrated off the emitted events."""
    tree = sc.LogicBlock('probe').add(at, module(polarity=polarity))
    return float(sc.moments(tree, order).get(module.axis, 0.0))


def _delta_m1(module):
    """The change in first moment between the toggles -- the quantity venc is defined by."""
    return _moment(module, 1, polarity=1) - _moment(module, 1, polarity=-1)


# ------------------------------------------------------------------- the two invariants
@pytest.mark.parametrize('venc_m_s', VENCS)
@pytest.mark.parametrize('polarity', (1, -1))
def test_stationary_spins_see_nothing(opts, venc_m_s: float, polarity: int) -> None:
    """
    **Half the contract.** Net area zero, or the encoding displaces stationary spins too.

    A bipolar pair with unequal lobes still encodes velocity, so this failure does not look like a
    failure -- it looks like a geometric shift in the image, which is the record's hazard (2).
    """
    module = encode(opts, venc_m_s=venc_m_s)

    assert _moment(module, 0, polarity) == pytest.approx(
        0.0, abs=module.lobe_area_1_per_m * RELATIVE)


@pytest.mark.parametrize('venc_m_s', VENCS)
def test_the_toggles_differ_by_the_first_moment_venc_asks_for(opts, venc_m_s: float) -> None:
    """
    **The other half, and the claim the module exists to make good on.**

    `venc_m_s` is defined by the *change* in first moment across the toggles, and that change is
    measured here rather than the design being trusted for it.  In SeqCraft's Hz/m gradients the
    relation carries no gamma at all: ``delta_m1 == 1 / (2 * venc)``.
    """
    module = encode(opts, venc_m_s=venc_m_s)

    assert _delta_m1(module) == pytest.approx(1.0 / (2.0 * venc_m_s), rel=1e-12)
    assert module.delta_m1_s_per_m == pytest.approx(_delta_m1(module), rel=1e-12)


@pytest.mark.parametrize('venc_m_s', VENCS)
def test_a_spin_at_venc_accumulates_exactly_pi_between_the_toggles(opts, venc_m_s: float) -> None:
    """
    The definition of VENC, evaluated on the waveform rather than on the formula that made it.

    A spin at ``x(t) = x0 + v t`` accumulates ``phi = 2 pi (m0 x0 + m1 v)``.  Both moments come
    from the emitted events, so this is the phase the sequence delivers -- and with `m0` measured
    rather than assumed, the `x0` term is carried through instead of being waved away.
    """
    module = encode(opts, venc_m_s=venc_m_s)
    x0_m = 0.037

    phases = [2.0 * np.pi * (_moment(module, 0, p) * x0_m
                             + _moment(module, 1, p) * venc_m_s) for p in (1, -1)]
    assert phases[0] - phases[1] == pytest.approx(np.pi, rel=1e-9)


# ------------------------------------------------------------------------ the factor of two
@pytest.mark.parametrize('venc_m_s', VENCS)
def test_one_toggle_carries_half_the_change(opts, venc_m_s: float) -> None:
    """
    **The mistake the public API exists to prevent**, stated as a test.

    A caller who gives one toggle a first moment of ``delta_m1`` rather than ``delta_m1 / 2`` is
    out by a factor of two, and gets a plausible velocity map at half the VENC they asked for.
    The module never takes a moment, so the only way to check the split is to measure it.
    """
    module = encode(opts, venc_m_s=venc_m_s)

    for polarity in (1, -1):
        assert _moment(module, 1, polarity) == pytest.approx(
            polarity * module.delta_m1_s_per_m / 2.0, rel=1e-12)
        assert module.first_moment_s_per_m(polarity) == pytest.approx(
            _moment(module, 1, polarity), rel=1e-12)


def test_polarity_is_the_sign_of_the_emitted_first_moment(opts) -> None:
    """
    `polarity` is defined on the waveform, which is what makes it checkable without a convention.

    Not the sign of a velocity and not the sign of a reconstructed phase -- those belong to the
    reconstruction, and nothing here has to agree with it to assert this.
    """
    module = encode(opts)

    assert _moment(module, 1, +1) > 0.0
    assert _moment(module, 1, -1) < 0.0
    # +1 plays the negative lobe first, and that ordering is the whole of the sign.
    first = [event for _, event, _ in sorted(flatten(module(polarity=+1)))][0]
    assert float(first.amplitude) < 0.0


# --------------------------------------------------------------- what else the toggles share
@pytest.mark.parametrize('venc_m_s', VENCS)
def test_the_toggles_are_one_waveform_at_opposite_polarity(opts, venc_m_s: float) -> None:
    """
    The two toggles are the same gradient played with the sign flipped, and **every** moment
    negates with it -- not only the first.

    Worth stating exactly, because it is tempting to write that the toggles "differ in the first
    moment and in nothing else".  They do not: negating a waveform negates its second moment too,
    so the honest invariant is that one is the other's negation, with ``m0 = 0`` making the
    zeroth the only moment that is equal.
    """
    module = encode(opts, venc_m_s=venc_m_s)
    scale = module.lobe_area_1_per_m * module.duration_s ** 2

    assert module(polarity=1).duration == module(polarity=-1).duration
    for order in (0, 1, 2, 3):
        plus, minus = _moment(module, order, 1), _moment(module, order, -1)
        assert plus == pytest.approx(-minus, abs=scale * RELATIVE)
    assert _moment(module, 2, 1) != pytest.approx(_moment(module, 2, -1), abs=scale * 1e-6)


@pytest.mark.parametrize('venc_m_s', VENCS)
def test_the_first_moment_does_not_depend_on_where_the_pair_is_placed(opts,
                                                                     venc_m_s: float) -> None:
    """
    Once ``m0 = 0``, shifting the time origin changes ``m1`` by ``m0 * dt``, which is zero.

    This is the property that makes the module a leaf: it can be placed anywhere in a repetition
    without being told where the echo is, and still deliver the first moment it promised.  A
    module that had failed to null `m0` would fail *this* test with a placement-dependent answer,
    which is a second, independent route to the same defect.
    """
    module = encode(opts, venc_m_s=venc_m_s)

    early = _moment(module, 1, 1, at=0.0)
    late = _moment(module, 1, 1, at=17e-3)
    assert early == pytest.approx(late, rel=1e-9)


# ---------------------------------------------------------------------------- metamorphic
def test_halving_venc_doubles_the_change_and_leaves_the_area_at_zero(opts) -> None:
    """VENC is a scale, not a redesign: the relation holds across it and `m0` stays nulled."""
    slow, fast = encode(opts, venc_m_s=0.5), encode(opts, venc_m_s=0.25)

    assert _delta_m1(fast) == pytest.approx(2.0 * _delta_m1(slow), rel=1e-12)
    for module in (slow, fast):
        assert _moment(module, 0) == pytest.approx(
            0.0, abs=module.lobe_area_1_per_m * RELATIVE)


def test_a_larger_venc_is_a_shorter_pair(opts) -> None:
    """
    Less first moment is less gradient-time area, so a larger VENC is a shorter pair.

    A monotonic relation, and deliberately not a minimality claim: the design rounds the
    continuous solution onto the raster rather than searching the lattice, so it is legal and
    exact but not necessarily the shortest legal pair -- see
    `test_the_design_is_not_claimed_to_be_the_shortest_on_the_raster`.
    """
    durations = [encode(opts, venc_m_s=v).duration_s for v in VENCS]

    assert durations == sorted(durations, reverse=True)


def test_the_axis_is_the_callers_and_nothing_leaks_onto_the_others(opts) -> None:
    """Each logical axis is encoded independently, so a 4D-flow study is three of these."""
    for axis in ('x', 'y', 'z'):
        module = encode(opts, axis=axis)

        assert set(sc.moments(module(polarity=1), 1)) == {axis}


# ------------------------------------------------------------------------ hardware and raster
@pytest.mark.parametrize('venc_m_s', VENCS)
def test_the_pair_stays_inside_the_gradient_system(opts, venc_m_s: float) -> None:
    """
    The times are rounded **up** onto the raster and the amplitude is re-solved against them.

    That ordering is what keeps the design legal: growing the times can only lower the amplitude
    the target needs, so a rounding can never push the gradient or the slew over the limit.
    """
    module = encode(opts, venc_m_s=venc_m_s)
    raster = float(opts.grad_raster_time)

    assert module.amplitude_hz_per_m <= float(opts.max_grad)
    assert module.amplitude_hz_per_m / module.rise_time_s <= float(opts.max_slew)
    for interval in (module.rise_time_s, module.flat_time_s, module.duration_s):
        assert round(interval / raster, 9) == pytest.approx(round(interval / raster), abs=1e-9)


@pytest.mark.parametrize('polarity', (1, -1))
def test_the_compiler_accepts_a_toggle_on_its_own(opts, polarity: int) -> None:
    """The contract is the emitted file, so a toggle has to survive compilation."""
    seq = sc.compile(encode(opts)(polarity=polarity), opts)

    assert len(seq.block_events) >= 2


def test_a_weaker_gradient_system_gives_a_longer_pair_with_the_same_moment(opts) -> None:
    """
    The design adapts to the hardware and the physical target does not move.

    A claim about `venc_m_s` that only held on one gradient system would not be a contract.
    """
    weak = copy.copy(opts)
    weak.max_grad, weak.max_slew = float(opts.max_grad) / 4.0, float(opts.max_slew) / 4.0

    strong_pair, weak_pair = encode(opts, venc_m_s=0.5), encode(weak, venc_m_s=0.5)

    assert weak_pair.duration_s > strong_pair.duration_s
    assert _delta_m1(weak_pair) == pytest.approx(_delta_m1(strong_pair), rel=1e-12)
    assert weak_pair.amplitude_hz_per_m <= weak.max_grad


# ------------------------------------------------------------------------------ the refusals
@pytest.mark.parametrize('polarity', (0, 2, -2, 0.5, '+', None))
def test_a_polarity_that_is_not_a_toggle_is_refused(opts, polarity: object) -> None:
    """Two toggles is what a phase-difference reconstruction subtracts; there is no third."""
    with pytest.raises(sc.ConfigurationError, match='polarity') as caught:
        encode(opts)(polarity=polarity)

    assert 'first moment' in str(caught.value)


@pytest.mark.parametrize('venc_m_s', (0.0, -1.5))
def test_a_venc_that_is_not_a_speed_is_refused(opts, venc_m_s: float) -> None:
    with pytest.raises(sc.ConfigurationError, match='venc_m_s'):
        encode(opts, venc_m_s=venc_m_s)


def test_an_unknown_axis_is_refused(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='axis'):
        encode(opts, axis='q')


def test_polarity_has_no_default_and_must_be_chosen(opts) -> None:
    """
    **A toggle nobody selected is how a velocity map gets made from two copies of one encoding.**

    One instance owns `venc_m_s`, the design and the pair invariant; which toggle a given
    acquisition plays is the caller's, and the pair is only a pair because someone chose both.
    So `polarity` is a required keyword, exactly as `line` is for `PhaseEncode` -- and this is the
    regression guard against a default quietly coming back.
    """
    module = encode(opts)

    with pytest.raises(TypeError, match='polarity'):
        module()
    with pytest.raises(TypeError, match='polarity'):
        module.build()
    with pytest.raises(TypeError, match='polarity'):
        module.first_moment_s_per_m()

    assert module(polarity=+1).duration > 0.0, 'chosen explicitly, it still builds'


def test_the_design_is_not_claimed_to_be_the_shortest_on_the_raster(opts) -> None:
    """
    **The design rounds the continuous optimum up; it does not search the discrete lattice.**

    Written as a test because the distinction is easy to lose in a docstring, and because the gap
    is real rather than theoretical: at `venc_m_s = 1.0` on the fixture's gradient system this
    design is 1.100 ms, and a legal raster point exists at 1.080 ms.  What the module *does*
    guarantee is asserted alongside it -- exact first moment, inside the hardware limits.

    If a later change ever adds a lattice search, this test is the one that should fail and be
    deleted deliberately, rather than a claim silently becoming true by accident.
    """
    module = encode(opts, venc_m_s=1.0)
    raster = float(opts.grad_raster_time)
    target = module.delta_m1_s_per_m / 2.0

    shortest = min(
        2.0 * (2.0 * rise + flat)
        for rise, flat in ((n * raster, m * raster) for n in range(1, 60) for m in range(200))
        if (amplitude := target / ((flat + rise) * (2.0 * rise + flat))) <= float(opts.max_grad)
        and amplitude / rise <= float(opts.max_slew)
    )

    assert shortest < module.duration_s, 'the lattice admits a shorter legal pair'
    assert _delta_m1(module) == pytest.approx(1.0 / (2.0 * 1.0), rel=1e-12)
    assert module.amplitude_hz_per_m <= float(opts.max_grad)
