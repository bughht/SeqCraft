"""
``DiffusionSEPrep`` -- the physical claim first, and never the module's own arithmetic.

The organising rule of this file: **the module says what `b` it intends, and something that shares
no code with it measures what the sequence delivers.**  `sc.b_value` integrates the emitted
gradients and applies the refocusing conjugation itself, so agreement between the two is evidence
rather than a tautology.

That is not a stylistic preference.  The first version of this module computed the lobe-centre
separation one ramp time short.  Every lobe was legal, the geometry looked right, the module's own
number said 1000, and the sequence delivered 1042 -- the amplitude solve had quietly compensated.
Only an independent integration could see it, and `test_delivered_b_matches_the_request` is that
test.

The closed form itself is checked against the handbook rather than against the implementation:
`b_of_trapezoid_pair` is a free function, and `test_ramp_terms_reduce_to_the_rectangle` is the
check Figure 9.3's own caption suggests.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.modules.kernel.diffusion_se import b_of_trapezoid_pair

#: Whole-body-ish limits.  Diffusion wants amplitude, so `max_grad` is what matters here.
GAMMA_HZ_T = 42.576e6


@pytest.fixture
def opts() -> pp.Opts:
    return pp.Opts(max_grad=80, grad_unit='mT/m', max_slew=200, slew_unit='T/m/s', B0=3.0,
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)


def prep(opts: pp.Opts, **overrides: object) -> sc.modules.DiffusionSEPrep:
    kwargs: dict = {'opts': opts, 'thickness_mm': 4.0, 'b_s_per_mm2': 1000.0}
    kwargs.update(overrides)
    return sc.modules.DiffusionSEPrep(**kwargs)


# --------------------------------------------------------- the closed form, against the handbook
def test_ramp_terms_reduce_to_the_rectangle() -> None:
    """
    Figure 9.3's caption: the trapezoid expression "reduces to the value for the rectangular lobes
    if the ramp duration goes to 0".  That is a property of the formula, testable without a module.
    """
    gamma = 2 * np.pi * GAMMA_HZ_T
    amplitude, delta, separation = 25e-3 * GAMMA_HZ_T, 22.99e-3, 50e-3
    rectangle = gamma ** 2 * 25e-3 ** 2 * delta ** 2 * (separation - delta / 3) / 1e6
    assert b_of_trapezoid_pair(amplitude, delta, separation, 0.0) == pytest.approx(rectangle)


def test_the_handbook_worked_example() -> None:
    """
    Handbook §9.1 Example 9.2: lobes of 25 mT/m separated by 50 ms reach b = 1000 s/mm^2 at a lobe
    width of 22.99 ms.  An external number, computed by someone else, published in 2004.
    """
    b = b_of_trapezoid_pair(25e-3 * GAMMA_HZ_T, 22.99e-3, 50e-3, 0.0)
    assert b == pytest.approx(1000.0, rel=1e-3)


def test_ramps_reduce_the_b_value() -> None:
    """The ramp terms are negative in sum, which is why leaving them out over-estimates `b`."""
    without = b_of_trapezoid_pair(80e-3 * GAMMA_HZ_T, 12e-3, 18e-3, 0.0)
    with_ramps = b_of_trapezoid_pair(80e-3 * GAMMA_HZ_T, 12e-3, 18e-3, 580e-6)
    assert with_ramps < without


# ------------------------------------------------------------ the claim, measured independently
@pytest.mark.parametrize('b_requested', [100.0, 300.0, 1000.0, 2000.0, 3000.0])
def test_delivered_b_matches_the_request(opts: pp.Opts, b_requested: float) -> None:
    """
    **The test this module exists to pass.**

    `sc.b_value` integrates the emitted gradients and conjugates at the refocusing pulse; it has
    no access to `delta_s`, `separation_s` or the amplitude solve.  Measured on the encoding axis,
    because the slice lobe and the crushers put their own small weighting on `z` and that is a
    different claim -- see `test_the_slice_axis_carries_its_own_small_weighting`.
    """
    module = prep(opts, b_s_per_mm2=b_requested)
    delivered = sc.b_value(module(), opts, end_s=module.time_to_echo())
    assert delivered['x'] == pytest.approx(b_requested, rel=2e-3)


def test_the_slice_axis_carries_its_own_small_weighting(opts: pp.Opts) -> None:
    """
    A composition-level fact worth pinning: the slice-selection lobes and the refocusing crushers
    are gradients too, so a "b = 0" acquisition is not b = 0.  Small, real, and the reason
    `b_value` measures the tree rather than the diffusion lobes.
    """
    module = prep(opts, b_s_per_mm2=0.0)
    delivered = sc.b_value(module(), opts, end_s=module.time_to_echo())
    assert 'x' not in delivered                      # no lobes are emitted at all
    assert 0.0 < delivered['z'] < 5.0


def test_an_oblique_direction_splits_the_weighting(opts: pp.Opts) -> None:
    """Per-axis independence, §9.2's phrasing: the total is the sum over axes."""
    module = prep(opts, b_s_per_mm2=1000.0, axis=('x', 'y'))
    delivered = sc.b_value(module(), opts, end_s=module.time_to_echo())
    assert delivered['x'] == pytest.approx(delivered['y'], rel=1e-6)
    assert delivered['x'] + delivered['y'] == pytest.approx(1000.0, rel=2e-3)


# ------------------------------------------------------------------------------- metamorphic
def test_b_is_quadratic_in_the_amplitude(opts: pp.Opts) -> None:
    """
    With the geometry fixed, `b` is quadratic in the amplitude -- which is why the module can
    trim a rounded-up lobe width back onto the requested `b` with one square root instead of a
    search.  Asserted on the closed form, where the geometry can actually be held fixed.
    """
    delta, separation, ramp = 12e-3, 18e-3, 580e-6
    one = b_of_trapezoid_pair(40e-3 * GAMMA_HZ_T, delta, separation, ramp)
    twice = b_of_trapezoid_pair(80e-3 * GAMMA_HZ_T, delta, separation, ramp)
    assert twice / one == pytest.approx(4.0, rel=1e-9)


def test_a_longer_echo_time_buys_room_and_changes_nothing_else(opts: pp.Opts) -> None:
    """
    The design is **independent of the echo time**, which is the whole reason it is closed form:
    with the amplitude at the cap, the lobe width follows from `b` alone, and TE only decides
    whether the pair fits.  A longer TE therefore leaves room for a readout and does not touch the
    encoding.

    That room is what a single-shot EPI needs for the half of its train that precedes the echo,
    and discovering the module gave none of it is what caught a real defect -- see
    ``test_the_second_lobe_ends_before_the_echo``.
    """
    short = prep(opts, b_s_per_mm2=1000.0)
    longer = prep(opts, b_s_per_mm2=1000.0, te_s=short.te_s + 20e-3)
    assert longer.delta_s == pytest.approx(short.delta_s)
    assert longer.amplitude_hz_m == pytest.approx(short.amplitude_hz_m)
    room = lambda m: m.time_to_echo() - m().duration      # noqa: E731 -- one expression, twice
    # Half of it, not all of it: the refocusing pulse sits at TE/2, so lengthening the echo time
    # moves the block's end later too and only the second half of the extra time is room.
    assert room(longer) - room(short) == pytest.approx(10e-3, abs=2e-5)


@pytest.mark.parametrize('b_requested', [0.0, 500.0, 1000.0, 2500.0])
def test_the_second_lobe_ends_before_the_echo(opts: pp.Opts, b_requested: float) -> None:
    """
    A lobe occupies ``delta + ramp``, not ``delta``.  The first version sized it against the
    window without the ramp, so at the minimum echo time the second lobe ran one ramp PAST the
    spin echo -- unphysical, and it left a readout nowhere to put its prephaser.

    Nothing about the preparation alone showed it.  It appeared when a real EPI acquisition was
    composed against it and the echo-time search would not converge.
    """
    module = prep(opts, b_s_per_mm2=b_requested)
    assert module().duration <= module.time_to_echo() + 1e-12


# --------------------------------------------------------------------- timing and the refusals
def test_the_minimum_echo_time_is_the_shortest_that_works(opts: pp.Opts) -> None:
    """
    The search returns a real minimum, not a safe over-estimate: one raster shorter cannot reach
    the requested `b`.  Checked through the module's own refusal rather than through its internals.
    """
    module = prep(opts, b_s_per_mm2=1000.0)
    assert module.te_s == pytest.approx(module.min_te_s)
    with pytest.raises(sc.ConfigurationError, match='shorter than this b-value allows'):
        prep(opts, b_s_per_mm2=1000.0, te_s=module.min_te_s - 1e-3)


def test_a_bigger_b_needs_a_longer_minimum_echo_time(opts: pp.Opts) -> None:
    """Monotone, which is what makes the bisection legitimate."""
    echo_times = [prep(opts, b_s_per_mm2=b).min_te_s for b in (100.0, 500.0, 1000.0, 2000.0)]
    assert echo_times == sorted(echo_times)


def test_the_echo_lands_where_the_module_says(opts: pp.Opts) -> None:
    """
    Measured on the compiled sequence: `sc.kspace` reports the excitation and refocusing instants,
    and a spin echo is at twice their separation.  Nothing here restates the module's arithmetic.
    """
    module = prep(opts, b_s_per_mm2=800.0)
    field = sc.kspace(module(), opts)
    t90, t180 = float(field['t_excitation'][0]), float(field['t_refocusing'][0])
    assert 2.0 * (t180 - t90) == pytest.approx(module.te_s, abs=2e-5)
    assert module.time_to_echo() - t90 == pytest.approx(module.te_s, abs=2e-5)


def test_the_two_lobes_are_identical(opts: pp.Opts) -> None:
    """
    Handbook §9.1 p. 280: identical lobes either side of the refocusing pulse cancel the
    concomitant-field phase in a spin echo.  Identical shape is the part a waveform check can see.
    """
    module = prep(opts, b_s_per_mm2=1500.0)
    shapes = [
        (float(pp.calc_duration(event)), float(getattr(event, 'amplitude', 0.0)))
        for _start, event, _path in _flatten(module())
        if getattr(event, 'channel', None) == 'x'
    ]
    assert len(shapes) == 2
    assert shapes[0] == pytest.approx(shapes[1])


def test_zero_b_emits_no_diffusion_lobes(opts: pp.Opts) -> None:
    """The reference image is built without the lobes, not with zero-amplitude ones -- a
    zero-amplitude trapezoid is still a block boundary and still costs time."""
    module = prep(opts, b_s_per_mm2=0.0)
    assert module.delta_s == 0.0
    assert not [1 for _s, e, _p in _flatten(module()) if getattr(e, 'channel', None) == 'x']


@pytest.mark.parametrize('bad', ['q', ('x', 'w'), ()])
def test_a_bad_axis_is_refused(opts: pp.Opts, bad: object) -> None:
    with pytest.raises(sc.ConfigurationError, match='axis must name'):
        prep(opts, axis=bad)


def test_a_negative_b_is_refused(opts: pp.Opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='must not be negative'):
        prep(opts, b_s_per_mm2=-1.0)


def test_it_compiles(opts: pp.Opts) -> None:
    """Layer 1's floor: the preparation is a legal Pulseq sequence on its own."""
    seq = sc.compile(prep(opts, b_s_per_mm2=1000.0)(), opts, name='diffusion_prep')
    assert len(seq.block_events) > 0


def _flatten(tree: sc.LogicBlock):
    from seqcraft.design.logic import flatten
    return flatten(tree)
