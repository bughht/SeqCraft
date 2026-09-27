"""
Velocity compensation on ``CartesianLine`` -- moments summed over the emitted events.

The organising rule comes from the flow-moments record: **correctness is a property of the whole
waveform on the axis, not of a compensating lobe in isolation.**  A lobe that nulls its own first
moment while ignoring what the readout has already accumulated is self-consistent and not
compensated, so every assertion here integrates every gradient event on the axis, truncated at the
echo, and sums them.  Nothing asks the module what it thinks it built.

The second rule is the semantic time.  A first moment is only defined once an origin is named, and
the origin here is the **echo** -- the instant ``k = 0`` is sampled.  ``time_to_echo`` is where
that instant comes from, and it is the same number the rest of the package places TE with.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.design.events import knots_of, pwl_moment
from seqcraft.design.logic import flatten

#: Three protocols that put very different areas and durations through the same construction.
PROTOCOLS = ((64, 500.0), (128, 250.0), (32, 1000.0))

#: The readout geometries the option is offered for.  The solve reads the lobe's measured
#: pre-echo shape, so each of these reaches it differently: partial Fourier moves the echo off
#: centre, and a train's first lobe is the one the winder is solved against.
GEOMETRIES = {
    'full echo': {},
    'partial fourier': {'partial_fourier': 0.75},
    'monopolar train': {'echoes': 4, 'polarity': 'monopolar'},
    'bipolar train': {'echoes': 4, 'polarity': 'bipolar'},
}


def line(opts, **overrides):
    kwargs = {'opts': opts, 'fov_mm': 220.0, 'matrix': 64, 'bandwidth_hz_px': 500.0}
    kwargs.update(overrides)
    return sc.modules.CartesianLine(**kwargs)


def moment_to_echo(module, order: int, echo_index: int = 0) -> float:
    """
    The `order`-th moment of **every** gradient event on the module's axis, at an echo.

    Truncated at that echo's time and referenced to it, integrated from each event's own
    piecewise-linear knots.  Summed across events, which is the part no per-event check reaches:
    the prephaser lobes and the readout's own ramp-up are three separate events whose moments only
    cancel together.
    """
    echo = module.te_s[echo_index]
    total = 0.0
    for start, event, _ in flatten(module()):
        if getattr(event, 'type', None) not in ('trap', 'grad'):
            continue
        if getattr(event, 'channel', None) != module.axis:
            continue
        times, amps = knots_of(event, start)
        if times.size < 2 or echo <= times[0]:
            continue
        if echo < times[-1]:
            cut = int(np.searchsorted(times, echo))
            edge = float(np.interp(echo, times, amps))
            times, amps = np.append(times[:cut], echo), np.append(amps[:cut], edge)
        total += pwl_moment(times - echo, amps, order)
    return float(total)


def tolerance(module, order: int) -> float:
    """
    What "zero" means for an `order`-th moment on this readout.

    Dimensioned rather than absolute: ``m0`` is an area in 1/m, ``m1`` an area times a time in
    s/m, so the scale each is compared against carries the same units it does.  The readout's own
    area to the echo and its echo time are the two quantities of the construction, so the natural
    scale is ``|area| * echo**order``, and the tolerance is a small fraction of it.
    """
    return abs(module.area_to_echo_per_m) * module.time_to_echo() ** order * 1e-9


# ------------------------------------------------------------------------ the moment conditions
@pytest.mark.parametrize(('matrix', 'bandwidth'), PROTOCOLS)
def test_the_first_moment_is_null_at_the_echo(opts, matrix: int, bandwidth: float) -> None:
    """
    **The condition the option exists to meet**, measured on the complete emitted waveform.

    A spin at ``x(t) = x0 + v t`` accumulates ``phi = 2 pi (m0 x0 + m1 v)`` over a gradient
    waveform, so ``m1 = 0`` at the echo removes the constant-velocity phase term there, just as
    ``m0 = 0`` removes the constant-position one.  Whether that reduces an artefact in an image is
    a separate question this does not reach.
    """
    module = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth, _null_moment_order=1)

    assert moment_to_echo(module, 1) == pytest.approx(0.0, abs=tolerance(module, 1))


@pytest.mark.parametrize(('matrix', 'bandwidth'), PROTOCOLS)
def test_the_zeroth_moment_is_still_null_at_the_echo(opts, matrix: int,
                                                     bandwidth: float) -> None:
    """
    The echo condition is not traded away to get the first moment.

    ``m0 = 0`` at the echo is what puts ``k = 0`` there at all; a compensated readout that lost it
    would image a shifted k-space and still look plausible.
    """
    module = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth, _null_moment_order=1)

    assert moment_to_echo(module, 0) == pytest.approx(0.0, abs=tolerance(module, 0))


@pytest.mark.parametrize(('matrix', 'bandwidth'), PROTOCOLS)
def test_an_uncompensated_readout_has_a_first_moment_to_null(opts, matrix: int,
                                                             bandwidth: float) -> None:
    """
    The control: without the option the first moment is large, so the test above can see it.

    An assertion that something is zero proves nothing unless the same measurement is non-zero
    when the thing being measured is switched off.  The ordinary prephaser nulls ``m0`` and leaves
    ``m1`` at tens to hundreds of times the tolerance used above.
    """
    module = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth)

    assert abs(moment_to_echo(module, 1)) > tolerance(module, 1) * 1e6
    assert moment_to_echo(module, 0) == pytest.approx(0.0, abs=tolerance(module, 0))


# ------------------------------------------------------------------------------ the construction
@pytest.mark.parametrize(('matrix', 'bandwidth'), PROTOCOLS)
def test_the_two_lobes_sum_to_the_area_the_echo_needs(opts, matrix: int,
                                                      bandwidth: float) -> None:
    """
    Two lobes of opposite sign, and it is their **total** that cancels the readout's area.

    The individual areas are neither ``-area_to_echo_per_m`` nor half of it, which is why
    `prephaser_area_per_m` is documented as the total and the lobes are exposed separately.
    """
    module = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth, _null_moment_order=1)
    areas = [float(lobe.area) for lobe in module.prephaser_lobes]

    assert len(areas) == 2
    assert areas[0] * areas[1] < 0.0, 'opposite signs'
    assert sum(areas) == pytest.approx(-module.area_to_echo_per_m, rel=1e-9)
    assert module.prephaser_area_per_m == pytest.approx(sum(areas), rel=1e-12)


def test_the_ordinary_prephaser_is_one_lobe(opts) -> None:
    """``_null_moment_order=0`` is exactly today's readout, and the default."""
    module = line(opts)

    assert module._null_moment_order == 0
    assert len(module.prephaser_lobes) == 1
    assert module.prephaser_area_per_m == pytest.approx(-module.area_to_echo_per_m, rel=1e-12)


@pytest.mark.parametrize(('matrix', 'bandwidth'), PROTOCOLS)
def test_compensation_costs_time_and_the_cost_is_readable(opts, matrix: int,
                                                          bandwidth: float) -> None:
    """
    Three lobes before the echo instead of two, so the echo is later -- and it is reported.

    The record asks for the cost to be reported rather than hidden; `prephaser_duration_s` and
    `time_to_echo` are both longer, and both are the numbers a composing kernel already places
    against.
    """
    plain = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth)
    compensated = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth, _null_moment_order=1)

    assert compensated.prephaser_duration_s > plain.prephaser_duration_s
    assert compensated.time_to_echo() - plain.time_to_echo() == pytest.approx(
        compensated.prephaser_duration_s - plain.prephaser_duration_s, abs=1e-12)


# ------------------------------------------------------------------------ hardware and raster
@pytest.mark.parametrize(('matrix', 'bandwidth'), PROTOCOLS)
def test_the_compensated_prephaser_stays_inside_the_gradient_system(opts, matrix: int,
                                                                    bandwidth: float) -> None:
    """A compensated waveform meets the same limits as the one it replaces."""
    module = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth, _null_moment_order=1)
    raster = float(opts.grad_raster_time)

    for lobe in module.prephaser_lobes:
        amplitude = abs(float(lobe.amplitude))
        assert amplitude <= float(opts.max_grad) * (1.0 + 1e-9)
        assert amplitude / float(lobe.rise_time) <= float(opts.max_slew) * (1.0 + 1e-9)
        duration = float(pp.calc_duration(lobe))
        assert round(duration / raster, 9) == pytest.approx(round(duration / raster), abs=1e-9)


@pytest.mark.parametrize(('matrix', 'bandwidth'), PROTOCOLS)
def test_the_pair_is_the_shortest_this_realisation_can_build(opts, matrix: int,
                                                             bandwidth: float) -> None:
    """
    Minimal **for this realisation**, tested on the lobes pypulseq actually builds.

    The realisation gives the two winders the same duration; whether letting them differ could
    be shorter is not established and is not searched for. What is checked is what the public
    number means: at `prephaser_duration_s` the built pair is legal, and one raster step shorter
    it is not -- so a bisection that stopped early, which would cost every repetition, fails
    here.

    The predicate is itself the built pair -- the lobes pypulseq emits, checked against the
    gradient and slew limits -- so there is one authority for legality rather than a smooth
    estimate and a correction after it.
    """
    module = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth, _null_moment_order=1)
    raster = float(opts.grad_raster_time)
    half = module.prephaser_duration_s / 2.0

    assert module._compensated_pair_fits(half)
    assert not module._compensated_pair_fits(half - raster)


@pytest.mark.parametrize(('matrix', 'bandwidth'), PROTOCOLS)
def test_the_feasibility_predicate_is_monotone_in_the_lobe_duration(opts, matrix: int,
                                                                    bandwidth: float) -> None:
    """
    **What licenses the bisection**, and the reason the predicate carries the limit areas.

    ``|A1|`` and ``|A2|`` each fall as `D` grows while the area a trapezoid of duration `D` can
    carry rises, so once the pair fits every longer pair fits.  Bisection on a predicate that was
    not monotone could return a feasible duration that is not the shortest one.

    Asserted here as one transition across a swept range, which is the property itself; the
    signs it follows from are checked separately.
    """
    module = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth, _null_moment_order=1)
    raster = float(opts.grad_raster_time)
    half = module.prephaser_duration_s / 2.0

    steps = [module._compensated_pair_fits(n * raster)
             for n in range(1, int(round(half / raster)) + 60)]
    transitions = sum(1 for a, b in zip(steps, steps[1:]) if a != b)
    assert transitions == 1, 'infeasible then feasible, once'
    assert steps[-1], 'the long end is the feasible one'


@pytest.mark.parametrize('geometry', list(GEOMETRIES))
def test_the_premises_the_monotonicity_argument_rests_on(opts, geometry: str) -> None:
    """
    The signs that make the real predicate monotone, checked rather than assumed.

    The argument is that ``|A1|`` and ``|A2|`` are each a positive constant plus ``C/D``, so both
    decrease as `D` grows while the area a trapezoid can carry increases.  That needs three
    things to be true of this construction, and each is measured here:

    * the readout has one sign before the first echo, so its area to the echo is positive;
    * ``C = m0_ro e - M`` is non-negative, because it is ``integral of g(t) t dt`` over an
      interval where both factors are;
    * and therefore ``A1 > 0`` and ``A2 < 0``, with both magnitudes falling as `D` grows.

    If a future geometry broke any of them the bisection would no longer be licensed, and this
    is where that would show.
    """
    module = line(opts, _null_moment_order=1, **GEOMETRIES[geometry])
    area, moment = (module._readout_moment_to_echo(order) for order in (0, 1))
    constant = area * module._echo_in_gx - -moment

    assert area > 0.0, 'one sign before the first echo'
    assert constant >= 0.0, 'C is the integral of a non-negative product'

    half = module.prephaser_duration_s / 2.0
    first, second = module._compensated_areas(half)
    assert first > 0.0 and second < 0.0
    longer = module._compensated_areas(half * 2.0)
    assert abs(longer[0]) < abs(first) and abs(longer[1]) < abs(second)


def test_the_compiler_accepts_a_compensated_readout(opts) -> None:
    """The contract is the emitted file, so the block has to survive compilation."""
    seq = sc.compile(line(opts, _null_moment_order=1)(), opts)

    assert len(seq.block_events) >= 3


@pytest.mark.parametrize('echoes', (1, 4))
def test_a_requested_duration_is_honoured_and_still_nulls_the_moment(opts, echoes: int) -> None:
    """
    Stretching the pair is what a composing kernel does when another axis needs longer.

    The areas are re-solved against the requested duration rather than scaled, so the first
    moment is still null -- which a design that had solved once and stretched afterwards would
    not give.
    """
    train = {'polarity': 'monopolar'} if echoes > 1 else {}
    shortest = line(opts, echoes=echoes, _null_moment_order=1, **train)
    stretched = line(opts, echoes=echoes, _null_moment_order=1, **train,
                     prephaser_duration_s=shortest.prephaser_duration_s + 200e-6)

    assert stretched.prephaser_duration_s > shortest.prephaser_duration_s
    assert moment_to_echo(stretched, 1) == pytest.approx(0.0, abs=tolerance(stretched, 1))
    assert moment_to_echo(stretched, 0) == pytest.approx(0.0, abs=tolerance(stretched, 0))


@pytest.mark.parametrize('polarity', ('monopolar', 'bipolar'))
@pytest.mark.parametrize('echoes', (2, 3, 5))
def test_every_acquired_echo_is_compensated(opts, polarity: str, echoes: int) -> None:
    """
    **The contract, measured at every echo rather than only the first.**

    The winder solves the first echo; the waveform between each pair of lobes carries the
    condition to the next one, by nulling the interval's zeroth *and* first moment.  Asserted at
    every echo of trains of three different lengths, because a design that served the first
    interval and drifted afterwards would pass a two-echo check.

    Both moments, because they are one condition: ``m0 = 0`` is what makes the instant an echo at
    all, and asserting ``m1`` at an instant that is not ``k = 0`` would measure nothing.
    """
    module = line(opts, echoes=echoes, polarity=polarity, _null_moment_order=1)

    assert len(module.te_s) == echoes
    for echo in range(echoes):
        assert moment_to_echo(module, 0, echo) == pytest.approx(0.0, abs=tolerance(module, 0))
        assert moment_to_echo(module, 1, echo) == pytest.approx(0.0, abs=tolerance(module, 1))


@pytest.mark.parametrize('polarity', ('monopolar', 'bipolar'))
def test_an_uncompensated_train_still_drifts(opts, polarity: str) -> None:
    """
    The control: without the option the later echoes are **not** compensated, as they never were.

    Worth keeping now that the compensated path passes, because a measurement that cannot fail
    proves nothing -- this is what says the assertion above is reading the design and not the
    tolerance.
    """
    module = line(opts, echoes=3, polarity=polarity)

    for echo in range(3):
        assert abs(moment_to_echo(module, 1, echo)) > tolerance(module, 1) * 1e6


@pytest.mark.parametrize('polarity', ('monopolar', 'bipolar'))
def test_the_transition_nulls_the_interval_and_not_itself(opts, polarity: str) -> None:
    """
    **What the waveform between the lobes is solved against.**

    The condition is on the interval ``[TE_i, TE_i+1]`` as a whole -- the post-echo tail of one
    lobe, the transition, and the pre-echo head of the next -- not on the transition's own first
    moment, which is a different and weaker statement.  So the transition's own ``m1`` about its
    own start is deliberately **not** zero, and this says so, because a future reader who
    "simplifies" it to the weaker condition would break every echo but the first.
    """
    module = line(opts, echoes=3, polarity=polarity, _null_moment_order=1)
    # Referenced to the interval's origin -- the echo -- which is where the window starts from,
    # not to the transition's own first knot.  The two differ by `window_start * m0`, which is
    # the whole reason a "the fly-back nulls its own first moment" shortcut is wrong.
    echo_s = module._echo_in_lobe_s(0)
    at = module._lobe_s
    own_m0 = own_m1 = 0.0
    for event in module.transition_events[0]:
        knots, amps = knots_of(event, at)
        own_m0 += float(pwl_moment(knots - echo_s, amps, 0))
        own_m1 += float(pwl_moment(knots - echo_s, amps, 1))
        at += float(pp.calc_duration(event))

    fixed_m0, fixed_m1 = module._interval_fixed(0, module.transition_duration_s)
    assert own_m0 == pytest.approx(-fixed_m0, rel=1e-9, abs=1e-9)
    assert own_m1 == pytest.approx(-fixed_m1, rel=1e-9)

    # And about its own start it is *not* zero, which is the weaker condition it is not solving.
    about_itself = own_m1 - (module._lobe_s - echo_s) * own_m0
    assert abs(about_itself) > tolerance(module, 1) * 1e6


@pytest.mark.parametrize('geometry', list(GEOMETRIES))
def test_the_first_echo_is_compensated_in_every_supported_geometry(opts, geometry: str) -> None:
    """
    **The contract holds wherever the option is accepted, not only for a centred single echo.**

    The solve uses the readout lobe's *measured* moments up to the echo, so partial Fourier --
    which moves the echo away from the lobe's midpoint -- and a reversed first lobe are handled
    by the same arithmetic rather than by a case each.  Worth asserting because "it is measured,
    so it must be right" is exactly the reasoning that hides a geometry nobody tried.
    """
    module = line(opts, _null_moment_order=1, **GEOMETRIES[geometry])

    assert moment_to_echo(module, 1, 0) == pytest.approx(0.0, abs=tolerance(module, 1))
    assert moment_to_echo(module, 0, 0) == pytest.approx(0.0, abs=tolerance(module, 0))


@pytest.mark.parametrize('geometry', list(GEOMETRIES))
def test_every_supported_geometry_stays_inside_the_gradient_system(opts, geometry: str) -> None:
    """A geometry that met the moment condition by exceeding the amplifier would not count."""
    module = line(opts, _null_moment_order=1, **GEOMETRIES[geometry])

    for lobe in module.prephaser_lobes:
        amplitude = abs(float(lobe.amplitude))
        assert amplitude <= float(opts.max_grad) * (1.0 + 1e-9)
        assert amplitude / float(lobe.rise_time) <= float(opts.max_slew) * (1.0 + 1e-9)


# ------------------------------------------------------------------------------ the refusals
def test_an_order_above_one_is_refused(opts) -> None:
    """
    Acceleration compensation needs a fourth lobe and is deferred, so it refuses.

    The message names the physics rather than the parameter: nulling this is not something a
    caller can ask for by keyword, so quoting one would send them looking for a spelling that
    does not exist.
    """
    with pytest.raises(sc.ConfigurationError, match='cannot null moments above the first') as got:
        line(opts, _null_moment_order=2)

    assert 'velocity compensation' in str(got.value)
    assert 'fourth lobe' in str(got.value)


def test_compensation_without_a_prephaser_is_refused(opts) -> None:
    """There is nothing to reshape: the compensation lives in the winder."""
    with pytest.raises(sc.ConfigurationError, match='no prephaser') as got:
        line(opts, prephase=False, _null_moment_order=1)

    assert 'prephase=True' in str(got.value)


def test_a_prephaser_too_short_for_the_compensated_areas_is_refused(opts) -> None:
    """Naming the floor, because "too short" without a number is a search for the caller."""
    module = line(opts, _null_moment_order=1)
    floor = module.prephaser_duration_s

    with pytest.raises(sc.ConfigurationError, match='prephaser_duration_s') as caught:
        line(opts, _null_moment_order=1, prephaser_duration_s=floor / 2.0)

    assert f'{floor:.6g}' in str(caught.value)
    assert line(opts, _null_moment_order=1, prephaser_duration_s=floor).prephaser_duration_s == (
        pytest.approx(floor, abs=1e-12))


# ---------------------------------------------------- the timing an all-echo train has to buy
@pytest.mark.parametrize('polarity', ('monopolar', 'bipolar'))
def test_the_automatic_spacing_is_the_shortest_that_works(opts, polarity: str) -> None:
    """
    ``echo_spacing_s=None`` is the shortest period this realisation family can compensate at.

    Shortest is asserted the only way it can be: one gradient raster step below it, the design
    refuses.  That is a statement about *this* family and these ``Opts``, not a claim that no
    waveform anywhere is shorter, which is why the refusal says so too.
    """
    module = line(opts, echoes=4, polarity=polarity, _null_moment_order=1)
    minimum = module.min_echo_spacing_s

    assert module.echo_spacing_s == pytest.approx(minimum, abs=1e-12)
    with pytest.raises(sc.ConfigurationError):
        line(opts, echoes=4, polarity=polarity, _null_moment_order=1,
             echo_spacing_s=minimum - float(opts.grad_raster_time))


@pytest.mark.parametrize('polarity', ('monopolar', 'bipolar'))
def test_the_automatic_spacing_lands_on_the_raster(opts, polarity: str) -> None:
    """An ESP off the gradient raster is a block the compiler refuses, however good the physics."""
    module = line(opts, echoes=4, polarity=polarity, _null_moment_order=1)
    raster = float(opts.grad_raster_time)

    assert module.echo_spacing_s / raster == pytest.approx(round(module.echo_spacing_s / raster))
    assert module.transition_duration_s / raster == pytest.approx(
        round(module.transition_duration_s / raster))


@pytest.mark.parametrize('polarity', ('monopolar', 'bipolar'))
def test_a_spacing_that_is_too_short_reports_the_minimum(opts, polarity: str) -> None:
    """
    A refusal that does not say what would have worked makes the caller guess.

    The number in the message is the one that works, so it can be pasted back in -- asserted by
    pasting it back in.
    """
    module = line(opts, echoes=4, polarity=polarity, _null_moment_order=1)

    with pytest.raises(sc.ConfigurationError) as raised:
        line(opts, echoes=4, polarity=polarity, _null_moment_order=1, echo_spacing_s=200e-6)
    assert 'min_echo_spacing_s' in str(raised.value)
    assert f'{module.min_echo_spacing_s:.6g}' in str(raised.value)
    assert line(opts, echoes=4, polarity=polarity, _null_moment_order=1,
                echo_spacing_s=module.min_echo_spacing_s) is not None


@pytest.mark.parametrize('polarity', ('monopolar', 'bipolar'))
def test_a_longer_requested_spacing_is_honoured_and_still_compensated(opts,
                                                                      polarity: str) -> None:
    """A hard timing request is a request, not a hint -- and it does not cost the physics."""
    module = line(opts, echoes=4, polarity=polarity, _null_moment_order=1)
    wanted = module.min_echo_spacing_s + 500e-6
    longer = line(opts, echoes=4, polarity=polarity, _null_moment_order=1, echo_spacing_s=wanted)

    assert longer.echo_spacing_s == pytest.approx(wanted, abs=1e-12)
    assert longer.transition_duration_s > module.transition_duration_s
    for echo in range(4):
        assert moment_to_echo(longer, 0, echo) == pytest.approx(0.0, abs=tolerance(longer, 0))
        assert moment_to_echo(longer, 1, echo) == pytest.approx(0.0, abs=tolerance(longer, 1))


@pytest.mark.parametrize('polarity', ('monopolar', 'bipolar'))
def test_compensation_costs_echo_spacing(opts, polarity: str) -> None:
    """
    The guarantee is paid for in time, and the payment is visible rather than silent.

    Both polarities pay; the bipolar train pays far more, because its lobes abut and the
    transition it needs has to be created rather than reshaped.
    """
    plain = line(opts, echoes=4, polarity=polarity)
    compensated = line(opts, echoes=4, polarity=polarity, _null_moment_order=1)

    assert compensated.echo_spacing_s > plain.echo_spacing_s
    assert compensated.transition_duration_s > plain.transition_duration_s


@pytest.mark.parametrize('polarity', ('monopolar', 'bipolar'))
def test_the_minimum_spacing_depends_on_the_hardware(opts, derated_opts, polarity: str) -> None:
    """
    A weaker amplifier needs longer to carry the same moment, and the design says so.

    This is the property that makes the timing *derived* rather than a formula: nothing in the
    module knows what a scanner can do except ``Opts``.
    """
    strong = line(opts, echoes=4, polarity=polarity, _null_moment_order=1)
    weak = line(derated_opts, echoes=4, polarity=polarity, _null_moment_order=1)

    assert weak.min_echo_spacing_s > strong.min_echo_spacing_s
    for echo in range(4):
        assert moment_to_echo(weak, 1, echo) == pytest.approx(0.0, abs=tolerance(weak, 1))


# ------------------------------------------------- what a reconstruction still has to be able to do
@pytest.mark.parametrize('polarity', ('monopolar', 'bipolar'))
def test_the_sampling_semantics_survive_compensation(opts, polarity: str) -> None:
    """
    The transition changes the gradients between echoes and **nothing a reconstruction reads**.

    Echo count, which sample is ``k = 0``, which lobes are reversed and the dwell all come back
    unchanged; only the period moves, and ``te_s`` moves with it.  A compensated train whose
    ``REV`` flags had quietly stopped matching its lobes would reconstruct mirrored.
    """
    plain = line(opts, echoes=4, polarity=polarity)
    compensated = line(opts, echoes=4, polarity=polarity, _null_moment_order=1)

    assert compensated.num_samples == plain.num_samples
    assert compensated.dwell_s == plain.dwell_s
    assert compensated.pre_echo_samples == plain.pre_echo_samples
    for echo in range(4):
        assert compensated.echo_sample(echo) == plain.echo_sample(echo)
        assert compensated.polarity_of(echo) == plain.polarity_of(echo)
    assert float(compensated.gx.amplitude) == pytest.approx(float(plain.gx.amplitude))
    assert float(compensated.gx.area) == pytest.approx(float(plain.gx.area))


@pytest.mark.parametrize('polarity', ('monopolar', 'bipolar'))
def test_the_echo_times_are_the_period_plus_the_echo_in_its_lobe(opts, polarity: str) -> None:
    """
    ``te_s`` stays the arithmetic it was, because one transition duration serves the whole train.

    A per-interval duration would have been the cheaper design and would have made the echo times
    a cumulative sum instead -- so this asserts the property that choice bought.
    """
    module = line(opts, echoes=5, polarity=polarity, _null_moment_order=1)
    spacing = np.diff(module.te_s)

    assert len(module.te_s) == 5
    if polarity == 'monopolar':
        assert spacing == pytest.approx(module.echo_spacing_s, abs=1e-12)
    else:
        # Still alternating by two dwells: that is the sample grid, which compensation does not
        # touch.  It is `te_s`'s whole reason for existing, and it survives.
        assert np.ptp(spacing) == pytest.approx(2 * module.dwell_s, abs=1e-9)


@pytest.mark.parametrize('polarity', ('monopolar', 'bipolar'))
def test_the_compensated_train_compiles(opts, polarity: str) -> None:
    """Legal physics that the compiler refuses is not a design, it is a proposal."""
    module = line(opts, echoes=4, polarity=polarity, _null_moment_order=1)

    sequence = sc.compile(sc.LogicBlock('train').add(0.0, module()), opts)

    assert sequence.block_events


@pytest.mark.parametrize('polarity', ('monopolar', 'bipolar'))
def test_every_readout_axis_gradient_is_still_a_trapezoid(opts, polarity: str) -> None:
    """
    A transition of two adjacent lobes merges into one extended gradient without a barrier.

    That is legal, simulates fine and is reported only as a merge warning -- and it changes the
    waveform the moments were solved for, so it is asserted rather than hoped for.
    """
    module = line(opts, echoes=4, polarity=polarity, _null_moment_order=1)
    sequence = sc.compile(sc.LogicBlock('train').add(0.0, module()), opts)

    for index in range(1, len(sequence.block_events) + 1):
        gradient = getattr(sequence.get_block(index), module.axis, None)
        if gradient is not None:
            assert gradient.type == 'trap', f'block {index} is a {gradient.type}'
