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
def test_the_first_echo_is_compensated_and_later_echoes_are_not(opts, polarity: str) -> None:
    """
    **The documented limitation, measured on both later echoes.**

    The winder is solved against the first echo, and later echoes of a train accumulate their own
    first moment from the lobes between them. So this option compensates the first echo and does
    not promise the rest of the train -- which is worth locking down, because a multi-echo
    protocol that assumed otherwise would look fine and be uncompensated everywhere but the
    start.
    """
    module = line(opts, echoes=3, polarity=polarity, _null_moment_order=1)

    assert moment_to_echo(module, 1, 0) == pytest.approx(0.0, abs=tolerance(module, 1))
    for echo in (1, 2):
        assert abs(moment_to_echo(module, 1, echo)) > tolerance(module, 1) * 1e6


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
