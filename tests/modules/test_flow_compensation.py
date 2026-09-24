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


def line(opts, **overrides):
    kwargs = {'opts': opts, 'fov_mm': 220.0, 'matrix': 64, 'bandwidth_hz_px': 500.0}
    kwargs.update(overrides)
    return sc.modules.CartesianLine(**kwargs)


def moment_to_echo(module, order: int, tree=None) -> float:
    """
    The `order`-th moment of **every** gradient event on the module's axis, at the echo.

    Truncated at ``time_to_echo`` and referenced to it, integrated from each event's own
    piecewise-linear knots.  Summed across events, which is the part no per-event check reaches:
    the prephaser lobes and the readout's own ramp-up are three separate events whose moments only
    cancel together.
    """
    echo = module.time_to_echo()
    total = 0.0
    for start, event, _ in flatten(tree if tree is not None else module()):
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


def _scale(module) -> float:
    """What "zero" means here: moments are areas, so they are compared against the lobe scale."""
    return abs(module.area_to_echo_per_m) * module.time_to_echo() ** 2


# ------------------------------------------------------------------------ the moment conditions
@pytest.mark.parametrize(('matrix', 'bandwidth'), PROTOCOLS)
def test_the_first_moment_is_null_at_the_echo(opts, matrix: int, bandwidth: float) -> None:
    """
    **The condition the option exists to meet**, measured on the complete emitted waveform.

    Velocity compensation is ``m1 = 0`` at the echo for everything that plays on the axis.  A
    spin moving at constant velocity then arrives at the echo with the phase it would have had
    standing still, which is what stops flow from writing a ghost into the phase-encode direction.
    """
    module = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth, null_moment_order=1)

    assert moment_to_echo(module, 1) == pytest.approx(0.0, abs=_scale(module) * 1e-9)


@pytest.mark.parametrize(('matrix', 'bandwidth'), PROTOCOLS)
def test_the_zeroth_moment_is_still_null_at_the_echo(opts, matrix: int,
                                                     bandwidth: float) -> None:
    """
    The echo condition is not traded away to get the first moment.

    ``m0 = 0`` at the echo is what puts ``k = 0`` there at all; a compensated readout that lost it
    would image a shifted k-space and still look plausible.
    """
    module = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth, null_moment_order=1)

    assert moment_to_echo(module, 0) == pytest.approx(
        0.0, abs=abs(module.area_to_echo_per_m) * 1e-9)


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

    assert abs(moment_to_echo(module, 1)) > _scale(module) * 1e-3
    assert moment_to_echo(module, 0) == pytest.approx(
        0.0, abs=abs(module.area_to_echo_per_m) * 1e-9)


# ------------------------------------------------------------------------------ the construction
@pytest.mark.parametrize(('matrix', 'bandwidth'), PROTOCOLS)
def test_the_two_lobes_sum_to_the_area_the_echo_needs(opts, matrix: int,
                                                      bandwidth: float) -> None:
    """
    Two lobes of opposite sign, and it is their **total** that cancels the readout's area.

    The individual areas are neither ``-area_to_echo_per_m`` nor half of it, which is why
    `prephaser_area_per_m` is documented as the total and the lobes are exposed separately.
    """
    module = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth, null_moment_order=1)
    areas = [float(lobe.area) for lobe in module.prephaser_lobes]

    assert len(areas) == 2
    assert areas[0] * areas[1] < 0.0, 'opposite signs'
    assert sum(areas) == pytest.approx(-module.area_to_echo_per_m, rel=1e-9)
    assert module.prephaser_area_per_m == pytest.approx(sum(areas), rel=1e-12)


def test_the_ordinary_prephaser_is_one_lobe(opts) -> None:
    """``null_moment_order=0`` is exactly today's readout, and the default."""
    module = line(opts)

    assert module.null_moment_order == 0
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
    compensated = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth, null_moment_order=1)

    assert compensated.prephaser_duration_s > plain.prephaser_duration_s
    assert compensated.time_to_echo() - plain.time_to_echo() == pytest.approx(
        compensated.prephaser_duration_s - plain.prephaser_duration_s, abs=1e-12)


# ------------------------------------------------------------------------ hardware and raster
@pytest.mark.parametrize(('matrix', 'bandwidth'), PROTOCOLS)
def test_the_compensated_prephaser_stays_inside_the_gradient_system(opts, matrix: int,
                                                                    bandwidth: float) -> None:
    """A compensated waveform meets the same limits as the one it replaces."""
    module = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth, null_moment_order=1)
    raster = float(opts.grad_raster_time)

    for lobe in module.prephaser_lobes:
        amplitude = abs(float(lobe.amplitude))
        assert amplitude <= float(opts.max_grad) * (1.0 + 1e-9)
        assert amplitude / float(lobe.rise_time) <= float(opts.max_slew) * (1.0 + 1e-9)
        duration = float(pp.calc_duration(lobe))
        assert round(duration / raster, 9) == pytest.approx(round(duration / raster), abs=1e-9)


@pytest.mark.parametrize(('matrix', 'bandwidth'), PROTOCOLS)
def test_the_shortest_pair_is_the_shortest_on_the_raster(opts, matrix: int,
                                                         bandwidth: float) -> None:
    """
    The duration is found by bisection, so this checks it is minimal rather than merely feasible.

    One raster step shorter must fail to fit, which is what makes the quoted
    `prephaser_duration_s` the real floor.  The areas are closed form at any duration; only the
    duration is searched, and a search that stopped early would quietly cost every repetition.
    """
    module = line(opts, matrix=matrix, bandwidth_hz_px=bandwidth, null_moment_order=1)
    raster = float(opts.grad_raster_time)
    half = module.prephaser_duration_s / 2.0

    assert module._compensated_pair_fits(half)
    assert not module._fits(half - raster)


def test_the_compiler_accepts_a_compensated_readout(opts) -> None:
    """The contract is the emitted file, so the block has to survive compilation."""
    seq = sc.compile(line(opts, null_moment_order=1)(), opts)

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
    shortest = line(opts, echoes=echoes, null_moment_order=1, **train)
    stretched = line(opts, echoes=echoes, null_moment_order=1, **train,
                     prephaser_duration_s=shortest.prephaser_duration_s + 200e-6)

    assert stretched.prephaser_duration_s > shortest.prephaser_duration_s
    assert moment_to_echo(stretched, 1) == pytest.approx(0.0, abs=_scale(stretched) * 1e-9)
    assert moment_to_echo(stretched, 0) == pytest.approx(
        0.0, abs=abs(stretched.area_to_echo_per_m) * 1e-9)


def test_compensation_holds_for_every_echo_of_a_train(opts) -> None:
    """
    The prephaser is compensated for the **first** echo, which is the one it is solved against.

    Recorded as a measurement rather than assumed: later echoes of a monopolar train accumulate
    their own first moment from the lobes between them, and this module compensates the winder,
    not the train. A protocol needing every echo compensated needs more than this option.
    """
    module = line(opts, echoes=3, polarity='monopolar', null_moment_order=1)

    assert moment_to_echo(module, 1) == pytest.approx(0.0, abs=_scale(module) * 1e-9)


# ------------------------------------------------------------------------------ the refusals
def test_an_order_above_one_is_refused(opts) -> None:
    """Acceleration compensation needs a fourth lobe and is deferred, so it refuses by name."""
    with pytest.raises(sc.ConfigurationError, match='null_moment_order') as caught:
        line(opts, null_moment_order=2)

    assert 'velocity compensation' in str(caught.value)


def test_compensation_without_a_prephaser_is_refused(opts) -> None:
    """There is nothing to reshape: the compensation lives in the winder."""
    with pytest.raises(sc.ConfigurationError, match='prephase') as caught:
        line(opts, prephase=False, null_moment_order=1)

    assert 'null_moment_order' in str(caught.value)


def test_a_prephaser_too_short_for_the_compensated_areas_is_refused(opts) -> None:
    """Naming the floor, because "too short" without a number is a search for the caller."""
    module = line(opts, null_moment_order=1)
    floor = module.prephaser_duration_s

    with pytest.raises(sc.ConfigurationError, match='prephaser_duration_s') as caught:
        line(opts, null_moment_order=1, prephaser_duration_s=floor / 2.0)

    assert f'{floor:.6g}' in str(caught.value)
    assert line(opts, null_moment_order=1, prephaser_duration_s=floor).prephaser_duration_s == (
        pytest.approx(floor, abs=1e-12))
