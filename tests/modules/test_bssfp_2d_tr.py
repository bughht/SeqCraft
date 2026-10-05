"""
``bSSFP2DTR`` -- the balance condition, measured where it is defined.

**Every balance measurement here integrates the compiled waveforms from one RF effective centre
to the next.**  Not the module's own reported areas, which would be checking its arithmetic
rather than its output, and not a ``LogicBlock``'s edges, which are a container boundary that
merely resembles the physical interval.

The interval ends inside the *next* repetition, so a single block cannot be measured at all:
these tests build ``N + 1`` compatible consecutive repetitions and measure the first ``N``
intervals.  **Compatible** is the module's own word -- the selection-gradient geometry at the
boundary has to match, while the flip angle and the phase-encode line are free to differ.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc

MATRIX = (64, 32)
#: The one protocol every test varies from, so a neighbour differs in exactly one field.
SPEC = dict(fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0, flip_deg=35.0,
            bandwidth_hz_px=800.0)
#: Areas here are hundreds of 1/m, so this is about twelve orders below the quantity measured.
ZERO = 1e-9


@pytest.fixture(scope='module')
def tr(opts):
    """One repetition of a small but complete 2D bSSFP, at the canonical symmetric default."""
    return sc.modules.bSSFP2DTR(opts=opts, **SPEC)


# ------------------------------------------------------------------ the measurement apparatus
def train(opts, *reps, lines=None, phases=None):
    """
    Compile a train of repetitions laid end to end, and return the sequence.

    `reps` are the modules, one per repetition, so a test can vary the flip angle between
    neighbours.  Lines default to consecutive, because the compiler rejects two readouts writing
    one k-space address and a repeated line is therefore not a sequence anybody can run.
    """
    lines = list(range(len(reps))) if lines is None else lines
    phases = [180.0 * (k % 2) for k in range(len(reps))] if phases is None else phases
    scan, at = sc.LogicBlock(), 0.0
    for rep, line, phase in zip(reps, lines, phases):
        scan.add(at, rep(line=line, phase_deg=phase))
        at += rep.tr_s
    return sc.compile(scan, opts=opts)


def rf_centres(seq) -> np.ndarray:
    """
    Absolute times of every RF effective centre, seconds.

    The effective centre is a property of the waveform rather than of the block, so it is read
    through ``pp.calc_rf_center`` exactly as :meth:`Excitation.time_to_center` does -- a test
    that assumed the midpoint would be measuring a different instant for a minimum-phase pulse.
    """
    out, t = [], 0.0
    for i in range(1, len(seq.block_events) + 1):
        block = seq.get_block(i)
        if getattr(block, 'rf', None) is not None:
            out.append(t + float(block.rf.delay) + float(pp.calc_rf_center(block.rf)[0]))
        t += float(seq.block_durations[i])
    return np.array(out)


def m0_between(seq, a: float, b: float) -> np.ndarray:
    """
    ``(M0x, M0y, M0z)`` over ``[a, b]``, 1/m, integrated from the emitted waveforms.

    The interval edges are added to the waveform's own knots rather than rounded onto them: an
    RF effective centre is not on any raster, and snapping it to one would move the boundary of
    the very interval under test.
    """
    waveforms = seq.waveforms_and_times()[0]
    out = []
    for axis in range(3):
        t, g = waveforms[axis][0], waveforms[axis][1]
        if len(t) == 0:
            out.append(0.0)
            continue
        grid = np.unique(np.concatenate([t[(t > a) & (t < b)], [a, b]]))
        out.append(float(np.trapezoid(np.interp(grid, t, g, left=0.0, right=0.0), grid)))
    return np.array(out)


def rf_to_rf_m0(seq) -> np.ndarray:
    """``M0`` over every RF-centre-to-RF-centre interval in `seq`: shape ``(N, 3)``."""
    centres = rf_centres(seq)
    return np.array([m0_between(seq, centres[k], centres[k + 1])
                     for k in range(len(centres) - 1)])


def echo_time_s(seq, rep) -> float:
    """Seconds from the first RF effective centre to the first echo sample."""
    t, adc_start = 0.0, None
    for i in range(1, len(seq.block_events) + 1):
        block = seq.get_block(i)
        if getattr(block, 'adc', None) is not None and adc_start is None:
            adc_start = t + float(block.adc.delay)
        t += float(seq.block_durations[i])
    echo = adc_start + (rep.ro.echo_sample(0) + 0.5) * float(rep.ro.adc.dwell)
    return echo - rf_centres(seq)[0]


# ------------------------------------------------------------------------- claim 1: balance
def test_the_rf_to_rf_interval_carries_no_net_area(opts, tr) -> None:
    """
    The defining condition, on the axis of every gradient that plays.

    Three repetitions give two intervals.  The third exists only to close the second: its own
    interval would end inside a fourth repetition that was never built, which is what it means
    for the contract to be compositional.
    """
    m0 = rf_to_rf_m0(train(opts, tr, tr, tr))

    assert len(m0) == 2
    assert np.abs(m0).max() < ZERO


def test_the_measured_interval_is_the_repetition_time(opts, tr) -> None:
    """RF centre to RF centre is TR, which is what lets a caller stack repetitions at ``n*TR``."""
    centres = rf_centres(train(opts, tr, tr, tr))

    assert np.diff(centres) == pytest.approx(tr.tr_s)


# ------------------------------------------------------------- claims 2 and 3: the echo itself
def test_the_requested_line_is_encoded_at_the_echo(opts, tr) -> None:
    """``k_y`` at the echo is the line that was asked for, and the other two axes are at zero."""
    for line in (0, MATRIX[1] // 2, MATRIX[1] - 1):
        seq = train(opts, tr, tr, lines=[line, (line + 1) % MATRIX[1]])
        centres = rf_centres(seq)
        k = m0_between(seq, centres[0], centres[0] + echo_time_s(seq, tr))

        assert k[1] == pytest.approx(tr.pe.k_per_m(line))
        assert abs(k[0]) < ZERO
        assert abs(k[2]) < ZERO


def test_the_echo_arrives_at_the_declared_te(opts, tr) -> None:
    """The declared TE is measured from the RF effective centre, not from the block."""
    seq = train(opts, tr, tr)

    assert echo_time_s(seq, tr) == pytest.approx(tr.te_s)


# ---------------------------------------------------------------- claim 4: the symmetric default
def test_passing_nothing_gives_the_canonical_symmetric_realisation(tr, opts) -> None:
    """
    ``TE = TR/2``, which is what ``te_s=None`` means here.

    Exactly, up to :attr:`bSSFP2DTR.symmetry_residual_s` -- the echo is the ADC sample where
    ``k = 0`` and that sample is half a dwell off the window's centre, while every fill has to
    be a whole number of rasters.
    """
    assert abs(tr.te_s - tr.tr_s / 2) <= opts.grad_raster_time / 2
    assert tr.symmetry_residual_s == pytest.approx(tr.te_s - tr.tr_s / 2)


def test_the_default_is_not_the_earliest_reachable_echo(tr) -> None:
    """
    The one place ``None`` means something different from what it means on ``GRE2DTR``.

    If this ever becomes an equality, the default has silently become "minimum TE" and every
    caller who passed nothing is getting an asymmetric acquisition they did not ask for.
    """
    assert tr.te_s > tr.min_te_s


def test_a_requested_tr_still_puts_the_echo_at_the_middle(opts, tr) -> None:
    """``te_s=None`` with an explicit ``tr_s`` is still the symmetric realisation."""
    longer = sc.modules.bSSFP2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                                  flip_deg=35.0, bandwidth_hz_px=800.0,
                                  tr_s=tr.min_tr_s + 2e-3)

    assert longer.tr_s == pytest.approx(tr.min_tr_s + 2e-3, abs=opts.grad_raster_time)
    assert abs(longer.symmetry_residual_s) <= opts.grad_raster_time


# -------------------------------------------------- claim 5: the definition, not the realisation
def test_an_asymmetric_te_is_still_balanced(opts, tr) -> None:
    """
    **The point of the whole contract.**  Balance is zero net area per axis, and nothing else.

    A repetition whose echo is nowhere near the midpoint is still bSSFP if this holds, so a test
    that only ever measured the symmetric case would be testing ``TE = TR/2`` and calling it
    balance.
    """
    asymmetric = sc.modules.bSSFP2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                                      flip_deg=35.0, bandwidth_hz_px=800.0,
                                      te_s=tr.min_te_s, tr_s=tr.tr_s)
    m0 = rf_to_rf_m0(train(opts, asymmetric, asymmetric, asymmetric))

    assert asymmetric.te_s < asymmetric.tr_s / 2 - 1e-4      # genuinely off the midpoint
    assert np.abs(m0).max() < ZERO


# ------------------------------------------------------------ claims 6 and 7: the two phases
def test_the_rf_carries_the_requested_phase(opts, tr) -> None:
    """The progression itself belongs to the caller; placing it on the RF belongs here."""
    seq = train(opts, tr, tr, tr, tr, phases=[0.0, 180.0, 0.0, 180.0])
    phases = [float(seq.get_block(i).rf.phase_offset)
              for i in range(1, len(seq.block_events) + 1)
              if getattr(seq.get_block(i), 'rf', None) is not None]

    assert np.allclose(np.mod(phases, 2 * np.pi), [0.0, np.pi, 0.0, np.pi])


def test_the_receiver_follows_the_transmitter(opts, tr) -> None:
    """
    One argument reaches both, because a phase on the RF alone writes itself into ``ky``.

    Checked against the RF of the *same* repetition rather than against the number the test
    passed in: that is the comparison that fails if the two ever come from different places.
    """
    seq = train(opts, tr, tr, tr, tr, phases=[0.0, 180.0, 90.0, 270.0])
    rf, adc = [], []
    for i in range(1, len(seq.block_events) + 1):
        block = seq.get_block(i)
        if getattr(block, 'rf', None) is not None:
            rf.append(float(block.rf.phase_offset))
        if getattr(block, 'adc', None) is not None:
            adc.append(float(block.adc.phase_offset))

    assert np.allclose(np.mod(rf, 2 * np.pi), np.mod(adc, 2 * np.pi))


# ------------------------------------------------------------------ claim 8: minimum timing
def test_the_reported_minimum_repetition_is_achievable(opts, tr) -> None:
    """A minimum that cannot be built is a refusal threshold, not a minimum."""
    shortest = sc.modules.bSSFP2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                                    flip_deg=35.0, bandwidth_hz_px=800.0, tr_s=tr.min_tr_s)

    assert shortest.tr_s == pytest.approx(tr.min_tr_s)
    assert np.abs(rf_to_rf_m0(train(opts, shortest, shortest, shortest))).max() < ZERO


@pytest.mark.parametrize(('field', 'attr'), [('te_s', 'min_te_s'), ('tr_s', 'min_tr_s')])
def test_a_request_below_the_minimum_names_the_achievable_value(opts, tr, field, attr) -> None:
    """The number a caller needs is in the message, not in a docstring they have to go find."""
    with pytest.raises(sc.errors.ConfigurationError) as err:
        sc.modules.bSSFP2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                             flip_deg=35.0, bandwidth_hz_px=800.0,
                             **{field: getattr(tr, attr) * 0.5})

    assert attr in str(err.value)


def test_a_symmetric_tr_is_refused_below_what_symmetry_costs(opts, tr) -> None:
    """
    ``min_tr_s`` is the shortest *symmetric* repetition, which is longer than the shortest one.

    A caller who wants the shorter one is asking for an asymmetric realisation and says so by
    passing `te_s`, which is exactly the distinction the refusal points at.
    """
    assert tr.min_tr_s > tr.min_te_s + tr.min_te_s

    relaxed = sc.modules.bSSFP2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                                   flip_deg=35.0, bandwidth_hz_px=800.0, te_s=tr.min_te_s)

    assert relaxed.tr_s < tr.min_tr_s


# --------------------------------------------------------------- claim 9: the block and hardware
def test_the_block_lasts_exactly_one_repetition_time(tr) -> None:
    """What lets the layer above stack repetitions at ``n * tr_s`` and nothing else."""
    assert tr(line=3).duration == pytest.approx(tr.tr_s)


def test_the_repetition_is_pure_and_compiles(opts, tr, component_checks) -> None:
    """The two checks the compiler cannot make, plus the ones it can."""
    component_checks.all(tr, line=5)


# ------------------------------------------------- claim 10: balanced is not in steady state
def test_balance_holds_on_the_very_first_interval(opts, tr) -> None:
    """
    Two repetitions, nothing in steady state, and the waveform condition already holds.

    This is the test that keeps *balanced* and *in steady state* separable. The first is true of
    the first interval a scanner ever plays; the second is a property of a train and of the
    tissue, and this module never claims it.
    """
    m0 = rf_to_rf_m0(train(opts, tr, tr))

    assert len(m0) == 1
    assert np.abs(m0).max() < ZERO


# ------------------------------------------------------ claim 11: balance across unlike neighbours
def test_balance_survives_neighbours_encoding_different_lines(opts, tr) -> None:
    """
    The y rewind belongs to the repetition that applied the blip, not to a block edge.

    Lines far apart rather than adjacent, so a rewind that was short by one line's worth would
    not be hiding inside the tolerance.
    """
    m0 = rf_to_rf_m0(train(opts, tr, tr, tr, lines=[0, MATRIX[1] - 1, MATRIX[1] // 2]))

    assert np.abs(m0).max() < ZERO


def test_balance_survives_a_start_up_flip_ramp(opts) -> None:
    """
    A neighbour with a different flip angle is **compatible**: the angle rides the RF, not Gz.

    This is what makes a start-up ramp something a caller can write as an ordinary loop. If the
    flip angle ever reached the selection gradient, this is the test that would say so.
    """
    ramp = [sc.modules.bSSFP2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                                 flip_deg=flip, bandwidth_hz_px=800.0)
            for flip in (8.75, 17.5, 26.25, 35.0)]

    assert len({rep.tr_s for rep in ramp}) == 1          # compatible: one timing, four angles
    assert np.abs(rf_to_rf_m0(train(opts, *ramp))).max() < ZERO


#: Heterogeneous neighbours, and what each one changes.  ``timing`` is whether the neighbour
#: leaves ``time_to_rf_center()`` alone, which is what simple block stacking needs.
NEIGHBOURS = [
    ('flip_deg', 10.0, True),
    ('bandwidth_hz_px', 400.0, True),
    ('thickness_mm', 8.0, False),
    ('rf_duration_s', 2e-3, False),
]


@pytest.mark.parametrize(('field', 'value'), [(f, v) for f, v, _ in NEIGHBOURS])
def test_balance_survives_a_neighbour_of_different_geometry(opts, tr, field, value) -> None:
    """
    The geometry at the boundary does **not** have to match, and this is why.

    Each repetition discharges its own selection halves -- the trailing lobe cancels ``A_post``
    and the leading winder cancels ``A_pre`` -- so the interval is two independent
    cancellations rather than one shared sum, and neither repetition is built against an
    assumption about the other.  A thicker slice, a longer pulse and a slower readout all change
    the selection gradient or the window, and none of them unbalances the interval.

    The lumped realisation would fail every one of these, which is the argument for this split.
    """
    other = sc.modules.bSSFP2DTR(**{**SPEC, 'opts': opts, field: value})
    seq = train(opts, tr, other, other, lines=[4, 5, 6])

    assert np.abs(rf_to_rf_m0(seq)).max() < ZERO


# ----------------------------------------------- balance-compatible is not timing-compatible
def test_the_rf_centre_is_where_te_is_measured_from(tr) -> None:
    """
    The identity that keeps the timing origin in one place.

    ``time_to_rf_center`` exists so that nothing downstream reconstructs the origin from a pulse
    duration or assumes it is halfway through the event -- which for a minimum-phase pulse it is
    not.
    """
    assert tr.time_to_echo() == pytest.approx(tr.time_to_rf_center() + tr.te_s, abs=1e-12)
    assert tr.time_to_rf_center() == pytest.approx(tr._lead_s + tr.exc.time_to_center())


def test_a_homogeneous_train_gets_the_declared_tr_from_block_stacking(opts, tr) -> None:
    """``block.duration == tr_s`` delivers the declared TR only because the offsets agree."""
    centres = rf_centres(train(opts, tr, tr, tr))

    assert np.diff(centres) == pytest.approx(tr.tr_s)
    assert tr(line=1).duration == pytest.approx(tr.tr_s)


@pytest.mark.parametrize(('field', 'value', 'timing_compatible'), NEIGHBOURS)
def test_timing_compatibility_is_measured_rather_than_inferred_from_balance(
    opts, tr, field, value, timing_compatible,
) -> None:
    """
    Balance-compatible does not imply timing-compatible, and this is where the two part.

    Stacking by block duration puts the next block at ``start + tr_s``, so the interval the
    magnetisation sees is ``block_duration(n) + c(n+1) - c(n)``.  A neighbour that moves the
    RF-centre offset therefore misses the declared TR by the difference, **while still
    balancing** -- the companion test above measures that half.

    The parametrisation carries the expected verdict rather than only the compatible cases, so
    a change that quietly made one of the incompatible ones agree would fail here instead of
    passing silently.
    """
    other = sc.modules.bSSFP2DTR(**{**SPEC, 'opts': opts, field: value})
    offset_delta = other.time_to_rf_center() - tr.time_to_rf_center()
    centres = rf_centres(train(opts, tr, other, other, lines=[4, 5, 6]))

    assert (abs(offset_delta) < 1e-12) is timing_compatible
    assert centres[1] - centres[0] == pytest.approx(tr.tr_s + offset_delta)
    if timing_compatible:
        assert centres[1] - centres[0] == pytest.approx(tr.tr_s)


def test_placing_from_the_offsets_restores_the_declared_tr(opts, tr) -> None:
    """
    The composition formula a heterogeneous train needs, asserted rather than only documented.

    ``start(n+1) = start(n) + tr_s + c(n) - c(n+1)``.  This module does not do the placement --
    that would be an acquisition framework, and this is a repetition -- but the contract it
    documents has to be one that works.
    """
    other = sc.modules.bSSFP2DTR(**{**SPEC, 'opts': opts, 'rf_duration_s': 2e-3})
    train_of = (tr, other, other)
    scan, at = sc.LogicBlock(), 0.0
    for k, rep in enumerate(train_of):
        if k:
            previous = train_of[k - 1]
            at += tr.tr_s + previous.time_to_rf_center() - rep.time_to_rf_center()
        scan.add(at, rep(line=4 + k))
    centres = rf_centres(sc.compile(scan, opts=opts))

    assert np.diff(centres) == pytest.approx(tr.tr_s)


def test_a_neighbour_that_keeps_no_leading_winder_is_outside_the_contract(opts, tr) -> None:
    """
    The other half of stating a contract: saying what it does **not** cover.

    The condition is on the contract the neighbour keeps, not on the protocol it was configured
    with.  A neighbour that plays no leading winder leaves its own pre-centre half uncancelled
    in the interval -- a spoiled gradient-echo repetition is exactly such a neighbour -- and the
    interval is then unbalanced by that half.  Nothing in this module is wrong when that
    happens; the train is.

    Built by removing the lobe rather than by composing a different class, so that what is
    being tested is the one term and not a second module's arithmetic.
    """
    from seqcraft.design.events import derive

    no_winder = sc.modules.bSSFP2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                                     flip_deg=35.0, bandwidth_hz_px=800.0)
    no_winder._z_lead = derive(pp.make_trapezoid(channel='z', area=0.0, system=opts,
                                                 duration=no_winder.winder_s))
    m0 = rf_to_rf_m0(train(opts, tr, no_winder, lines=[4, 5]))

    assert abs(m0[0][2]) == pytest.approx(abs(tr._a_pre), rel=1e-6)


# --------------------------------------------------------------------- what it does not own
def test_nothing_here_spoils(tr) -> None:
    """A spoiler is a deliberate non-zero net area, which is the condition this module negates."""
    assert not hasattr(tr, 'spoilers')
    assert 'spoil' not in str(sc.modules.bSSFP2DTR.__init__.__doc__ or '')


def test_a_start_up_repetition_loads_the_same_gradients(opts, tr) -> None:
    """
    ``acquire=False`` drops the ADC and the label and nothing else.

    What start-up *means* is the caller's; all this module offers is a repetition that does not
    sample, which is why the gradients have to be identical to an acquired one's.
    """
    def areas(block):
        return sorted(round(float(n.item.area), 9) for n in block
                      if getattr(n.item, 'type', '') == 'trap')

    assert areas(tr(line=7, acquire=False)) == areas(tr(line=7, acquire=True))
    assert np.abs(rf_to_rf_m0(train(opts, tr, tr, lines=[7, 8]))).max() < ZERO
