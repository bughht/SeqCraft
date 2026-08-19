"""
``Refocusing`` -- whether the block is symmetric about the conjugation instant, in time and area.

Two oracles sit under everything here.  The first is :func:`seqcraft.kspace`, which reaches
pypulseq's ``calculate_kspacePP`` and therefore applies the **conjugation** at ``use='refocusing'``
itself -- so a signed k at an ADC sample is a statement about physics rather than about this
module's algebra.  It is checked first, in five raw blocks, because an oracle nothing checks is an
assumption.

The second is the *failure mode*, asserted as its own test.  A module rewritten to be wrong in a new
way still passes a test that only asserts the fix, and the wrong version here is not exotic: it is
what every reference implementation writes -- the same crusher trapezoid twice.  Its symptom is a
``k_z`` that **alternates sign echo to echo**, which is the odd/even modulation an FSE is famous for
and which reads as a hardware fault rather than as a design error.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest
from pypulseq.opts import Opts

import seqcraft as sc
from seqcraft.modules._support import area_until, ceil_raster

#: The protocol every measured number in ``docs/se_api.md`` was taken on.
THICKNESS_MM, CRUSH_VOXEL_MM, CRUSH_CYCLES = 6.25, 5.0, 3.0


@pytest.fixture(scope='module')
def se_opts() -> Opts:
    """The reference spin-echo scanner: 32 mT/m, 130 T/m/s, dead time != ringdown."""
    return Opts(
        max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s', B0=3.0,
        rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
    )


@pytest.fixture(scope='module')
def symmetric_opts(se_opts: Opts) -> Opts:
    """The references' scanner: ``rf_dead_time == rf_ringdown_time``, so the plateau is symmetric.

    This is the configuration that makes ``writeTSE.m``'s two-identical-trapezoids design correct
    **by accident**, which is why the balance is asserted on both.
    """
    return sc.opts.derate(se_opts, rf_ringdown_time=100e-6)


def refocusing(opts: Opts, **kwargs) -> sc.modules.Refocusing:
    kwargs.setdefault('thickness_mm', THICKNESS_MM)
    kwargs.setdefault('crush_voxel_mm', CRUSH_VOXEL_MM)
    kwargs.setdefault('crush_cycles_per_voxel', CRUSH_CYCLES)
    return sc.modules.Refocusing(opts=opts, **kwargs)


# ------------------------------------------------------------------------------- the oracle
def test_kspace_conjugates_at_a_refocusing_pulse(se_opts: Opts) -> None:
    """
    Two **positive** lobes cancelling is the whole test: without conjugation the answer is ``2A``.

    Every physics assertion below reads k out of this function, so it is checked in five raw
    pypulseq blocks first -- and against the same tree with ``use='other'``, which is where it would
    show if pypulseq keyed the conjugation off something other than ``use``.
    """
    def tree(use: str) -> sc.LogicBlock:
        ninety = pp.make_block_pulse(flip_angle=np.pi / 2, duration=200e-6, delay=1e-4,
                                     use='excitation', system=se_opts)
        one_eighty = pp.make_block_pulse(flip_angle=np.pi, duration=400e-6, delay=1e-4, use=use,
                                         system=se_opts)
        lobe = pp.make_trapezoid('x', area=500.0, duration=1e-3, system=se_opts)
        adc = pp.make_adc(num_samples=4, dwell=1e-5, delay=1e-5, system=se_opts)
        out = sc.LogicBlock('probe').add(0.0, ninety)
        out.add(1e-3, lobe).add(3e-3, one_eighty)
        return out.add(5e-3, lobe).add(7e-3, adc)

    refocused = sc.kspace(tree('refocusing'), se_opts)
    unmarked = sc.kspace(tree('other'), se_opts)

    assert refocused['k_adc'][0, 0] == pytest.approx(0.0, abs=1e-6)
    assert unmarked['k_adc'][0, 0] == pytest.approx(1000.0, rel=1e-6), 'two lobes, not conjugated'
    assert len(refocused['t_refocusing']) == 1
    assert len(unmarked['t_refocusing']) == 0


# ------------------------------------------------------------------------------ the invariant
#: Three pulses.  The first is what hides the bug -- a symmetric plateau makes the naive design
#: correct -- and the second and third are what the references get wrong.
PULSES = {
    'symmetric sinc, dead != ringdown': ({}, 'se_opts'),
    'symmetric sinc, dead == ringdown': ({}, 'symmetric_opts'),
    'asymmetric, centre three-quarters through': ({'pulse_opts': {'center_pos': 0.75}}, 'se_opts'),
}


@pytest.mark.parametrize('case', list(PULSES))
def test_the_block_is_symmetric_in_time_and_in_area(request, case) -> None:
    """
    Both halves of the invariant, on all three pulses.  **They fail independently.**

    Balancing the areas makes k exact at the echo and says nothing about *when* the echo lands;
    symmetrising the plateau puts the echo at the midpoint between two pulses and says nothing
    about k.  A train needs both, because its stimulated-echo pathways only coincide with the
    primary echo when the pulses are uniformly spaced *and* the echo is midway between them.
    """
    kwargs, scanner = PULSES[case]
    opts = request.getfixturevalue(scanner)
    module = refocusing(opts, **kwargs)

    assert module.area_to_center_per_m == pytest.approx(module.area_from_center_per_m, abs=1e-9)
    assert module.time_to_center() == pytest.approx(module().duration / 2, abs=1e-12)


@pytest.mark.parametrize('case', list(PULSES))
def test_the_crusher_carries_the_area_it_was_asked_for(request, case) -> None:
    """
    Measured off the waveform with ``sc.moments``, not read back off a parameter.

    The total is twice the pre-centre area, because the block is symmetric about the centre -- so
    this checks the units (cycles across `crush_voxel_mm`) and the symmetry in one number.
    """
    kwargs, scanner = PULSES[case]
    opts = request.getfixturevalue(scanner)
    module = refocusing(opts, **kwargs)

    wanted = CRUSH_CYCLES / (CRUSH_VOXEL_MM / 1e3)

    assert module.area_to_center_per_m == pytest.approx(wanted, rel=1e-12)
    assert sc.moments(module())['z'] == pytest.approx(2 * wanted, rel=1e-9)


def test_time_to_center_is_where_the_rf_peak_actually_is(se_opts: Opts) -> None:
    """
    The package convention is only true if something checks it: every ``time_to_*`` is measured
    from the start of the block that module's ``build`` returns.  A method with a different origin
    would break every timeline built by addition, and nothing else would notice.
    """
    module = refocusing(se_opts)

    _grid, _grads, marks = sc.sample(module(), se_opts)
    rf_spans = [(start, end) for kind, start, end, _ in marks if kind == 'rf']

    assert len(rf_spans) == 1
    start, end = rf_spans[0]
    assert start < module.time_to_center() < end
    assert module.time_to_center() == pytest.approx((start + end) / 2, abs=6e-5), (
        'a symmetric sinc peaks at the middle of its own window, within the dead time'
    )


def test_time_to_crusher_is_where_the_trailing_window_begins(se_opts: Opts) -> None:
    """
    The number a caller drops the readout block at, and the reason it is worth 1.4 ms per echo.

    Everything after it is crusher, so the area from there to the end of the block is what a
    caller's own axes have to fit beside.
    """
    module = refocusing(se_opts)

    assert module.time_to_crusher() == pytest.approx(
        module().duration - module.crush_duration_s, abs=1e-12)
    assert module.time_to_crusher() > module.time_to_center()


# -------------------------------------------------------------------------- the failure mode
def _naive_refocusing(opts: Opts, *, window_s: float, duration_s: float = 4e-3):
    """
    The references' design: the same crusher trapezoid twice, and a plateau merely long enough.

    Two errors, and both compile.  The plateau covers the dead time and the ringdown without being
    *symmetrised* about the pulse, and the crusher area is taken from half the plateau's **total**
    area rather than from the area up to the effective centre.
    """
    rf, gz, _ = pp.make_sinc_pulse(
        flip_angle=np.pi, duration=duration_s, slice_thickness=THICKNESS_MM / 1e3,
        apodization=0.5, time_bw_product=4, delay=float(opts.rf_dead_time), use='refocusing',
        return_gz=True, system=opts,
    )
    bare_s = ceil_raster(float(opts.rf_dead_time) + duration_s + float(opts.rf_ringdown_time),
                        opts.grad_raster_time)
    plateau = pp.make_trapezoid('z', amplitude=float(gz.amplitude), flat_time=bare_s, system=opts)
    crusher = pp.make_trapezoid(
        'z', area=CRUSH_CYCLES / (CRUSH_VOXEL_MM / 1e3) - float(plateau.area) / 2,
        duration=window_s, system=opts,
    )
    plateau_s = float(pp.calc_duration(plateau))
    delay = window_s + float(plateau.rise_time) + float(opts.rf_dead_time)
    block = (sc.LogicBlock('naive')
             .add(0.0, sc.events.derive(rf, delay=delay))
             .add(0.0, crusher).add(window_s, plateau).add(window_s + plateau_s, crusher))
    return block, delay + float(pp.calc_rf_center(rf)[0])


def _k_z_after_each_pulse(opts: Opts, block: sc.LogicBlock, center_s: float, *,
                          echoes: int = 4, spacing_s: float = 20e-3) -> np.ndarray:
    """``k_z`` sampled in the quiet gap after each of `echoes` refocusing pulses."""
    tree = sc.LogicBlock('train').add(0.0, pp.make_block_pulse(
        flip_angle=np.pi / 2, duration=200e-6, delay=float(opts.rf_dead_time),
        use='excitation', system=opts))
    for n in range(echoes):
        start = 2e-3 + n * spacing_s
        tree.add(start, block)
        # In the gap, where nothing is on z, so the sample time is not delicate.
        tree.add(start + center_s + spacing_s / 2,
                 pp.make_adc(num_samples=4, dwell=1e-5, delay=1e-5, system=opts))
    return sc.kspace(tree, opts)['k_adc'][2].reshape(echoes, 4)[:, 0]


def test_two_identical_crushers_alternate_k_z_and_the_solve_does_not(se_opts: Opts) -> None:
    """
    The failure mode as its own test, because a test that only asserts the fix passes on a module
    rewritten to be wrong in a new way.

    The residual is ``a_sel * (dead - ringdown)`` -- 11.2 1/m here, 0.056 cycles across the slice.
    Small, and ``k_n = -k_(n-1) + delta`` turns it into ``delta, 0, delta, 0``: a modulation
    between odd and even echoes that no k-space *extent* check can see and that reads as hardware.
    """
    module = refocusing(se_opts)
    naive, naive_center = _naive_refocusing(se_opts, window_s=module.crush_duration_s)

    wrong = _k_z_after_each_pulse(se_opts, naive, naive_center)
    right = _k_z_after_each_pulse(se_opts, module(), module.time_to_center())

    residual = float(se_opts.max_grad) and float(module.gz.amplitude) * (
        float(se_opts.rf_dead_time) - float(se_opts.rf_ringdown_time))
    assert wrong[0] == pytest.approx(-residual, rel=1e-3)
    assert np.allclose(wrong[::2], wrong[0], atol=1e-6)
    assert np.allclose(wrong[1::2], 0.0, atol=1e-6)
    assert abs(wrong[0]) > 1.0, 'the residual has to be big enough to be worth catching'

    assert np.allclose(right, 0.0, atol=1e-6), 'the solved pair leaves nothing to alternate'


def test_a_symmetric_plateau_hides_it(symmetric_opts: Opts) -> None:
    """
    ``writeTSE.m`` sets its dead time and ringdown to the same 100 us, so it gets this right by
    accident.  Change either and the sequence is quietly wrong -- which is the argument for the
    module doing the arithmetic rather than the protocol happening to.
    """
    module = refocusing(symmetric_opts)
    naive, naive_center = _naive_refocusing(symmetric_opts, window_s=module.crush_duration_s)

    assert np.allclose(_k_z_after_each_pulse(symmetric_opts, naive, naive_center), 0.0, atol=1e-6)


# ------------------------------------------------------------------------------ the refusals
def test_a_two_millisecond_180_is_refused_with_the_duration_that_fixes_it(se_opts: Opts) -> None:
    """
    ``sc.compile`` checks gradient amplitude and slew and never looks at RF amplitude, and pypulseq
    warns and hands the pulse back -- so a 2 ms 180 reaches the console, which refuses it there.

    The message has to carry the floor, because "too big" without a number is a search.
    """
    with pytest.raises(sc.ConfigurationError, match='max_b1') as caught:
        refocusing(se_opts, duration_s=2e-3)

    assert '2.7' in str(caught.value), 'the floor at TBW 4 on a 20 uT system'
    assert '130' in str(caught.value), 'and how far over it the request was'


def test_the_same_shape_at_a_lower_flip_angle_is_allowed(se_opts: Opts) -> None:
    """Peak B1 scales with the flip angle, so the refusal is about the pair and not the duration."""
    assert refocusing(se_opts, duration_s=2e-3, flip_deg=120.0)


def test_a_crusher_window_below_the_area_it_carries_is_refused(se_opts: Opts) -> None:
    module = refocusing(se_opts)

    with pytest.raises(sc.ConfigurationError, match='min_crush_duration_s'):
        refocusing(se_opts, crush_duration_s=module.min_crush_duration_s / 2)


def test_its_own_minimum_is_accepted_back(se_opts: Opts) -> None:
    """A composite passes ``max(refoc.min_crush_duration_s, ...)`` down, which can land one ulp
    below where it came from after snapping onto the raster."""
    module = refocusing(se_opts)

    assert refocusing(se_opts, crush_duration_s=module.min_crush_duration_s).crush_duration_s == \
        pytest.approx(module.min_crush_duration_s)


def test_a_non_selective_pulse_needs_a_crusher_voxel(se_opts: Opts) -> None:
    """The same refusal ``IRPrep.spoil_voxel_mm`` makes, for the same reason: no length to count."""
    with pytest.raises(sc.ConfigurationError, match='crush_voxel_mm'):
        sc.modules.Refocusing(opts=se_opts, thickness_mm=None, pulse='block', duration_s=1e-3)


def test_axis_refuses_any_value_when_non_selective(se_opts: Opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='axis'):
        sc.modules.Refocusing(opts=se_opts, thickness_mm=None, pulse='block', duration_s=1e-3,
                              crush_voxel_mm=5.0, axis='z')


def test_position_mm_refuses_only_a_non_zero_value_when_non_selective(se_opts: Opts) -> None:
    """``None`` for an argument with no identity value, the identity value itself when there is
    one -- the rule ``Excitation`` and ``IRPrep`` both follow."""
    hard = sc.modules.Refocusing(opts=se_opts, thickness_mm=None, pulse='block', duration_s=1e-3,
                                 crush_voxel_mm=5.0)

    assert hard(position_mm=0.0)
    with pytest.raises(sc.ConfigurationError, match='non-selective'):
        hard(position_mm=20.0)


def test_a_crusher_below_the_plateaus_own_half_warns_with_the_number(se_opts: Opts) -> None:
    """
    Legal, and no longer anything a reader would recognise as a crusher: the leading lobe changes
    sign.  The balance still holds, which is why this warns rather than raises.
    """
    with pytest.warns(sc.SeqCraftWarning, match='changes sign'):
        module = refocusing(se_opts, crush_cycles_per_voxel=1.0)

    assert module.area_to_center_per_m == pytest.approx(module.area_from_center_per_m, abs=1e-9)


# ------------------------------------------------------------------- the component contract
def test_it_is_pure_and_compiles_alone(se_opts: Opts, component_checks) -> None:
    """
    The block starts and ends at zero, which is why compiling it on its own is a complete check --
    a module that only works when something happens to be beside it is not reusable.

    Purity is checked after *each* of two calls: the canonical bug is an involution, so comparing
    only before and after the second call finds it back where it started.
    """
    component_checks.all(refocusing(se_opts))
    component_checks.all(refocusing(se_opts, crush_axis='x', crush_voxel_mm=2.0))
    component_checks.all(sc.modules.Refocusing(opts=se_opts, thickness_mm=None, pulse='block',
                                               duration_s=1e-3, crush_voxel_mm=5.0))


def test_the_waveform_starts_and_ends_at_zero(se_opts: Opts) -> None:
    """What makes it splittable at any block boundary the compiler happens to choose."""
    module = refocusing(se_opts)

    for _, event in module.gradients:
        times, amps = sc.events.knots_of(event)
        assert amps[0] == pytest.approx(0.0)
        assert amps[-1] == pytest.approx(0.0)
        assert times[-1] == pytest.approx(float(pp.calc_duration(event)), abs=1e-12)


def test_a_phase_offset_moves_the_pulse_and_nothing_else(se_opts: Opts) -> None:
    """The CPMG relation belongs to whatever holds both pulses, so this is a build argument."""
    module = refocusing(se_opts)

    rf = next(n.item for n in module(phase_deg=90.0) if getattr(n.item, 'type', '') == 'rf')

    assert rf.phase_offset == pytest.approx(np.deg2rad(90.0))
    assert module.rf.phase_offset == 0.0, 'the design was not mutated'


@pytest.mark.parametrize('crush_axis', ['x', 'y'])
def test_a_crusher_on_another_axis_still_balances(se_opts: Opts, crush_axis) -> None:
    """
    The selection axis is the only one whose crusher can be *fused* with anything -- but the
    invariant is per axis, so an unfused pair has to hold it too, and on the selection axis at the
    same time.
    """
    module = refocusing(se_opts, crush_axis=crush_axis, crush_voxel_mm=2.0)

    assert module.area_to_center_per_m == pytest.approx(module.area_from_center_per_m, abs=1e-9)
    assert module.time_to_center() == pytest.approx(module().duration / 2, abs=1e-12)
    # ...and the selection gradient is symmetric about the centre on its own, with no solve.
    gz = next(event for _, event in module.gradients if event.channel == 'z')
    start = next(t for t, event in module.gradients if event is gz)
    half = area_until(gz, module.time_to_center() - start)
    assert half == pytest.approx(area_until(gz, float(pp.calc_duration(gz))) - half, abs=1e-9)


def test_a_fixed_rise_time_changes_the_amplitudes_and_not_the_balance(se_opts: Opts) -> None:
    """
    ``rise_time_s`` is slew headroom, not speed: measured, it changes the echo spacing by exactly
    zero, because the spacing depends on the crusher *window* and not on its ramps.
    """
    slewed = refocusing(se_opts)
    fixed = refocusing(se_opts, rise_time_s=250e-6)

    assert fixed().duration == pytest.approx(slewed().duration)
    assert fixed.area_to_center_per_m == pytest.approx(fixed.area_from_center_per_m, abs=1e-9)
    peak = [np.abs(sc.events.knots_of(event)[1]).max() for _, event in fixed.gradients]
    assert max(peak) > 0.0


def test_plateau_padding_is_the_price_of_the_time_symmetry(se_opts: Opts,
                                                           symmetric_opts: Opts) -> None:
    """Reported rather than hidden: ``|dead - ringdown|`` for a symmetric pulse, and twice the
    pulse's own asymmetry for one whose centre is not its midpoint."""
    assert refocusing(se_opts).plateau_padding_s == pytest.approx(
        abs(float(se_opts.rf_dead_time) - float(se_opts.rf_ringdown_time)), abs=1e-9)
    assert refocusing(symmetric_opts).plateau_padding_s == pytest.approx(0.0, abs=1e-9)
    asymmetric = refocusing(se_opts, pulse_opts={'center_pos': 0.75})
    assert asymmetric.plateau_padding_s > refocusing(se_opts).plateau_padding_s
