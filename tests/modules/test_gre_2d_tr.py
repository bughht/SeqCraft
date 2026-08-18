"""
``GRE2DTR`` -- one repetition, and the arithmetic no caller should have to repeat.

Everything here is a property of the *composition* rather than of any leaf: the winder match, TE
and TR, the rewinder actually rewinding, and the ``LIN`` label naming the line that was acquired.
Each is checked against the compiled sequence or against pypulseq's own k-space integration,
because those are the two independent oracles available.
"""

from __future__ import annotations

import numpy as np
import pytest

import seqcraft as sc

MATRIX = (64, 32)


@pytest.fixture(scope='module')
def tr(opts):
    """One repetition of a small but complete 2D GRE."""
    return sc.modules.GRE2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0)


def _labels(block) -> list[int]:
    """``LIN`` values sitting directly in a repetition."""
    return [n.item.value for n in block if getattr(n.item, 'type', '') == 'labelset']


# ------------------------------------------------------------------- the winder coupling
def test_the_winder_is_the_longest_of_three(opts, tr) -> None:
    """
    Slice rephaser on z, phase-encode blip on y, readout prephaser on x: all at once.

    Each reports its own minimum and accepts an override, the composite takes the maximum and
    passes it down, and no leaf knows the other two exist.
    """
    # A nanosecond of slack, because the maximum is snapped onto the gradient raster and can
    # land one ulp below the participant it came from.
    assert tr.pe(line=0).duration == pytest.approx(tr.winder_s)
    assert tr.ro.prephaser_duration_s == pytest.approx(tr.winder_s)
    assert tr.winder_s >= tr.pe.min_duration_s - 1e-9
    assert tr.winder_s >= tr.exc.rephaser_duration_s - 1e-9


def test_all_three_winder_axes_start_together(opts, tr) -> None:
    """
    Overlap on different axes is free, so nothing waits for anything.

    The rephaser is inside the excitation's block, so its instant is measured through
    ``time_to_rephaser`` rather than off a node -- which is the whole reason that method exists.
    """
    block = tr(line=17)
    starts: dict[str, list[float]] = {}
    for node in block:
        if isinstance(node.item, sc.LogicBlock):
            starts.setdefault(node.item.tag, []).append(node.start)

    winder_start = starts['CartesianLine'][0]
    # Two PhaseEncode blocks per repetition -- the blip and the rewinder -- so it is the first
    # that shares an instant with the readout.
    assert starts['PhaseEncode'][0] == winder_start
    assert starts['PhaseEncode'][1] > winder_start

    rephaser_start = starts['Excitation'][0] + tr.exc.time_to_rephaser()
    assert rephaser_start <= winder_start + 1e-12, 'the rephaser is not being waited for'
    assert rephaser_start + tr.exc.rephaser_duration_s <= winder_start + tr.winder_s + 1e-12, (
        'the readout gradient would start while the slice was still being rephased'
    )


def test_overlapping_the_rephaser_is_what_shortens_te(opts, tr) -> None:
    """
    The measurement behind the decision, so a future rearrangement has to beat it.

    Placing the winder after the whole excitation block instead of alongside its tail costs the
    rephaser's entire duration in TE -- here about a tenth of the echo time, given away by a
    placement decision rather than by any physics.
    """
    waited = tr.exc().duration - tr.exc.time_to_rephaser()

    assert waited > 0.4e-3, 'this protocol has a rephaser long enough for the test to mean something'
    assert tr.min_te_s + waited > tr.min_te_s * 1.05


def test_the_receiver_and_the_transmitter_share_one_phase(opts, tr) -> None:
    """
    The composite is the layer that holds both events, so it is where they are made to agree.

    A schedule applied to the RF alone writes its quadratic phase into ``ky``; on a single
    off-centre voxel that scattered a point source over the phase-encode direction and shifted
    its peak by thirteen pixels, with the readout direction perfectly correct beside it.
    """
    block = tr(line=17, phase_deg=117.0)
    events = {e.type: e for _, e, _ in sc.flatten(block) if getattr(e, 'type', '') in ('rf', 'adc')}

    assert events['rf'].phase_offset == pytest.approx(np.deg2rad(117.0))
    assert events['adc'].phase_offset == pytest.approx(events['rf'].phase_offset)


# ---------------------------------------------------------------------------- timing
def test_the_measured_echo_time_is_the_reported_one(opts, tr) -> None:
    """
    TE from the compiled sequence: the RF's effective centre to k = 0 at the ADC.

    ``calculate_kspacePP`` supplies both, so nothing in this assertion comes from the module's
    own arithmetic.
    """
    k = sc.kspace(sc.LogicBlock('probe').add(0.0, tr(line=17)), opts)

    echo_at = k['t_adc'][tr.ro.pre_echo_samples]

    assert echo_at - float(k['t_excitation'][0]) == pytest.approx(tr.te_s, abs=1e-9)
    assert echo_at == pytest.approx(tr.time_to_echo(), abs=1e-9)


def test_k_is_zero_on_every_axis_at_the_echo(opts, tr) -> None:
    """kx from the prephaser, ky from the line being 0, kz from the slice rewinder."""
    centre = sc.modules.GRE2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0)
    k = sc.kspace(sc.LogicBlock('probe').add(0.0, centre(line=centre.center_line)), opts)

    at_echo = k['k_adc'][:, centre.ro.pre_echo_samples]

    assert np.allclose(at_echo, 0.0, atol=1e-3)


def test_the_block_lasts_exactly_one_tr(opts, tr) -> None:
    """
    The fill is what makes the measured duration equal TR, which is what lets the layer above
    stack repetitions at ``n * tr_s`` and nothing else.
    """
    assert tr(line=0).duration == pytest.approx(tr.tr_s, abs=1e-12)
    assert tr(line=31).duration == pytest.approx(tr.tr_s, abs=1e-12)


def test_a_requested_te_lengthens_the_repetition(opts, tr) -> None:
    longer = sc.modules.GRE2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                                te_s=tr.min_te_s + 2e-3)

    assert longer.te_s == pytest.approx(tr.min_te_s + 2e-3, abs=1e-5)
    assert longer.min_tr_s > tr.min_tr_s, 'a later echo is a longer minimum TR'

    k = sc.kspace(sc.LogicBlock('probe').add(0.0, longer(line=17)), opts)
    measured = k['t_adc'][longer.ro.pre_echo_samples] - float(k['t_excitation'][0])
    assert measured == pytest.approx(longer.te_s, abs=1e-9)


def test_a_requested_tr_is_the_block_duration(opts, tr) -> None:
    longer = sc.modules.GRE2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                                tr_s=50e-3)

    assert longer.tr_s == pytest.approx(50e-3)
    assert longer(line=0).duration == pytest.approx(50e-3, abs=1e-12)


@pytest.mark.parametrize(('kwarg', 'attribute'), [('te_s', 'min_te_s'), ('tr_s', 'min_tr_s')])
def test_a_time_below_the_minimum_raises_naming_it(opts, tr, kwarg, attribute) -> None:
    """Never silently lengthened: the message says what is achievable, and the caller decides."""
    with pytest.raises(sc.ConfigurationError, match=attribute):
        sc.modules.GRE2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                           **{kwarg: getattr(tr, attribute) / 2})


# ----------------------------------------------------------------------- the rewinder
def test_the_phase_encode_moment_over_one_tr_is_zero(opts, tr) -> None:
    """
    The rewinder actually rewinds, for every line rather than for the one that was checked.

    Without it the phase encoding accumulates across repetitions, which for the first few lines
    looks like nothing at all.
    """
    for line in (0, 7, tr.center_line, 31):
        assert sc.moments(tr(line=line))['y'] == pytest.approx(0.0, abs=1e-9)


def test_the_slice_axis_ends_on_the_spoiler_rather_than_at_zero(opts, tr) -> None:
    """
    Deliberate, and worth stating so it is not read as the previous test failing.

    The z axis carries the spoiler, whose whole job is to leave dephasing behind: four turns
    across a 5 mm slice is 800 1/m.  The other term is the selection gradient's *leading* half,
    which the rephaser does not undo and should not -- spins are not transverse before the
    pulse's centre.  So a raw moment over one repetition is the sum of the two, and the sum is
    what a reader should expect to see here rather than zero.
    """
    spoiler_area = tr.spoil_cycles_per_voxel / (tr.voxel_mm('z') / 1e3)
    before_the_pulse = float(tr.exc.gz.area) + float(tr.exc.gzr.area)

    assert sc.moments(tr(line=0))['z'] == pytest.approx(
        spoiler_area + before_the_pulse, rel=1e-6,
    )


# ---------------------------------------------------------------------- the spoilers
def test_spoiling_is_more_than_one_axis_by_default(opts, tr) -> None:
    """
    ``('x', 'z')``, because a gradient dephases only along its own direction.

    Winding cycles on ``z`` alone leaves the residual perfectly coherent in ``x``, so the FID and
    the stimulated echoes that follow it reach the next readout carrying a ``ky`` that does not
    belong to it -- which is a ghost along the phase-encode direction rather than anything that
    reads as a spoiling problem.
    """
    assert tr.spoil_axis == ('x', 'z')
    assert sorted(tr.spoilers) == ['x', 'z']


def test_each_axis_winds_its_cycles_across_its_own_voxel(opts, tr) -> None:
    """
    The in-plane voxel and the slice differ by a factor of three here, so one shared area would
    be a different amount of spoiling on each axis while reading as the same number.
    """
    assert tr.voxel_mm('z') == tr.thickness_mm
    assert tr.voxel_mm('x') == pytest.approx(tr.fov_mm[0] / tr.matrix[0])
    assert tr.voxel_mm('y') == pytest.approx(tr.fov_mm[1] / tr.matrix[1])

    for axis, block in tr.spoilers.items():
        area = float(block.nodes[0].item.area)
        assert area == pytest.approx(tr.spoil_cycles_per_voxel / (tr.voxel_mm(axis) / 1e3))


def test_the_readout_axis_does_almost_no_spoiling_of_its_own(opts) -> None:
    """
    The measurement that makes ``'x'`` a default rather than an option.

    After the echo only half the readout's area is left on ``kx``, which is half a cycle across
    one voxel -- so the largest gradient in the repetition contributes almost nothing, and a
    reader who assumed otherwise would leave the axis unspoiled.
    """
    bare = sc.modules.GRE2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                              spoil_axis='z')

    voxel_m = bare.voxel_mm('x') / 1e3
    cycles = sc.moments(bare(line=0))['x'] * voxel_m

    assert 0.4 < cycles < 0.6, f'{cycles:.3f} cycles per voxel from the readout alone'


def test_naming_y_adds_a_spoiler_beside_the_rewinder_not_instead_of_it(opts) -> None:
    """
    So the net dephasing on ``y`` is the same every repetition, whatever line was acquired.

    Replacing the rewinder would make it line-dependent, which is a different sequence and a
    worse one: the point of spoiling is that every repetition starts from the same place.
    """
    spoiled = sc.modules.GRE2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                                 spoil_axis=('x', 'y', 'z'))

    expected = spoiled.spoil_cycles_per_voxel / (spoiled.voxel_mm('y') / 1e3)
    for line in (0, 7, spoiled.center_line, 31):
        assert sc.moments(spoiled(line=line))['y'] == pytest.approx(expected, rel=1e-6)


def test_one_axis_is_still_allowed_and_is_shorter(opts, tr) -> None:
    """The default costs tail time, so turning it off has to remain a one-word change."""
    bare = sc.modules.GRE2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                              spoil_axis='z')

    assert bare.spoil_axis == ('z',)
    assert bare.min_tr_s < tr.min_tr_s
    assert sc.moments(bare(line=0))['x'] < sc.moments(tr(line=0))['x']


def test_an_unknown_or_empty_spoil_axis_raises(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='spoil_axis'):
        sc.modules.GRE2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                           spoil_axis='xz')
    with pytest.raises(sc.ConfigurationError, match='nothing would be spoiled'):
        sc.modules.GRE2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0,
                           spoil_axis=())


# ------------------------------------------------------------------------- the label
def test_every_repetition_carries_the_line_it_acquired(opts, tr) -> None:
    for line in (0, 1, tr.center_line, 31):
        assert _labels(tr(line=line)) == [line]


def test_the_label_reaches_the_adc_in_the_compiled_sequence(opts, tr) -> None:
    """
    Emitted before the readout in the tree, and assigned to the block holding its ADC by the
    compiler -- which is the property that makes a label a k-space address rather than a
    position in a file.
    """
    lines = [3, tr.center_line, 30]
    scan = sc.LogicBlock('probe')
    for n, line in enumerate(lines):
        scan.add(n * tr.tr_s, tr(line=line))

    seq = sc.compile(scan, opts)
    seen = seq.evaluate_labels(evolution='adc')['LIN']

    assert [int(v) for v in np.atleast_1d(np.asarray(seen))] == lines


def test_the_centre_line_is_one_convention_with_three_consumers(opts, tr) -> None:
    """``LIN`` at the centre, ``center_line``, and ``kSpaceCenterLine`` in the file."""
    seq = sc.compile(
        sc.LogicBlock('probe').add(0.0, tr(line=tr.center_line)), opts,
        definitions={'kSpaceCenterLine': tr.center_line},
    )

    assert _labels(tr(line=tr.center_line)) == [tr.center_line]
    assert tr.center_line == MATRIX[1] // 2
    assert seq.definitions['kSpaceCenterLine'] == tr.center_line


# ------------------------------------------------------------------------ the dummy
def test_a_dummy_loads_identical_gradients(opts, tr) -> None:
    """
    Same duration, same waveforms, no ADC and no label.

    A dummy that does not load the gradients identically establishes a different steady state
    from the one that gets acquired, and nothing about the file looks wrong.
    """
    real, dummy = tr(line=17), tr(line=17, acquire=False)

    assert dummy.duration == pytest.approx(real.duration, abs=1e-12)
    assert _labels(dummy) == []
    assert not [e for _, e, _ in sc.flatten(dummy) if getattr(e, 'type', '') == 'adc']
    assert sc.moments(dummy) == pytest.approx(sc.moments(real))


# ------------------------------------------------------------------------ the geometry
def test_the_pair_is_split_at_this_level(opts) -> None:
    """
    ``matrix`` is ``(nx, ny)`` here and scalar in the leaves.

    The leaves each act on one axis and do not know the other exists, which is what lets
    ``PhaseEncode`` serve ``y`` today and ``z`` when 3D arrives.
    """
    rect = sc.modules.GRE2DTR(opts=opts, fov_mm=(250.0, 180.0), matrix=(128, 48),
                              thickness_mm=5.0)

    assert rect.ro.matrix == 128
    assert rect.pe.matrix == 48
    assert rect.ro.fov_mm == 250.0
    assert rect.pe.fov_mm == 180.0
    assert isinstance(rect.pe.matrix, int), 'the leaf never sees the pair'


def test_a_scalar_fov_means_square(opts) -> None:
    square = sc.modules.GRE2DTR(opts=opts, fov_mm=250.0, matrix=MATRIX, thickness_mm=5.0)

    assert square.fov_mm == (250.0, 250.0)


def test_a_three_element_matrix_raises(opts) -> None:
    with pytest.raises(sc.ConfigurationError, match='pair'):
        sc.modules.GRE2DTR(opts=opts, fov_mm=250.0, matrix=(64, 32, 16), thickness_mm=5.0)


# ------------------------------------------------------------------------ center_mm
def test_the_slice_offset_reaches_the_rf_and_the_readout_offset_the_adc(opts, tr) -> None:
    """Two unrelated mechanisms on two different events, held by the layer that holds both."""
    block = tr(line=17, center_mm=(30.0, 0.0, 20.0))
    events = {e.type: e for _, e, _ in sc.flatten(block) if getattr(e, 'type', '') in
              ('rf', 'adc')}

    assert events['rf'].freq_offset == pytest.approx(float(tr.exc.gz.amplitude) * 0.020)
    assert events['adc'].freq_offset == pytest.approx(float(tr.ro.gx.amplitude) * 0.030)


def test_a_phase_encode_offset_raises_rather_than_being_ignored(opts, tr) -> None:
    """
    Out of scope, and said so.

    A ``ky`` shift is a per-line ADC phase rather than a per-TR frequency, so it cannot be
    folded in with the other two; the three-tuple is what lets it be added later without
    changing a signature.
    """
    with pytest.raises(sc.ConfigurationError, match='phase-encode'):
        tr(line=17, center_mm=(0.0, 10.0, 0.0))


# ------------------------------------------------------------------------ provenance
def test_the_provenance_path_is_the_tree(opts, tr) -> None:
    """Four leaf tags under one kernel, and not one of them written by hand."""
    tree = sc.LogicBlock('scan').add(0.0, tr(line=17))

    assert {path for _, _, path in sc.flatten(tree)} == {
        ('scan', 'GRE2DTR'),                       # the LIN label and the TR fill
        ('scan', 'GRE2DTR', 'Excitation'),
        ('scan', 'GRE2DTR', 'PhaseEncode'),
        ('scan', 'GRE2DTR', 'CartesianLine'),
        ('scan', 'GRE2DTR', 'spoiler'),
    }


def test_it_is_pure_and_compiles_alone(opts, tr, component_checks) -> None:
    component_checks.all(tr, line=17)
