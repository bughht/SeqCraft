"""
``bSSFP3DTR`` -- the balance condition, measured where it is defined.

**Every balance measurement here integrates the compiled waveforms from one RF effective centre
to the next**, because that is the interval the contract is about.  A tree-level moment is not
the same quantity: it runs from block start to block end, which is neither of those instants, and
it reads zero for decompositions that are not balanced at all.

The test this file exists for is :func:`test_balance_survives_any_partition_transition`.  The
interval ends inside the *successor*, so a decomposition can hold ``M0 = 0`` for a train of
identical repetitions and still be wrong -- and a probe that repeats one partition cannot tell the
two apart.  Every neighbouring pair below changes **both** indices.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft import LogicBlock
from seqcraft.errors import ConfigurationError
from seqcraft.modules import bSSFP3DTR

#: The compiled waveforms balance to about 1e-12 1/m, set by float accumulation over the train
#: rather than by anything the design does.  A tenth of a nanometre of k.
ZERO = 1e-9

FOV = (200.0, 200.0, 160.0)
MATRIX = (64, 64, 16)
SLAB_MM = 100.0


@pytest.fixture(scope='module')
def opts() -> pp.Opts:
    return pp.Opts(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6, B0=3.0)


@pytest.fixture(scope='module')
def tr(opts: pp.Opts) -> bSSFP3DTR:
    return bSSFP3DTR(opts=opts, fov_mm=FOV, matrix=MATRIX, slab_thickness_mm=SLAB_MM)


# --------------------------------------------------------------------------- measuring helpers
def rf_centres(seq) -> list[float]:
    """Every RF **effective** centre in the compiled sequence, seconds from the start."""
    t, out = 0.0, []
    for i in range(1, len(seq.block_events) + 1):
        block = seq.get_block(i)
        rf = getattr(block, 'rf', None)
        if rf is not None:
            centre, _ = pp.calc_rf_center(rf)
            out.append(t + float(rf.delay) + centre)
        t += float(seq.block_durations[i])
    return out


def echo_times(seq, tr: bSSFP3DTR) -> list[float]:
    """Every ``k = 0`` sample time, seconds from the start of the compiled sequence."""
    t, out = 0.0, []
    for i in range(1, len(seq.block_events) + 1):
        block = seq.get_block(i)
        if getattr(block, 'adc', None) is not None:
            out.append(t + float(block.adc.delay)
                       + (tr.ro.echo_sample(0) + 0.5) * float(tr.ro.adc.dwell))
        t += float(seq.block_durations[i])
    return out


def moment(seq, a: float, b: float, axis: str, order: int = 0, ref: float | None = None) -> float:
    """``M_order`` of the compiled waveform on `axis` over ``[a, b]``, about `ref`."""
    times, amps = seq.waveforms_and_times()[0]['xyz'.index(axis)]
    t, g = np.asarray(times, float), np.asarray(amps, float)
    if t.size == 0:
        return 0.0
    grid = np.unique(np.concatenate([t[(t > a) & (t < b)], [a, b]]))
    y = np.interp(grid, t, g, left=0.0, right=0.0)
    if order == 1:
        y = y * (grid - (a if ref is None else ref))
    return float(np.trapezoid(y, grid))


def train(tr: bSSFP3DTR, views, *, acquire: bool = True):
    """A compiled train of `views`, stacked at the declared TR."""
    root, t = LogicBlock(), 0.0
    for line, partition in views:
        root.add(t, tr.build(line=line, partition=partition, acquire=acquire))
        t += tr.tr_s
    return sc.compile(root, tr.opts)


def violent_views(tr: bSSFP3DTR) -> list[tuple[int, int]]:
    """A view order where every neighbouring pair changes **both** indices."""
    ny, nz = tr.matrix[1], tr.matrix[2]
    cl, cp = tr.center_line, tr.center_partition
    return [(cl, 0), (0, nz - 1), (ny - 1, cp), (1, 1), (cl, nz - 1),
            (ny - 1, 0), (cp, cp), (0, nz - 2), (ny - 1, nz - 1), (cl, 2)]


# ------------------------------------------------------------------- claim 1: it is balanced
def test_balance_survives_any_partition_transition(tr: bSSFP3DTR) -> None:
    """
    The property the three-lobe decomposition exists for, and the one a fixed probe cannot see.

    ``lead`` carries no partition term and ``tail``/``rewind`` carry ``+K(p)`` and ``-K(p)`` in
    the *same* repetition, so the interval reduces to the selection halves and their balancing
    lobes whatever the successor acquires.  Give the encoding to the lobe before the next pulse
    instead and the interval carries ``K(p_n) - K(p_{n+1})``, which this train would catch and a
    train of identical repetitions would not.
    """
    views = violent_views(tr)
    seq = train(tr, views)
    centres = rf_centres(seq)

    assert len(centres) == len(views)
    for k in range(len(centres) - 1):
        for axis in 'xyz':
            m0 = moment(seq, centres[k], centres[k + 1], axis)
            assert abs(m0) < ZERO, (
                f'{views[k]} -> {views[k + 1]} left {m0:.3e} 1/m on {axis}')


def test_the_transitions_really_do_move_the_partition(tr: bSSFP3DTR) -> None:
    """Guards the guard: a view order that did not move would make the test above vacuous."""
    views = violent_views(tr)
    steps = [abs(tr.pe_z.k_per_m(b) - tr.pe_z.k_per_m(a))
             for (_, a), (_, b) in zip(views, views[1:])]

    assert min(steps) > 0.0, 'every neighbouring pair must change the partition'
    assert max(steps) > 0.9 * (tr.pe_z.k_per_m(tr.matrix[2] - 1)
                               - tr.pe_z.k_per_m(0)), 'and one pair must cross the table'


def test_a_homogeneous_train_is_balanced_too(tr: bSSFP3DTR) -> None:
    """The easy case still has to hold; it is simply not sufficient on its own."""
    seq = train(tr, [(tr.center_line, tr.center_partition)] * 3, acquire=False)
    centres = rf_centres(seq)

    for k in range(len(centres) - 1):
        for axis in 'xyz':
            assert abs(moment(seq, centres[k], centres[k + 1], axis)) < ZERO


# ------------------------------------------------------------------ claim 2: it encodes k
def test_the_echo_lands_on_the_requested_k_space_point(tr: bSSFP3DTR) -> None:
    """Balance is not enough: a repetition that rewinds the wrong partition is still wrong."""
    views = violent_views(tr)
    seq = train(tr, views)
    centres, echoes = rf_centres(seq), echo_times(seq, tr)

    for i, (line, partition) in enumerate(views):
        kx = moment(seq, centres[i], echoes[i], 'x')
        ky = moment(seq, centres[i], echoes[i], 'y')
        kz = moment(seq, centres[i], echoes[i], 'z')

        assert abs(kx) < ZERO, f'{(line, partition)}: readout is not at k=0 at the echo'
        assert ky == pytest.approx(tr.pe.k_per_m(line), abs=ZERO)
        assert kz == pytest.approx(tr.pe_z.k_per_m(partition), abs=ZERO)


def test_the_centre_partition_encodes_no_z_moment(tr: bSSFP3DTR) -> None:
    assert tr.pe_z.k_per_m(tr.center_partition) == pytest.approx(0.0, abs=1e-12)
    assert tr.center_partition == tr.matrix[2] // 2


# ------------------------------------------------------------------ claim 3: the timing is flat
def test_te_and_tr_do_not_depend_on_the_encoding(tr: bSSFP3DTR) -> None:
    """
    One shared window for the whole table is what buys this.

    Letting each partition take its shortest lobe would make TE a function of ``kz`` -- a
    contrast gradient across the volume that no k-space check would show.
    """
    views = violent_views(tr)
    seq = train(tr, views)
    centres, echoes = rf_centres(seq), echo_times(seq, tr)

    tes = [echoes[i] - centres[i] for i in range(len(views))]
    trs = [centres[i + 1] - centres[i] for i in range(len(centres) - 1)]

    assert max(tes) - min(tes) < 1e-12
    assert max(trs) - min(trs) < 1e-12
    assert tes[0] == pytest.approx(tr.te_s, abs=1e-9)
    assert trs[0] == pytest.approx(tr.tr_s, abs=1e-9)


def test_the_default_is_the_symmetric_realisation(tr: bSSFP3DTR) -> None:
    """``TE = TR/2`` to within the residual the module reports, not to within hand-waving."""
    assert tr.te_s == pytest.approx(tr.tr_s / 2.0, abs=abs(tr.symmetry_residual_s) + 1e-12)
    assert abs(tr.symmetry_residual_s) <= float(tr.opts.grad_raster_time) / 2
    assert tr.te_s == pytest.approx(tr.min_symmetric_te_s)
    assert tr.min_te_s <= tr.min_symmetric_te_s


# --------------------------------------------------- claim 4: the first moment, as realised
def test_the_centre_partition_keeps_the_symmetric_slice_limit(tr: bSSFP3DTR) -> None:
    """
    At ``kz = 0`` the z axis carries only slab-selection balancing, and the interval's first
    moment vanishes.

    **This is a property of the realisation, not of balance.**  Both selection-side lobes take
    one shared window, which is what places them symmetrically about the echo; balance itself
    would be satisfied by any pair of durations carrying the same areas.
    """
    seq = train(tr, [(tr.center_line, tr.center_partition),
                     (0, tr.center_partition + 5), (tr.matrix[1] - 1, tr.center_partition - 3)])
    centres, echoes = rf_centres(seq), echo_times(seq, tr)

    m1z = moment(seq, centres[0], centres[1], 'z', order=1, ref=echoes[0])

    assert abs(m1z) < 1e-9, f'centre partition left M1z = {m1z:.3e}'


def test_the_partition_encoding_contributes_its_own_first_moment(tr: bSSFP3DTR) -> None:
    """
    Away from the centre, ``z`` behaves like ``y``: an encoding-dependent first moment.

    Recorded rather than nulled.  The encode and its rewind sit on opposite sides of the echo,
    so their moments add rather than cancel -- which is what phase encoding does on any axis, and
    is why this is not a defect introduced by going to 3D.
    """
    nz = tr.matrix[2]
    out = {}
    for partition in (0, tr.center_partition, nz - 1):
        seq = train(tr, [(tr.center_line, partition), (0, (partition + 5) % nz),
                         (tr.matrix[1] - 1, (partition + 9) % nz)])
        centres, echoes = rf_centres(seq), echo_times(seq, tr)
        out[partition] = moment(seq, centres[0], centres[1], 'z', order=1, ref=echoes[0])

    assert abs(out[tr.center_partition]) < 1e-9
    assert abs(out[0]) > 1e-3, 'the low edge must carry a first moment'
    assert abs(out[nz - 1]) > 1e-3, 'and so must the high edge'
    assert out[0] * out[nz - 1] < 0, 'and they sit on opposite sides, as the signed table does'


# ------------------------------------------------------------------------ claim 5: the labels
def test_every_acquired_repetition_addresses_its_own_k_space_point(tr: bSSFP3DTR) -> None:
    """``LIN`` alone cannot address a volume: without ``PAR`` every partition writes one line."""
    views = violent_views(tr)
    seq = train(tr, views)

    labels = []
    for i in range(1, len(seq.block_events) + 1):
        block = seq.get_block(i)
        if getattr(block, 'label', None) is not None:
            labels.append({lab.label: int(lab.value) for lab in np.atleast_1d(block.label)})
    addresses = [(d.get('LIN'), d.get('PAR')) for d in labels if 'LIN' in d or 'PAR' in d]

    assert len(set(addresses)) == len(views), 'each repetition needs a distinct address'
    assert set(addresses) == set(views)


def test_a_non_acquiring_repetition_sets_no_labels_and_keeps_its_gradients(tr: bSSFP3DTR) -> None:
    """
    ``acquire=False`` drops the ADC and the labels and changes nothing else.

    A start-up repetition that loaded different gradients would not be preparing the steady
    state this module's train reaches.
    """
    acquired = tr.build(line=7, partition=3, acquire=True)
    quiet = tr.build(line=7, partition=3, acquire=False)

    assert quiet.duration == pytest.approx(acquired.duration)
    for axis in 'xyz':
        assert sc.moments(quiet).get(axis, 0.0) == pytest.approx(
            sc.moments(acquired).get(axis, 0.0), abs=1e-12)
    assert sc.moments(quiet, order=1).get('z', 0.0) == pytest.approx(
        sc.moments(acquired, order=1).get('z', 0.0), abs=1e-12)


# ------------------------------------------------------------------- claim 6: phase is caller's
def test_the_carrier_and_the_receiver_take_the_same_phase(tr: bSSFP3DTR) -> None:
    """One value reaches both, so a phase-cycled train cannot have them disagree."""
    for phase_deg in (0.0, 90.0, 180.0):
        events = [ev for _, ev, _ in sc.flatten(tr.build(line=4, partition=2,
                                                         phase_deg=phase_deg))]
        rf = next(ev for ev in events if getattr(ev, 'type', None) == 'rf')
        adc = next(ev for ev in events if getattr(ev, 'type', None) == 'adc')

        assert np.rad2deg(float(rf.phase_offset)) == pytest.approx(phase_deg, abs=1e-9)
        assert np.rad2deg(float(adc.phase_offset)) == pytest.approx(phase_deg, abs=1e-9)


# ------------------------------------------------------------- claim 7: the excitation mode
def _excitation(module: bSSFP3DTR) -> tuple[int, set[str]]:
    """``(RF sample count, gradient channels)`` of the excitation it actually plays."""
    events = [ev for _, ev, _ in sc.flatten(module.exc(phase_deg=0.0))]
    rf = next(ev for ev in events if getattr(ev, 'type', None) == 'rf')
    channels = {ev.channel for ev in events if getattr(ev, 'type', None) in {'grad', 'trap'}}
    return int(np.asarray(rf.signal).size), channels


def test_a_slab_and_no_slab_are_two_excitation_modes(opts: pp.Opts) -> None:
    """
    Not one mode with the selection gradient switched off: the RF family differs too.

    The rule is :func:`~seqcraft.modules._support.resolve_excitation_mode`, shared with
    ``GRE3DTR`` so the two 3D kernels cannot answer it differently.
    """
    slab_n, slab_ch = _excitation(
        bSSFP3DTR(opts=opts, fov_mm=FOV, matrix=MATRIX, slab_thickness_mm=SLAB_MM))
    hard_n, hard_ch = _excitation(
        bSSFP3DTR(opts=opts, fov_mm=FOV, matrix=MATRIX, slab_thickness_mm=None))

    assert slab_ch == {'z'} and hard_ch == set()
    assert slab_n > 100 and hard_n <= 4, 'the RF families differ, which is what the mode is'


def test_an_explicit_shaped_pulse_without_a_slab_is_still_available(opts: pp.Opts) -> None:
    """The expert path: spatially non-selective, spectrally shaped, asked for deliberately."""
    n, channels = _excitation(bSSFP3DTR(opts=opts, fov_mm=FOV, matrix=MATRIX,
                                        slab_thickness_mm=None, rf_pulse='sinc'))

    assert channels == set() and n > 100


def test_a_non_selective_repetition_is_still_balanced(opts: pp.Opts) -> None:
    """``A_pre`` and ``A_post`` are zero, so the z axis carries the encoding and nothing else."""
    module = bSSFP3DTR(opts=opts, fov_mm=FOV, matrix=MATRIX, slab_thickness_mm=None)
    seq = train(module, violent_views(module))
    centres = rf_centres(seq)

    for k in range(len(centres) - 1):
        for axis in 'xyz':
            assert abs(moment(seq, centres[k], centres[k + 1], axis)) < ZERO


# ---------------------------------------------------------------------------- claim 8: refusals
def test_a_partition_outside_the_table_is_refused(tr: bSSFP3DTR) -> None:
    for partition in (-1, MATRIX[2]):
        with pytest.raises(ConfigurationError, match='partition'):
            tr.build(line=0, partition=partition)


def test_a_two_component_protocol_is_refused(opts: pp.Opts) -> None:
    """The 2D protocol belongs to ``bSSFP2DTR``, and the refusal says so."""
    with pytest.raises(ConfigurationError, match='three components'):
        bSSFP3DTR(opts=opts, fov_mm=(200.0, 200.0), matrix=MATRIX)  # type: ignore[arg-type]


def test_a_time_bandwidth_product_is_refused_on_the_hard_default(opts: pp.Opts) -> None:
    """Inherited from the shared resolver rather than re-derived here."""
    with pytest.raises(ConfigurationError, match='block pulse'):
        bSSFP3DTR(opts=opts, fov_mm=FOV, matrix=MATRIX, slab_thickness_mm=None,
                  rf_time_bw_product=4.0)


def test_an_echo_time_below_the_minimum_names_what_is_reachable(tr: bSSFP3DTR) -> None:
    with pytest.raises(ConfigurationError, match='min_te_s|shorter than'):
        bSSFP3DTR(opts=tr.opts, fov_mm=FOV, matrix=MATRIX, slab_thickness_mm=SLAB_MM,
                  te_s=tr.min_te_s / 2)
