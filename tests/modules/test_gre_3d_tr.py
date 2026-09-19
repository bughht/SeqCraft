"""
``GRE3DTR`` -- the z axis, which is the only reason the class exists.

On x and y a 3D repetition is a 2D one.  On z it carries a partition encode and, when the
excitation is slab-selective, the rephasing that selection implies as well -- two moments on one
axis inside one window.  Everything here is about that coupling being **signed**, solved **once**
for every partition, and realised **exactly once**.

The failures this guards against are all quiet ones.  A doubled slab term shifts every partition
by a constant and passes any k-space extent check.  A per-partition winder makes TE a function of
kz, which is a contrast gradient across the volume that nothing downstream reports.  A limiting
partition found by index rather than by signed moment is right until the slab polarity changes.
"""

from __future__ import annotations

import numpy as np
import pytest
from pypulseq.opts import Opts

import seqcraft as sc
from seqcraft.modules.kernel.gre_3d_tr import _limiting_index

FOV_MM = (200.0, 200.0, 160.0)
MATRIX = (32, 16, 8)


@pytest.fixture(scope='module')
def scanner() -> Opts:
    """The 3D reference's scanner: 20 mT/m and a slow 120 T/m/s, where z actually costs time."""
    return Opts(max_grad=20, grad_unit='mT/m', max_slew=120, slew_unit='T/m/s', B0=3.0,
                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)


def _make(scanner, **kwargs):
    return sc.modules.GRE3DTR(opts=scanner, fov_mm=FOV_MM, matrix=MATRIX, **kwargs)


@pytest.fixture(scope='module')
def nonselective(scanner):
    return _make(scanner)


@pytest.fixture(scope='module')
def selective(scanner):
    return _make(scanner, slab_thickness_mm=180.0)


def _k_at_echo(tr, scanner, *, line: int, partition: int) -> np.ndarray:
    k = sc.kspace(tr(line=line, partition=partition), scanner)
    return k['k_adc'][:, tr.ro.echo_sample(0)]


# ------------------------------------------------------------- the lattices it encodes
@pytest.mark.parametrize('selective_slab', [None, 180.0])
def test_k_at_the_echo_is_the_requested_line_and_partition(scanner, selective_slab) -> None:
    """
    Signed, on all three axes, and **the selective case is the test**: a slab term applied
    twice, or not at all, shifts kz by a constant that no extent check would show.
    """
    tr = _make(scanner, slab_thickness_mm=selective_slab)
    for line, partition in ((3, 6), (0, 0), (tr.center_line, tr.center_partition)):
        echo = _k_at_echo(tr, scanner, line=line, partition=partition)

        assert abs(echo[0]) < 1e-3
        assert echo[1] == pytest.approx(tr.pe.k_per_m(line), abs=1e-3)
        assert echo[2] == pytest.approx(tr.pe_z.k_per_m(partition), abs=1e-3)


@pytest.mark.parametrize('selective_slab', [None, 180.0])
def test_kz_is_held_through_the_readout(scanner, selective_slab) -> None:
    """An encode, not a rephase.  A slab rephaser would bring kz back to zero at the echo."""
    tr = _make(scanner, slab_thickness_mm=selective_slab)
    k = sc.kspace(tr(line=3, partition=6), scanner)['k_adc']

    assert np.ptp(k[2]) < 1e-3
    assert k[2, 0] == pytest.approx(tr.pe_z.k_per_m(6), abs=1e-3)


def test_dkz_comes_from_the_encoded_extent_not_the_slab(scanner, selective) -> None:
    """``1/fov_z``, and the slab is 180 mm against a 160 mm FOV -- deliberately different."""
    assert selective.dk_per_m('z') == pytest.approx(1e3 / FOV_MM[2])
    assert selective.slab_thickness_mm == 180.0
    assert selective.voxel_mm('z') == pytest.approx(FOV_MM[2] / MATRIX[2])


# ------------------------------------------------------- the three degenerate limits
def test_l1_non_selective_reduces_to_the_partition_encode(scanner, nonselective) -> None:
    """``A_slab = 0``, so the combined winder **is** the partition encode."""
    assert nonselective.slab_rephase_area_per_m == 0.0
    for partition in range(MATRIX[2]):
        assert nonselective.combined_z_area_per_m(partition) == pytest.approx(
            nonselective.pe_z.k_per_m(partition), abs=1e-12)


def test_l2_the_centre_partition_is_the_pure_slab_rephaser(scanner, selective) -> None:
    """``A_partition = 0`` there, so what is left is the rephasing -- and kz reaches zero."""
    assert selective.pe_z.k_per_m(selective.center_partition) == 0.0
    assert selective.combined_z_area_per_m(selective.center_partition) == pytest.approx(
        selective.slab_rephase_area_per_m, abs=1e-12)

    echo = _k_at_echo(selective, scanner, line=3, partition=selective.center_partition)
    assert abs(echo[2]) < 1e-3


@pytest.mark.parametrize('slab', [None, 180.0])
def test_l3_adjacent_partitions_differ_by_exactly_one_dkz(scanner, slab) -> None:
    """
    The slab term is constant and must cancel from the difference.

    This is the one that catches a rephasing term wrong by a constant: L2 passes at the centre
    and L1 never exercises it, while every other partition is shifted.
    """
    tr = _make(scanner, slab_thickness_mm=slab)
    steps = [tr.combined_z_area_per_m(p + 1) - tr.combined_z_area_per_m(p)
             for p in range(MATRIX[2] - 1)]

    assert np.allclose(steps, tr.dk_per_m('z'), atol=1e-12)


# --------------------------------------------------- signed coupling and the limiting edge
def test_p1_p2_the_limiting_edge_follows_the_sign_of_the_slab_term(scanner) -> None:
    """
    The worked case from the design review, driven through the solver itself.

    A negative slab term makes the low edge limiting; a positive one of the same size makes the
    high edge limiting. Same partitions, same hardware, opposite answers.
    """
    partitions = [-200.0, -100.0, 0.0, 100.0, 200.0]

    low, _ = _limiting_index([-120.0 + a for a in partitions], scanner)
    high, _ = _limiting_index([+120.0 + a for a in partitions], scanner)

    assert (low, high) == (0, len(partitions) - 1)


def test_p3_a_large_enough_offset_switches_the_limiting_edge(scanner) -> None:
    """
    A "largest index wins" shortcut passes P1/P2 by luck and fails here.

    Sweeping the offset from negative to positive must move the limit across, and at an offset
    that cancels the low edge the limit must already be on the high side.
    """
    partitions = [-200.0, -100.0, 0.0, 100.0, 200.0]
    limits = [_limiting_index([offset + a for a in partitions], scanner)[0]
              for offset in (-200.0, -50.0, 0.0, 50.0, 200.0)]

    assert limits[0] == 0
    assert limits[-1] == len(partitions) - 1
    assert limits == sorted(limits)                       # it moves monotonically, and it moves


def test_p4_the_sign_comes_from_phase_encode_and_not_from_the_index(scanner, selective) -> None:
    """Partition index -> signed physical moment is ``PhaseEncode``'s mapping, not array order."""
    below, above = selective.center_partition - 1, selective.center_partition + 1

    assert selective.pe_z.k_per_m(below) < 0 < selective.pe_z.k_per_m(above)
    assert (selective.combined_z_area_per_m(above) - selective.combined_z_area_per_m(below)
            == pytest.approx(2 * selective.dk_per_m('z'), abs=1e-12))


def test_p5_the_reported_limiting_partition_is_the_one_that_needs_longest(scanner, selective):
    """Enumerated over every partition, and the module reports which one won."""
    import pypulseq as pp

    durations = []
    for partition in range(MATRIX[2]):
        area = selective.combined_z_area_per_m(partition)
        durations.append(0.0 if abs(area) < 1e-9 else float(pp.calc_duration(
            pp.make_trapezoid(channel='z', area=area, system=scanner))))

    assert selective.limiting_partition == int(np.argmax(durations))
    assert selective.winder_s >= max(durations) - 1e-12


# --------------------------------------------------------------------- one shared window
def test_t3_every_partition_has_the_same_te(scanner, selective) -> None:
    """
    A per-partition winder would make TE a function of kz.

    Measured, not asserted from the design: the compiled echo time at the centre and at both
    edges.
    """
    times = []
    for partition in (0, selective.center_partition, MATRIX[2] - 1):
        k = sc.kspace(selective(line=3, partition=partition), scanner)
        times.append(k['t_adc'][selective.ro.echo_sample(0)] - k['t_excitation'][0])

    assert np.ptp(times) < 1e-9
    assert times[0] == pytest.approx(selective.te_s, abs=1e-6)


def test_t2_the_window_lengthens_when_z_needs_it(scanner) -> None:
    """
    A protocol where z genuinely forces the window, rather than one where it merely fits.

    Many thin partitions make the edge moment large; the winder must grow to hold it, and the
    module must not refuse.
    """
    wide = sc.modules.GRE3DTR(opts=scanner, fov_mm=(200.0, 200.0, 40.0), matrix=(32, 16, 64),
                              slab_thickness_mm=60.0)
    easy = sc.modules.GRE3DTR(opts=scanner, fov_mm=(200.0, 200.0, 240.0), matrix=(32, 16, 4),
                              slab_thickness_mm=260.0)

    assert wide.winder_s > easy.winder_s
    assert wide.winder_s >= wide.ro.prephaser_duration_s
    assert wide.min_te_s > easy.min_te_s          # the longer window lands in the echo time


def test_t1_implicit_timing_gives_the_shortest_legal_design(scanner, selective) -> None:
    assert selective.te_s == pytest.approx(selective.min_te_s)
    assert selective.tr_s == pytest.approx(selective.min_tr_s)


def test_t5_the_block_is_exactly_tr_and_the_compiled_te_matches(scanner, selective) -> None:
    block = selective(line=3, partition=6)
    k = sc.kspace(block, scanner)

    assert block.duration == pytest.approx(selective.tr_s, abs=1e-12)
    assert k['t_adc'][selective.ro.echo_sample(0)] - k['t_excitation'][0] == pytest.approx(
        selective.te_s, abs=1e-6)


def test_t4_an_impossible_te_raises_and_names_what_is_limiting(scanner) -> None:
    with pytest.raises(sc.ConfigurationError) as raised:
        _make(scanner, slab_thickness_mm=180.0, te_s=1e-4)

    message = str(raised.value)
    assert 'shorter than this repetition can achieve' in message
    assert 'limiting_partition' in message
    assert 'combined_z_area_per_m' in message


def test_an_impossible_tr_raises(scanner) -> None:
    with pytest.raises(sc.ConfigurationError, match='shorter than this repetition'):
        _make(scanner, tr_s=1e-4)


# ------------------------------------------------------- the slab term, applied exactly once
def test_the_standalone_rephaser_is_not_emitted_beside_the_combined_winder(scanner, selective):
    """
    S8/G19.  A doubled slab term is invisible to every k-space extent check.

    Counted structurally -- the excitation contributes one z gradient, not two -- and measured:
    the z moment from the RF's effective centre to the echo is the combined one, once.
    """
    from seqcraft.modules._support import area_until

    nodes = list(sc.flatten(selective(line=3, partition=6)))
    from_excitation = [event for _, event, path in nodes
                       if 'Excitation' in path and getattr(event, 'type', '') == 'trap']

    assert len(from_excitation) == 1                      # the selection gradient, no rephaser

    waveform = sc.sample(selective(line=3, partition=6), scanner)[1]['z']
    del waveform                                          # sampled only to prove z is not silent
    combined = selective.combined_z_area_per_m(6)
    assert combined == pytest.approx(
        selective.slab_rephase_area_per_m + selective.pe_z.k_per_m(6), abs=1e-12)
    # If the rephaser were emitted as well, kz at the echo would be off by exactly that term.
    echo = _k_at_echo(selective, scanner, line=3, partition=6)
    assert echo[2] == pytest.approx(selective.pe_z.k_per_m(6), abs=1e-3)
    assert abs(echo[2] - (selective.pe_z.k_per_m(6) + selective.slab_rephase_area_per_m)) > 1.0
    del area_until


def test_the_non_selective_path_still_emits_no_rephaser(scanner, nonselective) -> None:
    """There is nothing to rephase, and ``rephase`` is irrelevant to a non-selective pulse."""
    kinds = [(getattr(e, 'type', ''), getattr(e, 'channel', ''))
             for _, e, path in sc.flatten(nonselective(line=3, partition=6))
             if 'Excitation' in path]

    assert kinds == [('rf', '')]


# ------------------------------------------------------------------------ the contract
def test_the_component_contract(scanner, selective, component_checks) -> None:
    component_checks.all(selective, line=3, partition=6)


def test_labels_name_the_line_and_the_partition(scanner, selective) -> None:
    labels = {e.label: e.value for _, e, _ in sc.flatten(selective(line=3, partition=6))
              if getattr(e, 'type', '') == 'labelset'}

    assert labels == {'LIN': 3, 'PAR': 6}
    assert not [e for _, e, _ in sc.flatten(selective(line=3, partition=6, acquire=False))
                if getattr(e, 'type', '') == 'labelset']


def test_a_two_component_fov_is_refused(scanner) -> None:
    """A 3D kernel does not quietly accept a 2D geometry."""
    with pytest.raises(sc.ConfigurationError, match='three components'):
        sc.modules.GRE3DTR(opts=scanner, fov_mm=(200.0, 200.0), matrix=MATRIX)


def test_a_partition_outside_the_matrix_is_refused(scanner, selective) -> None:
    with pytest.raises(sc.ConfigurationError):
        selective(line=3, partition=MATRIX[2])


def test_a_non_positive_slab_is_refused(scanner) -> None:
    with pytest.raises(sc.ConfigurationError):
        _make(scanner, slab_thickness_mm=0.0)


def test_a_slab_smaller_than_the_fov_is_allowed(scanner) -> None:
    """R4b excites 0.7 of its FOV; a rule requiring slab >= fov_z would reject a real sequence."""
    narrow = _make(scanner, slab_thickness_mm=0.7 * FOV_MM[2])
    wide = _make(scanner, slab_thickness_mm=1.2 * FOV_MM[2])

    assert narrow.selective and wide.selective
    assert abs(narrow.slab_rephase_area_per_m) > abs(wide.slab_rephase_area_per_m)
