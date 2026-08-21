"""
The two EPI composites that stay in their notebooks, asserted from outside them.

``GREEPI2D`` and ``SEEPI2D`` have one consumer each -- their own build notebook -- so neither ships
and neither can be imported here.  What CI can still do is **run the notebook and assert against
what it defined**, which is the bargain ``test_se_notebooks.py`` already strikes.

The claims worth a test rather than a cell are the ones whose failure is *silent* and whose symptom
is indistinguishable from the artefact the ``02`` notebooks exist to measure:

- **Every readout-axis gradient in every written file is a ``trap``.**  A centred blip straddles
  the only boundary the compiler can cut, and without ``sc.barrier`` at the seam it splits every
  readout lobe instead.  That version compiles, passes every k-space check and simulates
  correctly; the only report of it is a merge warning naming the readout axis.
- **k at the ``k = 0`` sample of every echo, signed, on both axes.**  A wrong dephaser sign is a
  legal sequence that acquires a different part of k-space, and ``|k|`` is symmetric, so a k-space
  *extent* check passes on the mirrored one.
- **The reverse echoes land on the forward grid.**  Off-centre by half a gradient raster and the
  two polarities sample different k grids -- an N/2 ghost the sequence made itself.
- **``ramp_sampling`` on and off traverse the same k extent.**  A ramp-sampled lobe reaching a
  different amount of k changes the field of view rather than the echo spacing, and the image
  comes back the right shape at the wrong scale.
- **The spin echo lands where the kernel says**, through ``sc.kspace``, which is the only oracle
  here that models the refocusing conjugation.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

import numpy as np
import pytest

import seqcraft as sc

nbformat = pytest.importorskip('nbformat', reason='needs seqcraft[dev]')

EXAMPLES = Path(__file__).resolve().parents[2] / 'examples'


def _run(notebook: Path, tmp_path_factory, stop_at: str) -> dict:
    """Execute a notebook's code cells up to and including the one containing `stop_at`."""
    if not notebook.exists():                                       # pragma: no cover
        pytest.skip(f'{notebook} is not present')
    cells = [cell.source for cell in nbformat.read(notebook, as_version=4).cells
             if cell.cell_type == 'code']
    end = next(i for i, source in enumerate(cells) if stop_at in source)
    sources = [s for s in cells[: end + 1] if 'plt.' not in s and 'sc.plot_block' not in s]

    namespace: dict = {'__name__': '__notebook__'}
    here = os.getcwd()
    os.chdir(tmp_path_factory.mktemp(notebook.parent.name))
    try:
        with warnings.catch_warnings():
            # The blip axis carries two half-blips per block by design, so the compiler sums them
            # and says so.  The notebooks pin that warning themselves; here it is noise.
            warnings.simplefilter('ignore', sc.SeqCraftWarning)
            for index, source in enumerate(sources):
                exec(compile(source, f'{notebook.name}:{index}', 'exec'), namespace)  # noqa: S102
    finally:
        os.chdir(here)
    return namespace


@pytest.fixture(scope='module')
def gre(tmp_path_factory) -> dict:
    return _run(EXAMPLES / 'gre_epi_2d' / '01_build.ipynb', tmp_path_factory,
                stop_at='class GREEPI2D(sc.Module):')


@pytest.fixture(scope='module')
def se(tmp_path_factory) -> dict:
    return _run(EXAMPLES / 'se_epi_2d' / '01_build.ipynb', tmp_path_factory,
                stop_at='class SEEPI2D(sc.Module):')


def _gre(notebook: dict, **kwargs):
    """A ``GREEPI2D`` on the notebook's own protocol."""
    base = dict(opts=notebook['opts'], fov_mm=notebook['FOV_MM'], matrix=notebook['MATRIX'],
                thickness_mm=notebook['THICKNESS_MM'], dwell_s=notebook['DWELL_S'],
                oversampling=notebook['OVERSAMPLING'], echoes=notebook['NY'],
                center_echo=notebook['CENTER_LINE'], tr_s=notebook['TR_S'])
    return notebook['GREEPI2D'](**{**base, **kwargs})


def _se(notebook: dict, **kwargs):
    base = dict(opts=notebook['opts'], fov_mm=notebook['FOV_MM'], matrix=notebook['MATRIX'],
                thickness_mm=notebook['THICKNESS_MM'], dwell_s=notebook['DWELL_S'],
                oversampling=notebook['OVERSAMPLING'], echoes=notebook['NY'],
                center_echo=notebook['CENTER_LINE'], te_s=notebook['TE_S'],
                tr_s=notebook['TR_S'])
    return notebook['SEEPI2D'](**{**base, **kwargs})


def _shaped(kernel, notebook, lines):
    """``sc.kspace`` for one shot, reshaped to ``(axis, readout, sample)``.

    Rows are counted from the array rather than from `lines`: a navigated shot carries
    `navigator_echoes` more readouts than it has lines, which is the whole point of them.
    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        k = sc.kspace(sc.LogicBlock('probe').add(0.0, kernel(lines=lines)), notebook['opts'])
    n = kernel.epi.num_samples
    rows = k['t_adc'].size // n
    assert rows == len(lines) + kernel.epi.navigator_echoes
    return (k['k_adc'].reshape(3, rows, n), k['t_adc'].reshape(rows, n),
            float(k['t_excitation'][0]))


def _gradient_kinds(seq, axis):
    return {getattr(g, 'type', '?') for i in range(1, len(seq.block_events) + 1)
            for g in (getattr(seq.get_block(i), axis, None),) if g is not None}


# --------------------------------------------------------------------------- the headline
@pytest.mark.parametrize('extra', [
    pytest.param({}, id='ramp'),
    pytest.param({'ramp_sampling': False}, id='flat-top'),
    pytest.param({'partial_fourier': 0.75}, id='partial-echo'),
    pytest.param({'navigator_echoes': 3}, id='navigators'),
])
def test_k_is_exact_at_the_echo_sample_of_every_echo(gre, extra) -> None:
    """
    Signed, on both axes, at the sample :meth:`echo_sample` names -- which is the assertion that
    the mirror is right.  ``1e-5`` 1/m against a ``dk`` of 4.5 is a factor of 400 000 below one
    sample spacing, and the measured worst case is 1.8e-7.
    """
    ny, centre = gre['NY'], gre['CENTER_LINE']
    kernel = _gre(gre, **extra)
    lines = list(range(ny))
    nav = kernel.epi.navigator_echoes

    k, _, _ = _shaped(kernel, gre, lines)
    k = k[:, nav:] if nav else k
    wanted = (np.asarray(lines) - centre) * (1e3 / gre['FOV_MM'][1])

    at_echo = np.array([[k[axis][n, kernel.epi.echo_sample(n)] for n in range(ny)]
                        for axis in (0, 1)])
    assert np.abs(at_echo[0]).max() < 1e-5, 'k_x: the prephaser and the alternating lobes'
    assert np.abs(at_echo[1] - wanted).max() < 1e-5, 'k_y: the blip-axis prephaser and the blips'


def test_the_reverse_echoes_land_on_the_forward_grid(gre) -> None:
    """
    Rule 1 as a statement about the *compiled* waveform: reversing a reverse echo's samples puts
    them on exactly the forward echo's positions.  Measured worst case 7.1e-12 1/m.

    This is the check the whole centred-window design exists to pass, and it is measured against
    what plays rather than against ``k_read_per_m``, which agreeing with itself would prove
    nothing.
    """
    kernel = _gre(gre)
    lines = list(range(gre['NY']))
    k, _, _ = _shaped(kernel, gre, lines)

    forward = [n for n in range(len(lines)) if kernel.epi.polarity(n) > 0]
    reverse = [n for n in range(len(lines)) if kernel.epi.polarity(n) < 0]
    assert forward and reverse

    grid = k[0][forward[0]]
    assert np.abs(k[0][reverse[0]][::-1] - grid).max() < 1e-4
    assert np.abs(np.array([k[0][n] for n in forward]) - grid).max() < 1e-4
    assert np.abs(np.array([k[0][n][::-1] for n in reverse]) - grid).max() < 1e-4


def test_ramp_sampling_changes_the_echo_spacing_and_not_the_k_extent(gre) -> None:
    """
    A ramp-sampled lobe that traversed a different amount of k would change the **field of view**,
    and the image would come back the right shape at the wrong scale -- which no other check here
    would notice.  Measured: 290.1 against 286.4 1/m, both within one ``dk`` of ``N * dk``.
    """
    ramped, flat = _gre(gre).epi, _gre(gre, ramp_sampling=False).epi
    dk = ramped.dk_read_per_m

    assert ramped.echo_spacing_s < flat.echo_spacing_s * 0.85
    assert abs(np.ptp(ramped.k_read_per_m) - np.ptp(flat.k_read_per_m)) < dk
    assert ramped.readout_amplitude_hz_m > flat.readout_amplitude_hz_m


@pytest.mark.parametrize('extra', [
    pytest.param({}, id='ramp'),
    pytest.param({'ramp_sampling': False}, id='flat-top'),
    pytest.param({'partial_fourier': 0.75}, id='partial-echo'),
    pytest.param({'blip_lines': 4}, id='multi-shot'),
    pytest.param({'navigator_echoes': 3}, id='navigators'),
])
def test_every_readout_gradient_survives_as_one_trapezoid(gre, extra) -> None:
    """
    **The barrier's test.**  Without ``sc.barrier`` at the seam the compiler cuts the midpoint of
    the boundary gap, which is inside the readout lobe's fall ramp, and 63 of 65 ``gx`` events
    become arbitrary waveforms.  Nothing else fails when that happens.
    """
    ny = gre['NY']
    kernel = _gre(gre, echoes=ny // 4 if 'blip_lines' in extra else ny,
                  center_echo=list(range(0, ny, 4)).index(gre['CENTER_LINE'])
                  if 'blip_lines' in extra else gre['CENTER_LINE'], **extra)
    lines = list(range(0, ny, 4)) if 'blip_lines' in extra else list(range(ny))

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', sc.SeqCraftWarning)
        seq = sc.compile(sc.LogicBlock('probe').add(0.0, kernel(lines=lines)), gre['opts'])

    assert _gradient_kinds(seq, 'gx') == {'trap'}
    merges = [str(w.message) for w in caught if 'gradient merge' in str(w.message)]
    axes = {part.split(')')[0] for m in merges for part in m.split('(axis ')[1:]}
    assert axes <= {kernel.epi.blip_axis}, f'a merge on {axes - {kernel.epi.blip_axis}}'


def test_time_to_echo_is_the_true_sample_time_across_the_whole_train(gre) -> None:
    """
    Navigators included, which is what negative echo indices are for: they are not simply
    ``echo * echo_spacing_s`` before imaging echo zero, because the blip-axis prephaser plays
    between the two.  Measured: 0 ps.
    """
    kernel = _gre(gre, navigator_echoes=3)
    lines = list(range(gre['NY']))
    nav = kernel.epi.navigator_echoes
    _, t_adc, _ = _shaped(kernel, gre, lines)

    worst = max(abs(t_adc[row, kernel.epi.echo_sample(echo)]
                    - kernel._train_start_s - kernel.epi.time_to_echo(echo))
                for row, echo in enumerate(range(-nav, len(lines))))
    assert worst < 1e-12
    assert [kernel.epi.polarity(e) for e in range(-nav, 1)] == [1, -1, 1, -1]


def test_te_is_measured_to_the_k_zero_sample_of_the_centre_echo(gre) -> None:
    """
    ``center_echo`` is what makes TE a design-time number, and the table is what decides it.  A
    reverse table reaches the centre one echo earlier *and* on the opposite polarity, so the two
    differ by one echo spacing plus one dwell rather than by one echo spacing.
    """
    ny, centre = gre['NY'], gre['CENTER_LINE']
    lines = list(range(ny))
    forward = _gre(gre)
    backward = _gre(gre, center_echo=lines[::-1].index(centre))

    _, t_adc, t_exc = _shaped(forward, gre, lines)
    sample = forward.epi.echo_sample(forward.center_echo)
    assert t_adc[forward.center_echo, sample] - t_exc == pytest.approx(forward.te_s, abs=1e-9)

    gap = forward.te_s - backward.te_s
    assert gap == pytest.approx(forward.epi.echo_spacing_s + forward.epi.dwell_s, abs=1e-12)


# ------------------------------------------------------------------------- the spin echo
def test_the_spin_echo_lands_where_the_kernel_says(se) -> None:
    """
    Through ``sc.kspace``, which applies pypulseq's own conjugation at ``use='refocusing'`` -- the
    only oracle here that models it.  The two numbers cannot be made equal: the gradient echo
    lives on the ADC raster and the spin echo on the gradient raster.  Measured 6.75 us apart.
    """
    kernel = _se(se)
    lines = list(range(se['NY']))
    k, t_adc, t_exc = _shaped(kernel, se, lines)
    sample = kernel.epi.echo_sample(kernel.center_echo)

    assert t_adc[kernel.center_echo, sample] - t_exc == pytest.approx(kernel.te_s, abs=1e-9)
    assert abs(kernel.spin_echo_s - kernel.te_s) < 2 * float(se['opts'].grad_raster_time)
    assert abs(k[0][kernel.center_echo, sample]) < 1e-5
    assert abs(k[1][kernel.center_echo, sample]) < 1e-5


def test_a_flipped_dephaser_sign_is_a_legal_sequence_that_is_wrong(se) -> None:
    """
    The reason the sign is asserted through ``sc.kspace`` rather than looked at: both wrong signs
    compile, meet every limit, and acquire a different part of k-space.  Against a symmetric
    phantom both still look like an image.
    """
    import pypulseq as pp

    kernel = _se(se)
    lines = list(range(se['NY']))
    sample = kernel.epi.echo_sample(kernel.center_echo)

    for sign_x, sign_y, axis in ((-1, -1, 0), (1, 1, 1)):
        probe = sc.LogicBlock().add([
            [0.0, kernel.exc(phase_deg=90.0)],
            [kernel._deph_start_s,
             sc.events.derive(pp.scale_grad(kernel._deph_x, sign_x)),
             sc.events.derive(pp.scale_grad(
                 kernel._deph_y, sign_y * kernel.epi.k_blip_per_m(lines[0])
                 / float(kernel._deph_y.area)))],
            [kernel._refoc_start_s, kernel.refoc(phase_deg=0.0)],
            [kernel._train_start_s, kernel.epi(lines=lines)],
        ])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', sc.SeqCraftWarning)
            k = sc.kspace(sc.LogicBlock('p').add(0.0, probe), se['opts'])['k_adc']
        k = k.reshape(3, len(lines), kernel.epi.num_samples)
        assert abs(k[axis][kernel.center_echo, sample]) > 100.0, 'the wrong sign was not wrong'


def test_both_notebooks_build_the_same_readout(gre, se) -> None:
    """
    ``EPI2D`` holds no RF and cannot tell a gradient echo from a spin echo, which is the whole
    argument for it being a leaf.  The two notebooks differ in what is played *before* the train
    and in nothing about the train itself.
    """
    one, two = _gre(gre).epi, _se(se).epi

    assert one.echo_spacing_s == two.echo_spacing_s
    assert one.num_samples == two.num_samples
    assert one.guard_s == two.guard_s
    assert np.array_equal(one.k_read_per_m, two.k_read_per_m)
    assert one.prephaser_duration_s > 0          # the gradient echo owns a prephaser
    with pytest.raises(sc.ConfigurationError):   # the spin echo's moved before the 180
        _ = two.prephaser_duration_s


def test_the_kernel_gives_one_carrier_phase_to_the_transmitter_and_the_receiver(gre) -> None:
    """`GREEPI2D` is what makes RF spoiling reach an EPI, and it has to reach both halves.

    Spoiling the transmitter alone is worse than not spoiling at all: the schedule's quadratic
    phase lands in the shot rather than in the carrier, which under segmentation is a different
    constant per shot -- a modulation of `ky` with the period of the shot plan.
    """
    kernel = _gre(gre, echoes=8, center_echo=4)
    wanted = float(np.deg2rad(117.0))

    events = list(sc.flatten(kernel(lines=range(8), phase_deg=117.0)))
    rf = [e for _, e, _ in events if getattr(e, 'type', '') == 'rf']
    adcs = [e for _, e, _ in events if getattr(e, 'type', '') == 'adc']

    assert len(rf) == 1
    assert float(rf[0].phase_offset) == pytest.approx(wanted, abs=1e-12)
    assert len(adcs) == 8
    assert [float(a.phase_offset) for a in adcs] == pytest.approx([wanted] * 8, abs=1e-12)


def test_the_spoiling_schedule_is_counted_from_the_first_dummy(gre) -> None:
    """Restarting it after the run-in is a discontinuity where the run-in exists to remove one."""
    schedule = [gre['spoil_phase'](n) for n in range(4)]

    assert schedule == pytest.approx([0.0, 117.0, 351.0, 702.0])
    assert gre['spoil_phase'](2, 0.0) == 0.0            # the naive file's fixed carrier
