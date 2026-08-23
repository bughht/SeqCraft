"""
``examples/megre_2d/01_build.ipynb``, asserted from outside it.

This directory's example defines **no class**, which no other example directory can say: a
multi-echo gradient echo is ``GRE2D`` with two more arguments, so there is nothing here for a test
to import and nothing for a notebook to have got wrong about a composition.  What there *is* to
check is the four ``.seq`` files it writes, and specifically the claims whose failure is **silent**
-- the ones that compile, pass every k-space extent check, and produce a plausible image:

- **Every readout-axis gradient in every file is a ``trap``.**  The train states its block
  boundaries with ``sc.barrier``; a split lobe is legal, correct, simulates fine, and is reported
  only by a merge warning.
- **Signed ``k`` at the ``k = 0`` sample of every echo.**  ``|k|`` is symmetric, so a k-space
  extent check passes on a mirrored echo -- and a mirrored echo against a symmetric phantom looks
  entirely correct.  The sample is the one ``EchoSamples`` names, read out of the file.
- **``ECO`` runs ``0 … 7`` and ``REV`` alternates or is constant**, read off the *compiled*
  sequence rather than off the tree.
- **The ``te_s`` in each file's ``[DEFINITIONS]`` matches the module's**, to a nanosecond.  ``02``
  fits against those numbers; if they are not the sequence's own echo times, every map it produces
  is scaled and nothing in either notebook would notice.
- **Flip, TR, matrix and slice are identical across the imaging files**, which is the pin that
  keeps ``02``'s comparison a comparison of readouts rather than of contrasts.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc

nbformat = pytest.importorskip('nbformat', reason='needs seqcraft[dev]')

EXAMPLES = Path(__file__).resolve().parents[2] / 'examples'
NOTEBOOK = EXAMPLES / 'megre_2d' / '01_build.ipynb'

#: The three files ``01`` §7 writes.  Two polarities and a partial echo, on one protocol.
FILES = ('megre_mono', 'megre_bipolar', 'megre_mono_pf075')


def _run(notebook: Path, tmp_path_factory, stop_at: str) -> tuple[dict, Path]:
    """Execute a notebook's code cells up to and including the one containing `stop_at`."""
    if not notebook.exists():                                       # pragma: no cover
        pytest.skip(f'{notebook} is not present')
    cells = [cell.source for cell in nbformat.read(notebook, as_version=4).cells
             if cell.cell_type == 'code']
    end = next(i for i, source in enumerate(cells) if stop_at in source)
    # Skip the cells that *draw*, not every cell that mentions matplotlib: `01`'s setup cell sets
    # `plt.rcParams` and also defines the protocol, so filtering on `plt.` alone would drop the
    # constants every later cell needs.  Every figure here opens with `plt.subplots`.
    sources = [s for s in cells[: end + 1]
               if 'plt.subplots(' not in s and 'sc.plot_block' not in s]

    namespace: dict = {'__name__': '__notebook__'}
    scratch = tmp_path_factory.mktemp('megre_2d')
    here = os.getcwd()
    os.chdir(scratch)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', sc.SeqCraftWarning)
            for index, source in enumerate(sources):
                exec(compile(source, f'{notebook.name}:{index}', 'exec'), namespace)  # noqa: S102
    finally:
        os.chdir(here)
    return namespace, scratch / 'seq'


@pytest.fixture(scope='module')
def built(tmp_path_factory) -> tuple[dict, Path]:
    """``01`` run to the end of section 7, and the directory it wrote its files into."""
    return _run(NOTEBOOK, tmp_path_factory,
                stop_at='WRITTEN = {name: write(name, module) for name, module in SCANS.items()}')


@pytest.fixture(scope='module')
def notebook(built) -> dict:
    return built[0]


@pytest.fixture(scope='module')
def seq_dir(built) -> Path:
    return built[1]


def _definitions(path: Path) -> dict:
    """The ``[DEFINITIONS]`` block, parsed the way ``02`` parses it and for the same reason."""
    out, inside = {}, False
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line.startswith('['):
            inside = line == '[DEFINITIONS]'
            continue
        if not inside or not line:
            continue
        key, _, rest = line.partition(' ')
        try:
            values = [float(p) for p in rest.split()]
        except ValueError:
            out[key] = rest.strip()
            continue
        out[key] = values[0] if len(values) == 1 else values
    return out


@pytest.fixture(scope='module')
def defs(seq_dir) -> dict:
    return {name: _definitions(seq_dir / f'{name}.seq') for name in FILES}


# ------------------------------------------------------------------ the files exist
def test_all_three_files_are_written(seq_dir) -> None:
    missing = [name for name in FILES if not (seq_dir / f'{name}.seq').exists()]
    assert not missing, missing


# ------------------------------------------------------------------ the trapezoids
@pytest.mark.parametrize('name', FILES)
def test_every_readout_gradient_survives_as_one_trapezoid(notebook, name) -> None:
    """
    A split lobe compiles, simulates and images correctly, and is reported only by a merge
    warning naming the readout axis.  ``01`` §3 measures that this train does not need the
    barriers to avoid it; this asserts that the files it actually wrote did not.
    """
    tree = notebook['SCANS'][name](lines=range(notebook['NY']), dummies=notebook['DUMMIES'])

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', sc.SeqCraftWarning)
        seq = sc.compile(tree, notebook['opts'], name=name)

    kinds = {getattr(g, 'type', '?') for i in range(1, len(seq.block_events) + 1)
             for g in (getattr(seq.get_block(i), 'gx', None),) if g is not None}
    assert kinds == {'trap'}, kinds
    merges = [str(w.message) for w in caught if 'merge' in str(w.message)]
    assert not merges, merges


# --------------------------------------------------------------------------- k = 0
def _k_at_echoes(readout, opts) -> np.ndarray:
    """Signed ``k`` along the readout axis at the ``k = 0`` sample of every echo."""
    axis = 'xyz'.index(readout.axis)
    k = sc.kspace(sc.LogicBlock('probe').add(0.0, readout()), opts)['k_adc'][axis]
    rows = k.reshape(readout.echoes, readout.num_samples)
    return np.array([rows[n][readout.echo_sample(n)] for n in range(readout.echoes)])


@pytest.mark.parametrize('name', FILES)
def test_k_is_zero_at_the_echo_sample_of_every_echo(notebook, name) -> None:
    """Signed, at the sample ``echo_sample`` names, on every echo of every imaging file."""
    readout = notebook['SCANS'][name].tr.ro
    worst = float(np.abs(_k_at_echoes(readout, notebook['opts'])).max())

    assert worst < 1e-9, f'{name}: {worst:.3e} 1/m'


@pytest.mark.parametrize('name', FILES)
def test_the_reverse_echoes_land_on_the_forward_grid(notebook, name) -> None:
    """Rule 1 on the compiled waveform.  Vacuous under monopolar, which is why it says so."""
    readout = notebook['SCANS'][name].tr.ro
    axis = 'xyz'.index(readout.axis)
    k = sc.kspace(sc.LogicBlock('probe').add(0.0, readout()), notebook['opts'])['k_adc'][axis]
    rows = k.reshape(readout.echoes, readout.num_samples)

    grid = rows[0]
    for echo in range(readout.echoes):
        got = rows[echo] if readout.polarity_of(echo) > 0 else rows[echo][::-1]
        assert np.abs(got - grid).max() < 1e-9, f'{name}: echo {echo} is on a different grid'


# ------------------------------------------------------------------------ the labels
@pytest.mark.parametrize('name', FILES)
def test_eco_and_rev_are_in_the_compiled_sequence(notebook, name) -> None:
    """
    Read off the compiled sequence rather than off the tree, because that is what a
    reconstruction sees.  ``ECO`` runs ``0 … 7`` inside every repetition; ``REV`` alternates under
    bipolar and is constantly zero under monopolar -- and the monopolar zeros are **not**
    omissible, because pulseq labels are stateful and an omitted ``SET`` inherits.
    """
    readout = notebook['SCANS'][name].tr.ro
    echoes = readout.echoes
    tree = notebook['SCANS'][name](lines=[notebook['CENTER_LINE']], dummies=0)
    seq = sc.compile(tree, notebook['opts'], name=name)

    seen = []
    for index in range(1, len(seq.block_events) + 1):
        label = getattr(seq.get_block(index), 'label', None)
        if label is None:
            continue
        events = label if isinstance(label, (list, tuple)) else [label]
        seen.append({e.label: e.value for e in events})

    assert [d['ECO'] for d in seen if 'ECO' in d] == list(range(echoes))
    assert [d['REV'] for d in seen if 'REV' in d] == [
        int(readout.polarity_of(n) < 0) for n in range(echoes)]
    assert [d['LIN'] for d in seen if 'LIN' in d] == [notebook['CENTER_LINE']], (
        'LIN belongs to the repetition and is emitted once, whatever the readout does'
    )


# ------------------------------------------------------------------- the definitions
@pytest.mark.parametrize('name', FILES)
def test_the_written_echo_times_are_the_modules_own(notebook, defs, name) -> None:
    """
    To a nanosecond.  ``02`` fits against these numbers rather than rebuilding the module, so if
    they are not the sequence's own echo times every map it produces is scaled -- and a scale
    error in a field map is exactly the kind of wrong no image inspection finds.
    """
    scan = notebook['SCANS'][name]
    offset = scan.tr.te_s - scan.tr.ro.te_s[0]
    wanted = np.asarray(scan.tr.ro.te_s) + offset
    written = np.asarray(defs[name]['TE'])

    assert written.shape == wanted.shape
    assert np.abs(written - wanted).max() < 1e-9, (
        f'{name}: worst {np.abs(written - wanted).max() * 1e9:.3f} ns'
    )


@pytest.mark.parametrize('name', FILES)
def test_the_definitions_describe_the_train(defs, name) -> None:
    """`echoes`, `polarity`, the samples and the polarities, so ``02`` never guesses a parity."""
    d = defs[name]
    echoes = int(d['Echoes'])

    assert d['Polarity'] in ('monopolar', 'bipolar')
    assert echoes == 8
    assert len(d['TE']) == echoes
    assert len(d['EchoSamples']) == echoes
    assert len(d['EchoPolarity']) == echoes
    assert set(d['EchoPolarity']) <= {1.0, -1.0}
    if d['Polarity'] == 'monopolar':
        assert set(d['EchoPolarity']) == {1.0}
        assert len(set(d['EchoSamples'])) == 1
    else:
        assert [int(v) for v in d['EchoPolarity']] == [1, -1] * (echoes // 2)


def test_the_files_share_one_protocol(defs) -> None:
    """
    The pin that keeps ``02``'s comparison honest.  Flip, TR, slice and FOV identical across every
    file, so what differs between them is the readout and nothing else.
    """
    shared = {name: (defs[name]['FlipAngle'], defs[name]['TR'], defs[name]['SliceThickness'],
                     tuple(defs[name]['FOV']), int(defs[name]['kSpaceCenterLine']))
              for name in FILES}

    assert len(set(shared.values())) == 1, shared


# --------------------------------------------------------------------- no new class
def test_the_notebook_defines_no_module_class() -> None:
    """
    The finding, asserted rather than described.  Every other example directory defines a
    :class:`~seqcraft.Module` subclass because its sequence is a composition the package does not
    have; a multi-echo GRE is ``GRE2D`` with two more arguments, so this one does not.
    """
    if not NOTEBOOK.exists():                                       # pragma: no cover
        pytest.skip(f'{NOTEBOOK} is not present')
    sources = '\n'.join(cell.source for cell in nbformat.read(NOTEBOOK, as_version=4).cells
                        if cell.cell_type == 'code')

    assert 'sc.Module' not in sources
    assert 'class ' not in sources
    assert 'MEGRE2D' not in sources
    assert not hasattr(sc.modules, 'MEGRE2D')


def test_the_protocol_is_one_constructor_call(notebook) -> None:
    """And the three files that are the protocol are three calls that differ in two arguments."""
    for name in FILES:
        scan = notebook['SCANS'][name]
        assert isinstance(scan, sc.modules.GRE2D)
        assert scan.tr.ro.echoes == 8
        assert scan.tr.ro.polarity in ('monopolar', 'bipolar')


def test_the_flyback_is_the_whole_lobe_in_the_written_protocol(notebook) -> None:
    """``01`` §2's number, checked against the module that wrote the file rather than the print."""
    from seqcraft.modules._support import area_until

    readout = notebook['SCANS']['megre_mono'].tr.ro
    total = area_until(readout.gx, float(pp.calc_duration(readout.gx)))

    assert readout.flyback_area_per_m == pytest.approx(-total, abs=1e-9)
    assert readout.flyback_area_per_m != pytest.approx(-readout.area_to_echo_per_m, rel=0.1)
