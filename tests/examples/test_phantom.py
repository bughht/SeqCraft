"""
The shared example phantom, tested because seven notebooks now depend on it agreeing with itself.

``examples/phantom.py`` is not part of the library -- it is the BrainWeb slab every 2D example
simulates against, prepared in one place instead of seven.  That makes it exactly the kind of thing
that breaks quietly: nothing imports it except notebooks, notebooks are not run by the fast test
tier, and a phantom that changes shape or orientation makes every notebook's numbers wrong without
making any of them fail.

What is worth pinning is the part the notebooks *rely on* and cannot check themselves:

- **the orientation convention**, because an image is the phantom transposed and a notebook that
  guesses draws every picture sideways -- which is what the first draft of the EPI examples did;
- **the slab is centred on isocentre**, because a slab left where it was in the volume is excited
  by nothing and the symptom is a dark image rather than an error;
- **``field_hz=None`` leaves the phantom's own B0 alone**, because every example simulates against
  it and silently replacing it would change what those notebooks measure;
- **``field_hz()`` returns that same map**, because the EPI notebooks draw it beside the image it
  distorted, and a map that does not match the simulation is worse than no map;
- **``slab_3d`` covers an extent where ``slab`` counts slices**, because the 3D examples assert
  their phantom against the sequence's encoded FOV, and those assertions are only worth making if
  the helper's physical extent really is independent of the matrix it is sampled on.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

EXAMPLES = Path(__file__).resolve().parents[2] / 'examples'


@pytest.fixture(scope='module')
def ph():
    """
    Import ``examples/phantom.py`` by path -- it is not on the package path, deliberately.

    The skip is here rather than at :func:`fetched` because ``phantom.py`` imports ``torch`` at
    module scope, so the import itself is what the fast tier cannot do.  ``MRzeroCore`` it imports
    lazily, inside the one function that needs it, which is why that guard stays where it is.
    """
    pytest.importorskip('torch', reason='needs seqcraft[sim]')
    spec = importlib.util.spec_from_file_location('example_phantom', EXAMPLES / 'phantom.py')
    if spec is None or spec.loader is None:            # pragma: no cover
        pytest.skip('examples/phantom.py is not present')
    module = importlib.util.module_from_spec(spec)
    sys.modules['example_phantom'] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope='module')
def fetched(ph):
    """Skip rather than download: the fast tier does not pull 19 MB."""
    pytest.importorskip('MRzeroCore', reason='needs seqcraft[sim]')
    if not ph.PHANTOM_PATH.exists():                   # pragma: no cover
        pytest.skip(f'{ph.PHANTOM_PATH} is not present -- run examples/gre_2d/02 once')
    return ph


def test_same_frame_is_the_one_transpose(ph) -> None:
    """
    An image is ``[ky, kx]`` and the phantom is ``[x, y]``, so drawing them together needs exactly
    one transpose -- named once here rather than remembered in seven notebooks.
    """
    grid = np.arange(12).reshape(3, 4)
    assert ph.same_frame(grid).shape == (4, 3)
    assert np.array_equal(ph.same_frame(ph.same_frame(grid)), grid)


@pytest.mark.parametrize(('nz', 'centre'), [(1, 64), (2, 64), (4, 64), (2, 50), (4, 58)])
def test_the_slab_is_centred_on_the_slice_it_names(ph, nz, centre) -> None:
    """One expression for the slice range, so no notebook writes a second that disagrees."""
    got = ph.slices_of(nz, centre)

    assert len(got) == nz
    assert got == list(range(got[0], got[0] + nz))
    assert got[0] <= centre <= got[-1] + 1
    assert abs(sum(got) / nz - centre) <= 0.5


def test_the_slab_is_moved_to_isocentre(fetched) -> None:
    """
    ``slices()`` leaves the slab where it was in the volume, tens of millimetres off centre, and a
    sequence exciting at isocentre then selects nothing.  The symptom is a dark image.
    """
    nz = 2
    slab = fetched.slab(matrix=64, nz=nz, n_coils=0)
    step = float(slab.affine[2, 2])

    assert float(slab.affine[2, 3]) == pytest.approx(-0.5 * nz * step)
    assert abs(float(slab.affine[2, 3])) < abs(step) * nz


def test_field_hz_none_leaves_the_phantom_s_own_b0_alone(fetched) -> None:
    """
    Every example simulates against the phantom's own B0.  ``None`` has to mean *do not touch*,
    because a shared helper that quietly replaced it would change their physics without changing a
    line of any of them.
    """
    native = fetched.slab(matrix=64, nz=2, n_coils=0)
    shimmed = fetched.slab(matrix=64, nz=2, n_coils=0, field_hz=0.0)

    assert float(native.B0.abs().max()) > 1.0, 'the phantom is supposed to carry a B0 map'
    assert float(shimmed.B0.abs().max()) == 0.0


def test_a_scalar_field_broadcasts_and_an_array_is_taken_as_given(fetched) -> None:
    imposed = np.zeros((64, 64, 2), np.float32)
    imposed[:, :, 0] = 25.0

    assert float(fetched.slab(matrix=64, nz=2, n_coils=0,
                              field_hz=17.0).B0.mean()) == pytest.approx(17.0)
    got = fetched.slab(matrix=64, nz=2, n_coils=0, field_hz=imposed).B0.numpy()
    assert np.array_equal(got, imposed)


@pytest.mark.parametrize(('nz', 'centre'), [(2, 64), (4, 64), (2, 50)])
def test_field_hz_is_exactly_what_the_simulation_sees(fetched, nz, centre) -> None:
    """
    The EPI notebooks draw this map beside the image it distorted, so it has to *be* the map the
    simulation ran against -- not a second one computed the same way.
    """
    drawn = fetched.field_hz(matrix=64, nz=nz, z_center=centre)
    simulated = fetched.slab(matrix=64, nz=nz, n_coils=0, z_center=centre).B0.numpy().real

    assert drawn.shape == (64, 64, nz)
    assert np.array_equal(drawn, simulated)


def test_the_phantom_s_b0_is_mrzero_s_own_and_is_already_centred(fetched) -> None:
    """
    Two properties the notebooks state and rely on, and neither is obvious from the file name.

    ``subject05.npz`` carries **no** ``B0_map`` -- MRzero's loader generates one from the proton
    density when it is absent, demeaned over that density.  So it arrives centred on zero, which is
    why nothing here shims it, and it is analytic rather than measured, which is why the notebook
    says so before showing a distortion figure.
    """
    field = fetched.field_hz(matrix=128, nz=4, z_center=64)
    tissue = fetched.mask(matrix=128, nz=4, z_center=64)

    assert abs(float(np.median(field[tissue]))) < 5.0, 'the map is supposed to arrive demeaned'
    assert field[tissue].max() > 10.0, 'a map with no structure would make section 6 pointless'

    with np.load(fetched.PHANTOM_PATH) as data:
        assert 'B0_map' not in data, 'the file now ships a B0 map -- the docstring needs updating'


# ------------------------------------------------------------------------ the 3D slab is separate
def test_slab_3d_covers_the_physical_extent_it_was_asked_for(fetched) -> None:
    """
    The property the helper exists for: `z_extent_mm` is physical, `matrix[2]` is sampling.

    ``interpolate`` rescales the affine along with the tensors, so resolution and extent are
    independent -- this pins that they stay independent, because every geometric assertion a 3D
    notebook makes about its phantom rests on it.
    """
    out = fetched.slab_3d(matrix=(48, 32, 8), z_extent_mm=96.0)

    assert fetched.extent_mm(out)[2] == pytest.approx(96.0, abs=fetched.NATIVE_MM)
    assert tuple(out.PD.shape) == (48, 32, 8)


def test_the_two_slab_helpers_are_not_interchangeable(fetched) -> None:
    """
    ``slab`` counts *native slices*; ``slab_3d`` covers an *extent*.  Reading one as the other is
    the mistake this pair of helpers exists to prevent.

    At ``nz = 8`` they differ by almost an order of magnitude in z: eight native slices is 12 mm of
    anatomy, which inside a sequence encoding 96 mm fills one partition and leaves seven empty.
    That looks like a reconstruction fault and is not one, so it is worth a test rather than a
    comment.
    """
    thin = fetched.extent_mm(fetched.slab(matrix=48, nz=8, n_coils=0))[2]
    thick = fetched.extent_mm(fetched.slab_3d(matrix=(48, 32, 8), z_extent_mm=96.0))[2]

    assert thin == pytest.approx(8 * fetched.NATIVE_MM, abs=0.5)
    assert thick > 7 * thin, 'slab_3d must not have become an alias for slab'


@pytest.mark.parametrize('nz', [4, 8, 16])
def test_the_z_sampling_does_not_move_the_phantom(fetched, nz: int) -> None:
    """Coarser sampling of the same slab is a coarser phantom, not a smaller one."""
    out = fetched.slab_3d(matrix=(32, 32, nz), z_extent_mm=96.0)

    assert fetched.extent_mm(out)[2] == pytest.approx(96.0, abs=fetched.NATIVE_MM)
    assert abs(float(out.affine[2, 2])) == pytest.approx(96.0 / nz, abs=fetched.NATIVE_MM / nz)


def test_slab_3d_is_centred_on_isocentre(fetched) -> None:
    """Same reason as :func:`test_the_slab_is_moved_to_isocentre`: an off-centre slab is dark."""
    nz = 8
    out = fetched.slab_3d(matrix=(32, 32, nz), z_extent_mm=96.0)
    step = float(out.affine[2, 2])

    assert float(out.affine[2, 3]) == pytest.approx(-0.5 * nz * step)


def test_the_in_plane_extent_is_the_head_and_not_the_matrix(fetched) -> None:
    """
    In-plane, too, the matrix is sampling.  A 3D notebook asserts its phantom fits inside the
    encoded FOV, and that assertion is only meaningful if this holds.
    """
    coarse = fetched.extent_mm(fetched.slab_3d(matrix=(32, 32, 8), z_extent_mm=96.0))
    fine = fetched.extent_mm(fetched.slab_3d(matrix=(96, 64, 8), z_extent_mm=96.0))

    assert coarse[:2] == pytest.approx(fine[:2], abs=1e-3)
    assert coarse[0] == pytest.approx(128 * fetched.NATIVE_MM, abs=1.0)
