"""
The BrainWeb phantom every 2D example simulates against, prepared **one way**.

Five notebooks had grown five copies of the same twenty lines -- load, interpolate, slice, move the
slab to isocentre, build a coil ring.  The five agreed, which is the good case and not the stable
one: the moment two more examples were added they came with a **different download URL**, and
nothing would have noticed until the phantom moved and one pair of notebooks fetched a different
object from the other.  A phantom that differs between notebooks makes their numbers incomparable,
and that is the one thing a shared reference object exists to prevent.

So the preparation lives here and the *parameters* stay in the notebooks, because those are real
choices with reasons attached: ``gre_2d/02`` wants four slices, ``fse_2d/02`` wants two because it
runs three acquisitions, and the EPI notebooks want two for the same reason.

Orientation, once, because it is the thing everyone gets wrong
--------------------------------------------------------------
The phantom is indexed ``[x, y, z]``: **x is left-right**, **y is anterior-posterior with the
frontal lobes at high y**, and z is the slice axis and the B0 direction.  Read that off the
phantom rather than assumed -- the mid-sagittal slice puts the cerebellum and the splenium of the
corpus callosum at *low* y and the genu at *high* y, and the inferior slices, which are cerebellum
and brainstem, sit a good ten voxels posterior of the whole brain's centroid.

A Cartesian reconstruction with the readout on ``x`` produces ``image[ky, kx]``: rows are the
phantom's ``y`` and columns its ``x``, verified by putting a spin at a known off-centre position
and finding where it lands.  So an image needs **no transform at all** -- draw it with
``imshow(image, origin='lower')`` and it is a conventional axial slice, anterior up.  A
phantom-indexed map is ``[x, y]`` and needs one transpose to match: :func:`same_frame`.

``origin='lower'`` is the half that is easy to leave out, and every example here uses it:
matplotlib's default puts row zero at the *top*, row zero is the most posterior line, and an
otherwise correct image then comes back with the frontal lobe at the bottom.

Off-resonance
-------------
:func:`field_hz` hands back **the phantom's own B0 map**, which is what every example simulates
against and what :func:`slab` leaves in place unless told otherwise.  It is worth knowing exactly
what that map is, because it is easy to assume it is a measurement:

``subject05.npz`` carries ``PD_map``, ``T1_map``, ``T2_map``, ``T2dash_map`` and ``D_map`` -- and
**no** ``B0_map``.  MRzero's loader notices that and calls its own ``generate_B0_B1(PD)``, whose
comment reads *"Generate a somewhat plausible B0 and B1 map.  Visually fitted to look similar to
the numerical_brain_cropped."*  It is two Lorentzians in the distance from a point near the front
of the head, ``7/(0.05 + d^2) - 45/(0.3 + d^2)``, demeaned over the proton density.

So it is smooth, it is already centred on zero -- nothing here needs to shim it -- and it spans
about -11 to +45 Hz.  What it is not is measured, or a susceptibility calculation.  For a spin
echo or a Cartesian gradient echo that hardly matters.  For an EPI it sets the whole distortion
figure: at 128 lines and a 700 us echo spacing it displaces the object by about **-1 to +4
pixels**, which is visible and mild.  Scale it if a stronger case is wanted -- and say so if you
do, because a scaled map is no longer the phantom's.
"""

from __future__ import annotations

import urllib.request
from pathlib import Path

import numpy as np
import torch

PHANTOM_URL = ('https://github.com/MRsources/MRzero-Core/raw/main/'
               'documentation/playground_mr0/subject05.npz')
#: Beside this file rather than beside the caller, so every example finds one copy.
PHANTOM_PATH = Path(__file__).resolve().parent / 'data' / 'subject05.npz'

_VOLUMES: dict[tuple[int, int], object] = {}


def download(*, quiet: bool = False) -> Path:
    """Fetch the ~19 MB phantom into ``examples/data/`` if it is not already there."""
    if not PHANTOM_PATH.exists():
        PHANTOM_PATH.parent.mkdir(parents=True, exist_ok=True)
        if not quiet:
            print(f'downloading {PHANTOM_URL} -> {PHANTOM_PATH} (~19 MB, once)')
        urllib.request.urlretrieve(PHANTOM_URL, PHANTOM_PATH)      # noqa: S310
    return PHANTOM_PATH


def volume(*, matrix: int = 128, nz: int = 128):
    """The whole head, interpolated in plane and cached.

    BrainWeb's own in-plane resolution is 128, so ``matrix=128`` resamples nothing at all.
    """
    import MRzeroCore as mr0

    key = (int(matrix), int(nz))
    if key not in _VOLUMES:
        _VOLUMES[key] = mr0.VoxelGridPhantom.load(str(download())).interpolate(
            int(matrix), int(matrix), int(nz))
    return _VOLUMES[key]


def coil_ring(nx: int, ny: int, nz: int, *, n_coils: int = 8) -> torch.Tensor:
    """A ring of `n_coils` surface coils around the object, as complex sensitivities."""
    gx, gy = np.meshgrid(np.linspace(-1, 1, nx), np.linspace(-1, 1, ny), indexing='ij')
    sens = np.stack([
        np.exp(-((gx - 1.5 * np.cos(a)) ** 2 + (gy - 1.5 * np.sin(a)) ** 2) / 1.5)
        * np.exp(1j * np.arctan2(gy - np.sin(a), gx - np.cos(a)))
        for a in np.arange(n_coils) * 2 * np.pi / n_coils
    ]).astype(np.complex64)
    return torch.tensor(np.repeat(sens[:, :, :, None], nz, axis=3))


def slices_of(nz: int, z_center: int = 64) -> list[int]:
    """The `nz` slice indices centred on `z_center`.  One expression, so nobody writes a second."""
    first = int(z_center) - int(nz) // 2
    return list(range(first, first + int(nz)))


def slab(*, matrix: int = 128, nz: int = 4, n_coils: int = 8, z_center: int = 64,
         field_hz: np.ndarray | float | None = None):
    """A BrainWeb slab at isocentre, with a receive ring and an optional off-resonance override.

    Parameters
    ----------
    matrix, nz
        In-plane resolution and the number of slices.
    n_coils
        Elements in the receive ring.  ``0`` leaves the phantom with no coil sensitivity at all.
    z_center
        Which slice of the 128 the slab is centred on.  64 is mid-brain.
    field_hz
        Off-resonance to impose, in Hz: an array shaped like the slab, or a scalar -- ``0.0`` for a
        perfectly shimmed magnet, which is the control an off-resonance measurement needs.
        **``None`` leaves the phantom's own B0 map alone**, which is what every example that is not
        deliberately overriding it wants.  See the module docstring for what that map is.

    Returns
    -------
    VoxelGridPhantom
        Call ``.build()`` on it, or use :func:`built` which does both.
    """
    out = volume(matrix=matrix).slices(slices_of(nz, z_center))
    # slices() leaves the slab where it was in the original volume.  Put it at isocentre, which is
    # where the sequence excites -- otherwise the slice profile selects the wrong tissue, or none
    # at all, and the symptom is a dark image rather than an error.
    out.affine[2, 3] = -0.5 * nz * float(out.affine[2, 2])
    if n_coils:
        out.coil_sens = coil_ring(matrix, matrix, nz, n_coils=n_coils)
    if field_hz is not None:
        out.B0 = torch.tensor(np.broadcast_to(np.asarray(field_hz, np.float32),
                                              tuple(out.B0.shape)).copy())
    return out


def built(**kwargs):
    """:func:`slab`, built into the ``SimData`` a simulation actually takes."""
    return slab(**kwargs).build()


def field_hz(*, matrix: int = 128, nz: int = 4, z_center: int = 64) -> np.ndarray:
    """The phantom's own B0 map over the imaged slab, in Hz, as ``(matrix, matrix, nz)``.

    The same array :func:`slab` leaves in place, handed back so a notebook can draw the field it
    is simulating against beside the image it distorted.  The module docstring says what the map
    is and, more to the point, what it is not.
    """
    return np.asarray(volume(matrix=matrix).B0.numpy().real)[..., slices_of(nz, z_center)]


def mask(*, matrix: int = 128, nz: int = 4, z_center: int = 64,
         threshold: float = 0.05) -> np.ndarray:
    """Where there is tissue, as ``(matrix, matrix, nz)`` of bool."""
    return volume(matrix=matrix).PD.numpy()[..., slices_of(nz, z_center)] > threshold


def same_frame(volume_map: np.ndarray) -> np.ndarray:
    """Put a phantom-indexed ``[x, y]`` map into an image's ``[ky, kx]`` frame.

    One transpose, in one place, so a field map and the image it distorted can be drawn side by
    side without anybody working out which way round each of them goes.  A reconstruction needs
    nothing -- it is already ``[y, x]`` -- so this is only ever applied to the *map*, and both are
    drawn with ``origin='lower'``.
    """
    return np.asarray(volume_map).T
