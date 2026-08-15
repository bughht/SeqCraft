"""
Gradient hardware models, for peripheral-nerve-stimulation prediction.

A **hardware model** is not a limit.  ``Opts`` says how strong and how fast the amplifier may be
driven; this describes how the *body* responds to being driven that way -- three exponential time
constants per axis, a stimulation threshold, and the forbidden acoustic-resonance bands.  Nothing
in the compile path reads it, which is why it is not on the ``Opts`` and not in ``core``:
:meth:`seqcraft.CompiledSequence.pns` takes it as an argument, on a sequence that has already been
compiled.

**No vendor hardware file is ever read from inside this repository.**  Siemens ``.asc`` gradient
descriptors carry proprietary PNS/CNS response coefficients and forbidden acoustic-resonance
bands.  :func:`load_hardware` resolves them through the ``SEQCRAFT_ASC_DIR`` environment variable
only, and :func:`synthetic_hardware` provides a vendor-free stand-in so PNS checks can run in CI.

Examples
--------
>>> import seqcraft as sc
>>> hw = sc.hardware.synthetic_hardware()
>>> hw.is_synthetic
True
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
from types import SimpleNamespace

from ..errors import ConfigurationError, format_error

__all__ = ['ASC_ENV_VAR', 'load_hardware', 'synthetic_hardware']

#: Environment variable naming a directory of vendor ``.asc`` files.  Never a path inside
#: this repository; ``.asc`` is gitignored.
ASC_ENV_VAR = 'SEQCRAFT_ASC_DIR'


def synthetic_hardware(name: str = 'synthetic_generic') -> SimpleNamespace:
    """
    Return a vendor-free PNS hardware model for tests and CI.

    Shaped like the output of ``pypulseq.utils.siemens.asc_to_hw`` so it can be handed
    straight to ``Sequence.calculate_pns``, but the coefficients are the illustrative
    values from pypulseq's own ``safe_pns_prediction.safe_example_hw()`` reference
    implementation, not measurements from any scanner.  **It is not a real scanner and
    must never be used to clear a sequence for human scanning** -- use
    :func:`load_hardware` with the site's own ``.asc`` for that.

    Examples
    --------
    >>> hw = synthetic_hardware()
    >>> hw.name
    'synthetic_generic'
    >>> hw.x.stim_limit
    30.0
    """

    def axis(stim_limit: float, stim_thresh: float, g_scale: float) -> SimpleNamespace:
        return SimpleNamespace(
            tau1=0.20,
            tau2=0.03,
            tau3=3.0,
            a1=0.4,
            a2=0.10,
            a3=0.50,
            stim_limit=stim_limit,
            stim_thresh=stim_thresh,
            g_scale=g_scale,
        )

    return SimpleNamespace(
        name=name,
        checkID=0,
        x=axis(30.0, 24.0, 0.4),
        y=axis(15.0, 12.0, 0.7),
        z=axis(25.0, 20.0, 0.3),
        acoustic_resonances=(),
        is_synthetic=True,
    )


def load_hardware(
    filename: str,
    *,
    cardiac_model: bool = False,
    directory: str | os.PathLike[str] | None = None,
) -> SimpleNamespace:
    """
    Load a vendor ``.asc`` gradient descriptor from outside the repository.

    Parameters
    ----------
    filename
        Bare file name, e.g. ``'CimaX.asc'``.  A path containing directory separators is
        rejected: the whole point is that the location comes from the environment.
    cardiac_model
        Pass ``True`` for the CNS (cardiac) response model instead of PNS.
    directory
        Overrides ``$SEQCRAFT_ASC_DIR`` for this call.

    Returns
    -------
    SimpleNamespace
        The response model, ready for :meth:`seqcraft.CompiledSequence.pns`.  It carries
        ``.source``, a provenance string of the form ``'<filename> sha256:<12 hex>'`` -- the
        file *name* and hash only, never the contents, which are vendor-confidential.

    Raises
    ------
    ConfigurationError
        If ``SEQCRAFT_ASC_DIR`` is unset, the file is missing, or `filename` is a path.

    Notes
    -----
    ``.asc`` files describe the gradient amplifier's PNS/CNS response and its forbidden
    acoustic-resonance frequencies.  They are site- and vendor-confidential, are excluded
    by ``.gitignore``, and are never bundled, copied into, or read from this repository.

    Only the response model is returned.  The acoustic-resonance bands are in the same file and
    nothing in seqcraft checks against them, so returning them would be inventing a consumer; read
    them with pypulseq's own ``asc_to_acoustic_resonances`` if you need them.
    """
    if os.sep in filename or '/' in filename:
        msg = format_error(
            f'load_hardware() takes a bare file name, got a path: {filename!r}.',
            {'reason': 'the directory must come from $' + ASC_ENV_VAR},
        )
        raise ConfigurationError(msg)

    root = directory or os.environ.get(ASC_ENV_VAR)
    if not root:
        msg = format_error(
            f'no gradient hardware directory configured, cannot load {filename!r}.',
            {
                'expected': f'${ASC_ENV_VAR} pointing at a directory of vendor .asc files',
                'why': 'vendor .asc files carry proprietary PNS and acoustic-resonance data '
                'and are never stored in this repository',
            },
            [
                f'set {ASC_ENV_VAR} to the directory holding your .asc files',
                'or use seqcraft.hardware.synthetic_hardware() for tests',
            ],
        )
        raise ConfigurationError(msg)

    path = Path(root) / filename
    if not path.is_file():
        msg = format_error(
            f'gradient hardware file not found: {filename!r}.',
            {'looked in': str(root)},
        )
        raise ConfigurationError(msg)

    from pypulseq.utils.siemens.asc_to_hw import asc_to_hw  # noqa: PLC0415
    from pypulseq.utils.siemens.readasc import readasc  # noqa: PLC0415

    asc, _extra = readasc(str(path))
    hardware = asc_to_hw(asc, cardiac_model=cardiac_model)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()[:12]
    hardware.source = f'{filename} sha256:{digest}'
    return hardware
