"""
Gradient hardware models, for peripheral-nerve-stimulation prediction.

A **hardware model** is not a limit.  ``Opts`` says how strong and how fast the amplifier may be
driven; this describes how the *body* responds to being driven that way -- three exponential time
constants per axis, a stimulation threshold, and the forbidden acoustic-resonance bands.  Nothing
in the compile path reads it, which is why it is not on the ``Opts``:
:func:`seqcraft.analysis.pns` takes it as a third argument.

**No vendor hardware file is ever read from inside this repository.**  Siemens ``.asc`` gradient
descriptors carry proprietary PNS/CNS response coefficients and forbidden acoustic-resonance
bands.  :func:`load_hardware` resolves them through the ``SEQCRAFT_ASC_DIR`` environment variable
only, and :func:`synthetic_hardware` provides a vendor-free **illustrative** stand-in, taken from
pypulseq's public example model, so PNS checks can run in CI.

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


class _SyntheticHardware(SimpleNamespace):
    """A PNS model that says what it is in its own ``repr``."""

    def __repr__(self) -> str:
        """Name the model and the caveat, so it cannot be mistaken for a measured one."""
        return (
            f'<synthetic PNS hardware {self.name!r} -- illustrative coefficients, '
            f'NOT a real scanner; never use it to clear a human scan>'
        )


def synthetic_hardware(name: str = 'synthetic_generic') -> SimpleNamespace:
    """
    Return a vendor-free PNS hardware model for tests, examples and CI.

    **Derived from pypulseq's public reference implementation**,
    ``safe_pns_prediction.safe_example_hw()``, rather than copied from it: every SAFE field
    ``Sequence.calculate_pns`` reads is taken from the upstream object, and this adds only
    SeqCraft's own metadata on top.  A second, hand-maintained table of the same coefficients had
    already drifted from the one it named -- ``y`` and ``z`` carried ``x``'s time constants and
    all three ``g_scale`` values differed -- which is exactly the failure delegating removes.
    :func:`tests.scanner.test_hardware` compares the two objects field by field so a future drift
    cannot be silent.

    .. warning::

       **This is not a real scanner.**  It is an *illustrative example* response model, so that
       PNS checks can run without a vendor file, and it **must never be used to clear a sequence
       for human scanning**.  Use :func:`load_hardware` with the site's own ``.asc`` for that.

       It is illustrative, not conservative: it is **not** an upper bound, a worst case or a
       safety margin, and it is not representative of any particular scanner.  A real model may
       differ from it in absolute peak, in which axis dominates, and in how it ranks one waveform
       against another.

       The object carries ``is_synthetic=True`` and says so in its ``repr``, so the caveat travels
       with the model rather than living in the docstring of whatever happens to consume it.

    Examples
    --------
    >>> hw = synthetic_hardware()
    >>> hw.name
    'synthetic_generic'
    >>> hw.is_synthetic
    True
    >>> hw.x.stim_limit
    30.0
    >>> hw
    <synthetic PNS hardware 'synthetic_generic' -- illustrative coefficients, NOT a real
    scanner; never use it to clear a human scan>
    """
    from pypulseq.utils.safe_pns_prediction import safe_example_hw  # noqa: PLC0415

    # Copied per axis rather than handed over: `safe_example_hw` builds a fresh object today, but
    # nothing promises that, and a caller mutating a shared upstream singleton would be a hard
    # bug to find.
    upstream = safe_example_hw()
    axes = {
        axis: SimpleNamespace(**vars(getattr(upstream, axis)))
        for axis in ('x', 'y', 'z')
    }
    return _SyntheticHardware(
        name=name,
        checkID=0,
        acoustic_resonances=(),
        is_synthetic=True,
        **axes,
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
        Bare file name, e.g. ``'scanner.asc'``.  A path containing directory separators is
        rejected: the whole point is that the location comes from the environment.
    cardiac_model
        Pass ``True`` for the CNS (cardiac) response model instead of PNS.
    directory
        Overrides ``$SEQCRAFT_ASC_DIR`` for this call.

    Returns
    -------
    SimpleNamespace
        The response model, ready for :func:`seqcraft.pns`.  It carries
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
