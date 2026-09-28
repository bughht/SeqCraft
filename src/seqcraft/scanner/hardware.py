"""
Gradient hardware models, for peripheral-nerve-stimulation prediction.

A **hardware model** is not a limit.  ``Opts`` says how strong and how fast the amplifier may be
driven; this describes how the *body* responds to being driven that way -- three exponential time
constants per axis, a stimulation threshold, and the forbidden acoustic-resonance bands.  Nothing
in the compile path reads it, which is why it is not on the ``Opts``:
:func:`seqcraft.analysis.pns` takes it as a third argument.

**Nothing here bundles, discovers or implicitly locates a vendor hardware file.**  Siemens
``.asc`` gradient descriptors carry proprietary PNS/CNS response coefficients and forbidden
acoustic-resonance bands, so :func:`load_hardware` reads one only when a caller explicitly hands
it a path -- an ordinary path, which seqcraft does not police -- and :func:`synthetic_hardware`
provides a vendor-free **illustrative** stand-in, taken from pypulseq's public example model, so
PNS checks can run without any file at all.

Where such a file *lives* is repository policy rather than a runtime rule: vendor descriptors
belong outside this tree and ``.asc`` is gitignored.

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

__all__ = ['load_hardware', 'synthetic_hardware']


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

    The response parameters are taken from pypulseq's public example hardware model,
    ``safe_pns_prediction.safe_example_hw()``: every SAFE field ``Sequence.calculate_pns`` reads
    comes from the object that function returns, and seqcraft adds only its own labelling and
    provenance semantics.  **SeqCraft derives these fields from the returned upstream object
    rather than maintaining an independently maintained coefficient table.**

    .. warning::

       **This is not a scanner-specific model.**  It is an *example* response model, so that PNS
       checks can run without a vendor file, and it **must never be used to clear a sequence for
       human scanning**.  Use :func:`load_hardware` with the site's own ``.asc`` for that.

       **No upper-bound, worst-case or conservatism claim is made.**  A site model may differ
       from it in absolute peak, in which axis dominates, and in how it ranks one waveform
       against another; that is two response models disagreeing, not either one being wrong.

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
    path: str | os.PathLike[str],
    *,
    cardiac_model: bool = False,
) -> SimpleNamespace:
    """
    Load a vendor ``.asc`` gradient descriptor the caller supplies.

    Parameters
    ----------
    path
        Where the file is, as a string or :class:`~pathlib.Path`.  An ordinary path, resolved the
        way Python resolves any other: a relative one is relative to the working directory.  The
        file is **read**, never copied, and nothing about its location is retained.
    cardiac_model
        Pass ``True`` for the CNS (cardiac) response model instead of PNS.

    Returns
    -------
    SimpleNamespace
        The response model, ready for :func:`seqcraft.pns`.  It carries ``.source``, a provenance
        string of the form ``'<basename> sha256:<12 hex>'`` -- the file's *name* and hash only,
        never its directory and never its contents, both of which are vendor-confidential.

    Raises
    ------
    ConfigurationError
        If there is no file at `path`.

    Examples
    --------
    >>> hw = load_hardware('/path/to/scanner.asc')          # doctest: +SKIP
    >>> sc.pns(tree, opts, hw)                              # doctest: +SKIP

    Notes
    -----
    Vendor descriptors carry proprietary response coefficients and forbidden acoustic-resonance
    bands, so one is an **explicit opt-in**: nothing in seqcraft looks for a descriptor on its
    own, and there is no configured location for them.  A caller who has one passes it; a caller
    who has not uses :func:`synthetic_hardware`, and can tell the two apart by reading their own
    code.
    """
    where = Path(path)
    if not where.is_file():
        msg = format_error(
            f'no gradient hardware file at {where.name!r}.',
            {'path': str(path)},
            [
                'pass the path to a vendor .asc descriptor, e.g. '
                "load_hardware('/path/to/scanner.asc')",
                'or use seqcraft.hardware.synthetic_hardware() for an illustrative model that '
                'needs no file',
            ],
        )
        raise ConfigurationError(msg)

    from pypulseq.utils.siemens.asc_to_hw import asc_to_hw  # noqa: PLC0415
    from pypulseq.utils.siemens.readasc import readasc  # noqa: PLC0415

    asc, _extra = readasc(str(where))
    hardware = asc_to_hw(asc, cardiac_model=cardiac_model)
    digest = hashlib.sha256(where.read_bytes()).hexdigest()[:12]
    # The basename only: a provenance string that carried the directory would put a site's
    # filesystem layout into every record that quoted it.
    hardware.source = f'{where.name} sha256:{digest}'
    return hardware
