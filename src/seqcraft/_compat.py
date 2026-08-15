"""
pypulseq feature probing.

The installed pypulseq is checked by **capability, never by version string**.  There is a
concrete reason: the build in use here reports ``__version__ == '1.5.1'`` while its
distribution metadata says ``1.5.0.post1``, so a version comparison would be wrong in one
direction or the other.  Probing what the module actually exposes is both simpler and
correct.

:func:`require` is called once from :mod:`seqcraft` at import time and raises a single clear
error listing everything missing, rather than letting the first module that needs a function
fail with an opaque ``AttributeError`` deep in a build.

That is the whole of this module's job.  It carries no optional-capability registry: a feature
seqcraft does not use needs no probe, and one it does use belongs in the required list, where a
missing entry is an error at import rather than a silently disabled branch.
"""

from __future__ import annotations

import pypulseq as pp

from .errors import ConfigurationError, format_error

__all__ = ['PYPULSEQ_VERSION', 'require']

#: What ``pypulseq.__version__`` reports.  Informational only -- never compared.
PYPULSEQ_VERSION: str = getattr(pp, '__version__', 'unknown')

#: Functions seqcraft calls that must exist on the top-level pypulseq module.
_REQUIRED_TOPLEVEL = (
    'Sequence',
    'Opts',
    'calc_duration',
    'calc_rf_bandwidth',
    'calc_rf_center',
    'make_adc',
    'make_delay',
    'make_label',
    'make_sinc_pulse',
    'make_trapezoid',
    'make_trigger',
    'scale_grad',
)

# These fields identify the Pulseq 1.5.1 compatibility build more reliably than its package
# version, which is identical to the incompatible PyPI wheel.
_REQUIRED_OPTS_FIELDS = (
    'max_b1',
    'max_freq_offset',
    'rf_samples_limit',
)

_PYPULSEQ_SOURCE = (
    'https://github.com/m-a-x-i-m-z/pypulseq-matlab-like/archive/'
    '22ef2db1f71ff38c8ce355c61913cfd8fceaac3b.zip'
)

def require() -> None:
    """
    Check that every function seqcraft depends on is present.

    Raises
    ------
    ConfigurationError
        Listing everything missing, with the install hint.
    """
    missing = [name for name in _REQUIRED_TOPLEVEL if not hasattr(pp, name)]
    if hasattr(pp, 'Opts'):
        opts = pp.Opts()
        missing.extend(f'Opts.{name}' for name in _REQUIRED_OPTS_FIELDS if not hasattr(opts, name))
    if missing:
        msg = format_error(
            'the installed pypulseq is missing capabilities seqcraft needs.',
            {
                'missing': ', '.join(missing),
                'found version': PYPULSEQ_VERSION,
                'module': getattr(pp, '__file__', 'unknown'),
            },
            [
                f'python -m pip install --force-reinstall "pypulseq @ {_PYPULSEQ_SOURCE}"',
                'seqcraft probes capabilities rather than version strings, so a partial or '
                'patched install shows up here',
            ],
        )
        raise ConfigurationError(msg)
