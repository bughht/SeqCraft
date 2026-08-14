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
"""

from __future__ import annotations

import importlib
from typing import Any

import pypulseq as pp

from .core.errors import ConfigurationError, format_error

__all__ = [
    'PYPULSEQ_VERSION',
    'has',
    'probe',
    'require',
    'rotate_3d',
    'supported_rf_uses',
]

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

#: Optional capabilities.  Absence disables a feature rather than failing the import.
_OPTIONAL = {
    'soft_delay': 'make_soft_delay',
    'rotation_extension': 'make_rotation',
    'adc_segments': 'calc_adc_segments',
    'b_tensor': None,          # Sequence.calc_moments_b_tensor, checked as a method
    'gradient_spectrum': None,  # Sequence.calculate_gradient_spectrum
}


def probe() -> dict[str, bool]:
    """
    Report which optional pypulseq capabilities are available.

    Returns
    -------
    dict
        ``capability -> present``.

    Examples
    --------
    >>> caps = probe()
    >>> caps['soft_delay'] in (True, False)
    True
    """
    out: dict[str, bool] = {}
    for name, attr in _OPTIONAL.items():
        if attr is not None:
            out[name] = hasattr(pp, attr)
    out['b_tensor'] = hasattr(pp.Sequence, 'calc_moments_b_tensor')
    out['gradient_spectrum'] = hasattr(pp.Sequence, 'calculate_gradient_spectrum')
    out['rotate_3d'] = _rotate_3d_available()
    return out


def has(capability: str) -> bool:
    """Report whether one optional capability is available."""
    return probe().get(capability, False)


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


def _rotate_3d_available() -> bool:
    """``pypulseq.rotate_3d`` is not re-exported at package level; import the submodule."""
    try:
        importlib.import_module('pypulseq.rotate_3d')
    except ImportError:
        return False
    return True


def rotate_3d(*args: Any, **kwargs: Any) -> Any:
    """
    Call ``pypulseq.rotate_3d.rotate_3d``, imported here and nowhere else.

    Kept behind one accessor because the function is **not** exported from the pypulseq
    package namespace, so the awkward import path is confined to a single place.
    """
    module = importlib.import_module('pypulseq.rotate_3d')
    return module.rotate_3d(*args, **kwargs)


def supported_rf_uses() -> tuple[str, ...]:
    """
    Return the RF ``use`` strings pypulseq accepts.

    Lives in ``pypulseq.supported_labels_rf_use`` rather than at package level, so the
    import path is confined here alongside the other non-exported helpers.
    """
    module = importlib.import_module('pypulseq.supported_labels_rf_use')
    return tuple(module.get_supported_rf_uses())
