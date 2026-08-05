"""
Module registry.

Every concrete :class:`~seqcraft.core.module.SeqModule` subclass registers itself, and the
contract test suite parametrises over the registry.  That is what makes "add a subclass
without touching the base" honest rather than aspirational: a new module inherits the full
purity / anchor / raster / limits / docstring suite the moment it is registered, with no
new test file.

Examples
--------
>>> from seqcraft.core.registry import registered, lookup
>>> 'phase_encode' in registered()
True
>>> lookup('phase_encode').__name__
'PhaseEncode'
"""

from __future__ import annotations

import difflib
from typing import TYPE_CHECKING, TypeVar

from .errors import ConfigurationError, format_error

if TYPE_CHECKING:
    from collections.abc import Mapping

    from .module import SeqModule

__all__ = ['lookup', 'register', 'registered']

_REGISTRY: dict[str, type] = {}

T = TypeVar('T', bound=type)


def register(name: str | None = None) -> object:
    """
    Class decorator registering a module class under a snake_case name.

    Parameters
    ----------
    name
        Registry key.  Defaults to the class name converted to snake_case
        (``PhaseEncode`` -> ``phase_encode``).

    Returns
    -------
    callable
        The decorator.

    Raises
    ------
    ConfigurationError
        If the name is already taken, which would make the contract suite silently skip
        one of the two classes.
    """

    def decorate(cls: T) -> T:
        key = name or _snake_case(cls.__name__)
        if key in _REGISTRY and _REGISTRY[key] is not cls:
            msg = format_error(
                f'registry name {key!r} is already taken.',
                {'existing': _REGISTRY[key].__qualname__, 'new': cls.__qualname__},
                ["pass an explicit name: @register('my_name')"],
            )
            raise ConfigurationError(msg)
        cls.registry_name = key  # type: ignore[attr-defined]
        _REGISTRY[key] = cls
        return cls

    return decorate


def registered() -> Mapping[str, type[SeqModule]]:
    """Return the full registry as an immutable view, keyed by snake_case name."""
    from types import MappingProxyType  # noqa: PLC0415

    return MappingProxyType(_REGISTRY)  # type: ignore[arg-type]


def lookup(name: str) -> type[SeqModule]:
    """
    Return the module class registered under `name`.

    Raises
    ------
    ConfigurationError
        If unknown, naming the closest matches.
    """
    try:
        return _REGISTRY[name]  # type: ignore[return-value]
    except KeyError:
        fields: dict[str, object] = {'known': ', '.join(sorted(_REGISTRY))}
        close = difflib.get_close_matches(name, list(_REGISTRY), n=3)
        if close:
            fields['did you mean'] = ', '.join(close)
        raise ConfigurationError(format_error(f'no module registered as {name!r}.', fields)) from None


def _snake_case(name: str) -> str:
    """
    Convert CamelCase to snake_case, keeping acronyms together.

    Examples
    --------
    >>> _snake_case('PhaseEncode')
    'phase_encode'
    >>> _snake_case('EPIReadout')
    'epi_readout'
    >>> _snake_case('SincExcitation')
    'sinc_excitation'
    """
    out: list[str] = []
    for i, ch in enumerate(name):
        if ch.isupper() and i:
            prev_lower = name[i - 1].islower()
            next_lower = i + 1 < len(name) and name[i + 1].islower()
            if prev_lower or next_lower:
                out.append('_')
        out.append(ch.lower())
    return ''.join(out)
