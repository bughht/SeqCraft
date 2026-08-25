"""Errors raised at the LogicBlock Tree Exchange boundary."""

from __future__ import annotations

from ..errors import ConfigurationError

__all__ = ['ExchangeError']


class ExchangeError(ConfigurationError):
    """An LBTX document is malformed, unsupported, or cannot preserve its input semantics."""

    def __init__(self, message: str, *, path: str = '$') -> None:
        self.path = path
        super().__init__(f'{path}: {message}')
