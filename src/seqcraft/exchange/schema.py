"""Validation for the small structure owned by LBTX."""

# ruff: noqa: TRY003

from __future__ import annotations

import copy
import json
import math
from importlib.resources import files
from typing import Any

from jsonschema import Draft202012Validator
from jsonschema.exceptions import best_match

from .errors import ExchangeError

__all__ = ['FORMAT_NAME', 'FORMAT_VERSION', 'exchange_schema', 'validate_document']

FORMAT_NAME = 'seqcraft.logicblock'
FORMAT_VERSION = 1


def exchange_schema() -> dict[str, Any]:
    """Return an independent copy of the LBTX 1 JSON Schema."""
    resource = files(__package__).joinpath('lbtx-1.schema.json')
    return copy.deepcopy(json.loads(resource.read_text(encoding='utf-8')))


_VALIDATOR = Draft202012Validator(exchange_schema())


def _json_path(parts: Any) -> str:
    path = '$'
    for part in parts:
        path += f'[{part}]' if isinstance(part, int) else f'.{part}'
    return path


def _require_finite(value: Any, path: str = '$') -> None:
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return
    if isinstance(value, (int, float)):
        if not math.isfinite(float(value)):
            raise ExchangeError('numeric values must be finite', path=path)
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            _require_finite(item, f'{path}[{index}]')
        return
    if isinstance(value, dict):
        for key, item in value.items():
            _require_finite(item, f'{path}.{key}')


def validate_document(document: Any) -> None:
    """Validate only the tree structure LBTX owns; Pulseq validates event semantics."""
    errors = list(_VALIDATOR.iter_errors(document))
    if errors:
        error = best_match(errors)
        assert error is not None
        raise ExchangeError(error.message, path=_json_path(error.absolute_path))
    _require_finite(document)
