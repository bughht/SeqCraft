"""Schema loading and validation for LBTX 0.1."""

# ruff: noqa: TRY003

from __future__ import annotations

import copy
import json
import math
from decimal import Decimal, InvalidOperation
from importlib.resources import files
from typing import Any

from jsonschema import Draft202012Validator
from jsonschema.exceptions import best_match

from .errors import ExchangeError

__all__ = ['SCHEMA_NAME', 'SCHEMA_VERSION', 'exchange_schema', 'validate_document']

SCHEMA_NAME = 'seqcraft.logicblock'
SCHEMA_VERSION = '0.1'


def exchange_schema() -> dict[str, Any]:
    """Return an independent copy of the LBTX 0.1 JSON Schema."""
    resource = files(__package__).joinpath('lbtx-0.1.schema.json')
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


def _decimal(value: str, path: str) -> Decimal:
    try:
        parsed = Decimal(value)
    except (InvalidOperation, ValueError) as err:
        raise ExchangeError('expected a finite decimal time string', path=path) from err
    if not parsed.is_finite():
        raise ExchangeError('time values must be finite', path=path)
    return parsed


def _validate_event(item: dict[str, Any], path: str) -> None:
    event_type = item['event_type']
    payload = item['payload']
    for name, value in payload.items():
        if name.endswith('_s') and isinstance(value, str):
            if _decimal(value, f'{path}.payload.{name}') < 0:
                raise ExchangeError('event time values must be non-negative',
                                    path=f'{path}.payload.{name}')
        elif name.endswith('_times_s'):
            previous: Decimal | None = None
            for index, text in enumerate(value):
                current = _decimal(text, f'{path}.payload.{name}[{index}]')
                if current < 0:
                    raise ExchangeError('sample times must be non-negative',
                                        path=f'{path}.payload.{name}[{index}]')
                if previous is not None and current <= previous:
                    raise ExchangeError('sample times must be strictly increasing',
                                        path=f'{path}.payload.{name}[{index}]')
                previous = current

    if event_type == 'grad':
        if len(payload['waveform_hz_per_m']) != len(payload['sample_times_s']):
            raise ExchangeError('waveform and sample_times must have equal lengths',
                                path=f'{path}.payload')
        if not payload['waveform_hz_per_m']:
            raise ExchangeError('gradient waveform must not be empty',
                                path=f'{path}.payload.waveform_hz_per_m')
    elif event_type == 'rf':
        signal = payload['signal_hz']
        sizes = (len(signal['real']), len(signal['imag']), len(payload['sample_times_s']))
        if len(set(sizes)) != 1 or sizes[0] == 0:
            raise ExchangeError('RF real, imaginary, and sample-time arrays must have equal '
                                'non-zero lengths', path=f'{path}.payload.signal_hz')
    elif event_type == 'adc':
        modulation = payload['phase_modulation_rad']
        if modulation and len(modulation) != payload['num_samples']:
            raise ExchangeError('ADC phase modulation must be empty or have num_samples entries',
                                path=f'{path}.payload.phase_modulation_rad')


def _validate_block(block: dict[str, Any], path: str) -> None:
    for index, node in enumerate(block['nodes']):
        node_path = f'{path}.nodes[{index}]'
        _decimal(node['start_s'], f'{node_path}.start_s')
        item = node['item']
        if item['kind'] == 'block':
            _validate_block(item['block'], f'{node_path}.item.block')
        else:
            _validate_event(item, f'{node_path}.item')


def validate_document(document: Any) -> None:
    """Validate one LBTX document, raising :class:`ExchangeError` with a JSON path."""
    errors = list(_VALIDATOR.iter_errors(document))
    if errors:
        error = best_match(errors)
        assert error is not None
        raise ExchangeError(error.message, path=_json_path(error.absolute_path))
    _require_finite(document)
    _validate_block(document['root'], '$.root')
