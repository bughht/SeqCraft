"""Read MATLAB-produced LBTX into the existing Python compiler model."""

# ruff: noqa: TRY003

from __future__ import annotations

import inspect
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
from pypulseq import Opts

from ..design.events import HANDLED_KINDS
from ..design.logic import BARRIER, LogicBlock, barrier
from .errors import ExchangeError
from .schema import validate_document

__all__ = ['read_logicblock']


def _event_value(value: Any, path: str) -> Any:
    if isinstance(value, dict):
        if set(value) == {'real', 'imag'}:
            real = np.asarray(value['real'], dtype=float)
            imag = np.asarray(value['imag'], dtype=float)
            if real.shape != imag.shape:
                raise ExchangeError('complex real and imaginary values must have equal shapes',
                                    path=path)
            return real + 1j * imag
        return {key: _event_value(item, f'{path}.{key}') for key, item in value.items()}
    if isinstance(value, list):
        decoded = [_event_value(item, f'{path}[{index}]') for index, item in enumerate(value)]
        if all(isinstance(item, (bool, int, float, complex, np.number)) for item in decoded):
            return np.asarray(decoded)
        return decoded
    return value


def _event_from_document(payload: dict[str, Any], path: str) -> SimpleNamespace:
    kind = str(payload['type'])
    if kind not in HANDLED_KINDS or kind == BARRIER:
        raise ExchangeError(f'unsupported Pulseq event type {kind!r}', path=f'{path}.type')
    return SimpleNamespace(**{
        key: _event_value(value, f'{path}.{key}') for key, value in payload.items()
    })


def _block_from_document(payload: dict[str, Any], path: str) -> LogicBlock:
    block = LogicBlock(str(payload['tag']))
    for index, node in enumerate(payload['nodes']):
        node_path = f'{path}.nodes[{index}]'
        if 'block' in node:
            item = _block_from_document(node['block'], f'{node_path}.block')
        elif 'barrier' in node:
            item = barrier(str(node['barrier']))
        else:
            item = _event_from_document(node['event'], f'{node_path}.event')
        block.add(float(node['start_s']), item)
    return block


def _opts_from_document(payload: dict[str, Any]) -> Opts:
    parameters = inspect.signature(Opts).parameters
    by_casefold = {name.casefold(): name for name in parameters}
    kwargs: dict[str, Any] = {}
    for field, value in payload.items():
        parameter = by_casefold.get(field.casefold())
        if parameter is None or parameter in {'set_as_default', 'reset_default'}:
            raise ExchangeError(f'unsupported mr.opts field {field!r}', path=f'$.opts.{field}')
        kwargs[parameter] = None if isinstance(value, list) and not value else value

    # MATLAB Pulseq stores these three quantities in its canonical internal units.
    kwargs.update(grad_unit='Hz/m', slew_unit='Hz/m/s', b1_unit='Hz')
    return Opts(**kwargs)


def read_logicblock(path: str | Path) -> tuple[LogicBlock, Opts, dict[str, Any]]:
    """Read and validate one MATLAB-produced LBTX document."""
    source = Path(path)
    try:
        document = json.loads(source.read_text(encoding='utf-8'))
    except json.JSONDecodeError as err:
        message = f'invalid JSON at line {err.lineno}, column {err.colno}: {err.msg}'
        raise ExchangeError(message) from err
    validate_document(document)
    root = _block_from_document(document['root'], '$.root')
    opts = _opts_from_document(document['opts'])
    return root, opts, document['definitions']
