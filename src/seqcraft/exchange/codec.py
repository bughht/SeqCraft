"""Lossless Python codec for LogicBlock Tree Exchange Format 0.1."""

# ruff: noqa: TRY003

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
from pypulseq.opts import Opts

from ..design.logic import LogicBlock
from .errors import ExchangeError
from .schema import SCHEMA_NAME, SCHEMA_VERSION, validate_document

__all__ = [
    'export_logicblock',
    'import_logicblock',
    'read_logicblock',
    'semantic_hash',
    'write_logicblock',
]

_OPTS_TO_WIRE = {
    'max_b1': 'max_b1_hz',
    'max_grad': 'max_grad_hz_per_m',
    'max_slew': 'max_slew_hz_per_m_per_s',
    'max_freq_offset': 'max_frequency_offset_hz',
    'rf_dead_time': 'rf_dead_time_s',
    'rf_ringdown_time': 'rf_ringdown_time_s',
    'adc_dead_time': 'adc_dead_time_s',
    'adc_raster_time': 'adc_raster_time_s',
    'rf_raster_time': 'rf_raster_time_s',
    'grad_raster_time': 'grad_raster_time_s',
    'block_duration_raster': 'block_duration_raster_s',
    'adc_samples_limit': 'adc_samples_limit',
    'rf_samples_limit': 'rf_samples_limit',
    'adc_samples_divisor': 'adc_samples_divisor',
    'flag_trid': 'trigger_identifier_supported',
    'gamma': 'gamma_hz_per_t',
    'B0': 'b0_t',
}
_TIME_OPTS = {
    'rf_dead_time', 'rf_ringdown_time', 'adc_dead_time', 'adc_raster_time',
    'rf_raster_time', 'grad_raster_time', 'block_duration_raster',
}


def _number(value: Any, path: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as err:
        raise ExchangeError('expected a numeric physical value', path=path) from err
    if not math.isfinite(result):
        raise ExchangeError('numeric physical values must be finite', path=path)
    return result


def _time(value: Any, path: str = '$') -> str:
    return format(_number(value, path), '.17g')


def _time_value(value: str, path: str) -> float:
    return _number(value, path)


def _array(values: Any, path: str) -> list[float]:
    array = np.asarray(values).reshape(-1)
    return [_number(value, f'{path}[{index}]') for index, value in enumerate(array)]


def _time_array(values: Any, path: str) -> list[str]:
    array = np.asarray(values).reshape(-1)
    return [_time(value, f'{path}[{index}]') for index, value in enumerate(array)]


def _json_value(value: Any, path: str = '$') -> Any:
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return _number(value, path)
    if isinstance(value, np.ndarray):
        return [_json_value(item, f'{path}[{index}]') for index, item in enumerate(value.tolist())]
    if isinstance(value, (list, tuple)):
        return [_json_value(item, f'{path}[{index}]') for index, item in enumerate(value)]
    if isinstance(value, Mapping):
        out: dict[str, Any] = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise ExchangeError('JSON object keys must be strings', path=path)
            out[key] = _json_value(item, f'{path}.{key}')
        return out
    raise ExchangeError(f'{type(value).__name__} is not representable as protocol JSON', path=path)


def _opts_to_document(opts: Opts) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for attribute, wire_name in _OPTS_TO_WIRE.items():
        if not hasattr(opts, attribute):
            raise ExchangeError(f'Opts has no required field {attribute!r}', path='$.opts')
        value = getattr(opts, attribute)
        if attribute in _TIME_OPTS:
            out[wire_name] = _time(value, f'$.opts.{wire_name}')
        elif attribute in {'adc_samples_limit', 'rf_samples_limit', 'adc_samples_divisor'}:
            out[wire_name] = int(value)
        elif attribute == 'flag_trid':
            out[wire_name] = bool(value)
        else:
            out[wire_name] = _number(value, f'$.opts.{wire_name}')
    return out


def _opts_from_document(payload: Mapping[str, Any]) -> Opts:
    values = {attribute: payload[wire] for attribute, wire in _OPTS_TO_WIRE.items()}
    for attribute in _TIME_OPTS:
        values[attribute] = _time_value(values[attribute], f'$.opts.{_OPTS_TO_WIRE[attribute]}')
    return Opts(
        max_b1=float(values['max_b1']), b1_unit='Hz',
        max_grad=float(values['max_grad']), grad_unit='Hz/m',
        max_slew=float(values['max_slew']), slew_unit='Hz/m/s',
        max_freq_offset=float(values['max_freq_offset']),
        rf_dead_time=float(values['rf_dead_time']),
        rf_ringdown_time=float(values['rf_ringdown_time']),
        adc_dead_time=float(values['adc_dead_time']),
        adc_raster_time=float(values['adc_raster_time']),
        rf_raster_time=float(values['rf_raster_time']),
        grad_raster_time=float(values['grad_raster_time']),
        block_duration_raster=float(values['block_duration_raster']),
        adc_samples_limit=int(values['adc_samples_limit']),
        rf_samples_limit=int(values['rf_samples_limit']),
        adc_samples_divisor=int(values['adc_samples_divisor']),
        flag_trid=bool(values['flag_trid']), gamma=float(values['gamma']), B0=float(values['B0']),
    )


def _event_to_document(event: Any, path: str) -> dict[str, Any]:  # noqa: C901, PLR0915
    kind = getattr(event, 'type', '')
    if kind == 'delay':
        payload = {'duration_s': _time(event.delay, f'{path}.duration_s')}
    elif kind == 'trap':
        payload = {
            'channel': event.channel,
            'amplitude_hz_per_m': _number(event.amplitude, path),
            'rise_time_s': _time(event.rise_time, path),
            'flat_time_s': _time(event.flat_time, path),
            'fall_time_s': _time(event.fall_time, path),
            'area_per_m': _number(event.area, path),
            'flat_area_per_m': _number(event.flat_area, path),
            'delay_s': _time(event.delay, path),
            'first_hz_per_m': _number(getattr(event, 'first', 0.0), path),
            'last_hz_per_m': _number(getattr(event, 'last', 0.0), path),
        }
    elif kind == 'grad':
        payload = {
            'channel': event.channel,
            'waveform_hz_per_m': _array(event.waveform, f'{path}.waveform_hz_per_m'),
            'sample_times_s': _time_array(event.tt, f'{path}.sample_times_s'),
            'delay_s': _time(event.delay, path),
            'shape_duration_s': _time(event.shape_dur, path),
            'area_per_m': _number(event.area, path),
            'first_hz_per_m': _number(event.first, path),
            'last_hz_per_m': _number(event.last, path),
        }
    elif kind == 'rf':
        signal = np.asarray(event.signal).reshape(-1)
        payload = {
            'signal_hz': {
                'real': _array(signal.real, f'{path}.signal_hz.real'),
                'imag': _array(signal.imag, f'{path}.signal_hz.imag'),
            },
            'sample_times_s': _time_array(event.t, f'{path}.sample_times_s'),
            'shape_duration_s': _time(event.shape_dur, path),
            'frequency_offset_hz': _number(event.freq_offset, path),
            'phase_offset_rad': _number(event.phase_offset, path),
            'frequency_offset_ppm': _number(getattr(event, 'freq_ppm', 0.0), path),
            'phase_offset_ppm': _number(getattr(event, 'phase_ppm', 0.0), path),
            'dead_time_s': _time(event.dead_time, path),
            'ringdown_time_s': _time(event.ringdown_time, path),
            'delay_s': _time(event.delay, path),
            'center_s': _time(event.center, path),
            'use': getattr(event, 'use', None),
        }
    elif kind == 'adc':
        payload = {
            'num_samples': int(event.num_samples),
            'dwell_s': _time(event.dwell, path),
            'delay_s': _time(event.delay, path),
            'frequency_offset_hz': _number(event.freq_offset, path),
            'phase_offset_rad': _number(event.phase_offset, path),
            'frequency_offset_ppm': _number(getattr(event, 'freq_ppm', 0.0), path),
            'phase_offset_ppm': _number(getattr(event, 'phase_ppm', 0.0), path),
            'dead_time_s': _time(event.dead_time, path),
            'phase_modulation_rad': _array(getattr(event, 'phase_modulation', []), path),
        }
    elif kind in {'labelset', 'labelinc'}:
        payload = {'label': str(event.label), 'value': int(event.value)}
    elif kind in {'trigger', 'output'}:
        payload = {
            'channel': str(event.channel), 'delay_s': _time(event.delay, path),
            'duration_s': _time(event.duration, path),
        }
    elif kind == 'seqcraft_barrier':
        payload = {'tag': str(getattr(event, 'tag', 'barrier'))}
    else:
        raise ExchangeError(f'unsupported event type {kind!r}', path=path)
    return {'kind': 'event', 'event_type': kind, 'payload': payload}


def _block_to_document(block: LogicBlock, path: str) -> dict[str, Any]:
    nodes = []
    for index, node in enumerate(block.nodes):
        node_path = f'{path}.nodes[{index}]'
        if isinstance(node.item, LogicBlock):
            item = {'kind': 'block', 'block': _block_to_document(node.item, f'{node_path}.item.block')}
        else:
            item = _event_to_document(node.item, f'{node_path}.item.payload')
        nodes.append({'start_s': _time(node.start, f'{node_path}.start_s'), 'item': item})
    return {'tag': block.tag, 'nodes': nodes}


def _event_from_document(kind: str, p: Mapping[str, Any], path: str) -> SimpleNamespace:  # noqa: C901
    def t(name: str) -> float:
        return _time_value(p[name], f'{path}.{name}')

    if kind == 'delay':
        return SimpleNamespace(type=kind, delay=t('duration_s'))
    if kind == 'trap':
        return SimpleNamespace(
            type=kind, channel=p['channel'], amplitude=float(p['amplitude_hz_per_m']),
            rise_time=t('rise_time_s'), flat_time=t('flat_time_s'), fall_time=t('fall_time_s'),
            area=float(p['area_per_m']), flat_area=float(p['flat_area_per_m']), delay=t('delay_s'),
            first=float(p['first_hz_per_m']), last=float(p['last_hz_per_m']),
        )
    if kind == 'grad':
        return SimpleNamespace(
            type=kind, channel=p['channel'], waveform=np.asarray(p['waveform_hz_per_m'], dtype=float),
            tt=np.asarray([_time_value(x, path) for x in p['sample_times_s']], dtype=float),
            delay=t('delay_s'), shape_dur=t('shape_duration_s'), area=float(p['area_per_m']),
            first=float(p['first_hz_per_m']), last=float(p['last_hz_per_m']),
        )
    if kind == 'rf':
        real = np.asarray(p['signal_hz']['real'], dtype=float)
        imag = np.asarray(p['signal_hz']['imag'], dtype=float)
        return SimpleNamespace(
            type=kind, signal=real + 1j * imag,
            t=np.asarray([_time_value(x, path) for x in p['sample_times_s']], dtype=float),
            shape_dur=t('shape_duration_s'), freq_offset=float(p['frequency_offset_hz']),
            phase_offset=float(p['phase_offset_rad']), freq_ppm=float(p['frequency_offset_ppm']),
            phase_ppm=float(p['phase_offset_ppm']), dead_time=t('dead_time_s'),
            ringdown_time=t('ringdown_time_s'), delay=t('delay_s'), center=t('center_s'), use=p['use'],
        )
    if kind == 'adc':
        return SimpleNamespace(
            type=kind, num_samples=int(p['num_samples']), dwell=t('dwell_s'), delay=t('delay_s'),
            freq_offset=float(p['frequency_offset_hz']), phase_offset=float(p['phase_offset_rad']),
            freq_ppm=float(p['frequency_offset_ppm']), phase_ppm=float(p['phase_offset_ppm']),
            dead_time=t('dead_time_s'), phase_modulation=list(p['phase_modulation_rad']),
        )
    if kind in {'labelset', 'labelinc'}:
        return SimpleNamespace(type=kind, label=p['label'], value=int(p['value']))
    if kind in {'trigger', 'output'}:
        return SimpleNamespace(type=kind, channel=p['channel'], delay=t('delay_s'),
                               duration=t('duration_s'))
    if kind == 'seqcraft_barrier':
        return SimpleNamespace(type=kind, tag=p['tag'], delay=0.0)
    raise ExchangeError(f'unsupported event type {kind!r}', path=path)


def _block_from_document(payload: Mapping[str, Any], path: str) -> LogicBlock:
    block = LogicBlock(str(payload['tag']))
    for index, node in enumerate(payload['nodes']):
        node_path = f'{path}.nodes[{index}]'
        item = node['item']
        if item['kind'] == 'block':
            value = _block_from_document(item['block'], f'{node_path}.item.block')
        else:
            value = _event_from_document(item['event_type'], item['payload'],
                                         f'{node_path}.item.payload')
        block.add(_time_value(node['start_s'], f'{node_path}.start_s'), value)
    return block


def export_logicblock(
    root: LogicBlock,
    opts: Opts,
    *,
    definitions: Mapping[str, Any] | None = None,
    provenance: Mapping[str, Any] | None = None,
    extensions: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return a validated LBTX document for one tree and scanner configuration."""
    if not isinstance(root, LogicBlock):
        raise ExchangeError(f'expected LogicBlock, got {type(root).__name__}', path='$.root')
    document = {
        'schema': SCHEMA_NAME,
        'version': SCHEMA_VERSION,
        'units': 'SI',
        'opts': _opts_to_document(opts),
        'definitions': _json_value(definitions or {}, '$.definitions'),
        'root': _block_to_document(root, '$.root'),
    }
    if provenance is not None:
        document['provenance'] = _json_value(provenance, '$.provenance')
    if extensions is not None:
        document['extensions'] = _json_value(extensions, '$.extensions')
    validate_document(document)
    return document


def import_logicblock(document: Mapping[str, Any]) -> tuple[LogicBlock, Opts, dict[str, Any]]:
    """Validate and reconstruct ``(root, opts, definitions)`` from an LBTX document."""
    validate_document(document)
    root = _block_from_document(document['root'], '$.root')
    opts = _opts_from_document(document['opts'])
    definitions = _json_value(document['definitions'], '$.definitions')
    return root, opts, definitions


def write_logicblock(
    path: str | Path,
    root: LogicBlock,
    opts: Opts,
    *,
    definitions: Mapping[str, Any] | None = None,
    provenance: Mapping[str, Any] | None = None,
    extensions: Mapping[str, Any] | None = None,
) -> Path:
    """Write one validated, reviewable LBTX JSON document and return its path."""
    target = Path(path)
    document = export_logicblock(root, opts, definitions=definitions, provenance=provenance,
                                 extensions=extensions)
    target.write_text(json.dumps(document, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    return target


def read_logicblock(path: str | Path) -> tuple[LogicBlock, Opts, dict[str, Any]]:
    """Read, validate, and reconstruct one LBTX JSON file."""
    source = Path(path)
    try:
        document = json.loads(source.read_text(encoding='utf-8'))
    except json.JSONDecodeError as err:
        raise ExchangeError(f'invalid JSON at line {err.lineno}, column {err.colno}: {err.msg}') from err
    return import_logicblock(document)


def semantic_hash(document: Mapping[str, Any]) -> str:
    """Return the SHA-256 of canonical semantic content, excluding provenance and extensions."""
    root, opts, definitions = import_logicblock(document)
    canonical = export_logicblock(root, opts, definitions=definitions)
    encoded = json.dumps(canonical, sort_keys=True, separators=(',', ':'), ensure_ascii=False)
    return hashlib.sha256(encoded.encode('utf-8')).hexdigest()
