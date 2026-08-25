"""The Python side reads MATLAB Pulseq structs without redefining their event types."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import seqcraft as sc
from seqcraft.exchange import ExchangeError, read_logicblock, validate_document

FIXTURE = Path(__file__).parents[1] / 'fixtures' / 'lbtx' / 'minimal.lbtx.json'


def test_handwritten_fixture_reads_and_compiles() -> None:
    root, opts, definitions = read_logicblock(FIXTURE)
    sequence = sc.compile(root, opts, definitions=definitions)

    assert root.tag == 'minimal'
    assert [node.start for node in root] == [0.0, 0.001, 0.002]
    assert root.nodes[1].item.tag == 'readout'
    assert root.nodes[1].item.nodes[0].item.type == 'trap'
    assert root.nodes[1].item.nodes[1].item.num_samples == 64
    assert opts.max_grad == 1703040
    assert opts.B0 == 3
    assert sequence.duration()[0] == pytest.approx(0.002)
    assert sequence.definitions['FOV'] == [0.22, 0.22, 0.005]


def test_generic_complex_and_array_values_decode(tmp_path) -> None:
    document = json.loads(FIXTURE.read_text(encoding='utf-8'))
    event = {
        'type': 'rf',
        'signal': {'real': [1.0, 3.0], 'imag': [2.0, 4.0]},
        't': [0.5e-6, 1.5e-6],
        'shape_dur': 2e-6,
        'freq_offset': 0,
        'phase_offset': 0,
        'freq_ppm': 0,
        'phase_ppm': 0,
        'dead_time': 0,
        'ringdown_time': 0,
        'delay': 0,
        'center': 1e-6,
        'use': 'excitation',
    }
    document['root']['nodes'] = [{'start_s': 0, 'event': event}]
    path = tmp_path / 'rf.lbtx.json'
    path.write_text(json.dumps(document), encoding='utf-8')

    root, _, _ = read_logicblock(path)
    restored = root.nodes[0].item

    assert np.array_equal(restored.signal, np.array([1 + 2j, 3 + 4j]))
    assert np.array_equal(restored.t, np.array([0.5e-6, 1.5e-6]))


@pytest.mark.parametrize(
    ('mutation', 'path'),
    [
        (lambda doc: doc.pop('version'), '$'),
        (lambda doc: doc['root'].update({'unknown': True}), '$.root'),
        (lambda doc: doc['root']['nodes'][0].update({'start_s': float('nan')}),
         '$.root.nodes[0].start_s'),
    ],
)
def test_malformed_documents_name_a_json_path(mutation, path) -> None:
    document = json.loads(FIXTURE.read_text(encoding='utf-8'))
    mutation(document)

    with pytest.raises(ExchangeError) as caught:
        validate_document(document)

    assert caught.value.path == path


def test_unknown_event_uses_existing_seqcraft_vocabulary(tmp_path) -> None:
    document = json.loads(FIXTURE.read_text(encoding='utf-8'))
    document['root']['nodes'] = [{'start_s': 0, 'event': {'type': 'not-pulseq'}}]
    path = tmp_path / 'unknown.lbtx.json'
    path.write_text(json.dumps(document), encoding='utf-8')

    with pytest.raises(ExchangeError, match='unsupported Pulseq event type'):
        read_logicblock(path)
