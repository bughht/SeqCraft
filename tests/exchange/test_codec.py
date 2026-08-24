"""LBTX preserves the tree semantics the compiler consumes."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc

FIXTURES = Path(__file__).parents[1] / 'fixtures' / 'lbtx'


def _event_tree(opts) -> sc.LogicBlock:
    rf = pp.make_block_pulse(flip_angle=0.5, duration=100e-6, delay=opts.rf_dead_time,
                             system=opts, use='excitation')
    grad = pp.make_arbitrary_grad('y', [0.0, 100.0, 0.0], first=0.0, last=0.0, system=opts)
    adc = pp.make_adc(num_samples=8, dwell=1e-6, delay=opts.adc_dead_time, system=opts)
    return sc.LogicBlock('all-events').add([
        [0.0, pp.make_delay(0.5e-3)],
        [0.5e-3, pp.make_trapezoid('x', area=20.0, duration=0.5e-3, system=opts)],
        [1.0e-3, grad],
        [1.1e-3, rf],
        [1.4e-3, pp.make_label('LIN', 'SET', 2), pp.make_label('REP', 'INC', 1), adc],
        [1.6e-3, pp.make_trigger('physio1', duration=20e-6, delay=10e-6)],
        [1.7e-3, pp.make_digital_output_pulse('osc0', duration=20e-6, delay=10e-6)],
        [1.8e-3, sc.barrier('cut')],
        [1.9e-3, sc.LogicBlock('nested').add(0.0, pp.make_delay(0.1e-3))],
    ])


def _semantic_document(document):
    return {key: value for key, value in document.items() if key not in {'provenance', 'extensions'}}


def test_handwritten_minimum_fixture_validates_and_compiles() -> None:
    document = json.loads((FIXTURES / 'minimal.lb.json').read_text(encoding='utf-8'))

    sc.validate_document(document)
    root, opts, definitions = sc.import_logicblock(document)
    sequence = sc.compile(root, opts, definitions=definitions)

    assert root.tag == 'minimal'
    assert [node.start for node in root] == [0.0, 0.001]
    assert sequence.duration()[0] == pytest.approx(0.002)
    assert sequence.definitions['FOV'] == [0.22, 0.22, 0.005]


def test_every_supported_event_round_trips_in_insertion_order(opts) -> None:
    root = _event_tree(opts)
    document = sc.export_logicblock(
        root, opts, definitions={'FOV': np.array([0.22, 0.22, 0.005])},
        provenance={'producer': 'pytest'}, extensions={'test:note': {'ignored': True}},
    )

    restored, restored_opts, definitions = sc.import_logicblock(document)
    canonical = sc.export_logicblock(restored, restored_opts, definitions=definitions)

    assert canonical == _semantic_document(document)
    assert [node.start for node in restored] == [node.start for node in root]
    assert [node.item.type if hasattr(node.item, 'type') else 'block' for node in restored] == [
        'delay', 'trap', 'grad', 'rf', 'labelset', 'labelinc', 'adc', 'trigger', 'output',
        'seqcraft_barrier', 'block',
    ]
    assert restored.nodes[-1].item.tag == 'nested'
    assert restored_opts.max_grad == opts.max_grad
    assert restored_opts.max_slew == opts.max_slew
    assert restored_opts.grad_raster_time == opts.grad_raster_time
    assert definitions['FOV'] == [0.22, 0.22, 0.005]


def test_round_tripped_all_event_tree_compiles_the_same(opts) -> None:
    root = _event_tree(opts)
    restored, restored_opts, _ = sc.import_logicblock(sc.export_logicblock(root, opts))

    original = sc.compile(root, opts)
    exchanged = sc.compile(restored, restored_opts)

    assert exchanged.duration()[0] == original.duration()[0]
    assert len(exchanged.block_events) == len(original.block_events)


def test_provenance_and_extensions_do_not_change_semantic_hash(opts) -> None:
    root = _event_tree(opts)
    first = sc.export_logicblock(root, opts, provenance={'producer': 'python'})
    second = sc.export_logicblock(
        root, opts, provenance={'producer': 'matlab', 'metadata': {'clock': 'tomorrow'}},
        extensions={'lab:private': {'anything': [1, 2, 3]}},
    )

    assert sc.semantic_hash(first) == sc.semantic_hash(second)


@pytest.mark.parametrize(
    ('mutation', 'path'),
    [
        (lambda doc: doc.pop('version'), '$'),
        (lambda doc: doc['root'].update({'unknown': True}), '$.root'),
        (lambda doc: doc['root']['nodes'][0].update({'start_s': 'NaN'}),
         '$.root.nodes[0].start_s'),
        (lambda doc: doc['root']['nodes'][1]['item']['block']['nodes'][0]['item'][
            'payload'].update({'rise_time_s': '-0.1'}),
         '$.root.nodes[1].item.block.nodes[0].item.payload.rise_time_s'),
    ],
)
def test_malformed_documents_name_a_json_path(mutation, path) -> None:
    document = json.loads((FIXTURES / 'minimal.lb.json').read_text(encoding='utf-8'))
    mutation(document)

    with pytest.raises(sc.ExchangeError) as caught:
        sc.validate_document(document)

    assert caught.value.path == path


def test_write_and_read_are_a_file_round_trip(tmp_path, opts) -> None:
    path = tmp_path / 'tree.lb.json'
    root = _event_tree(opts)

    returned = sc.write_logicblock(path, root, opts, definitions={'NameHint': 'roundtrip'})
    restored, restored_opts, definitions = sc.read_logicblock(path)

    assert returned == path
    assert path.read_text(encoding='utf-8').endswith('\n')
    assert sc.export_logicblock(restored, restored_opts, definitions=definitions) == json.loads(
        path.read_text(encoding='utf-8')
    )
