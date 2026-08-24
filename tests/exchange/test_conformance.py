"""Vertical LBTX slices compare physical semantics, not JSON or .seq bytes alone."""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc

FIXTURES = Path(__file__).parents[1] / 'fixtures' / 'lbtx'


@pytest.mark.parametrize(
    ('name', 'expected_hash', 'duration_s'),
    [
        ('cartesian.lb.json', '87dfb2a19cf14449d71316a48008fe4e6e16525b72d9128cdd3c4031690b8970',
         3.56e-3),
        ('spiral.lb.json', '4034cfe990b2be296fc7e4b2a0ce1687f31962be9ca4a82ed9c082008f879b5b',
         0.98e-3),
    ],
)
def test_versioned_vertical_golden_fixtures(name, expected_hash, duration_s) -> None:
    document = json.loads((FIXTURES / name).read_text(encoding='utf-8'))
    root, opts, definitions = sc.import_logicblock(document)

    sequence = _compiled(root, opts, definitions)

    assert sc.semantic_hash(document) == expected_hash
    assert sequence.duration()[0] == pytest.approx(duration_s)
    assert all(key in sequence.definitions for key in definitions)


def _roundtrip(root, opts):
    return sc.import_logicblock(sc.export_logicblock(root, opts))[0:2]


def _compiled(root, opts, definitions=None):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sc.SeqCraftWarning)
        return sc.compile(root, opts, definitions=definitions)


def test_cartesian_line_preserves_samples_echo_and_compiled_structure(opts) -> None:
    readout = sc.modules.CartesianLine(
        opts=opts, fov_mm=220.0, matrix=64, bandwidth_hz_px=312.5,
    )
    original = sc.LogicBlock('cartesian').add(0.0, readout())
    restored, restored_opts = _roundtrip(original, opts)

    k_original = sc.kspace(original, opts)['k_adc'][0]
    k_restored = sc.kspace(restored, restored_opts)['k_adc'][0]
    seq_original = _compiled(original, opts)
    seq_restored = _compiled(restored, restored_opts)

    assert np.array_equal(k_restored, k_original)
    assert k_restored[readout.pre_echo_samples] == pytest.approx(0.0, abs=1e-6)
    assert seq_restored.duration()[0] == seq_original.duration()[0]
    assert len(seq_restored.block_events) == len(seq_original.block_events)


def _spiral_tree(opts):
    count = 98
    angle = np.linspace(0.0, 6.0 * np.pi, count)
    envelope = np.sin(np.linspace(0.0, np.pi, count))
    amplitude = 0.015 * opts.max_grad
    gx = pp.make_arbitrary_grad(
        'x', amplitude * envelope * np.cos(angle), first=0.0, last=0.0, system=opts,
    )
    gy = pp.make_arbitrary_grad(
        'y', amplitude * envelope * np.sin(angle), first=0.0, last=0.0, system=opts,
    )
    adc = pp.make_adc(num_samples=count - 2, dwell=opts.grad_raster_time,
                      delay=opts.grad_raster_time, system=opts)
    return sc.LogicBlock('spiral').add(0.0, gx, gy, adc)


def test_spiral_preserves_waveforms_adc_times_moments_and_kspace(opts) -> None:
    original = _spiral_tree(opts)
    restored, restored_opts = _roundtrip(original, opts)

    original_events = [event for _, event, _ in sc.flatten(original)]
    restored_events = [event for _, event, _ in sc.flatten(restored)]
    original_gradients = [event for event in original_events if event.type == 'grad']
    restored_gradients = [event for event in restored_events if event.type == 'grad']

    for before, after in zip(original_gradients, restored_gradients, strict=True):
        assert np.array_equal(after.waveform, before.waveform)
        assert np.array_equal(after.tt, before.tt)
        assert after.first == before.first
        assert after.last == before.last
        assert after.area == before.area

    assert sc.moments(restored, order=0) == pytest.approx(sc.moments(original, order=0))
    assert sc.moments(restored, order=1) == pytest.approx(sc.moments(original, order=1))
    assert np.array_equal(
        sc.kspace(restored, restored_opts)['t_adc'], sc.kspace(original, opts)['t_adc'],
    )
    assert np.allclose(
        sc.kspace(restored, restored_opts)['k_adc'], sc.kspace(original, opts)['k_adc'],
        rtol=0.0, atol=1e-12,
    )
