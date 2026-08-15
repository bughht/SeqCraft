"""
Every event type is either emitted or rejected -- never dropped.

The bug this guards was silent: `rot3D`, `soft_delay` and `rf_shim` matched none of the
compiler's positive branches, so they were placed with a zero-width reservation, collected by
nothing, and vanished.  A tree carrying a rotation extension compiled clean, reported no
issues, and played unrotated.
"""

from __future__ import annotations

import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.core._compiler.model import HANDLED_KINDS
from seqcraft.core._compiler.placement import UNSUPPORTED_KINDS


def test_a_rotation_extension_is_rejected_not_dropped(opts) -> None:
    """The headline case: this used to compile clean and play unrotated."""
    rot = pp.make_rotation(np.pi / 4)
    tree = sc.LogicBlock('dwi').add(0.0, pp.make_trapezoid('x', area=100.0, system=opts))
    tree.add(0.0, rot)
    with pytest.raises(sc.CompileError, match='rotation extension'):
        sc.compile(tree, opts)


def test_the_rotation_error_says_how_to_get_the_rotation(opts) -> None:
    """An error that only says 'no' would leave a real sequence with nowhere to go."""
    tree = sc.LogicBlock('dwi').add(0.0, pp.make_rotation(np.pi / 4))
    tree.add(0.0, pp.make_trapezoid('x', area=100.0, system=opts))
    with pytest.raises(sc.CompileError) as err:
        sc.compile(tree, opts)
    text = str(err.value)
    assert 'bakes rotations' in text
    assert 'rotate_3d' in text or 'rotate()' in text
    assert 'dwi' in text, 'the message must name where in the tree it came from'


@pytest.mark.parametrize('kind', sorted(UNSUPPORTED_KINDS))
def test_every_unsupported_type_names_itself_and_offers_a_way_round(kind) -> None:
    """Each entry must carry a description and at least one hint, or the error is useless."""
    what, hints = UNSUPPORTED_KINDS[kind]
    assert what, f'{kind} has no description'
    assert not what.endswith('.'), 'the description is a phrase, not a sentence'
    assert hints, f'{kind} has no hints'
    assert all(h and h[0].islower() for h in hints), 'hints read as continuations, lowercase'


def test_a_soft_delay_is_rejected(opts) -> None:
    sd = pp.make_soft_delay(numID=0, hint='TE', default_duration=1e-3)
    tree = sc.LogicBlock('t').add(0.0, sd)
    with pytest.raises(sc.CompileError, match='soft-delay'):
        sc.compile(tree, opts)


def test_a_typo_type_is_rejected_with_the_handled_list(opts) -> None:
    """LogicBlock.add() accepts anything with a .type, so a hand-built namespace gets here."""
    bogus = SimpleNamespace(type='trapezoid', channel='x', delay=0.0)   # 'trap', not 'trapezoid'
    tree = sc.LogicBlock('t').add(0.0, bogus)
    with pytest.raises(sc.CompileError) as err:
        sc.compile(tree, opts)
    text = str(err.value)
    assert "unknown event type 'trapezoid'" in text
    assert 'trap' in text, 'the handled list must be shown so the fix is obvious'


def test_the_rejection_names_the_nested_tag_path(opts) -> None:
    inner = sc.LogicBlock('spoiler').add(0.0, pp.make_rotation(0.5))
    tree = sc.LogicBlock('tr').add(0.0, inner)
    tree.add(0.0, pp.make_trapezoid('z', area=100.0, system=opts))
    with pytest.raises(sc.CompileError, match=r'tr\.spoiler'):
        sc.compile(tree, opts)


# --------------------------------------------------------------------------------- coverage
def _pypulseq_event_types() -> set[str]:
    """
    Every ``.type`` literal pypulseq assigns, read out of its source.

    Introspective on purpose: when a pypulseq upgrade adds a thirteenth event type, this test
    fails and someone decides where it belongs.  Enumerating them by hand would silently rot.
    """
    root = Path(pp.__file__).parent
    found: set[str] = set()
    for path in root.rglob('*.py'):
        text = path.read_text(encoding='utf-8', errors='ignore')
        found.update(re.findall(r"\.type\s*=\s*'([A-Za-z_0-9]+)'", text))
    return found


def test_the_whitelist_accounts_for_every_pypulseq_event_type() -> None:
    """
    Each pypulseq type is either handled or explicitly unsupported.

    This is the test that would have caught the original bug, and it is the one that catches the
    next pypulseq release adding an event type seqcraft would otherwise drop.
    """
    produced = _pypulseq_event_types()
    expected_core = set(HANDLED_KINDS) - {sc.core.logic.BARRIER}
    assert expected_core <= produced, (
        f'introspection missed handled PyPulseq types {sorted(expected_core - produced)}; '
        'the source regex or installed package layout has changed'
    )
    unaccounted = produced - set(HANDLED_KINDS) - set(UNSUPPORTED_KINDS)
    assert not unaccounted, (
        f'pypulseq can produce {sorted(unaccounted)}, which the compiler neither emits nor '
        f'rejects by name -- they would fall through as "unknown event type". Add each to '
        f'HANDLED_KINDS in model.py or UNSUPPORTED_KINDS in placement.py.'
    )


def test_unsupported_and_handled_do_not_overlap() -> None:
    assert not (set(HANDLED_KINDS) & set(UNSUPPORTED_KINDS))


def test_all_handled_types_actually_compile(opts) -> None:
    """
    The other half of the whitelist's promise: everything on it works.

    A type listed as handled but broken would be worse than one rejected, so each gets a
    minimal tree.
    """
    rf = pp.make_sinc_pulse(flip_angle=0.5, duration=1e-3, system=opts, use='excitation',
                            slice_thickness=5e-3, apodization=0.5, time_bw_product=4)
    cases = {
        'trap': pp.make_trapezoid('x', area=100.0, system=opts),
        'grad': pp.make_extended_trapezoid(
            'y', times=np.array([0.0, 100e-6, 300e-6, 400e-6]),
            amplitudes=np.array([0.0, 1e4, 1e4, 0.0]), system=opts),
        'rf': rf,
        'adc': pp.make_adc(num_samples=64, dwell=10e-6, system=opts),
        'delay': pp.make_delay(1e-3),
        'labelset': pp.make_label('LIN', 'SET', 3),
        'labelinc': pp.make_label('LIN', 'INC', 1),
        'trigger': pp.make_trigger('physio1', duration=200e-6, system=opts),
        'output': pp.make_digital_output_pulse('osc0', duration=200e-6, system=opts),
        sc.core.logic.BARRIER: sc.barrier(),
    }
    assert set(cases) == set(HANDLED_KINDS), (
        f'HANDLED_KINDS and this table disagree: {set(HANDLED_KINDS) ^ set(cases)}'
    )
    for kind, event in cases.items():
        tree = sc.LogicBlock('t').add(0.0, pp.make_delay(2e-3)).add(0.0, event)
        out = sc.compile(tree, opts)          # must not raise
        assert out.n_blocks >= 1, kind
