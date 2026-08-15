"""
Labels address ADCs, so the compiler must place them by target, not by containment.

A pulseq label is a running register: the interpreter applies a block's labels on reaching that
block, and an ADC in the same block then samples the state.  So a label may live in any block
strictly after the previous ADC's and at or before its target ADC's -- all equivalent.  The
compiler therefore *chooses*, and choosing by containment was wrong: a boundary pushed later
than an ADC's reservation end (which the gap-midpoint fallback does whenever a gradient covers
the natural candidates) put a label into the previous readout's block and overwrote its k-space
address.

That failure was silent past the duplicate-address check whenever the shifted addresses stayed
unique, which is the normal case for a real acquisition.
"""

from __future__ import annotations

import warnings

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc


def _lin(compiled) -> list[int]:
    """LIN as seen by each ADC, in acquisition order."""
    labels = compiled.seq.evaluate_labels(evolution='adc')
    return [int(v) for v in np.atleast_1d(np.asarray(labels['LIN']))]


def _blocking_gradient(opts):
    """
    A gradient covering both natural boundary candidates between two ADCs.

    This is what forces the gap-midpoint fallback, which is what moved the boundary past the
    first ADC's reservation end and caused the corruption.
    """
    return pp.make_trapezoid('x', amplitude=0.15 * opts.max_grad, duration=4700e-6,
                             rise_time=200e-6, system=opts)


def test_a_label_between_two_adcs_reaches_the_later_one(opts) -> None:
    """
    The measured corruption: LIN was [7, 7] instead of [1, 7].

    The label sits at 1000 us, comfortably after the first ADC stops sampling at 650 us, yet the
    boundary lands at 2830 us -- so containment put it in the first ADC's block.
    """
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = (
        sc.LogicBlock('t')
        .add(0.0, adc, pp.make_label('LIN', 'SET', 1))
        .add(500e-6, _blocking_gradient(opts))
        .add(1000e-6, pp.make_label('LIN', 'SET', 7))
        .add(5000e-6, adc)
    )
    out = sc.compile(tree, opts)
    assert _lin(out) == [1, 7], (
        'the label at 1000 us addresses the second ADC, not the first'
    )


def test_the_result_does_not_depend_on_where_the_boundary_lands(opts) -> None:
    """
    The property that makes this correct rather than merely fixed.

    The same tree with and without the gradient that moves the boundary must label identically.
    Under containment these two disagreed, which is the definition of the bug.
    """
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)

    def build(*, blocker: bool):
        tree = sc.LogicBlock('t').add(0.0, adc, pp.make_label('LIN', 'SET', 1))
        if blocker:
            tree.add(500e-6, _blocking_gradient(opts))
        tree.add(1000e-6, pp.make_label('LIN', 'SET', 7))
        tree.add(5000e-6, adc)
        return sc.compile(tree, opts)

    assert _lin(build(blocker=False)) == _lin(build(blocker=True)) == [1, 7]


@pytest.mark.parametrize('label_t_us', [660, 1000, 2000, 2830, 3000, 4990])
def test_any_time_in_the_gap_gives_the_same_answer(opts, label_t_us) -> None:
    """
    Every time strictly between the two reservations is equivalent, so all must agree.

    660 us is the instant the first reservation ends and 2830 us is the chosen boundary -- the
    two places containment and target assignment used to differ most.
    """
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = (
        sc.LogicBlock('t')
        .add(0.0, adc, pp.make_label('LIN', 'SET', 1))
        .add(500e-6, _blocking_gradient(opts))
        .add(label_t_us * 1e-6, pp.make_label('LIN', 'SET', 7))
        .add(5000e-6, adc)
    )
    assert _lin(sc.compile(tree, opts)) == [1, 7]


def test_a_label_at_the_same_time_as_its_adc_still_reaches_it(opts) -> None:
    """The ordinary case, and the one the previous behaviour got right."""
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    tree = sc.LogicBlock('t').add(0.0, g).add(1e-3, adc, pp.make_label('LIN', 'SET', 5))
    assert _lin(sc.compile(tree, opts)) == [5]


def test_a_label_before_the_first_adc_reaches_it(opts) -> None:
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = (
        sc.LogicBlock('t')
        .add(0.0, pp.make_label('LIN', 'SET', 4))
        .add(0.0, pp.make_trapezoid('x', area=100.0, system=opts))
        .add(2e-3, adc)
    )
    assert _lin(sc.compile(tree, opts)) == [4]


def test_several_labelincs_on_one_key_commute(opts) -> None:
    """
    Addition commutes, so grouping two increments onto one block is safe.

    This matters because pypulseq gives no control over intra-block order (see
    :func:`test_order_dependent_labels_for_one_readout_are_rejected`), so "safe" has to mean
    "order-independent", not "emitted in the right order".
    """
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = (
        sc.LogicBlock('t')
        .add(0.0, adc, pp.make_label('LIN', 'SET', 0))
        .add(2000e-6, pp.make_label('LIN', 'INC', 3))
        .add(2500e-6, pp.make_label('LIN', 'INC', 4))
        .add(5000e-6, adc)
    )
    assert _lin(sc.compile(tree, opts)) == [0, 7]


def test_labels_on_different_keys_share_a_block_safely(opts) -> None:
    """Different keys do not interact, so order is irrelevant and grouping is free."""
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = (
        sc.LogicBlock('t')
        .add(0.0, adc)
        .add(2000e-6, pp.make_label('LIN', 'SET', 5))
        .add(2500e-6, pp.make_label('SLC', 'SET', 2))
        .add(5000e-6, adc)
    )
    out = sc.compile(tree, opts)
    labels = out.seq.evaluate_labels(evolution='adc')
    assert [int(v) for v in np.atleast_1d(np.asarray(labels['LIN']))] == [0, 5]
    assert [int(v) for v in np.atleast_1d(np.asarray(labels['SLC']))] == [0, 2]


def test_order_dependent_labels_for_one_readout_are_rejected(opts) -> None:
    """
    pypulseq discards intra-block label order, so this intent cannot be expressed.

    ``Sequence/block.py`` sorts a block's extensions by library reference id -- "we rely on the
    sorting of the extension IDs" -- so ``add_block(set, inc)`` and ``add_block(inc, set)``
    build the identical block.  Verified directly: both yield LIN=10, never 13.  Guessing would
    silently pick one of two different k-space addressings, so it is an error.
    """
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = (
        sc.LogicBlock('t')
        .add(0.0, adc, pp.make_label('LIN', 'SET', 0))
        .add(2000e-6, pp.make_label('LIN', 'SET', 10))
        .add(2500e-6, pp.make_label('LIN', 'INC', 3))
        .add(5000e-6, adc)
    )
    with pytest.raises(sc.CompileError) as err:
        sc.compile(tree, opts)
    text = str(err.value)
    assert 'order cannot be expressed' in text
    assert 'barrier' in text, 'the message must give a way to express it'


def test_a_barrier_makes_an_order_dependent_pair_expressible(opts) -> None:
    """
    The remedy the error suggests has to actually work.

    In separate blocks, block order determines application order -- which pulseq does define.
    """
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = (
        sc.LogicBlock('t')
        .add(0.0, adc, pp.make_label('LIN', 'SET', 0))
        .add(2000e-6, pp.make_label('LIN', 'SET', 10))
        .add(2250e-6, sc.barrier('between_labels'))
        .add(2500e-6, pp.make_label('LIN', 'INC', 3))
        .add(5000e-6, adc)
    )
    out = sc.compile(tree, opts)
    assert _lin(out) == [0, 13], 'SET 10 then INC 3, ordered by block'


def test_a_label_with_no_following_adc_is_emitted_and_reported(opts) -> None:
    """
    Nothing to address, so it keeps containment placement -- but it must not be lost, and the
    user has to be told, because containment can put it in the *preceding* readout's block and
    change that address.
    """
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = (
        sc.LogicBlock('t')
        .add(0.0, adc, pp.make_label('LIN', 'SET', 2))
        .add(3000e-6, pp.make_label('LIN', 'SET', 99))
        .add(3000e-6, pp.make_delay(1e-3))
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        out = sc.compile(tree, opts)
    emitted = [
        lab.value
        for i in sorted(out.seq.block_events)
        for lab in (getattr(out.seq.get_block(i), 'label', None) or [])
    ]
    assert 99 in emitted, 'a trailing label must still reach the file'
    orphan = [
        str(w.message) for w in caught
        if issubclass(w.category, sc.SeqCraftWarning) and 'no ADC after them' in str(w.message)
    ]
    assert orphan, f'expected a warning about the orphan label, got {[str(w.message) for w in caught]}'
    assert 'LIN' in orphan[0], 'the warning must name the label'


def test_a_multislice_address_survives_a_blocking_gradient(opts) -> None:
    """
    The realistic version: unique-but-shifted addresses pass the duplicate check.

    This is why the corruption was silent -- LIN 1,2,3 shifted by one collides with nothing, so
    the duplicate-address check that exists to catch mis-addressing never fires.  Labels are
    written before the readout they address, which is the correct idiom.
    """
    adc = pp.make_adc(num_samples=64, dwell=10e-6, system=opts)
    tree = sc.LogicBlock('t')
    period = 6000e-6
    for i in range(4):
        t0 = i * period
        tree.add(t0 + 1000e-6, pp.make_label('LIN', 'SET', i + 1))
        tree.add(t0 + 1500e-6, _blocking_gradient(opts))
        tree.add(t0 + 3000e-6, adc)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        out = sc.compile(tree, opts)
    assert _lin(out) == [1, 2, 3, 4], 'each readout gets the label written before it'
    assert not [w for w in caught if 'no ADC after them' in str(w.message)], (
        'no orphan warnings for a well-formed sequence'
    )
