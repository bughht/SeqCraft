"""
The data model: two attributes, one method, and a recursive walk.

These tests are short because the model is small.  That is the point of the design -- the
difficulty lives in the compiler, which has its own directory.
"""

from __future__ import annotations

import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.core.logic import flatten, span


# ------------------------------------------------------------------------------------- create
def test_add_accumulates_and_returns_self(opts) -> None:
    """`add` appends and returns the same block, so calls may chain."""
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    lb = sc.LogicBlock('t')
    assert lb.add(0.0, g) is lb
    lb.add(1e-3, g, g)
    assert len(lb) == 3
    assert [n.start for n in lb] == [0.0, 1e-3, 1e-3]


def test_add_takes_several_events_at_one_time(opts) -> None:
    """Several events at one start reads like ``add_block``, which is the point."""
    rf, gz, _ = pp.make_sinc_pulse(
        flip_angle=0.26, duration=1e-3, slice_thickness=5e-3, system=opts,
        return_gz=True, use='excitation',
    )
    lb = sc.LogicBlock('exc').add(0.0, rf, gz)
    assert [n.item.type for n in lb] == ['rf', 'trap']


def test_add_rejects_a_non_event() -> None:
    """A wrong type fails where it was added, not deep inside the compiler."""
    with pytest.raises(sc.ConfigurationError, match='takes pulseq events or LogicBlocks'):
        sc.LogicBlock('t').add(0.0, 42)


# --------------------------------------------------------------------------------------- read
def test_duration_is_measured_not_declared(opts) -> None:
    """A block cannot advertise a length that disagrees with what it plays."""
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    lb = sc.LogicBlock('t').add(0.0, g)
    assert lb.duration == pytest.approx(pp.calc_duration(g))
    lb.add(1e-3, g)
    assert lb.duration == pytest.approx(1e-3 + pp.calc_duration(g))


def test_duration_tracks_mutation(opts) -> None:
    """Editing a node's start changes the duration, with no cache to go stale."""
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    lb = sc.LogicBlock('t').add(0.0, g)
    before = lb.duration
    lb.nodes[0].start += 1e-3
    assert lb.duration == pytest.approx(before + 1e-3)
    del lb.nodes[0]
    assert lb.duration == 0.0


def test_a_delay_pads_a_block_with_no_gradients() -> None:
    """The b=0 diffusion case: no gradients, same length as an encoded lobe."""
    lb = sc.LogicBlock('b0').add(0.0, pp.make_delay(4.2e-3))
    assert lb.duration == pytest.approx(4.2e-3)
    assert not [n for n in lb if n.item.type in ('trap', 'grad')]


def test_nested_duration_is_recursive(opts) -> None:
    """A parent's duration accounts for where its children start."""
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    inner = sc.LogicBlock('inner').add(0.0, g)
    outer = sc.LogicBlock('outer').add(0.0, inner).add(2e-3, inner)
    assert outer.duration == pytest.approx(2e-3 + inner.duration)
    deeper = sc.LogicBlock('deeper').add(5e-3, outer)
    assert deeper.duration == pytest.approx(5e-3 + outer.duration)


def test_span_of_an_event_includes_its_own_delay(opts) -> None:
    """``calc_duration`` already covers the delay, so `span` must not add it twice."""
    g = pp.make_trapezoid('x', area=100.0, delay=200e-6, system=opts)
    assert span(g) == pytest.approx(pp.calc_duration(g))
    assert span(g) > 200e-6


def test_repr_is_one_readable_line(opts) -> None:
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    assert repr(sc.LogicBlock('ro').add(0.0, g)) == 'LogicBlock(ro, 1 node, 0.26 ms)'
    assert '2 nodes' in repr(sc.LogicBlock('ro').add(0.0, g, g))


def test_describe_renders_the_tree(opts) -> None:
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    inner = sc.LogicBlock('pre').add(0.0, g)
    text = sc.LogicBlock('ro').add(0.0, inner).add(1e-3, g).describe()
    assert 'ro' in text
    assert 'pre' in text
    assert 'trap x' in text


# ------------------------------------------------------------------------------------- update
def test_copy_shares_items_but_not_the_list(opts) -> None:
    """The answer to the one hazard mutability introduces."""
    g = pp.make_trapezoid('x', area=100.0, system=opts)
    a = sc.LogicBlock('a').add(0.0, g)
    b = a.copy()
    assert b.nodes is not a.nodes
    assert b.nodes[0].item is a.nodes[0].item
    b.add(1e-3, g)
    assert len(a) == 1


def test_the_same_block_placed_twice_is_shared_not_copied(opts) -> None:
    """1700 references to one spoiler, not 1700 copies -- the reason composition is O(1)."""
    g = pp.make_trapezoid('z', area=500.0, system=opts)
    spoil = sc.LogicBlock('spoil').add(0.0, g)
    root = sc.LogicBlock('seq')
    for i in range(50):
        root.add(i * 10e-3, spoil)
    assert len({id(n.item) for n in root}) == 1
    assert len(list(flatten(root))) == 50


# ------------------------------------------------------------------------------------ flatten
def test_flatten_sums_starts_through_three_levels(opts) -> None:
    """The compiler's only entry point: a recursive sum of node starts."""
    g = pp.make_trapezoid('x', area=10.0, system=opts)
    deep = sc.LogicBlock('a').add(
        0.0, sc.LogicBlock('b').add(1e-3, sc.LogicBlock('c').add(2e-3, g))
    )
    (start, event, path), = flatten(deep)
    assert start == pytest.approx(3e-3)
    assert event is g
    assert path == ('a', 'b', 'c')


def test_flatten_skips_untagged_levels(opts) -> None:
    """An untagged block contributes no path element, so tags stay optional."""
    g = pp.make_trapezoid('x', area=10.0, system=opts)
    tree = sc.LogicBlock('outer').add(0.0, sc.LogicBlock().add(0.0, g))
    (_, _, path), = flatten(tree)
    assert path == ('outer',)


def test_flatten_does_not_add_the_event_delay(opts) -> None:
    """The compiler adds it, because an RF's delay carries the transmit dead time."""
    g = pp.make_trapezoid('x', area=10.0, delay=500e-6, system=opts)
    (start, _, _), = flatten(sc.LogicBlock('t').add(1e-3, g))
    assert start == pytest.approx(1e-3)


def test_barrier_occupies_no_time() -> None:
    b = sc.barrier('mid')
    assert span(b) == 0.0
    assert sc.LogicBlock('t').add(5e-3, b).duration == pytest.approx(5e-3)


# ------------------------------------------------------------- what the model deliberately lacks
@pytest.mark.parametrize('name', ['marks', 'labels', 'then', 'over', 'chain', 'stack', 'scaled'])
def test_model_stays_minimal(name: str) -> None:
    """
    Guards the design decision, not just the code.

    Each of these was in an earlier draft and was removed for a reason: marks because a block
    cannot know where it sits, labels because a label is an event, and the combinators because
    ``add`` is shorter than any algebra.  If one comes back it should be a deliberate decision
    with this test deleted, not an accident.
    """
    assert not hasattr(sc.LogicBlock('t'), name)


def test_logic_block_has_exactly_two_attributes() -> None:
    """``__slots__`` states the whole model: a tag and a list of nodes."""
    assert set(sc.LogicBlock.__slots__) == {'tag', 'nodes'}
