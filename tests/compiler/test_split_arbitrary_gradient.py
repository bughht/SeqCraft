"""
An arbitrary gradient that changes sign, cut by a block boundary, compiles.

The compiler rebuilds each side of the cut.  It used not to round-trip exactly: where the
waveform's value ``v`` at the seam was below 1e-6 of ``max_grad``, ``axis_gradient`` set it to zero
on both sides, removing a triangle of height ``v`` out to the nearest sample, half a raster away, on
each side, and the emitted area moved by ``-v r / 2``.  When m1 was compared for the whole
sequence, that residual sat in the noise of its relative tolerance: six of ten legal sign-changing
gradients were refused until a floor proportional to peak ``|g|`` was added, and which six could not
be predicted from the sequence.

Only an end value float cannot tell from zero is zeroed now, so the seam value is kept and the area
does not move.  The same ten compile, and the same gradient played one raster late is refused.

Nothing here builds a module.  The property under test is *one arbitrary gradient that changes
sign inside itself, cut by a block boundary*, which spirals happen to meet first.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pypulseq as pp
import pytest
from recheck import compile_capturing, recheck

import seqcraft as sc
from seqcraft.compiler.verification import _sequence_moments
from seqcraft.design.events import knots_of, pwl_moment


def _opts(max_grad: float = 40.0) -> pp.Opts:
    return pp.Opts(max_grad=max_grad, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)


def _sign_changing(n: int, fraction: float, opts: pp.Opts):
    """One arbitrary gradient that runs out, crosses zero several times, and comes home."""
    r = np.linspace(0.0, 1.0, n)
    amp = fraction * opts.max_grad * r * np.sin(16.0 * n / 4000 * np.pi * r)
    return pp.make_arbitrary_grad(channel='x', waveform=amp, system=opts, delay=0.0,
                                  first=0.0, last=0.0)


def _split_tree(grad, n: int, opts: pp.Opts) -> sc.LogicBlock:
    """The same gradient, cut across a block boundary so the compiler has to rebuild it."""
    duration = n * opts.grad_raster_time
    return sc.LogicBlock('sign_changing').add(0.0, grad).add(duration / 2, sc.barrier())


@pytest.mark.parametrize('n', [1000, 2000, 4000, 8000, 16000])
@pytest.mark.parametrize('fraction', [0.3, 0.9])
def test_a_legal_sign_changing_gradient_compiles(n: int, fraction: float) -> None:
    """
    The regression, stated as the thing that was broken: these are all legal, and all were not.

    Six of these ten combinations were refused before the global m1 had a floor.
    """
    opts = _opts()
    sc.compile(_split_tree(_sign_changing(n, fraction, opts), n, opts), opts, name='legal')


@pytest.mark.parametrize('max_grad', [20.0, 40.0])
@pytest.mark.parametrize(('n', 'fraction'), [(2000, 0.2), (2000, 0.5), (8000, 0.2), (8000, 0.5)])
def test_a_small_genuine_seam_value_is_kept(max_grad: float, n: int, fraction: float) -> None:
    """
    Each of these seams sits below 1e-6 of ``max_grad``, where the old snap zeroed it.

    It is a genuine value, far above any float residue, so both blocks carry it and the emitted
    area moves by nothing beyond float -- not by the ``-v r / 2`` the snap used to remove.
    """
    opts = _opts(max_grad)
    grad = _sign_changing(n, fraction, opts)
    times, amps = knots_of(grad, 0.0)
    r = float(opts.grad_raster_time)
    # In ticks, as the compiler interpolates: at float seconds the value itself carries ~1e-10 Hz/m.
    ticks = np.rint(times * 1e12)
    v = float(np.interp(np.rint(n * r / 2 * 1e12), ticks, amps))
    assert 0.0 < abs(v) < 1e-6 * float(opts.max_grad), 'the old snap must have zeroed this one'

    seq, _ = compile_capturing(_split_tree(grad, n, opts), opts, check=False, name='kept')
    first, second = (seq.get_block(i).gx for i in sorted(seq.block_events))
    assert float(first.last) == pytest.approx(v, rel=1e-9)
    assert float(second.first) == pytest.approx(v, rel=1e-9)
    moved = _sequence_moments(seq, 0)['x'] - pwl_moment(times, amps, 0)
    assert abs(moved) < 1e-3 * abs(v) * r / 2


@pytest.mark.parametrize(('n', 'fraction'), [(1000, 0.3), (16000, 0.9)])
def test_a_split_gradient_played_one_raster_late_is_refused(n: int, fraction: float) -> None:
    """Modelling the snap does not blind the check: one raster late is far outside it."""
    opts = _opts()
    _, seen = compile_capturing(_split_tree(_sign_changing(n, fraction, opts), n, opts), opts,
                                name='late')
    dt = float(opts.grad_raster_time)
    late = [dataclasses.replace(p, node_t=p.node_t + dt, start=p.start + dt, end=p.end + dt,
                                res_start=p.res_start + dt, res_end=p.res_end + dt)
            if p.kind == 'grad' else p for p in seen['placed']]
    with pytest.raises(sc.CompilerContractError, match='compiled m1'):
        recheck(seen, placed=late)
