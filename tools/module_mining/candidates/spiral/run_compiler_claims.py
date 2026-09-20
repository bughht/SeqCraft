"""
Test PR23's two compiler claims with **no Spiral module present**.

The compiler-change escalation in the mining skill requires a candidate-free reproducer before a
compiler change is considered at all.  This is that reproducer.  The structural property both
claims depend on is not "spiral" -- it is one arbitrary gradient whose amplitude changes sign
inside itself, and a spiral is merely the first module in this library that has one.

Findings, at the time of writing, on current ``main``:

    m1  CONFIRMED.  The tolerance straddles the representation floor, and six of ten legal
        waveforms are refused unpredictably -- ``n=2000`` at 90 % passes with an error of 7e-15
        while ``n=4000`` at 90 % fails with 1.2e-7.
    m0  NOT REPRODUCED.  The comment/code mismatch is real, the failure is not.

    python tools/module_mining/candidates/spiral/run_compiler_claims.py
"""

from __future__ import annotations

import json
import logging
import warnings
from pathlib import Path

import numpy as np
import pypulseq as pp

import seqcraft as sc
import seqcraft.compiler as compiler_package
from seqcraft.compiler.verification import _sequence_moments
from seqcraft.design.events import knots_of, pwl_moment


#: PR23's proposed scale.  Reproduced here rather than imported: the point is to test the claim
#: against the compiler as it stands, not to run the branch's own code.
def traversed(times: np.ndarray, amps: np.ndarray) -> float:
    """``integral |g| dt`` over a piecewise-linear waveform, splitting at sign changes."""
    t0, t1 = np.asarray(times[:-1]), np.asarray(times[1:])
    g0, g1 = np.asarray(amps[:-1]), np.asarray(amps[1:])
    h = t1 - t0
    same = g0 * g1 >= 0.0
    span = np.abs(g0 - g1)
    crossing = np.divide(h * (g0**2 + g1**2), 2.0 * span, out=np.zeros_like(h),
                         where=~same & (span > 0.0))
    return float(np.sum(np.where(same, 0.5 * h * np.abs(g0 + g1), crossing)))


def _opts() -> pp.Opts:
    return pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)


def measure(n: int, fraction: float, *, splits: int, opts: pp.Opts) -> dict | None:
    """Return the m0/m1 discrepancy the compiler's own check would see, or None if illegal."""
    r = np.linspace(0.0, 1.0, n)
    amp = fraction * opts.max_grad * r * np.sin(16.0 * n / 4000 * np.pi * r)
    try:
        grad = pp.make_arbitrary_grad(channel='x', waveform=amp, system=opts, delay=0.0,
                                      first=0.0, last=0.0)
    except ValueError:
        return None                      # slew violation: not a legal input, so not evidence

    duration = n * opts.grad_raster_time
    knots = knots_of(grad, 0.0)
    tree_m0, tree_m1 = pwl_moment(*knots, 0), pwl_moment(*knots, 1)
    net, total = abs(tree_m0), traversed(*knots)

    tree = sc.LogicBlock('sign_changing').add(0.0, grad)
    for index in range(1, splits):
        tree.add(index * duration / splits, sc.barrier())

    # The check is what is under test, so it must not abort the measurement.
    original = compiler_package.verify_against_tree
    compiler_package.verify_against_tree = lambda *a, **k: None
    try:
        seq = sc.compile(tree, opts, name='repro')
    finally:
        compiler_package.verify_against_tree = original

    scale_net = max(net, 1.0)
    scale_total = max(total, 1.0)
    horizon = max(duration, 1e-3)
    return {
        'knots': n, 'amp_fraction': fraction, 'splits': splits,
        'net_per_m': net, 'traversed_per_m': total,
        'm0_error': abs(_sequence_moments(seq, 0)['x'] - tree_m0),
        'm1_error': abs(_sequence_moments(seq, 1)['x'] - tree_m1),
        'm0_tol_current': 1e-6 * scale_net,
        'm0_tol_traversed': 1e-6 * scale_total,
        'm1_tol_current': 1e-9 * max(scale_net * horizon, 1.0),
        'm1_tol_traversed_1e9': 1e-9 * max(scale_total * horizon, 1.0),
        'm1_tol_traversed_1e7': 1e-7 * max(scale_total * horizon, 1.0),
        'm1_one_raster_signal': total * opts.grad_raster_time,
    }


def main() -> None:
    """Sweep knot count and amplitude, and report which legal waveforms the compiler refuses."""
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    warnings.filterwarnings('ignore')
    opts = _opts()

    rows = [row for n in (1000, 2000, 4000, 8000, 16000) for fraction in (0.3, 0.9)
            if (row := measure(n, fraction, splits=2, opts=opts)) is not None]

    logging.info('%6s %5s %12s %14s %13s %13s %12s', 'knots', 'amp', 'm1 error',
                 'current tol', '1e-9 x trav', '1e-7 x trav', 'one raster')
    refused = 0
    for row in rows:
        bad = row['m1_error'] > row['m1_tol_current']
        refused += bad
        logging.info('%6d %4.0f%% %12.4g %14.4g %13.4g %13.4g %12.4g  %s',
                     row['knots'], row['amp_fraction'] * 100, row['m1_error'],
                     row['m1_tol_current'], row['m1_tol_traversed_1e9'],
                     row['m1_tol_traversed_1e7'], row['m1_one_raster_signal'],
                     'REFUSED' if bad else 'accepted')
    logging.info('')
    logging.info('m1: %d of %d legal waveforms refused by the current tolerance', refused,
                 len(rows))
    logging.info('m0: max error %.4g against a minimum tolerance of %.4g -- no failure observed',
                 max(r['m0_error'] for r in rows), min(r['m0_tol_current'] for r in rows))

    out = Path(__file__).with_name('compiler_claims_result.json')
    out.write_text(json.dumps(rows, indent=2))
    logging.info('-> %s', out)


if __name__ == '__main__':
    main()
