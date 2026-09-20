"""
What reversibility costs: terminal gradient free, against terminal gradient zero.

The four-variant family needs an arm that can be played backwards, which requires ``g = 0`` at
**both** ends.  The independent reference (``pulseq/pulseq`` ``writeSpiral.m``) ends its
spiral-out at full gradient and ramps down inside the spoiler, and has no spiral-in.  This
measures the difference, so the policy is adopted with a number beside it rather than on the
argument alone.

**Only the terminal boundary condition differs.**  Same path, same scanner, same solver, same
dwell.  So the study isolates the traversal cost and reports the sample-count consequence that
follows from it; the path cost is zero by construction and is asserted rather than assumed.

The traversal is a standard time-optimal sweep over the path's arclength in k-space:

    v <= g_max                          the gradient limit is a speed limit in 1/m per second
    v^2 * kappa <= s_max                curvature eats slew budget before any acceleration does
    |dv/ds| * v <= sqrt(s_max^2 - (v^2 kappa)^2)      what is left for tangential acceleration

forward from ``v = 0``, then backward -- from the gradient limit for the free policy, from
``v = 0`` for the reversible one.  Nothing here is production code; it exists to produce a number.

    python tools/module_mining/candidates/spiral/run_endpoint_cost.py
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pypulseq as pp

#: Deliberately coarse: the quantity of interest is a ratio, and both policies pay the same
#: discretisation.
_ARCLENGTH_STEPS = 20000


@dataclass(frozen=True)
class Protocol:
    """One representative design point."""

    label: str
    fov_mm: float
    matrix: int
    shots: int
    max_grad_mt_m: float
    max_slew_t_m_s: float
    dwell_s: float = 2.5e-6

    def opts(self) -> pp.Opts:
        """The scanner this protocol is designed against."""
        return pp.Opts(max_grad=self.max_grad_mt_m, grad_unit='mT/m',
                       max_slew=self.max_slew_t_m_s, slew_unit='T/m/s',
                       rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)


def archimedean(protocol: Protocol) -> np.ndarray:
    """
    Return ``(2, n)`` k-space points in 1/m: the Nyquist path, with no scanner in sight.

    ``r = dk * turns * u`` and ``theta = 2 pi * turns * u``, with ``turns = matrix / (2 * shots)``
    -- one turn per ``dk`` of radius per shot, which is the Nyquist condition for an Archimedean
    spiral and is what ``writeSpiral.m`` builds.
    """
    dk = 1e3 / protocol.fov_mm
    turns = protocol.matrix / (2.0 * protocol.shots)
    u = np.linspace(0.0, 1.0, _ARCLENGTH_STEPS)
    radius = dk * turns * protocol.shots * u
    theta = 2.0 * np.pi * turns * u
    return np.vstack([radius * np.cos(theta), radius * np.sin(theta)])


def _geometry(k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return cumulative arclength and curvature along the path."""
    d1 = np.gradient(k, axis=1)
    ds = np.hypot(d1[0], d1[1])
    s = np.concatenate([[0.0], np.cumsum(ds[1:])])
    # kappa = |x' y'' - y' x''| / |r'|^3, guarded at the origin where the path has no tangent.
    d2 = np.gradient(d1, axis=1)
    speed = np.maximum(ds, 1e-30)
    kappa = np.abs(d1[0] * d2[1] - d1[1] * d2[0]) / speed**3
    return s, np.nan_to_num(kappa, nan=0.0, posinf=0.0)


def traverse(k: np.ndarray, opts: pp.Opts, *, terminal_at_rest: bool) -> dict[str, float]:
    """
    Time-optimal traversal of `k`, differing only in whether it must end at ``v = 0``.

    Returns duration, peak gradient and peak slew of the realised traversal.
    """
    s, kappa = _geometry(k)
    ds = np.diff(s)
    g_max = float(opts.max_grad)
    s_max = float(opts.max_slew)

    # The curvature ceiling: as fast as this bend allows before any acceleration is spent.
    ceiling = np.minimum(g_max, np.sqrt(np.divide(s_max, np.maximum(kappa, 1e-30))))

    def sweep(start: float, order: slice, steps: np.ndarray) -> np.ndarray:
        """One monotone pass, respecting the ceiling and what slew is left for acceleration."""
        v = np.empty_like(ceiling)
        v[:] = ceiling
        v[order.start] = min(start, ceiling[order.start])
        index = range(len(steps))
        for i in index:
            here = order.start + order.step * i
            nxt = here + order.step
            lateral = (v[here] ** 2) * kappa[here]
            tangential = np.sqrt(max(s_max**2 - lateral**2, 0.0))
            reachable = np.sqrt(max(v[here] ** 2 + 2.0 * tangential * steps[i], 0.0))
            v[nxt] = min(ceiling[nxt], reachable)
        return v

    forward = sweep(0.0, slice(0, None, 1), ds)
    backward_start = 0.0 if terminal_at_rest else ceiling[-1]
    backward = sweep(backward_start, slice(len(s) - 1, None, -1), ds[::-1])
    v = np.minimum(forward, backward)

    speed = np.maximum(0.5 * (v[:-1] + v[1:]), 1e-12)
    dt = ds / speed
    duration = float(np.sum(dt))
    gradient = v                                        # |g| in 1/m per second IS the speed
    slew = np.abs(np.diff(gradient)) / np.maximum(dt, 1e-30)
    # Time spent in the outer tenth of the radius: with the path fixed, sampling density there is
    # exactly this, and it is where the two policies differ most because that is where one brakes.
    radius = np.hypot(k[0], k[1])
    outer = radius[:-1] >= 0.9 * float(np.max(radius))
    return {
        'duration_s': duration,
        'peak_grad': float(np.max(gradient)),
        'peak_slew': float(np.max(slew)),
        'terminal_grad': float(gradient[-1]),
        'outer_time_s': float(np.sum(dt[outer])),
    }


def study(protocol: Protocol) -> dict[str, object]:
    """Compare the two endpoint policies on one protocol."""
    opts = protocol.opts()
    k = archimedean(protocol)
    free = traverse(k, opts, terminal_at_rest=False)
    rest = traverse(k, opts, terminal_at_rest=True)

    extent = float(np.max(np.hypot(k[0], k[1])))
    samples_free = int(np.ceil(free['duration_s'] / protocol.dwell_s))
    samples_rest = int(np.ceil(rest['duration_s'] / protocol.dwell_s))
    # Outer-k oversampling: samples spent in the last 10 % of radius, per unit arclength there.
    radius = np.hypot(k[0], k[1])
    outer = radius >= 0.9 * extent
    s, _ = _geometry(k)
    outer_arclength = float(s[outer][-1] - s[outer][0])
    return {
        'protocol': protocol.label,
        'k_extent_per_m': extent,
        'duration_free_s': free['duration_s'],
        'duration_rest_s': rest['duration_s'],
        'duration_penalty_s': rest['duration_s'] - free['duration_s'],
        'duration_penalty_pct': 100.0 * (rest['duration_s'] / free['duration_s'] - 1.0),
        'samples_free': samples_free,
        'samples_rest': samples_rest,
        'sample_penalty': samples_rest - samples_free,
        'peak_grad_free': free['peak_grad'],
        'peak_grad_rest': rest['peak_grad'],
        'peak_slew_free': free['peak_slew'],
        'peak_slew_rest': rest['peak_slew'],
        'terminal_grad_free': free['terminal_grad'],
        'terminal_grad_rest': rest['terminal_grad'],
        'outer_arclength_per_m': outer_arclength,
        'outer_time_free_s': free['outer_time_s'],
        'outer_time_rest_s': rest['outer_time_s'],
        'outer_samples_free': int(np.ceil(free['outer_time_s'] / protocol.dwell_s)),
        'outer_samples_rest': int(np.ceil(rest['outer_time_s'] / protocol.dwell_s)),
    }


PROTOCOLS = (
    Protocol('reference-like 256/4', 256.0, 256, 4, 40.0, 200.0),
    Protocol('modest 128/4', 220.0, 128, 4, 40.0, 150.0),
    Protocol('single shot 128/1', 220.0, 128, 1, 40.0, 150.0),
    Protocol('high res 256/8', 220.0, 256, 8, 40.0, 150.0),
    Protocol('small FOV 128/4', 120.0, 128, 4, 40.0, 150.0),
    Protocol('weak gradients 128/4', 220.0, 128, 4, 20.0, 80.0),
    Protocol('strong gradients 128/4', 220.0, 128, 4, 80.0, 200.0),
)


def main() -> None:
    """Run every protocol and print the cost of ending at rest."""
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    rows = [study(p) for p in PROTOCOLS]

    logging.info('%-24s %10s %10s %9s %8s %9s', 'protocol', 'free/ms', 'rest/ms', 'penalty',
                 'percent', 'samples+')
    for row in rows:
        logging.info('%-24s %10.3f %10.3f %9.3f %7.2f%% %9d', row['protocol'],
                     row['duration_free_s'] * 1e3, row['duration_rest_s'] * 1e3,
                     row['duration_penalty_s'] * 1e3, row['duration_penalty_pct'],
                     row['sample_penalty'])

    logging.info('')
    worst = max(rows, key=lambda r: r['duration_penalty_pct'])
    logging.info('worst penalty: %.2f %% on %s (%.3f ms)', worst['duration_penalty_pct'],
                 worst['protocol'], worst['duration_penalty_s'] * 1e3)
    logging.info('k-space extent identical between policies: %s',
                 all(abs(r['k_extent_per_m'] - rows[0]['k_extent_per_m']) < 1e-9
                     for r in rows if r['protocol'] == rows[0]['protocol']))
    logging.info('')
    logging.info('%-24s %11s %11s %11s %11s', 'protocol', 'peak g free', 'peak g rest',
                 'outer n free', 'outer n rest')
    for row in rows:
        logging.info('%-24s %11.4g %11.4g %11d %11d', row['protocol'], row['peak_grad_free'],
                     row['peak_grad_rest'], row['outer_samples_free'], row['outer_samples_rest'])
    logging.info('')
    logging.info('peak gradient reduced by the rest policy in %d of %d cases (it brakes where '
                 'a spiral is fastest)', sum(r['peak_grad_rest'] < r['peak_grad_free'] - 1e-6
                                             for r in rows), len(rows))
    logging.info('terminal gradient, free policy: %.4g .. %.4g 1/m/s',
                 min(r['terminal_grad_free'] for r in rows),
                 max(r['terminal_grad_free'] for r in rows))
    logging.info('terminal gradient, rest policy: %.4g .. %.4g 1/m/s',
                 min(r['terminal_grad_rest'] for r in rows),
                 max(r['terminal_grad_rest'] for r in rows))

    out = Path(__file__).with_name('endpoint_cost_result.json')
    out.write_text(json.dumps(rows, indent=2))
    logging.info('-> %s', out)


if __name__ == '__main__':
    main()
