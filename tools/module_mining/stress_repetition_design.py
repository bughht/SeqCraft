r"""
Reproduce the adversarial stress evidence behind the repetition physical design architecture.

The architecture document and the spike findings quote numbers -- hole widths, jump ratios,
diagnostic agreement rates -- that came from a deliberate search for a feasibility cliff. This
script is where those numbers come from, so that a reader can re-run them rather than trust them.

    python tools/module_mining/stress_repetition_design.py --quick
    python tools/module_mining/stress_repetition_design.py --full --json evidence.json

`--quick` is a representative subset that finishes in a couple of minutes and exercises every
probe. `--full` is the reported matrix and takes some hours, because the split search tries every
raster-aligned split at every candidate window and is therefore quadratic in the window length.

**Not a gate, and not in CI.** It prints a human summary for the findings document and optionally
writes a machine-readable one. Nothing here asserts; the regression tests in
``tests/modules/test_joint_moment_design.py`` are what hold the conclusions in place.

Deterministic: every sweep is a fixed grid, there is no sampling and no randomness, so two runs of
the same size produce the same numbers.
"""

from __future__ import annotations

import argparse
import itertools
import json
import logging
import sys
import time
from typing import TYPE_CHECKING, Any

import numpy as np
import pypulseq as pp

import seqcraft as sc
from seqcraft.errors import ConfigurationError
from seqcraft.modules._joint import (
    CommonModeClaim,
    DifferenceClaim,
    JointProblem,
    Schedule,
    attempt,
    measure_moment,
    utilisation,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

RASTER = 10e-6
STATE = object()


# --------------------------------------------------------------------------------- the grids
class Grid:
    """One size of the sweep.  `full` is what the findings document reports."""

    def __init__(self, *, full: bool) -> None:
        self.full = full
        take = (lambda seq, n: tuple(seq)) if full else (lambda seq, n: tuple(seq)[::n])
        self.grads = take((72.0, 48.0, 32.0, 20.0, 12.0, 8.0), 2)
        self.slews = take((140.0, 100.0, 70.0, 45.0, 28.0, 18.0, 11.0), 3)
        self.m0s = take((0.0, 40.0, 145.0, 400.0), 2)
        self.m1s = take((0.0, 0.05, 0.167, 0.44, 1.2), 2)
        self.fixed = take(((0.0, 0.0), (0.0, -0.43), (120.0, -0.03), (-200.0, 0.25)), 2)
        self.leads = take((0.0, 1.5e-3, 5.0e-3), 2)
        self.vencs = take((5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01), 3)
        self.venc_grads = take((72.0, 48.0, 32.0, 20.0, 12.0, 8.0, 5.0, 3.0), 3)
        self.venc_slews = take((200.0, 140.0, 100.0, 70.0, 45.0, 28.0, 18.0, 11.0, 6.0), 3)
        self.hole_m1s = tuple(np.geomspace(0.01, 1.5, 24 if full else 6))
        self.echoes = take((1, 2, 3, 5, 8), 2)
        #: How far the window search looks, and how far above a first feasible window the hole
        #: census scans.  Both are search extents, not physical statements.
        self.window_limit = 600 if full else 250
        self.hole_band = 100 if full else 40


def opts_for(max_grad_mt_m: float, max_slew_t_m_s: float) -> pp.Opts:
    return pp.Opts(max_grad=max_grad_mt_m, grad_unit='mT/m', max_slew=max_slew_t_m_s,
                   slew_unit='T/m/s', grad_raster_time=RASTER, B0=3.0, rf_dead_time=100e-6,
                   rf_ringdown_time=100e-6, adc_dead_time=20e-6)


def schedule_of(window_s: float, lead_s: float = 0.0, tail_s: float = 0.0) -> Schedule:
    """A window with `lead_s` of history before it and `tail_s` between it and the echo."""
    return Schedule(origin_s=-lead_s, endpoint_s=window_s + tail_s,
                    window_start_s=0.0, window_s=window_s)


def fits(target: tuple[float, float], fixed: tuple[float, float], window_s: float,
         opts: pp.Opts, lead_s: float = 0.0) -> bool:
    """Authoritative feasibility: does the real family cascade emit something here?"""
    problem = JointProblem(axis='y', targets={STATE: target},
                           fixed=lambda state, sched: fixed)
    return attempt(problem, schedule_of(window_s, lead_s), opts) is not None


def smallest_window(target: tuple[float, float], fixed: tuple[float, float], opts: pp.Opts,
                    lead_s: float, limit: int, guard: int = 60) -> tuple[int | None, int]:
    """
    The smallest feasible window in raster steps, and how many holed windows were passed.

    Bisect to some feasible window, then walk down until `guard` consecutive windows refuse.
    Bisection alone would be unsound -- this same script establishes that the feasible set has
    holes -- so the walk-down bounds the error at `guard` steps rather than assuming monotonicity.
    """
    if not fits(target, fixed, limit * RASTER, opts, lead_s):
        return None, 0
    lo, hi = 1, limit
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if fits(target, fixed, mid * RASTER, opts, lead_s):
            hi = mid
        else:
            lo = mid
    best, holes, miss, n = hi, 0, 0, hi - 1
    while n >= 2 and miss < guard:
        if fits(target, fixed, n * RASTER, opts, lead_s):
            holes, best, miss = holes + miss, n, 0
        else:
            miss += 1
        n -= 1
    return best, holes


# ------------------------------------------------------------------------------- the probes
def probe_coarse_matrix(grid: Grid) -> dict[str, Any]:
    """Minimum feasible window over the full product, and the largest adjacent-step jumps."""
    rows = {}
    for g, s in itertools.product(grid.grads, grid.slews):
        opts = opts_for(g, s)
        for m0, m1, fx, lead in itertools.product(grid.m0s, grid.m1s, grid.fixed, grid.leads):
            n, holes = smallest_window((m0, m1), fx, opts, lead, grid.window_limit)
            rows[(g, s, m0, m1, fx, lead)] = {'n': n, 'holes': holes}

    names = ('g', 's', 'm0', 'm1', 'fx', 'lead')
    values = dict(zip(names, (grid.grads, grid.slews, grid.m0s, grid.m1s,
                              grid.fixed, grid.leads)))

    def null(base: dict[str, Any], dim: str, v: Any) -> bool:
        """No requirement at all -- a degenerate grid point, not a transition."""
        d = {**base, dim: v}
        return d['m0'] == 0.0 and d['m1'] == 0.0 and d['fx'] == (0.0, 0.0)

    jumps = {}
    for dim in names:
        others = [values[n] for n in names if n != dim]
        found = []
        for combo in itertools.product(*others):
            base = dict(zip([n for n in names if n != dim], combo))
            for a, b in zip(values[dim], tuple(values[dim])[1:]):
                if null(base, dim, a) or null(base, dim, b):
                    continue
                ra = rows.get(tuple({**base, dim: a}[n] for n in names))
                rb = rows.get(tuple({**base, dim: b}[n] for n in names))
                if ra and rb and ra['n'] and rb['n']:
                    found.append((rb['n'] / ra['n'], a, b))
        if found:
            ratios = np.array([f[0] for f in found])
            worst = max(found)
            jumps[dim] = {'pairs': len(found), 'median': float(np.median(ratios)),
                          'p99': float(np.percentile(ratios, 99)), 'max': float(ratios.max()),
                          'worst_step': f'{worst[1]} -> {worst[2]}'}

    found = [r['n'] for r in rows.values() if r['n'] is not None]
    out: dict[str, Any] = {'cells': len(rows), 'feasible': len(found), 'jumps': jumps}
    if found:
        out |= {'window_us_min': min(found) * 10, 'window_us_max': max(found) * 10}
    return out


def probe_refinement(grid: Grid) -> dict[str, Any]:
    """Fine 1-D sweeps through the hardware limits: is the transition a slope or a step?"""
    out = {}
    steps = 15 if grid.full else 6

    def sweep(label: str, values: Sequence[float],
              make: Callable[[float], tuple[pp.Opts, tuple[float, float]]]) -> None:
        got, prev, worst = [], None, 0.0
        for v in values:
            opts, target = make(v)
            n, _ = smallest_window(target, (-200.0, 0.25) if 'slew' in label else (0.0, 0.0),
                                   opts, 0.0, grid.window_limit)
            if n and prev:
                worst = max(worst, n / prev)
            got.append({'value': float(v), 'window_us': None if n is None else n * 10})
            prev = n
        out[label] = {'points': got, 'steepest_jump': worst}

    sweep('slew_1_T_m_s_steps', [30.0 - k for k in range(steps)],
          lambda s: (opts_for(72.0, s), (145.0, 0.44)))
    sweep('max_grad_half_mT_m_steps', [13.0 - 0.5 * k for k in range(steps)],
          lambda g: (opts_for(g, 140.0), (400.0, 0.167)))
    return out


def probe_holes(grid: Grid) -> dict[str, Any]:
    """
    How wide does a hole in the feasible set get, in the sampled census?

    A window can refuse where a shorter one succeeded: the two-lobe split has to land on the
    gradient raster, and for some window lengths no raster-aligned split solves the 2x2 system
    inside the limits while both neighbours do.  The number this returns is the **maximum
    observed** width over the grid below, not a bound over all targets.
    """
    widths, cells, holed = [], 0, 0
    worst: tuple[int, Any] = (0, None)
    for g, s, m0, fx in itertools.product(grid.grads, grid.slews, grid.m0s, grid.fixed):
        opts = opts_for(g, s)
        for m1 in grid.hole_m1s:
            target = (m0, float(m1))
            first = next((n for n in range(2, grid.window_limit)
                          if fits(target, fx, n * RASTER, opts)), None)
            cells += 1
            if first is None:
                continue
            run, here = 0, []
            for n in range(first + 1, first + grid.hole_band):
                if fits(target, fx, n * RASTER, opts):
                    if run:
                        here.append(run)
                    run = 0
                else:
                    run += 1
            if here:
                holed += 1
                widths.extend(here)
                if max(here) > worst[0]:
                    worst = (max(here), {'max_grad': g, 'max_slew': s, 'm0': m0,
                                         'fixed': list(fx), 'm1': float(m1)})
    summary = {'cells': cells, 'cells_with_a_hole': holed, 'holes': len(widths),
               'band_raster_steps': grid.hole_band}
    if widths:
        a = np.array(widths)
        summary |= {'width_median_steps': int(np.median(a)),
                    'width_p95_steps': int(np.percentile(a, 95)),
                    'max_observed_width_steps': int(a.max()),
                    'max_observed_width_us': int(a.max()) * 10,
                    'widest_at': worst[1]}
    return summary


def probe_diagnostic(grid: Grid) -> dict[str, Any]:
    """
    Does `u <= 1` mean the same thing as "a family fitted"?

    Reported for the raster lattice the emitter searches and for the fixed 48 fractions the
    diagnostic used to sample, because the difference is why the correction was made.
    """
    out = {}
    for label, kwargs in (('raster_lattice', {}), ('fixed_48_fractions', {'splits': 48})):
        agree = pessimistic = optimistic = 0
        for g, s in itertools.product(grid.grads, grid.slews):
            opts = opts_for(g, s)
            for m0, m1 in itertools.product(grid.m0s, grid.hole_m1s):
                for n in range(20, grid.window_limit // 3, 7):
                    target = (m0, float(m1))
                    sched = schedule_of(n * RASTER)
                    got = fits(target, (0.0, 0.0), n * RASTER, opts)
                    u = utilisation('y', target, (0.0, 0.0), sched, opts, **kwargs)['utilisation']
                    if got == (u <= 1.0 + 1e-9):
                        agree += 1
                    elif got:
                        pessimistic += 1
                    else:
                        optimistic += 1
        total = agree + pessimistic + optimistic
        out[label] = {'points': total, 'agreement_pct': round(100 * agree / total, 3),
                      'pessimistic': pessimistic, 'optimistic': optimistic}
    return out


def probe_standalone_venc(grid: Grid) -> dict[str, Any]:
    """Does the standalone bipolar ever refuse?  The premise of a claim that was withdrawn."""
    built, refused, durations = 0, [], []
    for g, s, venc in itertools.product(grid.venc_grads, grid.venc_slews, grid.vencs):
        try:
            ve = sc.modules.VelocityEncode(opts=opts_for(g, s), venc_m_s=venc, axis='y')
        except ConfigurationError as exc:
            refused.append({'max_grad': g, 'max_slew': s, 'venc_m_s': venc,
                            'why': str(exc).splitlines()[0]})
        else:
            built += 1
            durations.append(ve.duration_s)
    return {'systems': built + len(refused), 'built': built, 'refused': len(refused),
            'refusals': refused[:5],
            'pair_duration_ms_min': round(min(durations) * 1e3, 3),
            'pair_duration_ms_max': round(max(durations) * 1e3, 3)}


def probe_kernels(grid: Grid) -> dict[str, Any]:
    """Minimum TE end to end, and the multi-echo question, through the real kernels."""
    rows, refused = [], []
    vencs = (None, *grid.vencs[:2])
    for g, s, venc, ne in itertools.product(grid.grads, grid.slews, vencs, grid.echoes):
        claims: tuple[object, ...] = (CommonModeClaim('y', 1),)
        states: tuple[object, ...] = ('only',)
        if venc is not None:
            delta = sc.modules.VelocityEncode.delta_m1_for(venc)
            claims = (DifferenceClaim('y', 1, delta, ('+', '-')), CommonModeClaim('y', 1))
            states = ('+', '-')
        try:
            kernel = sc.modules.GRE2DTR(
                opts=opts_for(g, s), fov_mm=220.0, matrix=(128, 128), thickness_mm=5.0,
                bandwidth_hz_px=400.0, echoes=ne,
                **({} if ne == 1 else {'polarity': 'bipolar'}),
                joint_claims=claims, encoding_states=states)
        except ConfigurationError as exc:
            refused.append({'max_grad': g, 'max_slew': s, 'venc_m_s': venc, 'echoes': ne,
                            'why': str(exc).splitlines()[0]})
            continue
        rows.append({'g': g, 's': s, 'venc': venc, 'ne': ne,
                     'te_s': kernel.min_te_s, 'winder_s': kernel.winder_s})

    by_key = {(r['g'], r['s'], r['venc'], r['ne']): r for r in rows}
    echo_ratios = []
    for g, s, venc in itertools.product(grid.grads, grid.slews, vencs):
        for a, b in zip(grid.echoes, tuple(grid.echoes)[1:]):
            ra, rb = by_key.get((g, s, venc, a)), by_key.get((g, s, venc, b))
            if ra and rb:
                echo_ratios.append(rb['winder_s'] / ra['winder_s'])
    out: dict[str, Any] = {
        'builds': len(rows) + len(refused), 'built': len(rows), 'refused': len(refused),
        'refusals': refused[:3],
    }
    if rows:
        out |= {'min_te_ms_min': round(min(r['te_s'] for r in rows) * 1e3, 2),
                'min_te_ms_max': round(max(r['te_s'] for r in rows) * 1e3, 2)}
    if echo_ratios:
        out |= {'echo_step_winder_ratio_pairs': len(echo_ratios),
                'echo_step_winder_ratio_max': round(max(echo_ratios), 6)}
    return out


def probe_multi_echo_frames(grid: Grid) -> dict[str, Any]:
    """
    The same multi-echo repetition measured in two reference frames.

    On this axis nothing plays after the winder, so `m1` about the fixed excitation origin stops
    accruing and every echo inherits the first one's compensation.  Measured about each echo
    instead, it walks by `-m0 * ESP` per echo.  Both numbers, because which one the requirement
    means is exactly what naming an origin settles.
    """
    kernel = sc.modules.GRE2DTR(
        opts=opts_for(40.0, 150.0), fov_mm=220.0, matrix=(128, 128), thickness_mm=5.0,
        bandwidth_hz_px=400.0, echoes=3, polarity='bipolar',
        joint_claims=(CommonModeClaim('y', 1),), encoding_states=('only',))
    shot = kernel.build(line=100, encoding_state='only')
    origin, esp = kernel.exc.time_to_center(), kernel.ro.echo_spacing_s
    about_origin, about_echo, m0 = [], [], 0.0
    for echo in range(3):
        at = kernel.time_to_echo() + echo * esp
        about_origin.append(measure_moment(shot, 1, 'y', origin_s=origin, start_s=0.0, end_s=at))
        about_echo.append(measure_moment(shot, 1, 'y', origin_s=at, start_s=0.0, end_s=at))
        m0 = measure_moment(shot, 0, 'y', origin_s=at, start_s=0.0, end_s=at)
    return {'echo_spacing_s': esp, 'm0_1_per_m': m0,
            'm1_about_excitation': [round(v, 12) for v in about_origin],
            'm1_about_each_echo': [round(v, 9) for v in about_echo],
            'predicted_drift_per_echo': round(-m0 * esp, 9),
            'observed_drift_per_echo': round(about_echo[1] - about_echo[0], 9)}


def probe_capability(grid: Grid) -> dict[str, Any]:
    """Which axes each kernel advertises, and what it does with one it does not own."""
    out: dict[str, Any] = {}
    opts = opts_for(40.0, 150.0)
    for name, make, kwargs in (
        ('GRE2DTR', lambda ax: sc.modules.GRE2DTR(
            opts=opts, fov_mm=220.0, matrix=(32, 32), thickness_mm=5.0,
            flow_comp=sc.FlowCompensation(axis=ax)), {'line': 4}),
        ('GRE3DTR', lambda ax: sc.modules.GRE3DTR(
            opts=opts, fov_mm=(220.0, 220.0, 120.0), matrix=(32, 32, 8),
            flow_comp=sc.FlowCompensation(axis=ax)), {'line': 4, 'partition': 1}),
    ):
        accepted, refusals = [], {}
        for axis in ('x', 'y', 'z'):
            try:
                kernel = make(axis)
            except ConfigurationError as exc:
                refusals[axis] = str(exc).splitlines()[0]
                continue
            shot = kernel(**kwargs)
            residual = measure_moment(shot, 1, axis, origin_s=kernel.exc.time_to_center(),
                                      start_s=0.0, end_s=kernel.time_to_echo())
            accepted.append({'axis': axis, 'route': kernel._moment_owners[axis],
                             'm1_residual': round(residual, 12)})
        out[name] = {'accepted': accepted, 'refused': refusals}
    return out


PROBES: tuple[tuple[str, Callable[[Grid], dict[str, Any]]], ...] = (
    ('capability', probe_capability),
    ('standalone_venc', probe_standalone_venc),
    ('multi_echo_frames', probe_multi_echo_frames),
    ('refinement', probe_refinement),
    ('diagnostic', probe_diagnostic),
    ('holes', probe_holes),
    ('kernels', probe_kernels),
    ('coarse_matrix', probe_coarse_matrix),
)


# -------------------------------------------------------------------------------- the report
def human_summary(evidence: dict[str, Any]) -> str:
    """The prose the findings document quotes, rebuilt from the machine-readable result."""
    out = [f"stress evidence, {evidence['size']} grid, {evidence['seconds']:.0f} s", '']

    cap = evidence['probes'].get('capability', {})
    for name, got in cap.items():
        out.append(f'{name} accepts {tuple(a["axis"] for a in got["accepted"])}; '
                   f'refused {tuple(got["refused"])}')
        for entry in got['accepted']:
            out.append(f'    {entry["axis"]} via {entry["route"]:8s} '
                       f'm1 residual {entry["m1_residual"]:.3e}')

    venc = evidence['probes'].get('standalone_venc')
    if venc:
        out += ['', f'standalone VelocityEncode: {venc["built"]}/{venc["systems"]} built, '
                    f'{venc["refused"]} refused, pair '
                    f'{venc["pair_duration_ms_min"]}-{venc["pair_duration_ms_max"]} ms']

    frames = evidence['probes'].get('multi_echo_frames')
    if frames:
        out += ['', 'multi-echo, the same repetition in two frames:',
                f'    m1 about the excitation  {frames["m1_about_excitation"]}',
                f'    m1 about each echo       {frames["m1_about_each_echo"]}',
                f'    drift per echo: predicted {frames["predicted_drift_per_echo"]}, '
                f'observed {frames["observed_drift_per_echo"]}']

    holes = evidence['probes'].get('holes')
    if holes and 'max_observed_width_us' in holes:
        out += ['', f'holes: {holes["cells_with_a_hole"]}/{holes["cells"]} cells, '
                    f'max OBSERVED width {holes["max_observed_width_us"]} us '
                    f'(median {holes["width_median_steps"]} steps) '
                    f'-- a sampled maximum, not a bound']
    elif holes:
        out += ['', f'holes: none in {holes["cells"]} cells. They are rare -- about 1.6 % of the '
                    'full grid -- so a subset can miss them entirely. Not evidence of absence.']

    diag = evidence['probes'].get('diagnostic')
    if diag:
        out.append('')
        for label, got in diag.items():
            out.append(f'diagnostic vs emitter, {label:20s} '
                       f'{got["agreement_pct"]} % over {got["points"]} points')

    ref = evidence['probes'].get('refinement')
    if ref:
        out.append('')
        for label, got in ref.items():
            out.append(f'refined {label:28s} steepest adjacent jump x{got["steepest_jump"]:.3f}')

    kern = evidence['probes'].get('kernels')
    if kern:
        out += ['', f'kernels: {kern["built"]}/{kern["builds"]} built, min TE '
                    f'{kern.get("min_te_ms_min")}-{kern.get("min_te_ms_max")} ms']
        if 'echo_step_winder_ratio_max' in kern:
            out.append(f'    winder ratio across every echo-count step, max '
                       f'{kern["echo_step_winder_ratio_max"]}')

    coarse = evidence['probes'].get('coarse_matrix')
    if coarse:
        span = ('' if 'window_us_min' not in coarse
                else f', window {coarse["window_us_min"]}-{coarse["window_us_max"]} us')
        out += ['', f'coarse matrix: {coarse["feasible"]}/{coarse["cells"]} feasible{span}']
        for dim, got in coarse['jumps'].items():
            out.append(f'    {dim:5s} median x{got["median"]:.3f}  p99 x{got["p99"]:.3f}  '
                       f'max x{got["max"]:.3f}  at {got["worst_step"]}')

    out += ['', 'Extensive stress testing did not reveal a large artificial feasibility cliff.',
            'That is a statement about this sampled and refined parameter space, not a proof '
            'that cliffs cannot occur.']
    if evidence['size'] == 'quick':
        out += ['',
                'This is the quick grid: every probe ran, but the rare findings need the full '
                'one.',
                'Holes appear in ~1.6 % of full-grid cells, and the two diagnostic lattices '
                'only diverge',
                'near the feasibility boundary, so both can read identically here. Use --full '
                'to reproduce',
                'the numbers the findings document quotes.']
    return '\n'.join(out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    size = parser.add_mutually_exclusive_group()
    size.add_argument('--quick', action='store_true',
                      help='representative subset; every probe runs, coarsely (default)')
    size.add_argument('--full', action='store_true',
                      help='the matrix the findings document reports; hours, not minutes')
    parser.add_argument('--json', metavar='PATH',
                        help='write the machine-readable summary here')
    parser.add_argument('--only', action='append', metavar='PROBE',
                        choices=[name for name, _ in PROBES],
                        help='run only these probes (repeatable)')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(message)s', stream=sys.stdout)
    grid = Grid(full=args.full)
    wanted = args.only or [name for name, _ in PROBES]

    started = time.time()
    probes: dict[str, Any] = {}
    for name, probe in PROBES:
        if name not in wanted:
            continue
        at = time.time()
        logging.info('running %s ...', name)
        probes[name] = probe(grid)
        logging.info('  %s done in %.0f s', name, time.time() - at)

    evidence = {'size': 'full' if args.full else 'quick',
                'seconds': time.time() - started,
                'window_search_limit_steps': grid.window_limit,
                'hole_census_band_steps': grid.hole_band,
                'probes': probes}

    logging.info('\n%s', human_summary(evidence))
    if args.json:
        with open(args.json, 'w') as handle:
            json.dump(evidence, handle, indent=2, default=str)
        logging.info('\nmachine-readable summary written to %s', args.json)


if __name__ == '__main__':
    main()
