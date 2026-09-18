"""
A comparator stack, not a hash.

Each layer answers a different question, and each returns **measurements** so that a difference
can be read rather than merely detected:

``L0``  metadata sanity -- duration, event and sample counts, definitions.
``L2``  the physical waveform on its own exact knots, integrated where an area is what matters.
``L3``  acquisition semantics -- k at every ADC sample, and which sample is the echo.
``L4``  TSE physics -- refocusing spacing, echo placement, and the area balance about each pulse.

``L1`` (the exact SeqCraft event digest) is not here: it applies only between two SeqCraft trees
and ``tests/`` already owns that comparison.

Two rules this module exists to enforce:

* **Block count is never an equivalence criterion.**  It is reported at L0 for diagnostics and
  used by nothing.
* **Nothing is assumed about which sample carries the echo.**  The echo is located as the
  ``kx`` zero *crossing*, interpolated between samples, because an implementation is not obliged
  to put a sample there -- and the first thing measuring these references found is that the
  PyPulseq TSE demo does not: with 64 samples its k-space centre falls exactly between samples 31
  and 32.  Ranking by ``argmin|kx|`` there picks between two values half a k-step either side of
  zero, and which one wins is floating-point noise.  Measuring the crossing instead makes ``ky``
  and ``kz`` at the echo comparable across implementations that disagree about sampling, and
  reports the sampling difference separately as ``echo_offset_dwells``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    import pypulseq as pp

    from .reference import ReferenceSequence

#: Channel order used everywhere here, matching ``calculate_kspacePP``'s rows.
AXES = ('x', 'y', 'z')


# --------------------------------------------------------------------------- L0 / structure
def _readout_sample_counts(seq: pp.Sequence) -> list[int]:
    """Samples per ADC, in play order, read from the block table rather than from ``get_block``."""
    counts = []
    for block_id in sorted(seq.block_events):
        adc_id = int(seq.block_events[block_id][5])
        if adc_id:
            counts.append(int(round(float(seq.adc_library.data[adc_id][0]))))
    return counts


def _rf_count(seq: pp.Sequence) -> int:
    return sum(1 for block_id in seq.block_events if int(seq.block_events[block_id][1]))


# ------------------------------------------------------------------------------ L2 / areas
def _area(knots: np.ndarray, t0: float, t1: float) -> float:
    """
    Integrate one piecewise-linear gradient channel between two times, in 1/m.

    `knots` is a ``(2, n)`` array of times and amplitudes in Hz/m, as
    ``Sequence.waveforms_and_times`` returns it.  The endpoints are interpolated, so the result
    does not depend on where the knots happen to fall.

    Caveat worth knowing when a number looks wrong: intervals where a channel plays nothing carry
    no knots, and this interpolates straight across them.  That is correct while events begin and
    end at zero, which is what both the compiler and the references produce, and it is not correct
    for a waveform deliberately left at a non-zero amplitude across a gap.
    """
    t, g = np.asarray(knots[0], dtype=float), np.asarray(knots[1], dtype=float)
    if t.size < 2 or t1 <= t0:
        return 0.0
    inner = (t > t0) & (t < t1)
    times = np.concatenate(([t0], t[inner], [t1]))
    amps = np.concatenate((
        [float(np.interp(t0, t, g))], g[inner], [float(np.interp(t1, t, g))],
    ))
    return float(np.trapezoid(amps, times)) if hasattr(np, 'trapezoid') else float(
        np.trapz(amps, times))  # noqa: NPY201


def _echo_instant(times: np.ndarray, kx: np.ndarray) -> tuple[float, int, bool]:
    """
    Locate ``kx = 0`` within one readout: its time, the nearest sample, and whether it was found.

    Interpolated across the sign change rather than taken from the nearest sample, for the reason
    in the module docstring.  When ``kx`` never changes sign -- a readout that does not cross the
    centre of k-space at all -- the nearest sample is returned and ``crossed`` is ``False``, so a
    caller can tell a measurement from a fallback instead of being handed a plausible number.
    """
    sign_change = np.nonzero(np.diff(np.signbit(kx)))[0]
    if kx.size < 2 or sign_change.size == 0:
        nearest = int(np.argmin(np.abs(kx)))
        return float(times[nearest]), nearest, False
    i = int(sign_change[0])
    span = kx[i + 1] - kx[i]
    t_echo = float(times[i] - kx[i] * (times[i + 1] - times[i]) / span) if span else float(times[i])
    nearest = int(np.argmin(np.abs(times - t_echo)))
    return t_echo, nearest, True

# ------------------------------------------------------------------------------- measurement
def measure(ref: ReferenceSequence) -> dict[str, Any]:
    """
    Every deterministic measurement this comparator can take of one executed reference.

    Returns a nested dict rather than an object so that a report can be written out, diffed and
    read by a human without importing anything.
    """
    seq = ref.sequence
    k_adc, t_adc, _k, _t_k, t_exc, t_refoc = seq.calculate_kspacePP()[:6]
    k_adc, t_adc = np.asarray(k_adc), np.asarray(t_adc)
    t_exc, t_refoc = np.sort(np.asarray(t_exc)), np.sort(np.asarray(t_refoc))
    waveforms = seq.waveforms_and_times()[0]

    counts = _readout_sample_counts(seq)
    if sum(counts) != t_adc.size:
        msg = (f'{ref.name}: {sum(counts)} samples in the block table against {t_adc.size} ADC '
               f'sample times; the readout grouping would be wrong.')
        raise ValueError(msg)

    readouts, start = [], 0
    for index, n in enumerate(counts):
        span = slice(start, start + n)
        times, kx = t_adc[span], k_adc[0, span]
        dwell = float(np.median(np.diff(times))) if n > 1 else float('nan')
        t_echo, nearest, crossed = _echo_instant(times, kx)
        readouts.append({
            'index': index,
            'n_samples': n,
            'shot': int(np.searchsorted(t_exc, times[0], side='right') - 1),
            'sample_of_echo': nearest,
            'echo_is_a_crossing': crossed,
            # How far the nearest ADC sample sits from k = 0, in dwells.  0 means a sample carries
            # DC; 0.5 means the centre of k-space falls between two samples and the acquisition
            # has no DC sample at all.  A design difference, not an error -- and invisible to
            # every k-space check that only looks at the nearest sample.
            'echo_offset_dwells': float(abs(times[nearest] - t_echo) / dwell) if n > 1 else 0.0,
            't_echo_s': float(t_echo),
            't_first_sample_s': float(times[0]),
            'dwell_s': dwell,
            'k_at_echo_per_m': [float(np.interp(t_echo, times, k_adc[axis, span]))
                                for axis in range(3)],
            'k_at_nearest_sample_per_m': [float(k_adc[axis, span][nearest]) for axis in range(3)],
            'k_extent_per_m': [float(np.ptp(k_adc[axis, span])) for axis in range(3)],
        })
        start += n

    shots = []
    for shot, t0 in enumerate(t_exc):
        t1 = t_exc[shot + 1] if shot + 1 < t_exc.size else np.inf
        centres = t_refoc[(t_refoc >= t0) & (t_refoc < t1)]
        mine = [r for r in readouts if r['shot'] == shot]
        spacing = np.diff(centres)
        shots.append({
            'shot': shot,
            't_excitation_s': float(t0),
            'n_refocusing': int(centres.size),
            'n_readouts': len(mine),
            'refocusing_s': [float(v) for v in centres],
            'spacing_s': [float(v) for v in spacing],
            'spacing_spread_s': float(np.ptp(spacing)) if spacing.size else 0.0,
            'echo_spacing_s': float(np.median(spacing)) if spacing.size else float('nan'),
            'first_interval_s': float(centres[0] - t0) if centres.size else float('nan'),
        })

    return {
        'name': ref.name,
        'provenance': dict(ref.provenance),
        'parameters': dict(ref.parameters),
        'semantic_claimed': dict(ref.semantic),
        'l0': {
            'duration_s': float(seq.duration()[0]),
            'n_blocks': len(seq.block_events),          # diagnostics only -- never an equivalence
            'n_rf': _rf_count(seq),
            'n_adc': len(counts),
            'n_adc_samples': int(sum(counts)),
            'n_excitation': int(t_exc.size),
            'n_refocusing': int(t_refoc.size),
            'definitions': {k: v for k, v in (ref.definitions or {}).items()},
        },
        'readouts': readouts,
        'shots': shots,
        '_waveforms': waveforms,
    }


# ------------------------------------------------------------------------- invariant checks
def check_invariants(m: dict[str, Any], *, expected_ky_per_m: list[float] | None = None,
                     k_tol_per_m: float = 1e-3, spacing_tol_s: float = 1e-9,
                     midpoint_tol_s: float = 10e-6) -> list[dict[str, Any]]:
    """
    Measure the TSE/FSE invariant table against one reference's measurements.

    Returns one row per invariant, each carrying the **worst observed value** and where it
    occurred, so a failure names an axis, a shot and an echo rather than returning ``False``.
    Tolerances default to the table's (``docs/plans/module-mining/tse/invariant_table.md``).
    """
    rows: list[dict[str, Any]] = []
    readouts, shots = m['readouts'], m['shots']

    def worst(values: list[tuple[float, dict]]) -> tuple[float, dict]:
        return max(values, key=lambda pair: pair[0]) if values else (0.0, {})

    # I1 is not "kx = 0 at the echo", which the interpolation would make true by construction.
    # What is worth asserting is that every readout crosses the centre of k-space at all; where
    # the nearest sample falls relative to that crossing is reported, not judged.
    missing = [r['index'] for r in readouts if not r['echo_is_a_crossing']]
    offsets = [r['echo_offset_dwells'] for r in readouts]
    rows.append({'invariant': 'I1', 'what': 'every readout crosses kx = 0',
                 'worst': float(len(missing)), 'tolerance': 0.0, 'unit': 'readouts',
                 'where': {'readouts_without_a_crossing': missing[:8]},
                 'pass': not missing})
    rows.append({'invariant': 'I1-dc', 'what': 'distance from the nearest ADC sample to kx = 0',
                 'worst': float(max(offsets)) if offsets else 0.0, 'tolerance': None,
                 'unit': 'dwells', 'where': {'note': '0 = a sample carries DC; 0.5 = none does'},
                 'pass': None})

    value, where = worst([
        (abs(r['k_at_echo_per_m'][2]), {'shot': r['shot'], 'readout': r['index']})
        for r in readouts])
    rows.append({'invariant': 'I3', 'what': 'kz = 0 at the echo', 'worst': value,
                 'tolerance': k_tol_per_m, 'unit': '1/m', 'where': where,
                 'pass': value <= k_tol_per_m})

    if expected_ky_per_m is not None:
        if len(expected_ky_per_m) != len(readouts):
            msg = f'{len(expected_ky_per_m)} expected ky against {len(readouts)} readouts'
            raise ValueError(msg)
        value, where = worst([
            (abs(r['k_at_echo_per_m'][1] - want), {'shot': r['shot'], 'readout': r['index'],
                                                   'want_per_m': want,
                                                   'got_per_m': r['k_at_echo_per_m'][1]})
            for r, want in zip(readouts, expected_ky_per_m)])
        rows.append({'invariant': 'I2', 'what': 'ky = the requested line at the echo sample',
                     'worst': value, 'tolerance': k_tol_per_m, 'unit': '1/m', 'where': where,
                     'pass': value <= k_tol_per_m})

    value, where = worst([(s['spacing_spread_s'], {'shot': s['shot']}) for s in shots])
    rows.append({'invariant': 'I5', 'what': 'refocusing centres uniformly spaced within a shot',
                 'worst': value, 'tolerance': spacing_tol_s, 'unit': 's', 'where': where,
                 'pass': value <= spacing_tol_s})

    offsets, first = [], []
    for s in shots:
        centres = np.asarray(s['refocusing_s'])
        mine = sorted((r for r in readouts if r['shot'] == s['shot']), key=lambda r: r['t_echo_s'])
        for n in range(min(len(mine), centres.size) - 1):
            midpoint = (centres[n] + centres[n + 1]) / 2
            offsets.append((abs(mine[n]['t_echo_s'] - midpoint),
                            {'shot': s['shot'], 'echo': n, 'offset_s': mine[n]['t_echo_s'] - midpoint}))
        if centres.size > 1:
            first.append((abs(s['first_interval_s'] - s['echo_spacing_s'] / 2), {'shot': s['shot']}))

    value, where = worst(offsets)
    rows.append({'invariant': 'I6', 'what': 'echo at the midpoint between its two refocusing centres',
                 'worst': value, 'tolerance': midpoint_tol_s, 'unit': 's', 'where': where,
                 'pass': value <= midpoint_tol_s})
    value, where = worst(first)
    rows.append({'invariant': 'I7', 'what': 'first interval = half the echo spacing',
                 'worst': value, 'tolerance': midpoint_tol_s, 'unit': 's', 'where': where,
                 'pass': value <= midpoint_tol_s})

    rows.extend(_check_area_balance(m, k_tol_per_m))
    return rows


def _check_area_balance(m: dict[str, Any], k_tol_per_m: float) -> list[dict[str, Any]]:
    """
    I4, measured independently of I1-I3, and **not the same statement on every axis**.

    On the readout and selection axes the requirement is equal area either side of the refocusing
    centre: whatever is wound before the pulse is unwound after it, so the echo returns to the
    same k.  Running the same check on the phase axis is a category error, and measuring it was
    how that became obvious -- both references "fail" it by ~230 1/m, because a phase encode is
    *meant* to leave k somewhere different from where it started.  That is the whole job of an
    encode.

    What the phase axis must satisfy instead is ``ky = 0`` **at** each refocusing centre: the
    encode is rewound before the pulse and reapplied after it.  Both references build it that way
    -- blip after the 180, rewinder before the next -- and the SeqCraft notebook says why, which
    is that putting the blip before the pulse makes every later pulse flip the encode.

    Integrated from the emitted waveform so a failure names an axis and an echo, which ``k`` at
    the echo can only report as a total.
    """
    waveforms = m['_waveforms']
    rows = []
    for axis_index, axis in enumerate(AXES):
        balance, at_centre = (0.0, {}), (0.0, {})
        for s in m['shots']:
            centres = np.asarray(s['refocusing_s'])
            mine = sorted((r for r in m['readouts'] if r['shot'] == s['shot']),
                          key=lambda r: r['t_echo_s'])
            for n in range(min(len(mine) - 1, centres.size - 1)):
                centre = centres[n + 1]
                before = _area(waveforms[axis_index], mine[n]['t_echo_s'], centre)
                after = _area(waveforms[axis_index], centre, mine[n + 1]['t_echo_s'])
                where = {'shot': s['shot'], 'interval': n, 'axis': axis,
                         'area_before_per_m': before, 'area_after_per_m': after}
                if abs(before - after) > balance[0]:
                    balance = (abs(before - after), where)
                k_centre = mine[n]['k_at_echo_per_m'][axis_index] + before
                if abs(k_centre) > at_centre[0]:
                    at_centre = (abs(k_centre), where | {'k_at_centre_per_m': k_centre})

        if axis == 'y':
            rows.append({'invariant': 'I4y', 'what': 'ky = 0 at every refocusing centre',
                         'worst': at_centre[0], 'tolerance': k_tol_per_m, 'unit': '1/m',
                         'where': at_centre[1], 'pass': at_centre[0] <= k_tol_per_m})
        else:
            rows.append({'invariant': f'I4{axis}',
                         'what': f'equal {axis} area either side of each refocusing centre',
                         'worst': balance[0], 'tolerance': k_tol_per_m, 'unit': '1/m',
                         'where': balance[1], 'pass': balance[0] <= k_tol_per_m})
    return rows


# ------------------------------------------------------------------------------- comparison
def compare(a: dict[str, Any], b: dict[str, Any], *, k_tol_per_m: float = 1e-3,
            time_tol_s: float = 1e-9) -> dict[str, Any]:
    """
    Compare two references measured at the **same protocol**, and return a difference report.

    Refuses rather than guesses when the two do not have the same number of readouts or samples:
    a comparison across different protocols is not a physical statement, and silently truncating
    to the shorter one is how a comparator learns to pass.
    """
    report: dict[str, Any] = {'a': a['name'], 'b': b['name'], 'result': 'pass',
                              'l0': {}, 'l3': {}, 'first_failure': None}

    for key in ('duration_s', 'n_adc', 'n_adc_samples', 'n_rf', 'n_excitation', 'n_refocusing'):
        va, vb = a['l0'][key], b['l0'][key]
        report['l0'][key] = {'a': va, 'b': vb, 'delta': (vb - va)}
    report['l0']['n_blocks'] = {'a': a['l0']['n_blocks'], 'b': b['l0']['n_blocks'],
                                'note': 'diagnostic only; not an equivalence criterion'}

    if a['l0']['n_adc'] != b['l0']['n_adc'] or a['l0']['n_adc_samples'] != b['l0']['n_adc_samples']:
        report['result'] = 'incomparable'
        report['first_failure'] = {'check': 'readout inventory',
                                   'detail': 'different ADC count or sample count'}
        return report

    t_err = max(abs(ra['t_echo_s'] - a['shots'][ra['shot']]['t_excitation_s']
                    - (rb['t_echo_s'] - b['shots'][rb['shot']]['t_excitation_s']))
                for ra, rb in zip(a['readouts'], b['readouts']))
    # x and z are physics: both implementations must refocus them at the echo.  y is the
    # phase-encode **table**, which is caller policy and legitimately different between two
    # implementations of the same sequence -- so it is reported and never judged here.  Whether
    # each covers k-space once is `check_ky_coverage`'s question, not this one.
    k_err = [max(abs(ra['k_at_echo_per_m'][axis] - rb['k_at_echo_per_m'][axis])
                 for ra, rb in zip(a['readouts'], b['readouts'])) for axis in range(3)]
    sample_err = max(abs(ra['sample_of_echo'] - rb['sample_of_echo'])
                     for ra, rb in zip(a['readouts'], b['readouts']))

    physics_err = max(k_err[0], k_err[2])
    report['l3'] = {
        'te_from_excitation_max_error_s': t_err,
        'k_at_echo_max_error_per_m': dict(zip(AXES, k_err)),
        'ky_note': 'a table difference, reported and not judged; see check_ky_coverage',
        'echo_sample_index_max_delta': int(sample_err),
        'echo_offset_dwells': {'a': max(r['echo_offset_dwells'] for r in a['readouts']),
                               'b': max(r['echo_offset_dwells'] for r in b['readouts'])},
    }
    if physics_err > k_tol_per_m:
        report['result'] = 'fail'
        axis = 'x' if k_err[0] >= k_err[2] else 'z'
        report['first_failure'] = {'check': 'k at echo', 'axis': axis, 'error_per_m': physics_err}
    elif t_err > time_tol_s:
        report['result'] = 'fail'
        report['first_failure'] = {'check': 'echo time from excitation', 'error_s': t_err}
    return report


def check_ky_coverage(m: dict[str, Any], *, dky_per_m: float, n_lines: int,
                      tol_per_m: float = 1e-3) -> dict[str, Any]:
    """
    Does the acquired set of ``ky`` land on the expected lattice, once each?

    The implementation-independent form of I2, for a reference that never says which line each
    echo carries.  It cannot catch a permutation of the table -- which is why it does not replace
    passing ``expected_ky_per_m`` where the table is known -- but it does catch a missing line, a
    duplicated line, a wrong ``dk``, and a mirrored axis.
    """
    got = np.asarray([r['k_at_echo_per_m'][1] for r in m['readouts']])
    index = np.round(got / dky_per_m).astype(int)
    lattice_error = float(np.max(np.abs(got - index * dky_per_m))) if got.size else 0.0
    unique, counts = np.unique(index, return_counts=True)
    expected = set(range(-(n_lines // 2), n_lines - n_lines // 2))
    return {
        'invariant': 'I2-coverage',
        'what': 'acquired ky lands on the expected lattice, each line once',
        'lattice_max_error_per_m': lattice_error,
        'tolerance': tol_per_m,
        'n_acquired': int(got.size),
        'n_distinct': int(unique.size),
        'duplicated': [int(v) for v in unique[counts > 1]],
        'missing': sorted(expected - set(int(v) for v in unique)),
        'pass': lattice_error <= tol_per_m and bool(unique.size == got.size),
    }
