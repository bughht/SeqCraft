# The PyPulseq dependency

This page records which PyPulseq SeqCraft is pinned to and why. It covers what SeqCraft relies on
that official PyPulseq does not yet reproduce, one known defect in the pinned package and how
SeqCraft avoids it, and what to re-check before the pin moves. Update it whenever the pin changes.

## The pin

| | |
|---|---|
| package | `pypulseq @ https://github.com/m-a-x-i-m-z/pypulseq-matlab-like/archive/22ef2db1f71ff38c8ce355c61913cfd8fceaac3b.zip` |
| pinned in | `pyproject.toml`, `ci/constraints.txt`, and the install hint in `src/seqcraft/_compat.py` |
| what it is | a compatibility fork of PyPulseq that follows MATLAB Pulseq 1.5.1 behaviour |
| reports | `pypulseq.__version__ == '1.5.1'`; distribution metadata `1.5.0.post1` |

**Why the source URL is the contract.** The distribution metadata is the same `1.5.0.post1` as
the incompatible PyPI wheel, so no version specifier can tell them apart. `seqcraft._compat.require()`
therefore probes capabilities at import instead.

**History.**
- SeqCraft has pinned this commit since 2026-08-14 (`d505ed6`), after a fresh CI runner installed
  the PyPI wheel, which rejects `Opts(max_b1=...)`.
- `22ef2db` is the fork before it renamed its package to `pypulseq_matlab_like`. Its later commits
  are a separate project, not an upgrade path.

**One PyPulseq in the runtime.** SeqCraft builds every event and `Sequence` with the pinned
package. It does not install official PyPulseq beside it to borrow individual functions, because
events and sequences from two divergent implementations would then mix. SeqCraft also does not
maintain a fork of its own to repair individual methods. A patched copy was built and tested on
2026-11-06 for the trajectory defect below, and set aside in favour of the named function that
avoids it.

## Why official PyPulseq is not yet the runtime dependency

Measured against official [`pulseq/pypulseq`](https://github.com/pulseq/pypulseq) at `a8cb2a8`
(2026-10-29; metadata `1.5.0.post1`; Python ≥ 3.10), on 2026-11-05:

| behaviour | where SeqCraft relies on it | official at `a8cb2a8` |
|---|---|---|
| `Opts.max_b1` (with `b1_unit`, default 20 µT) | the transmit-limit refusal (`modules/_support.require_usable_max_b1`, the compiler's RF check); user code passes `pp.Opts(max_b1=...)` | absent: the keyword raises `TypeError` |
| `Opts.rf_samples_limit` | the compiler's RF sample-count check | absent |
| `Opts.max_freq_offset` | the import-time probe only | absent; the probe could drop it on a migration |
| `make_slr_pulse` | `pulse='slr'` in `Excitation`, `Refocusing`, `IRPrep`, `SaturationPrep` | absent. `sigpy_n_seq` exists without `dwell` and `return_delay` |
| `make_trapezoid(area=, duration=)` gives MATLAB Pulseq's shape: the lowest amplitude for the duration, with maximum-slew ramps | 17 call sites, including the moment-designed lobes of `design/joint.py` | takes the shortest trapezoid's ramps and stretches the flat top, so different waveforms, TE/TR minima, PNS and vector-norm load. It also refuses an exactly minimal duration through an exact float comparison |
| `make_extended_trapezoid(max_grad=inf, max_slew=inf, skip_check=True)` builds | the compiler builds merged gradients this way and reports an over-limit merge itself, as a structured error | `inf` means the system limit, and slew is always checked: a bare `ValueError` |
| `calculate_kspacePP()` returns `(k_adc, t_adc, k, t_k, t_excitation, t_refocusing, ...)` | `sc.sequence_kspace`, including the dense `k` and `t_k` | a deprecated alias of `calculate_kspace()`, with its other 5-value order. The dense trajectory's time base is not returned |
| `get_block()` gives labels and triggers as lists (`block.label`, `block.trig`) | tests and `tools/capture_compiler_baseline.py` (`src/` does not read them) | dicts, and the trigger is `block.trigger`. Both are still preserved through write/read |
| RF `use='other'` is preserved | a test control value only | becomes `'undefined'`, which counts as an excitation |
| the writer emits Pulseq **1.5.1** | the scanner interpreter | 1.5.0 (same sections) |
| shapes are deduplicated at six significant digits | the compiler's registration tolerance accounts for it | exact |

**Not a blocker, whichever package is used.** `RosetteReadout` passes `make_adc` a negative delay,
relying on the fork raising it to `adc_dead_time`. Official refuses it, so clamping it at zero in
SeqCraft is correct under both.

## One nominal trajectory

> In an ideal sequence there is one nominal `k(t_adc)`. A SeqCraft module's reported trajectory,
> the exact trajectory of the compiled `Sequence`, the trajectory of the written `.seq` file and an
> ideal simulator should agree within numerical and file-rounding tolerance.

SeqCraft's supported routes to it are:

| you have | call |
|---|---|
| a tree | `sc.kspace(tree, opts)` |
| a compiled `Sequence`, or one read back from a `.seq` file | `sc.sequence_kspace(seq)` |

Both return the same named fields: `k_adc` and `t_adc` at the true ADC sample times, the dense
`k` and `t_k`, and `t_excitation` and `t_refocusing`. `sc.kspace(tree, opts)` is
`sc.sequence_kspace(sc.compile(tree, opts))`. In the current pin, `sc.sequence_kspace` uses
`calculate_kspacePP()`, which is an exact piecewise-polynomial integration.

## Known defect in the pin: `Sequence.calculate_kspace()`

**`seq.calculate_kspace()` is not a supported route to the nominal trajectory.** SeqCraft returns a
bare `pypulseq.Sequence`, so the method is within reach, but examples and documentation do not use
it, and a reconstruction should not.

The pinned `calculate_kspace()` ports MATLAB Pulseq's legacy raster routine, which MATLAB keeps
only as an unused `calculateKspaceUnfunc`. It samples every gradient at raster edges, accumulates a
right-hand rectangle sum restarted at `round(t_exc / D)`, and interpolates linearly to the ADC
times. Its error at `t = jD + τ` is:
- `D/2 · G(t)`;
- plus `(G[j+1] − G[j]) · τ(D − τ) / (2D)` between raster edges;
- plus a constant per excitation from any gradient running through the RF centre.

That comes to up to 0.9 1/m on a Wave readout and 4 1/m on the axis of a 5 mm slice. A fixed
`trajectory_delay = −D/2` cancels the first term only when every sample is on a raster edge and
nothing runs through the RF centre, so it is not a repair.

`examples/bssfp_2d` once printed this as "|kx| at the echo = 0.224 of a k-space step" and explained
it as half a dwell of sampling. Read through `sc.sequence_kspace`, `kx` at the echo is 8e-11 1/m.
The half dwell is a timing quantity: `symmetry_residual_s`, how far the echo sits from TR/2.

**Where it is tested.** `tests/analysis/test_sequence_kspace.py` checks `sc.sequence_kspace`
against an independent exact integration of the block events. The cases are ramps, trapezoids
after a slice-selective excitation, arbitrary gradients, Wave readouts with odd and even cycle
counts, and a spin echo. Each is sampled on raster edges, half a raster off them and at an
incommensurate dwell, in memory and read back. One test there asserts that the pinned
`calculate_kspace()` is still biased. When the pin moves to a package whose `calculate_kspace()` is
exact, that test fails: delete it, and update this section.

## Moving to official PyPulseq

Official `calculate_kspace()` at `a8cb2a8` is already an exact piecewise-polynomial integration,
and official `calculate_kspacePP()` is a deprecated wrapper around it. Once official PyPulseq
becomes viable as the runtime dependency, `sc.sequence_kspace` switches its backend to official
`calculate_kspace()`, and its named contract does not change. The one gap is the dense trajectory's
time base, which official does not return (checklist item 6).

**Not as of the audit above.** Moving changes more than the trajectory API. The emitted gradient
shapes would change across nearly every sequence, and the `Opts` limits and SLR would need a new
home.

Whenever official `pulseq/pypulseq` or the compatibility fork changes materially, re-check each item
below on the candidate, and record the result here:

1. **RF and system-limit fields.** `Opts.max_b1`, `rf_samples_limit`, `max_freq_offset`,
   `b1_unit`, and their defaults.
2. **SLR.** `make_slr_pulse`, or a replacement producing the same waveforms for the four modules.
3. **Trapezoid realisation.** The shape `make_trapezoid(area=, duration=)` emits, its handling of an
   exactly minimal duration, and its rephasers inside `make_sinc_pulse` and similar.
4. **Arbitrary and extended gradients.** `make_extended_trapezoid` with infinite construction
   limits; first/last and raster-centre semantics; `restore_additional_shape_samples` on read.
5. **Exact `k(t_adc)`.** The candidate backend for `sc.sequence_kspace` against an independent
   integration, with a gradient through the RF centre and samples off the raster
   (`tests/analysis/test_sequence_kspace.py`).
6. **The trajectory API.** The return order of `calculate_kspace()` and `calculate_kspacePP()`, and
   whether the dense trajectory's time base is available.
7. **Labels, triggers and extensions.** Their `get_block` representation, round trip, and
   `evaluate_labels(evolution='adc')`; rotation extensions.
8. **Write/read.** The file-format version, `[DEFINITIONS]`, signatures, RF `use` (all seven
   values), and deduplication.
9. **Emitted physics on representative sequences.** Write every example's `.seq` under both and
   compare block by block: gradient knots, RF, ADC, labels, `k(t_adc)`, TE/TR minima, PNS.
10. **Full checks.** The test suite, doctests, the API reference, and the notebook smoke run
    (`docs/testing.md`), before and after any compatibility adjustment.

A detailed record of the 2026-11-05 comparison and the 2026-11-06 trajectory investigation, with
probes, per-test classifications and the patched copy, is kept outside this repository, in the
SeqCraft workspace at `docs/plans/pypulseq_dependency/README.md`.
