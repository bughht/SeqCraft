# Phase 0 baseline report

> **Historical.** This records the state at the time it was written.  The compiler was
> subsequently changed by [ADR-004](../adr/004-compile-returns-a-sequence.md): `compile` returns a
> bare `pypulseq.Sequence`, every legality failure raises, and `CompiledSequence`, `Report` and
> `Issue` no longer exist.  For the current shape see [`../compiler.md`](../compiler.md) and
> [`../architecture.md`](../architecture.md).
>
> Final disposition: the original artifact was re-captured after ADR-003/004 changed the public
> architecture and raw-PyPulseq fixtures replaced module-built recipes. It remains a regression
> baseline for current emitted structure, not proof across the original pre-refactor boundary.

## Baseline identity

- Source commit: `afe576bafe17e9cce3edfc768b343ad958994560`
- Remote CI: [run #7](https://github.com/bughht/SeqCraft/actions/runs/31848625433), seven jobs passed
- Local platform: macOS 15.7.1, arm64
- Python 3.12.13; SeqCraft 0.3.0; PyPulseq 1.5.1; NumPy 2.4.6
- pytest 8.3.3; Ruff 0.6.9; mypy 1.11.2
- Dependency source: `ci/constraints.txt`, including the immutable PyPulseq compatibility commit

The full machine-readable record is `tests/baselines/compiler_phase0.json`. Generate a reviewed
replacement with:

```bash
python tools/capture_compiler_baseline.py --iterations 3
```

Regeneration is not an automatic fix for a failing snapshot. Review the field-level difference and
update the artifact only for an approved behavior change or an intentional baseline recapture on a
documented runner.

## Representative sequence observations

| Recipe | Placed events | Blocks | Input gradients | Gradient splits | `.seq` bytes | Median compile time | Peak Python allocation |
|---|---:|---:|---:|---:|---:|---:|---:|
| GRE 2D | 640 | 383 | 382 | 0 | 30,663 | 0.119 s | 2.73 MiB |
| Spin echo 2D | 768 | 575 | 446 | 0 | 98,273 | 0.299 s | 10.62 MiB |
| Spiral DTI | 300 | 165 | 212 | 0 | 96,093 | 0.117 s | 3.14 MiB |
| EPI DWI | 1,294 | 487 | 212 | 644 | 96,283 | 0.300 s | 4.67 MiB |

Times are three local samples after environment import; allocation is a separate `tracemalloc` run
and therefore measures Python allocations, not total process RSS or native NumPy storage. File size,
performance, and the raw emitted-content digest are runner observations. The digest intentionally
retains exact floating-point bytes for same-runner comparisons, so it is not an exact cross-platform
assertion. Performance is compared on the same runner profile with a 20% tolerance.

The integration snapshot does make the following cross-platform exact assertions:

- total duration in integer picosecond ticks;
- emitted block count and block-duration digest;
- event counts and provenance-path digest;
- compile/check issue kinds and severities;
- whole-sequence m0, m1, and m2 summaries.

Full waveform equality is proved separately by `tests/compiler/fidelity.py`, which evaluates the
tree and emitted piecewise-linear waveforms on the union of their knots. The snapshot is not used as
a substitute for that independent oracle.

## Warning classification

| Warning class | Phase 0 disposition |
|---|---|
| RF delay increased to transmit dead time | Expected during current RF construction; counted in the artifact, not suppressed globally. |
| RF amplitude warning in conflict fixtures | Expected because those tests deliberately construct impossible overlaps; keep local to those fixtures. |
| Missing arbitrary-gradient first/last points | Test-fixture debt; explicit first/last values should be considered in an isolated cleanup. |
| PyPulseq shape-restoration deviation during spiral trajectory calculation | Known upstream/representation risk; physics and fidelity tests currently pass, so record and investigate separately. |
| Complex-cast warnings in optional reconstruction/example paths | Outside the build-only per-PR example tier; inspect in simulation/reconstruction validation. |
| SeqCraft `grad_merge`, `grad_resample`, raster, limit, and norm issues | Product report surface, not Python warnings; frozen by characterization tests and the artifact. |

## Regression gates

Every no-change compiler phase runs:

1. `ruff check .` and the scoped strict mypy command from `docs/testing.md`;
2. the non-heavy pytest tier with no new skip/xfail;
3. all source doctests;
4. compiler fidelity, invariants, errors, labels, triggers, and Phase 0 snapshots;
5. GRE, spin echo, spiral DTI, and EPI DWI integration tests;
6. the three build notebook smoke tests;
7. structural artifact comparison and same-runner performance comparison;
8. Linux/Windows and Python 3.11/3.12 GitHub Actions.

macOS is explicitly represented by this Phase 0 capture. Simulation/reconstruction notebooks need
MRzeroCore, Torch, and additional reconstruction dependencies; they remain a lab/nightly gate and
must run for phases that alter ADC, k-space, or bridge behavior.

## Local validation result

The Phase 0 candidate passed the following macOS/Python 3.12.13 gates on 2026-08-14:

| Gate | Result |
|---|---|
| Ruff correctness baseline | passed |
| Scoped strict mypy | 5 source files, passed |
| Non-heavy pytest/coverage tier | 658 passed, 2 skipped, 181 warnings, 82% coverage |
| Source doctests | 140 passed, 23 warnings |
| Build notebook smoke | 3 notebooks passed in an isolated temporary directory |
| Baseline capture | 4 recipes, 3 identical structural samples per recipe |

The two existing skips are zero-duration barrier cases in the generic module contract suite. The
warning classes are categorized above and are not globally suppressed. Simulation and reconstruction
notebooks were not run locally because the approved public/dev environment does not include the lab
simulation stack; they remain external Phase 0 evidence rather than a public PR gate.

## Exit status

Phase 0 local and remote gates are complete. GitHub Actions
[run #10](https://github.com/bughht/SeqCraft/actions/runs/31852265376) passed all seven jobs at
commit `7be79bf`: lint, types, examples, and Ubuntu/Windows tests on Python 3.11 and 3.12. At that
point the contract extraction was cleared to begin. Its final disposition and current paths are
recorded in [`phase1_contracts.md`](phase1_contracts.md).
