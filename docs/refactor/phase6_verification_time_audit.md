# Phase 6 verification, warning, and time-policy audit

- Status: Complete locally; hosted CI pending
- Date: 2026-08-16
- Branch: `refactor/compiler-phases`
- Baseline: `00fe41e`
- Policy change: None

## Decision

Keep every existing verification layer always on. The layers answer different questions, their
known failure modes are independently tested, and the realistic compile remains far inside the
existing performance budget. Do not add verification levels, a debug snapshot framework, unified
error codes, or a compiler-wide integer-time IR.

One concrete warning-contract defect was fixed: a non-empty transformation-note category without
an entry in `_WARNING_TEXT` used to fall back silently to its internal key, despite the documented
rule that every category must have explicit singular and plural text. It now raises
`CompilerContractError` before any warning is emitted.

## Verification ownership

| Check | Owner | Question | Why it is not redundant |
|---|---|---|---|
| Placed-event contract | placement boundary | Are absolute active/reservation intervals and provenance structurally valid? | Runs before scheduling and sees the placed IR, which PyPulseq never sees. |
| Ready-block contract | legalization boundary | Are block indices, intervals, exclusivity and provenance structurally valid across the complete immutable result? | Catches compiler-stage defects before PyPulseq mutation and names the private IR boundary. |
| Limit and required-duration checks | legalization | Are the combined waveforms within per-axis limits, and can every event fit its selected block? | Requires the post-superposition waveform and compiler-selected duration. |
| Event-size and duplicate-address checks | finished sequence | Can the interpreter hold each RF/ADC event, and are imaging addresses unique? | These facts exist only after PyPulseq registration and label folding. |
| `Sequence.check_timing()` | PyPulseq backend | Will the emitted backend representation pass Pulseq timing rules? | It validates PyPulseq's registered representation rather than SeqCraft IR. |
| Against-tree verification | final compiler stage | Did emission preserve duration, m0, m1 and the label address of every ADC? | It compares independently traversed source and emitted structures and catches legal-looking compiler corruption. |

The apparent overlaps are intentional:

- ready-block duration checks protect the legalization contract; PyPulseq timing checks the
  registered backend representation;
- duplicate-address checking catches collisions, while against-tree label checking catches a
  unique but shifted address sequence;
- m0 catches lost area, while m1 catches a waveform moved in time; and
- the explicit total-duration invariant is still required because PyPulseq's `TotalDuration`
  complaint contains a known float-equality artifact and does not compare against the source tree.

No contradictory owner or unreachable verification branch was found.

## Warning and error surface

The five warning-note categories have one owner and explicit display text:

| Category | Producer | Meaning |
|---|---|---|
| `snap` | orchestration | RF/ADC reservation moved to the block raster |
| `orphan_label` | boundary/label policy | label has no ADC after it |
| `merge` | legalization | same-axis gradients were superposed |
| `resample` | legalization | an otherwise unrepresentable mixed-knot waveform was resampled |
| `norm` | legalization | vector norm exceeds a per-axis limit but remains legal on real amplifiers |

Warnings remain aggregated once per category and are emitted only after every fatal and semantic
check has passed. Fatal user inputs remain `CompileError`, machine ceilings remain
`HardwareLimitError`, conflicting definitions remain `DefinitionConflict`, and broken private IR
or warning contracts remain `CompilerContractError`. There is no consumer that needs stable error
codes or debug snapshots, so adding either would create an API without a use case.

## Time-policy coverage

The current hybrid policy is deliberate: public and PyPulseq-facing times remain seconds as
floats, while quantization, interval subtraction, waveform-knot union and raster reconstruction
use integer picosecond ticks. A compiler-wide integer-time migration would duplicate PyPulseq's
float representation and add conversions at every boundary without a failing fixture to justify
it.

Existing tests cover the observed failure modes:

- an exact raster multiple must not be ceiled by one extra raster;
- `n * raster` reconstruction must return the canonical decimal float;
- repeated duration addition must not drift off-raster;
- negative nearest/floor/ceil behavior must be symmetric and ordered;
- non-Siemens and non-integer-microsecond rasters must compile;
- a gradient starting off its raster must fail instead of being silently moved;
- block ends must ceil rather than truncate RF ringdown or ADC dead time;
- in-block RF, ADC and gradient delays must land on their own rasters;
- arbitrary-gradient centre samples and edge knots must survive split/superposition;
- block seams must remain amplitude-continuous; and
- m1 must detect a whole-raster time displacement that m0 cannot see.

No new float/raster regression was found, so the time representation is unchanged.

## Performance evidence

A temporary profiler wrapped the existing stage functions without changing their results. Each
recipe reused one pre-built tree so the numbers measure compilation rather than fixture creation.
The realistic GRE contains 10,240 tree nodes (`matrix=128`, `n_slices=8`). Results are local
medians on the same Python 3.12 `seqcraft-dev` environment:

| Recipe | Iterations | Median compile | Against-tree share | Emitted-moment share | PyPulseq timing share | Label-fold share |
|---|---:|---:|---:|---:|---:|---:|
| GRE default | 7 | 0.118 s | 52.1% | 26.1% | 2.6% | 2.6% |
| Spin echo default | 5 | 0.284 s | 25.8% | 12.9% | 2.1% | 2.8% |
| GRE 128 × 8 slices | 3 | 1.909 s | 51.6% | 26.0% | 2.6% | 2.5% |

The semantic verifier has measurable cost, but the realistic compile is still over 30 times below
the 60-second integration budget. Every invariant has a test that deliberately makes it fail, so
the cost is buying demonstrated coverage rather than duplicate assertions.

Two consolidation prototypes were rejected before commit:

1. compute m0 and m1 during one emitted-block traversal; and
2. reuse one `evaluate_labels()` result for uniqueness and against-tree checks.

They reduced call counts but produced no repeatable wall-time improvement: moment integration,
not block reconstruction, dominates the first path, while PyPulseq already caches the label fold.
The realistic single-run candidate moved from 1.887 s to 1.940 s, within noise and in the wrong
direction. Keeping the simpler existing functions is therefore preferable.

## Local verification

| Gate | Result |
|---|---|
| Compiler plus time-policy tests | 170 passed |
| Non-heavy suite with coverage | 363 passed, 6 optional-dependency skips; 89% coverage |
| Source doctests | 44 passed |
| Ruff | passed |
| Strict mypy | 9 source modules passed |
| Executable API reference | 51 Python blocks and full `__all__` index passed |
| Getting-started notebook smoke | passed |
| Phase 0 structural baseline | exact match through the integration baseline tests |
| Phase 0 local capture | GRE: 383 blocks, 0.126 s median, 2,760,314 peak Python bytes |
| Phase 0 local capture | spin echo: 575 blocks, 0.594 s median, 10,978,457 peak Python bytes |

The six skips require the optional `pulseq_systems` dependency and are unchanged from the existing
local environment. Hosted Linux/Windows and Python 3.11/3.12 validation remains the final gate.

## Phase disposition

- Keep all verification always on; no verification-level option.
- Keep the current warning/error classes; no code registry or snapshot framework.
- Keep the current float/tick boundary; no compiler-wide time migration.
- Retain only the explicit warning-category contract fix and its regression test.
- Run the normal local and hosted gates; Phase 7 may then freeze the documented architecture.
