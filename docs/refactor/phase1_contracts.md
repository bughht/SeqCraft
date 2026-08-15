# Phase 1 compiler contracts

- Status: Implemented; hosted CI evidence pending
- Date: 2026-08-14
- Branch: `refactor/compiler-phases`
- Parent: `dd5ec6a`

## Scope delivered

- Added immutable `PlacedEvent` and `PulseqReadyBlock` contracts in `core/_compiler/model.py`.
- Reused `PlacedEvent` as the existing compiler's `_Placed` model instead of adding a parallel IR.
- Centralized event-kind groups and named time-comparison/duration helpers.
- Added structured placed-event and ready-block verifier skeletons.
- Adapted the authoritative assembly loop to construct and consume `PulseqReadyBlock` without
  moving boundary selection, gradient legalization, label targeting, or PyPulseq emission.
- Added strict type checking for both new private modules.
- Recorded the time and immutability decisions in ADR-001 and ADR-002.

No public import, compiler signature, event support policy, or unrelated `core` module changed.

## Verification evidence

Local macOS/Python 3.12.13 results:

| Gate | Result |
|---|---|
| Ruff correctness baseline | passed |
| Scoped strict mypy | 7 source files, passed |
| Focused contracts, invariants, and boundaries | 20 passed |
| Non-heavy pytest/coverage tier | 665 passed, 2 skipped, established warning classes, 82% coverage |
| Source doctests | 140 passed, 23 warnings |
| Build notebook smoke | 3 notebooks passed in an isolated temporary directory |
| Phase 0 structural comparison | all four recipes exactly equal |

The same-runner median compile-time ratios versus Phase 0 range from 0.966 to 1.042, below the 1.20
investigation threshold. Streaming ready-block verification keeps peak Python allocation ratios
between 1.001 and 1.025 rather than retaining the full intermediate block list.

The final phase commit also requires the complete Phase 0 regression gates and the hosted
Linux/Windows, Python 3.11/3.12 GitHub Actions matrix. Results are recorded here before the phase is
declared complete.

## Exit status

Pending the full local and hosted gates. Phase 2 must not begin until this section records their
results and the Phase 1 commit is pushed without changing `main`.
