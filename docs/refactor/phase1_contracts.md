# Phase 1 compiler contracts

> **Historical implementation record.** The contracts remain current, but their original
> `core/_compiler/` paths and migration alias do not. They now live in `seqcraft/compiler/model.py`,
> with no legacy adapter; see [ADR-002](../adr/002-compiler-ir-contracts.md) and the
> [architecture freeze](compiler_architecture_freeze.md).

- Status: Complete
- Date: 2026-08-14
- Branch: `refactor/compiler-phases`
- Parent: `dd5ec6a`
- Implementation commit: `a609fa4`

## Scope delivered

- Added immutable `PlacedEvent` and `PulseqReadyBlock` contracts in the then-current
  `core/_compiler/model.py` path.
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

Hosted [GitHub Actions run 31854985214](https://github.com/bughht/SeqCraft/actions/runs/31854985214)
passed all seven jobs: lint, types, examples, and Ubuntu/Windows tests on Python 3.11 and 3.12.

## Exit status

All Phase 1 local and hosted gates are complete. Phase 1 was merged through
[PR #4](https://github.com/bughht/SeqCraft/pull/4) at merge commit `824f380`; the long-lived
`refactor/compiler-phases` branch was fast-forwarded to the same commit for Phase 2.

The later structure revision moved the contracts to `seqcraft/compiler/model.py`; the final
legalization extraction added `LegalizationResult` and removed the temporary `_Placed` alias.
