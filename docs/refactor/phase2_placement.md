# Phase 2 placement extraction

- Status: Complete
- Date: 2026-08-14
- Branch: `refactor/compiler-phases`
- Parent: `824f380`
- Implementation: `1bef4ed`

## Scope delivered

- Moved the authoritative tree traversal and absolute-time resolution to
  `core/_compiler/placement.py`.
- Moved intrinsic event-duration interpretation and unsupported-kind rejection with that stage.
- Made the stage output an ordered immutable `tuple[PlacedEvent, ...]`.
- Kept `compiler._place` as a temporary identity alias so existing downstream instrumentation uses
  the extracted implementation rather than a parallel path.
- Kept negative-time, gradient-raster, exclusivity, boundary, label, gradient, and emission policy
  in the façade for their later phases.
- Added direct tests for nested offsets, insertion order, active/reservation intervals, provenance,
  input immutability, emission independence, and façade identity.
- Added machine-readable `stage` and `source_path` attributes to placement-created `CompileError`
  instances without changing their class or frozen human-facing text.
- Added the placement module to the strict mypy scope.

No public import, compiler signature, supported-kind policy, warning/error text, or unrelated
`core` module changed.

## Verification evidence

| Gate | Result |
|---|---|
| Placement, event-kind, contract, invariant, and boundary tests | 37 passed |
| Compiler test suite | 127 passed |
| Non-heavy test suite with coverage | 671 passed, 2 skipped; 82% coverage |
| Source doctests | 140 passed |
| Notebook smoke tests | 3 notebooks passed |
| Ruff | passed |
| Strict mypy for compiler helpers and typed utility modules | passed |
| Phase 0 structural baseline | exact match for all four recipes |

Same-runner performance ratios relative to the Phase 0 baseline:

| Recipe | Median wall-time ratio | Peak Python-memory ratio |
|---|---:|---:|
| `epi_dwi` | 1.087 | 1.000 |
| `gre_2d` | 1.015 | 1.029 |
| `se_2d` | 0.962 | 1.000 |
| `spiral_dti` | 1.009 | 1.002 |

All wall-time ratios remain below the 1.20 review threshold. GitHub Actions
[run 31856549261](https://github.com/bughht/SeqCraft/actions/runs/31856549261) passed all seven jobs:
lint, types, examples, and Ubuntu/Windows tests on Python 3.11/3.12.

## Exit status

All Phase 2 local and hosted gates are complete. The implementation remains in Draft
[PR #5](https://github.com/bughht/SeqCraft/pull/5); `main` is unchanged. The merge review passed and
the Phase 3 boundary was accepted in [`phase2_merge_review.md`](phase2_merge_review.md) and
[`phase3_boundary.md`](phase3_boundary.md).
