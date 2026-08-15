# Phase 2 placement extraction

- Status: Implemented; local evidence complete, hosted evidence pending
- Date: 2026-08-14
- Branch: `refactor/compiler-phases`
- Parent: `824f380`

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

All wall-time ratios remain below the 1.20 review threshold. The final candidate still requires the
hosted Linux/Windows Python 3.11/3.12 matrix.

## Exit status

Pending the hosted gates. Phase 3 must not begin until the evidence is complete and the Phase 2
merge recommendation has been reviewed.
