# Compiler architecture freeze

- Status: Stable
- Date: 2026-08-17
- Branch: `refactor/compiler-phases`
- Starting point: PR #11 merge commit `74359032b782fae4dcd54385d8989d49e094cb77`
- Semantic change: None

## Decision

Freeze the current logical two-layer contract and the single staged implementation:

```text
any producer -> LogicBlock -> sc.compile(tree, opts) -> pypulseq.Sequence

compiler implementation:
placement -> boundaries -> legalization -> emission -> verification
```

`LogicBlock` is the only SeqCraft data model a producer must return. The compiler is a deterministic
constrained-scheduling transform, not a framework for user-facing modules. `pypulseq.Opts` remains
the scanner input and `pypulseq.Sequence` remains the output.

The architecture may be reopened only for a reproducible severe defect, a measured performance
bottleneck, or an accepted ADR for a necessary contract change. A preferred directory shape,
speculative backend, or desire for a more general IR is not sufficient evidence.

## Public and private import inventory

| Surface | Status | Compatibility meaning |
|---|---|---|
| `seqcraft.compile` | Public | Preferred root-level compiler entry point. |
| `seqcraft.compile_sequence` | Public | Identity alias for callers avoiding the builtin name. |
| `seqcraft.compiler.compile_sequence` | Public | The sole export in `seqcraft.compiler.__all__`. |
| Root compiler exceptions | Public | Catchable as `seqcraft.CompileError`, `HardwareLimitError`, `DefinitionConflict`, and `CompilerContractError`. |
| `compiler.placement`, `boundaries`, `legalization`, `emission`, `verification`, `model`, `errors` | Internal | Physically importable for repository tests, but not a compatibility surface. |
| `PlacedEvent`, `PulseqReadyBlock`, `LegalizationResult`, and stage functions | Internal | May evolve with their owning stage; never re-exported from the package root or compiler `__all__`. |

`tests/test_layering.py` enforces the public facade, internal dependency graph, absence of legacy
paths, and the broader package dependency order from source. Physical importability is deliberately
not treated as public API.

## Legacy and duplication audit

- There is one compile orchestration function and no old/new feature flag.
- Tracked source contains no `src/seqcraft/core/`, `core/compiler.py`, `core/_compiler/`, or nested
  `compiler/_compiler/` package.
- Temporary `_Placed` and `_place` migration aliases are gone.
- Placement alone traverses the source tree; legalization produces complete immutable ready blocks;
  emission accepts only those blocks and cannot import policy-bearing stages.
- The compiler-stage source import graph is acyclic.
- Historical documents retain old names only where they explain a dated decision. Each such record
  is marked historical or superseded and links to `architecture.md`, `compiler.md`, or this record.

Ignored local `__pycache__` directories may still contain bytecode from old layouts. They are not
tracked, imported, or part of the package artifact.

## Maintenance ownership

The current owner map and the checklists for adding an event/extension or changing boundary policy
are in [`compiler_maintenance.md`](../compiler_maintenance.md). The non-negotiable boundaries are:

- placement owns tree traversal and reservation construction;
- boundaries owns block-edge and label-target policy;
- legalization owns scheduling, gradient transformation, limits, and complete ready-block output;
- emission is mechanical and makes no policy decision; and
- verification independently checks stage contracts, PyPulseq legality, and tree equivalence.

## Verification

Local macOS/Python 3.12.13 results:

| Gate | Result |
|---|---|
| New API/layering/legacy-path guards | 14 passed |
| Non-heavy suite with coverage | 366 passed, 6 optional `pulseq_systems` skips; 89% coverage |
| Source doctests | 44 passed |
| Ruff | passed |
| Strict mypy | 9 source modules passed |
| Executable API reference | 51 Python blocks and the complete `__all__` index passed |
| Getting-started notebook smoke | passed |
| Current structural baseline | exact match for GRE and spin echo |
| Realistic 128 × 8 GRE budget | passed the 60-second guardrail |

A three-iteration same-runner capture measured median compile times of 0.119 s for GRE and 0.284 s
for spin echo, with peak Python allocations of 2,752,224 and 10,984,439 bytes respectively. These
match the immediately preceding audit and show no performance regression from the documentation and
source-level test changes.

[GitHub Actions run 31994234063](https://github.com/bughht/SeqCraft/actions/runs/31994234063)
passed all seven jobs on PR #12: tests on Ubuntu and Windows with Python 3.11 and 3.12, Ruff,
strict mypy, and executable examples.
