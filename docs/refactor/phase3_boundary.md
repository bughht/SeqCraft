# Phase 3 legalization boundary

> **Historical.** This records the state at the time it was written.  The compiler was
> subsequently changed by [ADR-004](../adr/004-compile-returns-a-sequence.md): `compile` returns a
> bare `pypulseq.Sequence`, every legality failure raises, and `CompiledSequence`, `Report` and
> `Issue` no longer exist.  For the current shape see [`../compiler.md`](../compiler.md) and
> [`../architecture.md`](../architecture.md).

- Status: Accepted implementation boundary; implementation not started
- Date: 2026-08-14
- Applies after: Phase 2 merge
- Semantic-change policy: None

## Purpose

Phase 3 extracts the current legalization algorithm without changing any boundary, conflict,
gradient, label, warning, or error decision. The stage turns verified placed events into verified
PyPulseq-ready blocks. It may construct or derive individual PyPulseq event objects, but it must not
construct, mutate, or inspect a PyPulseq `Sequence`.

```text
LogicBlock
  -> placement -> tuple[PlacedEvent, ...]
  -> legalization -> immutable legalization result
       - tuple[PulseqReadyBlock, ...]
       - label assignment metadata
       - transformation warnings
       - current-policy boundary trace
  -> façade / later emission stage -> PyPulseq Sequence
```

## Input boundary

Legalization accepts only:

- an ordered immutable sequence of verified `PlacedEvent` objects;
- the minimum read-only hardware context needed for event and block rasters, dead times, limits,
  and maximum block duration;
- source diagnostic metadata only when required to preserve an existing error message.

It does not accept a `LogicBlock`, `Geometry`, sequence name, definitions, `CompiledSequence`, or a
partially built PyPulseq `Sequence`. Placement remains the only tree traversal owner.

## Output boundary

The stage returns one immutable result containing:

- contiguous `PulseqReadyBlock` objects covering the complete legalized duration;
- immutable label-target metadata required by the existing address verifier;
- warnings produced by legalization itself, including current raster, orphan-label, and gradient
  resampling warnings;
- an immutable trace of the current boundary decisions sufficient to explain which existing rule
  selected each edge.

The stage verifies the complete ready-block tuple before returning. It must not expose builder lists,
gradient sweep cursors, active-piece lists, or mutable issue accumulators.

## Current function ownership

| Current responsibility | Phase 3 owner | Notes |
|---|---|---|
| Negative reservation, off-raster gradient, and RF/ADC overlap checks | Legalization preflight | Preserve exception type, text, ordering, and first-offender selection. |
| Total duration and block-raster normalization | Legalization | Preserve ceiling and fixed-width maximum-block rules. |
| `_Spans`, `_covering`, barrier/gap helpers, `_boundaries` | Boundary selection within legalization | Extract mechanically; do not introduce Phase 4 boundary categories or scoring. |
| Label targeting, order conflicts, and orphan warnings | Scheduling within legalization | Preserve barrier behavior and current ADC assignment semantics. |
| Singleton scheduling and `_in_block_delay` | Ready-block assembly | Derived event objects are stage-owned; source events remain read-only. |
| `_superpose`, `_as_arbitrary`, `_axis_gradient`, `_resampled` | Gradient legalization | Refactor hidden issue-list mutation into explicit returned events and warnings. |
| `_adc_conflict` | Gradient legalization safety check | Preserve the protected ADC-window failure path. |
| `_required_duration`, `_common_path`, ready-block construction | Legalization exit validation | A returned ready block must already fit its interval and carry provenance. |
| `_limit_issues` | Façade until Phase 5 | Hardware checks are not part of the Phase 3 extraction. |
| `pp.Sequence`, `add_block`, and PyPulseq error translation | Façade until Phase 5 emission | Legalization must not call these paths. |
| Definition merge and `CompiledSequence` construction | Compiler façade | They are orchestration, not legalization. |
| `_expected_addresses` and `CompiledSequence._verify` | Existing semantic verification | Move or consolidate only in Phase 6. |

## Required extraction order

Phase 3 uses the existing `refactor/compiler-phases` branch and one phase PR with three reviewable
commits:

1. extract current boundary and conflict selection behind a pure input/output seam;
2. extract gradient split and superposition with explicit warning output;
3. assemble and verify the immutable legalization result, then remove the corresponding façade
   implementation.

Each commit must keep the compiler façade on the single new implementation. Temporary identity
aliases are allowed for internal tests, but parallel algorithms or old/new feature flags are not.

## Frozen behavior and non-goals

Phase 3 must not:

- change midpoint, edge, raster rounding, barrier, maximum-block, or tie-breaking decisions;
- introduce the hard/natural/soft/forbidden policy or a new boundary score;
- change block count, block edges, gradient samples, moments, label ownership, warnings, or errors;
- move mechanical sequence emission or final semantic verification;
- relocate unrelated `core` modules or change public imports.

Those changes belong to Phase 4 or later and require their own evidence and review.

## Definition of ready

- [x] Input and output ownership are explicit.
- [x] Existing functions have a stage owner.
- [x] Semantic changes are forbidden.
- [x] Hidden mutation to remove is identified.
- [x] Rollback is the Phase 2 merge commit.
- [x] Review units and required evidence are defined.

Phase 3 coding may begin only after Phase 2 is merged and the long-lived refactor branch is
synchronized to that merge commit.

## Exit evidence

Every Phase 3 commit and the final phase candidate require:

- direct synthetic-`PlacedEvent` tests for the extracted seam;
- exact Phase 0 structural comparison, including block edges, block count, split count, events,
  warnings, errors, labels, and waveform hashes;
- source-event immutability checks;
- compiler tests, non-heavy suite, source doctests, and notebook smoke tests;
- strict typing and Ruff;
- same-runner performance review against Phase 0;
- hosted lint, types, examples, and Ubuntu/Windows Python 3.11/3.12 checks.
