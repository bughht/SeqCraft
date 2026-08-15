# ADR-002: Compiler IR contracts and immutability

- Status: Accepted
- Date: 2026-08-14
- Applies from: Compiler refactor Phase 1

## Context

The compiler already had a frozen `_Placed` record, but its definition lived beside every compiler
algorithm and there was no typed boundary between legalization and PyPulseq emission. Creating a
second placement model would introduce two sources of truth; exposing mutable block assembly lists
would make later stages difficult to reason about and test independently.

## Decision

The private package defines two shallowly immutable contracts:

1. `PlacedEvent` records source-tree time, active interval, reservation interval, the read-only
   PyPulseq event reference, and its immutable source path. `compiler._Placed` is an import alias to
   this class during migration, so the authoritative placement algorithm produces one model.
2. `PulseqReadyBlock` records its contiguous index, absolute span, exact duration, tuple of events,
   all contributing source paths, and their common origin immediately before emission.

Both contracts provide concise summaries that do not expand waveform arrays. Their dataclasses are
frozen and compiler-owned collections are tuples. They do not deep-copy or freeze PyPulseq events;
all compiler stages must continue treating those user-provided objects as read-only.

`verify_placed_events` and `verify_ready_blocks` return structured internal contract violations.
The compatible `compile_sequence` path invokes them and raises a private internal error only if a
compiler stage produces a structurally impossible IR. Existing semantic verification of duration,
moments, labels, limits, and waveforms remains authoritative and unchanged.

These types live in `seqcraft.core._compiler` and are not re-exported by `seqcraft` or
`seqcraft.core`. They may evolve between releases while the public compiler façade remains stable.

## Consequences

Placement and emission now have independently constructible test boundaries without a second
compiler path. The adapters add one linear structural verification pass and shallow tuple
construction. Phase 0 performance guardrails determine whether that overhead is acceptable.

