# ADR-002: Compiler IR contracts and immutability

- Status: Accepted
- Date: 2026-08-14
- Applies from: Compiler refactor Phase 1
- Last amended: 2026-08-17 (final package location and completed stage contracts)

## Context

The compiler already had a frozen `_Placed` record, but its definition lived beside every compiler
algorithm and there was no typed boundary between legalization and PyPulseq emission. Creating a
second placement model would introduce two sources of truth; exposing mutable block assembly lists
would make later stages difficult to reason about and test independently.

## Decision

The private package defines three shallowly immutable contracts:

1. `PlacedEvent` records source-tree time, active interval, reservation interval, the read-only
   PyPulseq event reference, and its immutable source path.
2. `PulseqReadyBlock` records its contiguous index, absolute span, exact duration, tuple of events,
   all contributing source paths, and their common origin immediately before emission.
3. `LegalizationResult` carries the complete ready-block tuple and immutable transformation notes
   across the legalization/emission boundary.

Both contracts provide concise summaries that do not expand waveform arrays. Their dataclasses are
frozen and compiler-owned collections are tuples. They do not deep-copy or freeze PyPulseq events;
all compiler stages must continue treating those user-provided objects as read-only.

`verify_placed_events` and `verify_ready_blocks` return structured internal contract violations.
The single `compile_sequence` path invokes them and raises `CompilerContractError` only if a compiler
stage produces structurally impossible IR. Semantic verification of duration, moments, labels,
limits, and waveforms remains independently owned and always on.

These types live in `seqcraft.compiler.model`. They are not re-exported by `seqcraft` or included in
`seqcraft.compiler.__all__`; physical importability is not a compatibility promise. They may evolve
between releases while the public compiler façade remains stable.

## Consequences

Placement, legalization, and emission have independently constructible test boundaries without a
second compiler path. The final implementation has no migration aliases or legacy adapter. Its
always-on structural verification cost remains inside the recorded performance guardrails.
