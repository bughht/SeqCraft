# Phase 2 merge review

- Status: Passed; PR #5 merged as `3e4024b`
- Date: 2026-08-14
- Pull request: [#5](https://github.com/bughht/SeqCraft/pull/5)
- Reviewed head: `3e91507`
- Base: `824f380`

## Review result

No blocking correctness, compatibility, maintainability, or test-coverage finding remained. The
review cleared PR #5, which was subsequently merged. The boundary recorded in
[`phase3_boundary.md`](phase3_boundary.md) is a historical proposal; the delivered boundary is
[`phase3_legalization_emission.md`](phase3_legalization_emission.md).

## Scope reviewed

- Compared the extracted traversal, duration, reservation, and unsupported-kind branches with the
  former authoritative implementation in `core/compiler.py`.
- Traced every consumer of the placement output through boundary selection, label scheduling,
  gradient legalization, semantic verification, and compiler tests.
- Checked the private-module and public-import boundaries.
- Reviewed direct placement tests for ordering, nested offsets, timing, provenance, source-event
  immutability, error context, emission independence, and façade identity.
- Confirmed the structural baseline, local suites, and hosted Linux/Windows matrix recorded in
  [`phase2_placement.md`](phase2_placement.md).

## Compatibility assessment

The extraction preserves the event-kind policy, timing calculations, reservation spans, exception
classes, exception text, public compiler signature, and emitted sequence behavior. Two intentional
internal contract changes are approved:

1. placement returns `tuple[PlacedEvent, ...]` instead of the former private mutable list;
2. placement-created `CompileError` objects carry `stage` and `source_path` attributes while their
   existing human-readable messages remain unchanged.

At review time, `compiler._place` remained an identity alias to the single extracted
implementation. It was a private migration seam, not a second path, and has since been removed.

## Merge gate

- [x] Implementation has one authoritative placement traversal.
- [x] Placement does not construct or mutate a PyPulseq `Sequence`.
- [x] Source events remain read-only.
- [x] Phase 0 structural output is unchanged.
- [x] Same-runner performance remains within the 1.20 review threshold.
- [x] GitHub Actions passed all seven jobs at the reviewed head.
- [x] Phase 3 ownership and stage boundaries are recorded before merge.
