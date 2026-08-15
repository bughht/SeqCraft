# ADR-000: Compiler refactor scope and compatibility policy

- Status: Accepted
- Date: 2026-08-14
- Last amended: 2026-08-14 (private package location and core-audit boundary)
- Baseline: `afe576bafe17e9cce3edfc768b343ad958994560`

## Context

The compiler currently concentrates placement, label targeting, boundary selection, gradient
legalization, PyPulseq emission, reporting, and verification in one module and one assembly loop.
This keeps the difficult logic centralized, but makes a local scheduling change hard to review and
hard to test at its actual stage boundary.

The existing implementation already supports real EPI, spiral, diffusion, RF, ADC, label, trigger,
provenance, and round-trip workflows. Refactoring without first freezing those outputs risks
silently changing scanner waveforms while producing a sequence that still compiles.

## Decision

The refactor will proceed in numbered phases on the long-lived `refactor/compiler-phases` branch.
Each completed phase has a dedicated commit, complete exit evidence, and an explicit recommendation
about merging that phase to `main`.

Phases 0–3 and the extraction portions of Phase 5 are no-change work. They must preserve:

- public imports, signatures, result properties, and error classes;
- supported and explicitly unsupported event types;
- waveform, absolute timing, block structure, label state, provenance, and definitions;
- warning/error kinds, severities, and actionable messages;
- deterministic output and source-tree/event immutability;
- performance within the Phase 0 guardrails.

Phase 1 must account for the existing frozen `_Placed` dataclass rather than create a parallel
placement model. A ready-block contract may be introduced behind adapters, but the current
`compile_sequence` path remains authoritative until differential tests prove the replacement.

`src/seqcraft/core/compiler.py` remains the compatible façade. Private contracts and extracted
stages live under `src/seqcraft/core/_compiler/`; they are not placed at the top level of
`seqcraft`, and they are not public imports. Phases 1–6 may organize this private compiler package
but may not relocate unrelated core modules.

The broader question of whether every current core module belongs in `core` is handled by a
dependency- and cohesion-based audit in Phase 7. The audit records public/private imports, ownership,
test boundaries, and keep/split/move/merge recommendations. Except for changes required to finish
the single compiler path, its recommendations require a separate post-refactor plan and ADR before
implementation. The audit charter is `docs/refactor/core_package_boundary_audit.md`.

An intentional behavior change requires all of the following:

1. a separate ADR describing the old behavior, new behavior, and MRI/hardware rationale;
2. focused characterization or property tests demonstrating the difference;
3. waveform, timing, block, metadata, error, and performance baseline review;
4. explicit reviewer approval before updating the versioned artifact.

Known bugs are recorded during Phase 0 and are not fixed opportunistically in extraction commits.
Dependency upgrades, repository-wide formatting, and broad typing cleanup remain separate changes.

## Consequences

The refactor carries more adapters and differential tests temporarily, and some duplicated control
flow may exist between phases. In return, each extraction is reviewable, regressions identify the
stage that introduced them, and policy changes cannot hide inside mechanical moves.

The Phase 0 artifact is deliberately partly exact and partly observational. Stable structural,
timing, event-count, provenance, issue, and moment fields are cross-platform tests. The exact raw
floating-point content digest, file size, wall time, and Python allocation are observations tied to
a matching runner profile. A 20% compile-time regression requires investigation and approval, not
an automatic artifact rewrite.

Merging a phase does not require abandoning the long-lived branch. After a phase commit is merged,
the same branch continues with the next phase; the branch must retain the merged history and must not
rewrite already reviewed phase commits without an explicit reason.
