# Core package boundary audit charter

- Status: Charter accepted; findings deferred to Phase 7
- Scope: `src/seqcraft/core/`
- Implementation authority: None until a follow-up decision is accepted

## Purpose

The compiler refactor and a package-wide core cleanup solve different problems. Phases 1–6 turn the
compiler into an explicit, deterministic pipeline. This audit asks whether every other core module
has a coherent responsibility and dependency position after that pipeline exists.

The audit prevents two opposite mistakes: retaining unrelated code in `core` merely because it is
already there, and moving cohesive low-level code merely to make the directory look different.

## Fixed compiler boundary

The compiler refactor uses the following private layout:

```text
src/seqcraft/core/
├── compiler.py              compatible façade and orchestration
└── _compiler/               private stage contracts and implementations
    ├── model.py
    ├── placement.py
    ├── legalization.py
    ├── emission.py
    └── verification.py
```

Phases 1–6 may create and populate `core/_compiler`. They may not relocate other core modules.

## Audit scope

Phase 7 reviews `compiler`, `logic`, `events`, `system`, `geometry`, `timing`, `units`, `validate`,
`errors`, and `report`. For each module, record:

- its responsibility and owner;
- incoming and outgoing imports, including public import paths;
- whether it is required on the `LogicBlock → legal validated .seq` path;
- cohesion, test ownership, and circular-import risk;
- compatibility cost if its path changes;
- a recommendation: keep, split, move, or merge.

Recommendations require concrete evidence from the dependency graph and tests. File length or a
preference for a flatter or deeper directory is not sufficient evidence.

## Non-goals

This audit does not authorize:

- public import changes;
- broad file moves in compiler extraction commits;
- repository-wide typing, formatting, or naming cleanup;
- changes to `LogicBlock`, module, timing, waveform, or error semantics;
- speculative abstractions for future backends.

## Deliverables and decision rule

Phase 7 completes this document with a dependency graph, a module-by-module decision table, and
follow-up issues. Moves required to remove the legacy compiler path may remain in Phase 7. All other
approved structural changes require a separate post-refactor plan, compatibility review, ADR, and
PR sequence.

