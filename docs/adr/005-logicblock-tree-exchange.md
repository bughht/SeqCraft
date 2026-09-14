# ADR-005: LBTX is a thin, one-way MATLAB-to-Python boundary

- Status: Superseded by [ADR-006](006-the-matlab-frontend-is-withdrawn.md)
- Date: 2026-08-24

**This decision no longer holds.** The frontend it describes was removed on 2026-09-14; the
document is kept because it records what was built and what it cost, which is the argument
ADR-006 rests on. The implementation is parked whole in
[`salvage/matlab-frontend/`](../../salvage/matlab-frontend/), frozen at `d27217d`.

## Context

SeqCraft's established Python path is already small:

```text
Python Module -> Python LogicBlock -> sc.compile -> pypulseq.Sequence
```

MATLAB users need the same Module and LogicBlock composition model while retaining the official
MATLAB Pulseq objects they already use. The first LBTX implementation instead added a scanner
wrapper, event builders, a per-event wire type system, Python writer/round-trip, semantic hash,
provenance, and structured diagnostics. Those layers repeated definitions already owned by MATLAB
Pulseq, PyPulseq, or SeqCraft and made the exchange feature resemble a second Pulseq implementation.

## Decision

LBTX version 1 is the smallest representation needed for this path:

```text
mr.opts + mr.make* -> MATLAB LogicBlock -> LBTX -> Python reader -> sc.compile -> .seq
                                                                         |
                                                                         v
                                                                    mr.Sequence
```

Python's normal compilation path does not use LBTX. The existing Python `Module`, `LogicBlock`,
compiler stages, public compile signature, warning policy, and return type are unchanged.

### Ownership

LBTX owns only:

- `format`, integer `version`, and compatible Pulseq version;
- a normalized snapshot of official `mr.opts` fields;
- definitions passed to the compiler;
- block tags, ordered nodes, and relative `start_s` values;
- three node variants: nested `block`, official Pulseq `event`, and SeqCraft `barrier`.

An event is the public struct returned by `mr.make*`, with field names normalized from camelCase to
snake_case. It keeps its official `type` and public fields directly. LBTX does not define event
constructors, per-event payload schemas, derived-field rules, or physical validation. MATLAB Pulseq
creates and measures the event; the existing Python compiler validates what it receives.

JSON numbers carry numeric values. A complex MATLAB value, principally RF `signal`, is represented
as `{"real": [...], "imag": [...]}`. This is the only value-level encoding not handled directly by
JSON.

### MATLAB user contract

- scanner: official `mr.opts`;
- events: official `mr.make*` structs;
- `seqcraft.LogicBlock < handle`: `add` mutates and returns the same block, matching Python;
- `LogicBlock.duration`: measured with `mr.calcDuration`, never declared;
- `LogicBlock.copy`: new block/node container, shared child items, matching Python's shallow copy;
- `seqcraft.Module`: `opts`, `tag`, and one build lifecycle;
- `seqcraft.compile`: writes `.seq`, reads it through `mr.Sequence.read`, returns `mr.Sequence`;
- `seqcraft.writeLBTX`: optional public exchange-file API, secondary to `compile`.

MATLAB does not emulate Python callable objects. The user calls `module.build(...)`. A subclass
implements protected `buildImplicit(...)`; the base `build(...)` calls it, checks that it returned a
`LogicBlock`, and supplies the module tag when the block is unnamed. `buildImplicit` is a confirmed
language-adapter name, not a second lifecycle. A future rename should migrate only this hook and its
documentation.

### Python process boundary

The Python exchange package contains a generic reader and structural schema. It restores event
objects as `SimpleNamespace` values without calling per-event constructors and checks their `type`
against the existing `HANDLED_KINDS` vocabulary. Scanner fields are matched to the public `Opts`
signature rather than copied into another scanner class.

The CLI has one operation, `compile-tree`. It prints ordinary errors and warnings and uses its exit
status; it does not introduce a diagnostic result model. MATLAB turns non-zero status into an
exception and successful compiler output into a MATLAB warning.

## Deferred

Version 1 intentionally has no:

- Python LBTX writer or mandatory Python LBTX path;
- MATLAB LBTX reader or round-trip;
- semantic hash;
- provenance or extensions framework;
- validation command or structured diagnostics protocol;
- native MATLAB compiler.

Each requires independent user evidence and an approved design change.

## Consequences

There is one compiler and one set of Pulseq event constructors. New official Pulseq public fields
can cross the language boundary without adding an event-specific branch, subject to the declared
Pulseq compatibility version and the existing compiler vocabulary.

The tradeoff is deliberately asymmetric: version 1 solves MATLAB authoring and compilation, not
general serialization. A consumer needing Python export, offline MATLAB import, provenance, or
machine-readable diagnostics must propose that feature rather than assuming it is bundled with
LBTX.
