# ADR-001: Compiler internal time policy

- Status: Accepted
- Date: 2026-08-14
- Applies from: Compiler refactor Phase 1

## Context

The compiler receives seconds as binary floating-point values, but Pulseq legality is defined on
scanner rasters. Direct addition accumulates drift, direct equality is brittle, and applying one
blanket tolerance to arithmetic and ordering can hide a genuinely off-raster event.

The existing `core.timing` module already provides picosecond integer ticks, exact addition and
subtraction, and raster-aware quantization. The monolithic compiler also has established
`EPS = 1 ns` comparison behavior that is frozen by the Phase 0 baseline.

## Decision

Compiler stages use one policy:

- External and IR time fields remain seconds as `float` for PyPulseq compatibility.
- Duration addition/subtraction and raster quantization use `core.timing` integer-tick helpers.
- Ordering and interval containment use the named helpers in `core._compiler.model`, all based on
  the existing `EPS` tolerance.
- Active intervals (`start`, `end`) and hardware reservations (`res_start`, `res_end`) remain
  distinct. A comparison must name which interval it is using.
- A stage must not silently snap a source event to a raster. Existing explicit raster errors and
  established block-boundary quantization remain unchanged.

Phase 1 defines and tests the policy without mechanically replacing every established comparison
inside `compiler.py`. Later stage extraction adopts the helpers as code moves; differential tests
must accompany any replacement whose boundary behavior could differ.

## Consequences

New stage code has a single vocabulary for equality, ordering, containment, and exact duration.
The compiler continues to interoperate with PyPulseq without adding a public time type. During the
incremental refactor, some frozen legacy comparisons remain inline until their owning stage moves.

