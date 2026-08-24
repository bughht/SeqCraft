# ADR-005: LBTX is the versioned compiler-input exchange boundary

- Status: Accepted
- Date: 2026-08-24
- Protocol version: 0.1

## Context

Pulseq `.seq` is already portable between Python and MATLAB, but it is the output of the SeqCraft
compiler. It no longer contains the nested `LogicBlock` structure, relative starts, insertion-order
tie breaking, pre-legalization overlap, or tree provenance. A MATLAB frontend that builds SeqCraft
trees therefore needs a boundary before placement, not a `.seq` importer and not a second compiler.

Python objects are not that boundary. `LogicBlock`, `pypulseq.Opts`, and PyPulseq events contain
language-specific representation and mutable registration caches. Dumping their `__dict__` values
would turn implementation details into a protocol and make MATLAB depend on a particular PyPulseq
release.

## Decision

SeqCraft defines **LogicBlock Tree Exchange Format (LBTX)** as a versioned JSON document. Version
0.1 uses the schema identifier `seqcraft.logicblock`, SI units, and the provisional `.lb.json`
extension. The extension is a preview convention rather than a permanent promise.

The boundary is after a frontend or module has produced its tree and before compiler placement:

```text
Python or MATLAB builder -> LBTX -> Python importer -> existing sc.compile -> Pulseq .seq
```

The compiler remains the only authoritative implementation. It does not import the exchange
package, and no compiler stage or warning policy changes for LBTX.

### Time and numeric representation

Every time value is a base-10 string in seconds. Readers parse a finite decimal number; writers use
the shortest decimal that round-trips through binary64. This preserves the current in-memory model
while preventing a JSON parser from pre-rounding a timing value before the receiving language can
apply its raster policy. Other physical values are finite JSON numbers in their named SI units.

Integer ticks were rejected because one time base cannot express all scanner-defined rasters without
adding a second conversion policy. Rational pairs were rejected because they add substantial MATLAB
surface while the current tree and PyPulseq events are binary64 values. A future major version may
choose either if cross-language fixtures demonstrate a real loss.

### Schema and compatibility

- `schema` is exactly `seqcraft.logicblock`; `version` is independently versioned from SeqCraft.
- Version 0.1 rejects unknown semantic fields and unknown event types. Silent guessing is forbidden.
- Optional forward data belongs in `extensions`, whose names must be namespaced. Readers ignore it.
- A breaking semantic change requires a new major protocol version, an ADR, migration notes, and old
  fixtures. A minor reader may add optional fields only when their absence has an explicit default.
- Validation failures identify a JSON path. The CLI never requires MATLAB to parse a traceback.

### Scanner options

The schema lists the `Opts` values needed to reconstruct the compiler input and its physical
interpretation. Names include their units. Values are canonical (`Hz/m`, `Hz/m/s`, `Hz`, seconds),
not the display units originally passed to `Opts`. `rise_time` is intentionally absent: it is a
constructor rule for deriving `max_slew`, while the resulting limit is already explicit.

### Events

Version 0.1 covers delay, trapezoid and arbitrary/extended gradients, RF, ADC, label set/increment,
trigger, digital output, and SeqCraft barrier events. Each has a discriminated payload containing
physical fields only. Registration IDs, shape-library IDs, object identity, caches, and producer
classes are excluded.

RF complex samples use separate real and imaginary arrays. Gradient and RF sample times remain
explicit, so extended trapezoids and nonuniform RF shapes do not acquire a uniform-raster
assumption during exchange.

### Definitions, provenance, and semantic identity

`definitions` are compiler input and therefore semantic. `provenance` and `extensions` are not.
`semantic_hash(document)` validates and imports a document, re-exports its semantic content in
canonical form, and hashes sorted compact JSON. Producer version, git commit, wall-clock time, and
metadata consequently cannot change the physical identity.

Semantic equality means that canonical imported trees preserve node order, relative starts, nested
tags, event payloads, scanner options, and definitions. It does not mean byte-identical `.seq`
files; library IDs and writer metadata can differ without changing sequence meaning.

### Diagnostics and exit codes

`seqcraft validate-tree` and `seqcraft compile-tree` support human output and one JSON object on
stdout. JSON diagnostics contain severity, category, exception type, message, and an optional source
path. Exit codes are:

| Code | Meaning |
| --- | --- |
| 0 | Validated or compiled successfully |
| 2 | LBTX/schema/input error |
| 3 | SeqCraft compiler rejected the imported tree |
| 4 | File-system or unexpected adapter failure |

Compiler warnings are captured and returned as structured warning diagnostics without changing the
compiler's own warning contract.

## Consequences

MATLAB can use ordinary value objects and structs and invoke a stable CLI while Python owns all
placement and legalization. Protocol fixtures become the shared correctness source. Supporting a
new event requires an explicit schema/codec/conformance change, but cannot accidentally change the
compiler.

This does not provide compilation without Python. A native MATLAB or shared native compiler remains
deferred until deployment or measured performance evidence justifies reopening that decision.
