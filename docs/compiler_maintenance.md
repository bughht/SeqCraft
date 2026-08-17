# Maintaining the compiler

The compiler architecture is stable. Prefer changing a producer of `LogicBlock` objects when the
request is about sequence design; change the compiler only when the Pulseq transform itself needs a
new capability or a demonstrated correction.

## Before changing compiler code

1. Name the current owner: placement, boundaries, legalization, emission, or verification.
2. Write the input/output contract and invariant that changes.
3. Decide whether waveform, timing, block structure, labels, diagnostics, or the public API may
   change. If any may change, add an ADR and a focused characterization test first.
4. Keep one authoritative path. Do not add a compatibility adapter, feature-flagged compiler, or
   second IR for migration convenience.
5. Use raw PyPulseq events in compiler and integration fixtures. A user-facing module must not
   become a compiler-test dependency.

## Adding an event or extension

1. Classify the event in `design/events.py`: handled kinds are an explicit allowlist.
2. If support is intentionally deferred, add a named entry with actionable remedies to
   `compiler/placement.py::UNSUPPORTED_KINDS`. Unknown kinds must fail; none may disappear.
3. Define its intrinsic duration and reservation in `compiler/placement.py`. Keep active and
   reserved intervals distinct.
4. Decide its Pulseq block rules in `compiler/model.py`: whether it is indivisible and whether it
   is exclusive. These properties are independent; triggers are the canonical example.
5. Add scheduling or ownership policy to `compiler/boundaries.py`, transformation and limit policy
   to `compiler/legalization.py`, and only mechanical conversion to `compiler/emission.py`.
6. Add a minimal successful event fixture, rejection coverage, source immutability coverage, and
   the relevant label/trigger/fidelity/invariant tests. The PyPulseq source-inventory test must
   still account for every event type.
7. Run the complete gates in `testing.md`; compare the structural baseline and waveform oracle.

## Changing boundary policy

Boundary selection is a semantic policy, even when total duration and moments stay equal. Before
changing it:

1. Provide a concrete tree that the current midpoint/edge/barrier policy handles incorrectly or
   measurably poorly.
2. Record an ADR covering old and new block edges, tie-breaking, raster rounding, label targets,
   gradient splits, downstream block-index consumers, and rollback.
3. Add the failing fixture at the boundary stage and an integration characterization of block
   structure, waveform, labels, warnings, and errors.
4. Keep `emission.py` mechanical. A new boundary choice belongs in `boundaries.py`; any waveform
   consequence belongs in `legalization.py`.
5. Compare determinism, the current structural baseline, and same-runner performance before asking
   to merge.

## Review checklist

- The public entry remains `sc.compile` / `sc.compile_sequence`; stage functions and IR types are
  not added to `seqcraft.__all__` or `seqcraft.compiler.__all__`.
- Dependencies still run `errors -> design -> compiler -> analysis -> display`, with `scanner`
  independent of the compile path.
- Placement is the only tree traversal; legalization is the last policy-bearing stage; emission
  does not schedule or repair input.
- User-provided events and tree nodes remain read-only.
- Every fatal condition raises with stage/source context; every transformation warning has explicit
  singular and plural display text.
- Contract, PyPulseq, hardware, and against-tree verification remain always on unless measured
  evidence and an ADR approve a change.
- Documentation and tests identify the new maintenance owner without adding commentary that merely
  repeats clear code.

Reopen the architecture itself only for a reproducible severe bug, a measured performance
bottleneck, or an accepted ADR for a necessary contract change. Directory aesthetics and speculative
backends are not sufficient reasons.
