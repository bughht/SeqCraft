# Phase E — synthesis across C1, C2 and C3

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-22, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way.

**Date:** 2026-09-22 · **Skill:** v0, frozen and unmodified throughout the analysis
**Status:** recommendation for review. **Nothing implemented, nothing changed in
`.claude/skills/module-mining/`, no schema field added.**

Three fine scans, run against frozen v0 with the express purpose of seeing where it bends.

```text
C1  T2 preparation                 NEW_LEAF        / APPROVED_FOR_IMPLEMENTATION
C2  balanced SSFP                  NEW_KERNEL      / APPROVED_FOR_IMPLEMENTATION  (deferred pending C3)
C3  velocity encoding              NEW_LEAF        / APPROVED_FOR_IMPLEMENTATION
    flow compensation              EXTEND_EXISTING / APPROVED_FOR_IMPLEMENTATION
```

Four recommendations from three scans, none of them implemented. That ratio is the point: v0 was
used as a mining tool and it produced designs, not modules.

---

## The architecture these findings must not disturb

```text
Module  -->  build  -->  LogicBlock  -->  compose  -->  Tree  -->  Compiler  -->  Sequence
```

Everything below **refines what the Module layer contains**. Nothing below adds a layer between
`Module` and `LogicBlock`, or between `LogicBlock` and the compiler.

```text
modules know physics     tree knows times     compiler knows Pulseq
```

The rule that decides future cases:

> **"What waveform should physically exist?" → Module/kernel physical design.**
> **"How do these already-selected events become legal Pulseq blocks?" → compiler.**

The compiler must stay ignorant of VENC, m1 targets, flow compensation, b-value, RF isodelay
physics and steady-state physics. `LogicBlock` must stay a thin timed-event tree — `(t, event)`,
hierarchical, arbitrary overlap — with no moment requirements, optimisation constraints or
MRI-specific state.

---

## A. Should `evidence_class: domain-reference` become Skill vocabulary?

**Yes.** It is the one recommendation here that meets the Skill's own bar — *growth by reasoning
class rather than by sequence count* — because it changed the outcome three times **in three
different ways**:

| | corpus | what domain evidence did |
|---|---|---|
| **C1** T2Prep | **one** witness | **supplied a family** the corpus could not establish. Without it, YELLOW on witness count — which was the wrong standard |
| **C2** bSSFP | **three** independent witnesses | **corrected a contract they agreed on.** Two of three hard-wire `TE = TR/2`; the Handbook makes it the canonical realisation and not part of the balance invariant |
| **C3** flow | **one** witness, copyleft | **corrected an architectural inference.** The canonical joint solve is a quadratic with an explicit root, so a code-only reading would have built a solver layer for a closed-form problem |

C2 is the load-bearing one. If domain evidence only helped when the corpus was thin it would be a
fallback, and a fallback does not need vocabulary. C2 shows it corrects a *rich* corpus, because
implementations encode one set of choices and cannot say which are essential.

**Minimum change.** `evidence[].role` today describes a reference's relationship to *other
references* and every value presumes the reference is an implementation; `authoritative-external`
is the nearest fit and means something else. Add **one optional field**, `evidence_class`, with
three values and no new required field:

```text
domain-reference   handbook / authoritative review / classic paper
                   -> canonical concepts, semantic quantities, mature family boundaries
design-witness     an external implementation
                   -> corroborates a realisation, a timing convention, a mode
oracle             analytic calculation / simulation / independent measurement
                   -> tests whether OUR realisation produces the claimed physics
```

Plus one rule, which is the thing actually worth writing down:

> **The number of executable repositories is not a vote on whether a physical abstraction exists.**
> Code witnesses constrain implementation-specific claims; they do not determine the family.

And the curated-card workflow that already exists in `tools/module_mining/domain_evidence/` —
1–2 pages, provenance per claim, *establishes* and *does not establish* for every statement,
progressive escalation, and **report an unavailable source rather than reconstructing it**.

**Cost if wrong:** one optional field and a directory. Low.

## B. Does v0 need a concept for one scan producing several records?

**No. Prose and folder grouping are sufficient, and were.**

C3 produced two records because the candidate had two owners. `plans/flow_moments/` holds both
plus `findings.md` and `gap_note.md`, and a reader is not confused. A `scan_id` would be a field
invented from **one case**, which is the rule v0 sets for itself.

Recommend instead a sentence in `reference/artifacts.md`: *a fine scan may produce more than one
candidate record; group them in one folder with a shared findings document.* Documentation, not
schema.

**Revisit if** a third scan splits and the folder convention proves ambiguous.

## C. Placement kernel vs joint-realization kernel

**Clarify the existing definition. Do not create new runtime architecture.**

Every shipped kernel places: `GRE2DTR`, `GRE3DTR` and `TSEShot` decide offsets and hand each
leaf's events through unchanged. C2 and C3 both point at a second kind:

```text
PLACEMENT kernel            decides WHEN leaves happen; their waveforms are their own
JOINT-REALIZATION kernel    takes physical facts and requirements from several leaves and
                            designs the coupled waveform itself
```

Three things must be said precisely, and two of them are constraints:

1. **It works before materialisation.** A joint-realization kernel does not build child
   `LogicBlock`s, take them apart and rewrite their events. It holds the coupled physical design
   and emits the final events once. This is why it needs no new architectural layer — it is a
   Module that happens to consult several leaves' physics before building.
2. **It receives semantic facts and events, not cached derived quantities.** See E.
3. **It is still a Module.** `build()` still returns a `LogicBlock`. The tree and compiler never
   learn that anything unusual happened.

**Minimum change:** extend the ownership table's kernel row to name both kinds, and add the
before-materialisation constraint. **No new status, no new layer, no code.**

## D. Should the Skill require composition-level validation?

**Yes, and it is the second-strongest recommendation here.**

Flow compensation's central invariant is:

```text
the n-th moment of EVERY gradient on one axis,
from a semantic origin to the echo,
equals the target
```

That spans `Excitation`, `PhaseEncode`, `CartesianLine` and the compensation. Rules E and F are
both **per-module** — measure on the lattice *you* emit; validate occupancy after *your* events are
placed. Neither reaches a summed cross-module quantity, and `sc.moments` already takes a tree, so
the measurement exists while the rule that would require it does not.

Two clauses, both derived from cases rather than invented:

> **G. When a candidate's claim is a property of the composition, validate it on the composition.**
> A per-module check cannot see a quantity that several modules contribute to.

> **The validator integrates; it does not ask.** A module reporting its own `M1` and a validator
> checking that against the module's own target tests arithmetic, not physics.

**Cost if wrong:** one rule that occasionally applies to nothing. Low. **Cost if omitted:** a
GREEN flow-compensation module whose every constituent passes and whose axis sum is non-zero.

## E. The minimum architecture clarification justified by C2 + C3

Both scans found the same thing from opposite directions:

> **The most efficient physical waveform may span boundaries that SeqCraft currently assigns to
> separate leaf/module instances.**

Stated that way — not "one abstraction cannot own a waveform several leaves emit" — because in the
preferred resolution those leaves never emit the contested pieces. **The minimum clarification is
C's kernel distinction and nothing more.** Specifically **not**: `RequirementIR`,
`ConstraintSolver`, `AugmentationPlugin`, a generic `GradientDesignProblem`, or a backend protocol.
Those remain hypotheses, and C3 is a standing warning about designing from an inference.

Three supporting points, each a constraint on what gets built later.

**Keep three questions separate.** Answering one does not answer the others:

```text
WHAT must be true?    the physical requirement -- a target moment at a semantic instant
WHO owns it?          which abstraction holds the coupled design
HOW is it realised?   a closed-form construction, or a numerical backend
```

The same owner serves `total m1 at echo = 0` by analytic quadratic and a multi-constraint problem
by optimizer. Ownership and realization are independent axes.

**Moments are derived, not intrinsic.** A moment depends on placement time, semantic origin,
endpoint and order, so `M1` is not a property a `GradientEvent` can carry or a leaf can publish
truthfully. Moment analysis is an **internal helper**:

```text
(events on one axis, their placement times, semantic origin, endpoint, order)  ->  M0 / M1 / M2
```

If it ever becomes hot, cache primitive integrals (∫G dt, ∫tG dt, ∫t²G dt) *inside* that helper.
Do not put cache state in the event API without a measurement forcing it. And a leaf exposes its
**semantic instants and its events** — RF centre or isodelay, echo instant, the gradients it
intends — not a cached `m1` for one consumer.

**Physics-aware design is not compiler fusion.** C2 may be a case where two already-decided events
could be fused with *no change of physical semantics* — same integrated waveform, same RF and ADC
timing, same requirements. That would be a semantics-preserving compiler optimisation available to
every sequence. C3's merged waveform is *not*: its shape depends on VENC and a moment target.

```text
physics-aware joint realization       -> Module / kernel layer
semantics-preserving waveform fusion  -> possible future compiler optimisation
```

Which kind C2's is **has not been established.** Recovering its 220 µs by teaching the compiler
about bSSFP would be the second implemented as the first.

## F. Realization strategy

Record the principle; build nothing.

```text
physical design problem
        |
        +-- is there a supported canonical analytic construction?
              yes -> analytic fast path
              no  -> numerical backend
```

> **The physical problem definition belongs to SeqCraft. A realization backend does not define
> SeqCraft's public semantics.**

The Phase B inference and its correction, stated so a later reader cannot re-derive the error:

```text
WRONG       Open4DFlow uses GrOpt, therefore SeqCraft needs a general solver layer
CORRECT     the canonical cases do not justify introducing a general solver layer
ALSO WRONG  coupled gradient design never needs an optimizer
```

The Handbook gives closed forms for canonical constructions and **does not establish them for the
general problem** — several moment orders at once, several semantic windows, a minimum-TE search,
fixed waveform segments, non-zero endpoints, arbitrary pre-existing waveforms, eddy-current,
concomitant-field or PNS constraints. An optimizer is a legitimate backend there.

**GrOpt** stays evidence and a plausible future backend, not a dependency. When one is first
needed, write **the smallest adapter for that concrete problem** — not a universal interface
assuming every optimizer shares an input format. Evaluate a backend protocol only after a *second,
genuinely different* optimizer exists. Do not expose an optimizer's native structures as
SeqCraft's physical API.

**CI**, if an adapter is ever written: core CI stays SeqCraft-only, no GrOpt, canonical analytic
cases, blocking; an optional integration job may install and import GrOpt and run adapter tests.
Installing a GPL tool in a runner is not distribution. What needs separate licensing review is
vendoring it, making `pip install seqcraft` pull it, or shipping a combined wheel. **No CI change
is required now.**

**Testing an optimizer-backed realization**, when it exists: the optimizer proposes, SeqCraft's
independent validator measures — moments, amplitude, slew, fixed sections, timing windows,
semantic echo relationships. Analytic constructions serve as metamorphic references where they
exist. **Do not require sample-for-sample equality with optimizer output**, and distinguish
*feasible under the requested constraints* from *globally minimum TE* — a short solution is not a
proof of optimality. A complex waveform with no published reference is still testable: it
establishes that every requested constraint is satisfied, and does not establish global optimality,
equivalence to a canonical waveform, or scanner-level accuracy.

## G. Does anything require changing the existing architecture?

**No.** Confirmed against all three scans:

| | change needed? | |
|---|---|---|
| `LogicBlock` semantics | **no** | stays `(t, event)`, hierarchical, arbitrary overlap. No moments, constraints or MRI state |
| compiler responsibility | **no** | stays flatten, placement, boundaries, split/superpose/merge, legality, emit/verify. Learns nothing about VENC, m1, b-value, isodelay or steady state |
| `Module → LogicBlock → Tree → Compiler` | **no** | a joint-realization kernel is a Module whose `build()` returns a `LogicBlock` like any other |
| `GradientEvent` | **no** | moments stay derived; no `M0`/`M1` properties, no cache in the event API |

Every finding lands **inside** the Module layer:

```text
MODULE WORLD
  leaf modules          local MRI physical design
  kernels               placement, and -- new -- cross-leaf coupled physical design
  internal helpers      moment analysis; analytic waveform construction;
                        optional numerical realization adapters
        | build()
LOGICBLOCK WORLD        (t, event), hierarchical timed composition, arbitrary overlap
        |
COMPILER WORLD          flatten / placement / boundaries / split / superpose / merge /
                        Pulseq legality / emit / verify
```

---

## Recommended Skill v1, in priority order

| | change | evidence | cost if wrong |
|---|---|---|---|
| 1 | `evidence_class` (optional) + the not-a-vote rule + the curated-card workflow | C1, C2, C3 — three different mechanisms | one optional field |
| 2 | Rule G: validate composition-level claims on the composition; the validator integrates rather than asks | C3, with C2 adjacent | a rule that sometimes applies to nothing |
| 3 | Clarify the kernel row: placement vs joint realization, the latter designing **before** materialisation | C2 + C3 | a paragraph |
| 4 | Note that a fine scan may produce several records; group by folder | C3 | none |
| 5 | Record the analytic-first / optimizer-fallback principle and the ownership-vs-realization split | C3 | none — it forbids rather than builds |

**Not recommended:** `scan_id`; any requirement IR, constraint solver or plugin system; any change
to `LogicBlock`, the compiler or the top-level architecture; moments on events; a leaf publishing
cached moments; a backend protocol.

## What would change these conclusions

- **A**: a fine scan where domain evidence adds nothing a code witness did not, *and* the corpus is
  thin. Three for three so far, by three different mechanisms.
- **C**: discovering that a joint-realization kernel cannot in fact work before materialisation —
  that it genuinely needs to inspect built `LogicBlock`s. That would be a real architectural
  problem and is the assumption most worth attacking first.
- **D**: a composition-level claim that a per-module rule does catch, which would make rule G
  redundant rather than wrong.
- **F**: a case where the canonical construction is *not* closed form, which would move the
  optimizer from fallback to necessity for that family. **Diffusion is the obvious next test** —
  the MIT Pulseq witness inverts the b-factor analytically but its own comment marks trapezoid
  ramps `TODO`, and the copyleft witness holds the minimum-TE-for-target-b direction.

## What Phase E does not settle

Implementation of anything. `T2Prep`, `bSSFPTR`, `VelocityEncode` and the flow-compensation
extension are all approved designs and all unbuilt, deliberately. C2's 220 µs stays undecided
because deciding it needs the C2-vs-C3 fusion question answered first, and that is an ownership
question this document poses rather than closes.
