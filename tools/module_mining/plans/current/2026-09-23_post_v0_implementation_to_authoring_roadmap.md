# Post-v0 implementation → sequence-authoring roadmap

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-23, so that the pull request carries the direction behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../../README.md) for which one to edit.

**Status:** current direction after PR #36
**Date:** 2026-09-23

## Purpose

The post-v0 research phase has done enough architectural work. Module Mining v1, the C1/C2/C3 fine
scans, and the diffusion architecture stress test have established a usable boundary:

```text
reusable MRI physics / policy
        ↓
Module layer
    leaf / placement kernel / joint-realization kernel
        ↓
LogicBlock tree
        ↓
Compiler
        ↓
PyPulseq Sequence
```

The diffusion implementation demonstrated that a coupled physical design can remain inside the
existing Module layer and be completed **before materialisation**, without introducing a
`RequirementIR`, a general constraint-solver framework, new `LogicBlock` semantics, or
candidate-specific compiler behaviour.

The next milestone is therefore not another architecture-search phase. It is to **cash out the four
post-v0 designs already approved for implementation**, then freeze a useful package-level module
baseline and build the first `sequence-authoring` skill on top of it.

The intended order is:

```text
T2Prep
  ↓
VelocityEncode
  ↓
Flow compensation
  ↓
bSSFPTR
  ↓
sequence-authoring skill v0
  ↓
future module mining becomes primarily demand-driven
```

This document is an implementation roadmap, not a replacement for the candidate records. The
candidate YAMLs, findings, domain-evidence cards, and the module-mining Skill remain the detailed
sources for each physical contract.

---

## Starting point

The following work is already complete:

- Module Mining v0 was calibrated and evolved to v1 from observed evidence.
- Batch corpus mining and the T2Prep, bSSFP, velocity-encoding, and flow-compensation fine scans
  are complete.
- `DiffusionSEPrep` has been implemented as the architecture stress test.
- The diffusion result supports joint physical design in the Module/kernel layer before child
  `LogicBlock`s are materialised.
- Peak-B1 enforcement is now a package-wide contract for every RF-producing path, with a compiler
  backstop.
- Compiler and `LogicBlock` architecture are not waiting on another redesign before the approved
  queue can be implemented.

The approved but unbuilt queue is:

| capability | decision | intended layer |
|---|---|---|
| `T2Prep` | `NEW_LEAF / APPROVED_FOR_IMPLEMENTATION` | `preparation/` |
| `VelocityEncode` | `NEW_LEAF / APPROVED_FOR_IMPLEMENTATION` | `encoding/` |
| flow compensation | `EXTEND_EXISTING / APPROVED_FOR_IMPLEMENTATION` | waveform owner; kernel only for cross-leaf coupling |
| `bSSFPTR` | `NEW_KERNEL / APPROVED_FOR_IMPLEMENTATION` | `kernel/` |

The purpose of the next phase is to implement these contracts, not to maximise module count.

---

## Rules for the implementation phase

### 1. Treat approved candidate records as specifications, but audit them before coding

The fine scans contain the intended physical contract, ownership, exclusions, and validation
claims. They should not be re-mined from scratch.

However, an implementation PR must begin with a **spec-consistency pass**. Research records evolved
over several revisions and may contain stale wording. If two sections disagree, reconcile the
record first; do not let the implementation silently choose one interpretation.

A known example already exists in the T2Prep record: the decision and semantic-timing sections
define `prep_time_s` as effective-tip-down-rotation → effective-tip-up-rotation, reducing to
RF-centre → RF-centre for the initial hard/composite mode, while one acceptance bullet still
contains the older tip-down-end → tip-up-start wording. The record must be made internally
consistent before that acceptance test is encoded.

### 2. Do not repeat fine scans unless implementation exposes genuinely new evidence

The four capabilities are already approved. Implementation may reveal a local ambiguity or defect,
but the default action is:

```text
approved contract
    ↓
implement
    ↓
measure emitted sequence
    ↓
complete acquisition
    ↓
Layer 3 where the user-facing physical claim requires it
```

Do not reopen candidate discovery simply because code is now being written.

### 3. Do not invent generic infrastructure from one implementation

Still prohibited without new evidence:

- `RequirementIR`;
- a generic constraint-solver layer;
- an augmentation/plugin framework;
- new semantic responsibilities in `LogicBlock`;
- candidate-specific behaviour in the compiler;
- one wrapper class per named source sequence.

Closed-form canonical cases should remain closed form. A future numerical backend can be justified
by a real case that needs it.

### 4. Correctness is measured on emitted compositions

Use the existing module-mining validation rules:

- intrinsic module invariants are measured on emitted events;
- occupancy/timing is checked after placement;
- composition-level moment or balance claims are validated on the full composition;
- the validator integrates or measures the emitted waveform rather than trusting module-reported
  derived values;
- when a new Layer-1 instrument is introduced, Layer 3 should independently challenge it where
  practical.

### 5. Each capability should leave a complete teaching path

A shipped abstraction should normally leave:

```text
production implementation
        +
focused tests
        +
public API/docs
        +
a complete acquisition example
        +
Layer-3 evidence when the physical claim requires it
```

Examples teach the sequence and its physics. Development-history arguments stay in findings /
plans / PR text.

---

# Stage A — implement `T2Prep`

## Scope

Implement **only the approved `mlev4-hard` mode**.

The initial public abstraction is a nonselective T2-weighting preparation leaf:

```text
tip-down
    ↓
MLEV-4 composite refocusing train
    ↓
tip-up
    ↓
spoiler
```

The public target is `prep_time_s`, with semantic timing defined as:

> effective tip-down rotation instant → effective tip-up rotation instant.

For the initial hard/composite mode this reduces to centre-to-centre timing.

The approved mode includes:

- hard/non-adiabatic RF;
- four composite refocusing pulses at the approved fractional locations;
- the MLEV-4 `+ + - -` phase pattern;
- the witnessed composite tip-up;
- nonselective operation;
- spoiler after the tip-up, outside `prep_time_s`;
- caller-selectable spoiler axes/moment only to the extent already approved by the record.

## Explicitly out of scope

Do not add:

- adiabatic BIR-4/AHP modes;
- B0/B1-robustness claims;
- arbitrary refocusing-count APIs;
- `mlev2-hard`;
- the G1/G2 cross-boundary robustness variant;
- diffusion-preparation behaviour;
- a new generic composite-RF framework unless a second real consumer independently justifies it;
- compiler or `LogicBlock` changes.

## PR #36 implication: T2Prep is a new RF-producing path

Every RF event produced by `T2Prep` must obey the shared peak-B1 contract introduced in PR #36.

The implementation should use the same shared RF limit machinery as the existing RF/preparation
modules rather than creating a parallel check. The compiler remains the backstop, not the primary
place where a package Module discovers that it emitted an impossible RF pulse.

## Validation

### Layer 1

At minimum, independently measure from the emitted sequence:

- the semantic `prep_time_s`;
- refocusing locations and interval symmetry on the emitted lattice;
- MLEV-4 phase progression;
- pulse count / mode structure;
- spoiler placement outside the preparation interval;
- zero net gradient area inside the preparation period;
- metamorphic behaviour when `prep_time_s` changes;
- peak-B1 compliance for every RF path;
- short-preparation refusal with an actionable diagnostic.

Do not encode the stale edge-to-edge timing sentence as the invariant; fix the specification first.

### Layer 2

Build the simplest complete acquisition that makes T2Prep a real consumer rather than an isolated
waveform demonstration. Prefer an existing stable imaging path; do not introduce a new imaging
abstraction just to showcase the preparation.

### Layer 3

Use an independent Bloch-simulation path to test the claim that makes this abstraction worth
having:

```text
prepared longitudinal signal ∝ exp(-prep_time_s / T2)
```

Test over multiple preparation times and tissue T2 values, and include a control that
distinguishes T2 weighting from simple T2* decay where practical.

The initial mode does **not** claim quantitative B0/B1 robustness.

## Completion condition

`T2Prep` is shipped, documented, validated, represented in a complete acquisition, and the T2Prep
candidate/findings/current-roadmap records are updated from "approved, deliberately unbuilt" to the
implemented state.

Then move to VelocityEncode. Do not start sequence-authoring yet.

---

# Stage B — implement `VelocityEncode`

## Scope

Implement the approved **appended bipolar** velocity-encoding leaf only.

Its physical contract is self-contained:

```text
m0 = 0
Δm1 = π / (γ · VENC)
```

The public quantity is VENC; the important quantity is the **change in first moment between
toggles**, not the first moment of either toggle in isolation.

The leaf owns:

- its bipolar waveform;
- its own hardware-limited closed-form design;
- the two encoding toggles;
- its duration and moments.

The caller owns:

- which axis or axes are encoded;
- how toggles are scheduled across the acquisition;
- reconstruction / phase subtraction.

## Explicitly out of scope

Do not fold the merged flow-encoding/readout form into this leaf. The merged form changes a
neighbouring imaging waveform and belongs to the next stage.

Do not introduce a general moment-requirement IR or numerical optimizer for the canonical bipolar
case.

## Validation

Measure emitted moments independently:

- `m0 = 0`;
- toggle difference equals target `Δm1`;
- sign convention;
- VENC relation;
- origin independence once `m0 = 0`;
- hardware/raster behaviour.

Add a complete velocity-encoded acquisition or other composition that demonstrates the leaf under
real placement rather than only in isolation.

The independent moment-analysis machinery built here should be reusable as validation
infrastructure for the next stage, without becoming the definition of the module's physics.

---

# Stage C — implement flow compensation without assuming a `FlowComp` class

## Governing decision

The fine scan deliberately did **not** approve a new leaf named `FlowComp`.

Gradient moment nulling is a requirement on an imaging waveform:

> modify the waveform so that the required gradient moment of everything on that axis reaches its
> target at the semantic instant.

Therefore the implementation should follow ownership:

```text
self-contained waveform modification
        → extend the leaf that owns that waveform

cross-leaf coupled waveform
        → joint-realization kernel
```

Do not create an insertable compensation block simply to make the API look symmetric with
`VelocityEncode`.

## Implementation approach

Start from the canonical cases already supported by the evidence and the existing owners. Separate:

1. a self-contained case where one current leaf owns enough waveform to null the required moment;
   and
2. a cross-leaf case only where the physical solve genuinely needs information owned by more than
   one leaf.

For the cross-leaf case, follow the diffusion result:

```text
construct leaves / obtain construction-time physical facts
        ↓
solve the coupled physics
        ↓
materialise final events once
        ↓
return ordinary LogicBlock
```

No emitted child block should be built, taken apart, and rewritten.

## `wave-gre-flow-comp` / HarmonizedMRI

`HarmonizedMRI/wave-gre-flow-comp` is authored by the SeqCraft project owner and contains
substantial protocol-specific hard-coded flow-compensation logic.

Treat it as:

- a concrete real-sequence case;
- a source of edge cases;
- a regression/stress target;
- a useful example of what should become less hard-coded.

Do **not** treat it as an independent vote for the reusable abstraction boundary.

When this stage begins, analyse it explicitly:

```text
what is generic gradient-moment physics?
what is wave-encoding-specific?
what is GRE timing policy?
what is acquisition policy?
what is merely hard-coded implementation convenience?
```

Admit individual HarmonizedMRI repositories/files to the corpus only when they are needed and after
file-level provenance/license review. Do not bulk-admit the organization merely because it is
available.

## Completion condition

The canonical flow-compensation capabilities promised by the post-v0 queue are expressible through
the correct owners, with composition-level moment validation. The result may be extensions plus a
joint-realization kernel; it does not need to produce a public `FlowComp` class.

---

# Stage D — implement `bSSFPTR`

## Scope

Implement the approved RF-to-RF balanced-SSFP repetition kernel.

The kernel owns the per-repetition physical conditions it can state and verify:

- zero net gradient area on every axis over the balanced RF-to-RF interval;
- declared TE placement;
- RF/ADC phase convention/progression;
- coupled timing inside the repetition.

It does **not** own the claim that magnetisation has reached steady state. Startup/catalyzation and
repetition ordering remain composition/acquisition policy.

## The remaining 220 µs question

The fine scan found a measurable cost when a rewinder and the next repetition's prephaser are kept
as separate physical pieces.

Before encoding an API, settle which of these is true:

```text
A. the efficient waveform changes because of bSSFP physical requirements
   → physics-aware joint realization in the Module/kernel layer

B. two already-decided waveforms can be fused with identical physical semantics
   → candidate for a general semantics-preserving compiler optimisation
```

bSSFP physics must not be smuggled into the compiler.

A compiler optimisation, if eventually justified, should be general and semantics-preserving.
Correct bSSFP should not depend on candidate-specific compiler logic.

## Validation

Validate the balanced condition on the whole RF-to-RF composition, not by asking children whether
they are balanced. Include a complete acquisition and an independent physical/simulation check
appropriate to bSSFP contrast and phase cycling.

---

# Stage E — establish `sequence-authoring` skill v0

Only begin this stage once the four approved post-v0 capabilities above have been implemented or
have reached an explicit human-reviewed final disposition.

At that point SeqCraft has a much more stable package surface on which to teach an authoring agent.

## Responsibility split

`module-mining` asks:

> Is there a reusable physical abstraction that belongs in the SeqCraft package?

`sequence-authoring` asks:

> Given a requested MRI protocol, how do I correctly build the complete sequence using SeqCraft?

These are different optimization targets.

```text
module-mining
    optimize for reusable abstraction quality

sequence-authoring
    optimize for correct complete sequence construction
```

A successful sequence does not imply a new Module. `NO_NEW_MODULE` does not imply authoring
failure.

## Authoring workflow

```text
user protocol / sequence requirement
        ↓
identify acquisition family and physical requirements
        ↓
inspect shipped SeqCraft Modules + relevant examples
        ↓
design composition
        ↓
is a reusable abstraction apparently missing?
        │
        ├─ no → build the sequence
        │
        └─ yes → emit ModuleGap
                    ↓
                 module-mining
                    ↓
         independent reusable-boundary decision
                    ↓
      package Module / NO_NEW_MODULE / pending
                    ↓
                continue authoring
        ↓
compile
        ↓
physics / timing / hardware validation
        ↓
script / notebook / .seq
```

## Fallback hierarchy

The authoring skill should prefer:

```text
1. shipped SeqCraft Module
        ↓
2. file-domain local Module
        ↓
3. direct LogicBlock / raw event assembly
```

Use a local Module when the current sequence genuinely repeats a semantic component but evidence is
insufficient to put it in the package.

Use direct `LogicBlock` assembly for one-off glue or simple placement. Do not manufacture classes
merely to make generated code look modular.

## ModuleGap handoff

`sequence-authoring` may report a reusable-abstraction hypothesis, but it must not create the final
module-mining `Candidate` or command mining to produce a Module.

A thin handoff is sufficient:

```text
ModuleGap

requested capability:
why existing shipped modules are insufficient:
can current primitives express it:
candidate reusable physical boundary:
physics that appears reusable:
current protocol needs:
fallback if no package module:
evidence worth inspecting:
```

`module-mining` then evaluates the evidence independently. `NO_NEW_MODULE` remains a successful
outcome.

## Initial calibration prompts

Calibrate with at least three classes:

```text
A. fully covered
   should compose shipped modules without inventing abstractions

B. locally composable
   should choose a local Module or direct LogicBlock glue without package pollution

C. genuine reusable gap
   should produce a ModuleGap and continue with a fallback while module-mining evaluates it
```

Use several real complete acquisitions, including at least one that exercises preparation, one
motion/moment case, and one kernel-heavy sequence.

---

# Corpus policy after this milestone

Do not turn corpus expansion into a goal of its own.

HarmonizedMRI and other repositories should be admitted when they answer a concrete authoring or
mining question. Evaluate individual repositories/files for:

- licence and provenance;
- independence from existing witnesses;
- role: discovery / design-witness / oracle / concrete case;
- whether hard-coded sequence policy is being mistaken for reusable physics.

Once `sequence-authoring` exists, new module mining should become primarily **demand-driven**:

> What reusable physical gap does the authoring agent repeatedly encounter in real protocols?

rather than:

> What is the next named MRI sequence we can turn into a class?

---

# Definition of the milestone

The transition to `sequence-authoring` is ready when:

- `T2Prep` is shipped for the approved initial mode;
- `VelocityEncode` is shipped for the approved appended bipolar form;
- canonical flow-compensation behaviour is implemented through the correct owners, without forcing
  a `FlowComp` leaf;
- `bSSFPTR` is shipped with the balance/phase/timing contract validated on the composition;
- each capability has appropriate Layer 1/2/3 evidence;
- compiler and `LogicBlock` remain candidate-neutral;
- no unresolved post-v0 architecture question blocks normal sequence composition.

At that point the primary project loop becomes:

```text
stable reusable Module library
        ↓
sequence-authoring
        ↓
real protocol gaps
        ↓
ModuleGap
        ↓
module-mining
        ↓
library evolves only when evidence justifies it
```

That is the intended handoff from **library construction** to an **agentic sequence-programming
system**.
