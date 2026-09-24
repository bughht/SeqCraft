# Repetition-Level Joint Waveform Design

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-24, so that the pull request carries the direction behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../../README.md) for which one to edit.

**Status:** architecture-revisit reference; design direction, not a frozen public API  
**Project:** SeqCraft  
**Date:** 2026-09-24  
**Trigger:** `VelocityEncode` is shipped, the first `CartesianLine` flow-compensation realization is implemented in PR #39, and the earlier `repetition_augmentation_design_note.md` explicitly said to revisit shared repetition augmentation once two real consumers exist.

---

## 1. Decision summary

SeqCraft should **not** solve reusable velocity encoding and flow compensation by teaching every sequence kernel a separate implementation of each feature.

That creates an `N sequences × M augmentations` maintenance problem:

```text
GRE2DTR  × VelocityEncode
GRE2DTR  × FlowComp
GRE3DTR  × VelocityEncode
GRE3DTR  × FlowComp
EPI      × VelocityEncode
EPI      × FlowComp
...
```

The target architecture is instead:

```text
base repetition requirements
        +
augmentation requirements
        ↓
one repetition-level physical design problem
        ↓
shared joint waveform designer
        ↓
analytic realization where possible
optimizer-backed realization only when necessary
        ↓
ordinary final gradient events
        ↓
LogicBlock
```

The important reusable abstraction is therefore **requirements-before-realization**, not a universal numerical optimizer.

PR #39's `CartesianLine(null_moment_order=1)` remains useful as:

- a correct local canonical realization;
- a regression/reference case;
- potentially one realization primitive used by the shared designer.

It should **not yet be treated as the final Stage C architecture**. In particular, the public placement of `null_moment_order` on `CartesianLine` should remain reviewable until the repetition-level design is tested.

---

## 2. Why revisit now

The earlier non-binding design note deliberately deferred a shared repetition-augmentation mechanism until at least two real consumers existed.

That threshold is now crossed:

1. `VelocityEncode`
   - shipped;
   - carries reusable velocity-encoding physics;
   - currently has a standalone bipolar realization;
   - should not require each sequence family to reimplement VENC semantics.

2. Flow compensation
   - PR #39 demonstrates a local readout-axis realization;
   - the actual physical requirement is composition-level: total gradient moment at a semantic echo;
   - `wave-gre-flow-comp` is a concrete stress case where hard-coded per-sequence logic should become less hard-coded.

3. `DiffusionSEPrep`
   - adjacent evidence that joint physical design can live in the Module/kernel layer before `LogicBlock` materialization;
   - not necessarily a direct consumer of the same first implementation, but evidence that the architectural layer is viable.

This is no longer abstraction from a single hypothetical case.

---

## 3. The problem to solve

### 3.1 What must be reusable

A reusable augmentation should define its **physical requirement once**.

For velocity encoding, conceptually:

```text
axis = x
state = + / -
at semantic echo:
    M0 target = imaging requirement
    M1 target = +delta_m1/2 or -delta_m1/2
delta_m1 = 1 / (2 * venc)
```

For first-order flow compensation:

```text
at semantic echo:
    M1_total(axis) = 0
```

The augmentation should not need to know whether the base acquisition is:

- 2D GRE;
- 3D GRE;
- EPI;
- a future Cartesian kernel;
- another compatible repetition family.

Likewise, a new repetition kernel should not reimplement:

- the VENC formula;
- velocity polarity semantics;
- the definition of first-moment nulling;
- every future augmentation one by one.

### 3.2 What remains sequence-specific

A repetition kernel still owns the base sequence physics:

- RF and echo semantic instants;
- readout geometry;
- phase/partition encoding;
- which gradient regions are fixed;
- which gradient regions are free to redesign;
- TE/TR policy;
- line/partition/state dependence;
- shared timing across the acquisition.

The goal is not to make all sequences identical. The goal is to make augmentations act on a **shared physical design boundary** rather than on sequence-specific handwritten branches.

---

## 4. Architectural principle

Keep three questions separate:

```text
WHAT must be true?
    physical requirements at semantic instants

WHO owns the coupled design?
    repetition-level Module/kernel design

HOW is it realized?
    analytic construction first
    numerical optimization only when genuinely needed
```

Do not let the realization backend define public semantics.

The compiler remains unchanged in responsibility:

```text
Module/kernel physical design
        ↓
LogicBlock timed-event tree
        ↓
Compiler
        ↓
Pulseq legality / block emission
```

The compiler must not learn about:

- VENC;
- flow compensation;
- `M1` targets;
- RF isodelay physics;
- b-value;
- pathway physics.

---

## 5. Proposed internal model

The names below are provisional. The shape matters more than the names.

### 5.1 Repetition state

A build call supplies acquisition state:

```text
RepetitionState
    line / ky
    partition / kz
    velocity state
    echo selection where relevant
    other sequence-specific indices
```

This allows the physical problem to vary per repetition without duplicating augmentation logic.

Examples:

```text
GRE2D:
    line changes

GRE3D:
    line + partition change

phase-contrast GRE:
    line + partition + velocity state change
```

### 5.2 Semantic instants

The base repetition exposes named physical instants to the internal design problem:

```text
rf_center
rf_isodelay
echo[0]
echo[1]
...
refocusing_center[n]     when relevant
```

These are Module-layer semantics, not anchors added to `LogicBlock`.

### 5.3 Moment requirements

A minimal requirement is conceptually:

```text
MomentTarget
    axis
    order
    semantic endpoint
    target value
```

Examples:

```text
MomentTarget(x, order=0, at=echo[0], value=0)
MomentTarget(x, order=1, at=echo[0], value=0)

MomentTarget(x, order=1, at=echo[0], value=+delta_m1/2)
MomentTarget(x, order=1, at=echo[0], value=-delta_m1/2)
```

Do not generalize this immediately into an arbitrary symbolic constraint language.

### 5.4 Base waveform facts and degrees of freedom

The repetition design must know two different things.

**Fixed contributions** are physical waveform portions that are already determined by the base sequence.

**Adjustable regions** are waveform portions the joint designer is allowed to realize or reshape.

Conceptually:

```text
axis x:
    fixed readout lobe
    adjustable pre-echo winder window

axis y:
    line-dependent phase-encode requirement
    adjustable encoding / rewinding window

axis z:
    slab-selection contribution
    partition-dependent encoding requirement
    adjustable shared winder window
```

Augmentations should not receive private gradient event handles.

The base kernel supplies the design problem; the shared designer owns the actual realization of adjustable regions.

### 5.5 Hardware and timing

The problem also carries:

```text
max_grad
max_slew
gradient raster
available windows
required endpoints
shared timing constraints
```

For a full acquisition, the kernel may enumerate all supported states and choose a common timing envelope so TE/TR do not accidentally vary with line or partition.

`GRE3DTR` already provides a concrete precedent: it enumerates partition-dependent signed requirements, finds the limiting case, and uses one shared winder duration.

---

## 6. Augmentation interface

A reusable augmentation contributes **requirements**, not built events to be patched into a completed repetition.

Conceptually:

```text
augmentation
    ↓
reads semantic physical context
    ↓
adds requirements
```

An augmentation must not:

- inspect a built `LogicBlock` and rewrite events;
- reach into another module's private gradient object;
- contain sequence-specific branches such as `if isinstance(base, GRE3DTR)`;
- know the compiler.

### Velocity encoding

The existing `VelocityEncode` physics should remain single-source.

Its standalone bipolar `build(polarity=...)` can remain a canonical direct realization.

For repetition-level use, the same physical object should be able to contribute the corresponding moment target without each kernel reimplementing the VENC formula.

The exact API is not frozen. Possible shapes include an internal requirement-provider method or a wrapper/adaptor around the existing object.

### Flow compensation

The fine scan correctly rejected an **insertable `FlowComp` block**.

That does not prohibit a reusable non-block requirement object if one becomes useful.

For example, a future public or internal object could mean:

```text
first moment on selected axes shall be zero at the selected semantic echo
```

without itself emitting events.

The naming and public status should be decided only after the architecture spike.

---

## 7. Joint waveform designer

The shared designer consumes:

```text
base physical problem
+
augmentation requirements
+
repetition state
```

and returns the final physical waveform design for the adjustable regions.

### 7.1 Analytic backend first

Canonical cases should stay analytic when possible.

Examples already seen:

- appended bipolar velocity encoding;
- `CartesianLine` equal-duration two-winder first-moment nulling;
- slab rephase + partition encoding in `GRE3DTR`;
- published canonical flow-compensation constructions.

### 7.2 Numerical optimizer is a backend, not the abstraction

A numerical backend becomes justified when the actual problem requires it, for example:

- multiple moment orders;
- multiple semantic echoes;
- fixed waveform segments;
- non-zero start/end gradients;
- arbitrary existing waveforms;
- minimum-TE optimization with several free durations;
- eddy-current constraints;
- PNS constraints;
- concomitant-field constraints.

The architecture should permit this backend later without changing augmentation semantics.

Do **not** implement a general optimizer merely to satisfy this design document.

---

## 8. Why this avoids N × M duplication

Without the shared boundary:

```text
each sequence × each augmentation
    → handwritten integration logic
```

With the shared boundary:

```text
new repetition family
    → implement one base physical-problem adapter

new augmentation
    → implement one reusable requirement contributor

shared designer
    → combines compatible pairs
```

The maintenance shape becomes approximately:

```text
N sequence adapters + M augmentation contributors
```

instead of:

```text
N × M integrations
```

Compatibility is still explicit. Not every augmentation must work on every sequence.

A sequence can refuse an augmentation when it does not expose the required semantic instant or adjustable waveform freedom.

That refusal must be based on capability, not on a hard-coded list of class names.

---

## 9. Multi-state and acquisition-wide timing

A repetition-level designer does not need to own acquisition ordering, but it must accept the state supplied by acquisition policy.

For GRE:

```text
state = (ky, kz, velocity polarity, ...)
```

The imaging layer decides which states to acquire and in what order.

The repetition/kernel layer solves the waveform for a given state.

When a timing quantity must remain invariant across the acquisition:

```text
enumerate / bound supported states
        ↓
find limiting realization
        ↓
choose common timing window
        ↓
realize every state inside that window
```

This is the same pattern already used by `GRE3DTR` for partition-dependent z winders.

`wave-gre-flow-comp` should be used as a stress/regression case for this state-dependent design.

It is not independent evidence for the abstraction boundary.

---

## 10. Multi-echo

The first architecture should be able to represent multiple semantic checkpoints even if the first implementation only solves one.

For example:

```text
echo[0]: M1(x) = 0
echo[1]: M1(x) = 0
echo[2]: M1(x) = 0
```

or velocity-encoding targets per echo.

This is a strong candidate for a future optimizer-backed realization because the constraints span one train and may not admit the simple closed forms used by PR #39.

Do not make multi-echo support a blocker for the first repetition-level architecture.

---

## 11. Pathway-aware moments are deferred

This design concerns gradient-moment requirements for the intended sequence pathway and semantic intervals.

It does **not** introduce a general coherence-pathway model.

If pathway-aware analysis is later added, proceed in this order:

```text
first:
    pathway-aware zeroth-moment / k-space bookkeeping (M0)

then:
    pathway-aware first moment (M1)

later:
    higher moments as justified
```

Incomplete spoiling, stimulated echoes, and residual coherence across TRs therefore remain out of scope for this architecture revisit.

Do not build an EPG/pathway framework as part of Stage C.

---

## 12. PR #39 disposition

PR #39 currently establishes a useful local result:

```text
CartesianLine
    first echo
    readout axis
    equal-duration two-winder realization
    M0 = 0
    M1 = 0
```

The implementation and Layer-1/Layer-2 evidence should be preserved.

However, before merge, decide whether the public API:

```python
CartesianLine(..., null_moment_order=1)
```

is:

1. a useful standalone local option that should remain public;
2. an internal realization primitive used by the repetition-level designer;
3. or a temporary API that should be refactored before shipping.

Do not merge PR #39 under the claim that the overall Stage C architecture is complete until this question is answered.

The flow-compensation record may be GREEN for the **local readout-axis realization**, but the roadmap should continue to mark the repetition-augmentation architecture as unresolved until the shared design is demonstrated.

---

## 13. Architecture spike: required demonstrations

Before freezing an API, build the smallest prototype that proves or disproves the shared boundary.

### Demonstration A — same velocity-encoding semantics, two repetition families

Use the same `VelocityEncode` physical object/specification with:

- `GRE2DTR`;
- `GRE3DTR`.

Success means:

- the VENC formula is implemented once;
- polarity semantics are implemented once;
- neither kernel contains a `VelocityEncode`-specific physical derivation;
- whole-composition moment targets are independently validated.

### Demonstration B — same flow-compensation requirement across varying repetition state

Use a reusable first-moment-null requirement with at least:

- GRE2D line variation;
- GRE3D line/partition variation or the relevant `wave-gre-flow-comp` stress geometry.

Success means:

- compensation is not a `flow_comp=True` branch duplicated in every gradient leaf;
- line/partition-dependent existing moments are included automatically in the repetition problem;
- final whole-repetition moments meet the target at TE;
- timing cost is explicit;
- common timing is preserved where the acquisition requires it.

### Demonstration C — local realization survives as a backend/reference

Show whether PR #39's `CartesianLine` construction can be reused inside the shared design or remains independently valuable.

Do not preserve its public API merely to avoid changing the PR.

---

## 14. Failure conditions for the proposed architecture

Treat the spike as failed or incomplete if any of these are necessary:

- every sequence kernel branches on every augmentation type;
- every new augmentation requires editing all existing kernels;
- an augmentation rewrites an already-built `LogicBlock`;
- augmentations receive private gradient event handles from leaves;
- compiler or `LogicBlock` semantics are changed to carry MRI moment requirements;
- the shared layer is just a generic plugin bus with no physical contract;
- a numerical optimizer is introduced before a concrete problem needs it;
- correctness is checked only from solver outputs rather than from emitted waveform moments.

---

## 15. Validation rules

For every realized repetition:

```text
independent validator
    ↓
integrate final emitted waveform
    ↓
check moment targets at semantic instants
```

The validator must not ask the augmentation or solver for its own reported `M1`.

For state-dependent acquisitions, validate representative and limiting states, and enumerate all states when practical.

Keep separate:

```text
feasible under declared constraints
```

from:

```text
globally minimum TE
```

The latter needs an actual optimality argument.

---

## 16. Recommended near-term sequence

```text
PR #39
    preserve implementation/evidence
    keep open and unmerged during architecture revisit

architecture spike
    define minimal repetition physical-problem boundary
    test VelocityEncode + FlowComp across more than one GRE kernel
    use wave-gre-flow-comp as stress/regression
    no pathway framework
    no generic optimizer yet

decision
    keep/refactor PR #39 public API
    freeze the minimal shared boundary
    update roadmap + repetition_augmentation_design_note

then
    complete Stage C on repetition-level composition

then
    proceed to bSSFPTR
```

---

## 17. Core invariant to preserve

The architecture is successful only if all of these remain true:

```text
Module/kernel layer owns MRI physical design.
LogicBlock remains a thin timed-event tree.
Compiler owns Pulseq legality only.

Velocity encoding is defined once.
Flow-compensation requirement is defined once.
A new compatible sequence does not reimplement them.

Base requirements + augmentation requirements
are solved together before waveform materialization.

The final emitted composition, not the design variables,
is the source of truth for validation.
```
