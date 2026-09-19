# SeqCraft Module Mining — v0 Close-out and Next-Phase Roadmap

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../README.md) for which one to edit.

**Status:** v0 close-out active  
**Current scope:** finish and merge the current integrated PR cleanly  
**Important:** only the **v0 Close-out** section is actionable now. The later sections are intentionally stored in-repo as the roadmap for the next phase.

---

## 1. Why this document exists

The original post-TSE phase set a narrow success condition:

- `TSE/FSE = GREEN`
- `RadialReadout = GREEN`
- compare the two pilots before hardening the workflow into a durable AI-assisted module-mining process.

That condition has effectively been met, and the work has gone further:

- `TSEShot` / `FSE2D` are GREEN;
- `RadialReadout` is GREEN;
- the two-pilot retrospective and comparator known-answer tests exist;
- `GRE3DTR` has been mined as a third candidate;
- end-to-end build / simulate / reconstruct validation has been restored for FSE/GRE3D.

The remaining work is therefore **not another round of broad manual module extraction**. It is to close v0 cleanly, capture the lessons that were only visible after GRE3D, and then use those rules to build the next module-mining workflow.

---

# 2. v0 Close-out — actionable now

## 2.1 Treat physical modes as independent contracts

The first GRE3D implementation exposed an important failure mode.

The incorrect mental model was effectively:

```text
slab-selective
    = shaped RF + Gz

non-selective
    = the same thing with Gz removed
```

That produced a spatially non-selective sequence, but it still used the slab path's 3 ms sinc RF. The emitted sequence therefore did **not** match the intended conventional non-selective GRE3D protocol.

The correct rule for future mining is:

> **Alternative physical modes must be re-derived from their own physical contract. They must not be implemented by subtracting events from another mode unless independent reference evidence explicitly supports that equivalence.**

For every mode or major branch, the miner must state what the scanner actually plays:

- RF family / shape;
- RF duration and source of the default;
- spatially selective vs non-selective behavior;
- gradients that must exist;
- gradients that must not exist;
- intrinsic moment / rephasing requirements;
- ADC / readout semantics;
- timing consequences;
- reference evidence supporting those choices.

A mode name alone is not enough.

---

## 2.2 Add a Mode Contract Table before implementation

Any candidate with physical modes, variants, or defaults should have a table before code is written.

Example fields:

| Field | Mode A | Mode B |
|---|---|---|
| RF family | | |
| RF duration/default source | | |
| selection gradient | | |
| intrinsic rephasing requirement | | |
| encoding gradient role | | |
| ADC/readout semantics | | |
| timing consequence | | |
| reference evidence | | |
| intentionally user-overridable behavior | | |

The table is not an implementation plan. It is the physical contract the implementation must later realize.

The miner should also record a **differential contract**:

```text
mode A → mode B

changes:
    ...

must remain invariant:
    ...
```

This prevents "mode B is mode A minus one event" reasoning from becoming an implicit design rule.

---

## 2.3 Inspect the emitted sequence, not only derived measurements

A candidate can satisfy k-space and timing checks while still emitting the wrong physical pulse family.

For each promoted mode, validation must include at least one direct inspection of the compiled/emitted sequence:

```text
RF waveform/family
RF duration
selection-gradient presence/absence
rephaser realization/ownership
ADC/readout event
spoiler/rewinder structure where relevant
```

The question is:

> **If the emitted `.seq` were inspected without reading the Python, would its RF, gradients, ADC, and timing tell the same physical story as the API name and documentation?**

This is distinct from:

- analytic invariants;
- comparator output;
- simulation/reconstruction;
- notebook smoke execution.

All remain useful, but none replaces emitted-sequence inspection.

---

## 2.4 Reference comparisons must declare claim scope

"Validated against reference X" is too broad unless the comparison actually covers the whole protocol.

Every reference comparison should state what it proves.

Example:

```text
validated against official GRE3D reference for:
    k-space lattice
    semantic centre
    partition spacing
    requested encoding

not used to establish:
    slab-selective RF implementation
    alternative waveform decomposition
    signal equality when TE differs
```

A comparator result must not silently expand into a stronger physical claim than it measured.

---

## 2.5 Shared-leaf changes require a dependency-impact review

Not every Module change is isolated.

Two examples from v0:

```text
Excitation
    new semantic rephasing requirement
    build(rephase=False)

default behavior remains unchanged
→ existing consumers remain behaviorally stable unless they opt in
```

versus:

```text
CartesianLine
    ADC sample-count legality tightened

all callers inherit the stronger leaf contract
→ previously invalid protocols may now fail earlier
```

Before merging a shared-leaf change, explicitly map:

```text
changed leaf
    ↓
direct Module consumers
    ↓
transitive consumers
    ↓
affected examples / notebooks
```

For every affected consumer, classify the result:

```text
NO_BEHAVIOR_CHANGE
NEW_EARLY_REFUSAL_FOR_PREVIOUSLY_INVALID_INPUT
INTENTIONAL_BEHAVIOR_CHANGE
UNKNOWN / NEEDS_REVIEW
```

Passing the full suite is necessary but not sufficient; the impact should be stated explicitly.

---

## 2.6 Keep the validation ladder

The v0 validation ladder should remain:

```text
Layer 1
analytic invariants + regression tests
    → blocking

Layer 2
example notebook executes
    → blocking smoke

Layer 3
simulate + reconstruct + inspect
    → review evidence
```

Add the direct emitted-sequence inspection described above as a separate review activity rather than forcing it into one of those labels.

Each layer catches a different failure class.

---

## 2.7 Close the Radial example coverage gap

Current new/prompted Module coverage is approximately:

```text
TSEShot
    covered through FSE examples

FSE2D
    build + simulate/reconstruct + module API example

GRE3DTR
    build + simulate/reconstruct

RadialReadout
    package tests and mining evidence
    but no shipped user-facing example yet
```

Before v0 is closed, add a lightweight Radial example if practical.

Preferred minimum:

```text
examples/radial_gre/
    01_build.ipynb
```

It should compose:

```text
Excitation
+ RadialReadout
+ ordinary Python angle/spoke scheduling
```

and make visible:

- spoke orientation;
- `k=0` location;
- full-spoke / centre-out behavior when supported;
- emitted trajectory geometry;
- acquisition policy remaining outside the leaf.

Do **not** create a `RadialGRE` imaging class only to make the example convenient.

If the existing non-Cartesian simulation/reconstruction stack can be reused cheaply, add:

```text
02_simulate_and_reconstruct.ipynb
```

Otherwise do not create a new reconstruction architecture merely to satisfy the example checklist.

---

## 2.8 Final v0 close-out checklist

Before merging the current integrated PR:

```text
[ ] GRE3D non-selective default is the intended short hard/block RF
[ ] slab-selective GRE3D remains shaped/selective
[ ] explicit shaped-but-spatially-nonselective RF remains possible
[ ] GRE3D emitted `.seq` modes are inspected directly
[ ] physical-mode rules are written into the fine-scan / mining guidance
[ ] shared-leaf dependency-impact rule is documented
[ ] reference comparisons state claim scope
[ ] Radial example gap is closed, or explicitly deferred with rationale
[ ] full tests / doctests / lint / mypy / API checks pass from final combined tree
[ ] smoke notebooks pass from final combined tree
[ ] simulation notebooks are re-run where applicable
[ ] workspace/external planning docs are synchronized to the final v0 state
[ ] current PR is reviewed as the v0 close-out/integration PR
```

After this checklist, v0 should be considered closed.

---

# 3. What should *not* block v0 close-out

The following candidates were in the old recommended order, but they are **not required** to finish v0:

```text
SaturationPrep
DiffusionPrep / DiffusionEncoding
SpiralReadout / Spiral2D
```

They now serve a better role as future workflow-calibration and stress-test candidates.

Do not manually finish all of them before beginning skill formalization.

---

# 4. Next phase — Module Mining Skill v0

**Do not execute this section until the current v0 close-out is merged.**

The next goal is not immediate large-scale batch conversion.

First formalize a small module-mining workflow that preserves the successful parts of the pilots:

```text
candidate discovery
→ evidence collection
→ physical contract
→ ownership analysis
→ reference claim scope
→ invariant table
→ adapter / measurements
→ extraction
→ emitted-sequence inspection
→ parameter sweep / metamorphic tests
→ example/end-to-end validation where meaningful
→ GREEN / YELLOW / RED
→ human review
```

The skill should explicitly contain:

- physical-mode contract tables;
- dependency-impact maps for shared leaves;
- independent reference evidence;
- claim-scoped comparison;
- emitted `.seq` inspection;
- validator known-answer tests;
- no compiler change without an independent minimal reproducer;
- no automatic promotion when the physical boundary remains ambiguous.

The goal of this phase is a **supervised skill**, not autonomous batch merging.

---

# 5. Supervised calibration candidate: SaturationPrep

Saturation is useful because it tests a new family without extreme trajectory complexity.

Questions it should force the skill to answer:

```text
Is this preparation physics or excitation physics?
Does it legitimately reuse Excitation?
What does RF `use` mean here?
Who owns crushers/spoilers?
What timing belongs to the prep itself?
```

This candidate is intentionally comparatively simple.

Its purpose is to detect over-design:

> Can the skill extract one coherent physical concept without inventing unnecessary hierarchy?

---

# 6. Supervised calibration candidate: Spiral

Spiral should be the harder calibration case.

There is already an unmerged historical PR:

```text
PR #23
"Add Spiral2D, and the two tolerances that refused a legal spiral"
```

Do **not** treat that PR as the answer, and do **not** start from a blank page.

Treat it as an evidence corpus.

Useful material already present there includes:

- path vs traversal separation;
- variable-density trajectory design;
- spiral-out / spiral-in / in-out / out-in variants;
- `k_per_m` measured from emitted events;
- GRE spiral and SE spiral examples;
- non-Cartesian reconstruction;
- vector-gradient limit handling;
- proposed compiler-tolerance fixes.

The new mining workflow should independently ask:

```text
What is the reusable unit?
    SpiralReadout?
    SpiralInterleaf?
    Spiral2D?

Does one leaf own both path and traversal?

Who owns shot/interleaf rotation?

Are in/out/in-out/out-in genuinely one physical contract?

Is the proposed density vocabulary physically clear?

Do Radial + Spiral now justify a trajectory utility layer?

Are the compiler changes real compiler defects,
or is the candidate asking the compiler to rescue a design problem?
```

Then compare the new result with PR #23.

PR #23 is therefore:

```text
not the specification
not disposable work
a rich positive/negative evidence set
```

---

# 7. Compiler-change rule for future candidates

A candidate may expose a real compiler bug, but the burden of proof is high.

If a candidate claims the compiler is wrong:

```text
candidate failure
    ↓
remove the candidate
    ↓
construct minimal raw-PyPulseq / raw-LogicBlock reproducer
    ↓
show a physically/legal input is rejected or altered incorrectly
    ↓
characterize the numerical/error floor independently
    ↓
only then modify compiler behavior
```

This is especially relevant to the Spiral PR, which proposed changes to m0/m1 verification tolerances.

The compiler remains architecturally frozen unless such an independent reproducer demonstrates a genuine compiler defect.

---

# 8. After supervised calibration — batch scan

Only after at least the simple and hard supervised candidates have been run through the skill should broad batch scanning begin.

Batch scan should focus first on:

```text
inventory
candidate discovery
deduplication against existing modules
reference availability
reuse/coverage value
physical-boundary confidence
```

Batch scan should be allowed to conclude:

```text
NO_NEW_MODULE
YELLOW — NEEDS_HUMAN_PHYSICS_REVIEW
REFERENCE_INSUFFICIENT
WRAPPER_ONLY
DUPLICATES_EXISTING_MODULE
```

It should not be optimized for maximizing Module count.

---

# 9. Diffusion as a later architecture stress test

`DiffusionPrep / DiffusionEncoding` should remain later than Saturation and Spiral.

It combines:

- b-value / b-tensor semantics;
- refocusing context;
- gradient cross terms;
- timing coupling;
- direction tables;
- possible preparation-vs-encoding ownership ambiguity.

That makes it a useful stress test **after** the workflow is stable.

It is not a prerequisite for v0 close-out or initial skill formalization.

---

# 10. Proposed phase sequence

```text
CURRENT
v0 close-out
    ├── physical-mode guardrail
    ├── shared-leaf dependency impact
    ├── reference claim scope
    ├── emitted-sequence inspection
    ├── Radial example close-out
    └── merge current integrated PR

NEXT
Module Mining Skill v0
    └── supervised, human-gated

CALIBRATION 1
SaturationPrep
    └── simple new family / anti-overdesign case

CALIBRATION 2
Spiral
    └── use PR #23 as evidence corpus
        └── independently re-derive boundary and compiler claims

RETROSPECTIVE
compare supervised skill runs

THEN
batch scan

LATER
DiffusionPrep / DiffusionEncoding
    └── architecture stress test
```

---

# 11. Governing principle

The module-mining workflow should optimize for one thing:

> **Extract a reusable physical contract from evidence, then prove the implementation realizes that contract.**

Not:

```text
find similar code
→ wrap it
→ make tests pass
```

The GRE3D RF mistake is the clearest v0 reminder of why this distinction matters.
