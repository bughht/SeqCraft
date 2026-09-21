# SeqCraft Module Mining — v0 Close-out and Next-Phase Roadmap

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-20, so that the pull request carries the evidence behind the modules it adds.
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
examples/gre_radial_2d/
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

### 2.7.1 How this resolved, and what it does **not** claim

`examples/gre_radial_2d/01_build.ipynb` ships. `02_simulate_and_reconstruct.ipynb` is **deferred**, and
the deferral is a decision rather than an unfinished v0 task.

```text
RadialReadout

Layer 1 -- physical/reference validation      GREEN
Layer 2 -- complete acquisition build         GREEN
Layer 3 -- simulate + reconstruct             DEFERRED
```

`RadialReadout`'s GREEN claim is, exactly:

> the radial spoke leaf emits the correct waveform / ADC / trajectory semantics, independently
> validated against the official reference.

It is **not** the stronger claim:

> SeqCraft has an independently validated end-to-end radial imaging + non-Cartesian reconstruction
> pipeline.

Nothing in the repository has established the second one, and no wording anywhere should let a
reader slide from the first to the second. This is §2.4's rule applied to our own status table.

**Why Layer 3 is not cheap here.** A radial Layer 3 needs a non-Cartesian reconstruction path, and
that path carries its own conventions -- NUFFT coordinate scaling, density compensation, weighting,
trajectory timing -- each of which needs independent review before any image it produces can be
evidence for anything. An image that looks right through an unreviewed reconstruction is not a
validation of the sequence; it is two unvalidated things agreeing.

PR #23 contains useful non-Cartesian reconstruction infrastructure. Under this roadmap it is
**evidence for the Spiral supervised-calibration phase, not accepted architecture**, and it must
not be pulled into a radial example to complete a `01`/`02` pair.

### 2.7.2 The deferred item, and where it is picked up

During Spiral supervised calibration (§6):

```text
1. independently review PR #23's non-Cartesian reconstruction path;
2. determine what is genuinely trajectory-agnostic;
3. use Radial + Spiral as two real consumers of that contract;
4. if the shared reconstruction utility survives review,
   add gre_radial_2d/02_simulate_and_reconstruct.ipynb;
5. use the resulting end-to-end path as Layer-3 evidence
   for both families where appropriate.
```

Step 3 is the load-bearing one and is the reason this waits for Spiral rather than being done now:
a reconstruction contract with one consumer is that consumer's implementation with a more general
name. Radial alone cannot show which parts are trajectory-agnostic. Two real consumers can.

**Step 4 is a conditional, not a scheduled deliverable.** The notebook is added *only if* a
genuinely shared reconstruction utility survives review. If the review finds that PR #23's path is
spiral-shaped rather than trajectory-agnostic, the right outcome is that the radial `02` is never
written and Radial stays Layer-2 GREEN -- not that a radial-specific reconstruction gets built to
close the item. An item that can only be discharged one way is a deadline, not a decision.

Radial Layer 3 is therefore **not blocked on effort**; it is sequenced behind a contract that does
not exist yet and may not turn out to exist at all.

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

### 6.1 Radial Layer 3 is carried here

§2.7.2 defers `RadialReadout`'s simulate-and-reconstruct layer into this phase, and it is listed
here so the phase cannot start without it:

```text
1. independently review PR #23's non-Cartesian reconstruction path;
2. determine what is genuinely trajectory-agnostic;
3. use Radial + Spiral as two real consumers of that contract;
4. if the shared reconstruction utility survives review,
   add gre_radial_2d/02_simulate_and_reconstruct.ipynb;
5. use the resulting end-to-end path as Layer-3 evidence
   for both families where appropriate.
```

Radial is the **second consumer** that makes step 2 answerable, and it is a cheap one: the
trajectory is a straight line, so any coordinate, density-compensation or timing convention the
contract gets wrong shows up in a geometry simple enough to reason about by hand. Reviewing the
reconstruction path against radial first, then spiral, is the ordering that makes a failure
attributable.

Until that review happens, PR #23's reconstruction is **evidence, not accepted architecture**, and
nothing in `examples/` may depend on it. Steps 4 and 5 are conditional on the review's outcome:
no shared utility, no radial `02`, and Radial remains Layer-2 GREEN with the claim §2.7.1 states.

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

# 9b. Corrections from the MRzero coarse scan

Recorded here because they change the scan model, not any Module. Full reasoning in
`tools/module_mining/skill_v0/coarse_scan_mrzero.md`.

**Coverage is classified, not called "duplicate".** "SeqCraft can already build it" is not
"SeqCraft already ships the right reusable abstraction":

```text
DIRECT_SHIPPED_DUPLICATE            a public shipped abstraction already exists
DEGENERATE_CASE_OF_EXISTING_MODULE  a parameter limit of a shipped abstraction
COMPOSITION_COVERED                 modules express it; no equivalent public abstraction
NOTEBOOK_ONLY_EXISTING              the abstraction lives only in an example notebook
PRIMITIVE_COMPOSITION_COVERED       plain composition IS the intended API
```

**Spin echo is not a shipped duplicate.** SeqCraft ships no `SE2D`; the only one is notebook-local
in `examples/se_2d/`. `FSE2D(echoes=1)` is a conventional Cartesian 2D spin echo, so the status is
`DEGENERATE_CASE_OF_EXISTING_MODULE` + `ARCHITECTURE_REVISIT_CANDIDATE`. **No `SE2D` is to be
added now** — it would add a name and no new physical contract.

**And do not make `FSE2D` depend on a future `SE2D`.** FSE is not a complete SE scan repeated: the
crusher window, echo spacing and midpoint placement belong to the train. If a shared contract
exists it is *lower* than both, something like `SpinEchoCore` feeding single-SE and `TSEShot`
separately — the same sibling argument that kept `GRE3DTR` out of `GRE2DTR`. **`SpinEchoCore` is
not proposed**; the evidence supports only "it may exist", and what would settle it is comparing
Cartesian SE, SE-EPI and DWI SE-EPI for a shared refocusing contract that differs only in readout.

**FID is `PRIMITIVE_COMPOSITION_COVERED`**, not a missing module. A sequence can be fundamental and
common without deserving its own Module, and FID is a good future `NO_NEW_MODULE` calibration case.

**DREAM is reclassified from "likely `NO_NEW_MODULE`" to `YELLOW / DEFER`.** It is a B1/B0/TxRx
mapping method, and the question — does it retain an irreducible contract once the shared STEAM
physics is extracted — cannot be answered before the STEAM boundary is. DREAM is also **not** to be
grouped with BOLD: BOLD is a contrast/application, commonly GRE-EPI, expected to compose existing
physics plus time-series scheduling, which is caller policy.

**Skill fixtures are curated, not cumulative.** `skill_v0/migrated/` holds three calibration
fixtures because they cover three qualitatively different reasoning classes. A future candidate
joins them only if it exposes a reasoning class or failure mode the set does not already cover;
conversion alone is not a reason. Ordinary evidence stays in `plans/<candidate>/`.

---

# 9c. STEAM calibration result

`plans/steam/` holds the record. Closed as `UNDECIDED / YELLOW /
REFERENCE_NOT_INDEPENDENT`.

**What it established.** The Skill met a genuinely open candidate and stopped instead of forcing
an extraction — the first prospective stop, and the thing the calibration existed to test. It also
found that no permissively licensed executable stimulated-echo reference exists in PyPulseq or
Pulseq MATLAB, so the candidate is **evidence-blocked, not rejected**.

**What it did NOT establish.** It rejected the *broad, waveform-defined* grouping — "things that
look like the same three-pulse pattern" — on the grounds that the observed uses differ in RF
count, flip pattern, gradient realisation, storage and mixing timing, and readout coupling.

```text
broad waveform-defined STEAM abstraction      -> unsupported
physics-defined StimulatedEcho / STEAM        -> STILL UNRESOLVED
```

It is **not** a finding that stimulated echo should not be a Module. A future abstraction could
plausibly own the requested coherence pathway, storage and mixing intervals, a semantic
stimulated-echo time, required signed gradient-moment relationships, RF-centre timing
relationships and the suppression of unwanted pathways — if those parameterise while keeping clear
correctness conditions. If instead the general case requires the caller to specify nearly the whole
RF/gradient/timing/coherence program, the abstraction is a second sequence DSL and the physics
belongs in lower-level semantics plus family-specific kernels. **That boundary is deliberately not
resolved in this phase and does not block `SaturationPrep`.**

**Do not** keep searching for references merely to convert STEAM to GREEN. The calibration has produced
what it was for.

**DREAM.** The run observed a cross-component coupling — dephasing moment, readout waveform, and
the positions of the stimulated echo and FID inside one ADC window. That is **candidate evidence**
strengthening `DREAM -> DEFER / revisit later`. It does **not** justify `DREAM -> NEW_KERNEL`. A
later DREAM fine scan asks what remains once the shared stimulated-echo semantics are understood.

**Generic lessons promoted into the Skill** (and only these): a `status` value for "stopped before
deciding"; that a licence can block a candidate before its physics does, with the four uses of a
source distinguished; and that waveform or event-shape similarity is a discovery heuristic while
promotion requires a shared physical solve. The candidate record stays an ordinary mining record,
**not** a permanent calibration fixture.

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
    └── scan-source registry, incl. MRsources/MRzero-Core (discovery + oracle,
        NOT a design witness -- its playground builds with PyPulseq)

CALIBRATION 0  -- COMPLETE
StimulatedEcho / STEAM
    └── UNDECIDED / YELLOW / REFERENCE_NOT_INDEPENDENT
    └── demonstrated prospective STOP capability: a genuinely open candidate,
        and no extraction was forced
    └── rejected the BROAD waveform-defined candidate only; the physics-defined
        StimulatedEcho / STEAM abstraction remains INTENTIONALLY OPEN and does
        not block anything that follows

CALIBRATION 1  -- COMPLETE
SaturationPrep
    └── NEW_LEAF / GREEN
    └── simple extraction without overdesign; the scope refusals held

CALIBRATION 2  -- COMPLETE
Spiral
    └── NEW_LEAF / GREEN, after archaeology and implementation debugging
    └── PR #23 used as an evidence corpus and then closed as superseded
    └── the m1 compiler defect isolated WITHOUT Spiral and fixed in PR #30;
        the m0 claim did not reproduce and was not ported; the three helpers
        deferred for want of a second consumer

RETROSPECTIVE  -- COMPLETE
compare supervised skill runs
    └── tools/module_mining/skill_v0/prospective_calibration_retrospective.md
    └── the three outcomes differ on purpose: stop, extract small, take apart a
        large historical artifact and keep only what survives review

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
