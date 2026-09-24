# SeqCraft Fine Scan Playbook

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../../README.md) for which one to edit.

**Status:** Local execution guide  
**Date:** 2026-09-18  
**First recommended pilot:** TSE / FSE  
**Relationship to other documents:**
- [`archive/2026-09-18_original_mining_plan.md`](../archive/2026-09-18_original_mining_plan.md) — overall strategy and milestones.
- [`archive/2026-09-18_coarse_scan.md`](../archive/2026-09-18_coarse_scan.md) — coarse-scan evidence and candidate map.
- this file — how to execute a candidate fine scan locally.

---

# 1. Goal

A fine scan turns a coarse `ModuleCandidate` hypothesis into one of:

```text
NEW_LEAF
NEW_KERNEL
NEW_IMAGING
EXTEND_EXISTING
PROMOTE_NOTEBOOK
COMPOSITION_ONLY
NO_MODULE
ARCHITECTURE_REVIEW
MORE_EVIDENCE_REQUIRED
```

It does this using **source-code understanding plus executable comparison**.

The fine scan is not a source-to-source rewrite.

The desired chain is:

```text
multiple working references
        ↓
shared MRI contract
        ↓
candidate boundary
        ↓
SeqCraft extraction
        ↓
compile
        ↓
deterministic comparison
        ↓
parameter sweep
        ↓
promotion decision
```

---

# 2. Why the fine scan should be local

A fine scan requires more than GitHub browsing:

- full source and helper traversal;
- executable Python environments;
- MATLAB/Octave or pre-generated `.seq`;
- temporary sequence outputs;
- k-space calculations;
- waveform inspection;
- pytest;
- parameter sweeps;
- error injection;
- repeated edits to an extracted SeqCraft implementation.

Therefore the primary working environment should be a local SeqCraft checkout.

GitHub/web research remains useful for provenance and extra references.

---

# 3. Non-negotiable architectural rules

1. **Do not modify the compiler to make a candidate pass.**
2. **Do not add semantics to `LogicBlock`.**
3. **Do not create one class per source sequence.**
4. **Do not hide acquisition policy inside a lower-level module.**
5. **Do not trust source-code similarity as validation.**
6. **Do not trust one open-source implementation as ground truth.**
7. **Do not let the AI self-certify equivalence.**
8. **Do not copy incompatible licensed code into SeqCraft.**
9. **A fine scan may end with `NO_MODULE` or `MORE_EVIDENCE_REQUIRED`.**
10. **Promotion requires deterministic tests.**

---

# 4. Recommended local layout

Start outside the runtime package:

```text
SeqCraft/
  src/seqcraft/
  tests/
  examples/

  tools/
    module_mining/
      candidates/
      references/
      adapters/
      fingerprints/
      reports/
      sweeps/
```

Candidate-specific example:

```text
tools/module_mining/candidates/tse/
  candidate.yaml
  findings.md

  references/
    seqcraft_fse/
    pypulseq_tse/
    pulseq_matlab_tse/
    openmrf_tse/

  adapters/
    run_seqcraft.py
    run_pypulseq.py
    load_matlab_seq.py

  reports/
    invariant_table.md
    comparison_report.md
    sweep_report.md
```

Keep this tooling experimental. Do not treat these paths as public API.

---

# 5. Fine-scan workflow

## Pass 0 — freeze provenance

Before analysis, record:

```yaml
repo:
commit:
path:
license:
language:
dependency_versions:
entrypoint:
reference_parameters:
```

If the repository is GPL/AGPL or otherwise incompatible with SeqCraft's intended licensing, use it as **behavioral/physics evidence**, not as code to transplant.

---

## Pass 1 — code archaeology

Read the source and answer:

### 1. What does the user configure?

Examples:

- FOV;
- matrix;
- thickness;
- flip;
- echo spacing;
- turbo factor;
- b-value;
- trajectory angle;
- TI/TE/TR;
- bandwidth/dwell.

### 2. Which calculations encode MRI physics?

Examples:

- RF-center timing;
- gradient moment balance;
- k-space centre placement;
- crusher balance;
- b-value;
- echo spacing;
- radial/spiral rotation;
- balanced moments.

### 3. Which code only exists because Pulseq is flat?

Examples:

- manual split gradients;
- manual block stitching;
- repeated add/merge operations;
- event fragments whose only purpose is satisfying one-gradient-per-axis block constraints.

Those are usually **not** candidate API concepts; SeqCraft's compiler should absorb them.

### 4. Which choices are acquisition policy?

Examples:

- ky ordering;
- acceleration pattern;
- interleave ordering;
- dummy scan count;
- slice ordering;
- cardiac segment schedule.

Keep these caller-controlled unless the physics requires otherwise.

### 5. What is the natural reusable unit?

Possible answers:

```text
leaf event bundle
readout
preparation
repeating shot/TR kernel
complete imaging acquisition
no reusable unit
```

Write the answer before implementation.

---

# 6. Build a semantic decomposition

Produce one short diagram per reference.

Example TSE:

```text
Excitation
    ↓
initial dephase
    ↓
┌──────────────────────────────────────┐
│ Refocusing + crusher                 │
│ phase encode                         │
│ spin-echo Cartesian readout          │  × N echoes
│ phase rewind                         │
└──────────────────────────────────────┘
    ↓
tail / spoiler / TR fill
```

Then mark each component:

```text
EXISTING_SEQCRAFT
MISSING_PRIMITIVE
CANDIDATE_OWNED
CALLER_POLICY
REFERENCE_ARTIFACT
```

---

# 7. Build the invariant table before the candidate API

An invariant is a physical/timing property the reference is trying to maintain.

Do not derive the candidate API first.

For every invariant record:

```yaml
name:
physical_meaning:
source_evidence:
measurement:
tolerance:
family_specific: true|false
```

Example TSE invariants:

- RF effective centres follow the intended CPMG spacing.
- Spin echo lies at the intended location relative to refocusing centres.
- `kx = 0` at the semantic echo sample.
- `ky` equals the requested line at that sample.
- crusher moments are balanced around refocusing RF centre.
- effective TE is determined by the echo carrying central k-space.
- echo spacing remains constant across the train.
- shot duration and TR constraints are respected.

The invariant table becomes the specification for validation.

---

# 8. Reference adapters

The adapter layer should be thin and experimental.

Conceptually:

```python
class ReferenceSequence:
    sequence: pypulseq.Sequence
    parameters: dict
    definitions: dict
    provenance: dict
    semantic: dict
```

Do not freeze this type yet.

## Python references

Preferred route:

```text
reference function/script
    ↓
pypulseq.Sequence
    ↓
ReferenceSequence
```

When possible, refactor the *test harness* to call the source's builder without editing the source implementation itself.

## SeqCraft reference

Use the notebook-local implementation or extracted helper with the same parameters:

```text
LogicBlock
    ↓
sc.compile
    ↓
pypulseq.Sequence
```

Keep access to the original LogicBlock as well because SeqCraft-specific exact event digest tests are useful.

## MATLAB Pulseq

Initial route:

```text
MATLAB source
    ↓ execute or use reproducibly generated output
.seq
    ↓
pypulseq.Sequence.read
    ↓
ReferenceSequence
```

The MATLAB source is still inspected for semantics.

The `.seq` provides the executable cross-language artifact.

Do not build a MATLAB parser before this simple route proves insufficient.

---

# 9. SequenceFingerprint should be a comparator stack

Avoid one opaque hash.

Use progressively stronger layers.

## L0 — metadata sanity

- duration;
- RF count;
- ADC count;
- sample count;
- gradient event count;
- sequence definitions;
- label inventory.

## L1 — exact SeqCraft event digest

When the reference is a SeqCraft/raw-pypulseq extraction where exact identity is expected:

```text
absolute event time
+
event content hash
```

This is the strongest refactor regression test.

## L2 — canonical physical waveform

Compare:

```text
RF(t)
Gx(t)
Gy(t)
Gz(t)
ADC sample times
```

Prefer exact event knots / piecewise-linear evaluation over unnecessarily dense arbitrary resampling.

Block count must not be an equivalence criterion.

## L3 — acquisition semantics

Compare:

- k-space position at every ADC sample;
- echo sample location;
- echo polarity;
- phase/partition labels where meaningful;
- RF/receiver phase relationship;
- semantic TE/TR/TI values.

## L4 — family-specific physics

Examples:

### TSE
- CPMG/refocusing spacing;
- crusher balance;
- echo-centre symmetry;
- effective TE.

### Diffusion
- b-value;
- b-matrix/tensor;
- direction;
- cross terms;
- TE.

### bSSFP
- zeroth moments balanced per TR;
- RF phase cycling;
- TE≈TR/2 where intended.

### Spiral
- trajectory extent;
- centre timing;
- interleaf rotation;
- gradient/slew;
- ADC segmentation.

### Velocity encoding
- M1 target;
- VENC;
- moment compensation.

No candidate becomes GREEN without the relevant L4 checks.

---

# 10. Difference reports

Return measurements, not only pass/fail.

Example:

```yaml
result: fail

metadata:
  duration_error_s: 0
  adc_count_delta: 0

waveform:
  gx_max_error_hz_m: 0

acquisition:
  k_adc_max_error_per_m: 4.52
  adc_time_max_error_s: 0

physics:
  effective_te_error_s: 0
  crusher_balance_error_per_m: 11.3

first_failure:
  echo: 4
  check: crusher_balance
```

The AI can reason about this report, but deterministic code produces it.

---

# 11. Candidate extraction

Only after the shared invariant table exists should the SeqCraft API be proposed.

Ask in this order:

1. Can an existing module already express it?
2. Can an existing module be extended without making its API incoherent?
3. Is this one reusable physics concept?
4. Does it shorten multiple reference implementations?
5. Can it be tested independently?
6. Does it preserve caller choice?
7. Is the natural output one `LogicBlock`?
8. Is the natural level leaf, preparation, kernel, or imaging?

Candidate implementation should follow existing SeqCraft conventions:

- `opts` explicit;
- `__init__` designs;
- `build` assembles;
- calls pure;
- do not mutate stored events;
- use semantic timing methods where the tree cannot know the meaning;
- return a LogicBlock;
- keep compiler unaware of producer type.

---

# 12. Parameter sweep

A single reference protocol is not sufficient.

Choose a small structured sweep that varies the candidate's meaningful axes.

Example TSE:

```text
echoes:          1, 8, 16
matrix:          64, 128
echo spacing:    minimum, minimum + margin
partial Fourier: 1.0, 0.75
ordering:        linear, interleaved, centric
refocus flip:    180, reduced flip if supported
```

Do not explode into a Cartesian product initially.

Use pairwise/high-information cases.

---

# 13. Metamorphic tests

Some truths do not require another implementation.

Examples:

### TSE

- turbo factor 1 should reduce to a single-echo spin-echo family.
- changing ky ordering must not change the readout waveform.
- the echo assigned to central ky must determine effective TE.
- increasing turbo factor should preserve earlier echo timing when other timing parameters are unchanged.

### Radial

- rotating by angle `φ` rotates k-space by `φ`.
- FOV changes `dk` as expected.
- reversing angle/sign should produce the corresponding transformed trajectory.

### Diffusion

- zero b-value should remove diffusion weighting.
- rotating direction should rotate the gradient vector while preserving target b-value.
- negating a direction preserves scalar b-value.

Metamorphic tests are particularly valuable when independent open-source references disagree.

---

# 14. GREEN / YELLOW / RED

## GREEN

- executable reference exists;
- candidate boundary is clear;
- comparator passes required levels;
- parameter sweep passes;
- relevant metamorphic tests pass;
- SeqCraft package tests pass;
- no compiler changes;
- provenance/license recorded.

Result:

```text
candidate may enter PR review
```

## YELLOW

Examples:

- one semantic invariant remains ambiguous;
- MATLAB reference unavailable but Python evidence is strong;
- two references disagree in a way that may be representational or historical;
- API boundary not fully settled.

Result:

```text
human MRI review
```

## RED

Examples:

- reference cannot be reproduced;
- physical mismatch persists;
- candidate needs source-specific compiler behavior;
- abstraction only wraps code without removing meaningful duplication;
- licensing/provenance unresolved.

Result:

```text
do not promote
```

---

# 15. First pilot: TSE / FSE

## Reference set

Minimum:

```text
SeqCraft examples/fse_2d
official PyPulseq write_tse.py
```

Then add:

```text
official Pulseq MATLAB writeTSE.m
OpenMRF TSE
```

## Main architecture hypothesis

The coarse scan suggests:

```text
kernel:
    TSE/FSE shot

imaging:
    FSE2D
```

with:

```text
HASTE = FSE/TSE configuration
```

This is only a hypothesis.

The fine scan should prove or reject it.

## Specific questions

1. Is one full echo train the reusable unit?
2. Does the initial excitation belong in that unit?
3. Which dephasing/rewinding arithmetic belongs to `CartesianLine` versus the TSE kernel?
4. Should phase encoding be passed as a list of line indices per shot?
5. Can the same shot kernel support:
   - SE (`echoes=1`);
   - FSE;
   - HASTE;
   - future variable-flip TSE?
6. Are crusher semantics fully owned by `Refocusing`, or does the train need additional moment coordination?
7. How should `time_to_echo(n)` and effective TE be exposed?
8. Which labels are kernel-owned versus imaging-owned?

## Initial acceptance measurements

At minimum:

```text
RF centres
refocusing intervals
ADC times
k_adc
kx at echo
ky at echo
echo spacing
effective TE
gradient moments around refocusing
shot duration
hardware legality
```

## Initial deliverables

```text
candidate.yaml
findings.md
reference adapters
comparison report
parameter sweep report
candidate boundary decision
```

Do not require final package promotion in the first pilot.

The workflow itself is the first deliverable.

---

# 16. After TSE/FSE

Recommended next candidates:

```text
RadialReadout
GRE3DTR
SaturationPrep
DiffusionPrep/Encoding
SpiralReadout + trajectory utility
BalancedSSFPTR
```

The order may change based on what the fine-scan tooling exposes.

---

# 17. Definition of a successful fine-scan system

The system is ready for batch use when a local session can repeatedly do:

```text
choose candidate
    ↓
collect references
    ↓
produce invariant table
    ↓
run references
    ↓
generate/refine SeqCraft extraction
    ↓
run comparator
    ↓
run sweep
    ↓
emit GREEN/YELLOW/RED report
```

without inventing a new validation framework for every candidate.

Family-specific physics checks are expected.

The orchestration and evidence format should become reusable.

---

# 18. Immediate local-session task

Start with TSE/FSE.

Do **not** start by coding `TSEShot`.

Start by:

```text
1. inspect current SeqCraft FSE notebook;
2. inspect official PyPulseq write_tse.py;
3. write side-by-side semantic decomposition;
4. list shared invariants and disagreements;
5. create the smallest runnable Python reference adapters;
6. compare the existing SeqCraft FSE implementation numerically;
7. implement only the comparator checks needed to explain real differences;
8. decide the candidate boundary;
9. then extract/promote code.
```

This order is important.

The fine scan is first an evidence-building process, and only second a code-generation process.

---

# 19. Four rules the third candidate added

These came out of `GRE3DTR`, and the first one came out of a real mistake. They apply to every
candidate from here on, and they are the part of this playbook a batch process would most need.

## 19.1 Derive every physical mode from its own contract

The first `GRE3DTR` reasoned, in effect:

```text
slab-selective  = shaped RF + Gz
non-selective   = the same thing with Gz removed
```

The result was spatially non-selective and it still played the slab path's **3 ms shaped sinc**.
It satisfied every k-space and timing invariant, compiled legally, and put three milliseconds into
every echo time for a pulse that selected nothing. Every official Pulseq 3D reference uses a short
hard block pulse instead.

> **Alternative physical modes must be re-derived from their own physical contract. Do not
> implement one mode by subtracting events from another unless independent reference evidence
> establishes that equivalence.**

A mode *name* is not a contract. Before implementing, state what the scanner actually plays:

```text
RF family / shape
RF duration, and where its default comes from
selection gradient: present or absent
intrinsic rephasing / moment requirement
encoding gradients
ADC / readout semantics
spoiler / rewinder structure
timing consequence
the reference evidence for each of the above
what is deliberately left overridable
```

### The mode contract table

Any candidate with modes, variants or mode-dependent defaults gets this table **before** code:

| | mode A | mode B |
|---|---|---|
| RF family / shape | | |
| RF duration, default source | | |
| selection gradient | | |
| intrinsic rephasing requirement | | |
| encoding role | | |
| ADC / readout semantics | | |
| timing consequence | | |
| reference evidence | | |
| intentionally user-overridable | | |

and a differential contract beside it:

```text
mode A -> mode B
    changes:            ...
    must stay invariant: ...
```

The table is not an implementation plan. It is the physical contract the implementation has to be
shown to realise, and writing it is what stops "mode B is mode A minus an event" becoming a design
rule by accident.

## 19.2 Inspect the emitted sequence, not only what was measured from it

Correct k-space, a legal compile and a passing notebook do not prove the emitted experiment is the
intended one — the RF mistake above had all three.

For at least one representative protocol of **every promoted mode**, read the compiled `.seq`
and ask:

> If I inspected this file without reading the Python, would its RF, gradients, ADC and timing
> tell the same physical story as the API name and the docstring?

`tools/module_mining/inspect_emitted.py` prints what a file says about itself — pulse family,
duration, which gradients accompany the RF, whether it is selective at all. It is a **review aid
and prints rather than asserts**: it belongs beside the analytic tests, the comparator and the
simulation, not inside any of them.

## 19.3 A reference comparison proves what it measured, and no more

"Validated against reference X" is not a statement until it says against *what*. Record the scope
with the result:

```text
validated against the official GRE3D reference for:
    k-space lattice, semantic centre, partition spacing, requested encoding

not used to establish:
    the slab-selective RF implementation
    the waveform decomposition
    equal signal when TE differs
```

The GRE3D lattice comparison was correct and it established the encoding geometry. It said nothing
about the RF family, and for a while the RF family was wrong. A comparator result must not quietly
grow into a stronger physical claim than the measurement supports.

## 19.4 A shared-leaf change needs a dependency-impact map

Two leaf changes in this phase, with opposite consequences:

```text
Excitation   + rephaser_area_per_m, + build(rephase=False)
             additive, default unchanged -> no consumer changes unless it opts in

CartesianLine  ADC sample-count legality tightened
               a contract change -> every consumer inherits it
```

So before merging a shared-leaf change, map it and classify every consumer:

```text
changed leaf -> direct consumers -> transitive consumers -> examples and notebooks
```

```text
NO_BEHAVIOR_CHANGE
NEW_EARLY_REFUSAL_FOR_PREVIOUSLY_INVALID_INPUT
INTENTIONAL_BEHAVIOR_CHANGE
UNKNOWN / NEEDS_REVIEW
```

**Read the imports rather than the names.** `EPI2D` sounds like a readout that would compose
`CartesianLine` and does not — it owns its own train and ADC geometry, so a `CartesianLine`
contract change does not reach it. A map built from naming would have got that wrong in both
directions.

A green suite is necessary and is not this. The suite says nothing broke; the map says who could
have.

# 20. Realisation is not contract

**Added 2026-09-23, from the T2Prep implementation.**

A fine scan reads a small number of implementations, and one of them is usually the one whose
licence lets us look closely. Everything that implementation does is then in front of us in
detail, and the temptation is to write all of it down as the contract.

> **Implementation witnesses tell us how someone realised the physics. They do not define the
> reusable physical contract.**

Do not go looking for a reason to keep the detail. **Classify it**, and let the classification say
where it goes:

```text
Does domain evidence define this detail as part of the family itself?
    yes  ->  contract invariant

Do independent implementations converge on it because it is intrinsic
to the same physical claim?
    yes  ->  evidence for a contract invariant

Does the detail instead solve an additional, separately named physical problem?
    yes  ->  realisation / mode / option, carrying its own claim

otherwise
         ->  implementation detail; do not promote it
```

The third branch is the one that is easy to get backwards, because it is the branch where the
detail has an obvious good reason behind it. **That reason is not a promotion.** A detail that
exists to solve a *separate* named problem is evidence that it belongs in a realisation, a mode or
an option with a claim of its own — not evidence that it belongs in the parent family's contract.
"It is there for a good reason" is a reason to name it and to state what it claims; it is not a
reason to make it an invariant of a family that does not make that claim.

The second branch has a narrower catch: convergence counts only when it is convergence on the
*same* physical claim, and only between witnesses that are actually independent. Two ports of one
original are one witness, and the registry's `independence_note` is what decides that.

### What it looked like in practice

The T2Prep scan's only permissively licensed executable witness realises the tip-up as a composite
`+270x` then `−360x`. The scan recorded that as part of the `mlev4-hard` mode, and the first
implementation shipped it as the only behaviour.

Walk the branches and it never reaches the first two.

**Does domain evidence define it as part of the family?** No. The domain reference describes the
family as `+90`, one or more refocusing pulses, `−90`, spoiler, and names composite pulses
separately as an *optional* robustness choice. MLEV-4 defines the *refocusing* train's phase
cycling and says nothing about the tip-up at all.

**Do independent implementations converge on it?** No. The second, copyleft witness uses adiabatic
tip pulses instead, so on this point the witnesses do not agree.

**Does it solve an additional, separately named problem?** Yes — reduced sensitivity to B0 and B1
error is exactly what a composite pulse is reached for, and the domain reference names it in those
terms. So the detail lands in the third branch: a realisation with a claim of its own. And that is
the part the first implementation got backwards. Having a good reason behind it is what makes it a
*named option*, not what makes it an invariant — especially here, where the claim that reason
refers to is one this module does not make.

The cost of getting it the wrong way round is not academic: the composite is longer in RF, so it
pays more relaxation and more RF energy, and it carried a timing convention — "the centre of a
composite" — that the domain evidence does not define and that no measurement we had could settle.

And when the two realisations were finally measured against each other under a B1 error, **no
consistent advantage for the composite appeared** — at one point in the sweep it was the worse of
the two. The robustness intuition that made the composite look contract-worthy was not there to
find. That is not a finding about which part of the preparation dominates B1 sensitivity, and the
record does not claim one; it is only enough to say that the option must not be presented as a
robustness feature.

### The smell

Ask of any invariant in a candidate record: **which branch put this here?** If the answer is "the
witness did it", it is an implementation detail that has not been labelled as one yet. If the
answer is "it is there for a good reason", that is the third branch, and the record owes a name
for the option and a statement of what it claims.

The corollary is about evidence classes rather than rules. The canonical form of a family is a
*domain* question, and the registry now carries an `educational-synthesis` role for sources that
answer it well — a curated secondary account is often the fastest way to see what practitioners
consider canonical, and it never outranks a handbook or a primary paper.
