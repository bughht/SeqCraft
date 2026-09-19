# SeqCraft Fine Scan Playbook

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../README.md) for which one to edit.

**Status:** Local execution guide  
**Date:** 2026-09-18  
**First recommended pilot:** TSE / FSE  
**Relationship to other documents:**
- [`module_mining_plan.md`](module_mining_plan.md) — overall strategy and milestones.
- [`coarse_scan_2026-09-18.md`](coarse_scan_2026-09-18.md) — coarse-scan evidence and candidate map.
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
