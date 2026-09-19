# SeqCraft AI Module Mining & Batch Conversion Plan

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../README.md) for which one to edit.

> **Status:** Active implementation plan — coarse scan completed; the TSE/FSE fine-scan pilot is
> **COMPLETE — GREEN** (`TSEShot`, `FSE2D` and the `CartesianLine` echo-geometry extension are
> extracted and verified); the next candidate is `RadialReadout`, planned in
> [`next_phase_plan.md`](next_phase_plan.md)  
> **Primary goal:** Build an AI-assisted pipeline that mines reusable MRI sequence concepts from mainstream open-source Pulseq/PyPulseq sequences, converts high-confidence concepts into SeqCraft modules, and verifies that the conversion preserves the reference sequence's physical behavior.  
> **Architectural constraint:** Do **not** change the SeqCraft compiler to accommodate individual source sequences. The compiler architecture is treated as stable. This project operates above the `LogicBlock` boundary.

---

## 1. Executive summary

SeqCraft already has a clean architectural split:

```text
Reusable MRI modules / arbitrary Python
                ↓
            LogicBlock
                ↓
        SeqCraft compiler
 placement → boundaries → legalization → emission → verification
                ↓
        pypulseq.Sequence
```

The current bottleneck is no longer how to represent or compile overlapping sequence events. The bottleneck is **module coverage**: the current module library covers only a subset of the MRI building blocks needed to construct vendor-like families of sequences.

Writing modules manually one by one is expensive because each module requires:

- understanding the source sequence's MRI physics;
- deciding what is genuinely reusable versus sequence-specific;
- identifying the right parameters and invariants;
- preserving timing and waveform semantics;
- checking corner cases;
- validating the extracted implementation against a known-working reference sequence.

The proposed solution is an AI-assisted **module mining and reference-preserving extraction pipeline**, not a simple source-to-source converter.

The key principle is:

> **AI proposes abstractions; deterministic tooling proves preservation.**

The pipeline should therefore have three logically separate jobs:

```text
1. Discover
   Open-source sequence corpus
       ↓
   Identify recurring MRI concepts and missing SeqCraft capabilities

2. Extract
   Working reference sequence
       ↓
   SeqCraft Module / extension of existing Module

3. Verify
   Reference implementation
       ↕
   Canonical physical fingerprint
       ↕
   SeqCraft implementation after compilation
```

The immediate first step should **not** be implementing the full converter. The immediate first step should be a **coarse corpus scan** that produces a provisional `ModuleCandidate` inventory and a SeqCraft module coverage map. This can begin now, before the final `ModuleCandidate` schema or `SequenceFingerprint` comparator is implemented.

---

# 2. Core architectural principles

The project should inherit and preserve the architectural rules already established in SeqCraft.

## 2.1 Compiler is out of scope

The compiler remains responsible for:

```text
LogicBlock tree
    ↓ placement
absolute-time events
    ↓ boundaries
Pulseq block boundaries
    ↓ legalization
split / sum / limit checks
    ↓ emission
pypulseq.Sequence
    ↓ verification
legal output
```

The module-mining system must not introduce special compiler cases such as:

- "spiral mode" in the compiler;
- "diffusion mode" in the compiler;
- source-repository-specific legalization behavior;
- module-specific block boundary rules;
- compiler knowledge of `Module` subclasses.

If a converted module only works after adding special behavior to the compiler, that is a strong signal that the abstraction is wrong or the compiler has a separate reproducible defect that should be addressed independently.

## 2.2 `LogicBlock` remains the only runtime interface

A generated reusable component may be:

- an `sc.Module` subclass;
- a plain function;
- a domain-specific class with multiple named outputs.

But what enters the compiler remains a `LogicBlock`.

## 2.3 Modules own MRI physics; the tree owns timing relationships

A module may know:

- where its effective RF center is;
- where k-space center occurs;
- its echo times;
- its gradient moments;
- its intended b-value relationship;
- its trajectory geometry;
- whether a particular echo is reversed;
- the physical meaning of a parameter.

A `LogicBlock` should not learn those concepts.

## 2.4 Extraction, not invention

The most important rule for AI-generated modules is:

> **A module should be extracted from a working sequence, not designed in the abstract.**

The working reference implementation acts as the behavioral specification.

A stronger rule for promotion into the shared module library is:

> **A reusable abstraction should ideally be supported by multiple independent consumers or source sequences.**

One source may justify a candidate. Multiple sources justify generalization.

## 2.5a Package modules are Python-only, for now

*Added 2026-09-18, with the TSE/FSE boundary decision.*

Package-level SeqCraft `Module` implementations are **Python-only for now**. MATLAB Pulseq
references remain useful as physics and architecture evidence, and as optional `.seq` validation
evidence, but they are **not an implementation target**.

Do not design a MATLAB Module API, a cross-language Module abstraction or a MATLAB implementation
layer. A one-off MATLAB run to settle a `.seq` question is fine; it must not expand the module
architecture.

The practical consequence for mining: where a candidate has both a Python and a MATLAB reference,
the **Python** one is the priority to obtain and execute.

## 2.5 No class-per-sequence explosion

The mining system must be able to conclude:

```text
NO_NEW_MODULE
```

or:

```text
EXTEND_EXISTING_MODULE
```

rather than always creating a new class.

For example, a multi-echo GRE sequence may indicate that the correct abstraction is an extension of an existing Cartesian readout rather than a new `MultiEchoGRE` wrapper.

This is essential. Without this rule, batch mining will recreate the complexity SeqCraft intentionally removed.

---

# 3. Proposed end-to-end pipeline

The target system is:

```text
                     ┌───────────────────────────┐
                     │ Open-source Pulseq corpus │
                     └─────────────┬─────────────┘
                                   ↓
                     ┌───────────────────────────┐
                     │  Stage A: Corpus scanner  │
                     │ metadata + source mapping │
                     └─────────────┬─────────────┘
                                   ↓
                     ┌───────────────────────────┐
                     │ Stage B: Coarse AI mining │
                     │ concepts, families, gaps  │
                     └─────────────┬─────────────┘
                                   ↓
                     ┌───────────────────────────┐
                     │ Provisional candidates    │
                     │ + coverage map            │
                     └─────────────┬─────────────┘
                                   ↓
                     ┌───────────────────────────┐
                     │ Stage C: Candidate merge  │
                     │ deduplicate / generalize  │
                     └─────────────┬─────────────┘
                                   ↓
                     ┌───────────────────────────┐
                     │ Stage D: Fine scan        │
                     │ static + dynamic analysis │
                     └─────────────┬─────────────┘
                                   ↓
                     ┌───────────────────────────┐
                     │ Stage E: Extraction       │
                     │ Module / extension code   │
                     └─────────────┬─────────────┘
                                   ↓
             ┌─────────────────────┴─────────────────────┐
             ↓                                           ↓
   Reference sequence path                    SeqCraft converted path
             ↓                                           ↓
        canonical                             LogicBlock → compile
        fingerprint                                      ↓
             └─────────────────────┬─────────────────────┘
                                   ↓
                     ┌───────────────────────────┐
                     │ Stage F: Differential     │
                     │ verification              │
                     └─────────────┬─────────────┘
                                   ↓
                     ┌───────────────────────────┐
                     │ Stage G: Parameter sweep  │
                     │ + metamorphic tests       │
                     └─────────────┬─────────────┘
                                   ↓
                         GREEN / YELLOW / RED
                                   ↓
                        reviewed batch PRs
```

---

# 4. Important answer: does `ModuleCandidate` need code first?

## Short answer

**No.**

The first version of `ModuleCandidate` should deliberately **not** be implemented as a rigid Python class or runtime schema.

Instead, use a provisional Markdown/YAML record during the coarse scan.

Recommended sequence:

```text
coarse scan now
    ↓
provisional candidate records
    ↓
observe real variation across sequence families
    ↓
revise candidate fields
    ↓
second/fine scan
    ↓
freeze stable schema
    ↓
only then implement Pydantic/dataclass validation if useful
```

This is preferable because the first 20–50 sequences will teach us what information is actually common across:

- Cartesian readouts;
- EPI;
- spiral;
- radial / UTE;
- diffusion;
- TSE / HASTE;
- balanced SSFP;
- inversion / saturation / T2 preparation;
- spectroscopy;
- CEST / MT;
- flow encoding.

Prematurely encoding the schema in Python risks designing it around GRE/EPI terminology and then forcing fundamentally different modules into the wrong model.

## Recommended schema lifecycle

### Level 0 — free-form research notes

For the first few sequences, allow the agent to write short structured notes.

Purpose:

- identify terminology;
- identify repeated concepts;
- discover fields we did not anticipate.

### Level 1 — provisional YAML template

After roughly 5–10 sequences, begin using a consistent but non-enforced YAML format.

Missing fields are allowed. New fields are allowed.

### Level 2 — stabilized YAML schema

After a broad coarse scan, normalize candidate fields and aliases.

At this stage, each candidate should be machine-readable enough for ranking and deduplication.

### Level 3 — code-backed schema

Only after the field set is stable should we consider:

- `dataclasses.dataclass`;
- `TypedDict`;
- Pydantic;
- JSON Schema.

This code would belong to the **module-mining tooling**, not `seqcraft.design` or the SeqCraft runtime API.

## Suggested schema freeze criterion

Do not freeze `ModuleCandidate` until all of the following are approximately true:

- at least **20–40 reference sequences** have been scanned;
- at least **5–8 distinct sequence families** are represented;
- leaf, kernel, preparation, encoding, and readout concepts have all appeared;
- fewer than roughly **10% of new candidates require adding a new top-level schema field**;
- candidate merge/deduplication can be performed without repeatedly reading the raw source again.

The exact numbers are not sacred; the principle is to freeze only after the schema is learned from the corpus.

---

# 5. Phase 1 — Coarse corpus scan: start immediately

This phase can begin now and should require little or no new SeqCraft code.

## 5.1 Goal

Produce:

1. a sequence corpus inventory;
2. a preliminary module coverage map;
3. a set of `ModuleCandidate` records;
4. evidence showing which source sequences support each candidate;
5. an initial estimate of which missing abstractions unlock the largest number of sequence families.

This phase is about **discovery and prioritization**, not code generation.

## 5.2 Recommended initial corpus

### Tier A — authoritative / baseline references

Start with sources closest to Pulseq itself:

- `pulseq/pypulseq` example scripts;
- `pulseq/pulseq` MATLAB demo sequences;
- official or author-maintained Pulseq extensions/examples.

These are useful because:

- APIs are conventional;
- sequence intent is usually clear;
- dependencies are limited;
- they provide common reference implementations;
- they cover many mainstream sequence families.

### Tier B — mature research repositories

Then add representative projects for missing families, for example:

- spiral implementations;
- diffusion sequences;
- CEST / MT;
- flow / phase contrast;
- advanced radial / UTE / ZTE;
- spectroscopy;
- pTx if relevant to SeqCraft's intended scope.

### Tier C — exploratory / noisy repositories

Only later include:

- one-off lab scripts;
- repositories with vendor/private dependencies;
- incomplete projects;
- old Pulseq versions;
- code without reproducible parameter sets.

These may still reveal useful abstractions, but should not define the initial schema.

## 5.3 Corpus manifest

Each scanned sequence should receive a lightweight record such as:

```yaml
source_id: pulseq-pypulseq/write_tse.py
repository: pulseq/pypulseq
commit: <sha>
language: python
license: MIT
sequence_family: TSE
reference_quality: authoritative
entrypoint: examples/scripts/write_tse.py

features:
  rf: [excitation, refocusing]
  readout: cartesian
  phase_encoding: true
  multi_echo: true
  labels: false
  arbitrary_gradients: false
  preparation: []

execution:
  runnable: unknown
  external_data_required: false
  scanner_specific_dependency: false

mining_status: coarse_scanned
```

This manifest can initially be written by the AI skill itself.

No Python schema is required yet.

---

# 6. Coarse scan methodology

The coarse scan should intentionally be inexpensive.

It should answer:

> What MRI concepts occur in this sequence, and which of them are absent or insufficiently generalized in SeqCraft today?

It should **not** yet prove waveform equivalence.

## 6.1 Inputs examined during coarse scan

For each source sequence:

- README / documentation around the sequence;
- main sequence construction script;
- helper functions called by the main script;
- obvious parameter definitions;
- loops that construct repetitions, echoes, phase encodes, shots, or segments;
- existing SeqCraft module inventory.

Runtime execution is optional at this phase.

## 6.2 Questions the AI should answer

For every source sequence:

### A. What are the conceptual components?

Examples:

```text
slice-selective excitation
Cartesian readout
phase encode
readout rewinder
RF spoiling
multi-echo train
refocusing pulse
crusher
inversion preparation
spiral interleaf
diffusion lobe pair
velocity encoding pair
fat saturation
```

### B. Which components already exist in SeqCraft?

Classify each as:

```text
EXISTING_SUFFICIENT
EXISTING_NEEDS_EXTENSION
NEW_CANDIDATE
SEQUENCE_SPECIFIC
UNCLEAR
```

### C. What is repeated?

Look for repeated code patterns across:

- echoes;
- phase-encode lines;
- TRs;
- shots;
- slices;
- directions;
- interleaves;
- preparation cycles.

Repetition is often evidence of a reusable component.

### D. What arithmetic encodes MRI meaning?

Examples:

- gradient area derived from `1/FOV`;
- prephaser area derived from half a readout area;
- timing derived from RF center to k-space center;
- diffusion timing derived from b-value;
- TSE symmetry around a refocusing pulse;
- spiral trajectory rotation;
- RF phase cycling;
- balanced gradient moment constraints.

This arithmetic is often the real reason a module should exist.

### E. What should remain a caller choice?

Examples:

- phase-encode ordering;
- acceleration pattern;
- number and arrangement of segments;
- slice ordering;
- sequence-level dummy scans;
- contrast-specific scheduling.

The AI should avoid hiding these choices inside a module unless the physical abstraction requires it.

---

# 7. Provisional `ModuleCandidate` format

The first coarse-scan schema should be explicitly labeled **provisional**.

Example:

```yaml
candidate_id: diffusion_encoding_v0
name: DiffusionEncoding
status: provisional
category: encoding

summary: >
  Bipolar or monopolar diffusion-weighting gradients arranged around
  one or more refocusing pulses, parameterized by b-value and direction.

source_evidence:
  - repository: pulseq/pulseq
    path: matlab/demoSeq/writeEpiDiffusionRS.m
    role: primary
  - repository: some/research-repo
    path: dwi_sequence.py
    role: supporting

seqcraft_relationship:
  decision: NEW_CANDIDATE
  related_existing:
    - Refocusing
    - EPI2D

physics:
  invariant_concepts:
    - gradient moment relationship across refocusing
    - target diffusion weighting
    - vector direction projection onto physical gradient axes

parameters_observed:
  - b_value
  - direction
  - max_grad
  - max_slew
  - timing

semantic_times:
  - refocusing_center

possible_interface:
  shape: multi_output_or_block_with_gap
  notes: >
    May be more honest as one LogicBlock containing a central gap, or as
    a domain object exposing named pre/post pieces. Do not decide during
    coarse scan.

reuse_evidence:
  source_sequence_count: 2
  sequence_families:
    - SE-EPI-DWI
    - spiral-DWI

risk:
  level: high
  reasons:
    - b-value verification required
    - direction scaling must preserve vector magnitude
    - concomitant timing interactions

open_questions:
  - Is b-value owned by the module or by a higher-level diffusion kernel?
  - Should crusher contributions be included or composed externally?

confidence: medium
```

## Why this should not be code-backed yet

Several fields above are intentionally ambiguous:

- `possible_interface` may evolve;
- not every module has semantic times;
- some components naturally have multiple outputs;
- some candidates should become an extension to an existing module;
- some concepts are functions rather than classes;
- some sequence-specific compositions should not be promoted at all.

The coarse scan should reveal the stable common denominator.

---

# 8. Candidate aggregation and deduplication

After scanning individual sequences, the next job is not conversion. It is **cross-sequence synthesis**.

## 8.1 Merge equivalent candidates

Examples:

```text
RadialReadout
RadialSpoke
RadialLine
GoldenAngleSpoke
```

may represent one physical primitive plus caller-controlled angle ordering.

Likewise:

```text
DiffusionLobes
StejskalTannerEncoding
DWIEncoding
```

may or may not be the same abstraction depending on what the source code reveals.

The AI should merge only after comparing:

- invariant physics;
- inputs;
- outputs;
- semantic timing requirements;
- intended composition boundaries.

## 8.2 Candidate outcome types

Every aggregated candidate must end in one of five decisions:

```text
A. NEW_MODULE
B. EXTEND_EXISTING_MODULE
C. COMPOSE_EXISTING_MODULES
D. KEEP_SEQUENCE_LOCAL
E. NEEDS_MORE_EVIDENCE
```

This decision is more important than the proposed class name.

## 8.3 Promotion heuristic

A candidate is more attractive when:

- it appears in multiple independent sequences;
- it carries nontrivial MRI arithmetic;
- multiple users would otherwise reimplement that arithmetic;
- it composes with several sequence families;
- its caller-facing parameters are physically meaningful;
- it can compile and validate in isolation;
- it removes sequence boilerplate without stealing sequence-programming choices from the caller.

---

# 9. Produce a SeqCraft Module Coverage Map

The most valuable deliverable from the first coarse scan may be the coverage map itself.

Example structure:

```text
RF
├── Excitation                     ✓ existing
├── Refocusing                     ✓ existing
├── Inversion                      ✓ via IRPrep
├── Saturation                     △ missing/generalize
├── FatSat                         ? candidate
└── adiabatic preparation          ? candidate

Preparation
├── IR                             ✓
├── T2Prep                         ?
├── MT / CEST                      ?
└── diffusion preparation          ?

Encoding
├── Cartesian phase encoding       ✓
├── 3D partition encoding          ?
├── diffusion encoding             ?
├── velocity encoding              ?
└── trajectory rotation            ?

Readout
├── Cartesian line                 ✓
├── multi-echo Cartesian           △ extension
├── EPI                            ✓
├── radial spoke                   ?
├── spiral interleaf               ?
├── UTE radial                     ?
└── ZTE / PETRA                    ?

Kernel
├── GRE TR                         ✓
├── spin-echo kernel               ?
├── TSE echo train                 ?
├── bSSFP TR                       ?
└── diffusion SE kernel            ?
```

This map should be regenerated from candidate data rather than maintained manually once the tooling matures.

---

# 10. Ranking: maximize sequence coverage, not module count

The goal should not be "how many modules can AI generate?"

The goal should be:

> **Which missing primitive unlocks the largest amount of useful sequence space?**

A simple priority model can eventually use:

```text
priority ≈
    reuse_count
  × sequence_family_count
  × downstream_unlock_count
  × confidence
  ÷ implementation_risk
```

This is only a planning metric, not a scientific score.

Examples of high-leverage primitives may include:

### Spiral interleaf

```text
Excitation + SpiralInterleaf
        ↓
GRE spiral
SE spiral
spiral fMRI
spiral MRF
spiral diffusion
```

### Diffusion encoding

```text
DiffusionEncoding + Refocusing + EPI2D
        ↓
SE-EPI DWI

DiffusionEncoding + Refocusing + SpiralInterleaf
        ↓
spiral DWI
```

### TSE echo train

```text
Excitation + Refocusing + Cartesian echo train
        ↓
SE
TSE / FSE
HASTE
variants with different echo ordering
```

### Radial spoke

```text
Excitation + RadialSpoke
        ↓
radial GRE
radial SE
UTE family
potentially stack-of-stars variants
```

This coverage-oriented ranking should determine conversion order.

---

# 11. Phase 2 — Fine scan

Once a candidate is selected for serious conversion, perform a fine scan.

The fine scan should combine **static semantic analysis** with **dynamic tracing**.

## 11.1 Why source code alone is insufficient

Source code tells us why calculations are being made, but AI can misinterpret local variables or miss hidden helper behavior.

For example:

```text
gx_pre
gx_read
gx_spoil
```

are syntactically similar gradients but semantically different.

The source is needed to discover:

- FOV relationships;
- TE/TR definitions;
- loop semantics;
- sequence ordering;
- intended echo center;
- physical parameter meanings.

## 11.2 Why final `.seq` alone is insufficient

The final sequence tells us what plays, but not why.

It cannot reliably tell us:

- which external variable was FOV;
- which loop variable was the phase-encode index;
- which delay was specifically chosen to satisfy TE;
- which waveform relationship represents a reusable invariant;
- what the intended user-facing API should be.

## 11.3 Combined model

```text
source code
    ↓
static semantic interpretation

runtime / generated sequence
    ↓
dynamic event trace

        ↓ combine

semantic sequence model
```

This combined model is the correct input to module extraction.

---


# 11A. Local fine-scan execution protocol

The fine scan should be performed **directly against source code in a local working environment**.

This is the point where the project changes from corpus research into executable sequence engineering.

The local environment is preferable because a fine scan needs to do all of the following together:

- read the full source and helper functions;
- trace parameters into RF / gradient / ADC construction;
- execute the reference implementation;
- generate and inspect `.seq` files;
- calculate k-space and semantic timing;
- compare multiple implementations under identical parameters;
- run parameter sweeps;
- intentionally inject errors to test the comparator;
- iterate on a SeqCraft extraction without committing premature runtime API decisions.

A browser/GitHub scan remains useful for provenance and comparison, but it is no longer sufficient as the primary working surface.

## 11A.1 Unit of work

A fine scan is organized around **one candidate**, supported by **multiple reference implementations**.

Example:

```text
candidate: TSE/FSE shot

references:
  - SeqCraft examples/fse_2d
  - official PyPulseq write_tse.py
  - official Pulseq MATLAB writeTSE.m
  - OpenMRF TSE
```

The goal is not to translate each reference independently.

The goal is to discover the **shared physical contract** that survives across those implementations.

## 11A.2 Recommended local directory

Do not place exploratory mining code directly inside `src/seqcraft`.

Use a tooling area first:

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

A candidate may then have:

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
    seqcraft_reference.py
    pypulseq_reference.py
    matlab_seq_reference.py

  reports/
    reference_inventory.md
    invariant_table.md
    comparison_report.md
```

The exact filesystem layout is not part of SeqCraft's public architecture. It can change freely while the workflow matures.

## 11A.3 Three passes of a fine scan

### Pass 1 — code archaeology

Read the candidate source implementations and build a semantic map.

For each source, identify:

```text
inputs
  ↓
waveform design calculations
  ↓
timing calculations
  ↓
event placement
  ↓
loops / ordering / acquisition policy
  ↓
outputs / definitions / labels
```

Explicitly separate:

- **invariant physics** — should belong to the reusable abstraction;
- **parameter variation** — should become module parameters;
- **sequence policy** — should remain with the caller;
- **implementation artifacts** — split gradients, temporary variables, block layouts that exist only because Pulseq is flat.

Do not write the SeqCraft candidate yet unless the boundary is already obvious.

### Pass 2 — executable reference normalization

Create a small adapter that converts each runnable reference into a common record.

Provisional conceptual interface:

```python
ReferenceSequence(
    sequence=...,
    parameters=...,
    definitions=...,
    source=...,
    semantic_times=...,
)
```

The exact Python type does not need to be frozen yet.

The adapter should make it possible for later code to ask the same questions regardless of whether the source was:

- PyPulseq;
- SeqCraft;
- MATLAB Pulseq exported to `.seq`.

For Python references, prefer direct execution and direct access to `pypulseq.Sequence`.

For MATLAB references, the first practical target is:

```text
MATLAB / Octave source
      ↓
generated .seq
      ↓
PyPulseq Sequence.read(...)
      ↓
common reference analysis
```

This gives a useful cross-language boundary without requiring a MATLAB AST or tracing framework.

### Pass 3 — extraction and differential validation

Only after the shared physical contract is understood:

```text
reference implementations
        ↓
shared invariants
        ↓
candidate API
        ↓
SeqCraft implementation
        ↓
compile
        ↓
differential comparison
        ↓
parameter sweep
```

The first extraction should be deliberately conservative.

If the candidate can be expressed by extending an existing module, prefer that over a sibling class.

If the references disagree materially, do not average them into a new abstraction. Record the disagreement and keep the candidate YELLOW until it is understood.

## 11A.4 What the fine scan must record

Every fine scan should produce at least five artifacts.

### A. Reference inventory

For each reference:

- repository;
- commit SHA;
- license;
- path;
- language;
- dependency/version assumptions;
- runnable/not runnable;
- parameter set used;
- output `.seq` if generated;
- known limitations.

### B. Semantic decomposition

A short diagram, for example:

```text
Excitation
    ↓
dephase / winder
    ↓
[Refocusing + crusher
 + phase encode
 + readout
 + rewind] × N
    ↓
tail / spoiler
```

### C. Invariant table

For TSE/FSE, an example table might include:

| Invariant | Meaning | How measured |
|---|---|---|
| refocusing spacing | CPMG timing | RF effective centres |
| echo midpoint | spin echo centered between refocus pulses | ADC semantic sample vs RF centres |
| `kx=0` at echo | readout refocusing | k-space at ADC |
| requested `ky` at echo | phase encoding | k-space at ADC |
| crusher balance | unwanted pathway suppression | integrated moments |
| effective TE | contrast-defining echo | echo carrying central ky |
| train duration | shot timing | event timeline |

### D. Candidate boundary decision

One of:

```text
NEW_LEAF
NEW_KERNEL
NEW_IMAGING
EXTEND_EXISTING
PROMOTE_NOTEBOOK
COMPOSITION_ONLY
NO_MODULE
ARCHITECTURE_REVIEW
```

with explicit reasons.

### E. Validation report

The report should say which levels were passed:

```text
L0 metadata
L1 exact event digest
L2 canonical waveform
L3 acquisition semantics
L4 family-specific physical invariants
parameter sweep
metamorphic tests
```

## 11A.5 Reference adapters should stay thin

A reference adapter is not a new abstraction layer for users.

Its job is only to make external sequence implementations measurable.

It should not:

- redesign the source;
- repair suspicious source behavior silently;
- normalize away real physical differences;
- infer missing parameters from desired output;
- make an incompatible reference look equivalent.

If the source has a bug or ambiguous behavior, preserve it as reference evidence and record the issue.

## 11A.6 Start with Python references

The first fine-scan pilot should use Python sources first:

```text
SeqCraft notebook
+
official PyPulseq
+
Python research repository when useful
```

This minimizes environment complexity while the workflow itself is still experimental.

Once the adapter/fingerprint path works, add MATLAB references through generated `.seq` files.

The MATLAB source is still required for semantic interpretation; the `.seq` is the executable physical artifact.

## 11A.7 Do not collapse the comparator into one hash

A single `SequenceFingerprint` hash is too blunt.

Use a structured comparator with layers:

```text
representation identity
        ↓
waveform equivalence
        ↓
acquisition equivalence
        ↓
family-specific physics
```

This matters because two correct implementations may:

- split blocks differently;
- encode one gradient as a trapezoid and another as an arbitrary waveform;
- place harmless labels differently;
- use different temporary event decomposition while producing the same physical timeline.

Likewise, two sequences can have nearly identical waveforms but differ in a crucial semantic quantity such as b-value or echo assignment.

`SequenceFingerprint` should therefore mean a **collection of deterministic measurements and comparisons**, not necessarily one cryptographic digest.

## 11A.8 Fine-scan stopping rules

A candidate should stop and remain YELLOW when:

- independent references disagree on the physical contract;
- the proposed API hides a sequence-policy choice;
- validation requires changing the compiler;
- the only reusable code is a thin wrapper;
- the candidate cannot be tested independently;
- a source license prevents safe code reuse and the physics cannot yet be re-derived independently;
- a family-specific invariant is not understood well enough to test.

A fine scan may legitimately conclude:

```text
NO_NEW_MODULE
```

or:

```text
MORE_EVIDENCE_REQUIRED
```

That is a successful result.

## 11A.9 TSE/FSE as the first pilot

TSE/FSE is the recommended first local fine scan because it has four useful forms of evidence:

```text
SeqCraft working notebook
+
official PyPulseq implementation
+
official Pulseq MATLAB implementation
+
independent modular implementation
```

The pilot should first answer:

1. What is the smallest reusable unit — refocused echo, full echo train, or one complete shot?
2. Which arithmetic belongs in the kernel rather than in `Refocusing`, `CartesianLine`, or `PhaseEncode`?
3. Which ordering choices must remain caller data?
4. Can HASTE be expressed entirely as TSE/FSE configuration?
5. Does the current notebook-local `FSE2D` naturally decompose into:

```text
TSE/FSE shot kernel
        ↓
FSE2D imaging composition
```

6. Which exact measurements are required for `SequenceFingerprint` to prove preservation?

The first pilot should not optimize for a perfect class name. It should optimize for a clear physical contract and a trustworthy validation path.

---

# 12. Dynamic tracing strategy

## 12.1 PyPulseq references

For Python sources, possible instrumentation includes:

- wrapping `Sequence.add_block`;
- observing `make_*` event creation calls;
- capturing final `Sequence` blocks;
- capturing definitions;
- recording parameter values;
- recording loop indices or execution context where practical.

The simplest useful first implementation may only require executing the reference script and parsing the resulting `pypulseq.Sequence`.

## 12.2 MATLAB Pulseq references

For MATLAB references, do not require a full MATLAB semantic parser initially.

Use:

```text
MATLAB source code
      +
reference .seq output
```

The source provides intent and parameter relationships.

The `.seq` file provides the executable physical reference.

If a reproducible MATLAB/Octave environment is available later, execution tracing can be added.

## 12.3 Non-executable references

If a repository cannot be executed because of:

- proprietary files;
- unavailable hardware dependencies;
- vendor APIs;
- missing data;
- obsolete packages;

then the candidate may still contribute to coarse mining, but must not receive a GREEN conversion status based solely on AI reasoning.

---

# 13. `SequenceFingerprint` comparator

Unlike `ModuleCandidate`, the comparator **does need code** once we move into actual conversion and verification.

This is the deterministic heart of the project.

## 13.1 Purpose

Given two implementations that may use different source structure or Pulseq block boundaries, determine whether they represent the same intended sequence behavior to an appropriate tolerance.

The comparator should avoid relying solely on:

- source code similarity;
- class structure;
- Pulseq block count;
- serialized `.seq` byte identity;
- event object identity;
- event content hashes alone.

## 13.2 Why the existing event digest is still valuable

SeqCraft already has a strong exact-equivalence pattern for internal extraction:

```text
absolute event time + event content hash
```

This is ideal when:

- reference and extracted implementation use the same PyPulseq event representation;
- the goal is to prove a refactor did not move or change any event.

Keep this comparator as **Level 1 exact structural equivalence**.

But external source conversions need a more canonical physical comparator.

## 13.3 Proposed fingerprint levels

### Level 0 — metadata sanity

Compare:

- total duration;
- block count for diagnostics only;
- definitions;
- event counts;
- number of ADC samples;
- RF count;
- label inventory.

Failure here often indicates a major mismatch.

### Level 1 — exact SeqCraft event digest

Where applicable:

```text
[(absolute_time, event_content_hash), ...]
```

This remains the strongest internal regression check.

### Level 2 — canonical waveform timeline

Represent each physical channel in a canonical form:

```text
RF magnitude / phase / frequency versus time
Gx versus time
Gy versus time
Gz versus time
ADC sample times
triggers / labels / extensions
```

Use exact knots whenever possible rather than arbitrary high-rate resampling.

### Level 3 — acquisition semantics

Compare quantities reconstruction and sequence timing depend on:

- ADC timestamps;
- k-space coordinates at each ADC sample;
- echo centers;
- echo polarity;
- echo number labels;
- phase-encode / partition labels;
- relevant sequence definitions.

### Level 4 — physical invariants

Depending on module family, compare:

- gradient moments;
- zeroth/first moments;
- RF effective center;
- integrated flip behavior where meaningful;
- b-value / b-matrix;
- trajectory extent;
- trajectory rotation;
- balanced moments;
- crusher moments;
- TE/TR timing relationships.

Different module families may register additional domain checks.

## 13.4 Comparator output

Do not return only `True/False`.

Return a structured difference report:

```yaml
result: fail

summary:
  total_duration_error_s: 0.0
  adc_time_max_error_s: 0.0
  kspace_max_error_per_m: 12.4
  gradient_waveform_max_error_hz_m: 0.0

failures:
  - check: kspace_at_adc
    axis: x
    echo: 3
    sample: 0
    reference: -145.2
    candidate: -132.8

interpretation:
  - likely reversed/prephaser mismatch before echo 3
```

AI may interpret this report, but deterministic code must generate the measurements.

---

# 14. Comparator canonicalization rules

This area deserves careful design because physical equivalence is not identical to source representation equivalence.

## 14.1 Gradients

Two gradients may be physically identical even when one is represented as:

- one trapezoid;
- several split trapezoids;
- an arbitrary gradient waveform;
- a waveform produced after overlap summation.

Therefore compare the continuous piecewise-linear waveform or exact knot representation whenever possible.

## 14.2 RF

Compare at minimum:

- RF waveform timing;
- amplitude/phase samples;
- frequency offset;
- phase offset;
- event delay;
- effective center used by timing logic.

For some references, exact RF representation may be impossible across toolchains; such cases should use explicit tolerances and domain-specific checks.

## 14.3 ADC

ADC timing is especially important.

Compare:

- number of samples;
- sample dwell;
- first sample time;
- every sample time when practical;
- frequency and phase offsets where relevant.

## 14.4 Labels and extensions

Labels are stateful and therefore cannot be treated as decorative metadata.

Compare the effective label state at relevant acquisition points, not only the literal number of label events.

## 14.5 Block boundaries

Block boundaries should normally be ignored as an equivalence criterion.

Different legal block decompositions may represent the same physical sequence.

However, block type can still be checked when a specific module contract intentionally requires something stronger, for example preserving an unsplit trapezoid because downstream tooling expects it.

---

# 15. Phase 3 — Module extraction

Only after a candidate is selected and a reference execution path exists should AI generate code.

## 15.1 Generated artifact set

Each conversion should ideally generate:

```text
src/seqcraft/modules/.../<candidate>.py

tests/modules/test_<candidate>.py

tests/reference/<source_id>/...

module_mining/candidates/<candidate>.yaml

module_mining/reports/<candidate>_verification.md
```

The exact repository layout can be adjusted, but module-mining metadata should remain outside the runtime package where possible.

## 15.2 Generated module rules

AI must follow the existing SeqCraft module contract:

### Design in `__init__`

Waveforms that do not change per call should be created once.

### Assemble in `build`

`build` should primarily:

- derive/copy stored events;
- scale or rotate as required;
- place components in a `LogicBlock`;
- return the block.

### Calls must be pure

Calling the module repeatedly with the same inputs must not mutate stored state or previous results.

### Pass `opts` explicitly

Never rely on process-global PyPulseq defaults.

### Module does not compile itself

The module emits a `LogicBlock`; compilation belongs outside.

### Semantic timing belongs to the module

Examples:

- `time_to_echo()`;
- `time_to_center()`;
- echo-specific timings;
- effective inversion center;
- trajectory sample indices.

Only expose what callers actually need.

### Do not recreate Pulseq block constraints inside the module

If two independently meaningful gradients overlap, express the intended timing and let the compiler legalize them.

---

# 16. Decide module versus extension versus local composition

Before generating a new class, AI must answer this checklist.

## 16.1 New module is justified when

- there is reusable MRI physics;
- the arithmetic is easy to get subtly wrong;
- multiple sequence families need it;
- it has a stable user-facing parameter vocabulary;
- it can be tested in isolation;
- extracting it significantly simplifies consumers.

## 16.2 Extend an existing module when

The new sequence differs by a capability that belongs naturally to an existing primitive.

Examples might include:

- one versus multiple Cartesian echoes;
- monopolar versus bipolar readout polarity;
- optional partial Fourier behavior that belongs to the readout itself;
- rotation/interleaf parameters that are intrinsic to a trajectory primitive.

## 16.3 Compose existing modules when

The source sequence is merely a new arrangement of already-reusable concepts.

This should remain sequence code rather than become a library class unless repeated consumers emerge.

## 16.4 Keep sequence-local when

- only one consumer exists;
- parameters are protocol choices rather than reusable physics;
- abstraction does not shorten or clarify the source;
- the proposed module would merely wrap a handful of calls without owning meaningful arithmetic.

---

# 17. Phase 4 — Differential verification

The generated module should never be accepted because "the code looks right."

It must be executed against the reference.

## 17.1 Reference-preserving test pattern

```text
reference implementation
        ↓
reference SequenceFingerprint

SeqCraft module
        ↓
LogicBlock
        ↓
sc.compile
        ↓
candidate SequenceFingerprint

        ↓ compare
```

## 17.2 Required baseline checks

Every converted candidate should test:

- isolated compilation succeeds;
- sequence timing check succeeds;
- gradient amplitude/slew limits are respected;
- relevant RF / ADC dead-time rules are respected;
- reference duration matches;
- relevant gradient waveforms match;
- ADC sample times match;
- k-space at ADC samples matches;
- semantic timing values match the actual generated sequence;
- labels/extensions match where reconstruction depends on them.

## 17.3 Exact versus tolerant verification

Verification should record which checks are:

```text
EXACT
NUMERICAL_TOLERANCE
SEMANTIC_INVARIANT
NOT_APPLICABLE
```

Do not silently downgrade an exact mismatch to "close enough."

---

# 18. Phase 5 — Automated parameter sweeps

A module that reproduces one nominal protocol may still be wrong as a reusable abstraction.

Therefore generated tests should explore the parameter space.

## 18.1 Generic sweep dimensions

Where applicable:

- matrix sizes: 32 / 64 / 128;
- small / medium / large FOV;
- minimum allowed timing;
- longer timing;
- low/high flip angle;
- several gradient axes;
- several scanner gradient limits;
- several slew limits;
- odd/even matrix sizes;
- center and edge phase-encode positions.

## 18.2 Readout-specific sweeps

Examples:

- echoes = 1 / 2 / 8;
- monopolar / bipolar;
- full / partial Fourier;
- reversed polarity;
- number of spiral interleaves;
- trajectory rotation angle;
- radial spoke angle.

## 18.3 Diffusion-specific sweeps

Examples:

- b = 0;
- low b;
- high b near hardware limit;
- x/y/z unit directions;
- oblique normalized direction;
- sign reversal;
- multiple refocusing arrangements if supported.

---

# 19. Metamorphic tests

Metamorphic tests are particularly useful because many MRI errors still produce legal sequences and plausible images.

Examples:

## Phase encoding

```text
center line → phase-encode moment = 0
line offset changes sign → moment changes sign
|line offset| equal → |moment| equal
```

## Cartesian readout

```text
reverse polarity
    → k-space sample order reverses
    → magnitude extent remains unchanged
```

## FOV

```text
increase FOV
    → dk decreases as 1/FOV
```

## Multi-echo readout

```text
increase echo count
    → existing earlier echo locations remain unchanged
```

## Diffusion

```text
rotate normalized diffusion direction
    → intended b-value magnitude remains constant
    → physical gradient-axis projections change
```

## Balanced sequence component

```text
one TR
    → zeroth gradient moment returns to the required balance condition
```

These tests should be generated from the candidate's physical invariants where possible.

---

# 20. Simulation as a secondary verification layer

Waveform/fingerprint equivalence should be the primary automated contract because it is precise and fast.

Simulation is still valuable for selected cases.

Use simulation to detect errors that survive structural checks or to validate a novel abstraction at integration level.

Potential checks include:

- expected k-space coverage;
- image orientation;
- ghosting sensitivity;
- T2* decay behavior across multi-echo readouts;
- expected contrast trends;
- diffusion attenuation consistency;
- off-resonance behavior for EPI/spiral.

Simulation should not replace the physical comparator.

A wrong sequence can sometimes produce a plausible-looking image.

---

# 21. Confidence gates

Every candidate should receive an automated outcome.

## GREEN

Requirements should approximately include:

- reference implementation executes;
- source provenance is known;
- license permits the intended use;
- reference fingerprint can be produced;
- candidate compiles independently;
- baseline fingerprint checks pass;
- required semantic invariants pass;
- parameter sweep passes;
- SeqCraft test suite passes;
- no compiler modifications were required.

Outcome:

> Candidate may be proposed as a PR for human review.

## YELLOW

Examples:

- waveform comparison passes but semantic interpretation is uncertain;
- reference code contains ambiguous conventions;
- only one source supports a proposed abstraction;
- external toolchain prevents complete dynamic tracing;
- some physical checks require expert review.

Outcome:

> Generate code/report if useful, but require MRI expert review before promotion.

## RED

Examples:

- reference cannot be reproduced;
- candidate changes acquisition semantics;
- unexplained fingerprint mismatch;
- converter requires compiler special cases;
- proprietary dependency prevents verification;
- sequence intent cannot be reconstructed reliably.

Outcome:

> Keep the mining report; do not generate/promote package code.

---

# 22. Human review should focus on abstraction, not bookkeeping

The entire point of this project is to shift human effort away from repetitive waveform comparison.

Human reviewers should mainly answer:

- Is this the correct reusable MRI concept?
- Are the module boundaries natural?
- Are the parameters expressed in a useful physical vocabulary?
- Does the module hide a choice that should remain with the sequence author?
- Should this extend an existing module instead?
- Is the abstraction supported by enough independent usage?

The tooling should answer:

- Did any event move unexpectedly?
- Did ADC timing change?
- Did k-space change?
- Did moments change?
- Did the labels change?
- Did hardware legality change?
- Did the converted implementation remain deterministic and pure?

---

# 23. Recommended skill behavior

The eventual AI skill should not be one monolithic prompt.

Conceptually it should support several modes.

## `scan-corpus`

Input:

- repository list or search scope.

Output:

- corpus manifest;
- sequence family inventory;
- coarse candidate records;
- coverage map.

## `analyze-sequence`

Input:

- one source sequence.

Output:

- conceptual decomposition;
- existing SeqCraft coverage;
- candidate/extension/local-composition decisions;
- open questions;
- risk assessment.

## `merge-candidates`

Input:

- candidate collection.

Output:

- deduplicated abstraction set;
- evidence graph;
- proposed generalizations;
- priority order.

## `convert-candidate`

Input:

- selected candidate;
- source references;
- stable SeqCraft module conventions.

Output:

- proposed module code;
- reference adapter;
- generated tests;
- verification configuration.

## `verify-candidate`

Input:

- reference implementation;
- converted implementation.

Output:

- fingerprint comparison;
- parameter-sweep report;
- confidence gate.

Keeping these modes separate makes failures easier to reason about and prevents one hallucinated interpretation from flowing silently through the whole pipeline.

---

# 24. Suggested repository/tooling layout

Avoid placing mining infrastructure in SeqCraft's runtime layers.

One possible layout is:

```text
SeqCraft/
├── src/seqcraft/
│   └── modules/                 # promoted, reviewed runtime modules only
│
├── tools/
│   └── module_miner/
│       ├── README.md
│       ├── corpus.py
│       ├── fingerprint.py
│       ├── compare.py
│       ├── reference_runner.py
│       ├── candidate_schema.py   # only after schema stabilizes
│       └── adapters/
│           ├── pypulseq.py
│           └── pulseq_seqfile.py
│
├── module_mining/
│   ├── corpus.yaml
│   ├── coverage.md
│   ├── candidates/
│   ├── reports/
│   └── references/
│
└── tests/
    └── modules/
```

An alternative is to keep the entire mining system in a separate repository if the experiment becomes large. The important constraint is that it should not become part of the compiler dependency graph.

---

# 25. Provenance and licensing

Every mined candidate must preserve evidence of where the idea came from.

At minimum record:

- repository URL;
- repository commit SHA;
- source path;
- license;
- whether code was copied, adapted, or independently reimplemented from behavior;
- source paper/citation if available.

The desired default should be:

> Extract the physical abstraction and independently implement it using SeqCraft conventions, rather than mechanically copying large source fragments.

The original source remains the reference implementation and evidence.

Do not promote code when license compatibility is unclear.

---

# 26. Failure modes to explicitly guard against

## 26.1 Over-abstraction

Symptom:

- every sequence becomes a new class;
- modules encode acquisition ordering or protocol choices;
- many modules have one consumer.

Countermeasure:

- require `NEW_MODULE / EXTEND / COMPOSE / LOCAL / MORE_EVIDENCE` decision.

## 26.2 Under-abstraction

Symptom:

- the same difficult arithmetic appears in many generated sequence scripts;
- AI repeatedly rewrites the same trajectory or preparation calculations.

Countermeasure:

- cross-corpus candidate merging;
- repeated-code and repeated-physics evidence.

## 26.3 Testing only nominal parameters

Symptom:

- module passes the source protocol but breaks at different matrix/FOV/TE.

Countermeasure:

- parameter sweeps and metamorphic tests.

## 26.4 Comparing serialized Pulseq files byte-for-byte

Symptom:

- physically equivalent representations fail because block decomposition differs.

Countermeasure:

- canonical physical fingerprint.

## 26.5 Comparing only k-space extent

Symptom:

- mirrored trajectories or mislabeled echoes pass.

Countermeasure:

- signed k-space at every ADC sample;
- echo polarity and labels;
- semantic timing checks.

## 26.6 Letting AI self-certify

Symptom:

- model writes code and claims it is equivalent based on inspection.

Countermeasure:

- deterministic comparator is the acceptance authority.

## 26.7 Modifying the compiler to make a candidate pass

Symptom:

- source-specific compiler branches accumulate.

Countermeasure:

- compiler frozen for this project;
- compiler defects require an independent reproduction and separate review.

---

# 27. Recommended implementation sequence

The project should be implemented in the following order.

## Milestone 0 — Freeze scope and rules

Deliverables:

- this plan;
- explicit "compiler out of scope" rule;
- existing module inventory;
- source/provenance policy.

Success criterion:

- team agrees on the distinction between module mining, conversion, and verification.

## Milestone 1 — Coarse scan pilot

Scan approximately 10–15 representative sequences across at least:

- GRE;
- multi-echo GRE;
- EPI;
- SE-EPI;
- TSE/HASTE;
- MPRAGE;
- radial;
- UTE;
- spiral;
- diffusion.

Deliverables:

- corpus manifest;
- provisional candidate files;
- first coverage map;
- list of schema fields that changed during the scan.

**No converter implementation is required for this milestone.**

## Milestone 2 — Broad coarse scan

Expand to roughly 20–40 sequences and several additional families.

Deliverables:

- deduplicated candidate set;
- evidence graph;
- candidate outcome classification;
- module priority ranking;
- revised provisional schema.

Success criterion:

- new sequences rarely require new top-level candidate fields.

## Milestone 3 — Stabilize `ModuleCandidate`

Only now decide whether code-backed validation is useful.

Deliverables:

- documented schema v1;
- optional Pydantic/dataclass/JSON Schema implementation;
- migration script for earlier candidate files if necessary.

Important:

- this schema belongs to tooling, not SeqCraft runtime.

## Milestone 4 — Local fine-scan pilot + `SequenceFingerprint` v1

Run the first fine scan directly against source code in a local SeqCraft checkout.

Recommended pilot:

```text
TSE / FSE
```

Reference set:

- SeqCraft `examples/fse_2d`;
- official PyPulseq `write_tse.py`;
- official Pulseq MATLAB `writeTSE.m`;
- one independent modular/research implementation when useful.

Work in three passes:

```text
code archaeology
    ↓
reference adapters
    ↓
extraction + differential validation
```

Deliverables:

- candidate `findings.md`;
- refined `candidate.yaml`;
- runnable Python reference adapter(s);
- MATLAB `.seq` ingestion path;
- exact SeqCraft event-digest comparator;
- duration / event / ADC metadata comparator;
- gradient canonicalization;
- ADC timestamp comparator;
- k-space-at-ADC comparator;
- TSE-specific semantic checks;
- structured difference report;
- first parameter sweep.

Success criterion:

- the workflow can compare the notebook-local SeqCraft FSE implementation against at least one independent reference;
- harmless block-splitting differences do not automatically fail physical equivalence;
- seeded timing/moment/echo-order errors are detected;
- the candidate boundary can be stated clearly as `NEW_KERNEL`, `PROMOTE_NOTEBOOK`, `EXTEND_EXISTING`, or another explicit outcome.

This milestone is the point where `SequenceFingerprint` becomes code.

It is also where the provisional `ModuleCandidate` schema is tested against real executable evidence.

## Milestone 5 — Convert/promote one medium-complexity candidate

Good choices:

- radial spoke;
- multi-echo Cartesian readout extension;
- a simple saturation/preparation module.

Avoid choosing diffusion as the first comparator validation target because it adds b-value complexity before the basic machinery is proven.

Deliverables:

- one complete reference-preserving extraction;
- module code;
- tests;
- fingerprint report;
- parameter sweep;
- human review notes.

## Milestone 6 — Convert one high-complexity candidate

Suggested candidates:

- spiral interleaf;
- diffusion encoding;
- TSE echo train.

This validates whether the architecture handles genuinely difficult physics rather than only simple event packaging.

## Milestone 7 — Batch conversion

Only after the previous milestones succeed should conversion become batch-oriented.

Batch behavior:

```text
candidate queue
    ↓
reference runner
    ↓
AI extraction
    ↓
automatic tests
    ↓
GREEN → PR queue
YELLOW → expert review queue
RED → report only
```

---

# 28. What can be done immediately, before any new code

A useful amount of work can begin now.

## Immediate task A — inventory existing SeqCraft modules

Record for each existing module:

- category;
- physical responsibility;
- parameters;
- semantic timing queries;
- known consumers;
- extension points.

This becomes the deduplication reference for mining.

## Immediate task B — coarse scan official Pulseq/PyPulseq sequences

No runtime instrumentation is required initially.

For each sequence produce:

- one sequence summary;
- conceptual decomposition;
- existing module matches;
- candidate list;
- likely extension opportunities;
- source evidence.

## Immediate task C — build candidate index manually/YAML-first

Do not implement schema validation yet.

Use candidate records as living research documents.

## Immediate task D — produce coverage map

After the first 10–15 sequences, aggregate all findings into a coverage map.

## Immediate task E — identify first comparator target

Select a case where:

- reference runs easily;
- SeqCraft already has a close equivalent;
- physical equivalence can be checked independently.

The existing GRE notebook/package equivalence tests are a useful calibration point.

---

# 29. Suggested first coarse-scan batch

A practical first batch should deliberately span different abstraction types.

For example:

```text
1. PyPulseq GRE
2. PyPulseq EPI
3. PyPulseq SE-EPI
4. PyPulseq TSE
5. PyPulseq radial GRE
6. PyPulseq UTE
7. PyPulseq MPRAGE
8. Pulseq multi-echo GRE
9. Pulseq spiral
10. Pulseq diffusion EPI
11. Pulseq HASTE
12. Pulseq bSSFP / TrueFISP
13. Pulseq PRESS
14. Pulseq semi-LASER
15. Pulseq ZTE/PETRA
```

The goal is not to convert all 15. The goal is to discover the abstraction vocabulary.

After this batch, we should be able to answer much more confidently:

- what constitutes a leaf module;
- what constitutes a kernel;
- where readout variants should be parameters versus classes;
- which semantic timing fields recur;
- whether preparation modules need additional conventions;
- whether multi-output components are common enough to deserve explicit tooling support.

---

# 30. Expected output of the first coarse-scan milestone

A useful first report might look like:

```text
Scanned sequences: 15
Sequence families: 11

Existing modules sufficient: 6 concepts
Existing modules needing extension: 4 concepts
New candidate primitives: 9 concepts
Sequence-local compositions: 7
Unclear / needs evidence: 3

Highest-leverage missing candidates:
- SpiralInterleaf
- RadialSpoke
- DiffusionEncoding
- TSEEchoTrain
- SaturationPrep

Likely extensions:
- CartesianLine → multi-echo
- PhaseEncode → 3D partition/generalized encoding? [needs evidence]

Likely non-modules:
- MPRAGE scan ordering
- GRAPPA sampling pattern
- centric line ordering
```

The precise result should emerge from the scan rather than being assumed in advance.

---

# 31. Decision on `ModuleCandidate` after coarse scan

At the end of the broad coarse scan, explicitly review the schema.

Ask:

1. Which fields were present in almost every candidate?
2. Which fields were specific to only one family?
3. Which fields were repeatedly ambiguous?
4. Which fields are needed for ranking?
5. Which fields are needed for conversion?
6. Which fields are needed only for human explanation?
7. Which data can be computed later instead of stored?

Then split the schema into:

```text
Core candidate fields
    stable across all concepts

Optional family-specific metadata
    readout / RF / diffusion / preparation / etc.

Derived fields
    reuse count, coverage count, priority
```

Only after this step should a code-backed schema be considered stable enough to be worth maintaining.

---

# 32. Relationship between `ModuleCandidate` and `SequenceFingerprint`

These two concepts solve different problems and should not be conflated.

## `ModuleCandidate`

Question:

> **What reusable abstraction might exist here?**

Nature:

- semantic;
- exploratory;
- AI-heavy;
- initially flexible;
- evolves during corpus mining.

Therefore:

> **Start without code. Stabilize later.**

## `SequenceFingerprint`

Question:

> **Did the extracted implementation preserve the reference behavior?**

Nature:

- numerical;
- deterministic;
- test infrastructure;
- acceptance-critical.

Therefore:

> **Implement as code before trusting batch conversion.**

This distinction gives the project a natural two-speed workflow:

```text
Discovery can move fast and remain flexible.
Verification must move slowly and remain strict.
```

---

# 33. Definition of success

The project succeeds when SeqCraft module growth changes from:

```text
human finds one missing module
        ↓
human studies one sequence
        ↓
human implements module
        ↓
human manually validates it
        ↓
repeat
```

into:

```text
AI continuously scans sequence corpus
        ↓
coverage map shows missing reusable physics
        ↓
team selects high-leverage candidate
        ↓
AI extracts candidate from working references
        ↓
deterministic tooling proves preservation
        ↓
human reviews abstraction and MRI meaning
        ↓
promote module
```

The desired outcome is **not zero human review**.

The desired outcome is that human effort is spent on the hard scientific question — "is this the right abstraction?" — rather than repetitive code translation and waveform bookkeeping.

---

# 34. Recommended next action

The coarse scan has now been completed and has produced:

- a module coverage map;
- a provisional candidate taxonomy;
- candidate layer/promotion recommendations;
- two architecture findings:
  - composite preparations should remain under `preparation/`;
  - non-Cartesian trajectory design may deserve a non-Module utility layer.

> **Done, 2026-09-18.** Milestone 4 completed GREEN; see
> [`tse/boundary_decision.md`](tse/boundary_decision.md) for the outcome and
> [`next_phase_plan.md`](next_phase_plan.md) for what follows. The rest of this section is
> kept as the reasoning that selected TSE/FSE.

The next action is therefore **Milestone 4: local fine-scan pilot**.

Recommended first target:

```text
TSE / FSE
```

Recommended execution order:

```text
1. Read SeqCraft FSE notebook + official PyPulseq TSE source side by side.
2. Write a semantic decomposition and invariant table.
3. Create thin runnable reference adapters.
4. Compare the existing SeqCraft implementation to the reference numerically.
5. Add `SequenceFingerprint` checks only as required by real observed differences.
6. Bring in the MATLAB Pulseq `.seq` reference after the Python path works.
7. Decide whether the reusable abstraction is:
      TSE/FSE shot kernel
      + FSE2D imaging composition.
8. Run a small parameter sweep.
9. Record GREEN / YELLOW / RED and the final candidate action.
```

Do not start batch conversion yet.

Do not freeze the `ModuleCandidate` schema in Python yet.

Do not add trajectory infrastructure, diffusion taxonomy, or new module folders merely because the coarse scan suggested them. Those remain hypotheses until a fine scan produces executable evidence.

The purpose of the first fine scan is to prove the **workflow**, not merely to add the first new class.
