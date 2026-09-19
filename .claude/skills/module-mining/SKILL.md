---
name: module-mining
description: Supervised, human-gated workflow for mining a reusable physical contract out of a corpus of MR sequence implementations and proving a SeqCraft Module realises it. Use when asked to evaluate a sequence family as a module candidate, run a fine scan, extract a Module from reference implementations, or decide GREEN/YELLOW/RED on one. Not for ordinary module edits.
---

# Module mining

This workflow exists because the obvious one is wrong. The obvious one is:

```text
find similar code -> wrap it -> make tests pass
```

That produces a wrapper with a confident name, and it produces it just as readily when there is no
reusable physics to find. The governing rule instead:

> **Extract a reusable physical contract from evidence first, then prove the implementation
> realises that contract.**

Evidence before API. Contract before code. Claim before comparison.

**You do not decide.** Every terminal outcome is a recommendation for a human, and
`NO_NEW_MODULE` is a success. Do not optimise for producing Modules.

## The pipeline

```text
source corpus
    -> candidate discovery
    -> reference evidence
    -> physical contract
    -> ownership / boundary analysis
    -> claim-scoped acceptance criteria
    -> reference adapter / measurements
    -> only then extraction
    -> validation
    -> GREEN / YELLOW / RED
    -> human review
```

Stages do not overlap. In particular, **no extraction before the acceptance claim is written
down** — an acceptance criterion invented after the code exists is a description of the code.

Depth for every stage: `tools/module_mining/plans/fine_scan_playbook.md`. This skill is the
operational layer over it and does not repeat it.

## Ownership, which is what most of the argument turns out to be

| layer | owns |
|---|---|
| leaf (`rf/ preparation/ encoding/ readout/`) | intrinsic physics and geometry, determinable from its own parameters |
| kernel | coupling **between** leaves that no leaf can see |
| imaging | acquisition policy — how many, which order |
| compiler | Pulseq legality, and nothing about any candidate |
| validation tooling | measures claims; never defines physics |

The test that settles most disputes:

> If a component cannot determine the correct value using only its own physical contract and
> parameters, it should not own that event.

Two corollaries that have each been needed: **reuse of semantics does not oblige reuse of an
emitted event**, and a composite that merely sequences leaves in a fixed order is composition, not
a kernel.

## The four rules a candidate must satisfy before promotion

Each came from a real failure. `reference/rules.md` has the full form, the failure behind it, and
the template; load it when a candidate reaches that stage.

**A. Physical mode contract.** Any candidate with meaningful modes or defaults publishes a mode
table (RF family/shape, RF duration and source, selection gradient, intrinsic rephasing, encoding
role, ADC semantics, timing consequence, reference evidence, user-overridable choices) plus an
explicit `mode A -> mode B: changes / must remain invariant` differential.

> A mode may not be implemented as "another mode minus one event" unless reference evidence
> establishes that equivalence.

**B. Emitted-sequence inspection.** At least one representative `.seq` per promoted physical mode,
read back with `python tools/module_mining/inspect_emitted.py <file.seq>`. The question is: *if
this file were inspected without reading the Python, would it tell the same physical story as the
API name and the documentation?*

**C. Claim-scoped reference validation.** Every comparison states both halves:

```text
this comparison establishes: ...
this comparison does NOT establish: ...
```

"Validated against X" may never mean more than what was measured.

**D. Shared-leaf dependency impact.** A candidate that changes an existing leaf maps
`changed leaf -> direct consumers -> transitive consumers -> examples/notebooks`, and classifies
each as `NO_BEHAVIOR_CHANGE`, `NEW_EARLY_REFUSAL_FOR_PREVIOUSLY_INVALID_INPUT`,
`INTENTIONAL_BEHAVIOR_CHANGE` or `UNKNOWN / NEEDS_REVIEW`. Grep for real composition, not for the
class name — a docstring mention is not a consumer.

## Acceptance claims are family-specific, on purpose

There is **no universal equality rule** and `SequenceFingerprint` is not one. The stable model is:

```text
generic measurement primitives  ->  family-specific check set  ->  structured DifferenceReport
```

What the three pilots accepted, none of which would have worked for the others:

| | claim that fit |
|---|---|
| TSE/FSE | event identity — it was promoting an existing implementation |
| Radial | trajectory / sample agreement — two legitimate constructions share no events |
| GRE3D | k-space lattice + signed z physics + timing — the selective path has no executable reference at all |

Your job is to help **choose and state** the claim, not to run every family through one comparator.

## Where candidates come from

`tools/module_mining/sources.yaml` is the scan-source registry: every corpus the workflow may
scan, with the two judgements that must be made **before** mining anything from it.

```text
discovery       tells us a family exists and is worth looking at
design-witness  an implementation whose construction can corroborate another's
oracle          an independent measurement path, for validating rather than designing
```

A corpus is never automatically a design-witness because it contains working sequences.
`MRsources/MRzero-Core` is registered `discovery` + `oracle` with `design-witness` **excluded**:
its playground builds with PyPulseq and simulates with MRzeroCore, so agreeing with a PyPulseq
example is one design agreeing with itself.

The registry's `role` is a property of the corpus and does **not** set `evidence[].role` for a
candidate. Re-ask, per candidate: independence relative to the *other* references cited here;
claim scope; shared implementation assumptions with PyPulseq/Pulseq; discovery or acceptance
evidence; and licence at the **file** level, not only the repository level.

Copyleft sources (AGPL/GPL, and MRzero additionally carries a EULA) are **physics and invariants
only** — no implementation text.

## Artifacts

Six structured records, in `reference/artifacts.md`, with a template at
`templates/candidate.yaml`. Validate with:

```bash
python tools/module_mining/schema.py path/to/candidate.yaml
```

`Candidate`, `Evidence`, `PhysicalContract`, `AcceptanceClaim`, `ValidationPlan`, `Decision`.

The schema is deliberately thin. Only 13 field paths survived all three pilots, and `validation`
survived as a key with **no** common sub-key — the evidence says family-specific validation must
stay free-form, so the schema requires that the section exist and does not prescribe its shape.
Do not add required fields that fewer than all three pilots needed.

## Coverage classification, at coarse-scan time

Before any decision, a coarse scan says what already covers a family. **"SeqCraft can already
build it" is not "SeqCraft already ships the right reusable abstraction"**, and the word
"duplicate" hides the difference. Vocabulary, with criteria in `reference/outcomes.md`:

```text
DIRECT_SHIPPED_DUPLICATE            a public shipped abstraction already exists
DEGENERATE_CASE_OF_EXISTING_MODULE  a parameter limit of a shipped abstraction
COMPOSITION_COVERED                 modules express it; no equivalent public abstraction
NOTEBOOK_ONLY_EXISTING              the abstraction lives only in an example notebook
PRIMITIVE_COMPOSITION_COVERED       plain composition IS the intended API
```

`DEGENERATE_CASE_OF_EXISTING_MODULE` may carry `ARCHITECTURE_REVISIT_CANDIDATE`: the parameter
limit works and the naming or layering may still be wrong.

## Outcomes, including the ones that stop

Full criteria in `reference/outcomes.md`. The catalogue:

```text
GREEN                                           may enter PR review
NO_NEW_MODULE                                   the corpus has no reusable contract
YELLOW - PHYSICAL_BOUNDARY_UNCLEAR
YELLOW - REFERENCES_DISAGREE
YELLOW - REFERENCE_NOT_INDEPENDENT
YELLOW - VALIDATION_ORACLE_UNTRUSTED
YELLOW - SHARED_LEAF_IMPACT_UNRESOLVED
YELLOW - COMPILER_CHANGE_CLAIM_NEEDS_MINIMAL_REPRODUCER
RED - WRAPPER_ONLY
RED - DUPLICATES_EXISTING_MODULE
RED - EXTRACTION_CHANGES_THE_PHYSICS
```

A corpus containing no new Module is a correct answer. Reaching `NO_NEW_MODULE` early is cheaper
than reaching `RED - WRAPPER_ONLY` after an extraction.

## Compiler changes are exceptional

Never change compiler behaviour because a candidate will not compile. The escalation, in order:

```text
candidate exposes apparent compiler problem
    -> remove the candidate implementation from the reproducer
    -> minimal raw PyPulseq / LogicBlock case
    -> prove a legal physical input is rejected or altered
    -> characterise the numerical / error floor independently
    -> only then consider a compiler change
```

Without a candidate-free reproducer the outcome is
`YELLOW - COMPILER_CHANGE_CLAIM_NEEDS_MINIMAL_REPRODUCER`, not a patch. A candidate that needs
source-specific compiler behaviour is `RED`.

## Validation layers, and deferring one honestly

```text
Layer 1  analytic / reference invariants          blocking
Layer 2  complete acquisition build (01 notebook) blocking
Layer 3  simulate + reconstruct                   review evidence
```

State status **per layer**. A GREEN that covers Layers 1 and 2 does not license an end-to-end
claim, and the difference must be written where a reader will hit it.

Layer 3 may be recorded as:

```text
Layer 3: DEFERRED - the validation oracle itself has not been independently promoted
```

An image produced through an unreviewed reconstruction is two unvalidated things agreeing, not
evidence about the sequence. **A missing `02_simulate_and_reconstruct.ipynb` is not automatically
a task** — see `RadialReadout`, whose Layer 3 is deferred by decision until the non-Cartesian
reconstruction path is independently reviewed with two real consumers.

## What is part of the skill, and what is not

The skill is small and stays small. Three tiers, and things do not drift between them:

```text
.claude/skills/module-mining/
    SKILL.md                workflow, reasoning rules, stop rules -- stable
    reference/, templates/  the candidate data contract

tools/module_mining/skill_v0/migrated/
    tse.yaml radial.yaml gre3d.yaml
                            CALIBRATION FIXTURES ONLY

tools/module_mining/plans/<candidate>/
    ordinary evidence and decisions -- NOT promoted into the skill
```

The three fixtures earn their place by proving the schema can represent three **qualitatively
different** accepted pilots: a promotion judged by event identity, a reimplementation judged by
trajectory agreement, and a reimplementation with no executable reference for one of its modes.

**A new candidate does not become a fixture because it was converted.** It becomes one only if it
exposes a **reasoning class or failure mode the existing suite does not cover** — a new acceptance
criterion, a stop outcome never exercised, a licence or independence shape that breaks the
`Evidence` record. Otherwise its evidence lives in `plans/<candidate>/` like any other candidate's
and the skill does not grow.

> **The skill grows by new reasoning classes and failure modes, not by the number of sequences
> processed.**

A calibration suite that accumulates one entry per conversion stops being a calibration suite and
becomes a second copy of the plans directory.

### What may and may not be written here

This file holds **generic** workflow, evidence, ownership and stop conditions. Rules of this shape
belong here:

```text
if no independent physical contract exists
    -> NO_NEW_MODULE / RED - WRAPPER_ONLY may be appropriate
if the physical boundary is unresolved
    -> YELLOW, rather than a forced extraction
if existing Modules already express the physics
    -> distinguish composition coverage from a missing abstraction
if evidence is non-independent or the oracle is untrusted
    -> do not promote merely because outputs agree
```

Mappings of this shape do **not**, and must be removed on sight:

```text
<family> -> reject
<family> -> defer
<family> -> NEW_KERNEL
```

A settled pilot may be cited as a **worked example** of a reasoning class — that is what the
fixtures are for. A family that has not yet had its fine scan may not have an expected verdict
recorded anywhere the skill reads.

### The prospective-run guardrail

> **Candidate-specific expected outcomes from coarse scans must not be encoded into Skill decision
> logic before a prospective run. The Skill may receive the evidence corpus and the open
> questions, but not the expected terminal classification.**

A coarse scan legitimately forms expectations; that is its job. Those expectations live in the
coarse-scan report, and a prospective run must be able to contradict them. A run that begins
knowing the answer is a rehearsal, and it demonstrates nothing about refusal.

This also grades the evidential value of a result. A run that reaches `NO_NEW_MODULE` on a family
whose coarse scan already expected it is a **regression control** — useful for checking the
workflow still stops, worthless as proof that it can. Only a stop reached on the run's own
evidence, where the outcome was genuinely open beforehand, demonstrates refusal capability.

## A candidate can be blocked by evidence rather than by physics

The pipeline puts evidence before contract for a reason, and the reason is not only ordering: **a
candidate can fail at the evidence stage while its physics remains entirely open.**

Check early, before any contract analysis:

```text
is there an executable reference at all?
is it permissively licensed, or copyleft / EULA-bound?
if copyleft: it can establish that a family exists and what its physics is.
             it cannot serve as acceptance evidence for an implementation you may not copy from
             and whose construction shares the arithmetic you would be testing.
are the executable references independent of each other, or one group's house style?
would simulating against the same oracle the reference used prove anything?
```

When the answer is that no usable acceptance evidence exists, the outcome is
`YELLOW - REFERENCE_NOT_INDEPENDENT` (or `VALIDATION_ORACLE_UNTRUSTED` where the instrument is the
problem). Record it as **evidence-blocked, not rejected** — the physics question stays open and
the candidate can be re-run if a permissive reference appears. Quietly dropping it loses the work;
recording it as `NO_NEW_MODULE` claims evidence nobody has.

`status: UNDECIDED` exists for exactly this record: a YELLOW still needs a status, and every other
value asserts an action nobody took.

## Beware a family defined by its waveform silhouette

A grouping criterion of the form *"the same N pulses keep appearing"* is a **discovery heuristic,
not a boundary.** Pulse counts and flip patterns are what a coarse scan can see cheaply, which is
why families get grouped that way and why the grouping must be re-tested before it is trusted.

The test: for each use, what value could somebody get wrong, and is it the *same* value? If the
uses differ in what is stored, what the intervening gradients are for, or how the preparation
couples to what follows, then what they share is a **coherence pathway or a pulse count**, and
neither is something a module can own. A module owns a value; spin physics is not a value.

Every boundary that has held so far was grouped by a shared **physical solve** — a crusher window
satisfying three axes, a trajectory, a signed moment coupling — not by a waveform shape.

## Standing constraints

- Package-level Module implementations are **Python-only**. MATLAB is reference and validation
  evidence.
- Do not restructure the compiler or `LogicBlock`.
- Do not create a wrapper because a source sequence has a name.
- Do not treat one reference implementation as unquestionable truth.
- **Do not claim physical equivalence from source inspection alone.**
- Verify the comparator before trusting it against a candidate. Its record so far is four wrong to
  the candidates' zero; every one was caught by measuring something whose answer was already known.
