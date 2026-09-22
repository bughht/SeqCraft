# The six artifacts, and the evidence for each field

The schema is small because it was derived, not designed. Three pilots have completed —
`tse/candidate.yaml`, `radial/candidate.yaml`, `gre3d/candidate.yaml` in
`tools/module_mining/plans/` — and the field sets were intersected rather than imagined.

## What actually survived all three

Thirteen field paths appear in every pilot:

```text
id
name
capability          capability.summary   capability.family
status
traffic_light
target              target.layer         target.untouched
evidence
validation
open_questions
```

Everything else is one or two pilots' vocabulary: `disagreements`, `priority`,
`existing_seqcraft` and `architecture` are TSE's; `known_limitations` and
`resolved_during_review` are Radial's; `hazards`, `independence`, `measured_so_far`,
`working_hypothesis` and `found_during_implementation` are GRE3D's.

**The sharpest result is about `validation`.** It is present in all three and *not one of its
sub-keys is*. TSE recorded `levels_run` / `reference_executable` / `experiments`; Radial recorded
`against_reference` / `sweep` / `metamorphic`; GRE3D recorded `non_selective` / `slab_selective` /
`polarity` / `timing`. Three disjoint shapes for three acceptance claims. The schema therefore
requires that `validation` exist and says nothing about its interior. That is §4 of the phase
brief falling out of the data rather than being asserted.

Two defects the intersection also exposed, both fixed below:

- **`physics` vs `working_hypothesis`** — TSE and Radial called the contract `physics`, GRE3D
  called the same thing `working_hypothesis`. One concept, two names, so it looked family-specific
  when it is not. The schema names it once, as `contract`.
- **`role` drifted into prose.** TSE used `primary` / `authoritative-external` / `deferred`;
  GRE3D wrote `"primary, slab-selective -- the only witness to that path"`. A field carrying both
  a role and an independence judgement and a scope qualifier cannot be checked or counted. Split
  into `role` (closed vocabulary), `covers` (which mode or aspect) and `independence`.

---

## 1. Candidate

```yaml
id: gre-3d-tr                 # kebab-case, stable
name: GRE3DTR                 # working name; the boundary decision may still reject the class
capability:
  summary: >                  # one paragraph of physics, no API nouns
  family: kernel              # readout | rf | preparation | encoding | kernel | imaging
target:
  layer: kernel               # same vocabulary; a list when the candidate spans layers
  path: src/seqcraft/modules/kernel/gre_3d_tr.py
  composes: [...]             # existing modules reused
  extended: [...]             # existing files this candidate changes -> triggers rule D
  untouched: [...]            # REQUIRED: what is deliberately not touched
  not_created: [...]          # classes considered and rejected, with the reason elsewhere
reuse_rationale: >            # why this is one reusable contract and not two, or none
```

`target.untouched` survived all three pilots, which is worth noticing: in practice the useful
sentence is as often *what the candidate does not own* as what it does. `Refocusing` appears in
all three untouched lists.

## 2. Evidence

```yaml
evidence:
  - repo: pulseq/pulseq
    commit: 2fd6ab6           # REQUIRED for an executable reference
    path: matlab/demoSeq/writeTSE.m
    license: MIT
    role: authoritative-external
    executable: documentary   # executed | documentary
    covers: slab-selective    # which mode or aspect this witnesses; omit if all
    independence: >           # relationship to the OTHER references, not to the candidate
      shares a house style and identical encoding arithmetic with R1; counts as one witness
```

`role` vocabulary: `primary`, `authoritative-external`, `independent-formulation`,
`independent-implementation`, `supporting`, `architecture-evidence`, `deferred`.

**`independence` is the field the pilots most needed and least consistently had.** Only GRE3D
recorded it, and only in prose at the top level — yet it is what decides whether three references
agreeing is three witnesses or one witness copied twice. Two ported implementations of the same
design are one witness. Required whenever more than one reference is cited.

### `evidence_class` — what kind of evidence this is

Optional, added in v1. `role` says how a reference relates to the *other* references cited for a
candidate, and every value presumes the reference is an implementation. A handbook is none of them.

```text
domain-reference   handbook / authoritative review / classic paper
design-witness     an external implementation
oracle             analytic calculation / simulation / independent measurement
```

A `domain-reference` carries a `citation` and a `curated_card` in place of a `repo`, and the card
lives in `tools/module_mining/domain_evidence/`. The schema enforces the citation and warns about a
missing card, because a reference with neither makes the next reader re-read the source.

**Why it exists:** across three post-v0 fine scans, domain literature changed the outcome three
times by three different mechanisms — supplying a family a single code witness could not establish,
correcting a contract three independent witnesses agreed on, and correcting an architectural
inference. The second is the one that made it vocabulary rather than a note: it means domain
evidence is not a fallback for a thin corpus.

## 3. PhysicalContract

```yaml
contract:
  intrinsic: >                # what this candidate determines from its own parameters
  invariants: [...]           # the measurable statements; tolerances live in invariant_table.md
  varying_parameters: [...]
  modes:                      # only when the candidate has meaningful modes -- rule A
    - name: non-selective
      rf_family: block
      ...
  excluded_caller_policy: [...]   # REQUIRED: what stays with the caller
  excluded_pulseq_artifacts: [...] # block splits, event counts, label usage -- representation
```

The two `excluded_*` lists are the artifact form of the ownership rule. A contract that cannot
name what it excludes has not drawn a boundary. `excluded_pulseq_artifacts` exists because block
count differed 185 vs 350 between two TSE implementations of identical physics — a
representation difference that an equality rule would have called a failure.

## 4. AcceptanceClaim

Written **before** extraction. If it is written after, it describes the code.

```yaml
acceptance:
  claim: >                    # the equivalence or preservation actually asserted
  establishes: [...]          # rule C, first half
  does_not_establish: [...]   # rule C, second half -- REQUIRED, never empty
  criterion: trajectory-agreement   # event-identity | trajectory-agreement |
                                    # lattice-and-semantics | invariant-only
  tolerance: >
  defined_before_extraction: true
```

`does_not_establish` being non-empty is enforced by the validator. A claim with nothing out of
scope has not been scoped.

## 5. ValidationPlan

```yaml
validation:                   # interior deliberately unconstrained -- see above
  generic_measurements: [...]
  family_checks: [...]
  metamorphic: [...]
  parameter_sweep: >
  emitted_inspection:         # rule B -- one entry per promoted physical mode
    - mode: non-selective
      file: examples/gre_3d/seq/gre_3d_nonselective.seq
      reads_as: 2 samples, hard/block, 0.200 ms, gradients on no axis
  layers:
    layer_1: GREEN
    layer_2: GREEN
    layer_3: DEFERRED
    layer_3_reason: >         # REQUIRED when layer_3 is DEFERRED or NOT_APPLICABLE
  compiler_changes: none
```

`layers` is the only part of `validation` the schema requires, because a per-layer status is what
stops a Layer-1 GREEN being read as an end-to-end claim.

## 5b. EvidenceState — what is owed, and what would make it actionable

Required only where something is actually owed: a deferred layer, or a non-GREEN light. A
candidate that owes nothing does not carry an empty block for symmetry.

```yaml
evidence_state:
  human_review: DONE | NEEDED | NOT_NEEDED
  experiment_design_review: DONE | NEEDED | NOT_NEEDED
  scanner_validation: NOT_REQUIRED_FOR_CLAIM | RECOMMENDED | REQUIRED | FUTURE | DONE
  deferred:
    - what: ...
      why: ...
      revisit_trigger: ...        # all three required
```

**`established` and `not_established` are deliberately absent** — they are
`acceptance.establishes` and `acceptance.does_not_establish`, and a second copy would drift. What
this block adds is the part no field carried: **each deferral's trigger**.

> A deferral without a trigger is an omission with better prose.

### Three kinds of human evidence, not one bucket

They differ in who can do them and what they cost, so they are separate fields:

| | |
|---|---|
| `human_review` | inspect emitted events, timing, trajectory, a reconstructed image — a workstation afternoon |
| `experiment_design_review` | decide what phantom, protocol and observable would actually test the claim |
| `scanner_validation` | run the `.seq` on hardware and inspect the physical outcome |

### Scanner work is claim-driven

**There is no rule that every GREEN module is eventually scanner-tested**, and one must not be
added. The question is whether hardware is needed for *this* claim. For many leaves an
independent reference comparison is more direct evidence than an image — a sample-for-sample
trajectory match, or a carrier offset read off the compiled file, beats a picture that would
confirm the same thing indirectly through a phantom.

Hardware becomes valuable when the uncertainty is eddy currents, gradient delay, RF hardware
behaviour, B0/B1 imperfection, SAR or duty cycle, receiver behaviour, actual spoiling efficiency
or a contrast outcome — or when no sufficiently independent reference exists, or simulation cannot
observe the failure, or a legal sequence can still produce the wrong physical signal.

### Where it lives when the record is closed

A candidate whose own record is immutable takes an **addendum**: a file with `id`, `applies_to`
and `evidence_state`, beside the record it points at. The validator recognises it by the absence
of `capability` and checks only that block.

### The cross-candidate view is generated

```bash
python tools/module_mining/schema.py --inventory tools/module_mining/plans/*/candidate.yaml \
    tools/module_mining/plans/*/evidence_state.yaml
```

Read from the records, never maintained beside them — a hand-written inventory is a second copy
of the same facts, and the copy is the one that goes stale. Validation records grow by meaningful
unresolved claims, not by duplicating test results.

## 6. Decision

```yaml
status: NEW_KERNEL            # NEW_LEAF | NEW_KERNEL | NEW_IMAGING | PROMOTE_NOTEBOOK |
                              # EXTEND_EXISTING | NO_NEW_MODULE | UNDECIDED
traffic_light: GREEN          # GREEN | APPROVED_FOR_IMPLEMENTATION | YELLOW | RED
                              # (+ reason code -- reference/outcomes.md)
reason_code: null             # required unless GREEN
decided_by: human             # ALWAYS human; the skill recommends
open_questions: [...]         # survived all three pilots; an empty list is suspicious
```

**`UNDECIDED` is a real status, not a placeholder.** A YELLOW record still needs one, and every
other value asserts an action nobody took -- pairing a provisional `NEW_LEAF` with a YELLOW puts a
decision in the record that was never made.

`status` and `traffic_light` are separate because they answered different questions: TSE was
`NEW_KERNEL + PROMOTE_NOTEBOOK` / `GREEN`, and earlier in its life `MORE_EVIDENCE_REQUIRED` /
`YELLOW` with the physics already settled and only the boundary open.

---

## One scan may produce more than one record

A fine scan asks one question; the answer may be that the candidate has more than one owner. The
flow scan asked whether velocity encoding and flow compensation are one contract, and found one
requirement model with two owners — `NEW_LEAF` for the encoding, `EXTEND_EXISTING` for the
compensation. A record holds one `status`, so forcing them together would have meant choosing
which half to misreport.

**Group them in one folder with a shared findings document.** That is the whole convention. There
is deliberately no `scan_id` field: it would be invented from one case, which is the rule the
schema sets for itself everywhere else.

## What the schema deliberately does not do

- It is **not** a Pydantic model or a dataclass hierarchy. It is a documented YAML shape with a
  validator that reports rather than constructs. Freezing these into classes before a fourth
  candidate would harden the vocabulary of exactly three.
- It does not constrain `validation`'s interior.
- It does not require `hazards`, `measured_so_far`, `disagreements` or `priority`, all of which
  earned their place in one pilot and would be noise in another. Use them freely; the validator
  does not mind extra keys.
