# Retrospective conformance: Skill v0 against the three known cases

Skill v0 was run backwards over `TSE/FSE`, `RadialReadout` and `GRE3DTR` before being given
anything new. The question is not whether it reproduces the old prose — it is whether the
workflow **recovers the same reasoning structure and the same safety boundaries**.

No production Module was modified and no reference was re-validated. Two things were executed:

```bash
python tools/module_mining/schema.py tools/module_mining/plans/{tse,radial,gre3d}/candidate.yaml
python tools/module_mining/schema.py tools/module_mining/skill_v0/radial_migrated.yaml
```

---

## 1. The nine questions

`Y` = the workflow forces the question and the pilot's answer is recoverable from it.
`P` = partially — the workflow asks, but something outside it supplied the answer.
`n/a` = the question does not arise for that candidate.

| | TSE/FSE | Radial | GRE3D |
|---|---|---|---|
| recovers the important questions | Y | Y | Y |
| identifies the correct layer/boundary | Y | Y | Y |
| states the physical contract understandably | Y | Y | Y |
| chooses a compatible acceptance claim | Y | Y | Y |
| notices non-independent evidence | Y | **Y** | P |
| states what the comparator does not prove | Y | Y | Y |
| identifies shared-leaf impact | **Y** | **Y** | **Y** |
| requires emitted-sequence inspection | n/a | n/a | **Y** |
| allows a deferred Layer 3 | n/a | **Y** | n/a |

Bold marks a cell where the workflow produces something the pilot **did not** have at the time.

### Boundary and layer

All three land where the humans landed, and for the reason the humans gave. The deciding test —
*if a component cannot determine the correct value using only its own physical contract and
parameters, it should not own that event* — separates every disputed case in the record:
`Refocusing` declining the readout-axis lobe pair, the three-axis window solve staying kernel-local,
`RadialGRETR` not existing because an angle schedule is caller policy, `GRE3DTR` being a sibling of
`GRE2DTR` rather than a wrapper around it.

The workflow also reproduces the *shape* of TSE's history, which is the harder test. TSE sat at
`MORE_EVIDENCE_REQUIRED / YELLOW` with its physics fully settled and only its boundary open.
Skill v0 represents that directly — `status` and `traffic_light` are separate fields, and
`YELLOW - PHYSICAL_BOUNDARY_UNCLEAR` is exactly that state — rather than forcing a single verdict
that would have had to be either wrong or premature.

### Acceptance claims

Each pilot's claim is expressible and none of them collapses into another:

```text
TSE       event-identity          promoting an existing implementation
Radial    trajectory-agreement    two legitimate constructions share no events
GRE3D     lattice-and-semantics   the reference has no DC sample where CartesianLine puts one
```

`criterion` is a per-candidate field rather than a tool-level setting, which is what keeps §4's
rule enforceable: the skill helps *choose* the claim and cannot impose one.

GRE3D exposes the case that matters most. Its slab-selective path has **no executable reference at
all**; the claim is carried by signed-moment arithmetic and three degenerate limits. A workflow
that required reference agreement would have had to either refuse the candidate or accept a
reference that does not exist. Skill v0 records `evidence[].covers`, so "primary, slab-selective —
the only witness to that path" becomes a structured fact rather than a sentence inside a role
string.

---

## 2. What the workflow catches that the pilots did not have

### Shared-leaf impact — all three, and one of them is a real gap

Rule D fires whenever `target.extended` is non-empty. It is non-empty for all three:

| | extended | had a dependency map? |
|---|---|---|
| TSE | `cartesian_line.py` (echo geometry) | no |
| Radial | `cartesian_line.py` (`adc_samples_divisor` refusal) | **not even declared** |
| GRE3D | `excitation.py` (`rephaser_area_per_m`, `build(rephase=)`) | prose, in `mode_contract.md` |

Radial is the interesting one. Its `candidate.yaml` records the `CartesianLine` change under
`resolved_during_review` as narrative and never declares `target.extended`, so nothing could have
fired. The migrated record declares it, and the resulting map is in
`skill_v0/radial_migrated.yaml`: four direct consumers, two transitive, one notebook class.

The map also reproduces the trap that rule D exists for. Searching for the *name* `CartesianLine`
finds ten source files; searching for **construction** finds four. `epi_2d.py`, `excitation.py`,
`refocusing.py` and `_support.py` mention it only in docstrings. The inflated list is worse than
useless — it would have sent a reviewer to check a module that composes nothing.

And the classification is not cosmetic. Every genuine consumer is
`NEW_EARLY_REFUSAL_FOR_PREVIOUSLY_INVALID_INPUT`, which is only safe to assert after checking
whether any of them had been *relying* on passing an illegal sample count through to the compiler.
`examples/fse_2d` runs HASTE at `partial_fourier=0.625` precisely because 0.6 of 128 is 77 — a fact
the notebook stated in prose and nothing enforced. That is the sentence rule D is designed to make
somebody write down.

### Emitted-sequence inspection — GRE3D

Rule B fires for GRE3D because `contract.modes` is non-empty, and the validator refuses a GREEN
with declared modes and no `emitted_inspection` entry. Applied to the original mistake:

> non-selective mode, inherited 3 ms shaped sinc, selection gradient removed

the mode table's **RF family** row differs between modes while the **reference evidence** row says
every official Pulseq 3D reference uses a block pulse. The contradiction is visible in the table
before any code exists, and rule A's prohibition — *a mode may not be "another mode minus one
event" unless reference evidence establishes that equivalence* — names the exact reasoning that
produced it.

**Honest limit:** this is a claim about a workflow, made after the fact, by the same process that
made the mistake. What can be said without hindsight is narrower and still worth something: the
emitted file reads as `3000 samples, shaped / soft, 2.999 ms, gradients on no axis`, and no
invariant in the comparator looks at pulse shape. Rule B routes around the whole class of error by
reading the artifact through a different path than the one that built it, whether or not rule A
would have caught this instance.

### Deferred Layer 3 — Radial

`validation.layers` is required, `layer_3_reason` is required when Layer 3 is deferred, and
`VALIDATION_ORACLE_UNTRUSTED` is a first-class YELLOW. So the Radial position is representable
exactly as decided: Layers 1 and 2 GREEN, Layer 3 `DEFERRED` because the oracle has not been
independently promoted.

This is the one place the workflow actively prevents a plausible wrong move. A pipeline that
treated a missing `02_simulate_and_reconstruct.ipynb` as an open task would generate work whose
correct outcome may be that it is never done.

---

## 3. What the migration exercise cost

`skill_v0/radial_migrated.yaml` re-expresses the Radial pilot in the new shape. It validates
clean, and **nothing the pilot knew was lost**. Three things had to move rather than copy, and each
move is a finding:

1. **`resolved_during_review` → `target.extended` + `dependency_impact`.** The narrative said what
   happened; the structured form says who is affected. Only the second can be checked.
2. **`role` split into `role` + `covers` + `independence`.** Radial's
   `role: decides full-spoke vs centre-out` is a scope note wearing a role's clothes. GRE3D was
   worse — `"primary, slab-selective -- the only witness to that path"` packs a role, a scope and
   an independence judgement into one string that nothing can count.
3. **`note` → `independence`.** All three pilots *made* the independence judgement and none had a
   field for it. TSE: "byte-identical to this file modulo the import name." Radial: "the same
   design as write_radial_gre.py, ported — not independent evidence." GRE3D put it at record level.
   The judgement was always there; it was never checkable, and the number of witnesses is exactly
   the thing a reader will otherwise over-count.

## 4. Validator output on the historical records

Run as-is, unmodified — this is the gap list, not a defect list. These records predate the rules.

```text
tse/candidate.yaml       3 errors:  evidence independence (3 of 4 silent)
                                    validation.layers missing
                                    dependency_impact missing (extends cartesian_line.py)
radial/candidate.yaml    2 errors:  evidence independence (4 of 5 silent)
                                    validation.layers missing
gre3d/candidate.yaml     2 errors:  validation.layers missing
                                    dependency_impact missing (extends excitation.py)
                         + role-vocabulary warnings on all 5 references
```

GRE3D passes the independence check that the other two fail, because it has a record-level
`independence` block — the only pilot that wrote one. That is the correct relative result and it
is a small piece of evidence that the check measures something real rather than penalising age.

**The historical records were deliberately not migrated.** Rewriting them would destroy the
evidence that the schema was derived from the pilots rather than imposed on them, and the whole
argument for the field set rests on that derivation.

---

## 5. Where the workflow would still not have saved us

Stated plainly, because a conformance report that finds only successes has not been run properly.

- **The comparator was wrong four times and the candidates zero.** Skill v0 carries this as a
  standing constraint and points at `tests/module_mining/`, but nothing in the workflow *derives*
  the need to verify an instrument against known answers. It is transmitted as a rule, which is
  weaker than a rule that falls out of the process.
- **TE short by exactly the winder** was caught by comparing the compiled echo time against the
  reported one. The workflow asks for family-specific checks; it does not tell you that
  *design-against-itself* is the failure mode, so a candidate could satisfy every listed step with
  a check that compares the module to its own arithmetic.
- **The z spoiler colliding with the partition rewinder** was caught by the compiler's slew limit,
  i.e. by SeqCraft, not by this workflow. That is fine — but it means the ladder's Layer-1 coverage
  depends on the package being strict, and a family whose failures are legal-but-wrong gets less
  protection than this one did.
- **Rule A's retrospective success is unfalsifiable** and is marked as such above.

None of these is a reason to widen the schema. Two of them are reasons to keep `open_questions`
mandatory and non-empty, which the validator warns about.
