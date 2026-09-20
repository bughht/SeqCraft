# Migrated records: the same evidence, in the current shape

Two things are kept, deliberately, and neither replaces the other:

```text
historical record   plans/{tse,radial,gre3d}/candidate.yaml
                    = what the workflow actually captured at the time

migrated record     skill_v0/migrated/{tse,radial,gre3d}.yaml
                    = how the same evidence is represented by Skill v0 now
```

**The historical records are unchanged and must stay that way.** The argument for the Skill v0
field set is that it was derived by intersecting those three files; editing them would destroy the
evidence for the thing they justify. The conformance errors they produce are not defects in any
production Module and are not a reason to touch PR #29.

All three migrated records validate:

```bash
python tools/module_mining/schema.py tools/module_mining/skill_v0/migrated/*.yaml
```

```text
gre3d.yaml    ok
radial.yaml   ok
tse.yaml      ok, 1 warning  (writeHASTE.m has no licence recorded -- see §4)
```

## 0. These are calibration fixtures, and the set is curated

The three files here are **not** a knowledge base and must not grow with every converted
sequence. They exist to prove one thing: that the current schema can represent three
qualitatively different accepted pilots.

| fixture | the reasoning class it covers |
|---|---|
| `tse.yaml` | a **promotion**, judged by event identity, with a superseded YELLOW in its history |
| `radial.yaml` | a **reimplementation**, judged by trajectory agreement, with a deliberately deferred Layer 3 |
| `gre3d.yaml` | a reimplementation with **no executable reference for one of its modes**, and a declared mode set that fires rule B |

A future candidate is added here **only if it exposes a reasoning class or failure mode this set
does not already cover** — a new acceptance criterion, a stop outcome never exercised (all of them
currently), or an evidence shape the `Evidence` record cannot carry. Conversion alone is not a
reason. Ordinary candidates keep their evidence in `plans/<candidate>/`.

## 1. Provenance convention

Every field the schema requires carries one of three markers:

| marker | meaning |
|---|---|
| `[VERBATIM]` | copied from the historical `candidate.yaml`, value unchanged |
| `[RECONSTRUCTED]` | assembled from prose elsewhere in that candidate's evidence set; **the source document and section are named inline** |
| `[NOT RECORDED]` | the pilot did not capture this. Stated as absent. Never invented. |

The rule this enforces: a migrated record may **re-shape** what the pilot knew and may **not add**
to it. Where a currently-required judgement exists only as prose, it is marked `[RECONSTRUCTED]`
with its source so a reader can check the reconstruction rather than trust it.

## 2. Field mapping

### Renames — one concept that had three names

| historical | normalized | note |
|---|---|---|
| `physics` (TSE, Radial) | `contract` | same concept |
| `working_hypothesis` (GRE3D) | `contract` | same concept, third name |
| `physics.caller_policy` | `contract.excluded_caller_policy` | now explicitly an *exclusion* |
| `existing_seqcraft.reused` (TSE) | `target.composes` | |
| `target.not_proposed` (GRE3D) | `target.not_created` | |

### Splits — one field carrying three jobs

`role` was free text and had absorbed scope and independence judgements:

| historical `role` | `role` | `covers` | `independence` |
|---|---|---|---|
| `primary, non-selective` | `primary` | `non-selective` | (from the record-level block) |
| `primary, slab-selective -- the only witness to that path` | `primary` | `slab-selective` | "the ONLY witness to the slab-selective path" |
| `architecture evidence (2D)` | `architecture-evidence` | `2D architecture only` | |
| `decides full-spoke vs centre-out` (Radial) | `supporting` | `full-spoke vs centre-out` | |

### Promotions — prose to structure

| was | now | which record |
|---|---|---|
| `evidence[].note` "ported — not independent evidence" | `evidence[].independence` | Radial |
| `evidence[].note` "byte-identical modulo the import name" | `evidence[].independence` | TSE |
| record-level `independence.note` | per-item `evidence[].independence` | GRE3D |
| `resolved_during_review` narrative | `target.extended` + `dependency_impact` | **Radial** |
| `acceptance_criterion.md` (a whole document) | `acceptance` block | GRE3D |
| `validation.reference_executable` R1/R2/R3 | `evidence[].executable` | TSE |
| the validation ladder, recorded in the PR body | `validation.layers` | all three |

## 3. The finding this exercise was for

**Radial's extraction changed `CartesianLine`'s contract and the record never said so.**

The `adc_samples_divisor` refusal — a new, *earlier* refusal of a sample count the receiver cannot
digitise — was recorded under `resolved_during_review` as a paragraph of narrative. `target.extended`
was left empty. So rule D had nothing to fire on: mechanically, that candidate looked like a pure
addition that touched no shared leaf.

`migrated/radial.yaml` declares it, and the map that follows is what would have been produced:

```text
CartesianLine
  -> gre_2d_tr.py      NEW_EARLY_REFUSAL_FOR_PREVIOUSLY_INVALID_INPUT
  -> tse_shot.py       NEW_EARLY_REFUSAL_FOR_PREVIOUSLY_INVALID_INPUT
  -> gre_3d_tr.py      NEW_EARLY_REFUSAL_FOR_PREVIOUSLY_INVALID_INPUT
  -> radial_readout.py INTENTIONAL_BEHAVIOR_CHANGE        (the candidate itself)
  -> gre_2d.py         NO_BEHAVIOR_CHANGE                 (transitive via GRE2DTR)
  -> fse_2d.py         NO_BEHAVIOR_CHANGE                 (transitive via TSEShot)
  -> examples/         NO_BEHAVIOR_CHANGE
```

Consumers were found by searching for **construction**, not for the class name. The name appears
in ten source files; four construct it. `epi_2d.py`, `excitation.py`, `refocusing.py` and
`_support.py` mention it only in docstrings.

The same search was run for `Excitation` when migrating GRE3D, and corrected a wrong assumption in
the first draft of that file: `IRPrep` does **not** compose `Excitation` — it builds its own
inversion pulse. Three construction sites, not four.

## 4. What could not be represented without loss

Reported rather than fixed by changing the schema.

**1. `writeHASTE.m` has no licence and no independence assessment.** It was deferred before any
comparison, so neither was ever established. The record carries `license: null` and an explicit
"not assessed", and the validator warns. The repository-level MIT signal is recorded in the coarse
scan, but that is a repository fact and this field is a file-level one; filling it in would be
inventing evidence of exactly the kind the convention forbids. **This warning is correct and
should stay.**

**2. GRE3D has no metamorphic test.** None was defined. The three degenerate limits play a related
role — the non-selective limit is the selective solve at zero slab rephasing — but they were framed
as invariants, not as transformations, and reclassifying them now would be a claim the pilot never
made. `metamorphic: []`, and a new open question records the possibility.

**3. TSE's dependency map covers the current tree, not the tree at the time.** `gre_3d_tr.py` and
`radial_readout.py` did not exist when TSE was extracted. They are listed, with a note saying so,
because the alternative — silently omitting them — would misrepresent what the map covers. A
retrospective blast radius is a blast radius **now**.

**4. The historical layer statuses are reconstructed, not recorded.** The validation ladder was
defined during the GRE3D pilot; TSE and Radial predate it. The `layers` blocks say what the
recorded evidence actually covers, marked `[RECONSTRUCTED]`. The one that matters — Radial's
Layer 3 `DEFERRED` — was an explicit decision and is recorded as one.

**5. TSE's `status` is compound.** `NEW_KERNEL + PROMOTE_NOTEBOOK` describes two acts: a new
kernel, and a notebook implementation promoted into the package. The validator splits on `+` and
checks each part rather than forcing one label. Flattening it would lose the fact that the pilot
did two distinguishable things.

**6. TSE's superseded YELLOW has no recorded reason code**, because the catalogue did not exist.
It is reconstructed as `PHYSICAL_BOUNDARY_UNCLEAR` — the historical note says the physics was
settled and the boundary was not, which is that code's definition — and marked as reconstructed.
It is kept rather than dropped: it is the more instructive half of the record.
