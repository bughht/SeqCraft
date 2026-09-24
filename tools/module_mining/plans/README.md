# Module mining plans — what is current, what is evidence, where to look

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-20, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../README.md) for which one to edit.

```text
plans/
    current/
        2026-09-24_repetition_physical_design_architecture.md
                                                THE GOVERNING ARCHITECTURE for the rest of
                                                Stage C.  Where MRI timing design lives, what
                                                stays local, what becomes shared, and the
                                                schedule-outside / realization-inside split
        2026-09-24_repetition_level_joint_waveform_design.md
                                                the direction it settles.  Base requirements plus
                                                augmentation requirements, solved together at the
                                                repetition level, so that velocity encoding and
                                                flow compensation are each defined once instead
                                                of per sequence family.  Opened the revisit that
                                                process/repetition_augmentation_design_note.md
                                                asked for
        2026-09-23_post_v0_implementation_to_authoring_roadmap.md
                                                the implementation roadmap it interrupts.  Still
                                                current for T2Prep, VelocityEncode and bSSFPTR;
                                                its Stage C is what the design above revisits
        roadmap.md                              the post-v0 roadmap it continues.  Still current
                                                for everything before the queue

    process/
        fine_scan_playbook.md                   how a candidate is scanned, in depth
        repetition_augmentation_design_note.md  a direction for joint-solve augmentations,
                                                deliberately not an API

    archive/                                    completed records, date-prefixed so the order
                                                is visible from the listing
        2026-09-18_original_mining_plan.md      the plan the project started from
        2026-09-18_coarse_scan.md               the original corpus scan -- evidence, not a plan
        2026-09-18_radial_fine_scan_proposal.md the Radial fine scan, as planned before its answer
        2026-09-19_two_pilot_retrospective.md   TSE against Radial, before the prospective runs
        2026-09-19_v0_closeout.md               v0 close-out; its FORWARD plan is superseded
        2026-09-20_pr23_reconciliation.md       one disposition per PR #23 capability
        2026-09-21_batch_corpus_scan.md         the batch corpus scan, against merged main
        2026-09-22_skill_v1_synthesis.md        what C1-C3 justified changing in the Skill.  Its
                                                five recommendations are APPLIED -- the Skill is
                                                v1, and diffusion/ is the first run under it

    one folder per candidate                    the record and its reasoning
        tse/  radial/  gre3d/  steam/  saturation/  spiral/  t2prep/  bssfp/
        flow_moments/  -- two records from one scan, plus a gap note
        diffusion/     -- the only one that did not stop at a decision: DiffusionSEPrep exists
```

Three directories and nothing else: **`current/`** is what to follow, **`process/`** is how the
work is done, **`archive/`** is what was decided and when. A document moves to `archive/` when it
is complete, not when it is wrong — several of them are still cited as evidence, and the audit
below says which.

Domain evidence — the curated literature cards, one per fine-scan candidate — is a separate
evidence class from the code corpus and lives in `../domain_evidence/`, beside `../sources.yaml`.

The Skill itself is `.claude/skills/module-mining/`. Its reports — the conformance run, the
generalization report, the MRzero coarse scan and the three-case prospective retrospective —
are in [`../skill_v0/`](../skill_v0/).

---

## Audit

| document | class | note |
|---|---|---|
| `current/2026-09-23_post_v0_implementation_to_authoring_roadmap.md` | **CURRENT** | the active direction: cash out the approved queue, then `sequence-authoring` v0. Not a replacement for the candidate records, which stay the detailed physical specs |
| `current/roadmap.md` | **CURRENT** | the post-v0 roadmap. Its *approved queue* section is superseded by the document above; everything before it still holds |
| `archive/2026-09-22_skill_v1_synthesis.md` | **CURRENT** | the Phase E synthesis across C1, C2 and C3. Its five recommendations were applied in Skill v1; the document stands as the reasoning behind them, and every change it proposed is still costed here |
| `archive/2026-09-21_batch_corpus_scan.md` | **CURRENT** | the phase B scan: what the corpus holds after PR #31 and #32, what SeqCraft already covers, and the candidate backlog. Current *evidence*, not a second roadmap — it recommends a Phase C shortlist and the roadmap records the decision |
| `process/fine_scan_playbook.md` | **DURABLE_REFERENCE** | the process in depth; the Skill is the operational layer over it |
| `process/repetition_augmentation_design_note.md` | **DURABLE_REFERENCE** | open architecture, deliberately not an API (a Chinese translation sits beside it in the workspace mirror) |
| `archive/2026-09-20_pr23_reconciliation.md` | HISTORICAL_UNIQUE | the disposition record; nothing else carries it |
| `archive/2026-09-19_v0_closeout.md` | HISTORICAL_UNIQUE | the close-out reasoning is unique; **its forward roadmap is superseded**, and the file says so at the top |
| `archive/2026-09-18_coarse_scan.md` | HISTORICAL_UNIQUE | **scan evidence**, not an outdated plan — the corpus, the licences and the per-family findings are cited by later candidates |
| `archive/2026-09-18_original_mining_plan.md` | HISTORICAL_UNIQUE | where the project started. Not rewritten: the audit found no reason to |
| `archive/2026-09-18_radial_fine_scan_proposal.md` | HISTORICAL_UNIQUE | the Radial questions *as posed before the answer*, which is part of what makes that scan prospective. Superseded as a plan, retained as evidence |
| `archive/2026-09-19_two_pilot_retrospective.md` | HISTORICAL_UNIQUE | not superseded by the three-case retrospective, which covers different candidates, asks a different question and cites this one |
| `handover_2026-09-18.md` | **SUPERSEDED — removed** | a session handover: start a Radial fine scan, its open questions, and a candidate order. Every question is answered in `radial/`, and the order is superseded by `current/roadmap.md`. No unique reasoning; recoverable from git |

**Nothing was removed for being old.** One file was removed for being wholly answered elsewhere,
and the two documents whose forward-looking halves are superseded say so at the top rather than
being folded into the live roadmap — close-out reasoning and a current plan are different kinds of
document, and merging one into the other loses both.
