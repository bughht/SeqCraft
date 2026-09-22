# Module mining plans — what is current, what is evidence, where to look

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-20, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../README.md) for which one to edit.

```text
CURRENT ROADMAP
    post_v0_roadmap.md                       the only active forward plan
    post_v0_batch_scan.md                    the phase B corpus scan, against merged main
    phase_e_synthesis.md                     what C1-C3 justify changing in the Skill.  Its
                                             five recommendations are APPLIED -- the Skill is
                                             v1, and diffusion/ is the first run under it

DOMAIN EVIDENCE                              ../domain_evidence/ -- curated literature cards,
                                             one per fine-scan candidate.  A separate evidence
                                             class from the code corpus in ../sources.yaml

DURABLE PROCESS / OPEN ARCHITECTURE
    fine_scan_playbook.md                    how a candidate is scanned, in depth
    repetition_augmentation_design_note.md   a direction for joint-solve augmentations, no API

CANDIDATE EVIDENCE                           one folder per candidate: the record and its reasoning
    tse/  radial/  gre3d/  steam/  saturation/  spiral/  t2prep/  bssfp/
    flow_moments/  -- two records from one scan, plus a gap note
    diffusion/     -- the only one that did not stop at a decision: DiffusionSEPrep exists

HISTORICAL / RECONCILIATION RECORDS
    pr23_reconciliation.md                   one disposition per PR #23 capability
    v0_closeout_and_roadmap.md               v0 close-out reasoning; its FORWARD plan is superseded
    coarse_scan_2026-09-18.md                the original corpus scan -- evidence, not a plan
    module_mining_plan.md                    the plan the project started from
    next_phase_plan.md                       the Radial fine scan, as planned before its answer
    two_pilot_retrospective.md               TSE against Radial, before the three prospective runs
```

The Skill itself is `.claude/skills/module-mining/`. Its reports — the conformance run, the
generalization report, the MRzero coarse scan and the three-case prospective retrospective —
are in [`../skill_v0/`](../skill_v0/).

---

## Audit

| document | class | note |
|---|---|---|
| `post_v0_roadmap.md` | **CURRENT** | the only forward roadmap. Everything else that points forward is superseded by it |
| `phase_e_synthesis.md` | **CURRENT** | the Phase E synthesis across C1, C2 and C3. Its five recommendations were applied in Skill v1; the document stands as the reasoning behind them, and every change it proposed is still costed here |
| `post_v0_batch_scan.md` | **CURRENT** | the phase B scan: what the corpus holds after PR #31 and #32, what SeqCraft already covers, and the candidate backlog. Current *evidence*, not a second roadmap — it recommends a Phase C shortlist and the roadmap records the decision |
| `fine_scan_playbook.md` | **DURABLE_REFERENCE** | the process in depth; the Skill is the operational layer over it |
| `repetition_augmentation_design_note.md` | **DURABLE_REFERENCE** | open architecture, deliberately not an API (a Chinese translation sits beside it in the workspace mirror) |
| `pr23_reconciliation.md` | HISTORICAL_UNIQUE | the disposition record; nothing else carries it |
| `v0_closeout_and_roadmap.md` | HISTORICAL_UNIQUE | the close-out reasoning is unique; **its forward roadmap is superseded**, and the file says so at the top |
| `coarse_scan_2026-09-18.md` | HISTORICAL_UNIQUE | **scan evidence**, not an outdated plan — the corpus, the licences and the per-family findings are cited by later candidates |
| `module_mining_plan.md` | HISTORICAL_UNIQUE | where the project started. Not rewritten: the audit found no reason to |
| `next_phase_plan.md` | HISTORICAL_UNIQUE | the Radial questions *as posed before the answer*, which is part of what makes that scan prospective. Superseded as a plan, retained as evidence |
| `two_pilot_retrospective.md` | HISTORICAL_UNIQUE | not superseded by the three-case retrospective, which covers different candidates, asks a different question and cites this one |
| `handover_2026-09-18.md` | **SUPERSEDED — removed** | a session handover: start a Radial fine scan, its open questions, and a candidate order. Every question is answered in `radial/`, and the order is superseded by `post_v0_roadmap.md`. No unique reasoning; recoverable from git |

**Nothing was removed for being old.** One file was removed for being wholly answered elsewhere,
and the two documents whose forward-looking halves are superseded say so at the top rather than
being folded into the live roadmap — close-out reasoning and a current plan are different kinds of
document, and merging one into the other loses both.
