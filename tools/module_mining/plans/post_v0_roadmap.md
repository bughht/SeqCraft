# Post-v0 roadmap

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-20, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../README.md) for which one to edit.

**v0 is a baseline, not a final theory.** Module Mining Skill v0 merged as PR #31 and is frozen:
the workflow, the schema, the vocabulary and the three prospective calibrations. Everything below
starts from that merged `main`.

The governing constraint on all of it:

> Future Skill changes are allowed and expected, but must be **evidence-driven**. Do not pre-invent
> a requirement IR, an augmentation plugin system or a generic constraint solver before real
> joint-solve cases require one.

---

## Phase A — validation-debt cleanup

**Goal:** recover what is still useful from PR #23, independently validate the shared
non-Cartesian path against **two** real consumers, and complete Radial's and Spiral's Layer 3.

```text
PR23 capability reconciliation        pr23_reconciliation.md -- one disposition per capability
non-Cartesian reconstruction review   reduce to the smallest contract two consumers justify
Radial / Spiral Layer 3               four notebooks, each with a different job
```

Architectural boundaries held throughout: no compiler change, no `LogicBlock` change, no
`GRESpiral2D`/`GRERadial2D`, **no reconstruction in `src/`**, no helper promotion without a second
real consumer, no Spiral `echoes > 1`.

The division of labour that decides where reconstruction lives:

```text
SeqCraft      builds pulse sequences
MRzeroCore    simulates them
SigPy         reconstructs the resulting data
```

Reconstruction exists here to teach how an emitted sequence is consumed, to show the expected
physical consequence, and to supply Layer-3 review evidence. It is not part of the
sequence-programming API, so it stays under `examples/`.

**SigPy is the implementation backend, not proof that our adapter is correct.** What this project
owns is the conversion — emitted trajectory → physical k and sample timing → SigPy's coordinate
convention — and that adapter needs validation that does not call SigPy.

## Phase B — batch mining with the frozen v0 Skill

Use v0 **unchanged** as a tool. A broader coarse scan over the registered corpora in
`sources.yaml`, producing a **candidate map**, not implementations.

Classify with the existing coverage vocabulary: `DIRECT_SHIPPED_DUPLICATE`,
`DEGENERATE_CASE_OF_EXISTING_MODULE`, `COMPOSITION_COVERED`, `NOTEBOOK_ONLY_EXISTING`,
`PRIMITIVE_COMPOSITION_COVERED`, or unresolved.

Answers: what families exist; what SeqCraft already expresses; where a reusable physical contract
plausibly remains; which candidates expose new reasoning classes. **No expected fine-scan outcome
may be encoded into skill logic** — the prospective-run guardrail applies to batch output too.

## Phase C — fine scans chosen from evidence

Select from the backlog for **learning and architectural value, not module count**. Run the v0
pipeline: evidence → physical contract → ownership → acceptance claim → implementation only when
justified. `NO_NEW_MODULE`, `YELLOW` and `RED` remain valid results.

## Phase D — architecture-stress candidates

Chosen because they do not all look like ordinary leaf extraction.

| | why it is interesting |
|---|---|
| **bSSFP** | already identified as `BalancedSSFPTR`; balanced moments, RF phase cycling, steady state, sequence-wide timing ownership. The conventional case |
| **DREAM / stimulated echo** | coherence-pathway semantics and cross-component timing — does the abstraction stay a physical contract or become a second sequence DSL? |
| **Flow compensation** | the first real **joint-solve**: base moment requirements + an M1 target + a timing window + hardware limits, solved together |
| **Diffusion encoding** | a deliberate stress test: b-tensor target coupled to TE, refocusing timing, imaging moments, crushers, available windows |
| **Velocity / motion encoding** | a second joint moment-solve, which is what makes the class generalisable rather than one case |

## Phase E — Skill evolution, when evidence requires it

v0 reasons well about reusable contracts, leaf/kernel/composition, reference evidence, validation
and compiler escalation. **It has not been calibrated on a family where correctness requires base
requirements and augmentation requirements to be solved together before waveform realisation.**
Flow compensation, diffusion and velocity encoding may expose exactly that.

When a real candidate shows v0 cannot express the correct abstraction:

```text
1. record how v0 behaves
2. document why the ownership / candidate model is insufficient
3. compare more than one real joint-solve case where possible
4. derive the SMALLEST generic rule or vocabulary extension
5. update the Skill in a separate follow-up PR
```

Do not pre-commit to `RequirementIR`, `AugmentationPlugin` or `ConstraintSolver`. The durable
principle is unchanged and is what any such abstraction must be reconciled against:

> Leaves know intrinsic physics, kernels own cross-leaf coupling, imaging owns acquisition policy,
> and the compiler owns Pulseq legality.

See also `repetition_augmentation_design_note.md`, which records the direction without proposing
an API.

---

## PR decomposition

The roadmap is one document; implementation stays reviewable.

```text
Post-v0 PR A     Phase A, as PR #32.  A1 (contract + adapter validation + reconciliation)
                 is pushed and open for early review; A2 (the four notebooks and the
                 Layer-3 evidence updates) lands on the same branch and the PR does not
                 merge until both are in.
Post-v0 PR B     batch scan / candidate backlog
Post-v0 PR C+    selected fine scans
Skill evolution  separate PR, only when new evidence requires it
```

Each phase starts from merged `main` on its own branch. If future evidence challenges a v0 rule,
that rule changes in a follow-up PR rather than by reopening the v0 calibration.
