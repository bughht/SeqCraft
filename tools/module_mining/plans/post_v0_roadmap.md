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

## Phase A — validation-debt cleanup — **COMPLETE** (PR #32, 2026-09-20)

```text
A1  the shared non-Cartesian contract, validated against a dense-DFT oracle, plus the PR23
    reconciliation                                                                     DONE
A2  four notebooks -- gre_radial_2d/02, gre_spiral_2d/02, se_spiral_2d/01 and /02 --
    plus the Radial and Spiral Layer-3 evidence updates and the reconciliation closure  DONE
```

Radial Layer 3 and Spiral Layer 3 are both **GREEN**, claim-scoped in
`radial/evidence_state.yaml` and `spiral/candidate.yaml`. The reconstruction stayed under
`examples/`: two consumers did not trigger promotion, which is the decision rather than an
omission. Two sequence-side defects were found by composing on top of `SpiralReadout` after its
own suite was green, and both were fixed with regressions before the work that found them
landed — `spiral/findings.md` §9 and §10.

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

## Phase B — batch mining with the frozen v0 Skill — **COMPLETE** (2026-09-21)

The scan report is [`post_v0_batch_scan.md`](post_v0_batch_scan.md). Nine registered corpora read,
the coverage map redrawn against merged `main`, and a candidate backlog with reasoning-class
annotations. **No Module was implemented and the Skill was not edited**, including where v0's
vocabulary was awkward — those are recorded as Phase E evidence in the report's §6.

The result that matters is not a candidate. A published 4D-flow implementation does not build its
velocity-encoding gradients at all: it states a residual **moment requirement**, relative to what
the imaging gradients already contribute, and calls an external constrained optimiser. The
requirement-before-realisation regime is now observed rather than hypothesised, and no SeqCraft
Module can express it.

Proposed Phase C shortlist, chosen to be complementary rather than likely GREEN: **T2
preparation** (the conventional control), **bSSFP** (two independent witnesses that disagree about
the boundary), and **velocity encoding with flow compensation scanned together** (can a
requirement be expressed at all). Diffusion is deliberately held back behind the third.

The report also recommends one ordering change — a bounded Skill pass between C1/C2 and C3 — and
says plainly why the stricter alternative is defensible. That is a decision to take, not one taken.

## Phase C and E — **fine scans COMPLETE, synthesis READY FOR REVIEW** (2026-09-22)

Three fine scans against frozen v0, four recommendations, **nothing implemented**:

```text
C1  T2 preparation      plans/t2prep/        NEW_LEAF        / APPROVED_FOR_IMPLEMENTATION
C2  balanced SSFP       plans/bssfp/         NEW_KERNEL      / APPROVED_FOR_IMPLEMENTATION
C3  velocity encoding   plans/flow_moments/  NEW_LEAF        / APPROVED_FOR_IMPLEMENTATION
    flow compensation   plans/flow_moments/  EXTEND_EXISTING / APPROVED_FOR_IMPLEMENTATION
```

Domain literature became a distinct evidence class during C1 and changed the outcome of all three
scans by three different mechanisms — supplying a family, correcting a contract three witnesses
agreed on, and correcting an architectural inference. Cards in `../domain_evidence/`.

The Phase E recommendation is [`phase_e_synthesis.md`](phase_e_synthesis.md). Its headline: the
findings **refine the meaning of the Module layer** and require no change to `LogicBlock`, the
compiler, or the `Module → LogicBlock → Tree → Compiler` architecture. C2's 220 µs stays
undecided pending one ownership question the synthesis poses.

## The approved queue — designed, reviewed, **deliberately unbuilt**

Four designs carry `APPROVED_FOR_IMPLEMENTATION` and none is implemented. That is a decision, not
a backlog that was forgotten.

| | record | status | why it is not built yet |
|---|---|---|---|
| **T2Prep** | `plans/t2prep/` | `NEW_LEAF` | **could be built now.** Held because it would not demonstrate the architecture the Phase C work found, and the next PR is about demonstrating it |
| **VelocityEncode** | `plans/flow_moments/velocity_encoding.yaml` | `NEW_LEAF` | the appended canonical form is implementation-ready; same reason as above |
| **FlowComp** | `plans/flow_moments/flow_compensation.yaml` | `EXTEND_EXISTING` | depends on how joint realization is finally expressed |
| **bSSFPTR** | `plans/bssfp/` | `NEW_KERNEL` | same dependency, plus the open question of whether its 220 µs is joint physical design or semantics-preserving compiler fusion |

Building all four now would either dilute the architectural point or produce code likely to be
reshaped immediately afterwards.

**The next PR is `Diffusion: architecture stress test + proof-by-implementation`** — a targeted
fine scan, an architecture decision from it, production implementation, and a Layer 1/2/3
showcase. It is explicitly not a fourth research-only scan. It tests two Phase E hypotheses on a
new physical family: *analytic first with an optimizer fallback*, where the MIT Pulseq witness
inverts the b-factor analytically but leaves trapezoid ramps `TODO` while the copyleft witness
holds the minimum-TE direction; and *joint realization before materialisation*, by asking whether
a joint owner can get everything it needs before child `LogicBlock`s are built.

Afterwards the queue is revisited, in the likely order T2Prep, VelocityEncode, FlowComp, bSSFPTR —
**and not forced through one mechanism if the physical evidence says otherwise.**

## PR decomposition

The roadmap is one document; implementation stays reviewable.

```text
Post-v0 PR A     Phase A, as PR #32.  A1 (contract + adapter validation + reconciliation)
                 is pushed and open for early review; A2 (the four notebooks and the
                 Layer-3 evidence updates) lands on the same branch and the PR does not
                 merge until both are in.
Post-v0 PR B     batch scan / candidate backlog.  COMPLETE -- post_v0_batch_scan.md.
Post-v0 PR C+    selected fine scans
Skill evolution  separate PR, only when new evidence requires it
```

Each phase starts from merged `main` on its own branch. If future evidence challenges a v0 rule,
that rule changes in a follow-up PR rather than by reopening the v0 calibration.
