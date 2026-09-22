# SeqCraft Module Mining — Two-Pilot Retrospective

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../README.md) for which one to edit.

**Date:** 2026-09-19
**Pilots:** TSE/FSE — COMPLETE, GREEN · `RadialReadout` — COMPLETE, GREEN
**Purpose:** decide what hardens into the module-mining workflow and what stays experimental,
**before** a third candidate makes the question harder to answer.

The milestone this closes is not "two modules converted". It is *two deliberately different
conversions succeeded, so we can tell which parts of the process were the process and which were
the candidate.*

---

## 1. The two pilots were different on purpose

| | TSE/FSE | `RadialReadout` |
|---|---|---|
| what it was | **promotion** of a working SeqCraft notebook class | **independent implementation** from external references |
| layer | kernel + imaging | one readout leaf |
| trajectory | Cartesian | non-Cartesian |
| central difficulty | cross-leaf coupling — one crusher window shared by three axes | trajectory geometry and orientation |
| strongest proof | **event identity**: same absolute time, same content hash, at turbo 1, 8, 16 and HASTE | **trajectory identity**: k at every ADC sample, 4.96e-12 1/m against the official reference |
| could it have used the other's proof? | yes, and did | **no** — two legitimate constructions of one trajectory share no events |
| new kernel? | yes, `TSEShot` | **no**, and that was the finding |
| references genuinely independent? | yes — two Python, one MATLAB | **no** — the two official demos are the same design ported; the weight fell on OpenMRF and the fast variant |

That last row is worth keeping. We nearly counted PyPulseq's `write_radial_gre.py` and Pulseq's
`writeRadialGradientEcho.m` as two witnesses. They are one, ported. Checking whether references
are independent is now part of Pass 0, not an afterthought.

---

## 2. What was identical in both, and is therefore the workflow

Eleven steps ran the same way twice, on candidates that had almost nothing else in common:

```text
0  freeze provenance -- repo, commit, licence, and whether the references are independent
1  code archaeology, reading for intent
2  semantic decomposition, with each component classified
3  the invariant table, written BEFORE any API
4  thin adapters: execute the reference, never reimplement it
5  measure the references, and reconcile their protocols
6  state the ownership questions -> HUMAN DECISION
7  extract
8  validate, by a criterion chosen for this candidate
9  parameter sweep
10 metamorphic tests
11 GREEN / YELLOW / RED
```

Three of these earned their place rather than merely occupying it.

**Step 3, the invariant table before the API.** In both pilots it made "did this change anything?"
answerable by something other than the new code's own tests. In Radial it went further: the
rotation-equivariance test was written and passed *against the official reference* before
`RadialReadout` existed, so it cannot encode the module's behaviour. That ordering is now a rule,
not a habit.

**Step 5, protocol reconciliation, produces findings before any comparison runs.** TSE: the
reference's 2 ms refocusing pulse peaks at 116 % of its own `max_b1` and warns, where SeqCraft
refuses; its readout time is fixed in the script body; its derived crusher window is assumed
raster-legal and is not at a 30 µs ringdown. Radial: matching geometry exposed that the two
official demos are one design. None of this needed a comparator.

**Step 6, the human checkpoint, is at the boundary and not at the end.** Both pilots went
*evidence → YELLOW → human ownership decision → extraction → GREEN*. Neither reached GREEN by
measurement alone, and neither needed a human to read waveforms. That is the split the project set
out to get, and it held twice.

---

## 3. What was necessarily different

**The acceptance criterion.** TSE was a promotion, so nothing was allowed to move and event
identity was exactly right. Radial was a reimplementation, so event identity was meaningless and
the trajectory was the contract. Both are valid GREENs. **This is the strongest single argument
against one universal equivalence rule**, and §5 acts on it.

Also candidate-specific: the L4 physics checks, the metamorphic tests, the semantic quantities a
module exposes, the sweep dimensions, and what "physically equivalent" even means. A framework
that pretends these are generic will either be wrong or empty.

---

## 4. The comparator was wrong three times. The candidates were wrong none.

| | what the comparator claimed | what was actually true |
|---|---|---|
| TSE | `kx` at the echo alternating ±1.947 1/m — a textbook odd/even error | the reference has **no DC sample**; `argmin \|kx\|` was choosing between two values half a k-step apart on floating-point noise |
| TSE | both implementations violate the equal-area rule by ~230 1/m | the rule does not apply to the **phase axis** — an encode is *meant* to move k |
| Radial | rotation equivariance out by **5.8e+02** 1/m | the harness stacked bare spokes in one sequence; `k` integrates from the start and resets only at an excitation, so each spoke began where the last ended |

Three for three. Every one looked exactly like a real physics bug, and every one was caught the
same way: **the failing quantity had a known answer**, and the known answer was used to audit the
measurement instead of the candidate.

Two rules follow, and they belong in the workflow rather than in a postmortem:

> **A failed comparator result is evidence to investigate, not proof that the candidate is wrong.**

> **Validate the validator against something whose answer is already known** — a metamorphic
> identity, an analytic value, or a reference implementation — before believing it about something
> whose answer is not.

There is a structural reason this keeps happening, and §7 states it: a comparator that decides
*where the echo is*, or *which axes a rule applies to*, has started modelling the sequence. It
should be measuring claims, not making them.

---

## 4a. Three layers, and what each one catches

The comparator work made it easy to forget that waveform arithmetic is not the whole of
validation. The ladder the hand-built modules were originally developed through is back, and the
distinction is worth keeping explicit:

| | question | status |
|---|---|---|
| **Layer 1** analytic invariants and regressions | does the design satisfy the stated contract? | **blocking CI** |
| **Layer 2** the example compiles and executes | did the emitted sequence preserve it? | **blocking CI** (notebook smoke) |
| **Layer 3** simulate and reconstruct | does the whole acquisition behave like the experiment intended? | review evidence, **not** a gate |

Layer 3 is where axis swaps, sign reversals, orientation mistakes and coverage errors stop being
numbers. It is also where three separate mistakes surfaced during the GRE3D work — a TR that could
not hold a longer arrangement, a comparison missing a spoiler, and a phantom passed unbuilt — none
of which any Layer 1 check would have asked about.

**A rendered image is never the criterion.** `assert image == golden_png` would be a gate that
fails for reasons nobody can read. A small number of numerical assertions is reasonable —
reconstructed shape, energy inside the object, a sequential-versus-combined difference — and the
rest is for a human to look at.

## 5. Category A — stable enough to formalise

These survived two very different candidates unchanged, and should become the durable part:

| | why it is ready |
|---|---|
| **the eleven-step lifecycle** (§2) | ran identically on a kernel promotion and a leaf reimplementation |
| **the provenance record** — repo, commit, path, licence, runnable, role, **independence** | the same eight fields served both; the independence field was earned by Radial |
| **evidence → boundary decision → extraction → validation → GREEN/YELLOW/RED** | both pilots, including the YELLOW that preceded the human decision |
| **the human checkpoint at the boundary question** | twice, neither decidable by measurement |
| **the thin-adapter contract** — execute the reference, never edit the source; the *harness* may intervene for one call | `pypulseq_tse.overridden_opts` patched a hard-coded `Opts`; the radial adapters call `main()` untouched |
| **structured difference reports** — worst value, its location, and the tolerance beside it | "pass" would have hidden all three comparator bugs |
| **the invariant table written before the API, and metamorphic tests validated on a reference first** | the ordering is what makes the later result mean anything |
| **the fine-scan document set** — inventory, findings, invariant table, comparison/sweep report, boundary decision, candidate record | produced twice, read by a human twice |

## 6. Category B — stable frame, family-specific content

The frame is the L0/L2/L3/L4 stack and the report shape. What goes *in* L4 is the family's, and
the framework must let it be:

- **L4 physics checks.** TSE: refocusing spacing, echo-at-midpoint, crusher balance, effective TE.
  Radial: angle, extent, Δk, straightness, in-plane, rotation equivariance. No overlap at all.
- **How the semantic centre of a readout is located.** Already a parameter —
  `measure(centre='kx-zero' | 'k-min')` — because a rotated spoke's `kx` need not cross zero.
  Adding a rule per family is cheap; guessing is silent.
- **Metamorphic tests, semantic quantities, sweep dimensions, and the meaning of equivalence.**

## 7. Category C — do not freeze on two examples

- the exact `ReferenceSequence` Python API;
- the exact `ModuleCandidate` schema — still research YAML, and it should stay there;
- **any single equivalence algorithm**, mandatory event-identity check, or universal L4 check;
- automatic acceptance thresholds;
- batch-conversion policy;
- a shared helper for the "coupled moment window" solve. `GRE2DTR` and `TSEShot` both do it and
  it is still two small pure functions in the kernel that needed them. Two is not a contract.

### `SequenceFingerprint`, redefined

If the name survives, it means:

> a structured collection of deterministic measurements relevant to a **sequence family** and a
> **validation goal**

and specifically **not**:

> a canonical hash proving sequence equivalence in general.

The two pilots' strongest proofs — event identity and trajectory identity — are not the same
measurement, are not reducible to one another, and are each wrong for the other candidate.

---

## 8. The ownership model, extended

Both pilots used it and it held. The last line is new, and it is the one the three comparator bugs
paid for:

```text
leaf       intrinsic waveform/trajectory physics, and the semantic facts derivable from that design
kernel     quantities that require several leaves to be solved together
imaging    acquisition scheduling and policy
compiler   Pulseq legality
tooling    measures the claims the others make -- and does not define their physics
```

Operationally, unchanged:

> If a component cannot determine the correct value from its own contract and parameters, it
> should not own that event.

Two corollaries the pilots added:

**Validation follows the same rule as events.** The layer that *determines* a quantity validates
the constraints on it: `CartesianLine` computes the ADC sample count, so it now refuses a count
the receiver cannot digitise, and `RadialReadout` inherits that by composing the line. The
compiler stays the final backstop; the caller gets the useful error first.

**Arithmetic moves down more often than classes move up.** TSE pushed the readout's echo
asymmetry into `CartesianLine`; Radial reused `partial_fourier` and the sample arithmetic whole,
and pushed the divisor check down too. Two of three extractions *shrank* the new code by making an
existing leaf state what it already knew.

> **A watch item, not a conclusion.** `CartesianLine` is becoming the shared geometric core —
> echo geometry, sample counts, partial Fourier, now consumed by a second readout. That is the
> right direction and it has a failure mode. If a third candidate pushes something down into it
> that a *line* would not naturally know, that is the signal to split the geometry out rather than
> keep adding.

---

## 9. What this changes for the next candidate

`GRE3DTR` remains next, and this retrospective is what it should be started from — not restarted
from scratch. Concretely, it should:

- run the eleven steps in §2 without inventing a twelfth;
- **choose its acceptance criterion explicitly, in writing, before extracting.** It is another
  promotion-shaped candidate near `GRE2DTR`, so event identity is probably right — but "probably"
  is the reason to say so first;
- check reference independence in Pass 0;
- expect the coupled-window question to return, and **still** not extract a shared helper unless
  a third instance defines the contract;
- watch whether 3D phase encoding wants something from `PhaseEncode` that a 2D encode would not
  know — the `CartesianLine` watch item, on a different leaf.

## 10. Open questions for the team

1. Should the fine-scan document set become a template, or stay hand-written? It has been
   hand-written twice with the same six documents.
2. Where should mined evidence live long term? It is currently duplicated: the workspace planning
   repo holds the working copy and `tools/module_mining/plans/` a mirror taken at PR #26.
3. Does the `tools/module_mining/` comparator belong in CI, or stay a local instrument? It has no
   gate today, and it has been wrong three times — which is an argument both ways.
4. How many candidates before `ModuleCandidate` is worth code? The original plan said 20–40
   scanned and 5–8 families; the coarse scan met that, but only two have been *converted*.
