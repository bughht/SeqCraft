# v0 prospective calibration — what three real candidates forced the Skill to learn

Three candidates were run **prospectively**: evidence first, outcome unknown at the start, with
the guardrail that no expected terminal classification could be encoded into skill logic
beforehand. This is not a summary of three implementations. It answers one question:

> What did real prospective use force Module Mining Skill v0 to learn, change, or validate?

---

## 1. The questions v0 started with

After three *retrospective* pilots (TSE/FSE, Radial, GRE3D — all GREEN, all already decided), the
open questions were:

```text
can the workflow stop?                    no stop outcome had ever been exercised
can it extract without over-designing?    every pilot so far was a real module
can it separate a candidate's problems
    from the compiler's?                  never tested against a real compiler claim
does the schema survive a fourth case?    derived from exactly three
```

Each prospective run was chosen to press on a different one.

## 2. STEAM — `UNDECIDED / YELLOW / REFERENCE_NOT_INDEPENDENT`

**The workflow stopped.** It met a genuinely open candidate and declined to extract, which is the
first stop this project has produced on its own evidence.

What it taught:

| | |
|---|---|
| **waveform similarity is a discovery heuristic** | "the same three pulses keep appearing" grouped three uses with three pulse counts, three flip patterns and two incompatible readout couplings. Promotion needs a shared *physical solve*, not a silhouette. |
| **reference independence is load-bearing** | 21 of 26 MRzero playground notebooks import PyPulseq, and the workflow pins `pypulseq==1.4.2post1`. Agreement with them is one design agreeing with itself. |
| **"not enough evidence yet" is a result** | recorded as *evidence-blocked, not rejected*, with a revisit trigger — because the physics question stays open. |
| **a licence can block before physics does** | the four prior candidates all had permissive evidence, so the case had never arisen. |

It also exposed a vocabulary gap: a YELLOW record still needs a `status`, and every value asserted
an action nobody took. `UNDECIDED` was added.

**Why this matters more than a fourth GREEN:** it is the evidence that the Skill is not a
module-producing machine. A workflow that had only ever produced modules would be indistinguishable
from one that always does.

## 3. SaturationPrep — `NEW_LEAF / GREEN`

**A simple concept extracted without over-designing it.** The pressure here was the opposite of
STEAM's: the temptation to widen the candidate until it covered CEST trains and regional
saturation slabs, which `findings.md` §7 refuses explicitly.

What it taught:

| | |
|---|---|
| **reference disagreement on protocol parameters is not abstraction disagreement** | two independent designs differ on flip angle, duration and bandwidth derivation — and agree on structure, sign convention and what must follow what. That is what a *caller-owned parameter* looks like, and escalating to `REFERENCES_DISAGREE` would have overstated the uncertainty. Promoted as a generic rule. |
| **reuse of waveform machinery ≠ reuse of the public abstraction** | `Excitation(use='saturation')` was a live alternative and was rejected: a semantic-purpose flag turns a module with one physical contract into a generic RF pulse whose *label* carries the meaning. |
| **claim scope and validation debt must be explicit** | the `evidence_state` block came from here — deferrals had been documented in five places and none of them carried a machine-readable trigger. |
| **scanner validation is claim-driven** | there is deliberately no rule that every GREEN module is eventually scanner-tested. A carrier offset read off the compiled file is *more* direct than an image confirming the same thing through a phantom. |

It also produced a finding that was not about the workflow at all: **pypulseq's `Opts` defaults
`B0` to 1.5 T in silence**, so a 2.89 T protocol missing `B0=2.89` puts fat at −220 Hz instead of
−424. Legal, plausible, wrong, invisible. Pinned by a test.

## 4. Spiral — `NEW_LEAF / GREEN`, after archaeology and debugging

The hardest case, and the one that pressed on separation. It arrived with a large historical PR
(#23) proposing a module, three shared helpers and two compiler changes, and the workflow had to
take it apart.

**What it separated, and how each part came out:**

| | |
|---|---|
| **m1 compiler defect** | **real.** Reproduced with no Spiral module — six of ten legal sign-changing gradients refused, unpredictably. Fixed in standalone **PR #30**, with the correction *claim-driven* (a floor on peak `\|g\|`) rather than calibrated to one family. |
| **m0 compiler claim** | **not reproduced.** The comment/code mismatch is real; the failure is not. Not ported. |
| **three promoted helpers** | **deferred.** One consumer each on current `main`, including after PR #27's EPI work. |
| **four-mode family** | **retained**, after independent physical re-derivation from the coupling test — not because PR23 had it. |
| **`writeSpiral.m`** | **changed evidence strength, not abstraction scope.** It witnesses spiral-out only, so `in`/`in-out`/`out-in` rest on analytic and metamorphic evidence. That is a claim-scope statement, not a reason to ship a smaller module. |
| **endpoint policy** | **deliberately differs from the reference**, which ends at full gradient and therefore has no spiral-in. Cost quantified at **0.20–0.78 %** of readout duration before adoption. |
| **a later m0-looking symptom** | **isolated to this module.** Splitting a bare arbitrary gradient across ADC-driven block boundaries is exact to 1e-13, so the compiler was proved correct and the defect was mine. |

That last row is the strongest evidence for the compiler-escalation rule in the whole project.
The same symptom that would have "confirmed" PR23's m0 claim appeared twice; the candidate-free
reproducer is what stopped it being claimed both times.

**Two generic implementation lessons** came out of the debugging and are now rules E and F:

> Measure a module on the lattice it actually emits, not on a more convenient representation.

> Validate shared temporal occupancy after events are jointly placed; individually correct events
> can become physically wrong when one extends the compiled block beyond another waveform's
> designed support.

The concrete instances — a rewinder sized against a design end-point, slew measured on
raster-spaced differences when the emitted edges are half a raster, an acquisition outlasting its
gradient by 10 µs — stay in the candidate record. Only the generic form is in the Skill.

## 5. Rules introduced because reality required them

None of these was invented in advance. Each exists because a run needed it.

| rule | forced by |
|---|---|
| waveform similarity is a discovery heuristic, not boundary evidence | STEAM |
| a licence can block a candidate before its physics does; four uses of a source distinguished | STEAM |
| `status: UNDECIDED` — stopped before deciding | STEAM |
| reference disagreement on protocol parameters ≠ abstraction disagreement | SaturationPrep |
| `evidence_state` with a required `revisit_trigger` per deferral | SaturationPrep |
| three kinds of human evidence kept apart; scanner work claim-driven | SaturationPrep |
| reference coverage sets evidence *tier*, not abstraction *boundary* | Spiral |
| internal lowering boundary ≠ public Module boundary | Spiral |
| measure on the emitted lattice (rule E) | Spiral |
| validate occupancy after joint placement (rule F) | Spiral |
| `APPROVED_FOR_IMPLEMENTATION` and `IMPLEMENTATION_DEFECT_UNRESOLVED` | Spiral |
| coverage classification replacing the word "duplicate" | MRzero coarse scan |

## 6. Rules that survived all three

| | |
|---|---|
| **the pipeline order**, and especially *no extraction before the acceptance claim is written* | a criterion written after the code describes the code |
| **the ownership test** — if a component cannot determine the value from its own contract and parameters, it should not own that event | separated every disputed case across five families without a tiebreaker |
| **claim scope stated in both halves** | every candidate; the half that gets dropped is always `does_not_establish` |
| **family-specific acceptance criteria** | four candidates, four criteria: event identity, trajectory agreement, lattice-and-semantics, analytic-plus-metamorphic |
| **one consumer does not justify promotion** | `RadialGRETR`, `seqcraft/trajectory/`, three Spiral helpers, a reconstruction package API — all correctly not created |
| **compiler changes need a candidate-free reproducer** | vindicated twice in one candidate |
| **`NO_NEW_MODULE` is a success** | STEAM stopped; nothing was forced |

## 7. Known limits left beyond v0

Deliberate, each with a trigger recorded in the relevant `evidence_state`:

```text
Spiral echoes > 1              an extension of the same leaf, not a gap in it
Radial Layer 3                 waiting on an independently reviewed reconstruction path
Spiral Layer 3                 same, and the trajectory claims are carried more directly anyway
STEAM                          evidence-blocked; re-run if a permissive reference appears
DREAM                          revisit after the STEAM boundary is understood
reconstruction promotion       a shared contract exists and is small; promoting it is a
                               separate decision, and example-only remains valid
the m0 comment/code mismatch   real, unproven as a defect, not ported
```

## 8. Why v0 is ready to close

**v0 success does not mean every prospective calibration produces a new GREEN Module.** The three
outcomes are useful *because they differ*:

```text
STEAM             the evidence says: do not decide yet
SaturationPrep    the evidence says: extract a simple preparation leaf, and do not overdesign it
Spiral            the evidence says: retain a sophisticated family, reject unsupported historical
                  changes, isolate the real compiler defect from the imagined one, and debug the
                  realization until the declared claims actually pass
```

A workflow that produced three GREENs would have told us only that it can produce GREENs. One that
stops, one that extracts something small, and one that takes a large historical artifact apart and
keeps only what survives independent review — between them those cover the regimes v0 was built
for, and each was reached without the answer being written down first.

The schema survived a fourth, fifth and sixth candidate, gaining four fields and three vocabulary
values, all traceable to a run that needed them. That is growth by reasoning class rather than by
sequence count, which is the rule the Skill sets for itself.

**What would change this conclusion:** a fourth prospective run reaching a regime none of these
covers — a `RED`, or a candidate whose physics is adequate but whose abstraction is not. Those
remain untested, and are named here rather than left to be discovered.
