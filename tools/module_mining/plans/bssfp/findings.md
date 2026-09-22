# balanced SSFP — Phase C second fine scan

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-22, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way.

**Date:** 2026-09-22 · **Skill:** v0, frozen and unmodified
**Result:** `NEW_KERNEL / APPROVED_FOR_IMPLEMENTATION` + `ARCHITECTURE_REVISIT_CANDIDATE`,
for the **RF-to-RF repetition** inside a layered answer — with one quantified cost for a reviewer
to decide.

Nothing implemented. Record: [`candidate.yaml`](candidate.yaml). Domain card:
[`../../domain_evidence/bssfp.md`](../../domain_evidence/bssfp.md).

## 1. Domain sources

*Handbook of MRI Pulse Sequences*, **Ch. 14 §14.1 "Gradient Echo"**, subsection **"Balanced SSFP
(True FISP)"**, book **pp. 592–597**; **Figures 14.12 and 14.13** (pp. 600–601); **Table 14.2**
(p. 597); vendor names in Table 14.1 (p. 584). Only those pages were read.

## 2. Witnesses, and this is the best evidence base in the project

| | licence | dim | TR | TE | catalyzation |
|---|---|---|---|---|---|
| `writeTrufi.m` (Pulseq) | MIT | 2D | **output** — "you just specify the ADC time" | `TR/2` | **α/2 + TR/2** |
| OpenMRF `BSSFP` | MIT | **3D** | input, fill split in halves | `TR/2` by construction | **flip-angle ramp** |
| `mrseq` `t1_molli_bssfp` | Apache-2.0 | 2D, gated | **input with a floor**, `None` → minimum | **independent, with a floor** | flip-angle ramp, counted |
| MRzero bSSFP notebook | AGPL+EULA | — | — | — | **oracle only**, not used this pass |

Three independent groups, two languages, 2D and 3D. They **agree on every physical condition** and
disagree only on parameterisation.

## 3. Canonical contract

**Two conditions, and they are not the same one.** "To establish SSFP, the gradient area on any
axis **must not vary among** the TR intervals" (p. 592). Balanced imposes "a **further
condition** … that the gradient area on any axis **is zero during each** TR interval" (p. 593).
Strictly stronger. **No implementation states this distinction**, because each implements one side
of it.

**The interval is RF-to-RF.** Figure 14.12's caption states the balance over "one TR interval" and
then draws the *next* repetition's RF and slice gradient "to emphasize the symmetry of the
waveform" — the domain source has to reach into the next repetition to show it. `writeTrufi`'s own
comment independently names it: the ringdown correction is "needed to center the ADC and the
gradient echo in the center of **RF-RF period**".

**TE = TR/2 is a contrast condition, not an invariant.** T2 rather than T2\* appears in the
exponent "**if** the balanced SSFP signal is rephased in the center of the TR interval"; away from
centre, T2′ weighting enters (p. 593), and Figure 14.13 gives partial-echo bSSFP at TE < TR/2 as a
named realisation. Two of three witnesses hard-wire TR/2.

**Phase progression is a frame convention.** Sign alternation "is **equivalent to** a pulse
sequence … that has constant precession … by φ = 180° in each TR interval" (p. 593). That is also
what a band *is* — a region already precessing 180°/TR loses the benefit. And the Handbook notes
alternation "is analogous to the driven equilibrium methods discussed in §17.4" — a direct link
back to C1.

## 4. Is RF-to-RF the natural boundary? Yes, and it costs something

It is the only one of the four hypotheses that can state a correctness condition **it can also
check on a single instance**: net area zero per axis over its own duration, echo at its declared
TE, no reference to neighbours.

But a module's events cannot continue into the next instance of itself, so the boundary must fall
where every axis is at zero. **Measured:**

```text
readout axis, the lobe pair that straddles the boundary
  rewinder of TR n, alone           380 us
  prephaser of TR n+1, alone        380 us
  two separate trapezoids           760 us
  ONE trapezoid, same moment        540 us
  cost of keeping them separate     220 us per repetition per axis   (~5.5 % of a 4 ms TR)
```

The slice axis carries the same split between the excitation rephaser and the next selection lobe.
`writeTrufi` reshuffles exactly this away — splitting the slice gradient at the RF and recombining
it with the *previous* repetition's rephaser — and it can only do that because it is a script, not
a composable unit.

**SeqCraft is better placed than Pulseq here, but not enough.** I measured that a gradient spanning
both an RF and an ADC stays continuous at constant amplitude across both seams inside one module —
two blocks, one waveform, never forced to zero. So the compiler already does *inside* a module what
`writeTrufi` does by hand. What it cannot do is carry a waveform across the repetition boundary.

## 5. Catalyzation ownership — outside the kernel

The Handbook settles it: the approach to steady state takes "four to five times T1", making
catalyzation "an important practical problem", and it names **at least three** schemes — α/2 with
TR/2, a linear flip-angle ramp (Le Roux 2003; Hargreaves 2001), and ramp-down storage on the
longitudinal axis for an interposed preparation. **An abstraction with three named solutions and
no criterion for choosing is policy.**

The witnesses agree structurally: all three put catalyzation outside the repetition, and `mrseq`
counts its startup pulses in units of the *same* TR — so a catalyzing repetition is the ordinary
repetition at a different flip angle. `GRE2DTR` sets the shipped precedent, with dummies in the
example's loop.

## 6. 2D versus 3D — segmentation is not repetition physics

OpenMRF is 3D and its repetition is structurally identical to `writeTrufi`'s 2D one; `loop_y` and
`loop_z` select encode values and nothing else. `mrseq` is 2D and cardiac-gated, and what its
segmentation changes is how many repetitions share one catalyzation and how often the steady state
is interrupted. §14.1 is dimension-agnostic.

**So a segment boundary carries steady-state physics, not bSSFP physics** — it is where
catalyzation is needed again, or where magnetisation must be stored and restored (p. 595). One
kernel serves 2D and 3D.

## 7. The four hypotheses

| | can it state a condition it can check? | needs from outside | verdict |
|---|---|---|---|
| **A** RF-to-RF kernel | **yes, on one instance** | repetition index, encodes, flip angle | **the kernel** |
| **B** balanced requirement | stateable, not checkable alone | every other event on the axis, budget, limits | **not expressible today** |
| **C** steady-state segment | **no** — its claim is about magnetisation | T1, T2, a simulation | policy, not an abstraction |
| **D** layered A + composition | yes, as A | as A | **selected** |

## 8. Recommendation

**Implement A as the kernel, in the layered form D.** Expose TE with `TR/2` as the default and the
T2′ consequence documented — the faithful choice, per the domain source and `mrseq`.

**The one thing for a reviewer to decide first:** accept 220 µs per axis per repetition, or wait
for the capability that removes it. Accepting ships a correct composable kernel a few per cent
slower than a hand-written script. Declining waits on C3.

**That convergence is this scan's most useful result.** bSSFP and velocity encoding reach the same
missing capability — a requirement realised jointly with its surroundings — from opposite
directions, and neither was looking for it.

## 9. v0

**Handled well.** The ownership test did the whole boundary analysis: "can this state a correctness
condition it can check" separated all four hypotheses without a tiebreaker. `NEW_KERNEL` and the
kernel-layer definition fit exactly. `ARCHITECTURE_REVISIT_CANDIDATE` carried the cost finding.
Rule E shaped the acceptance claim. The Layer 1/2/3 split is what forced the sharpest sentence in
the record — a balanced legal repetition is Layer 1 and Layer 2 green and says nothing about
steady state, which only Layer 3 reaches.

**Reasoned outside v0.** Two things. There is still no `evidence_class`, recorded by hand again.
And v0's lights have no way to say *"the boundary is settled and a quantified cost needs a product
decision"* — `APPROVED_FOR_IMPLEMENTATION` is right about the design and silent about the trade;
`YELLOW – PHYSICAL_BOUNDARY_UNCLEAR` would be false. The trade sits in `open_questions` and the
decision summary, which is adequate but is not in the light.

## 10. Did domain evidence change the outcome again? Yes — and differently

In C1 it established a family one code witness could not. **Here the corpus was rich — three
independent implementations — and would still have taught the wrong contract.** Two of three
hard-wire TE = TR/2, so a boundary derived from code alone would have made it an invariant; the
Handbook says it is a contrast condition with a named realisation on the other side. And the
SSFP-versus-balanced distinction appears in no implementation, because each implements one side.

So domain evidence is not a fallback for thin corpora. In C1 it supplied what was missing; here it
**corrected what three witnesses agreed on.** Only the second was a surprise.

Two for two. Holding the formalisation question until C3, as instructed.
