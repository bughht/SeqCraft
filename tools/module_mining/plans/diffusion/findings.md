# Diffusion — architecture stress test and proof by implementation

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-22, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way.

**Date:** 2026-09-22 · **Skill:** v1
**Result:** `DiffusionSEPrep` implemented, a complete acquisition showcased, both Phase-E
hypotheses stress-tested and **both survived** — and independent validation caught four defects,
two of them in the validator itself.

```text
1. targeted fine scan          section 1-2
2. architecture decision       section 3
3. production implementation   section 4
4. complete acquisition        section 5
5. Layer 1 / 2 / 3 validation  section 6
```

This candidate was chosen to try to **break** two claims Phase E made. Diffusion is the hardest
case the batch scan surfaced for both of them, and the honest report is that neither broke. What
did break was something else, and section 6 is about that.

---

## 1. What the scan asked

Phase E left two hypotheses standing, and both were reached by reasoning about C1, C2 and C3
rather than by building anything:

> **H1 — realization is analytic first, with a numerical optimizer as a fallback.** The canonical
> physical constructions have closed forms; an optimizer buys generality, not the canonical
> answer. It may be a realization strategy and may not define public physical semantics.

> **H2 — joint realization can happen before LogicBlock materialisation.** A kernel that owns a
> waveform spanning several leaves can complete its design from construction-time properties, so
> nothing needs to inspect or rewrite emitted events.

Diffusion is the adversarial case for both:

| | why it should break the hypothesis |
|---|---|
| **H1** | the copyleft corpus contains an explicit *minimum-TE optimizer* for exactly this problem (`opt_TE_bv_SE`), and a 4D-flow pipeline hands the same class of problem to a constrained solver. If any canonical case needs an optimizer, this is it |
| **H2** | the encoding is defined *around a pulse another module owns*. The lobe width depends on the refocusing pulse's duration, its position depends on that pulse's centre, and the echo time depends on both. If anything needs to look at emitted events, this is it |

The scan was targeted rather than broad: the goal was a decision and then an implementation, not
a survey. Sources, licences and what each establishes are in
[`candidate.yaml`](candidate.yaml); the literature card is
[`../../domain_evidence/diffusion.md`](../../domain_evidence/diffusion.md).

---

## 2. What the evidence situation actually is

The asymmetry here is the most instructive thing the scan found, and it is not about diffusion.

```text
permissive witness (MIT)       solves the EASY direction   with INCOMPLETE physics
copyleft witness (AGPL)        solves the HARD direction   with complete physics
authoritative literature       gives the complete physics  and no algorithm
```

`writeEpiDiffusionRS.m` takes TE as an input and inverts the analytic b-factor for the amplitude —
the right shape — but its own helper is documented "for rect gradients" with "for trapezoid
gradients: **TODO**", so the shipped physics ignores the ramps. `PulseqDiffusion` has both the
exact b-value and a minimum-TE search, and is discovery-only under AGPL.

**Only the domain reference closes the gap.** The Handbook's Figure 9.3 gives the ramp-inclusive
closed form that the permissive implementation marks TODO, which is what allowed this module to be
built from a permissive basis at all. That is the third time in Phase C that domain literature
changed an outcome, and by a third mechanism:

| | what the literature did |
|---|---|
| **C1 T2Prep** | supplied a family the corpus did not contain |
| **C2 bSSFP** | corrected a contract three witnesses agreed on |
| **C3 flow** | corrected an architectural inference drawn from one implementation |
| **Diffusion** | **supplied physics the only permissive implementation is missing** |

> **The counting rule stands and is worth restating.** The number of executable repositories is
> not a vote on whether a physical abstraction exists — and here it is not a vote on whether the
> physics is *right*, either. One MIT witness and a textbook beat one MIT witness.

---

## 3. The ownership decision

Applying the ownership test —

> If a component cannot determine the correct value using only its own physical contract and
> parameters, it should not own that event.

— to each candidate owner:

| owner | can it determine the lobe width? | |
|---|---|---|
| `Excitation` | no | it knows nothing about a refocusing pulse |
| `Refocusing` | no | it knows nothing about a b-value |
| a new `DiffusionLobe` leaf | **no** | Δ is measured between lobe centres *across* the 180, so the width depends on that pulse's duration |
| the acquisition | yes, but | it would have to re-derive the physics, and every acquisition would repeat it |
| **a kernel over all three** | **yes** | and it is the smallest thing that can |

So: a **kernel**, `modules/kernel/diffusion_se.py`, owning the excitation, the refocusing pulse
and the two lobes. This is the Phase-E "joint-realization kernel" in its first production form.

```text
Excitation, Refocusing   own   their own waveforms, and report their own timing
DiffusionSEPrep          owns  the lobes, the echo time, and where the 180 goes
the acquisition          owns  the echo time across b -- which no module can see
```

**What was NOT introduced**, deliberately: no `RequirementIR`, no `ConstraintSolver`, no
`GradientDesignSpec`, no design context of any kind. No change to `LogicBlock`, to the compiler,
or to what a `GradientEvent` knows. The kernel is an ordinary `sc.Module`.

---

## 4. H2, tested: no new interface was needed

The kernel needs six facts from the two RF leaves:

```text
Excitation.time_to_center()         Refocusing.time_to_center()
Excitation.time_to_rephaser()       Refocusing.time_to_crusher()
Excitation.rephaser_duration_s      Refocusing.crush_duration_s
```

**Every one of them was already a construction-time property.** The kernel constructs its leaves,
reads those six numbers, completes the coupled design, and only then calls `build()`. Nothing is
inspected after materialisation, nothing is rewritten, and no leaf was modified.

> **H2 survives, and more cheaply than Phase E expected.** The prediction was that a joint owner
> would need *some* interface for semantic facts. For this candidate the interface already
> existed, because the leaves were already written to report their own semantic instants — which
> is what `se_2d/`, `fse_2d/` and `se_spiral_2d/` needed them for.

This is one candidate and not a proof. What it does establish is that "the joint owner must see
the events" was not true here, and that the burden is on the next candidate to show otherwise
rather than on this architecture to pre-build a mechanism.

---

## 5. H1, tested: the canonical case is entirely closed form

With the amplitude at the cap and Δ = C + δ + ε (C the refocusing block, ε the ramp), the
Handbook's Figure 9.3 expression becomes a **cubic in δ**:

```text
(2/3) d^3 + (C + e) d^2 - (e^2/6) d + (e^3/30) = b / (gamma^2 G^2)
```

whose smallest positive real root is the lobe width. This is the Handbook's own Example 9.2
calculation, which states the cubic has "a single physically meaningful solution".

The echo time does not enter the lobe design at all. It only decides whether the pair *fits*, and
both windows — excitation-to-first-lobe and second-lobe-to-echo — grow as TE/2, so the shortest
echo time is a **maximum of two linear expressions**. Rastering the width up and trimming the
amplitude back needs one square root, because b is quadratic in amplitude.

**No optimizer, no search, no iteration.** Measured: the delivered b matches the request to
0.003 % from 100 to 3000 s/mm², and a third implementation that integrates |k(t)|² directly agrees
with the closed form to 1e-9.

### Where the line falls, stated rather than implied

H1 is confirmed *for this mode* and the correction from C3 must not over-shoot here either. What
is **not** established:

- **Twice-refocused / eddy-current-compensated designs.** Four lobes, two refocusing pulses, and
  cancellation conditions the single-refocused case does not have. The Handbook names them,
  motivates them and defers the design. This is the best test of whether an optimizer is ever
  needed here, and it should be scanned in its own right.
- **A genuine minimum-TE search over the readout as well as the encoding.** What this module
  returns is the shortest echo time for *this* geometry — full amplitude, symmetric about the
  refocusing pulse. Another arrangement might be shorter and no optimality claim is made.
- **Multi-constraint problems** — several moment orders, fixed waveform segments, non-zero
  endpoints, eddy-current or PNS constraints combined. A numerical backend stays a legitimate
  realization strategy there.

> `salvage/bvalue.py`, the b-value solve from the pre-reform module set, says "a closed form does
> not exist because `Delta` depends on `delta`" and steps up the raster to find the width. That is
> true only while the two are left independent; substituting Δ = δ + gap is what makes it a cubic.
> **The earlier code was not wrong about the physics — it stopped one substitution short.**

---

## 6. Validation, and the four defects it caught

The validation ladder was run in full. The interesting part is not that it passed.

| layer | what | outcome |
|---|---|---|
| **1** | 31 tests: the Handbook closed forms, delivered b measured independently, metamorphic scaling, timing against `sc.kspace`, identical lobes, the lobe-before-echo invariant, the refusals, and agreement with a third implementation | GREEN |
| **2** | `examples/dwi_se_epi_2d/01_build.ipynb` — a complete diffusion-weighted spin-echo EPI, compiled, with the readout's centre-line echo on the spin echo | GREEN |
| **3** | `examples/dwi_se_epi_2d/02_simulate_and_reconstruct.ipynb` — MRzeroCore on the emitted `.seq` | GREEN |

### 6.1 Two defects in the module, caught by Layer 1 and Layer 2

**The lobe-centre separation was one ramp short.** Every lobe was legal, the geometry looked
right, the module's own number said 1000 — and the sequence delivered 1042, because the amplitude
solve had quietly compensated for the wrong Δ. Only an integration of the *emitted* gradients
could see it.

**The lobe was sized without its own ramp**, so the second lobe ran one ramp past the spin echo
and left the readout nowhere to go. This one appeared **only when a real acquisition was composed
against it** — the echo-time fit loop would not converge. That is Rule G's argument in one
sentence: the composition-level claim has to be validated on the composition.

### 6.2 Two defects in the *validator*, caught by Layer 3

This is the result worth carrying forward.

`sc.b_value` is the independent analyser Rule G asks for: it integrates the emitted gradients,
conjugates at every refocusing pulse, shares no code with the module, and measures the whole tree.
It caught 6.1's first defect. Layer 3 then caught **it**.

Simulating an unweighted acquisition and reading off the attenuation gives the b-value the
magnetisation actually experienced. The two disagreed:

```text
sc.b_value said   0.639 s/mm^2  (and 0.479 for the same encoding at a shorter TE)
the simulator     0.074 s/mm^2  (and 0.074)
```

The *second* symptom is the diagnostic. b is an integral of k(t)², and k is zero while nothing is
playing, so a b that grows with **dead time** means the integral believes k is some nonzero
constant during the gaps. There were two independent reasons:

| defect | consequence |
|---|---|
| the phase origin was the start of the tree, not the excitation's centre | a slice-select rephaser undoes the lobe area *after* the RF centre, not half the lobe — so half a lobe, 510 1/m, was left as a constant offset |
| refocusing pulses were located by `pp.calc_rf_center` **without the event's own `delay`** | the conjugation was applied 700 µs early, leaving a second constant of 240 1/m after the 180 |

Both are fixed and pinned by three regression tests in `tests/analysis/`. After the fix the
analyser reports 0.0749 against the simulator's 0.0741 — 1 %, which is raster integration against
state propagation — and no longer grows with echo time.

**This record is where that history lives.** The notebook shows the corrected comparison and the
current agreement; it is a tutorial, and replaying a fixed software defect is not what it is for.
See `examples/README.md` → *Writing example notebooks*.

> **Neither defect was visible on the diffusion axis.** There, k is zero before the encoding and
> flat across the whole refocusing block, so conjugating 700 µs early gives *exactly* the same
> answer, and there is no excitation gradient to leave an offset. Both were visible only on the
> slice axis — which the module does not own, does not design, and which no test in the
> repository was pointed at.

### 6.3 The rule this suggests

Rule G says the validator integrates rather than asks. That held. What this candidate adds is
narrower and sharper:

> **An independent validator is only independent of what it was pointed at.** `sc.b_value` was
> genuinely independent of `DiffusionSEPrep` and caught a real error in it. It was not independent
> of the *assumption both of them shared* — that the axis under test is the axis that matters —
> and nothing inside the repository could see that, because everything inside the repository was
> checking the same axis.
>
> The escape was a validator from a different **layer**: a Bloch simulator, which does not know
> what a b-value is and computes the attenuation from the magnetisation. That is what Layer 3 is
> for, and this is the first time in this project that it corrected a Layer 1 instrument rather
> than confirming one.

This is not an argument for more layers everywhere. It is an argument for **Layer 3 on any
candidate whose Layer 1 instrument is itself new**, which is a much smaller set. `sc.b_value` was
written for this candidate; that is exactly the condition.

---

## 7. What the acquisition owns, which no module can see

`01` §5 measures a fact that belongs to neither the kernel nor the readout:

| b | encoding alone | with the readout | set by |
|---|---|---|---|
| 0 | 9.00 ms | 51.46 ms | readout |
| 500 | 29.04 ms | 71.50 ms | readout |
| 1000 | 34.94 ms | 77.40 ms | readout |

Half the EPI train comes before its own echo, and that half has to fit between the second
diffusion lobe and the spin echo. The composition takes the larger of the two requirements — and
then a **third** constraint appears that is neither module's either: every b must be acquired at
*one* echo time, or the attenuation ratio carries T2 as well as D. A mismatch adds
(TE_b − TE_0)/(T2·b) to every reported ADC — 0.43 µm²/ms at a 60 ms T2 across the echo times this
protocol produces, measured against the closed form to four decimal places and **additive**, so it
does not calibrate away.

> This is the Phase-E point about acquisition-level ownership, in numbers. No module can decide
> it, no module should try, and writing it down is the deliverable.

---

## 8. What was deliberately not done

- **No twice-refocused mode.** A different physical contract, not a parameter away.
- **No direction schemes, b-matrices or tensor encoding.** The module encodes along the axes it
  is given; DTI is its own candidate.
- **No optimizer, and no optimizer hook.** GrOpt remains a possible realization backend for a
  future candidate that needs one and is not a SeqCraft dependency; nothing here creates a seat
  for it, because nothing here needed one.
- **No change to the top-level architecture.** `Module → build → LogicBlock → Tree → Compiler →
  PyPulseq Sequence` is untouched, as is what each layer knows.

## 9. For the approved queue

The queue after Phase E was T2Prep, then bSSFP, then the flow/moment pair. This candidate changes
the *evidence* behind how they should be approached, not their order:

1. **T2Prep** (NEW_LEAF, approved) — a leaf, and H2 is not in question for it.
2. **bSSFP** (C2) and **velocity encoding / flow compensation** (C3) converge on one missing
   capability: a joint owner for a waveform currently split across leaves. **That capability now
   exists and has shipped once**, as a kernel that reads construction-time properties. The next
   one should be attempted the same way before any new mechanism is proposed.
3. The open question C3 recorded — whether a merged design needs a solver — is now narrower.
   Diffusion's canonical case did not, and the reason was a substitution, not an algorithm. The
   C2/C3 merged forms should be attacked analytically first and the result reported either way.
