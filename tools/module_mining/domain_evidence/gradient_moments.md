# Domain evidence — gradient moments: flow encoding and moment nulling

Curated card for the C3 fine scan, covering **both** families because the question was whether
they are one. **Not a literature review.**

**Primary source.** Bernstein, King & Zhou, *Handbook of MRI Pulse Sequences*, Elsevier 2004.
**Source ID:** `handbook_mri_pulse_sequences_2004` — see [`README.md`](README.md) for resolution.

Two sections, and **the source puts them in different chapters**, which is the first finding:

```text
Ch. 9  Motion-Sensitizing Gradients
  9.2  Flow-Encoding Gradients          (M. A. Bernstein)   book pp. 281-291
Ch. 10 Correction Gradients
  10.4 Gradient Moment Nulling          (M. A. Bernstein)   book pp. 331-348
```

Same author, different chapters: one **sensitizes** (puts information in), one **corrects**
(takes an artefact out). Read: §9.2 pp. 281–283 and 287–290; §10.4 pp. 331–333 and 343.

---

## The shared object

| claim | section / page | establishes | does not establish |
|---|---|---|---|
| "With FE gradients, **each of the three logical axes** (readout, phase encoding, and slice selection) **is treated independently.** Therefore FE gradients can be added to a single axis, to any two, or to all three." | 9.2, p. 281 | per-axis independence — there is no vector coupling in the requirement | anything about how axes interact in hardware limits |
| "With GMN, **each of the three logical axes … is treated independently.** Therefore, GMN can be performed on a single axis, any two, or all three axes at once." | 10.4, p. 331 | **the identical sentence for the other family** — strong evidence they share a model | that the two are one abstraction |
| When m0 = 0, "the first moments would be **properly nulled at the end of the waveform regardless of where the time origin was selected** because the net gradient area … is zero" (Simonetti et al. 1991) | 10.4, p. 343 | that m1 is **origin-independent once m0 = 0** — so a moment requirement need not fix a time origin when the zeroth moment is already nulled | origin-independence when m0 ≠ 0 |

## Flow encoding — §9.2

| claim | page | establishes | does not establish |
|---|---|---|---|
| A flow-encoding gradient "**encodes** information about macroscopic flow or coherent motion **into the phase** of the MR signal" | 281 | the purpose: put information in | |
| "The most common FE gradient is a **bipolar** velocity-encoding gradient, which comprises **two lobes of equal area and opposite polarity**. Because the **net area … is 0**, it produces no net phase for stationary spins … For spins moving along the gradient, [it] produces a phase accumulation **linearly proportional to their velocity**." | 281 | the canonical realisation, and that the requirement is a **pair**: m0 = 0 and m1 at a target | that a bipolar pair is the only realisation |
| `VENC = π / (γ |Δm1|)` (Eq. 9.14) — "when the velocity component along the gradient direction is equal to ±VENC, the resulting phase difference is ±π" | 287 | **the public physical quantity is Δm1, the *change* in first moment on toggling** — not m1 itself | a preferred Δm1 for a given application |
| For a toggled bipolar, `Δm1 = m1 − (−m1) = 2AΔT` (Eq. 9.16), with a worked example giving VENC = 32.6 cm/s | 287–288 | the closed-form link from lobe area and separation to VENC | |
| "VENC is usually selected to be slightly higher (e.g. 30 %) than the expected peak vessel velocity" | 287 | that VENC is a protocol choice driven by the subject | |

## Moment nulling — §10.4

| claim | page | establishes | does not establish |
|---|---|---|---|
| "GMN … is the process of **modifying a gradient waveform** in order to make a pulse sequence more immune to image artifacts arising from motion. The modifications can include **changing the amplitude, duration, shape, or number of gradient lobes** that make up the waveform." | 331 | **the purpose: take an artefact out — and it does so by reshaping an existing waveform, not by inserting a block** | that a separate inserted block cannot also work |
| Binomial patterns: `1,−1` nulls m0 but leaves m1 and m2 non-zero; adding a third lobe in the `1,−2,1` pattern nulls m1 as well | 333–335 | a canonical closed-form construction family, and that **order costs lobes** | an optimal order |
| "a **velocity-compensated frequency-encoding waveform can be constructed from a minimum of three gradient lobes**, whereas acceleration-compensating that same waveform requires a minimum of four" | 333 | a minimum lobe count per moment order | |
| "GMN beyond a certain order will invariably become **counterproductive**. The optimal order … depends on the specifics of the application and on the performance of the gradient hardware." | 333 | that the order is a protocol-and-hardware choice, not a physical constant | |
| **The velocity-compensated slice-selection construction.** `ms` is "the absolute value of the **first moment of the slice-select lobe from the RF isodelay point** to the end of the slice-selection gradient lobe"; `ms` is given in closed form from the lobe's amplitude, rise time and plateau (Eq. 10.70); trapezoid width relates to area through the hardware limits (Eq. 10.71); substituting "yields **a quadratic equation for the unknown lobe area**" (10.72) whose "**physically significant root**" is given explicitly (10.73). | 343 | **that the canonical joint solve — target moment, moment already accumulated by a neighbouring lobe measured from a semantic instant, and hardware limits — has a CLOSED-FORM solution.** No numerical optimiser is required for this case | that arbitrary orders, arbitrary waveforms, or added constraints such as eddy currents are closed-form |

## How the two families compose — §9.2.3, p. 290

> "When the bipolar velocity-encoding gradient is **appended** to a velocity-compensated
> slice-selection or frequency-encoding waveform, the resulting **five lobes can be combined into
> three**. This combination **always reduces the minimum TE** and often reduces higher-gradient
> moments as well."

So the two coexist on one axis, and the efficient form **merges** them. The source names the
drawback: the concomitant-field phase error "is not necessarily equal for the two toggled
settings" once merged, so it no longer subtracts out in the phase-difference reconstruction —
though it "can be calculated exactly and then removed".

**Three named design choices** for the merged waveform, distinguished by where the two toggled
first moments sit — design details in Bernstein et al. (1992):

```text
two-sided      m1 straddles zero symmetrically for the two toggles
one-sided      one toggle is velocity-compensated (m1 = 0), the other carries the full dm1
minimum-TE     the most compact waveform, irrespective of where the individual m1 lie
```

"the choice among them is often determined by **how the accompanying magnitude images are
reconstructed**" — a reconstruction-driven choice, above the waveform layer.

## What the source does not determine

- Any general algorithm. It gives closed forms for the canonical per-axis cases and cites
  Bernstein et al. (1992) for the merged designs.
- Eddy-current or concomitant-field constraints as part of the moment solve; §10.1 and §10.3 treat
  those separately.
- How to choose among the three merged designs beyond "how the magnitude images are reconstructed".
- Anything about a requirement *interface* — moments are computed and waveforms are constructed in
  the same breath throughout.

## Provenance note

Located by a full-text search for "flow-encoding gradients", "moment nulling" and "VENC", then by
the chapter table of contents. Only the pages listed at the top were read.
