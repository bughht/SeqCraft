# Domain evidence — balanced SSFP (True FISP / FIESTA / balanced FFE)

Curated card. **Not a literature review.**

**Primary source.** Bernstein, King & Zhou, *Handbook of MRI Pulse Sequences*, Elsevier 2004 —
**Ch. 14 "Basic Pulse Sequences", §14.1 "Gradient Echo"**, subsection **"Balanced SSFP (True
FISP)"**, book **pp. 592–597**, with Figures **14.12** and **14.13** (pp. 600–601) and **Table
14.2** (p. 597). Vendor-name cross-reference in Table 14.1 (p. 584). The PDF is not in the
repository.

Only pp. 592–597 and the two figure captions were read.

---

## The two conditions, and they are not the same one

| claim | page | establishes | does **not** establish |
|---|---|---|---|
| "To establish SSFP, **the gradient area on any axis must not vary among the TR intervals**." | 592 | the **SSFP** condition — constancy across repetitions, not zero | that the area must be zero |
| "If a **further condition** is imposed that **the gradient area on any axis is zero during each TR interval**, a very different steady state results." | 593 | **this is what "balanced" means** — zero net area per axis per TR, and it is a *strictly stronger* condition than SSFP | the interval's endpoints, which Fig. 14.12 supplies |
| Fig. 14.12 caption: "On any of the three gradient axes, the net gradient area is zero during one TR interval… **The RF pulse and part of the slice-selection gradient for the next TR interval are shown (dashed lines), to emphasize the symmetry of the waveform.**" | 600 | that the natural interval is **RF-to-RF**, and that showing the symmetry *requires reaching into the next repetition* | that a repetition is an independent object |

The distinction in the first two rows is load-bearing: a sequence can satisfy SSFP and not be
balanced. Spoiled GRE, SSFP-FID and SSFP-echo differ from balanced SSFP **only** in the gradient
waveform class (Table 14.2, p. 597).

## Consequences the source states

| claim | page | establishes | does not establish |
|---|---|---|---|
| The SSFP-FID and SSFP-echo peaks "coalesce; that is, they rephase at the same time TE. Therefore the balanced SSFP signal is the **coherent sum** of the two" | 593 | why balance changes the contrast regime rather than just the gradients | |
| Signal **with** sign alternation (14.23) exceeds that without (14.24), "so in practice sign alternation is used" | 593 | RF phase alternation is the canonical progression, on signal grounds | that it is required |
| "the effect of the sign alternation is **analogous to the driven equilibrium methods** discussed in §17.4" | 593 | a direct physical link to the T2Prep family (C1) | |
| "The effect of sign alternation is **equivalent to** a pulse sequence (without sign alternation) that has constant precession … by φ = 180° in each TR interval" (Hinshaw 1976) | 593 | that the phase progression is a **frame convention**, interchangeable with an off-resonance condition — so "alternation" and "180° per TR" are one thing | a preferred implementation |
| T2 rather than T2\* appears in the exponent "**if** the balanced SSFP signal is rephased in the center of the TR interval (i.e. **TE = TR/2**) … As the peak … is moved away from the center, T2′ weighting is introduced" | 593 | **TE = TR/2 is the condition for the stated T2/T1 contrast, not a legality constraint** | that TE must equal TR/2 |
| Fig. 14.13, partial-echo balanced SSFP: "The net area on the readout gradient waveform remains zero. **Because TE < TR/2**, a slight amount of susceptibility weighting is introduced." | 601 | that **TE ≠ TR/2 is a named, legitimate realisation** whose consequence is contrast, not illegality | |
| Banding: a region where φ ≈ 180° in a sign-alternated sequence, which "effectively removes the sign alternation"; worsens as B0 rises | 593–594 | off-resonance sensitivity as intrinsic to the method | a correction |

## Steady state, catalyzation and interruption — named as a separate problem

| claim | page | establishes | does not establish |
|---|---|---|---|
| "The approach to steady state … can take **on the order of four to five times T1**. This is a long time to wait relative to TR, so accelerating or **catalyzing** the approach to steady state is an important practical problem." | 595 | **a single repetition is waveform-correct and physically meaningless**; the steady-state claim is a train-level property | how many repetitions suffice |
| "A simple and relatively effective catalyzing method … is to apply a **half-flip angle pulse with a half-TR interval**, that is, (θx/2) − (TR/2) − θ₋ₓ − TR − θx …" | 595 | the α/2 scheme, canonically | that it is the only one |
| "More sophisticated catalyzing methods such as **linearly ramping up flip angle** are described in Le Roux (2003) and Hargreaves et al. (2001)." | 595 | **that both schemes the code witnesses use are canonical named alternatives** | which is preferred |
| "once the steady state is established it **often has to be interrupted**, for example for … chemical saturation pulses or to resynchronize … with a detected cardiac trigger"; the SSFP "can be temporarily **stored on the longitudinal axis**" by ramping down, then restored | 595 | interruption and store/restore as first-class, and a named RF recipe for it | |
| CISS / multi-acquisition: two acquisitions, one sign-alternated and one not, separately reconstructed and combined, to shift the bands | 595–596 | that multiple phase-cycle acquisitions are an **acquisition-and-reconstruction** policy above the repetition | |

## What the source does not determine

- Whether TR is an input or an output. It treats TR as a protocol parameter throughout and never
  discusses its minimum.
- Where a *software* repetition boundary should fall. It states the balance over "one TR interval"
  and draws the neighbouring RF to show the waveform is symmetric across that boundary.
- 2D against 3D, or segmentation. §14.1 is dimension-agnostic; nothing in the balanced-SSFP
  subsection makes segmentation a property of the repetition.
- How many catalyzation pulses, or how to verify that steady state was reached.

## Provenance note

Located by searching the PDF text for "balanced SSFP", "TrueFISP" and "FIESTA", which pointed at
§14.1; then the two figure captions and Table 14.2 by figure number. Nothing else in the Handbook
was consulted for this card.
