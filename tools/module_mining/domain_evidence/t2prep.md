# Domain evidence — T2 preparation (driven equilibrium)

A curated card, so a later reader does not have to re-read the source. **Not a literature review.**

**Primary source.** Bernstein, King & Zhou, *Handbook of MRI Pulse Sequences*, Elsevier 2004 —
**Ch. 17 "Advanced Pulse Sequence Techniques", §17.4 "Driven Equilibrium"** (author K. F. King),
book **pp. 888–895**. §17.4.2 *Gradient-Echo-Driven Equilibrium Preparation* is **pp. 893–895**.
The PDF is not in the repository and must not be committed; it lives outside the tree.

The Handbook's own term is **driven equilibrium (DE)**, also DEFT and fast recovery (FR). It names
the gradient-echo-preparation use explicitly:

> "This use of driven equilibrium is also called **T2 preparation (T2-prep)**." — p. 888

---

## What the source establishes

| claim | page | what it establishes | what it does **not** |
|---|---|---|---|
| The canonical RF pattern is **nonselective** `90x–τ–180y–τ–90₋ₓ`, or `90x–τ–180x–τ–90₋ₓ` less commonly | 888 | the family's base structure, the π/2 phase offset of the refocusing pulse relative to the excitation, and that the tip pulses are **not** slice-selective | any particular number of refocusing pulses, or any pulse shape |
| `Mz·exp(−Tprep/T2)` is the longitudinal magnetisation after the second 90°, assuming perfect flip angles | 893 | **the semantic quantity is one preparation time `Tprep`, and it alone sets the T2 weighting** | the flip-angle-error case, which is the next row |
| "the longitudinal magnetisation prior to the first 90° **must be calculated by taking into account steady-state conditions — it is not simply the equilibrium longitudinal magnetisation**" | 893 | that the preparation's *input* magnetisation is a sequence-wide property the preparation cannot know. It is `Mz` in the formula, supplied, not computed | how to compute it — the source defers to Parrish & Hu (1994) |
| "Typically Tprep ≈ 50 ms" | 893 | the protocol scale | a limit |
| "Greater T2 weighting can be obtained by using **multiple 180° pulses to prolong Tprep**" | 893 | **the refocusing train is the canonical means of reaching a longer `Tprep`** — so the train serves the preparation time and is subordinate to it | which count is standard |
| "Sometimes **composite pulses** are also used for the 90° or 180° pulses to **decrease sensitivity to B0 and B1 inhomogeneity**" | 893 | composite pulses are a canonical, optional robustness *realisation*, not part of the definition | that composite pulses are required, or which composite |
| DE preparation RF is "followed by **gradient spoiling lobes on one or more axes** to dephase residual transverse magnetisation" | 893 | the spoiler is part of the preparation, it follows the tip-up, and the axis count is free | a specific moment or axis set |
| Degraded T2 contrast from T1 regrowth during the imaging period is reduced "by acquiring the center of k-space before substantial regrowth occurs" — **centric or elliptical-centric ordering** | 893 | **k-space ordering is downstream acquisition policy**, outside the preparation, and is the canonical mitigation for the preparation's main weakness | that any particular ordering is required |
| B1 inhomogeneity leaves longitudinal magnetisation "not directly determined by T2 contrast"; mitigated by a dephasing lobe **G1 after the first 90°** and a rephasing lobe **G2 after each RF pulse in the subsequent imaging sequence**, where "in practice G2 is combined with the slice-selection rephaser" | 893–894, Fig. 17.36 | **a canonical variant whose correctness spans the preparation/acquisition boundary** — the preparation emits G1 and the *readout* must emit G2. Cost: reduced SNR, slightly longer minimum TE and TR | that this variant is standard rather than optional |
| Gradient lobes placed **between** the 90° and 180° make it a **diffusion preparation** | 895, Fig. 17.37 | diffusion preparation is a variant of *this same family*, not a separate one | the diffusion physics, which is §9.1 |
| Phase-cycling the tip-up — `90x–τ–180y–τ–90₋ₓ` then `90x–τ–180y–τ–90y`, combined as root-sum-of-squares — recovers signal lost to incomplete rephasing | 895 | that a **two-acquisition** scheme exists and that combining the two is a *reconstruction* step | that one acquisition is insufficient in general |
| "DE preparation can also be used along with other preparation pulses such as **fat saturation** and **adiabatic inversion**" | 895 | composability with sibling preparations — evidence for a preparation leaf rather than an imaging kernel | ordering constraints |
| DE is "a pulse sequence **module** that can be appended to" a sequence | 576 (§11.x cross-ref) | the source itself treats it as a modular unit | |

## Canonical uses named by the source

T2-weighted head imaging with gradient echo (Mugler 1991); cardiac imaging with susceptibility
contrast agents (Sakuma 1994); **blood–myocardial contrast in cardiac gradient echo** (Brittain
1995; Wielopolski 1995) — the last being the coronary-angiography application that the concrete
implementations cite.

## What the source leaves open

- the refocusing **count** — "multiple 180° pulses", no number;
- the **bookkeeping** of `Tprep` against real pulse durations. The physics is stated for ideal
  instantaneous pulses, where "the time between the DE 90° pulses" is unambiguous;
- any spoiler moment;
- whether the tip-up is a simple `−90` or a composite.

Those are exactly the questions the implementation witnesses answer differently, and therefore the
ones a mode contract has to fix — see `plans/t2prep/candidate.yaml`.

## Provenance note

Facts above were located by searching the PDF's text for "T2 preparation" and "driven
equilibrium", which pointed at the index entry (book p. 890 index → pp. 888–890) and then at
§17.4. Only pp. 888–889 and 893–895 were read. Nothing else in the Handbook was consulted for this
card.
