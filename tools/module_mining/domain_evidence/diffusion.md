# Domain evidence — diffusion weighting

Curated card. **Not a literature review.**

**Primary source.** Bernstein, King & Zhou, *Handbook of MRI Pulse Sequences*, Elsevier 2004 —
**Ch. 9 "Motion-Sensitizing Gradients", §9.1 "Diffusion-Weighting Gradients"** (X. J. Zhou), book
**pp. 274–281**, with **Figures 9.3 and 9.4** (pp. 278–279) and **Example 9.2** (p. 279).
**Source ID:** `handbook_mri_pulse_sequences_2004` — see [`README.md`](README.md) for resolution.

Read: pp. 274, 278–280. Companion card for the neighbouring section:
[`gradient_moments.md`](gradient_moments.md).

---

## The canonical structure

| claim | page | establishes | does not establish |
|---|---|---|---|
| "A diffusion-weighting gradient typically consists of **two lobes with equal area**. In pulse sequences based on spin echoes, the two lobes have the **same polarity** and are placed at either side of a refocusing RF pulse." | 274 | the PGSE / Stejskal–Tanner structure, and that **the encoding is defined around a pulse that belongs to another component** | that two lobes are the only option — the same page allows more |
| "In gradient-echo-based sequences, however, the two lobes must have **opposite polarity** and are often concatenated." | 274 | that the polarity follows from whether a refocusing pulse intervenes — the 180 does the sign flip | |
| "sometimes referred to as a bipolar gradient (or **Stejskal–Tanner** gradient; Stejskal and Tanner 1965), even though it is **unipolar in the spin-echo case**" | 274 | the naming, and a trap: "bipolar" describes the gradient-echo form | |
| "Waveforms with **more than two lobes** can also be used … as long as they follow the guidelines given in this section." | 274 | that the family is open, and that the b-value integral — not the lobe count — is the contract | which multi-lobe designs are worthwhile |
| "The **amplitude** of the diffusion-weighting gradient is typically **the maximum allowed by the system**, and its pulse width is considerably longer than that of most imaging gradients." | 274 | the realisation convention: pin the amplitude, solve the width | |

## The b-value, and it is an integral over everything

$$b = (2\pi)^2 \int_0^{TE} k(t)\cdot k(t)\,dt, \qquad k(t) = \frac{\gamma}{2\pi}\int_0^t G(t')\,dt'$$

| claim | page | establishes | does not establish |
|---|---|---|---|
| Eqs. 9.5–9.7, above; "used to calculate the b-value for **any gradient waveform in the pulse sequence**" | 278 | that `b` is a property of the **whole waveform**, not of the diffusion lobes — so slice-select lobes, crushers and readout prephasers all contribute. **This is why validation must integrate the tree.** | a decomposition into per-component contributions, which does not exist: `b` is quadratic in `k` |
| Fig. 9.3, spin echo, rectangular lobes: $b = \gamma^2G^2\delta^2(\Delta-\delta/3)$ | 278 | the canonical closed form, with $\delta$ the lobe width and $\Delta$ the **lobe-centre separation** | |
| Fig. 9.3, **trapezoidal** lobes: $b = \gamma^2G^2\left[\delta^2(\Delta-\delta/3) + \varepsilon^3/30 - \delta\varepsilon^2/6\right]$, and the caption: "the expression … **reduces to the value for the rectangular lobes if the ramp duration goes to 0**" | 278 | **the finite-ramp correction, in closed form** — the term the reference implementations omit. It is *negative* in sum, so ignoring it over-estimates `b` | higher-order shapes beyond those figured |
| Fig. 9.4 gives the gradient-echo forms, e.g. rectangular $b = 2\gamma^2G^2\delta^3/3$ | 279 | that the spin-echo and gradient-echo cases have genuinely different expressions | |
| **Example 9.2**: lobes of 25 mT/m separated by 50 ms, to reach b = 1000 s/mm²; "solving the **cubic equation** with respect to the gradient lobe width δ, we obtain **a single physically meaningful solution**, δ = 22.99 ms" | 279 | **that the canonical design problem is a cubic with an explicit physical root** — algebra, not search. An external, published number to test against | that every diffusion design is closed form |

## Practical constraints, and one that is a design default

| claim | page | establishes | does not establish |
|---|---|---|---|
| "To minimize TE, the **maximum gradient amplitude** available on the MRI scanner is usually employed." | 280 | pin amplitude, solve width | |
| "Although using the maximum gradient **slew rate** can also reduce the TE value, the improvement is usually **much less** … because the plateau duration typically greatly exceeds the ramp width. In addition, employing the maximum slew rate contributes to the overall dB/dt, which may induce unwanted **peripheral nerve stimulation** or excessive **eddy currents**. For these reasons, **maximum slew rate is often not used** in diffusion-weighting gradients, especially for imaging human subjects." | 280 | **a physical default with a source**: design the lobes below the slew limit. Not a safety margin someone guessed | what fraction — the source says "often not used", not how much |
| "The large diffusion-weighting gradient amplitude can induce substantial **eddy currents** … An effective way to address this problem is to **break the diffusion-weighting gradient waveform into multiple lobes** so that eddy currents produced by different lobes can partially cancel one another using **two refocusing pulses**" (Reese et al. 2003, §17.2) | 280 | that the twice-refocused design exists and **why** — eddy-current cancellation, not b-value efficiency | its design equations, which are in §17.2 |
| "If the diffusion gradient lobes before and after a refocusing pulse are **identical**, their phase accumulation from the **concomitant field is completely canceled** in spin-echo pulse sequences." | 280 | a reason to make the two lobes identical, and a real advantage of the spin-echo form | the gradient-echo case, where it "cannot be easily cancelled" |
| Typical b "on the order of 1000 s/mm²", with lobes "several tens of milliseconds", which "inevitably leads to a rather long TE", reducing SNR and adding T2 contrast — "**T2 shine-through**" | 279–280 | the protocol scale and the cost | |

## What the source does not determine

- **Minimum-TE design.** It says a long TE is the consequence and how to shorten it (amplitude), but
  gives no minimum-TE construction. That is where the copyleft `PulseqDiffusion` witness works.
- **Twice-refocused designs.** Named, motivated, and deferred to §17.2 and Reese et al. (2003).
- Direction schemes, b-matrices for cross terms, or anything about tensor encoding.
- How the diffusion lobes interact with the imaging gradients' own contribution to `b` — the
  integral in Eq. 9.7 includes them and the figures' closed forms do not.

## Provenance note

Located by a full-text search for "b-value" and "diffusion-weighting gradients", which pointed at
the chapter contents and then §9.1. Only the pages listed at the top were read.
