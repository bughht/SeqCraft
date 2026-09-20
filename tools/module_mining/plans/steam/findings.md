# StimulatedEcho / STEAM — Boundary Analysis

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Stage:** prospective supervised run. Outcome open at the start and at the time each section
below was written.

---

## 1. What the three uses actually are

Read off the emitted block structure, not off the notebook prose.

### 1.1 `mr0_STE_3pulses_5echoes_seq` — a coherence-pathway demonstration

```text
rf1 (90) -> adc -> delay 5 ms -> rf2 (90) -> adc -> rf3 (180) -> adc -> adc
```

**No gradients are played at all.** `gx`, `gx_pre` and `gspoil` are constructed in the notebook
and never added to a block. There is no spatial encoding, no storage interval that anything is
designed around, and no readout in the imaging sense. Its purpose is to show the five echoes three
pulses produce.

Note the flips: **90 / 90 / 180**. This is not a classical STEAM (90/90/90).

### 1.2 `mr0_diffusion_prep_STEAM_2D_seq` — a prepared-and-stored magnetisation module

```text
rf90 + Gz(slice)
Gdiff(y) + Gz_rephase
rf180  (phi = 90)
Gdiff(y) + Gx_store
rf90   (phi = 180)          <- tips back to +z
Gz_spoil                    <- destroys residual transverse
then: a train of 5 deg GRE readouts (prephase + phase encode, ADC, spoil + rewind)
```

Flips: **90 / 180 / 90**. The pulse in the middle is a refocusing pulse, not a second storage
pulse, so what is stored is a *refocused, diffusion-weighted* state.

One detail carries real architectural weight: `Gx_store` has area
`-(gx.area/2 + gx.area * spoil_factor)` — the store gradient **also carries the readout's
prephaser**. The preparation and the readout are not independent.

### 1.3 `mr0_DREAM_STE_seq` — two preparation pulses and a two-pathway readout

```text
rf1 (alpha, block, non-selective)
Gx_m2                       <- dephasing moment, area -gx_read_amp*2*t0, ~11*t0 long
rf2 (alpha, block, phi=180)
Gx_spoil
then: a train of alpha_readout sinc pulses, each followed by ONE extended-trapezoid x waveform
      that prephases, reads and spoils, with an ADC of 2*Nread samples over 4*t0
```

Flips: **α / α / (readout train)**. Only **two** preparation pulses. The readout ADC is twice the
usual length and the x waveform is a single shaped gradient, because **both the stimulated echo
and the FID must land inside one acquisition window at known positions** — that is what makes the
flip-angle map computable from their ratio.

---

## 2. The finding: waveform similarity is not contract similarity

The coarse scan grouped these as "the same three-pulse structure appearing three times". Measured,
that grouping does not survive.

| | pulses before readout | flips | gradient during the interval | what is stored | readout coupling |
|---|---|---|---|---|---|
| 1.1 | 3 | 90 / 90 / 180 | **none** | nothing; all pathways observed | no readout |
| 1.2 | 3 | 90 / 180 / 90 | diffusion, both sides of the 180 | a refocused, diffusion-weighted state | store gradient **carries the readout prephaser** |
| 1.3 | 2 | α / α | one dephasing moment | a modulated longitudinal state | dephasing moment **sets echo separation inside one ADC window** |

Three different pulse counts, three different flip patterns, three different gradient roles and
three different couplings to what follows. What they share is *Bloch/EPG behaviour* — that three
RF pulses generate a stimulated echo — which is a fact about spin physics, not a contract a module
can own.

> **The waveform grouping does not establish a boundary.**
>
> A module owns a value somebody could get wrong. Here the values that could be got wrong are
> different in each case: a diffusion moment balanced across a refocusing pulse; a dephasing
> moment that positions one echo relative to another inside a shared ADC window; and, in the
> first case, nothing at all.

This is the reasoning class this run exposes, and it generalises beyond this candidate: **a family
identified by its waveform pattern is not thereby a candidate.** The coarse scan's grouping
criterion — "the same three pulses keep appearing" — is a discovery heuristic, and discovery
heuristics are not boundary evidence. Every prior pilot was grouped by a *physical solve* (a
shared crusher window, a trajectory, a signed moment coupling), not by a waveform silhouette.

### What this does **not** settle

It rejects the **broad, waveform-defined** candidate. It says nothing about whether
stimulated-echo physics can be expressed as a reusable abstraction with a clear correctness
contract, and that question stays open:

```text
broad waveform-defined STEAM abstraction      -> unsupported by this evidence
physics-defined StimulatedEcho / STEAM        -> STILL UNRESOLVED
```

A future abstraction could plausibly own semantics such as the requested coherence pathway, the
storage and mixing intervals, a semantic stimulated-echo time, the required signed gradient-moment
relationships, RF-centre timing relationships, and the suppression of unwanted pathways. Whether
those can be parameterised **while retaining clear correctness conditions** is the real question.

If supporting the general case instead requires the caller to specify nearly the whole
RF / gradient / timing / coherence program, the abstraction has become a **second sequence DSL**,
and the physics is better served by lower-level semantics and helpers plus family-specific
kernels. That is the test a later scan should apply.

**This phase does not resolve it**, and nothing here should be read as "stimulated echo should not
be a Module".

## 3. Applying the ownership test

> If a component cannot determine the correct value using only its own physical contract and
> parameters, it should not own that event.

**A hypothetical `STEAMPrep` leaf** — two pulses, an interval, a dephasing moment, a third pulse —
can determine its own flip angles and its own dephasing area. It cannot determine `Gx_store`'s
readout-prephaser term (1.2) or the dephasing moment that fixes echo separation in the readout
window (1.3), because both depend on the readout. So a leaf covering all three uses would either
own values it cannot compute, or exclude precisely the values that make each use work.

**A hypothetical `STEAMKernel`** — preparation plus readout solved together — could own those
couplings. But the two couplings are *different couplings*: one is a prephaser folded into a store
gradient, the other is an echo-separation constraint inside a doubled ADC window. A kernel
covering both would be a switch between two unrelated solves, which is the
`GRE2DTR(balanced=True, spoiled=False, ...)` shape already ruled out for bSSFP.

**Against `IRPrep` (`EXTEND_EXISTING`).** `IRPrep` inverts and waits. The structural similarity —
prepare, wait, read — is real, but `IRPrep` owns no coupling to the readout, and every candidate
here is interesting *because* of its readout coupling. Generalising `IRPrep` to cover them would
mean adding readout-dependent terms to a module that currently has none.

## 4. What is genuinely reusable, and where it already lives

Stripping the application-specific parts of 1.2 and 1.3 leaves: an excitation, a refocusing or a
second small-flip pulse, a delay, a gradient whose area is set by something outside the
preparation, a store pulse, and a spoiler. SeqCraft already ships `Excitation`, `Refocusing`,
`spoiler` and `PhaseEncode`, and a delay is a delay. The composition is not currently shortened by
any class that could be extracted without also owning the readout coupling.

## 5. What this run did **not** establish

- **That no STEAM module could ever exist.** Only that the three uses in this corpus do not share
  one. A narrower candidate — the 90/180/90 diffusion-prepared store of 1.2 alone — may well be a
  real contract. It has **one consumer**, which is the same objection recorded for a
  single-consumer reconstruction contract: one consumer's implementation under a more general name.
- **Anything about DREAM's own contract.** §1.3 observed a cross-component coupling — the
  dephasing moment and the readout waveform appear to be solved together so that two pathways land
  in one window at known positions. **That is candidate evidence, not a verdict.** It was read off
  one notebook's block structure; whether the same solve is reused across DREAM variants (`STE`
  and `STID` differ in which pathway is read) is unexamined, and nothing here establishes that the
  coupling survives once shared stimulated-echo semantics are understood.

  So it strengthens `DREAM -> DEFER / revisit later`. It does **not** justify
  `DREAM -> NEW_KERNEL`. A later DREAM fine scan asks what remains after the shared semantics are
  known.
- **That the physics is correctly described.** Everything above is read from block structure and
  gradient areas in one non-independent corpus. No measurement was performed, no reference was
  executed by us, and no simulation was run.
