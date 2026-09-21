# Post-v0 batch scan — the corpus after PR #31 and PR #32

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-21, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../README.md) for which one to edit.

**Date:** 2026-09-21 · **Phase:** B · **Skill:** v0, frozen and unmodified

A batch mining pass with the v0 Module Mining Skill used as a *tool* rather than as the thing
under calibration. It produces a **candidate map**, not implementations. No Module was written,
no production code was changed, and `.claude/skills/module-mining/` was not touched — including
where v0 was awkward, which is recorded in §6 as Phase E evidence instead.

`coarse_scan_2026-09-18.md` is the previous scan and stays as historical evidence. This one
replaces nothing; it re-reads the corpus against **merged `main`**, which now ships
`RadialReadout`, `SpiralReadout`, `SaturationPrep` and `GRE3DTR`.

**Headline: very few genuinely new leaves, and one architectural finding that is worth more than
a Module.** A published 4D-flow implementation does not construct its velocity-encoding gradients
at all — it states a *moment requirement* and calls an external constrained optimiser. That is the
regime SeqCraft cannot currently express, and it is now observed rather than hypothesised.

---

## 1. Corpora actually scanned

| corpus | commit | licence | role used | scanned |
|---|---|---|---|---|
| `pypulseq` | `f2c582b` | MIT | discovery, design-witness | `examples/scripts/` — 19 scripts |
| `pulseq` (MATLAB) | `2fd6ab6` | MIT | discovery, design-witness | `matlab/demoSeq/` — 42 scripts |
| `MRzero-Core` | `5732d86` | AGPL + EULA | discovery, **oracle only** | `documentation/playground_mr0/` — 26 notebooks |
| `openmrf-core-matlab` | `fd03604` | MIT | discovery, design-witness | preparations, readouts, `main_sequences/` |
| `fmrifrey/lps` | `552d87a` | MIT | design-witness | 3D Cartesian GRE |
| `pulseq-cest` | `0756666` | MIT | discovery | `seq-generation/` |
| `PulseqDiffusion` | `0517ae1` | AGPL-3.0 | **discovery only** | file names and solver *structure*; no implementation text read for adaptation |
| `rtspiral_pypulseq` | `5d32a13` | GPL-3.0 | **discovery only** | not pursued — spiral shipped |
| `Open4DFlow` | `9e752d8` | GPL-3.0 | **discovery only** | physics and requirement structure only |

`openmrf-core-matlab` and `lps` matched the commits the registry already recorded. The other six
were cloned for this scan; `pypulseq`, `pulseq` and `MRzero-Core` were already checked out at the
registry's recorded commits.

**Verified rather than assumed:** 21 of 26 MRzero playground notebooks import `pypulseq`, exactly
as `sources.yaml` states. Its `design-witness` exclusion holds.

---

## 2. Coverage map against merged `main`

Merged `main` ships fifteen public items: `CartesianLine`, `EPI2D`, `Excitation`, `FSE2D`,
`GRE2D`, `GRE2DTR`, `GRE3DTR`, `IRPrep`, `PhaseEncode`, `RadialReadout`, `Refocusing`,
`SaturationPrep`, `SpiralReadout`, `TSEShot`, `spoiler`. Five more abstractions live **only** in
example notebooks — `MPRAGE2D`, `MP2RAGE2D`, `SE2D`, `GREEPI2D`, `SEEPI2D` — and two notebooks
(`fse_2d/01`, `gre_2d/01`) define a class that has since shipped, which is where `FSE2D` and
`GRE2DTR` came from rather than a second copy.

| family | sources | physical role | current coverage | classification |
|---|---|---|---|---|
| GRE / FLASH | pypulseq, pulseq, OpenMRF, MRzero | kernel + imaging | `GRE2DTR`, `GRE2D` | `DIRECT_SHIPPED_DUPLICATE` |
| 3D Cartesian GRE | pulseq, lps | kernel | `GRE3DTR` | `DIRECT_SHIPPED_DUPLICATE` |
| Multi-echo GRE | pypulseq, pulseq | readout policy | `megre_2d/`, which defines **no class** | `COMPOSITION_COVERED` |
| SE / CPMG | pulseq, MRzero | kernel | `Refocusing` + `SE2D` in notebook | `NOTEBOOK_ONLY_EXISTING` |
| TSE / FSE / HASTE / RARE | all four | kernel | `TSEShot`, `FSE2D` | `DIRECT_SHIPPED_DUPLICATE` |
| EPI (GRE and SE) | all four | readout | `EPI2D` | `DIRECT_SHIPPED_DUPLICATE` |
| MPRAGE / IR-FLASH | pypulseq, pulseq | imaging policy | `IRPrep` + `mprage_2d/`, `mp2rage_2d/` | `COMPOSITION_COVERED` (unchanged from v0) |
| Radial | pypulseq, pulseq, OpenMRF | readout | `RadialReadout` | `DIRECT_SHIPPED_DUPLICATE` |
| Spiral | pulseq, OpenMRF, rtspiral | readout | `SpiralReadout` | `DIRECT_SHIPPED_DUPLICATE` |
| Fat saturation | pulseq, OpenMRF | preparation | `SaturationPrep` | `DIRECT_SHIPPED_DUPLICATE` |
| Inversion | OpenMRF, pypulseq | preparation | `IRPrep` | `DIRECT_SHIPPED_DUPLICATE` |
| Crusher / spoiler | OpenMRF, all | primitive | `spoiler` | `DEGENERATE_CASE_OF_EXISTING_MODULE` |
| FID | pulseq, MRzero | — | plain composition | `PRIMITIVE_COMPOSITION_COVERED` (unchanged from v0) |
| Selective RF | pulseq | — | `Excitation` | `DEGENERATE_CASE_OF_EXISTING_MODULE` |
| **Cine / cardiac gating** | pulseq | imaging policy | trigger and digital-output **events already compile** | `PRIMITIVE_COMPOSITION_COVERED` |
| GRAPPA / calibration sampling | pulseq, pypulseq | imaging policy | `gre_epi_2d/` §7 | `COMPOSITION_COVERED` |
| Labels (`LIN`/`REV`/`NAV`…) | pulseq, pypulseq | annotation | `EPI2D` emits them | `PRIMITIVE_COMPOSITION_COVERED` |
| UTE | pulseq, pypulseq, OpenMRF | readout | none | **UNRESOLVED** |
| ZTE / PETRA | pulseq | readout + encoding | none | **UNRESOLVED** |
| PROPELLER / TSEprop | pulseq | encoding policy | none | **UNRESOLVED** |
| Rosette / cones | OpenMRF | readout | none | **UNRESOLVED** |
| PRESS / semi-LASER | pulseq, OpenMRF | localisation kernel | none | **UNRESOLVED** |
| **bSSFP / TrueFISP** | pulseq, OpenMRF, MRzero | kernel? | none | **UNRESOLVED** |
| **T2 preparation** | OpenMRF | preparation | none | **UNRESOLVED** |
| **Spin-lock / T1ρ** | OpenMRF | preparation | none | **UNRESOLVED** |
| CEST | OpenMRF, pulseq-cest | preparation | none | **UNRESOLVED** |
| Three-axis crusher policy | OpenMRF | primitive | `spoiler`, one axis at a time | `DEGENERATE_CASE_OF_EXISTING_MODULE` + `ARCHITECTURE_REVISIT_CANDIDATE` |
| **Diffusion weighting** | pulseq (MIT), PulseqDiffusion (AGPL), MRzero | encoding + kernel coupling | none | **UNRESOLVED** |
| **Velocity / flow encoding** | Open4DFlow (GPL) | encoding | none | **UNRESOLVED** |
| Flow compensation | Open4DFlow (GPL) | encoding | none | **UNRESOLVED** |
| MRF | OpenMRF | sequence-wide schedule | none | **UNRESOLVED** |
| STEAM / stimulated echo | MRzero, OpenMRF | unresolved | none | **unchanged from v0 — evidence-blocked** |
| DREAM | MRzero | unresolved | none | **unchanged from v0 — deferred behind STEAM** |

**Three classifications were measured rather than assumed.**

Multi-echo GRE: `megre_2d/01_build.ipynb` defines no class at all, so there is no abstraction to
be notebook-only. Composition covered.

`spoiler` / crusher: OpenMRF's `CRUSH_x_y_z` takes per-axis twist counts, which is what `spoiler`
already expresses one axis at a time — so, degenerate rather than duplicate. It carries
`ARCHITECTURE_REVISIT_CANDIDATE`, and the reason is sharper than "it is one object there":

```text
OpenMRF derates every crusher axis to  maxGrad * lim_grad  with lim_grad = 1/sqrt(3)
and STAGGERS the three ramps in time, each delayed by the previous axis's rise time.
```

Both mechanisms exist to keep the **vector** magnitude within the per-axis limit when three axes
crush at once. That is the same question Phase A left open for `SpiralReadout`'s prephaser and
rewinder, recorded there as a design-policy revisit with no evidence either way — and this is
evidence: an independently developed implementation has a deliberate policy for it, and two
distinct mechanisms to implement it. **The Phase A item should be updated to say a witness
exists.** Not done here; Phase B does not edit other records.

Cine: I expected a gap, because `sc.trigger` does not exist. It is not a gap —
`design/events.py` has `POINT_KINDS = {'labelset', 'labelinc', 'trigger', 'output'}`, `EPI2D`
already emits labels, and a trigger plus a digital output plus an excitation compiles to four
blocks. Plain composition is the intended API, and that is the answer rather than a missing
Module.

---

## 3. Unresolved candidate backlog, with reasoning classes

Reasoning classes are **scan annotations**, not Skill vocabulary, and are not added to the schema.

| candidate | reasoning class | shared contract hypothesis | the ownership question | witnesses | fine-scan value |
|---|---|---|---|---|---|
| **Velocity / flow encoding** | *joint waveform solve before realisation* | moments at semantic times, solved under hardware limits | who owns M1, when the imaging gradients have already contributed some of it | 1, GPL, discovery only | **highest — architectural** |
| **Flow compensation** | *joint waveform solve before realisation* | same requirement model as above with M1 → 0 | whether it is the same contract as velocity encoding or a coincidence of bipolar shapes | 1, GPL, discovery only | high, jointly with the above |
| **Diffusion weighting** | *joint moment solve* + *cross-leaf kernel* | b-value, timing windows and TE are mutually determined | who owns TE | 2 — one MIT, one AGPL, **and they solve opposite directions** | **high** |
| **bSSFP** | *steady-state / sequence-wide semantic coupling* | balanced moments per TR, TE = TR/2, RF phase cycling, catalyzation | the rewinder of TR *n* and the prephaser of TR *n+1* are one waveform — the repetition boundary cuts through it | 2 independent (pulseq, OpenMRF) + MRzero as oracle | **high** |
| **T2 preparation** | *ordinary reusable leaf* | tip-down, refocusing train, tip-up, crush — magnetisation returns to +z weighted by exp(−T<sub>prep</sub>/T2) | whether the refocusing train is the leaf's or a caller's | 1 (OpenMRF); pulseq-cest adjacent | **medium — the conventional control** |
| **Spin-lock / T1ρ** | *coherence-pathway problem* | lock along B1 for a defined time, with a B0/B1-compensating pattern | six named variants (SSL/hSSL/qSSL/RESL/CSL/PSCSL) — one physical solve or six? | 1 (OpenMRF) | medium |
| **CEST** | *ordinary reusable leaf*, probably | a saturation train at an offset, then readout | relation to shipped `SaturationPrep` — extension or different contract | 2 (OpenMRF, pulseq-cest) | medium |
| **UTE** | *ordinary reusable leaf* | centre-out readout with ramp sampling and half-pulse excitation | whether `RadialReadout(partial_fourier=0.5)` already is it | 3 | medium |
| **ZTE / PETRA** | *possible second sequence DSL* | encoding begins before excitation; k-space filled by two different mechanisms | whether one abstraction covers the shell and the Cartesian centre | 1 | low now |
| **PROPELLER** | *composition / imaging policy only* | rotated Cartesian blades | likely a schedule, like radial angles | 1 | low |
| **Rosette / cones** | *ordinary reusable leaf* | another parametrised non-Cartesian path | whether `SpiralReadout`'s path/traversal split generalises | 1 | low — but cheap, and it tests a v0 claim |
| **MRF** | *steady-state / sequence-wide semantic coupling* | per-repetition flip and TR schedules, hundreds long | whether a schedule is data (like radial angles) or a contract | 1 | medium |
| **PRESS / semi-LASER** | *cross-leaf kernel* | three orthogonal selective pulses intersecting in one voxel | localisation as a joint property of three pulses | 2 | low now |
| **STEAM / DREAM** | *evidence-blocked physical boundary* | unchanged | unchanged | unchanged | **do not re-run** |

**STEAM is unchanged and was not reopened.** DREAM's two playground notebooks differ in which
stimulated-echo pathway is read first (`STE` against `STID`), which is consistent with the v0
finding — the variation is in coherence-pathway ordering — and adds no new witness, because both
import PyPulseq. Nothing here changes the v0 verdict.

---

## 4. The architectural finding

`Open4DFlow/sequences/4Dflow.py` computes M0 and M1 analytically for the imaging gradients that
are already there, forms a **residual moment requirement** per axis — the desired M1 minus what
the readout and slice gradients already contribute — and then calls `gropt`, an external
constrained gradient optimiser, with a parameter block of moment targets, tolerances, `gmax`,
`smax` and a time budget. It caches the result, keyed by the moment pair, because the solve is
expensive.

```text
base sequence moment state   (M0, M1 already accumulated)
+ desired condition          (venc -> target M1; flow compensation -> M1 = 0)
+ available timing           (TE)
+ hardware limits            (gmax, smax)
        |
        v  numerical solve
    waveform
```

This is the shape `post_v0_roadmap.md` refused to build speculatively, observed in a published
implementation rather than hypothesised. **No SeqCraft Module can express it today**: a Module
computes its events in closed form from its own parameters, and this requirement is stated
*relative to what the rest of the sequence already does*.

Diffusion is the same class seen from the other side. The MIT Pulseq demo inverts the analytic
b-formula for amplitude **given** timing windows that TE has already fixed, and its own comment
says the formula is for rectangular gradients with trapezoids marked `TODO` — so the shipped
physics ignores the ramps. `PulseqDiffusion` holds the harder direction, `opt_TE_bv_SE` and
`opt_TE_bv_TRSE`: minimum TE for a target b-value, single- and twice-refocused. **That direction
is AGPL and discovery-only.**

---

## 5. Reference and licensing gaps

```text
velocity / flow      ONE witness, GPL, discovery only.  Any candidate would rest on
                     analytic re-derivation plus a permissive second implementation that
                     does not currently exist in the registry.

diffusion            TWO witnesses and they are NOT interchangeable.  pulseq's
                     writeEpiDiffusionRS.m is MIT and usable as a design-witness -- the
                     registry lists diffusion only under the AGPL PulseqDiffusion entry,
                     which understates what is available.  REGISTRY CORRECTION, recorded
                     here rather than applied, because sources.yaml is not Phase B's to edit
                     mid-scan.

gropt                NOT in the registry, present locally, GPL-3.0.  It is the solver a
                     real flow implementation delegates to, so it is discovery evidence
                     about the shape of the problem.  Adding it is a registry decision.

T2prep, spin-lock,   ONE witness each (OpenMRF).  Independent of the Pulseq line, which
MRF, rosette         makes it a genuine witness -- but a single one, and v0's rule is that
                     one witness sets the evidence tier, not the abstraction boundary.

bSSFP                TWO genuinely independent witnesses, and they DISAGREE on
                     catalyzation (alpha/2 prep against a flip-angle ramp) and on whether
                     TR is an input or an output.  That disagreement is the most useful
                     evidence in the scan: it is about the boundary, not about a protocol
                     parameter.
```

---

## 6. Where v0 was insufficient — Phase E evidence, not Phase B fixes

Recorded as observations. The Skill was not edited.

**1. The coverage vocabulary has no cell for "the events exist and the composition is the API,
but the *design problem* is unsolved."** Velocity encoding is not
`PRIMITIVE_COMPOSITION_COVERED` — a user cannot write it, because the waveform has to be solved
for. It is not `UNRESOLVED / NEW CANDIDATE` in the ordinary sense either, because the missing
thing is not a Module. v0 made me choose between two labels that both mislead.

**2. The ownership test asks "can this component determine the value from its own contract and
parameters?" and assumes the answer decides *placement*.** For a moment requirement the answer is
no for every component, including the kernel — the value depends on hardware limits and the time
budget, which belong to the compiler. v0's test correctly says "not this leaf" and has nothing to
say about what to do next.

**3. There is no vocabulary for a candidate whose contract is a *requirement* rather than an
event.** Rules A–F are all about emitted events: the mode contract, inspecting what is emitted,
measuring on the emitted lattice, occupancy after placement. A candidate that produces a
constraint to be solved later passes none of them vacuously and fails none of them meaningfully.

**4. `DEGENERATE_CASE_OF_EXISTING_MODULE` + `ARCHITECTURE_REVISIT_CANDIDATE` worked exactly as
designed** on the three-axis crusher, which is worth recording as a success rather than a gap.

**5. What a human reasoned about outside the Skill:** whether two families that look alike
(velocity encoding and flow compensation) share a *requirement model* rather than a waveform
shape. v0 says waveform similarity is discovery only — correctly — but offers no positive test for
"same requirement, different target value."

---

## 7. Proposed Phase C shortlist

Three, chosen to be complementary rather than likely to go GREEN.

| | candidate | what it tests | why this one |
|---|---|---|---|
| **C1** | **T2 preparation** | *does the ordinary path still work?* | The conventional control. A preparation leaf beside `IRPrep` and `SaturationPrep`, in a folder whose contract already names the third seat. If v0 cannot produce a clean GREEN here, that is a much louder signal than any of the hard cases. One witness, so the evidence tier is capped and the claim scope must say so — which is itself a v0 rule worth exercising again. |
| **C2** | **bSSFP** | *where does the kernel boundary actually fall?* | Two independent witnesses that agree on the physics and disagree on the boundary. The rewinder/prephaser pair spans the repetition boundary, so a `…TR` kernel shaped like `GRE2DTR` may be the wrong unit — and `post_v0_roadmap.md` names `BalancedSSFPTR` without evidence. Expected outcome is genuinely open, including `ARCHITECTURE_REVISIT_CANDIDATE` against the existing kernel layer. |
| **C3** | **Velocity encoding + flow compensation, scanned together** | *can SeqCraft express a requirement at all?* | The architectural stress case, and the two are scanned as one fine scan precisely to test whether they share a requirement model or merely a silhouette. A likely outcome is "no Module, and here is what would have to exist first", which is the most valuable result available. Licence-constrained to analytic re-derivation, which must be stated in the claim scope up front. |

**Diffusion is deliberately not in the shortlist**, despite being high value. It couples to EPI,
TE ownership, crushers and hardware limits simultaneously, and C3 tests the same *reasoning class*
with fewer confounds. It becomes the obvious next candidate once C3 has said whether requirements
can be expressed at all. Recorded so the omission is a decision rather than an oversight.

---

## 8. Recommended phase ordering

One change, and the reason is C3.

```text
current    B (scan) -> C (fine scans) -> D (?) -> E (Skill evolution)

proposed   B -> C1 + C2 -> E-lite -> C3 -> E
```

Run **C1 and C2 first**, then a *small* Skill evolution pass on what they actually produce, then
C3 — because C3 is the candidate most likely to be misrecorded by v0's vocabulary, and running it
against an unrevised Skill risks the finding being lost in a record that cannot hold it. §6 items
1–3 are already enough to justify a bounded vocabulary extension, and they came from a scan rather
than from design.

**This is a recommendation, not a decision.** The alternative — run all three fine scans against
frozen v0 and revise once, with more evidence — is defensible and keeps the "observe before
changing" discipline stricter. The trade is between cleaner evidence and a record that can hold it.
