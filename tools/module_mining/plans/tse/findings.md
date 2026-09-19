# TSE/FSE Fine Scan — Code Archaeology and Semantic Decomposition

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** TSE/FSE shot
**Date:** 2026-09-18
**Pass:** 1 — code archaeology; 2 — semantic decomposition
**Status:** evidence only. No candidate API is proposed here, and no boundary decision is taken.

References and provenance: [`reference_inventory.md`](reference_inventory.md).
Measurable contract: [`invariant_table.md`](invariant_table.md).

---

## 1. What each reference lets the user configure

| | R1 SeqCraft `FSE2D` | R2 `write_tse.py` | R3 `writeTSE.m` |
|---|---|---|---|
| geometry | `fov_mm`, `matrix`, `thickness_mm` | `fov`, `n_x`, `n_y`, `slice_thickness` | `fov`, `Nx`, `Ny`, `sliceThickness` |
| train | `echoes` | `n_echo` | `necho` |
| contrast timing | `echo_spacing_s` (optional) | `te` (= first echo spacing) | `TE1`, `TEeff` |
| repetition | `tr_s` | `tr`, `n_slices` | `TR`, `Nslices` |
| RF | `flip_deg`, `refocus_flip_deg`, durations, `refocus_thickness_factor` | `rf_flip_deg` (scalar or per echo), durations fixed in body | `rflip` (vector), durations in body |
| readout | `bandwidth_hz_px`, `partial_fourier` | `sampling_time` in body | `samplingTime` in body |
| spoiling | `crush_cycles_slice`, `crush_cycles_readout`, `spoil_cycles_per_voxel`, `spoil_axis` | `fsp_r`, `fsp_s` in body | `fspR`, `fspS` in body |
| ordering | **`segments`, a list of line lists, passed to `build`** | derived inside from `n_echo`/`n_y` | derived inside, then `circshift` by `TEeff/TE1` |

The first structural difference is visible here and it is the largest one: **R2 and R3 own the
phase-encode table; R1 takes it as data.** R3's own comment calls its reordering
`% TSE echo time magic`, and the shift it applies exists only to move the central line to a
chosen echo index. That is acquisition policy expressed as arithmetic, which is exactly the class
of thing the playbook says to leave with the caller.

## 2. Which calculations encode MRI physics

Identical in intent across all three, differing only in direction of derivation:

1. **Readout amplitude and area** from `1/FOV`: `deltak = 1/fov`, `kWidth = Nx*deltak`
   (R2/R3); R1 derives the same from `fov_mm`/`matrix`/`bandwidth_hz_px` inside `CartesianLine`.
2. **Dephasing before the first echo** — half the readout area, plus whatever readout-axis
   crusher is used:
   - R2/R3: `agr_preph = gr_acq.area/2 + agr_spr`.
   - R1: `dephase_area_per_m = lobe_in.area + ro.area_to_echo_per_m`, which with
     `crush_cycles_readout = 0` reduces to `gx.area / 2`.
   Same identity, written from opposite ends — R2/R3 from the designed trapezoid, R1 from the
   readout's own integrated area to the echo sample, which is what makes partial Fourier and the
   half-dwell sample offset fall out of it.
3. **Crusher balance about the refocusing pulse.** All three place equal gradient area on either
   side of the refocusing RF centre, on every axis. R2/R3 achieve it by construction from
   mirrored extended trapezoids; R1 achieves it by integrating the waveform it emits and solving
   the two crusher amplitudes (`Refocusing._design_fused`).
4. **Slice rephasing.** R2/R3 have **no** excitation rephaser: the excitation's post-centre
   plateau area is cancelled by the refocusing pulse's own trailing crusher, which is why
   `GSspr.area = AGSex*(1+fspS)` while `GSspex.area = AGSex*fspS` — the extra `AGSex` *is* the
   rephasing. R1's `Excitation` always emits its own rephaser (`gzr`), so `k_z` is already zero
   when the train starts and the crushers are symmetric.
   **Both reach `k_z = 0` at the echo; they are not the same waveform and they do not crush the
   same amount.** See §6, disagreement D1.
5. **Echo spacing.** R2/R3 invert it: the user gives `TE1`, and the crusher window follows as
   `tSp = 0.5*(TE1 - readoutTime - tRefwd)`. R1 goes forward: three axes each state a minimum
   window (`crush_min_s`), the maximum wins, and `echo_spacing_s` is the result — with an optional
   requested spacing inverted back onto the window (`+d` of window buys exactly `+2d` of spacing).
   Nothing in R2/R3 checks that its derived `tSp` is long enough for the phase-encode area or for
   the crusher it asked for; a short `TE1` silently produces an over-amplitude trapezoid or a
   `makeTrapezoid` failure.
6. **Echo at the midpoint between its two refocusing pulses.** R2/R3 get this from the same
   mirrored construction. R1 spends a named quantity on it, `shift_s`, because its echo spacing
   can be set by a bound *other* than the train (see §4), and without the shift the echo drifts
   off the midpoint — the notebook reports 1.23 ms for HASTE.

## 3. Which code only exists because Pulseq is flat

R3's `GS1 ... GS7` and `GR3 ... GR7`, and R2's identical `gs1 ... gs7`, `gr3 ... gr7`: nine
extended trapezoids hand-cut so that each pulseq block holds one gradient per axis, with
`amplitudes=[... GSex.amplitude, GSspex.amplitude, ..., GSref.amplitude]` carrying continuity
across the seams by hand. Roughly a third of both scripts is this.

R1 emits continuous waveforms and `sc.compile` cuts them, carrying `first`/`last` across each
seam. `Refocusing`'s module docstring names this explicitly as the reason the module emits one
event rather than seven.

**Classification: `REFERENCE_ARTIFACT`.** None of `GS1..GS7` is a candidate API concept, and
block count must not be an equivalence criterion in the comparator (playbook §9 L2).

Two more artifacts of the same kind:

- R2/R3's fixed `dG = 250e-6` ramp, commented "makes sequence structure much simpler" — a slew
  limit frozen at one system, which R1 derives.
- R2/R3's `gr_spex`/`GRspex` is computed and never used (R3 line `GRspex = ...`, then `GR3=GRpreph`).

## 4. Which choices are acquisition policy

- **the phase-encode table** — `pe_steps`, the even-`necho` `circshift`, and the `TEeff` shift
  (R2/R3 internal; R1 `segments`);
- **ordering within the table** — linear / interleaved / centric (R1 demonstrates three from one
  instance; R2/R3 have one);
- **dummy shots** — all three use exactly one, differently (D3 below);
- **slice ordering / `n_slices`** — R2/R3 loop slices inside the shot loop and set
  `rf.freq_offset` per slice; R1 takes `center_mm` per shot and leaves the loop to the caller;
- **which echo carries the centre line**, and therefore the effective TE.

R1 additionally exposes two *derived* policy diagnostics rather than choices: `echo_bands` and
the warning about an echo index scattered across `ky`. That is a property of the table, computed
by the module, owned by the caller.

## 5. Semantic decomposition

The hypothesis in the playbook survives, with one correction: the initial dephaser is not a
separate stage, it lives **inside the first refocusing pulse's leading crusher window**.

```text
Excitation  (90, +y phase; R1 also rephases the slice here)
      |
      |  x-axis dephaser, area = readout area to echo (+ readout crusher)
      |  placed to END where the refocusing plateau starts
      v
 +--------------------------------------------------+
 |  Refocusing (180, +x phase) + crusher pair,       |
 |     balanced in area about the RF effective centre|
 |  phase encode blip            (after the pulse)   |   x echoes
 |  spin-echo Cartesian readout  (no prephaser)      |
 |  phase rewind                 (before next pulse) |
 +--------------------------------------------------+
      |
      v
 tail: spoilers on x and z, then TR fill
```

Per-component classification:

| component | R1 | class |
|---|---|---|
| 90 excitation + slice rephase | `Excitation` | `EXISTING_SEQCRAFT` |
| initial x dephaser | `pp.make_trapezoid` in `FSE2D.__init__` | `CANDIDATE_OWNED` — it is the train's first-interval arithmetic |
| 180 + crusher pair, balanced about effective centre | `Refocusing` | `EXISTING_SEQCRAFT` |
| phase encode / rewind | `PhaseEncode(line=..)` / `(rewind=True)` | `EXISTING_SEQCRAFT` |
| spin-echo readout | `CartesianLine(prephase=False)` | `EXISTING_SEQCRAFT` |
| readout-axis lobe pair (`lobe_in`/`lobe_out`) | `pp.make_trapezoid` in `FSE2D` | `CANDIDATE_OWNED` — the `skew` identity is train physics |
| one crusher window shared by z, x, y | `crush_min_s` / `crush_s` | `CANDIDATE_OWNED` |
| echo spacing, `shift_s`, `t_refoc_s` | `FSE2D.__init__` | `CANDIDATE_OWNED` |
| tail spoiler | `sc.modules.spoiler` | `EXISTING_SEQCRAFT` |
| TR fill, shot stacking, dummies | `FSE2D.build` | `CALLER_POLICY` (imaging level) |
| `segments`, `te_eff_s`, `echo_bands` | `FSE2D` methods taking the table | `CALLER_POLICY` + derived diagnostics |
| `GS1..GS7`, `GR3..GR7`, fixed `dG` | — | `REFERENCE_ARTIFACT` |

**The natural reusable unit, provisionally: one complete shot — excitation plus echo train.**
The evidence is that the excitation cannot be separated from the train without losing arithmetic:
the first interval's length, the dephaser's area *and* its placement, and `shift_s` all couple the
excitation's effective centre to the first refocusing centre. All three references agree on this
coupling. Whether that unit belongs in `kernel/` and what `FSE2D` retains is deferred to the
boundary decision; the reusable-unit question is not the same as the folder question.

## 6. Where the references disagree

Recorded, not resolved.

**D1 — where the slice rephasing lives.** R2/R3 fold it into the refocusing crusher
(`AGSex*(1+fspS)` vs `AGSex*fspS`); R1 rephases in `Excitation` and crushes symmetrically.
Consequence: at `fspS = 0.5` the reference's *effective* crusher area on z differs before and
after the pulse in **meaning** even though it balances in **area**, and the amount of FID
spoiling per echo is not the same number in the two designs. Both satisfy `k_z = 0` at the echo.
This is a genuine design difference, not a representation difference, and it must not be
averaged away.

**D2 — readout-axis crusher.** `fspR = 1.0` in R2/R3, `crush_cycles_readout = 0.0` in R1. R1's
`Refocusing` docstring argues the readout-axis crusher has *no* k-space consequence (the
conjugation cancels what the two lobes share) and is paid for in echo spacing. R2/R3 pay ~2 x
`gr_acq.area` of window per echo for it. Neither is wrong; the difference is a contrast/spoiling
choice with a measurable timing cost, and it is a parameter, not a boundary question.

**D3 — what a dummy shot plays.** R2/R3 run the `kex=0`/`i_excitation=0` shot with
`phase_area = 0` on every echo; R1 replays segment 0's actual lines without acquiring. The
gradient history entering shot 1 is therefore different. R1's choice makes the dummy's eddy
current and diffusion history match a real shot; R2/R3's keeps `ky` at zero throughout. Worth a
measurement, not a rewrite.

**D4 — labels.** R2 emits none. R3 calls `seq.autoLabel(...)` after the fact, with
`mirrorFourier` and slice sorting. R1 emits `LIN` per echo and `SEG` per shot during
construction and deliberately omits `ECO`, on the argument that every echo of an FSE shot
belongs to one image. Label state is therefore **not** comparable like-for-like across R1 and R3.

**D5 — the plateau symmetry assumption.** This is the disagreement with teeth, and it is a
disagreement between R1 and *both* references:

R2/R3 place the RF inside a selection plateau of length `tExwd = tEx + ringdown + dead` with the
pulse delayed by `dead`. The RF effective centre is then at `dead + tEx/2`, and the plateau
midpoint at `(tEx + ringdown + dead)/2`. **These coincide if and only if
`rf_dead_time == rf_ringdown_time`,** which both references set to 100 us. When they differ, the
area before the centre and the area after it differ by `a_sel * (dead - ringdown) / 2`, and
because `k_n = -k_(n-1) + delta` that residual can alternate sign echo to echo — the odd/even
modulation an FSE is famous for.

> **Measured, it does not alternate here.** E1 found a *constant* `kz = -9.6 1/m` at all sixteen
> echoes. `k_n = -k_(n-1) + delta` solves to `delta/2 + (-1)^n (k_0 - delta/2)`, and the
> alternating term vanishes when the first interval carries the same imbalance as every later one
> — which this construction makes it do. See [`comparison_report.md`](comparison_report.md) §2.

R1's `Refocusing` removes the assumption in two independent steps: it symmetrises the plateau
about the measured effective centre
(`t_plateau = 2*max(dead + centre, duration - centre + ringdown)`), and it then solves the two
crusher amplitudes by integrating the emitted waveform. Its module docstring states the case
directly: *"The pulseq and pypulseq TSE demos get it right by accident."*

By hand, for R3: the first-interval z balance requires
`a_ex*(tExwd + dG)/2 == GSex.area/2`, which holds exactly; and the x balance requires
`gr_acq.area/2 + gr_spr.area == gr_spr.area + a_acq*dG/2 + kWidth/2`, which also holds exactly.
**The reference construction is exact — conditionally.** Both identities are built on the RF
centre lying at the plateau midpoint.

R1's own protocol uses `rf_ringdown_time = 30e-6`, i.e. precisely the case where the reference
construction breaks. This gives the fine scan a real experiment rather than a seeded fault
(experiment E1, below).

## 7. What the comparison must therefore do

Planned experiments, in order. Each produces a measurement, not a verdict.

- **E0 — common protocol.** *Done; the reconciliation itself produced findings.* Reconcile the three protocols (inventory §2) into one the three can
  all build: 256 mm FOV, 128x128, 16 echoes, TE1 fixed, `fsp_r`/`fsp_s` matched to R1's crusher
  cycles or R1's crushers matched to them. Expect to need R1's 4 ms refocusing (B1) or an explicit
  `max_b1`; that constraint is itself a finding about R2/R3.
- **E1 — the ringdown experiment.** *Done; see [`comparison_report.md`](comparison_report.md) §2.*
  The prediction held in substance and failed in one detail: R2's residual is constant rather than
  alternating, and at 30 us R2 cannot be built at all, because its derived window falls off the
  gradient raster. R1 is unchanged at every ringdown.
- **E2 — L0/L2/L3 comparison at the common protocol.** *Done; the two implementations agree.* block-count-blind: duration, ADC counts
  and times, gradient waveform on exact knots, k at every ADC sample, echo sample location.
- **E3 — L4 TSE physics.** *Done; every invariant holds on both.* refocusing-centre spacing spread, echo-vs-midpoint offset, crusher
  balance about each centre, effective TE from the table.
- **E4 — MATLAB path:** one run of `writeTSE.m` to `tse.seq`, read back with
  `pypulseq.Sequence.read`, then the same L0-L4 measurements. Deferred until E1-E3 run.

## 8. Open questions carried into the boundary decision

1. Is the reusable unit the shot (excitation + train), and does `FSE2D` then keep only the table,
   dummies and TR stacking?
2. Does the `skew` identity — `skew = (gx.area - 2*area_to_echo)/2` — belong to `CartesianLine`
   rather than to the caller? It is derivable entirely from the readout's own geometry, and every
   spin-echo consumer of `CartesianLine(prephase=False)` needs it. Moving it down is the first
   thing to test against "does this make the leaf's responsibility incoherent?".
3. Does the three-axis crusher window (`crush_min_s` -> `crush_s`) have a home, or is the
   coupling itself the kernel's whole job? `GRE2DTR` already does the same thing under the name
   "winder coupling", which is evidence that this is a *recurring kernel pattern* rather than a
   TSE one.
4. Should `Refocusing` own the readout-axis lobe pair? Its docstring explicitly declines
   ("they live in two different readout blocks, so the caller places those"). The train is the
   first caller that holds both. This is the one place where the existing leaf's stated boundary
   is under real pressure.
5. Is HASTE a configuration? R1 already answers yes by construction (`echoes=72`,
   `partial_fourier=0.625`, one shot) but this has not been compared against `writeHASTE.m`.
