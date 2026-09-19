# TSE/FSE Fine Scan — Invariant Table

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** TSE/FSE shot
**Date:** 2026-09-18
**Pass:** 3 — the physical contract, written before any candidate API

This table is the specification the comparator is built against, not a description of an API.
Provenance: [`reference_inventory.md`](reference_inventory.md). Evidence:
[`findings.md`](findings.md).

## How the "R1 achieved" column was obtained

Notebook `examples/fse_2d/01_build.ipynb` executed up to the `FSE2D` definition, then measured
through `sc.kspace` (which compiles and calls pypulseq's `calculate_kspacePP`, so ADC times are
true sample times). Protocol: the notebook's own — 256 mm FOV, 128x128, 5 mm slice, 200 Hz/px,
`echoes` in {1, 8, 16}, interleaved table, no dummy shots. Worst case over all echoes and shots.

Context measured at the same time, because several tolerances only mean something beside it:
echo spacing **10.700 ms** (train-bound at every turbo factor), TE1 **10.702 ms**, crusher window
**700 us**, chosen as the max of the three axis minima **z 700 / y 440 / x 40 us** — the slice
crusher sets it, confirming §2 of the findings. `shift_s` is 10 us here, one gradient raster.

---

## The invariants

| # | invariant | physical meaning | R1 | R2 | R3 | measurement | tolerance | R1 achieved |
|---|---|---|---|---|---|---|---|---|
| I1 | every readout crosses `kx = 0` | the echo exists: the readout's pre- and post-areas cancel somewhere inside the ADC window | dephaser `gx.area/2` + `lobe_in`/`lobe_out` pair | `agr_preph = gr_acq.area/2 + agr_spr`, mirrored `gr5`/`gr7` | same | sign change of `k_adc[0]` within each readout, the crossing located by interpolation | every readout | **all**, both implementations |
| I1-dc | how far the nearest ADC sample sits from `kx = 0` | whether the acquisition **has a DC sample**.  Reported, never judged: a design choice, and invisible to every other check | a sample at `k = 0` | **no sample at `k = 0`** -- 64 samples straddle it | same | `abs(t_sample - t_crossing) / dwell` | none -- measured, not required | seqcraft **0.000**, pypulseq **0.500** dwells |
| I2 | `ky = (line - centre) * dky` at the echo sample | the requested line, signed | `PhaseEncode(line=..)` after the 180 | `phase_areas[i_echo, i_ex]` | same | `k_adc[1]` at the echo sample vs the table | 1e-3 1/m | **2.0e-7 1/m** |
| I3 | `kz = 0` at the echo sample | slice rephased; through-slice phase not wound | `Excitation` rephaser + symmetric crushers | crusher asymmetry `AGSex*(1+fspS)` vs `AGSex*fspS` | same | `k_adc[2]` at the echo sample | 1e-3 1/m | **2.7e-8 1/m** |
| I4x / I4z | equal gradient area either side of every refocusing effective centre | the CPMG conjugation condition on the readout and selection axes.  **Not** `k = 0` at the centre -- k there is the accumulated crusher area; what must match is the area before and the area after | solved by integration in `Refocusing` | by mirrored construction | same | integrate each axis waveform from echo `n` to `t_refocusing[n+1]` and from there to echo `n+1`, on exact knots | 1e-3 1/m | **5.8e-10 / 3.2e-10 1/m** (pypulseq 3.0e-9 / 5.7e-10) |
| I4y | `ky = 0` **at** every refocusing centre | the phase axis is the one place the equal-area rule does **not** apply -- an encode is meant to move k.  What must hold is that it is rewound before the pulse and reapplied after it, or every later pulse flips the encode | blip after the 180, rewinder before the next | same | same | `k_y(echo n) + area_y(echo n -> centre)` | 1e-3 1/m | **4.2e-11 1/m** (pypulseq 4.3e-14) |
| I5 | refocusing centres uniformly spaced within a shot | any drift mixes primary and stimulated pathways with different phase histories | `t_refoc_s + n*echo_spacing_s` | `tref` per echo, identical blocks | same | `ptp(diff(t_refocusing))` per shot | 1 ns | **1.8e-15 s** |
| I6 | the echo lies at the midpoint between its two refocusing centres | primary and stimulated echoes coincide **only** there | `shift_s` exists for this | by mirrored construction | same | echo sample time vs mean of the two adjacent `t_refocusing` | 1 gradient raster (10 us) — the echo is on the ADC raster, the gradients are not | **1.95 us** |
| I7 | first interval = half the echo spacing | the same condition for the first echo, measured from the excitation's effective centre | `t_refoc_s` derived from `exc.time_to_center()` | `t_ex` fill | same | `t_refocusing[0] - t_excitation` vs `esp/2` | 1 ns | **9.2e-16 s** |
| I8 | the RF effective centre is the balance point, not the plateau midpoint | with `rf_dead_time != rf_ringdown_time` the two differ, and the residual **alternates sign echo to echo** | plateau symmetrised about `calc_rf_center`; crusher amplitudes solved by integration | assumes `dead == ringdown` (both 100 us) | same | area on the crush axis before vs after `time_to_center()`; and I3 as a function of `rf_ringdown_time` | 1e-3 1/m | holds at `ringdown = 30 us`, which is R1's protocol — **see experiment E1**, which predicts R2/R3 fail here |
| I9 | echo spacing constant across the train | contrast model and every T2 fit assume it | one `echo_spacing_s` | one `tref` | same | `diff(t_refocusing)` (I5), and echo-to-echo times | 1 ns | covered by I5 |
| I10 | effective TE = the TE of the echo carrying the centre line | the contrast the protocol claims | `te_eff_s(segments)` | `TEeff/TE1` `circshift` of `PEorder` | same | find `ky == 0` sample, report its time minus `t_excitation` | 1 gradient raster | to be measured in E3 |
| I11 | the shot fits in TR and shots are uniformly stacked | the next excitation must not land in the train | `min_tr_s` refusal + block-raster `tr_s` | `tr_delay` from `te_train`, clamped to 1 ms with a warning | same | shot start times; total duration | block raster | R1 refuses; R2/R3 warn and continue |
| I12 | hardware legality | the sequence must run | `sc.compile` verification | `seq.check_timing()` | `seq.checkTiming` | timing check + grad/slew limits | exact | R1 compiles; R2/R3 self-report |
| I13 | `z` is quiet while any ADC is open | through-slice dephasing during sampling is signal loss, and it is invisible to every k-space check | asserted in `tests/modules/test_se_notebooks.py` | `gs` is zero during `gr6` | same | sample the z waveform over each ADC window | 0 Hz/m | asserted in CI today |
| I14 | CPMG phase relation: excitation and refocusing 90 degrees apart, receiver at 0 | the fixed point of the phase map a refocusing pulse applies | `phase_deg=90` / `phase_deg=0` | `rf_ex_phase = pi/2`, `rf_ref_phase = 0` | same | RF phase offsets of the two pulse families | exact | all three agree |

Family-specific (`family_specific: true`): I4x/I4y/I4z, I5, I6, I7, I8, I10. The rest are general Cartesian
acquisition invariants that any 2D readout candidate would reuse.

---

## Notes on tolerance choices

- **1e-3 1/m for k.** `dky` here is 3.906 1/m, so this is 2.6e-4 of one k-space step — far below
  anything with an image consequence, and four orders above what R1 achieves. It is set loose
  enough that a *representation* difference between a MATLAB-written `.seq` and an in-process
  `pypulseq.Sequence` (shape compression, 1e-6 rounding in the file format) cannot fail it, and
  tight enough that a real sign or balance error cannot pass.
- **1 ns for pulse spacing.** R1 achieves 1.8e-15 s because it computes the times; a `.seq`
  round trip quantises to the raster, so 1 ns is the honest cross-implementation figure.
- **1 gradient raster for the echo-to-midpoint offset.** This is not slack: the echo instant lives
  on the ADC raster and gradient starts live on the gradient raster, so the two cannot coincide in
  general. R1's 1.95 us is well inside it. A comparator that demanded zero here would be measuring
  the raster, not the physics.
- **Exact for I12/I14.** These are either true or the sequence is wrong.

## What this table deliberately does not contain

- **Block count, block boundaries, or event decomposition.** R2/R3 spend `GS1..GS7` on what R1
  emits as one waveform; both are legal and the comparator must be blind to it (playbook §9 L2).
- **Total gradient area of the crushers.** R1 and R2/R3 disagree on it by design (findings D1/D2)
  and the disagreement is a parameter, not a violated invariant.
- **Label state.** R3 applies `autoLabel` after construction; R1 emits labels during it and omits
  `ECO` on purpose (findings D4). Not comparable like-for-like, so it is excluded until the
  boundary decision says who owns labels.


---

## Revisions

**2026-09-18, after the first comparison run.** Two rows above were wrong when written and were
corrected by measurement rather than by review; both are recorded in
[`comparison_report.md`](comparison_report.md) §1.

- **I1** said "`kx = 0` at the echo sample", which presumes a sample lands there.  `write_tse.py`
  has no DC sample, and ranking its samples by `|kx|` chooses between two values half a k-step
  apart on floating-point noise.  I1 is now the existence of the crossing, and the sampling
  difference is measured separately as I1-dc.
- **I4** was one statement across three axes.  Applied to the phase axis it "failed" for both
  implementations by ~230 1/m, because a phase encode is supposed to move k.  It is now the
  equal-area rule on x and z, and `ky = 0` at the pulse on y.
