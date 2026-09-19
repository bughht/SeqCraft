# RadialReadout Fine Scan — Invariant Table

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** `RadialReadout`
**Date:** 2026-09-18
**Pass:** 3 — the physical contract, written before any API

The specification the comparator is built against. Provenance:
[`reference_inventory.md`](reference_inventory.md). Evidence: [`findings.md`](findings.md).

## How the measured column was obtained

`python tools/module_mining/candidates/radial/run_archaeology.py` — the official
`write_radial_gre.py` executed unmodified at 8 spokes of 64 samples, FOV 260 mm, 2 dummies, then
measured through `calculate_kspacePP`. Worst case over every acquired spoke. Δk = 3.8462 1/m and
`k_max` = 123.08 1/m at this protocol.

Spoke geometry is measured, not assumed: the direction is the **principal axis of the sampled
points**, so a centre-out spoke — which does not start at `-k_max` — reports the same quantity as
a full one.

---

| # | invariant | physical meaning | how measured | tolerance | R1 measured |
|---|---|---|---|---|---|
| R1 | a sample lands on `k = 0` | every spoke shares the centre of k-space, and it is the most densely sampled point in the acquisition. A spoke whose samples straddle it puts each spoke's centre in a different place | `min |k|` over the spoke's ADC samples | 1e-3 1/m | **2.2e-12** |
| R2 | the samples lie on a straight line through the origin | it is a *spoke*. A curved or offset one still produces a plausible image through a gridding reconstruction | max distance of the samples from the fitted line | 1e-3 1/m | **1.9e-13** |
| R3 | the sample spacing is uniform along the spoke | non-uniform Δk is a density-compensation error, not an obvious artefact | `ptp(diff(k along the spoke))` | 1e-3 1/m | **5.6e-12** |
| R4 | every spoke has the same extent and spacing | otherwise the sampling density varies with angle and the image is shaded | `ptp` of extent and of Δk across spokes | 1e-3 1/m | **2.7e-12** |
| R5 | the spoke points where it was asked to | a wrong angle is a rotated image, or with a wrong *increment*, streaks | fitted angle against the requested one, modulo 180 deg | 1e-6 deg | **2.8e-14** |
| R6 | the spoke stays in plane | `kz` at the samples; through-slice winding is signal loss, invisible to every in-plane check | `max |kz|` | 1e-3 1/m | **1.5e-11** |
| R7 | **the spoke at `phi` is the spoke at 0, rotated** | the defining property of a radial acquisition, and the one no single spoke can exhibit | rotate spoke 0's samples by `phi`, compare to the measured spoke | 1e-3 1/m | **3.3e-12** |
| R8 | Δk follows `1/FOV` | the resolution the protocol asked for | Δk against `1/fov` (and against the asymmetry coupling where present) | 1e-9 relative | *sweep, not yet run* |
| R9 | hardware legality | the sequence must run | timing check, gradient and slew limits | exact | R1 self-reports pass |

R7 is the metamorphic test the next-phase plan asked for, and it already passes **on the
reference**. That is deliberate: a metamorphic test written after a candidate exists tends to
encode the candidate's own behaviour, and one validated against an independent implementation
first does not.

## The centre sample is not the middle sample

Measured on R1: 64 samples, centre at index **32**, k running −32Δk … **+31**Δk. A full spoke
with an even sample count is asymmetric by one sample. Under partial echo it moves further —
`write_ute.py` walks it from 64/128 to 0/128 as `readout_asymmetry` goes 0 → 1 (see
[`findings.md`](findings.md) §6 Q1).

Any candidate must therefore **report** its centre sample rather than let a caller assume
`n/2`, exactly as `CartesianLine.echo_sample` does.

## What this table deliberately does not contain

- **Block structure and how the rotation was achieved.** Four of the five references call a
  pulseq `rotate` helper and one scales amplitudes onto two axes; the trajectories are what must
  match.
- **The angle schedule.** Equal-increment, golden-angle and randomised orderings are policy
  (findings §6 Q4). R5 checks that a spoke points where it was *asked* to, not that the sequence
  of requests was wise.
- **The prephaser's formula.** Three references derive it three ways and one of them measures it
  instead; R1 is the statement of what they are all trying to achieve, and it is the only thing
  worth asserting.

## Not yet measured

- **R8**, which needs a sweep over FOV, matrix and bandwidth.
- Anything on a MATLAB-written `.seq`. Under the Python-only constraint these are evidence rather
  than a gate, and none has been run for this candidate.
- Centre-out beyond `write_ute.py`'s own protocol: the asymmetry sweep in findings §6 was measured
  at one FOV and one matrix.
