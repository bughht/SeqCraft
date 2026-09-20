# SaturationPrep — Reference Inventory

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** working name `SaturationPrep` (`preparation/`)
**Date:** 2026-09-19
**Stage:** prospective supervised run, evidence phase. **No extraction attempted, no Module written.**

No terminal outcome was fixed before this inventory was assembled.

---

## 1. Executable, permissively licensed references

Unlike the previous candidate, this family is well covered by MIT-licensed sources.

### R1 — PyPulseq `write_epi_se_rs.py`

`imr-framework/pypulseq` @ `f2c582b`, MIT, executable Python.

```text
B0        = 2.89
sat_ppm   = -3.45
sat_freq  = sat_ppm * 1e-6 * B0 * gamma          ->  -424.50 Hz
rf_fs     = make_gauss_pulse(flip=110 deg, duration=8 ms,
                             bandwidth=|sat_freq|, freq_offset=sat_freq,
                             use='saturation')
gz_fs     = make_trapezoid('z', delay=duration(rf_fs), area=1/1e-4)   # spoil to 0.1 mm
```

### R2 — Pulseq MATLAB `writeEpiSpinEchoRS.m`

`pulseq/pulseq` @ `2fd6ab6`, MIT, documentary (ships `epi_se_rs.seq`).

Identical design to R1 — same 110°, same 8 ms, same `bandwidth=abs(sat_freq)`, same spoil to
0.1 mm — **plus one line R1 does not have**:

```matlab
rf_fs.phaseOffset = -2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs);
% compensate for the frequency-offset induced phase
```

**Independence: R1 and R2 are one witness.** Same design, one ported from the other. Their
disagreement is recorded in `findings.md` §2 as a within-design difference, not as corroboration.

### R3 — `fmrifrey/lps` `gre3d/gre3d_write_seq.m`

`fmrifrey/lps` @ `552d87a`, MIT, documentary. **Independent of the Pulseq house style** — the same
judgement recorded for this source during the GRE3D pilot.

```text
fatOffres = gamma * B0 * 3.5e-6                  # note: 3.5 ppm, and the SIGN is applied later
fatBW     = 200 Hz            (fixed, not derived from the offset)
fatPW     = 12 ms
rf_fat    = makeGaussPulse(pi/2, Duration=fatPW, freqOffset=-fatOffres,
                           apodization=0.42, timeBwProduct=fatBW*fatPW,
                           use='saturation')
then gzSpoil
```

### R4 — SeqCraft itself

`IRPrep` (`preparation/ir_prep.py`) is the structural precedent: a preparation pulse plus a
crusher, owning `rf.use='inversion'`. Not external evidence — it is the architectural pattern any
candidate here would follow or deliberately depart from.

## 2. Negative and deferred

- **`writeMPRAGE_grappa.m` is not evidence.** A grep for `saturation` hits line 171, whose comment
  reads `% fat-saturation excitation` — but the block added is `rf180` followed by a delay and a
  spoiler. That is an **inversion**, already covered by `IRPrep`. The comment is stale and the
  match is spurious. Recorded so the next scan does not re-find it.
- **OpenMRF `FAT`/`SAT` preparation library** — MIT, and the most likely second independent
  witness. **No local checkout; uninspected.** Open item, not evidence.
- **`kherz/pulseq-cest`** — MIT, saturation periods for CEST. Uninspected. A CEST saturation train
  is a different physical problem (long trains, duty cycle, B1rms) and is **out of scope** for this
  candidate; recorded so the boundary question is explicit rather than assumed.

## 3. Evidence position

```text
executable, permissive, independent designs:   TWO   (R1/R2 as one witness, R3 as another)
they AGREE on:      structure, purpose, spectral selectivity, spoil-and-do-not-rephase
they DISAGREE on:   flip angle, duration, how bandwidth is chosen, and the sign path
```

This is a materially stronger position than the previous candidate and comparable to Radial's.
