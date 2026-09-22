# TSE/FSE Fine Scan — Reference Inventory

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** TSE/FSE shot
**Date:** 2026-09-18
**Pass:** 0 — freeze provenance

Recorded before any analysis, per
[`../process/fine_scan_playbook.md`](../process/fine_scan_playbook.md) §5 Pass 0.

---

## 1. References

### R1 — SeqCraft notebook-local `FSE2D`

```yaml
repo:         bughht/SeqCraft (nested checkout at SeqCraft/)
commit:       0260f8a  (main, clean worktree at scan time)
path:         examples/fse_2d/01_build.ipynb
license:      MIT (SeqCraft LICENSE)
language:     python
dependencies: seqcraft (this checkout), pypulseq 1.5.1, numpy
entrypoint:   class FSE2D(sc.Module) defined in cell 3; protocol in cell 1
runnable:     yes — executed in CI today by tests/modules/test_se_notebooks.py
role:         primary (the implementation under evaluation)
```

Reference parameters (notebook cell 1): FOV 256x256 mm, matrix 128x128, slice 5 mm,
refocus thickness factor 1.25, excitation 2.5 ms, refocusing 4 ms, bandwidth 200 Hz/px,
slice crusher 3.0 cycles/voxel, readout crusher 0.0, spoiler 4.0 cycles on x and z,
`echoes=16`, TR 2 s, 1 dummy shot. Also builds `echoes=1`, `8`, `72`
(HASTE, `partial_fourier=0.625`).

Scanner options differ from both external references in one load-bearing way:
`rf_dead_time=100e-6` with `rf_ringdown_time=30e-6`. Both external references set the two
equal. See invariant I8.

### R2 — PyPulseq `write_tse.py`

```yaml
repo:         m-a-x-i-m-z/pypulseq-matlab-like   # NOT pulseq/pypulseq
commit:       afcd4bb
path:         examples/scripts/write_tse.py
license:      MIT ("PyPulseq Contributors")
language:     python
dependencies: pypulseq_matlab_like (the fork's own package), numpy
entrypoint:   main(...) — parameterised, returns pp.Sequence
runnable:     not yet attempted; imports `pypulseq_matlab_like`, not the installed `pypulseq`
role:         supporting
```

**Provenance caveat.** The official `pulseq/pypulseq` repository is not checked out on this
machine; only the installed wheel (1.5.1, MIT), which ships no examples. The script above is
structurally the official `write_tse.py` — same variable names, same block order, same
`pe_steps`/`pe_order` arithmetic — but it has **not** been diffed against `pulseq/pypulseq`
HEAD, and it is a fork whose stated purpose is MATLAB-compatible behaviour. Treat it as a
*second Python witness to the same design*, not as the authoritative PyPulseq example, until
the official file is obtained.

Reference parameters (defaults): FOV 256 mm, 64x64, `n_echo=16`, 1 slice, refocusing flip 180,
slice 5 mm, TE1 12 ms, TR 2 s, sampling 6.4 ms, tEx 2.5 ms, tRef 2 ms, `fsp_r=1`, `fsp_s=0.5`,
`dG=250e-6`, max_grad 32 mT/m, max_slew 130 T/m/s, rf dead = ringdown = 100 us.

### R3 — Pulseq MATLAB `writeTSE.m`

```yaml
repo:         pulseq/pulseq
commit:       2fd6ab6   (worktree dirty in one unrelated file: matlab/+mr/@Sequence/gradSpectrum.m)
file commit:  7df3671   (last change to writeTSE.m)
path:         matlab/demoSeq/writeTSE.m
license:      MIT
language:     matlab
authors:      Juergen Hennig, Maxim Zaitsev
entrypoint:   script (no function wrapper); writes tse.seq
runnable:     not yet attempted (MATLAB required)
role:         authoritative external reference
```

Reference parameters: FOV 256 mm, 256x256, `necho=16`, 1 slice, TE1 12 ms, `TEeff=100e-3`,
TR 2 s, sampling 6.4 ms, tEx 2.5 ms, tRef 2 ms, `fspR=1.0`, `fspS=0.5`, `dG=250e-6`,
max_grad 30 mT/m, max_slew 170 T/m/s, rf dead = ringdown = 100 us.

R3 also calls `seq.autoLabel('mirrorFourier',true,'sortSlices','descending')` and
`seq.testReport`, neither of which R2 does.

### R4 — deferred

`writeHASTE.m`, `write_haste.py`, `writeTSEprop.m` (all present locally) and OpenMRF TSE are
**not** part of the primary comparison. They are the second wave, after the R1/R2/R3 contract
is understood.

---

## 2. Protocol comparison

The three references do not share a protocol, so no comparison can be run at their defaults.
This table is what a common protocol has to reconcile.

| | R1 SeqCraft | R2 PyPulseq | R3 MATLAB |
|---|---|---|---|
| max_grad | 32 mT/m | 32 mT/m | 30 mT/m |
| max_slew | 130 T/m/s | 130 T/m/s | 170 T/m/s |
| rf_dead_time | 100 us | 100 us | 100 us |
| rf_ringdown_time | **30 us** | 100 us | 100 us |
| adc_dead_time | 10 us | 10 us | 10 us |
| FOV | 256 mm | 256 mm | 256 mm |
| matrix | 128 x 128 | 64 x 64 | 256 x 256 |
| echoes | 16 | 16 | 16 |
| excitation | 2.5 ms sinc, TBW 4, apod 0.5 | same | same |
| refocusing | **4 ms** sinc | 2 ms sinc | 2 ms sinc |
| refocus slab | **1.25 x** slice | 1.0 x | 1.0 x |
| readout | 200 Hz/px -> 5.0 ms sampling | 6.4 ms sampling | 6.4 ms sampling |
| echo spacing | **derived** (min over 3 axes) | **derived from TE1 = 12 ms** | derived from TE1 = 12 ms |
| slice crusher | 3.0 cycles/voxel | `fsp_s=0.5` x half-plateau area | `fspS=0.5` x half-plateau area |
| readout crusher | **0.0** | `fsp_r=1.0` x readout area | `fspR=1.0` x readout area |
| ramp time | derived from slew | fixed 250 us | fixed 250 us |
| dummy shots | 1 (replays segment 0's lines) | 1 (`i_excitation=0`, **phase area 0**) | 1 (`kex=0`, phase area 0) |
| labels | LIN per echo, SEG per shot, no ECO | none | `autoLabel` after the fact |

R1's 4 ms / 1.25x refocusing is not cosmetic: `Refocusing` refuses a 2 ms 180 above `max_b1`
(refocusing.py module docstring — a 2 ms TBW-4 sinc 180 asks 130 % of a 20 uT limit). Any
common protocol must either raise `max_b1` or keep the longer pulse; this is itself a finding
about what the references leave unchecked.

---

## 3. Executable path and its risks

```text
R1  LogicBlock -> sc.compile -> pypulseq.Sequence          (direct, in-process)
R2  main() -> pp.Sequence                                   (needs the fork's package importable)
R3  MATLAB script -> tse.seq -> pypulseq.Sequence.read      (needs MATLAB once)
```

Verified today, because it is the assumption the whole MATLAB path rests on:

- `pulseq` at 2fd6ab6 writes file format **1.5.0**, and its RF section header is
  `# id ampl. mag_id phase_id time_shape_id center delay freqPPM phasePPM freq phase use`
  (`matlab/+mr/@Sequence/write.m:74`). The `use` field is therefore **in the file**.
- pypulseq 1.5.1 `read_seq.py:133` only guesses `use` for files below 1.5.0, and warns if asked
  to guess for 1.5.0+.

This matters because `calculate_kspacePP` conjugates k at pulses whose `use` is `refocusing`.
Had `use` been lost in the round trip, every k-space invariant measured on the MATLAB reference
would have been silently wrong. It is not lost.

Open risks:

- R2 imports `pypulseq_matlab_like`; the adapter must put that package on the path rather than
  edit the script. Not yet attempted.
- R3 needs one MATLAB run to produce `tse.seq`. Deferred until the Python path works, per
  playbook §8.
- R3's `autoLabel` mutates label state after the fact; label comparison against R1 is therefore
  not like-for-like and is excluded from L3 for now.
