# TSE/FSE Fine Scan — Comparison Report

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-18, so that the pull request extracting `TSEShot` and `FSE2D` carries the evidence
> that justified it.  The two copies are identical today and there is nothing keeping them
> that way; see [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** TSE/FSE shot
**Date:** 2026-09-18
**Passes:** 4 — reference adapters; 5 — numerical comparison; 6 — comparator checks driven by
observed differences

Tooling: `tools/module_mining/`. Every number below is regenerable with the two scripts named,
which write `ringdown_result.json` and `comparison_result.json` beside themselves -- untracked,
because `.gitignore` excludes scan output.

---

## 1. What the comparator had to be corrected about

Both corrections came from measuring, not from reading, which is the point of the exercise.

**The echo is not "the sample where `|kx|` is smallest".** The PyPulseq demo acquires 64 samples
whose k-space centre falls exactly *between* samples 31 and 32 — it has **no DC sample**. Ranking
by `|kx|` there chooses between two values half a k-step either side of zero, and which one wins
is floating-point noise: the first run reported `kx` at the echo alternating `-1.947, +1.947,
-1.947 …`, which looks exactly like a real odd/even error and is not one. The comparator now
locates the echo as the **interpolated `kx` zero crossing** and reports the distance from the
nearest ADC sample separately, as `I1-dc`, in dwells.

**The equal-area check about a refocusing pulse does not apply to the phase axis.** Run on all
three axes it "failed" for *both* implementations by ~230 1/m — because a phase encode is meant to
leave k somewhere different from where it started. The phase axis's invariant is `ky = 0` **at**
the refocusing centre (blip after the pulse, rewinder before the next), which both references
build and which both satisfy. `invariant_table.md` I4 is now three statements, not one.

Neither correction would have been found by reading the sources side by side. Both were found by
a comparator that reports measurements rather than a verdict.

---

## 2. Experiment E1 — does the references' construction assume `dead == ringdown`?

`python tools/module_mining/candidates/tse/run_ringdown.py`

Prediction from [`findings.md`](findings.md) §6 D5: R2 places the RF inside a plateau of length
`tEx + ringdown + dead` delayed by `dead`, so its effective centre coincides with the plateau
midpoint only when the two times are equal. 16 echoes, the reference's own protocol otherwise.

| reference | `I3` kz at echo (1/m) | `I6` echo vs midpoint (s) | `I4z` area balance (1/m) | `I1-dc` (dwells) |
|---|---|---|---|---|
| seqcraft @ 100 us | 2.4e-9 | 1.95e-6 | 8.0e-11 | 0.000 |
| pypulseq @ 100 us | 9.2e-9 | **2.5e-14** | 5.7e-10 | 0.500 |
| seqcraft @ 40 us | 3.7e-9 | 1.95e-6 | 8.0e-11 | 0.000 |
| pypulseq @ 40 us | **9.60** | **3.0e-5** | **19.2** | 0.500 |
| seqcraft @ 30 us | 3.7e-9 | 1.95e-6 | 8.0e-11 | 0.000 |
| pypulseq @ 30 us | *not buildable* | | | |

Three results, one of which was not predicted:

1. **At its own assumption the reference is exact**, and on the echo-to-midpoint check it is
   *better* than SeqCraft: 2.5e-14 s against 1.95e-6 s. SeqCraft's 1.95 us is raster quantisation,
   documented in the notebook and allowed by its own test at one gradient raster. This is worth
   stating plainly because the point of the exercise is not to find the reference wanting.

2. **Break the assumption by 60 us** (ringdown 40 us against dead time 100 us) and the reference
   develops a **19.2 1/m per-pulse slice-area imbalance**, and an echo displaced **30.0 us** from
   the midpoint — exactly `(dead - ringdown)/2`, the timing half of the residual. SeqCraft is
   unmoved at 3.7e-9 and 1.95 us: its `Refocusing` symmetrises the plateau about the measured
   effective centre and then solves the crusher amplitudes by integration.

3. **At SeqCraft's actual 30 us the reference cannot be built at all.** `t_sp = 0.5*(te -
   readout_time - t_refwd)` lands at 1.725 ms, off the 10 us gradient raster, and
   `make_extended_trapezoid` refuses. A *second* hidden assumption: that the derived crusher window
   happens to be raster-legal. This one fails loudly, which is the better of the two failure modes.

### The correction to the prediction

`findings.md` said the residual "alternates sign echo to echo". **Measured, it does not** — `kz`
at the echo is a constant `-9.6 1/m` at all sixteen echoes.

The recursion is `k_n = -k_(n-1) + delta`, whose general solution is
`k_n = delta/2 + (-1)^n (k_0 - delta/2)`. The alternating term vanishes when the first interval
carries the same imbalance as every later one, which in this construction it does, so the
particular solution `delta/2` is all that survives. A constant through-slice offset of 9.6 1/m is
0.048 cycles across a 5 mm slice: a slice-profile phase and a small signal loss, **not** the
odd/even ghost the alternating form would produce. The alternating case is the one where the
excitation-side and refocusing-side imbalances differ, which this protocol does not produce.

---

## 3. Experiments E0, E2, E3 — both Python references at one protocol

`python tools/module_mining/candidates/tse/run_comparison.py`

**E0 is a result, not a setting.** The three references cannot run at each other's defaults, and
reconciling them produced facts:

- `write_tse.py` fixes the readout at 6.4 ms of sampling in its body, so SeqCraft is given
  `bandwidth_hz_px = 1/6.4e-3`;
- its 2 ms refocusing pulse peaks at **116 % of its own `max_b1`** — it warns and continues.
  SeqCraft's `Refocusing` **refuses** the same pulse at 130 % of the notebook's limit. Rather than
  raise the limit to force agreement, the comparison keeps SeqCraft's 4 ms pulse and matches the
  **echo spacing**, which is what the physics depends on;
- SeqCraft's shortest train-bound spacing at this geometry is **12.120 ms**, so that — not the
  reference's nominal 12 ms — is the common `te`.

Common protocol: 256 mm FOV, 64 x 64, 16 echoes, 5 mm slice, 6.4 ms sampling, echo spacing
12.120 ms, TR 2 s, one dummy shot each.

### L0 — inventory

| | seqcraft_fse | pypulseq_tse | delta |
|---|---|---|---|
| duration | 10.000000 s | 10.000000 s | 3.2e-14 s |
| ADC events | 64 | 64 | 0 |
| ADC samples | 4096 | 4096 | 0 |
| RF events | 85 | 85 | 0 |
| excitations / refocusing | 5 / 80 | 5 / 80 | 0 / 0 |
| pulseq blocks | **185** | **350** | *diagnostic only* |

The same acquisition in **47 % fewer blocks**, which is the `GS1…GS7` / `GR3…GR7` hand-splitting
of [`findings.md`](findings.md) §3 measured rather than described. It is reported and used by
nothing.

### L3 — acquisition semantics

| measurement | value |
|---|---|
| TE from excitation, max error | 4.6e-14 s |
| k at echo, max error, x | 4.1e-11 1/m |
| k at echo, max error, z | 1.3e-8 1/m |
| k at echo, y | 242 1/m — **the table**, reported and not judged |
| DC sample offset | seqcraft 0.0 dwells, pypulseq 0.5 dwells |
| ky coverage | both: 64 lines, 64 distinct, none missing, none duplicated |

`comparison.result = pass`.

### L4 — TSE physics, each measured on its own sequence

| invariant | seqcraft_fse | pypulseq_tse | unit |
|---|---|---|---|
| I1 readouts without a kx crossing | 0 | 0 | readouts |
| I3 kz at the echo | 3.4e-9 | 1.4e-8 | 1/m |
| I4x equal area either side of the pulse | 5.8e-10 | 3.0e-9 | 1/m |
| I4y `ky = 0` at the pulse | 4.2e-11 | 4.3e-14 | 1/m |
| I4z equal area either side of the pulse | 3.2e-10 | 5.7e-10 | 1/m |
| I5 refocusing spacing spread | 2.8e-17 | 4.2e-17 | s |
| I6 echo vs midpoint | 7.1e-15 | 4.3e-14 | s |
| I7 first interval vs half the spacing | 1.8e-15 | 8.9e-16 | s |

**At a protocol both can express, the two implementations are the same sequence.** Every physical
invariant holds on both, to between four and eleven orders inside the table's tolerances.

---

## 4. What survives as a real difference

| | difference | kind |
|---|---|---|
| D1 | slice rephasing folded into the refocusing crusher vs rephased in `Excitation` | design; both reach `kz = 0` |
| D2 | readout-axis crusher `fspR = 1.0` vs `0.0` | parameter, paid for in echo spacing |
| D3 | dummy shot plays `ky = 0` vs segment 0's real lines | design; not yet measured |
| D4 | labels: none / `autoLabel` after the fact / `LIN`+`SEG` during build | ownership question |
| **D6** | **no DC sample vs a sample at `k = 0`** | **design; new, found by the comparator** |
| D7 | 2 ms refocusing warns at 116 % of `max_b1` vs refused at 130 % | safety policy |
| D8 | a derived window off the gradient raster refuses at build time vs never arises | robustness |

D6 is worth a sentence because it is invisible to everything else: a reconstruction that assumes
sample `n/2` carries DC is half a k-step wrong against `write_tse.py` and exactly right against
SeqCraft. Neither is a bug; the difference simply has to be known.

---

## 5. Still outstanding

- **E4, the MATLAB path.** `writeTSE.m` has not been run; `matlab_seq.py` is written and the
  `use`-field question that would have invalidated it is settled (inventory §3), but no `.seq` has
  been produced. Until then R3 contributes source-level evidence only.
- **R2's provenance.** Still the fork, not `pulseq/pypulseq` (inventory §1).
- **D3**, the dummy-shot difference, is recorded and unmeasured.
