# TSE/FSE Fine Scan — Sweep, Metamorphic Tests and Status

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-18, so that the pull request extracting `TSEShot` and `FSE2D` carries the evidence
> that justified it.  The two copies are identical today and there is nothing keeping them
> that way; see [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** TSE/FSE shot
**Date:** 2026-09-18
**Passes:** 8 — parameter sweep; 9 — metamorphic tests; 10 — status

Raw output: `tools/module_mining/candidates/tse/{sweep,metamorphic}_result.json`.

---

## 1. Parameter sweep

`python tools/module_mining/candidates/tse/run_sweep.py` — the SeqCraft implementation at eight
cases chosen so that each can fail on its own account, rather than a Cartesian product.

Worst value per invariant, over every shot and echo of each case:

| case | I2 ky (1/m) | I3 kz (1/m) | I4x (1/m) | I4y (1/m) | I4z (1/m) | I5 (s) | I6 (s) | I7 (s) | esp (ms) | bound |
|---|---|---|---|---|---|---|---|---|---|---|
| echoes=1 | 2.0e-7 | 2.7e-8 | — | — | — | — | — | — | 10.700 | train |
| echoes=8 | 2.0e-7 | 6.8e-9 | 3.4e-9 | 2.0e-7 | 4.4e-10 | 1.8e-15 | 1.95e-6 | 8.9e-16 | 10.700 | train |
| echoes=16 | 2.0e-7 | 3.7e-9 | 1.1e-9 | 2.0e-7 | 8.0e-11 | 1.8e-15 | 1.95e-6 | 8.9e-16 | 10.700 | train |
| matrix=64 | 9.5e-8 | 2.6e-9 | 2.4e-10 | 9.5e-8 | 4.0e-11 | 2.8e-17 | 8.25e-6 | 8.9e-16 | 10.700 | train |
| linear | 2.0e-7 | 3.7e-9 | 1.1e-9 | 2.0e-7 | 8.0e-11 | 1.8e-15 | 1.95e-6 | 8.9e-16 | 10.700 | train |
| centric | 2.0e-10 | 3.7e-9 | 1.1e-9 | 8.4e-11 | 8.0e-11 | 1.8e-15 | 1.95e-6 | 8.9e-16 | 10.700 | train |
| esp +2 ms | 7.4e-8 | 3.7e-9 | 2.1e-9 | 7.4e-8 | 1.4e-10 | 4.4e-16 | 1.95e-6 | 8.9e-16 | 12.720 | train |
| haste | 2.5e-8 | 1.9e-10 | 2.5e-10 | 2.5e-8 | 2.3e-11 | 1.1e-16 | 5.15e-6 | 8.7e-18 | 10.620 | **adc clearance** |

**8 of 8 cases pass every invariant.** Tolerances are the table's: 1e-3 1/m for k, 1 ns for pulse
spacing, one gradient raster (10 us) for echo placement.

Read with care where it matters:

- **`echoes=1` reports zero for I4-I7 because those invariants are vacuous there** — one
  refocusing pulse has no spacing, no midpoint and no interval to the next. Shown as `—` above
  rather than as a passing zero.
- **The largest echo-to-midpoint offset is 8.25 us, at `matrix=64`**, still inside one gradient
  raster. It grows with the dwell, which is the expected behaviour: the echo instant lives on the
  ADC raster while gradient starts live on the gradient raster.
- **HASTE is `adc clearance`-bound**, not train-bound — a different one of the three minima wins,
  which is exactly the case where the echo would drift off the midpoint without the shift the
  implementation applies. It does not drift: 5.15 us.
- The **requested** spacing of 12.700 ms is delivered as 12.720 ms, the next gradient raster.

## 2. Metamorphic tests

`python tools/module_mining/candidates/tse/run_metamorphic.py` — three truths that need no second
implementation, which matters because the two references agree so thoroughly that agreement has
stopped being informative.

| | claim | result |
|---|---|---|
| M1 | turbo factor 1 is a spin echo | **cited, not re-run** — SeqCraft CI already pins `FSE2D(echoes=1)` against `examples/se_2d/01` event for event, in `tests/modules/test_se_notebooks.py` |
| M2 | reordering the table leaves the readout waveform untouched | **pass** — `gx` differs by **exactly 0 Hz/m** between interleaved and linear, while `gy` differs by 6.6e5 Hz/m. Both halves matter: a table that changed nothing would mean it was never applied |
| M3 | the claimed effective TE equals the measured time to the `ky = 0` sample | **pass** — claimed 96.301950 ms, measured 96.301950 ms |
| M4 | lengthening the train does not move the echoes already there | **pass** — the first eight echo times are identical between `echoes=8` and `echoes=16`, drift **exactly 0 s** |

M3 is the one worth keeping. Every other check compares a construction against another
construction; M3 compares what the module **says** against what its sequence **does**, which is
the failure no amount of k-space agreement can catch — a module can compute a perfectly correct
sequence and report the wrong TE for it, and only the report reaches the protocol.

## 3. Status

```text
YELLOW
```

Not RED: nothing failed. Not GREEN: the playbook's GREEN requires that *the candidate boundary is
clear*, and it is not, because no extraction has been attempted.

What is settled, and what a reviewer should not need to re-derive:

- **The physics is proven.** At a protocol both can express, the notebook-local `FSE2D` and the
  PyPulseq TSE demo are the same sequence — identical L0 inventory, TE agreeing to 4.6e-14 s, k at
  the echo agreeing to 4e-11 (x) and 1.3e-8 (z) 1/m — in 47 % fewer pulseq blocks.
- **The SeqCraft implementation is strictly more robust in three specific ways**, each measured:
  it does not assume `rf_dead_time == rf_ringdown_time` (E1), it does not assume its derived
  crusher window lands on the gradient raster (E1, where the reference cannot be built at all), and
  it refuses rather than warns at a 180 over `max_b1`.
- **The reusable unit is one complete shot**, excitation included: the first interval's length,
  the dephaser's area *and* placement, and the midpoint shift all couple the excitation's
  effective centre to the first refocusing centre, in all three references.
- **Eight sweep cases and three metamorphic tests pass**, including a single-shot 72-echo HASTE
  that is bound by a different minimum than the others.

What is not settled, and is a design judgement rather than a measurement:

1. **Where the shot ends and the scan begins.** The coarse scan's hypothesis — `kernel/TSEShot`
   plus `imaging/FSE2D` — is consistent with everything measured, but so is promoting `FSE2D`
   whole. Nothing in the numbers chooses.
2. **Whether the readout skew identity belongs in `CartesianLine`.**
   `skew = (gx.area - 2*area_to_echo)/2` is derivable entirely from the readout's own geometry,
   and every spin-echo consumer of `CartesianLine(prephase=False)` needs it. This is the strongest
   candidate for moving logic *down* rather than extracting it sideways.
3. **Whether `Refocusing` should own the readout-axis lobe pair** it currently declines on the
   grounds that "they live in two different readout blocks". The train is the first caller that
   holds both.
4. **Whether the three-axis crusher-window coupling has a home.** `GRE2DTR` already does the same
   thing under the name "winder coupling", which is evidence that this is a recurring *kernel
   pattern* rather than a TSE one — and therefore evidence against giving it a TSE-specific name.
5. **Who owns labels** (D4), which the three references answer three different ways.

## 4. Recommended next actions

- **Ask for a human decision on questions 1-4 before writing `TSEShot`.** They are what the
  playbook reserves for MRI review, and the fine scan has done its job by isolating them.
- **Run `writeTSE.m` once** to close E4. The adapter exists and the `use`-field risk is settled;
  what is missing is one MATLAB session.
- **Obtain the official `pulseq/pypulseq` `write_tse.py`** and diff it against the fork, so R2's
  provenance caveat can be dropped.
- **Measure D3**, the dummy-shot difference, which is the one recorded disagreement with no
  number against it.

## 5. What the pilot established about the workflow itself

The stated goal was to prove the *workflow*, not to ship a class. Two things are worth carrying
into the next candidate:

- **The comparator was wrong twice, and measuring is what corrected it** — once about where the
  echo is, once about which axes the equal-area rule applies to. A comparator that returned
  `True/False` would have reported a false failure on both occasions and a plausible story would
  have been written around it.
- **Reconciling the protocols (E0) produced findings before any comparison ran**: the B1 overrun,
  the fixed sampling time, and the raster-legality assumption all surfaced while trying to make
  two references describe the same experiment.
