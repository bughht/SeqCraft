# Spiral — ownership map, boundary, helpers and reconstruction

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-20, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

Derived from the physics and from current `main`, not from PR23's class list.

---

## 1. The coupling test, applied

> If A can change while B stays physically valid without redesign, they are probably not one
> indivisible contract. If changing A forces recomputing B, that coupling is evidence they belong
> together.

| change | what must be recomputed | coupled? |
|---|---|---|
| `density` (path geometry) | the whole traversal — every gradient sample, the duration, the sample count, the ADC split | **yes, totally** |
| slew or amplitude limit | the traversal, the duration, the ADC split — but **not** the path | **one-way** |
| `angle_rad` (interleaf rotation) | nothing; the arm is rotated after it is designed | **no** |
| `variant` (out / in / in-out) | nothing, *given* `v = 0` at both ends — that is what D2 buys | **no** |
| `echoes` | fly-back presence, total duration, ADC split | yes |
| `shots` | the path (Nyquist depends on it), hence everything | **yes** |

**The path→traversal coupling is one-way**, and that is the finding. A path is a meaningful
geometric object without any scanner: it is the Nyquist condition `dθ/dk_r = 2π·FOV(k_r)/shots`
and nothing else. But no traversal is meaningful without its path. So they are **one contract with
an internal lowering boundary**, not two public modules — the caller has no use for a path they
cannot play, and no reason to swap traversal policies on a fixed path.

> **PR23's path/traversal separation is right as an internal design and wrong as a public API.**
> Verified independently: multiple traversal policies *could* realise one path, but nothing in the
> corpus or in SeqCraft's callers wants to choose one, and exposing the choice would export a
> solver's internals as a contract.

## 2. Proposed ownership

```text
SpiralReadout  (readout/, leaf)
    owns: the path, its traversal, the emitted arbitrary gradients, the ADC and its
          segmentation, the prephaser, the fly-backs, the rewinder, the join correction,
          and k_per_m / t_adc_s measured off the built events

    does NOT own: the angle schedule, the shot count as policy, the excitation,
                  the spin echo, the ordering, the reconstruction
```

Rotation is a **per-call argument**, exactly as `RadialReadout` settled it — the precedent is
directly applicable and was decided on its own evidence.

**Not proposed:** `SpiralPath`, `SpiralTraversal` (internal, §1), `SpiralInterleaf` (an interleaf
is a call, not an object — Radial's precedent), `GRESpiral2D` / `SESpiral2D` (composition until a
second consumer appears; the GRE3D precedent declined `GRE3D` on exactly this basis).

**Open, and the reason this is not yet a settled boundary:** whether `variant` and `echoes`
belong to the readout or to acquisition policy. `in-out`'s seam *is* `k=0`, which is a property of
the readout; but "read two echoes per excitation" is the same kind of statement as
`CartesianLine`'s `echoes`, which is owned by the readout there. The precedent points one way and
the mode-flag risk points the other, and four variants × echoes is where a flag explosion would
start if it starts anywhere.

## 3. Trajectory truth — PR23's strongest idea, independently assessed

Five candidate trajectories, and they are not interchangeable:

| | authoritative for |
|---|---|
| design-pass `k` | **nothing.** It is an input to a solver, not a description of the machine |
| built Module events | the module's own claims — what it *intends* to play |
| compiled sequence | what the interpreter will play, after splits, merges and any resample |
| `sc.kspace` | an independent measurement of the compiled sequence; shares no code with the module |
| sidecar consumed by reconstruction | must equal the compiled one, or the image is wrong |

**The claim survives review**: the trajectory handed to reconstruction must describe the waveform
that plays. Reporting the design trajectory is failure **C**, and it is silent.

**But exact equality is not achievable and should not be demanded.** PR23 measures 1.3–4.7e-5 1/m
between its `k_per_m` and `sc.kspace` — 3–10e-6 of a Nyquist step. My own independent measurement
of the compiler's round-trip (`compiler_claims.md`) puts the representation floor at ~1.2e-7 s/m
in m1 and ~6e-6 1/m in m0 for a split arbitrary gradient, which is the same order. **A tolerance
is needed, and it must be characterised on the floor rather than inherited from PR23's number.**

## 4. The three helpers — one consumer each

| helper | equivalent on main? | consumers today | assessment |
|---|---|---|---|
| `oblique_trapezoid` | no, but `RadialReadout._scaled` solves a *related* problem differently — it scales a canonical spoke's amplitude, area and flat area rather than building a trapezoid of area `hypot(ax, ay)` split by direction cosines | **one** (Spiral) | **do not promote yet.** Radial is a candidate second consumer only if refactored, and refactoring a shipped module to justify a helper is the wrong direction |
| `sample_quantum` | no; `dwell_quantum` is the same *idea* on the dwell rather than the count | **one** | defer |
| `segment_samples` | no. `check_event_sizes` already exists on main — added because of this very spiral — but nothing *segments*; EPI's `segment` is a shot index, an unrelated meaning | **one** | defer |

PR #27's EPI timing work did **not** create a second consumer for any of the three: it added
`shot_delay_s` / `shot_offset_s`, which is interleaf timing, not ADC segmentation.

> Two consumers justify examining a shared helper. One does not justify promotion. All three
> should live inside the Spiral module until something else needs them.

## 5. Reconstruction — audit it separately

`examples/noncartesian_recon.py` (470 lines) covers trajectory convention, a NUFFT operator, DCF,
a CG solve, coil sensitivities, off-resonance time segmentation and preconditioning.

**Its placement decision is correct and should be kept**: `examples/`, not `src/`, because a
reconstruction in the package would put an ESPIRiT dependency on the compile path. That is
PR23's own argument (D18) and it holds.

**The deferred Radial Layer-3 question now has its second consumer** — this is the trigger
recorded in `radial/evidence_state.yaml`. But the trigger says *inspect whether a minimal shared
contract exists*, not *promote one*. Which pieces are genuinely trajectory-agnostic is
**unassessed**: I have read the file list, not the file. That assessment is its own task and
should not be folded into the readout boundary, because reconstruction needs must not dictate the
Module API.

**A valid outcome remains that the utility stays example-only.**

---

# 6. `variant` and `echoes` — settled by ownership, not by combination count

Applying the test the brief sets: does the parameter change the physical solve the readout owns?

## `echoes` — **belongs to the readout**

Not "another echo repeats a completed readout". Under an arm that begins and ends at rest, a
second echo either **flies back** (one-arm variants) or **reverses** (two-arm variants), and
either way it changes:

```text
gradient continuity across the join       ADC regions -- one ADC may span joined arms
where k = 0 falls in each echo            echo spacing feasibility under the slew limit
total waveform duration
```

Every one of those is a correctness condition of the emitted waveform. The precedent agrees:
`CartesianLine.echoes` owns exactly this, and its monopolar/bipolar trade is the same trade at a
Cartesian scale. **Readout, not acquisition policy.**

## `variant` — **belongs to the readout, as a declared mode set**

All four share one arm, one traversal and one set of joins; the variant selects direction and
count, not a different machine. But it is not free either — each variant moves k = 0:

```text
out      k = 0 is the FIRST sample
in       k = 0 is the LAST sample
in-out   k = 0 is the SEAM, and must fall INSIDE a sampling window, not in an ADC guard
```

That last line is a correctness condition that differs per variant, which is what makes this a
**mode set** rather than a parameter — so **rule A fires** and a mode-contract table is required
before implementation, with the differential stated explicitly.

**The trap rule A exists for is live here.** "in is out reversed" is structurally the same
reasoning as "mode B is mode A minus Gz", which is what produced the GRE3D RF defect. And
`reference_review.md` §3 shows there is **no reference evidence** for in, in-out or out-in — the
independent witness ends at full gradient and is spiral-out only.

**That does not block them.** GRE3D set the precedent: its slab-selective path had no executable
reference either and was carried by signed-moment arithmetic that depends on no implementation.
The same applies — reversing a piecewise-linear waveform and integrating it is checkable
analytically. But the claim scope must say so:

```text
spiral-out            reference-supported, comparable against writeSpiral.m
in / in-out / out-in  ANALYTIC ONLY -- no independent witness exists
```

## Conclusion

One public module, `SpiralReadout`, owning `variant` and `echoes`, with a mode contract for
`variant` and a **two-tier acceptance claim**. A split is not justified: the variants are not
different physical machines with unrelated conditional behaviour, they are one arm traversed
differently — which is exactly what `v = 0` at both ends was chosen to buy.

# 7. Reconstruction audit — `examples/noncartesian_recon.py`

470 lines, read on the PR branch. Its `SpiralReadout` is a **data container**, and this is the
finding:

```python
k_per_m: (shots, 2, samples)     t_adc_s: (samples,)
fov_m: float                     matrix: int          te_s: tuple[float, ...]
```

Nothing in that shape is spiral-specific. **A radial acquisition fits it unchanged** — spokes are
shots, and `te_s` is where k = 0 falls.

| concern | trajectory-agnostic? |
|---|---|
| the trajectory + timing container | **yes** — only its *name* is spiral |
| `encoding_operator` (sigpy NUFFT), `_coord` normalisation | **yes** — takes coordinates |
| `dcf(method='pipe')` — Pipe–Menon | **yes**, and it needs no analytic density |
| `dcf(method='radial')`, weights by `\|k\|` | **yes**, and it is literally named for the other consumer |
| CG solve, row-sum preconditioner, spectral norm | **yes** |
| coil sensitivities | **yes** |
| off-resonance time segmentation (`segments_for`, `echo_relative_times`) | **yes** — needs `t_adc_s` and `te_s`, which any trajectory has |
| `delayed()` gradient-delay model | **partly** — shifting sample times is generic; "on a spiral this is a rotation *and* a radial shift" is the spiral reading of it |
| `dcf` analytic Jacobian `\|dk/dt\| / FOV(\|k\|)` | **no** — needs `density`; genuinely spiral-specific |

**So a shared contract plausibly exists and it is small**: the container, plus operators that take
coordinates. One method of one function is spiral-specific.

**This is an inspection result, not a promotion.** The Radial trigger says *inspect whether a
minimal shared contract exists*, and the answer is "yes, and it is smaller than the file". What to
do about it is a separate decision, and **keeping the utility example-only remains valid** — the
placement argument (an ESPIRiT dependency must not reach the compile path) is untouched by this.

---

# 8. Implementation findings

Written during the extraction, from what went wrong rather than from the plan.

## 8.1 Three silent failures the implementation hit

Each compiled, looked right, and was wrong. Two of them are PR23's failures arriving on schedule.

| | what looked reasonable | what it cost | how it was found |
|---|---|---|---|
| **rewinder sized against the design end-point** | the design says the arm ends at `k_max`, so rewind `-k_max` | **0.015 `dk`** left behind every TR — PR23's failure **E**, which it measured at 0.024 `dk` | total m0 of the assembled block |
| **limits measured on raster-spaced differences** | the waveform is on the raster, so difference it by the raster | the emitted first and last intervals are **half** a raster, so real slew is **twice** what that reports; the design passed its own check and `make_arbitrary_grad` refused it | the factory refusing a waveform the module had certified |
| **acquisition outlasting its gradient** | size the ADC by duration over dwell | each segment spends two dead times, so the ADC ran 8 µs past the arm; the compiler holds the block open and pads the waveform with area nobody designed | compiler m0 contract |

The second is the sharpest, and it generalises past spirals: **a module that measures its own
output on a more convenient lattice than the one it emits will certify something it cannot build.**
`limits()` now reads the emitted knots, which is what §9 of the acceptance plan asked for and is
why it is a measurement rather than a restatement.

## 8.2 A compiler defect that was not one

A multi-segment readout lost 0.02–0.08 1/m of m0, with the same signature as PR23's m0 claim —
tree sum ≈ 0, compiled visibly different, "a split or a merge lost area".

**The escalation rule caught it.** Removing the module and splitting a bare sign-changing
arbitrary gradient across ADC-driven block boundaries is exact to **1e-13**, across 1, 2, 4, 8 and
16 events. So the compiler splits correctly and the loss is this module's.

Recorded because the near-miss is the point: the same evidence that would have "confirmed" PR23's
m0 claim was available, and the reproducer is what stopped it being claimed.

## 8.3 The defect, minimized and fixed

**Minimal failing case:** `matrix=112, shots=1, variant='out'`, one protocol wide.

```text
m111  seg 2  n 12476  dur 49950 us  acq_end 49950 us   ok
m112  seg 2  n 13112  dur 52490 us  acq_end 52500 us   FAIL  m0 off by -0.0446 1/m
m113  seg 2  n 13284  dur 53190 us  acq_end 53190 us   ok
```

Both neighbours pass. Segment count is not the discriminator -- `m96` and `m104` pass with two
segments. **The discriminator is `acq_end > dur`**, and `m112` was the only case in the range
where the acquisition ran past its gradient.

**Where m0 first diverges**, traced through the pipeline:

| stage | m0 |
|---|---|
| designed path / traversal | correct |
| emitted gradient knots | `254.47896` 1/m against a design end of `254.54545` -- the expected representation difference, already absorbed by the rewinder |
| ADC segmentation / acquisition sizing | **the fault is introduced here, and is invisible in m0** |
| assembled `LogicBlock` | **exactly `0.0`** on both axes |
| compiled Pulseq blocks | **`-0.0446`** |

So the tree is right and the compile is wrong -- but not because the compiler is wrong. Block
durations are `26250 + 26250 + 400` us against a gradient of `52490` us: the two ADC blocks total
**10 us more** than the arm, the compiler holds the last block open, and the waveform is padded
with area nobody designed.

**Root cause.** Each ADC event spends a lead delay and a trailing dead time **and** its span is
rounded up to the gradient raster. So the overhead per segment is up to `lead + trail + raster`,
while the sample budget subtracted only `lead + trail`. One raster of under-estimate per segment
is enough.

**The fix**, and the smallest one that is actually justified: size the budget **against the
placement** rather than against an estimate of it. `budget()` now computes where the events would
land and shrinks the sample count by one divisor at a time until the last one finishes inside the
gradient. No contract changed, no refusal added, no compiler touched.

**Result:** 96 protocols -- four matrices by three shot counts by four variants -- all compile,
with `acquisition_end_s <= duration_s` in every one. The strict xfails are gone, replaced by a
regression test that asserts the invariant the defect violated.

**Why it hid.** The under-estimate is silent unless the rounding happens to push past the
gradient, which needs the sample count, the dwell, the dead times and the raster to line up
against a duration that is itself an output of the traversal. `m112` was that alignment; `m111`
and `m113` were not.

## 8.4 What PR23 got right, confirmed by re-deriving it

The path/traversal split, `v = 0` at both ends, measuring `k` off built events, correcting the
join on the assembled waveform, and keeping reconstruction out of `src/` — all five were reached
independently here, and three of them only after the implementation failed the way PR23 said it
would. **That is the calibration result**: the ideas survive independent re-derivation, and the
parts that did not survive are the m0 compiler change and the helper promotions.

---

# 9. A sequence-side defect found by the reconstruction work

Found while building the Phase-A reconstruction contract, and **isolated rather than fixed there**
— the rule is that reconstruction work does not change sequence code.

## The defect

```text
variant='in'   reports ZERO origin crossings, at every matrix tried
```

| matrix | variant | crossings | min \|k\| reached | half-step |
|---|---|---|---|---|
| 16 | `in` | **0** | 3.245 | 2.083 |
| 32 | `in` | **0** | 3.259 | 2.083 |
| 64 | `in` | **0** | 2.321 | 2.083 |
| 128 | `in` | **0** | 3.643 | 2.083 |

`out` and `in-out` are unaffected. `out-in` reports **2** crossings only at `matrix=128` and
**1** below it, which is the same cause at the tail.

## Why

`_crossings` searches the **acquired samples** for one within half a Nyquist step of the origin.
For `in`, the trajectory reaches the origin at the *end of the arm* — but the acquisition is sized
to finish **inside** the gradient, so sampling stops while the arm is still braking toward the
centre, several steps short. No sample is near enough, so the module truthfully reports none.

Truthfully, and uselessly: `time_to_echo(0)` then raises, and a caller has no echo for a variant
whose mode contract says it has one, at the last sample.

## Why it was not caught

Every existing test asks *"is each reported crossing near the origin?"* — which is vacuously true
of an empty tuple — or exercises `out` / `in-out`, which have crossings. Nothing asserted that a
variant **has** the number of crossings its mode contract promises. The contract says so in prose
and no test read it.

## Disposition

```text
status          OPEN, isolated, not fixed in the reconstruction phase
scope           variant='in' entirely; variant='out-in' loses its second crossing below m=128
workaround      none needed for Phase A -- the reconstruction tests use 'in-out', whose echo is
                mid-acquisition and makes the same point about not inferring an echo
revisit trigger BEFORE the Phase A PR merges.  It contradicts the mode contract, which is a
                v0 deliverable, so it should not sit open behind a merged post-v0 PR.
```

The fix is likely small — the acquisition should either extend to cover the arm's approach to the
origin, or the crossing search should locate the trajectory's closest approach rather than require
a sample within half a step — but choosing between those is a contract question about what
`origin_crossing_samples` promises, not a debugging question, and it belongs in its own change.
