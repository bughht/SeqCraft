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

# 9. A sequence-side defect found by the reconstruction work — found, then fixed

Found while building the Phase-A reconstruction contract and **isolated rather than fixed there**
— the rule is that reconstruction work does not change sequence code. Fixed in its own change
afterwards, before the Phase-A PR merged.

## The defect, as first seen

```text
variant='in'   reports ZERO origin crossings
```

| fov | matrix | shots | variant | crossings | min \|k\| | half-step |
|---|---|---|---|---|---|---|
| 240 | 16 | 1 | `in` | **0** | 3.245 | 2.083 |
| 240 | 32 | 1 | `in` | **0** | 3.259 | 2.083 |
| 240 | 64 | 1 | `in` | **0** | 2.321 | 2.083 |
| 240 | 128 | 1 | `in` | **0** | 3.643 | 2.083 |
| 220 | 128 | 4 | `in` | **0** | 2.281 | 2.273 |
| 220 | 128 | 4 | `out-in` | **1** of 2 | 2.550 | 2.273 |

`time_to_echo(0)` then raises, and a caller has no echo for a variant whose mode contract says it
has one, at the last sample.

## Why

`in` reaches the origin at the *end* of the arm, where the gradient is braking to rest, so the
last `dk` of the trajectory lives in the last few microseconds. `_plan_adc.budget` sized the
acquisition by **shrinking a deliberately conservative estimate and never growing it back**:

```python
samples = int((duration - count * (lead + trail + raster)) / dwell)   # conservative
samples -= samples % divisor
while placed_end(...) > duration:                                     # shrink only
    samples -= divisor
```

`lead + trail + raster` over-reserves, and the `divisor` floor rounds further down, so the
acquisition finished 32–36 µs before the gradient where the raster admitted 12–20. Those 20 µs
are the last `dk`.

**This is the measurement that names it.** The two ends of the arm should be symmetric — the same
reservation at each — and they were not. Measured as a fraction of the half-step, which is the
quantity that actually decides the outcome, over 24 protocols (3 FOVs × 4 matrices × 2 shot
counts):

```text
                        |k| at the sample nearest the origin, / half-step
out   (first sample)                   0.130 -- 0.299
in    (last sample, before)            0.453 -- 1.498      -- over 1.0 is a missed crossing,
                                                              and 6 of the 24 were
in    (last sample, after)             0.204 -- 0.901
```

The fix restores the symmetry rather than adjusting a threshold: `budget` now also grows while
`placed_end` still fits. `placed_end` is non-decreasing in the sample count, so the loop
terminates on the largest acquisition the raster admits.

It also **acquires more data**, which was the other cost of the bug and was invisible: every
variant was discarding up to five samples per segment for no reason.

**The margin is smaller at the tail than at the head, and that is not fixed.** `out` sits at a
third of the half-step at worst; `in` sits at nine tenths. The asymmetry is structural — the ADC
reserves a lead delay at the head and a dead time *plus* a raster rounding at the tail — and it
means a protocol harsher than any tested here (higher slew, larger FOV, coarser raster) could
still put `in`'s last sample outside the half-step. That is recorded as a margin to watch, not
closed, and not papered over with a runtime refusal on evidence this thin.

## Why it was not caught — and one correction to the first reading

The first reading of this was that *"nothing asserted that a variant has the number of crossings
its mode contract promises."* **That was wrong**, and worth recording because the wrong reading
would have produced the wrong fix. `test_origin_crossing_count` asserts exactly that, and has
since the module was written.

It ran at **one** protocol — `matrix=64, shots=4, fov=220` — which sits on the lucky side of a
knife edge. At that protocol the last sample lands at 1.53 against a half-step of 2.273 and the
test passes. One matrix step away it lands at 2.281 and the test would have failed.

So the gap was not a missing assertion but an **unstressed margin**: the quantity the assertion
depends on is not one the default protocol varies. `test_origin_crossing_count` is now
parametrized over six protocols spanning matrix 16–128, 1 and 4 shots, and FOV 220 and 400 mm —
FOV because it halves `dk` and so halves the margin the last sample must land inside. Two of the
twenty-four combinations fail without the fix.

A sweep of 96 combinations (3 FOVs × 4 matrices × 2 shot counts × 4 variants) now agrees with the
mode contract in every case; before the fix it did not.

## Disposition

```text
status          RESOLVED
fix             src/seqcraft/modules/readout/spiral_readout.py -- _plan_adc.budget grows as
                well as shrinks
test            tests/modules/test_spiral_readout.py::test_origin_crossing_count, parametrized
                over six protocols; fails on two of them without the fix
scope           no contract change.  origin_crossing_samples still promises a sample within half
                a Nyquist step of the origin, and now the acquisition actually reaches one
margin          variant='in' lands at 0.20--0.90 of the half-step against 'out' at 0.13--0.30.
                Structural: the tail reserves a dead time and a raster rounding, the head only a
                lead delay
revisit trigger if any protocol is found where 'in' or 'out-in' reports fewer crossings than its
                mode contract promises.  The parametrized test is the detector; widen its grid
                before widening the module's claimed protocol range
```

## What this says about the process

Two things, neither of them about spirals.

**A margin that no test varies is not tested.** The assertion was right, the fixture was one
point, and a one-point fixture cannot find a knife edge. This is the same shape as rule E —
measure on the lattice you actually emit — applied to the *parameter* space rather than the time
lattice.

**Reconstruction found a sequence defect that sequence tests did not.** Phase A's stated purpose
was to validate the reconstruction adapter; it also exercised the modules at protocols their own
suites did not, because a reconstruction has to have data at k=0 and therefore cares about a
quantity the sequence tests only asserted. That is an argument for Layer 3 as evidence, not only
as demonstration — and it is recorded here rather than promoted to a rule, because one instance
is one instance.

---

# 10. A second sequence-side defect, found by composing a spin echo

Found while building `examples/se_spiral_2d/01_build.ipynb`, which is the first thing in the
repository to place a *refocused* echo against a reported origin crossing. Fixed in its own change
before the notebook was written, rather than worked around in it.

## The defect

`SpiralReadout` reported **two different clocks under one docstring.**

```text
time_to_echo(0) docstring   "the time from the readout block's start to origin crossing index"
what it returned            the time from the start of the ARM
```

Those differ by the prephaser, and only for the variants that have one:

| variant | prephaser | `time_to_echo(0)` reported | measured on the compiled block |
|---|---|---|---|
| `out` | 0 µs | 12.0 µs | 12.0 µs |
| `out-in` | 0 µs | 12.0 µs | 12.0 µs |
| `in` | 300 µs | 2744.0 µs | **3044.0 µs** |
| `in-out` | 300 µs | 2772.0 µs | **3072.0 µs** |

`sample_times_s()` carried the same error, so the whole reported trajectory was 300 µs early on
its own block.

## Why it is the dangerous kind

Everything about the sequence is *correct*. The trajectory is right, the k positions are right to
1e-12, the file compiles, the arm plays exactly as designed. The only thing wrong is the **label
on the clock** — and the consequence is that a kernel composing `TE` from
`start + readout.time_to_echo(0) - excitation.time_to_center()` puts the echo 300 µs off, in a
sequence in which nothing looks wrong.

It is also the precise failure `time_to_echo`'s own docstring warns about, one level down: *"a
spin echo placed against the wrong crossing is a legal sequence whose echo time is wrong"*. The
module warned about the caller choosing the wrong crossing and then mis-stated when its own
crossings were.

## Why nothing caught it

Three things had to line up, and did.

```text
gre_spiral_2d/01 uses variant='out'      the only shipped consumer, and it has no prephaser
the k checks compare POSITIONS           sc.kspace agreement to 1e-12 says nothing about instants
RadialReadout was already right          so the convention existed and nothing compared them
```

`RadialReadout.time_to_center()` has always included its own 320 µs prephaser and is verified
against the compiled block. The two sibling readouts documented the same contract and implemented
different ones, and no test looked at both.

## The fix

The internal clock stays on the arm — `_knots` starts there, and `duration_s` and
`acquisition_end_s` are statements *about* the arm's gradient, so "the acquisition finishes inside
its gradient" stays a comparison of two arm-clock quantities. What moves is the **reporting**:

```text
_arm_sample_times()      new, private, the integration clock
sample_times_s()         public, block-relative -- the arm clock plus the prephaser
origin_crossing_times    block-relative, so time_to_echo and echo_spacing_s follow
prephaser_duration_s     new, public: the bridge, and what build(prephase=False) needs subtracted
```

A first attempt also moved `acquisition_end_s` onto the block clock. Two tests failed immediately,
and they were right to: that quantity is compared against `duration_s`, so moving one side of a
comparison broke the invariant it encodes. The two clocks are now named in the docstrings rather
than unified, because they measure different things.

## Disposition

```text
status          RESOLVED
fix             src/seqcraft/modules/readout/spiral_readout.py
tests           test_reported_times_are_on_the_block_the_module_builds (all four variants,
                against sc.kspace) and test_the_prephaser_is_the_offset_between_the_two_clocks.
                Both fail on the old code; the first fails only on 'in' and 'in-out'
verified        24 configurations -- four variants, three shot counts, two matrices -- agree
                with the compiled sequence on crossing time, sample times and k
scope           no change to any emitted waveform.  The 'out' variant, and therefore every
                shipped .seq before this change, is bit-identical
```

## What this says about the process

**Two modules documenting the same contract is not evidence that they implement it.** The
ownership rule put `time_to_center` on `RadialReadout` and `time_to_echo` on `SpiralReadout` for
the same reason, and the identical wording made it look settled. Nothing compared the two, and a
convention held by prose in two files is held by nothing.

**A check on positions is not a check on instants.** The Spiral candidate's Layer-1 evidence is
unusually strong on *where* the samples are — agreement with an independent measurement to
1e-12 — and that strength is exactly what made the timing gap invisible. Rule E says to measure on
the lattice actually emitted; this is the same rule in the time axis, and the candidate record's
acceptance criterion did not name it.

**The defect appeared at the composition, not in the leaf.** Both defects in this record
(section 9 and this one) were found by building something *on top of* the module, after its own
suite was green. That is an argument for the validation ladder having a Layer 3 at all, and it is
now two instances rather than one.
