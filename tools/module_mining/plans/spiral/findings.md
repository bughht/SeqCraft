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
