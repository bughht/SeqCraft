# RadialReadout Fine Scan — Boundary Decision, Extraction and Status

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** `RadialReadout`
**Date:** 2026-09-18
**Passes:** 7 — boundary; extraction; validation; sweep; status
**Decided by:** human review of the four questions in [`findings.md`](findings.md) §6

---

## 1. The decision

```text
NEW_LEAF   readout/RadialReadout
```

and, explicitly, **nothing else**:

| | |
|---|---|
| `RadialGRETR` | **not created.** Radial GRE today is `Excitation` → `RadialReadout` → spoiler → delay, which is composition. A kernel is justified by a shared physical or timing solve that no leaf can do alone — what `GRE2DTR` and `TSEShot` exist for — not by a recurring arrangement |
| `seqcraft/trajectory/` | **not created.** The numerical design fits inside the leaf. Revisit only if Spiral independently needs the same contract |
| `Refocusing`, the compiler, `LogicBlock` | untouched |

| question | decision |
|---|---|
| Q1 full spoke vs centre-out | **one leaf**, spanned by `partial_fourier` in `[0.5, 1.0]` |
| Q2 angle / orientation | **the leaf's**, as a per-call `build(angle_rad=…)` argument, emitting already-oriented events |
| Q3 semantic `k = 0` | **the leaf's**, exposed rather than reverse-engineered |
| Q4 spoke ordering / golden angle | **caller's**, above the leaf |

### The rule, restated

```text
leaf     intrinsic physics / geometry
kernel   coupling between leaves
imaging  acquisition policy
compiler Pulseq legality
```

> Repeated composition is not enough to justify a kernel. Repeated cross-leaf coupling that must
> be solved jointly may justify one.

## 2. What was built

SeqCraft branch `feat/module-mining-radial`:

| | |
|---|---|
| `src/seqcraft/modules/readout/radial_readout.py` | `RadialReadout` |
| `tests/modules/test_radial_readout.py` | 29 tests |
| `docs/api_reference.md`, `CHANGELOG.md` | index, a worked example, the argument |
| `tools/module_mining/candidates/radial/` | adapters, trajectory checks, archaeology and validation runs |

**It composes `CartesianLine`.** A 2D spoke is a Cartesian line pointing somewhere other than
along an axis: the prephaser, the dwell arithmetic, the ADC and which sample carries `k = 0` are
the same problem, already solved and already tested. The new module adds the rotation, the
trajectory geometry and a contract stated as a spoke. This is the architectural-economy rule in
`AGENTS.md` applied literally — an existing definition stays the single source of truth.

A consequence worth stating: `partial_fourier` was not chosen to *resemble* the existing
parameter, it **is** it. `CartesianLine._resolve_samples` already computes
`pre_echo = round(pf*matrix) - (matrix - matrix//2)`, which at `matrix=64` gives centre sample 32
for `pf=1.0` and 0 for `pf=0.5` — the reference's full spoke and a centre-out spoke, with no new
arithmetic.

## 3. Validation

### Against the official reference, sample for sample

Matched geometry (260 mm, 64 samples, 20 µs dwell), 8 spokes, **not** a matched construction —
the reference hand-builds trapezoids and calls `pp.rotate`; the candidate derives rotated copies
of a `CartesianLine`:

```text
worst k-space difference over all 8 spokes:  4.96e-12  1/m
centre sample:                               32 in both, every spoke
```

Two implementations sharing no code and agreeing on every sample. Per §10 of the phase plan,
event-for-event identity was **not** the criterion; the trajectory is.

### The invariant table, and the sweep

`python tools/module_mining/candidates/radial/run_validation.py` — **8/8 cases pass every
invariant** (a sample on the centre, straightness, uniform spacing, consistency across spokes, the
requested angle, staying in plane, and rotation equivariance):

| case | samples | centre | Δk 1/m | rotation error |
|---|---|---|---|---|
| reference protocol | 64 | 32 | 3.8462 | 4.3e-14 |
| matrix 128 | 128 | 64 | 3.8462 | 1.0e-13 |
| matrix 96 | 96 | 48 | 3.8462 | 5.7e-14 |
| fov 200 mm | 64 | 32 | 5.0000 | 8.5e-14 |
| dwell 10 µs | 64 | 32 | 3.8462 | 4.3e-14 |
| dwell 40 µs | 64 | 32 | 3.8462 | 3.6e-14 |
| partial 0.75 | 48 | 16 | 3.8462 | 2.8e-14 |
| centre-out 0.5 | 32 | 0 | 3.8462 | 1.4e-14 |

In every case the module's **claimed** `k_first_per_m`, `dk_per_m` and `center_sample` match the
compiled trajectory.

**The rotation test is the one written before the module existed.** It passed against the
official PyPulseq reference during the archaeology pass, unchanged, which is what stops it
encoding the candidate's own behaviour. It is now also a package test.

### Repository gates

`ruff check .` clean · **938 passed / 19 skipped** (902 before, so 36 new tests and nothing
broken) · 73 doctests · strict mypy subset clean · `check_api_reference.py` clean with the new
example executing · all nine build notebooks execute. **No compiler change.**

## 4. Where the code pushed back

**4.1 `require_range` is half-open and 0.5 had to be legal.** The shared helper checks
`low < value <= high`, deliberately, because the fractions it bounds are degenerate at their
lower end. 0.5 is not degenerate for a spoke — it is centre-out. Rather than add an `inclusive`
flag to a helper with one such caller, the check is written in the module, where the message can
say what 0.5 *means*.

**4.2 The first harness measurement was wrong, and the module was fine.** Stacking bare spokes in
one tree reported a rotation error of **5.8e+02 1/m**. `k` is integrated from the start of a
sequence and reset only by an excitation, so each spoke began wherever the last ended. The
adapter now compiles one spoke per sequence and says why. A harness that sums its own
trajectories is worse than no harness — and this is the second time in this project that the
comparator, not the candidate, was the thing that was wrong.

**4.3 An odd matrix did not compile. Resolved by review: refuse early, in the layer that
computes the count.** `matrix=65` gives 65 ADC samples, which is not a multiple of the scanner's
`adc_samples_divisor`; the failure used to arrive from the compiler, naming a block index rather
than the arguments responsible.

`CartesianLine` now refuses it at construction, and `RadialReadout` inherits that by composing the
line rather than keeping a second copy of the rule. Two points settled with it:

- **the rule is on the count, not on the parity of the matrix.** An even matrix produces an
  illegal count just as easily — 0.6 of 128 is 77, which is why `examples/fse_2d` runs HASTE at
  0.625, a fact the notebook stated in prose and nothing enforced until now;
- **nothing is rounded.** Moving 65 to 64 would quietly change the partial Fourier fraction, the
  k-space extent, the semantic centre sample and the sampling density. The error names both
  arguments and suggests only combinations the module would itself accept — tested, because a fix
  that walks the caller into the next refusal is worse than none.

This is the ownership rule applied to validation rather than to events: *the layer that determines
a quantity validates the constraints on it*, with the compiler still the final backstop.

**4.4 Rotating a trapezoid scales its stored `area`.** Both references that extend a readout's
flat time warn in a comment that their event's area becomes wrong. Scaling amplitude alone would
reproduce that defect in a module whose whole purpose is to be asked what it plays, so `area` and
`flat_area` are scaled too, and a test asserts it.

## 5. Status

```text
GREEN
```

Against §11 of the phase plan: the boundary is clear; two executable Python references exist and
run; provenance and licences are recorded; waveform and trajectory invariants pass; the rotation
metamorphic test passes; the sweep passes; the package suite passes; no compiler change was
required.

The one rough edge found during the sweep — an odd matrix failing at the compiler rather than at
the module — was resolved rather than recorded; see §4.3.

Not done, and not gates under the Python-only constraint: no MATLAB `.seq` was generated for this
candidate, and 3D radial, ramp sampling and `rot3D` are deliberately out of scope.

## 5b. The example, added at v0 close-out

The candidate shipped with tests but with no example, which was the one coverage gap the v0
close-out found: every other module in `sc.modules` has a notebook where a reader can see what it
does, and `RadialReadout` did not.

`examples/gre_radial_2d/01_build.ipynb` closes it, and it is a **build and trajectory-visualisation
notebook only**. What it shows:

| | |
|---|---|
| the repetition is assembled **in the notebook** | `Excitation` + the spoke + a spoiler + TR fill, four lines, no new class |
| the angle schedule is **ordinary Python** | `equal_increment` and `golden_angle` are two list comprehensions over the same readout instance |
| every spoke has a sample **exactly on `k = 0`** | `1.98e-08` 1/m measured off the compiled sequence, which is integration noise |
| the spokes point where they were asked to | `1.99e-13` degrees maximum error across 64 spokes |
| `partial_fourier` spans the family | `1.0` full spoke, `0.5` centre-out, `Δk` unchanged, `0.25` refused |

This re-confirms §1's decision from the other direction. The notebook was written to see whether a
`RadialGRE` would have earned its place, and the composition it would have wrapped is four lines
of which three are shared with every other GRE — so what the class would actually own is the angle
list, which is a sequence choice and belongs to the caller.

### Why there is no `02_simulate_and_reconstruct.ipynb`

**The status is layered, and §5's `GREEN` belongs to the first two layers only:**

```text
Layer 1 -- physical/reference validation      GREEN
Layer 2 -- complete acquisition build         GREEN
Layer 3 -- simulate + reconstruct             DEFERRED
```

So the claim this document establishes is exactly:

> the radial spoke leaf emits the correct waveform / ADC / trajectory semantics, independently
> validated against the official reference.

and it is **not**:

> SeqCraft has an independently validated end-to-end radial imaging + non-Cartesian reconstruction
> pipeline.

The deferral is a decision, not an unfinished v0 task, and it should not be closed for symmetry
with FSE and GRE3D. A radial Layer 3 needs a non-Cartesian reconstruction path whose own
conventions -- NUFFT coordinates, density compensation, weighting, trajectory timing -- each need
independent review first; an image produced through an unreviewed reconstruction is two unvalidated
things agreeing, not evidence about the sequence. PR #23 holds usable non-Cartesian reconstruction
infrastructure, but under the roadmap it is **evidence for the Spiral supervised-calibration phase,
not accepted architecture**, and nothing in `examples/` may depend on it until that review happens.
It is picked up in `v0_closeout_and_roadmap.md` §2.7.2 and §6.1, where the notebook is a
**conditional** on that review rather than a scheduled deliverable: if no genuinely shared
reconstruction utility survives, the radial `02` is never written and this module stays Layer-2
GREEN with the claim above.

Deferred, then, with the reason recorded rather than left implicit. A radial image needs a non-Cartesian
reconstruction; the example suite has none — every existing `02` is an FFT on a Cartesian grid, and
a search for NUFFT or gridding machinery across `examples/` found only the word in prose. Building
a gridding or NUFFT pipeline here would be **new reconstruction infrastructure created to complete
a pair**, which is the kind of scope the fine-scan playbook exists to refuse.

What that costs is worth stating plainly: the geometric claims above are measured on the compiled
trajectory, so **no image has been reconstructed from a SeqCraft radial acquisition.** The
sample-for-sample agreement with the official reference in §3 is the stronger evidence and it
stands on its own, but it is agreement with another sequence, not a picture. If a non-Cartesian
reconstruction arrives for some other reason, this is the notebook to pair with it.

## 6. Before the next candidate

The phase plan asks for the two pilots to be compared before anything is formalised. The short
version, for that conversation:

- **What repeated:** provenance first; an invariant table before the API; the ownership question
  *which layer can determine this value*; a metamorphic test defined against a reference before
  the candidate existed; and a comparator that returns measurements rather than verdicts.
- **What differed:** TSE needed event-for-event identity because it was promoting an existing
  implementation; Radial could not use it, because two legitimate constructions of the same
  trajectory share no events. The acceptance criterion is candidate-specific and should not be
  frozen into the tooling.
- **What was wrong both times:** the comparator. Once about where the echo is, once about which
  axes an area rule applies to, and once about summing trajectories across spokes. Every one was
  caught by measuring something whose answer was already known.
