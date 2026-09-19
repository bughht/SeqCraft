# GRE3DTR — Extraction Report and Status

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** `GRE3DTR`
**Date:** 2026-09-19
**Passes:** extraction; validation; sweep; status
**Acceptance criterion:** written **before** any of this, in
[`acceptance_criterion.md`](acceptance_criterion.md)

---

## 1. What was built

| | |
|---|---|
| `src/seqcraft/modules/kernel/gre_3d_tr.py` | `GRE3DTR` — one repetition, both excitation modes |
| `src/seqcraft/modules/rf/excitation.py` | `+ rephaser_area_per_m`, `+ build(rephase=…)` |
| `tests/modules/test_gre_3d_tr.py` | 27 tests |
| `tests/modules/test_excitation.py` | 5 added, including the backward-compatibility pin |
| `tools/module_mining/candidates/gre3d/` | validation against `gre3d.seq`, and the sweep |
| `docs/api_reference.md`, `CHANGELOG.md` | the index, a worked example, the argument |

`GRE2DTR` is untouched. So is the compiler, `LogicBlock`, `PhaseEncode` and every other module.

## 2. The ownership split, as built

```text
Excitation      states the rephasing its own slab implies   (rephaser_area_per_m)
                and can decline to realise it               (build(rephase=False))
PhaseEncode     states the signed moment a partition wants  (k_per_m)
GRE3DTR         holds both, so it owns their combined realisation
compiler        legalises the resulting timeline, and changes no timing
```

`Excitation`'s new property is the same shape `CartesianLine.area_to_echo_per_m` already had, and
for the same reason: a composite combining two axes' requirements should read a **physical
moment**, not reach into a pypulseq event and reinterpret it. It is derived by integrating the
selection gradient from the RF's effective centre rather than by reading `gzr`, and the two agree
to **6e-14** — which is evidence, because they are not the same computation.

## 3. Validation

### Non-selective path, against `pulseq/pulseq`'s shipped `gre3d.seq`

No MATLAB run. Compared on **lattice and semantic centre**, not sample position, because the
reference prephases with `-gx.area/2` and has no sample at the centre of k-space where
`CartesianLine` puts one.

| | reference | `GRE3DTR` |
|---|---|---|
| `Δk` x / y / z (1/m) | 5.0000 / 5.0000 / 6.2500 | **5.0000 / 5.0000 / 6.2500** |
| encoded index at the probed lines | `ky` −32 … +31 | −32, 0, +31 as asked |
| encoded index at the probed partitions | `kz` −31 … +32 *(mirrored: its own vendor sign flip)* | −32, 0, +31 as asked |
| `kz` spread through the readout | held | **0.0e+00** |
| `|kx|` at the semantic echo | ±2.500 — *no DC sample* | **2.7e-07** — a sample at DC |

### Both paths, in the package suite

All of L1/L2/L3, P1–P5 and T1–T5 pass as written:

| | result |
|---|---|
| **L1** non-selective limit | `A_slab = 0`, and the combined winder equals the partition encode for every partition |
| **L2** centre partition | combined moment is the pure rephasing; measured `kz` at the echo is zero |
| **L3** adjacent partitions | differ by exactly `Δkz`, both modes, the slab term cancelling |
| **P1/P2** | `A_slab = ∓120` moves the limiting index between the two ends — your worked example, reproduced through the solver |
| **P3** | sweeping the offset moves the limit monotonically across; a "largest index wins" shortcut fails this and passes P1/P2 |
| **P4** | the sign comes from `PhaseEncode.k_per_m`, not from array order |
| **P5** | the reported limiting partition is the one whose duration is actually largest, enumerated |
| **T1** | implicit timing returns the shortest legal design |
| **T2** | a thin slab with many partitions grows the winder **260 → 1110 µs** and TE with it, without refusing |
| **T3** | TE is identical at the centre and both edges, to < 1 ns |
| **T4** | an impossible TE raises, naming the limiting partition and its combined moment |
| **T5** | the block is exactly `tr_s`, and the compiled TE equals the reported one |
| **S8/G19** | the excitation contributes one z gradient, not two, and `kz` at the echo is the partition's — a doubled slab term would shift it by 12.5 1/m here |

### Sweep

| case | winder | TE | `A_slab` |
|---|---|---|---|
| reference geometry 64³ | 410 µs | 3.991 ms | 0 |
| unequal `ny`/`nz` | 260 µs | 3.863 ms | 0 |
| **thin slab, 64 partitions** | **1110 µs** | **4.713 ms** | 0 |
| slab = `fov_z` | 260 µs | 3.863 ms | −12.54 |
| slab < `fov_z` (0.7×) | 260 µs | 3.863 ms | −17.92 |
| slab > `fov_z` (1.25×) | 260 µs | 3.863 ms | −10.03 |
| sharper slab, TBW 8 | 260 µs | 3.863 ms | −25.08 |

The third row is the one that matters: z genuinely forces the window there, so T2 is exercised
rather than merely present. The last row shows the RF forwarding working — doubling the
time-bandwidth product doubles the rephasing moment.

### Repository gates, from the integrated tree

`ruff check .` clean · **985 passed / 19 skipped** · 76 doctests · strict mypy subset clean ·
`check_api_reference.py` clean with the new example executing · all ten notebooks execute.
**No compiler change.**

## 4. What the implementation found rather than assumed

**4.1 TE was short by exactly the winder duration.** `min_te_s` added `winder_s` on top of
`ro.time_to_echo()`, which already carries the prephaser it was designed with. Caught by T5
comparing the *compiled* echo time against the reported one — a test that would have passed had
it only compared the design against itself.

**4.2 The z spoiler collides with the partition rewinder.** Inheriting `GRE2DTR`'s
`spoil_axis=('x', 'z')` put a spoiler on the axis that already carries the rewinder, and the two
summed past the slew limit — the compiler refused, correctly. The default is now `('x',)`, which
is what `writeGradientEcho3D.m` does, its tail being `gyReph, gzReph, gxSpoil`. Asking for `'z'`
is still allowed and may need a derated design.

**4.3 The signed solver is a pure function.** `_limiting_index(areas, opts)` exists apart from the
class because that is the part worth testing on its own: a test hands it a slab term of either
sign and watches the answer move between the ends of the list, which is precisely what a
shortcut would get wrong.

## 5. Status

```text
GREEN
```

Boundary decided by review before extraction; the acceptance criterion written before
implementation and met on both paths; an executable reference for the non-selective path and
analytic limits for the selective one; sweep passing including a z-limited protocol; package
suite clean; no compiler change.

The selective path still has **no executable reference** — `fmrifrey/lps` is MATLAB and ships no
`.seq`. That is by design of the criterion rather than an outstanding gap: its claims are carried
by signed-moment arithmetic and the three degenerate limits, which do not depend on any external
implementation. One MATLAB run would corroborate G10/G11 and remains optional.

## 6. The end-to-end layer

The analytic and compiled-sequence checks above answer *does the design satisfy the contract* and
*did the emitted sequence preserve it*. They do not answer *does the whole acquisition behave like
the MRI experiment intended*, and the two bugs this extraction actually had — a reported TE short
by one winder, and a z spoiler colliding with the partition rewinder — are a reminder that local
agreement is not the same as a working scan.

So the ladder is restored for this candidate:

| layer | where | gate? |
|---|---|---|
| analytic invariants and regressions | `tests/modules/test_gre_3d_tr.py` | **blocking CI** |
| the example executes | `examples/gre_3d/01_build.ipynb`, in the notebook smoke tier | **blocking CI** |
| the reconstruction looks like an MRI experiment | `examples/gre_3d/02_simulate_and_reconstruct.ipynb` | review evidence, not a gate |

`examples/gre_3d/01_build.ipynb` builds a complete acquisition from `GRE3DTR` and two `for` loops
— which is itself the argument in §7 — and prints the signed coupling partition by partition.
`02_simulate_and_reconstruct.ipynb` reconstructs the volume with a 3D FFT and shows it in three
planes, which is where a reversed `kz` or a `ky`/`kz` swap stops being a number.

**No rendered image is a correctness criterion.** The two assertions in that notebook are the
reconstructed shape and that most of the energy lies inside the phantom; everything else is for a
human to look at.

### What the visual layer added, concretely

Two things the analytic checks had not covered, both found while writing it:

- **The sequential arrangement is inherently longer than the combined one**, so it cannot be
  stacked at the kernel's own minimum TR. A first attempt did, every repetition overlapped the
  next by exactly the 100 µs being measured, and the compiler's first-moment check refused it.
  Both scans now run at one explicit TR.
- **The comparison is only fair with the spoiler present.** Without it the steady state, not the
  z realisation, was the difference — a first run reported a 63 % discrepancy that had nothing to
  do with the thing under test.

Neither is a defect in `GRE3DTR`. Both are exactly the class of mistake an end-to-end layer exists
to surface.

### The measured claim

At one TR, one geometry and one phantom, the **combined** z winder and a **sequential**
slab-rephase-then-encode arrangement reconstruct the same volume to **0.0014 of the peak**
(0.0001 mean), while the combined one reaches the echo **100 µs earlier per repetition**.

> SeqCraft changed how the required z moment is realised, not what acquisition is performed.

## 7. `GRE3D` (imaging) is **not** proposed

The same rule that produced this kernel withholds the layer above it, and
`examples/gre_3d/01_build.ipynb` is the evidence: a complete 3D acquisition there is the kernel
plus two ordinary `for` loops. Partition and line
ordering, dummies and the RF-spoiling schedule are policy and would belong there — but there is
no consumer yet, no second arrangement to generalise from, and no evidence that the ordering
question is harder in 3D than `GRE2D` already answers. A kernel existing is not a reason for an
imaging class to exist.
