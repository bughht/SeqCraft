# GRE3DTR — Acceptance Criterion

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** `GRE3DTR`
**Date:** 2026-09-19
**Written before extraction, deliberately.** Nothing in this file was chosen after seeing what an
implementation happened to produce.

> **On the labels.** The bracketed IDs are cross-references for this plan and nothing more. The
> name is the concept; the number is a locator. Tests, package documentation and review notes use
> the names.

---

## 1. The claim this candidate makes

> `GRE3DTR` is an **independent implementation** of a 3D Cartesian gradient-echo repetition,
> reproducing the physical encoding of the Pulseq references without sharing their event
> structure.

**It is not a promotion**, and that is a correction to the two-pilot retrospective, which guessed
this candidate would be "promotion-shaped … so event identity is probably right". The evidence
says otherwise: **SeqCraft contains no 3D implementation to preserve.** `examples/mprage_2d/` and
`examples/mp2rage_2d/` are 2D GRE trains under `IRPrep`. There is no notebook-local `GRE3DTR`, no
`GRE3D`, and nothing whose events could be held fixed.

So this candidate is **Radial-shaped, not TSE-shaped**, and the retrospective's own reasoning
applies: event identity is the right criterion for a promotion and the wrong one for a
reimplementation.

## 2. Therefore: what preservation means here

**Event identity is explicitly NOT the criterion.** Neither is block count, nor event
decomposition, nor agreement with the references' hand-built trapezoid arrays.

Nor, and this one is specific to these references, is **raw sample-position agreement on x**: R1's
prephaser is `-gx.area/2` with no half-sample correction, so its samples straddle `k = 0` by
Δkx/2, where SeqCraft's `CartesianLine` puts a sample there. That is a design difference already
recorded twice in this project. The comparison is therefore on the **lattice and the semantic
centre**, not on where sample *n* lands.

### The contract, stated as checks

| | check | tolerance |
|---|---|---|
| **Readout refocused at the echo** (A1) | `kx = 0` at the module's own `echo_sample`, every repetition | 1e-3 1/m |
| **Requested line encoded** (A2) | `ky` at the echo is `(line - matrix//2) * 1/fov_y`, exactly the requested line | 1e-3 1/m |
| **Requested partition encoded** (A3) | `kz` at the echo is `(partition - matrix_z//2) * 1/fov_z`, exactly the requested partition | 1e-3 1/m |
| **kz held through the readout** (A4) | **`kz` is held constant from the first sample to the last** — it is an encode, not a rephase | 1e-3 1/m |
| **Both encodes rewound** (A5) | both encodes are rewound before the next excitation: `ky` and `kz` return to 0 | 1e-3 1/m |
| **Winder axes start together** (A6) | the three winder axes **start together**, and the window is the widest of the three minima | 1 ns |
| **Reported TE and TR achieved** (A7) | claimed `te_s` / `tr_s` match the compiled sequence | 1 gradient raster |
| **Block is exactly TR** (A8) | the block is exactly `tr_s` long | 1 ps |
| **Hardware legality** (A9) | hardware legality — timing check, gradient and slew limits | exact |
| **Lattice agrees with the reference** (A10) | against R1's `gre3d.seq` at matched geometry: same `Δk` per axis, same index range per axis, same number of distinct indices, each acquired once | 1e-3 1/m |

### Metamorphic tests, to be written against the *reference* first where possible

| | claim |
|---|---|
| **Partition index maps linearly to kz** (M1) | the partition index maps linearly to `kz`, with the sign the caller asked for |
| **Centre partition reaches kz = 0** (M2) | the single partition at `matrix_z // 2` gives `kz = 0` at the echo |
| **More partitions move none of the existing ones** (M3) | increasing `n_z` does not move the `kz` of partitions that already existed |
| **Swapping the y and z FOV swaps their dk** (M4) | swapping the y and z FOV components swaps `Δky` and `Δkz` and nothing else |

The first two are measurable on `gre3d.seq` before any candidate exists, and should be, for the
reason the radial rotation test was: a metamorphic test written after the candidate tends to
encode the candidate.

### The three degenerate limits that validate the combined-z solve

These are the strongest tests this candidate has, because each one collapses the combined
arithmetic to something whose answer is known without it. **They are central invariants, not
extras**, and they are what makes the selective path acceptable without an external
implementation that happens to share SeqCraft's construction.

| limit | | required result |
|---|---|---|
| **Non-selective limit** (L1) | `slab_thickness_mm = None` | `A_slab_rephase = 0`, so the combined z winder **is** the pure partition encode — which ties the selective implementation to the official non-selective references, measured against R1 |
| **Centre-partition limit** (L2) | `partition = matrix_z // 2` | `A_partition = 0`, so the combined z winder **is** the pure slab rephaser, and `kz = 0` at the echo |
| **Adjacent-partition step** (L3) | `partition p -> p + 1` | the combined z moment changes by **exactly `Δkz`**, whatever the slab rephasing contributes |

The adjacent-partition step is the one that catches an error the other two cannot: a slab-rephase term that is wrong by a
constant passes the centre-partition limit and still shifts every other partition, and the
non-selective limit never exercises it at all. Together the three pin the affine relation `A_z(p) = A_rephase + (p - centre)·Δkz`
completely.

### Polarity and the limiting partition

Both z terms are signed, so which partition is limiting is a *result* of the coupling rather than
a property of the index. These tests exist to ensure the solver is built on signed physical
moments and not on an accidental array convention.

| | claim |
|---|---|
| **Limiting edge, negative slab term** (P1) | with the slab rephasing negative, the limiting partition is one `kz` edge |
| **Limiting edge, positive slab term** (P2) | with it positive, the **opposite** edge becomes limiting — same geometry, same partitions |
| **Limiting edge switches with the offset** (P3) | a slab design change that moves the constant offset far enough **switches** the limiting edge, and the solver follows it |
| **Sign from semantics, not index order** (P4) | reversing the partition sign convention does not change which *physical* `kz` each partition reaches |
| **Reported limiting partition is the real one** (P5) | the limiting partition the module reports is the one whose `t_z_min` is actually largest, enumerated over every supported partition |

The switching case is the one that would catch a solver right for the wrong reason: a "largest
index wins" shortcut passes the two fixed-polarity cases if the lattice happens to be signed that
way, and fails this one.

### Timing behaviour, which is a stated contract and not an emergent one

| | claim |
|---|---|
| **Implicit timing is the shortest legal design** (T1) | with `te_s=None` and `tr_s=None`, the module returns the **shortest legal** design and reports what it achieved |
| **Window lengthens automatically** (T2) | when the worst partition needs a longer z winder than x or y, the shared window is **lengthened automatically** — the module does not refuse, and the compiler is not asked to rescue it |
| **TE independent of kz** (T3) | **every partition uses the same winder duration**, so TE does not depend on `kz`. Measured across the centre and both edges |
| **Impossible timing refused with the reason** (T4) | an explicit `te_s` or `tr_s` below the achievable minimum **raises**, naming the limiting partition, its combined moment and the minimum z-winder duration |
| **Compiled-TE agreement** (T5) | nothing downstream of the kernel changes its timing: the compiled TE equals the TE the module reported |

TE independence is not cosmetic. A per-partition window would make TE a function of `kz` — a contrast gradient
across the volume that no reconstruction expects and that no k-space check would show.

### Feasibility, which the sweep has to reach

The sweep must include a protocol where the combined moment at the worst partition does *not* fit
the window the other two axes need, so that the automatic lengthening is exercised rather than
assumed.

### Parameter sweep

```text
n_y, n_z            two sizes each, including unequal
fov_z               two slab FOVs -- dkz comes from the z FOV component, not the in-plane one
slab_thickness_mm   None; equal to fov_z; less than fov_z (R4b's 0.7); greater than fov_z
slab polarity       both signs of the rephasing term, so the limiting edge switches
rf_time_bw_product  the factory default, and one sharper value
te_s, tr_s          minimum and longer
partition sign      both, since the reference's negation is policy
partition index     centre, both edges
```

## 3. The two paths have different evidence, and therefore different criteria

They are one kernel, but what can be measured about each differs, and pretending otherwise would
either weaken the non-selective path or invent evidence for the selective one.

### Non-selective path — validated against an executable official reference

R1's `gre3d.seq` is read directly, so every check above applies: the three lattices, the semantic
centres, ADC timing, TE/TR, rewinding and the shared winder window, compared on **lattice and
semantic centre** rather than on sample positions.

### Slab-selective path — validated analytically, with R4b as corroboration

R4b is MATLAB and ships no `.seq`, so it is read for physics rather than measured. The criterion
is therefore:

| | check |
|---|---|
| **RF effective centre** (S1) | the RF effective centre is where the module says, and TE is measured from it |
| **Combined z moment from the RF centre** (S2) | the z moment accumulated from the RF centre to the echo is `A_slab_rephase + A_partition(p)`, measured by integrating the emitted waveform |
| **kz at the echo, held** (S3) | `kz` at the echo is the requested partition's, and is held through the readout |
| **Partition spacing** (S4) | `Δkz` between adjacent partitions is `1/fov_z` |
| **Degenerate limits** (S5) | limits L1, L2, L3 above |
| **Legality at the edge partitions** (S6) | gradient amplitude and slew legality at the edge partitions, where the combined moment is largest |
| **Minimum TE not longer than sequential** (S7) | the minimum TE is **not longer** than a sequential construction would need — the claim the combined solve is making |
| **Slab rephasing realised once** (S8) | the slab rephaser is realised **once**: the z moment from the RF's effective centre to the echo is `A_slab_rephase + A_partition(p)` and not `2·A_slab_rephase + …` |

The echo-time claim is the one worth stating out loud: R4b's own `te_min` pays for `calcDuration(gzReph)` *and*
`calcDuration(gzPre)`. If the combined solve does not remove that term, it has bought nothing and
the simpler sequential construction should win.

**Event identity is not used on either path**, and neither is agreement with R4b's two-event z
decomposition. External references validate physical claims; they do not dictate SeqCraft's event
decomposition.

## 4. What would make this YELLOW or RED instead

- **YELLOW** if the combined solve turns out to be infeasible at protocols users would reasonably
  ask for, so that the kernel refuses more often than it helps.
- **YELLOW** if `GRE3DTR` needs nothing `GRE2DTR` does not already own once the selective path is
  set aside — the case for a sibling kernel would then rest on a second encode axis alone.
- **RED** if the shared-window solve cannot be expressed without reaching into `CartesianLine`'s
  or `PhaseEncode`'s already-designed events — that would mean the nesting boundary is wrong,
  which is the warning sign the architecture review named.
- **RED** if the non-selective limit fails: a selective implementation that does not reduce to the non-selective one at
  zero slab rephasing is not the same kernel.
- **RED** if the standalone slab rephaser is emitted alongside the combined winder. It is a
  doubled moment, it is invisible in a k-space extent check, and it is the specific failure the
  architecture has been arranged to prevent.

## 5. What this criterion costs

Every check above is measurable with the existing comparator plus two additions: a
**partition axis**, which the Cartesian checks currently treat as "the axis that should be quiet",
and **moment integration from an RF effective centre**, which already exists for TSE's crusher
balance and needs pointing at z. Both are family-specific content in the sense the retrospective
defined — the frame is reused. Neither is a new comparator.

The three degenerate limits need no comparator at all: they are relations between two builds of the same
module, which is the cheapest and most durable evidence this candidate has.
