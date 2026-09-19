# GRE3DTR Fine Scan — Invariant Table

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** `GRE3DTR`
**Date:** 2026-09-19
**Pass:** 3 — the physical contract, written before any API

Provenance: [`reference_inventory.md`](reference_inventory.md). Evidence:
[`findings.md`](findings.md). What acceptance will require:
[`acceptance_criterion.md`](acceptance_criterion.md).

> **On the labels.** The `G` numbers are cross-references for this plan. The invariant is the
> concept; the number is a locator.

## How the measured column was obtained

`matlab/demoSeq/gre3d.seq` — shipped by `pulseq/pulseq`, written by `writeGradientEcho3D.m`, read
with `pypulseq.Sequence.read` and measured through `calculate_kspacePP`. **No MATLAB run.**
FOV 200 × 200 × 160 mm, 64³, so `Δkx = Δky = 5.000` and `Δkz = 6.250` 1/m. 4096 readouts of 64
samples.

---

| invariant | ref | physical meaning | how measured | tolerance | R1 measured |
|---|---|---|---|---|---|
| **`kx = 0` at the semantic echo** | G1 | the readout is refocused | k at the module's own centre sample | 1e-3 1/m | **no DC sample**: ±2.500 1/m at samples 31/32, exactly ±Δkx/2 |
| **`ky` on the `1/fov_y` lattice, each index once** | G2 | in-plane phase encoding | index = `round(ky/Δky)`, and its residual | 1e-3 1/m | 64 distinct, **−32 … +31**, each 64 times, off-lattice by **1.5e-4** |
| **`kz` on the `1/fov_z` lattice, each index once** | G3 | partition encoding; `Δkz` comes from the **slab** | as the row above | 1e-3 1/m | 64 distinct, **−31 … +32**, each 64 times, off-lattice by **9.6e-5** |
| **`kz` is held through the readout** | G4 | z is an *encode*. A slab rephaser would return it to zero at the echo; this is the measurement that tells the two apart | `kz` at the first sample, the echo, the last | 1e-3 1/m | **200.00 → 200.00 → 200.00** 1/m |
| **both encodes rewound before the next excitation** | G5 | otherwise k accumulates across TRs | k at the start of the next repetition | 1e-3 1/m | *to measure* |
| **the winder axes start together, in one window** | G6 | the coupling the kernel would own | event start times on x, y, z | 1 ns | R1 fixes `Tpre = 3 ms` by hand on all three |
| **TE measured from the pulse's **effective centre** | G7 | a block pulse's centre is not the end of its block | echo time minus excitation time | 1 gradient raster | R1 computes it that way; *to measure* |
| **the repetition is exactly `TR`** | G8 | a caller stacks them | block duration | 1 ps | *to measure* |
| **hardware legality** | G9 | it must run | timing check, gradient and slew limits | exact | R1's `.seq` reads and computes cleanly |

## The slab-selective rows

`fmrifrey/lps` (R4b) is the only witness, it is MATLAB and it ships no `.seq`, so these are
**stated from its source and from arithmetic**, not measured. They are what the selective path's
acceptance criterion tests.

| invariant | ref | physical meaning | how measured | tolerance |
|---|---|---|---|---|
| **the slab-selection gradient leaves a rephasing requirement of `-gz.area/2` from the RF's effective centre** | G10 | a selective pulse dephases across the slab and must be undone | integrate z from `time_to_center()` to the echo | 1e-3 1/m |
| **`A_z(p) = A_slab_rephase + A_partition(p)`** — the two requirements land on one axis in one window** | G11 | the coupling that makes this a kernel | integrate the emitted z waveform between the RF centre and the echo, per partition | 1e-3 1/m |
| **`A_z(p+1) - A_z(p) = Δkz` exactly, whatever the rephasing contributes** | G12 | the affine relation, which a constant error in the rephasing term would break | two builds, one subtraction | 1e-3 1/m |
| **at `slab_thickness_mm = None` the z winder is the pure partition encode** | G13 | ties the selective implementation to the non-selective references | compare against the same module built non-selective | exact |
| **at the centre partition the z winder is the pure slab rephaser, and `kz = 0` at the echo** | G14 | the other degenerate limit | as the combined-moment row | 1e-3 1/m |
| **the combined window is **no longer** than the partition encode alone would need** | G15 | the claim the combined solve makes; R4b pays for two windows and its `te_min` shows it | minimum TE, both constructions | — |
| **amplitude and slew legal at **every** partition** | G16 | the largest combined moment is wherever the signed sum puts it, not at a known index | limits check over all supported partitions | exact |
| **the limiting partition follows the **sign** of the rephasing term** | G17 | reversing the slab gradient's polarity moves the limit to the opposite `kz` edge | worst `t_z_min` over all partitions, both polarities | — |
| **one winder duration serves every partition** | G18 | otherwise TE is a function of `kz`, which is a contrast gradient across the volume | compiled TE at the centre and both edges | 1 ns |
| **the slab rephasing is realised **once** | G19 | a doubled moment passes every extent check | z moment from the RF effective centre to the echo | 1e-3 1/m |

`slab_thickness_mm` is deliberately **not** `fov_z` in any of these, and neither bounds the other:
R4b excites 0.7 of its FOV, a full-volume protocol commonly excites at least the FOV and often
more, and the only validation is that a selective slab is positive. A row asserting
`slab >= fov_z` would fail against the only selective reference there is.

## The two conventions worth naming

**The z lattice is mirrored relative to y** — indices −31…+32 against −32…+31 — because R1
negates `areaZ` for a vendor reconstruction, with the source comment saying so. `Δk`, the index
count and the extent are all identical; only the sign differs. that row therefore checks the lattice and
**not** the sign, and the sign is caller policy (findings §6).

**`Δkz` is `1/slab`, not `1/FOV`.** 6.250 against 5.000 1/m at R1's protocol, because the FOV is
a vector and its z component is the slab thickness. Anything that assumes a scalar FOV gets this
wrong by 25 % here and produces a plausible image.

## What this table deliberately does not contain

- **Event structure, block count, or the references' per-line trapezoid arrays.** R1 pre-builds
  one trapezoid per phase-encode line to halve its own construction time; that is a workaround
  for a flat event list.
- **The partition ordering.** Which partitions are acquired, in what order, with how many
  dummies, and with which sign, is acquisition policy.
- **How the z moment is decomposed into events.** R4b plays a rephaser and then a scaled encode;
  SeqCraft intends one combined winder. the combined-moment row asks what the z axis *achieves* between the RF centre
  and the echo, which both constructions can answer.

## Not yet measured

The rewinding, RF-centre TE and block-duration rows need a candidate or a second pass over R1,
and every slab-selective row needs a candidate. R4b could corroborate the slab-rephasing and
combined-moment rows with one MATLAB run, which is optional evidence rather than a gate. They are
all listed so the criterion is complete rather than convenient.
