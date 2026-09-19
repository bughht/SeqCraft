# GRE3DTR Fine Scan — Code Archaeology and Semantic Decomposition

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** `GRE3DTR`
**Date:** 2026-09-19
**Passes:** 1 — code archaeology; 2 — semantic decomposition
**Status:** evidence only. **No class has been written, and the boundary is not decided.**

References: [`reference_inventory.md`](reference_inventory.md). Contract:
[`invariant_table.md`](invariant_table.md). What preservation will mean:
[`acceptance_criterion.md`](acceptance_criterion.md).

---

## 1. What the user configures

| | R1 3D GRE | R2 3D MPRAGE |
|---|---|---|
| geometry | `fov` as `[x y z]`, `Nx = Ny = Nz` | `fov`, `n_x`, `n_y`, `n_z` |
| contrast | `TE` (vector allowed), `TR`, `alpha` | `te`, `tr`, `ti`, `t_recovery`, flip |
| readout | `Tread`, `Tpre` | `readout_duration`, oversampling |
| ordering | `iZ` outer, `iY` inner | partition/line loops |
| spoiling | areas in the body | in the body |
| dummies | `Ndummy = 50` | dummy train |

## 2. Which calculations encode MRI physics

1. **Per-axis sample spacing** `deltak = 1./fov` — a *vector*, and the z component uses the slab
   thickness rather than an in-plane FOV. `Δkz = 1/160 mm = 6.25 1/m` at R1's protocol, measured.
2. **The two encode lattices**, both `((0:N-1) - N/2) * dk` — which is exactly
   `PhaseEncode`'s `center_line = matrix // 2` convention, asymmetric by one at an even matrix.
   Measured on R1: 64 distinct indices per axis, each hit 64 times, on-lattice to 1.5e-4 1/m.
3. **The readout prephaser** `-gx.area/2` — with **no** half-sample correction. Measured: at the
   two central samples `kx = -2.50` and `+2.50` 1/m, exactly ±Δkx/2. **This reference has no DC
   sample**, like the PyPulseq TSE demo and unlike the radial one, which adds `-deltak/2`.
4. **TE from the block pulse's effective centre**, not from the end of an RF block:
   `TE - calcDuration(rf) + calcRfCenter(rf) + rf.delay - calcDuration(gxPre) - calcDuration(gx)/2`.
5. **One shared winder window.** `Tpre = 3 ms` is given once and used for `gxPre`, every `gyPre`
   and every `gzPre` — the references fix it by hand rather than solving for the minimum, which is
   precisely the arithmetic `GRE2DTR` owns and calls the winder coupling.

## 3. Which code exists only because Pulseq is flat

- `seq.registerGradEvent` / `registerRfEvent` pre-registration, "not necessary, but accelerates
  the sequence creation by up to a factor of 2" — a performance workaround for a flat event list;
- the pre-built `gyPre(iY)` / `gyReph(iY)` arrays, one trapezoid per line, for the same reason;
- R2's extended-trapezoid ramp sampling, which is a readout concern rather than a 3D one.

None of it is a candidate API concept.

## 4. The central finding: the official corpus excites non-selectively

The architecture review's hypothesis was that GRE3D's z axis needs **slab rephasing plus partition
encoding in one window**, and that solving them jointly is what justifies a separate kernel.

**No *official Pulseq* reference does that, because none of them excites selectively** — and the
one implementation that does, R4b, was found only after this section was first written, which is
why the candidate now has two paths rather than one.

| reference | excitation | slab rephaser |
|---|---|---|
| `writeGradientEcho3D.m` | `mr.makeBlockPulse(8°, 0.2 ms)` | none |
| `writeMPRAGE.m` | `mr.makeBlockPulse(alpha, rfLen)` | none |
| `writeMPRAGE_grappa.m` | `mr.makeBlockPulse(...)` | none |
| `writeMPRAGE_4ge.m` | `mr.makeBlockPulse(...)` | none |
| `write_3Dt1_mprage.py` | `pp.make_block_pulse(...)` | none |
| **`fmrifrey/lps` `gre3d_write_seq.m`** | **`mr.makeSincPulse`, slab = 0.7 x FOV** | **yes, `-gz.area/2`** |

Measured on R1, which settles it independently of reading the source: **`kz` is held at the
partition value right through the readout** — `200.00 → 200.00 → 200.00` 1/m from the first
sample, through the echo, to the last. A slab rephaser would bring `kz` to zero at the echo; a
partition encode holds it. On this evidence z is an **encode**, full stop.

### What that changes

**It splits the candidate into two paths rather than settling it.** Review has since decided that
`GRE3DTR` should support both, as one kernel — see §9 — so the finding is not "z is only ever an
encode" but "the official corpus only demonstrates the encode-only case, and the selective case
has exactly one witness (R4b)".

For the **non-selective** path the coupling is the **same kind** as `GRE2DTR`'s, with a different
participant on z:

```text
GRE2DTR      x readout prephaser | y phase encode | z slice rephaser
GRE3DTR      x readout prephaser | y phase encode | z partition encode
```

One window, three axis minima, the widest wins — not a combined-moment solve. That is still
genuine cross-leaf coupling, and by the standing rule (*a kernel owns arithmetic that requires
information from several leaves*) it still justifies a sibling kernel. But **it justifies it for a
different reason than the one hypothesised**, and the difference matters because it changes what
the kernel has to be able to do.

For the **slab-selective** path the coupling is the one the review predicted, and R4b shows it
exists in practice — see §9 for what the evidence says about how to solve it.

## 5. Semantic decomposition

**Non-selective** (R1, R2, R3):

```text
Excitation (RF only -- no gradient, no rephaser)
      |
      |  x readout prephaser  }
      |  y phase encode       }  one shared winder window
      |  z partition encode   }
      v
   TE delay  ->  [ readout + ADC ]  ->  [ y rewind | z rewind | x spoiler ]  ->  TR fill
```

**Slab-selective** (R4b), with the one structural difference this candidate is about:

```text
Excitation (slab-selective: RF + slab-select gradient, leaving a rephasing requirement)
      |
      |  x readout prephaser      }
      |  y phase encode           }  one shared winder window
      |  z slab rephase           }  <- R4b spends a separate window here...
      |    + partition encode     }  <- ...and another one here
      v
   TE delay  ->  [ readout + ADC ]  ->  [ y rewind | z rewind | x spoiler ]  ->  TR fill
```

| component | class |
|---|---|
| non-selective excitation | `EXISTING_SEQCRAFT` — `Excitation(thickness_mm=None)` |
| readout + prephaser + ADC | `EXISTING_SEQCRAFT` — `CartesianLine` |
| y phase encode + rewind | `EXISTING_SEQCRAFT` — `PhaseEncode(axis='y')` |
| z partition encode + rewind | `EXISTING_SEQCRAFT` — `PhaseEncode(axis='z')`, **axis-generic already** |
| **z slab rephase + partition encode as one moment** | **`CANDIDATE_OWNED`** — selective mode only; see §9 |
| the shared window across three axes | **`CANDIDATE_OWNED`** |
| TE/TR from the pulse's effective centre | **`CANDIDATE_OWNED`** |
| x spoiler | `EXISTING_SEQCRAFT` — `spoiler()`; R1's choice of a full readout area is a parameter |
| partition and line ordering, dummies, RF spoiling schedule | `CALLER_POLICY` — imaging layer |
| event pre-registration, per-line trapezoid arrays | `REFERENCE_ARTIFACT` |

**Every leaf already exists.** What the candidate adds is the coupling and the timing.

## 6. A second finding: the z polarity is a vendor convention, not physics

R1 negates the z lattice:

```matlab
areaY =  ((0:Ny-1)-Ny/2)*deltak(2);
areaZ = -((0:Nz-1)-Nz/2)*deltak(3); % on-line reconstruction on Siemens with autoLabel
                                    % doesn't work otherwise (e.g. without -)...
```

Measured consequence: `ky` runs over indices **−32 … +31** and `kz` over **−31 … +32**, mirrored,
and the `.seq` definitions record `kSpaceCenterLine 31` against `kSpaceCenterPartition 32`.

The sign is chosen to satisfy a vendor reconstruction, which makes it **acquisition policy**. It
must not be baked into `PhaseEncode` or into the kernel: the imaging layer chooses the partition
list, and a caller who needs the reference's convention should be able to express it by passing
the indices they want. This is the same conclusion the earlier pilots reached about `segments` and
about spoke schedules, arrived at from a third direction.

## 7. The decided shape

Settled by architecture review, recorded here so the rest of the scan tests it rather than
rediscovers it:

```text
GRE3DTR

x   Cartesian readout prephaser
y   phase encoding
z   non-selective mode   -> partition encode
    slab-selective mode  -> slab rephase + partition encode, one kernel-owned combined solve

RF  slab_thickness_mm = None  -> non-selective
    slab_thickness_mm = float -> slab-selective
    duration and time-bandwidth remain user-adjustable RF design parameters

policy above the kernel: partition ordering, acquisition ordering, scan-level schedules
```

**One kernel, two excitation contracts.** Not `GRE3DTRNonSelective` and
`GRE3DTRSlabSelective`: the acquisition is the same and only the excitation differs.

**The mode is chosen by a physical parameter, not a flag.** `slab_thickness_mm=None` means
non-selective, exactly as `Excitation(thickness_mm=None)` already does. A boolean beside a
thickness is two arguments that can contradict each other.

**`slab_thickness_mm` is not `fov_z`, and neither bounds the other.** `fov_z` is the nominal
encoded and reconstructed extent; `slab_thickness_mm` is what the RF actually excites. All three
relations are legitimate protocols:

| | |
|---|---|
| `slab < fov_z` | calibration, inner-volume or restricted-support acquisitions — R4b's `slabfrac = 0.7` |
| `slab ≈ fov_z` | a simple full-volume acquisition |
| `slab > fov_z` | slab oversampling, transition-band margin, alias suppression |

So the **only** validation is `slab_thickness_mm > 0` when selective excitation is requested.
`slab_thickness_mm >= fov_z` is **not** a rule and must not be enforced; R4b's 0.7 would fail it,
and R4b is a working sequence. What the documentation should say instead is that a conventional
full-volume GRE3D commonly chooses a slab at least as large as the nominal `fov_z`, and often
somewhat larger, so that the profile's transition region does not alias into the reconstructed
volume.

## 8. Which RF parameters the kernel forwards

Minimum set, in `Excitation`'s own vocabulary rather than a GRE3D dialect:

| | |
|---|---|
| `slab_thickness_mm` | `None` for non-selective; a thickness selects the slab and the mode |
| `rf_duration_s` | as `GRE2DTR` already forwards it |
| `rf_time_bw_product` | **default `None`**, deferring to the RF factory's own default |

`None` rather than a number, and specifically not `4.0`: `Excitation` and the pypulseq factory
below it already own that default, and a second copy here is two defaults that drift. What the
documentation should say instead is that a sharper slab profile is a deliberate choice —
*increase `rf_time_bw_product` explicitly when the transition band matters* — rather than a rule
that 3D imaging implies a particular value.

`pulse` and `pulse_opts` are **not** forwarded in the first extraction. They are added only if the
fine scan finds a reference that needs them; R4b uses a plain apodised sinc.

## 9. The combined-z hypothesis, stated so it can be tested

For the slab-selective path the two requirements land on **one axis inside one timing window**:

```text
A_z(p) = A_slab_rephase + A_partition(p)
```

and `GRE3DTR` designs **one legal z winder** for that total.

Neither leaf can do it: `Excitation` knows the rephasing its own slab selection implies and
nothing about partitions; `PhaseEncode(axis='z')` knows the moment a partition index wants and
nothing about a slab. The kernel is the first layer holding both.

**R4b does it sequentially**, `gzReph` then a scaled `gzPre`, and its own `te_min` pays for both
windows. So the combined solve makes a claim that can be measured rather than argued:

> the same final `kz` state, in one window instead of two, and therefore a shorter minimum TE.

### 9.1 Both terms are signed, and the worst partition is not the largest index

This is the part that is easy to get wrong, and it cannot be recovered from later. With

```text
A_slab_rephase = -120        partitions:  low edge -200 | centre 0 | high edge +200
A_z:                                      low -320      | centre -120 | high  +80
```

the **low** edge is limiting. Reverse the slab gradient's polarity:

```text
A_slab_rephase = +120        A_z:         low  -80      | centre +120 | high +320
```

and the **high** edge is. The limiting partition is a property of the *combined signed* moment,
so nothing may assume "the last partition is worst", or "the largest `|A_partition|` before
combining is worst", or that a larger index means a more positive `kz` — the corpus already
reverses the z lattice relative to y for a vendor reconstruction (§6).

The sign must therefore come from `PhaseEncode`'s own semantics — *partition index → signed
physical moment* — and never from array order.

### 9.2 The solve, in the order that makes it correct

```text
for each supported partition p:
    A_partition(p)   signed, from PhaseEncode semantics
    A_z(p)         = A_slab_rephase + A_partition(p)
    t_z_min(p)     = shortest legal duration for A_z(p) under this scanner's grad/slew limits

z_worst_s        = max over p of t_z_min(p)
shared_winder_s  = max(x_prephaser_min_s, y_encode_worst_s, z_worst_s)
```

then every participating gradient is designed to that one duration.

**Enumerate every partition.** `A_z(p)` is affine in `p`, so for a contiguous interval on
symmetric hardware the extreme is at one of the two edges — but the shortcut buys nothing at
realistic matrix sizes and stops being true for partial `kz`, non-contiguous partition sets or a
supplied partition table. Clarity over premature algebra.

### 9.3 One winder duration for every partition

Not each partition's own minimum. A per-partition window would make **TE depend on `kz`**, which
is a contrast gradient across the volume that no reconstruction expects and nothing downstream
would report. The worst case is found at design time and every partition uses it.

### 9.4 Lengthening the window is design, not legalization

If the worst partition needs more time than the other axes, the kernel **lengthens the shared
window** rather than refusing. That is not compiler work creeping upward:

```text
kernel      given a physical target and this scanner's limits, what legal waveform
            and duration implement it?
compiler    given already-designed events at already-chosen times, how are they
            split, summed and emitted as legal Pulseq blocks?
```

`GRE3DTR` already receives `opts` with `max_grad`, `max_slew` and the raster, so "this moment
needs 420 µs rather than 300 µs" is a design question it is equipped to answer. The compiler must
*not* rescue an infeasible module by stretching a gradient or moving a readout: that would change
TE, and MRI timing semantics are not the compiler's to change.

With `te_s=None` and `tr_s=None` the module therefore returns **the shortest legal design** and
reports what it achieved. With an explicit request below that minimum it **refuses**, naming the
limiting partition and its combined moment — the same rule `GRE2DTR` already applies to TE and TR.

### 9.5 The slab rephaser must not be emitted twice

If the kernel realises `A_slab_rephase + A_partition(p)` on z, then `Excitation` must not also
emit its own standalone rephaser, or the slab term is applied twice.

The precedent is exact: `CartesianLine(prephase=False)` drops the prephaser **event** while
`area_to_echo_per_m` keeps stating the intrinsic **requirement**, which is what lets `TSEShot` own
the compensation. The analogous shape here is for `Excitation` to keep reporting the rephasing its
slab selection implies while not realising it — and today it has no such property: `gz` and `gzr`
are public events, but there is no `rephaser_area_per_m`-style number, and `build` always emits
`gzr` when selective.

**Not settled here.** Whether this is a constructor option mirroring `prephase=False`, a private
mechanism, or something the fine scan finds cleaner, is for the extraction to decide with evidence.
What is settled is the failure it must avoid.

## 10. What external references are for here

R1 and R2 validate the non-selective path's k-space semantics. R4b validates that the selective
ingredients and their k-space consequences are right.

**Neither dictates event decomposition.** R4b plays two z events where SeqCraft may legitimately
play one; what has to agree is the physics — the `kz` reached, the partition spacing, the semantic
centre, TE, TR, rewinding and hardware legality. Requiring an external implementation that happens
to use SeqCraft's construction would be looking for a coding pattern rather than for evidence.

## 11. Still open

1. **How is the intrinsic rephasing requirement exposed without emitting it?** §9.5 — the
   precedent exists but the property does not.
2. **Does `GRE3D` (imaging) follow?** Partition/line ordering, dummies and RF spoiling are policy
   and would belong there, as `GRE2D` holds them today. Not proposed yet.
3. **Is the three-axis window solve now worth a shared private helper?** This would be its third
   instance after `GRE2DTR` and `TSEShot._shared_window_s`. The standing rule was "two is not a
   contract"; three might be.
4. **The z polarity** (§6) stays caller policy. The kernel must not bake in R1's sign.

## 12. Not yet done

- No sweep, and no candidate to sweep.
- R2 measured only near its default protocol, for the slew reason in the inventory.
- R4b not executed; no `.seq` ships with it.
- No comparison run yet: the acceptance criteria are written first, deliberately, in
  [`acceptance_criterion.md`](acceptance_criterion.md).
