# GRE3DTR Fine Scan — Reference Inventory

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** `GRE3DTR`
**Date:** 2026-09-19
**Pass:** 0 — freeze provenance

---

## 1. References

### R1 — Pulseq MATLAB 3D gradient echo *(the cleanest 3D GRE, and it ships a `.seq`)*

```yaml
repo:     pulseq/pulseq
commit:   2fd6ab6   (file last changed accf5f1)
path:     matlab/demoSeq/writeGradientEcho3D.m
license:  MIT
language: matlab
role:     primary
runnable: yes, without MATLAB -- the repository also ships matlab/demoSeq/gre3d.seq
```

Protocol: FOV 200 x 200 x 160 mm, 64³, TE 10 ms, TR 40 ms, read 3.2 ms, `Tpre` 3 ms, rise time
400 µs fixed, 50 dummies, RF spoiling 84°, **non-selective 8° block pulse of 0.2 ms**.

**`gre3d.seq` is the find that makes this candidate executable.** File format 1.5.1, 24 876
blocks, 165.9 s, written by the MATLAB toolbox, and its definitions carry
`kSpaceCenterLine 31`, `kSpaceCenterPartition 32`, `kSpaceCenterSample 31`. Reading it with
`pypulseq.Sequence.read` gives a measurable 3D reference with **no MATLAB run required** — which
matters under the Python-only constraint, where a MATLAB round trip is evidence rather than a
gate.

### R2 — official PyPulseq 3D MPRAGE

```yaml
repo:     imr-framework/pypulseq
commit:   f2c582b
path:     examples/scripts/write_3Dt1_mprage.py
license:  MIT
role:     supporting -- executable, but not a plain 3D GRE
runnable: yes, at or near its default protocol
```

A 3D GRE *train* under an inversion, with ramp-sampled extended-trapezoid readouts. Two
consequences: it is the only executable **Python** 3D reference, and it is a harder comparison
than R1 because the prep and the ramp sampling are not part of the candidate. It also refuses to
build at small `n_x` — `n_x=32` raises a 386 % slew violation inside
`make_extended_trapezoid` — so any sweep over it must keep `n_x` near its default.

### R3 — Pulseq MATLAB MPRAGE family

```yaml
repo:    pulseq/pulseq
commit:  2fd6ab6
paths:   matlab/demoSeq/writeMPRAGE.m, writeMPRAGE_grappa.m, writeMPRAGE_4ge.m
license: MIT
role:    architecture evidence
```

### R4 — SeqCraft's own 2D composites

`examples/mprage_2d/`, `examples/mp2rage_2d/` — `IRPrep` + a 2D GRE train. **Read for what they
are not:** there is no 3D implementation anywhere in SeqCraft, which is what decides this
candidate's acceptance criterion. See [`acceptance_criterion.md`](acceptance_criterion.md).

### R4b — `fmrifrey/lps` 3D GRE *(the slab-selective reference)*

```yaml
repo:     fmrifrey/lps
commit:   552d87a   (2026-08-12)
path:     gre3d/gre3d_write_seq.m
license:  MIT (Copyright 2025 David Frey)
language: matlab
role:     primary for the slab-selective path
runnable: not attempted -- read for physics; no .seq ships with it
```

**The only slab-selective Cartesian 3D GRE found.** It matters because the rest of the corpus
excites non-selectively, so without it the selective path would have no evidence at all.

```matlab
[rf, gz] = mr.makeSincPulse(fa*pi/180, 'Duration', 3e-3, ...
    'SliceThickness', slabfrac*fov*1e-2, 'timeBwProduct', 4, 'use', 'excitation', ...);
gzReph   = mr.makeTrapezoid('z', 'Area', -gz.area/2, ...);
gzPre    = mr.makeTrapezoid('z', 'Area', N/(fov*1e-2)/2, 'Duration', calcDuration(gxPre), ...);
```

Three facts worth recording before any design decision:

1. **The slab is a fraction of the FOV, not the FOV**: `slabfrac = 0.7` by default, so the excited
   slab is *narrower* than the nominal FOV here. Slab thickness and `fov_z` are separate
   quantities in a real implementation, which settles that question with evidence.
2. **It plays the two z moments sequentially** — `gzReph` in one window, then the scaled `gzPre`
   in the next — and **its own TE arithmetic pays for both**:
   `te_min = calcDuration(rf)/2 + calcDuration(gzReph) + calcDuration(gzPre) + calcDuration(gx)/2`.
   A combined solve would remove `calcDuration(gzReph)` from the minimum TE, which makes the
   SeqCraft hypothesis a *measurable* claim rather than a stylistic preference.
3. The partition encode is a **scaled** copy of one designed trapezoid
   (`mr.scaleGrad(gzPre, pescz)`), with `+ (pescz == 0)*eps` to avoid a zero-amplitude gradient —
   a pulseq artefact, not a candidate concept.

### R5 — OpenMRF `GRE`

```yaml
repo:    OpenMRF/openmrf-core-matlab
commit:  fd03604
path:    include_pulseq_toolbox/src_readouts/GRE/
license: MIT for original OpenMRF work
role:    independent implementation -- but 2D
```

Slice-selective 2D, and worth one line of evidence anyway: its `gz_spoil` **fuses** the slice
rephasing into the spoiler, `Area = spoil_nTwist/dz - gz.area/2`. That is a cross-leaf fusion on
z at the *tail* rather than in the winder, and it is the nearest thing in the corpus to the
coupling this candidate was hypothesised to need.

### Independence check

*Added to Pass 0 after the radial scan, where two references turned out to be one design ported.*

R1 and R3 share an author and a house style but are **not** ports of each other; R2 is an
independent Python implementation of the MPRAGE variant; R4b is independent of all of them, by a
different author and institution; R5 is independent and 2D. The 3D
Cartesian encoding arithmetic is nonetheless nearly identical across R1, R2 and R3 — all three
compute `((0:N-1) - N/2) * dk` per axis — so agreement on *that* is weak evidence, and it is
treated as one witness in the invariant table.

---

## 2. Protocol comparison

| | R1 3D GRE | R2 3D MPRAGE | R3 MPRAGE | **R4b lps** |
|---|---|---|---|---|
| excitation | **non-selective block**, 0.2 ms, 8° | **non-selective block** | **non-selective block** | **slab-selective sinc**, 3 ms, TBW 4 |
| slab rephaser on z | none | none | none | **yes**, `-gz.area/2`, its own window |
| slab vs FOV | — | — | — | `slabfrac = 0.7` of the FOV |
| z encoding | `gzPre(areaZ)` … `gzReph(-areaZ)` | partition encode per TR | partition encode per TR |
| readout | 3.2 ms flat, `gxPre = -gx.area/2` | ramp-sampled extended trapezoid | flat |
| DC sample on x | **no** — samples straddle `k = 0` by Δk/2 | — | — |
| x spoiler | `area = gx.area` (a full readout area) | — | — |
| winder window | one `Tpre = 3 ms` for x, y and z | shared | shared |
| RF spoiling | 84°, quadratic | 117°/quadratic | 84° |
| ramp time | fixed 400 µs via `mr.opts('riseTime', …)` | derived | derived |

**Every *official Pulseq* Cartesian 3D reference excites non-selectively**, and the one
slab-selective implementation found is outside that corpus. That split is the central finding of
the scan, and it is why the candidate now has two paths with different evidence behind each;
[`findings.md`](findings.md) §4 works through it.

## 3. Executable path

```text
R1   gre3d.seq -> pypulseq.Sequence.read          measured, no MATLAB
R2   main(...) -> pp.Sequence                     measured, near default protocol only
R3   source only
R4b  source only -- no .seq ships; a one-off MATLAB run would make it measurable
R5   source only
```

R4b is the one place a single MATLAB run would buy real evidence, since it is the only witness to
the selective path. It is **not** a gate: under the Python-only constraint the selective path is
validated by analytic moment checks and metamorphic limits
([`acceptance_criterion.md`](acceptance_criterion.md) §3), with R4b as corroboration.
