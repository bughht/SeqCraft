# RadialReadout Fine Scan — Reference Inventory

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** `RadialReadout`
**Date:** 2026-09-18
**Pass:** 0 — freeze provenance

---

## 1. References

### R1 — official PyPulseq radial GRE *(primary executable reference)*

```yaml
repo:     imr-framework/pypulseq
commit:   f2c582b
path:     examples/scripts/write_radial_gre.py
license:  MIT
language: python
entrypoint: main(...) -- parameterised, returns pp.Sequence
runnable: yes -- runs against the installed pypulseq 1.5.1, unmodified
role:     primary
```

Defaults: FOV 260 mm, `n_x=64`, flip 10 deg, slice 3 mm, TE 8 ms, TR 20 ms, 60 spokes,
20 dummies, angle increment `pi/n_spokes`, RF spoiling 117 deg.

### R2 — Pulseq MATLAB radial GRE

```yaml
repo:    pulseq/pulseq
commit:  2fd6ab6   (file last changed a4455ef)
path:    matlab/demoSeq/writeRadialGradientEcho.m
license: MIT
role:    authoritative external
```

**Finding, recorded before any analysis: R1 and R2 are the same design.** Same
`-gx.area/2 - deltak/2` prephaser, same `rotate('z', phi, ...)` per block, same timing
arithmetic; they differ only in protocol numbers (`Nx` 320 vs 64, `Nr` 256 vs 60, RF spoiling
84 vs 117 deg) and a multi-TE loop. R1 is a port of R2.

The consequence matters for the whole scan: **agreement between R1 and R2 is not independent
evidence**, and the weight therefore falls on R2b and R3 below. This is the opposite of the
TSE pilot, where the two Python witnesses and the MATLAB one were genuinely separate.

### R2b — Pulseq MATLAB *fast* radial GRE *(independent formulation)*

```yaml
repo:    pulseq/pulseq
commit:  2fd6ab6
path:    matlab/demoSeq/writeFastRadialGradientEcho.m
license: MIT
role:    independent-formulation evidence
```

Same family, materially different construction: readout oversampling (`ro_os=2`), a prephaser
written from **sample positions** rather than from half an area, readout spoiling implemented by
*extending the readout's flat time*, and the slice spoiler fused into the selection gradient with
`makeExtendedTrapezoidArea` + `splitGradientAt`.

### R3 — OpenMRF `RAD` *(independent implementation)*

```yaml
repo:    OpenMRF/openmrf-core-matlab
commit:  fd03604
path:    include_pulseq_toolbox/src_readouts/RAD/{RAD_init,RAD_add,RAD_calc_adc}.m
license: MIT for original OpenMRF work (third-party components covered separately)
role:    independent implementation
```

The genuinely independent one, and it does something none of the others does: it **measures**
the k-space error instead of deriving it. `RAD_init` compiles a two-block probe sequence, calls
`calculateKspacePP`, reads `kx` at the centre sample, and subtracts that residual from the
prephaser area. Its optional rewinder is found the same way, from the trajectory's final point.

It also separates the **angle schedule** from the readout design: `phi_mode` in
`{EqualPi, Equal2Pi, RandomPi, Random2Pi, GoldenAngle, RandomGoldenAngle}` produces a list of
angles, and `RAD_add` plays one of them per repetition.

### R4 — official PyPulseq UTE *(the centre-out evidence)*

```yaml
repo:     imr-framework/pypulseq
commit:   f2c582b
path:     examples/scripts/write_ute.py
license:  MIT
runnable: yes
role:     decides the full-spoke vs centre-out question
```

Carries a `readout_asymmetry` argument. At `0.0` it is a full spoke through the centre of
k-space; at `1.0` it is centre-out. One script, one formula, one continuum — see
[`findings.md`](findings.md) §5.

### Deferred

`writeFastRadialGradientEcho_rot3D.m`, `writeRadialGradientEcho_rotExt.m` (rotation extensions),
`writeUTE.m`, `writeUTE_rs.m` (ramp sampling) and OpenMRF `UTE`. Second wave.

---

## 2. Protocol comparison

| | R1 PyPulseq | R2 MATLAB | R2b fast | R3 OpenMRF | R4 UTE |
|---|---|---|---|---|---|
| max_grad / max_slew | 28 / 120 | 28 / 120 | 28 / 120 | config | 28 / 100 |
| FOV | 260 mm | 260 mm | 240 mm | config | 250 mm |
| samples | 64 | 320 | 240 x 2 | `Nxy x ro_os` | 64 x 2 |
| spokes | 60 | 256 | 256 | config | 32 |
| angle increment | `pi/n` | `pi/Nr` | `pi/Nr` | six modes incl. golden | `pi/n` |
| readout | 1.28 ms flat | 1.28 ms flat | 1.2 ms + spoil extension | `adcTime` | 2.56 ms flat |
| prephaser | `-gx.area/2 - dk/2` | same | sample-centre formula | half areas **+ measured correction** | asymmetry-dependent |
| RF spoiling | 117 deg | 84 deg | 84 deg | configurable, lin or quad | 117 deg |
| rotation | `pp.rotate` per block | `mr.rotate` per block | `mr.rotate` | `mr.rotate` | none — UTE is `rotate`-free, angles are in `gx`/`gy` amplitudes |

R4 is the exception in the last row and it is worth noting early: it builds each spoke by
**scaling a designed gradient onto two axes** rather than by calling a rotate helper. Pulseq's
`rotate` is a block-level helper with no SeqCraft equivalent, and the compiler is out of scope for
this project, so R4's approach is the one a SeqCraft module would have to take.

## 3. Executable path

```text
R1, R4  main(...) -> pp.Sequence              in process, no edits
R2, R2b MATLAB -> .seq -> Sequence.read       deferred; evidence only, not an implementation target
R3      MATLAB -> .seq -> Sequence.read       deferred; read for architecture, not executed
```

Under the project constraint recorded in the plan (`package modules are Python-only for now`),
R1 and R4 are the executable references that matter, and R2/R2b/R3 are physics and architecture
evidence. No MATLAB run has been made for this candidate.
