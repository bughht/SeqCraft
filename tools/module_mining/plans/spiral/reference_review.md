# `writeSpiral.m` — the independent witness, read on its own terms

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-20, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

`pulseq/pulseq` @ `2fd6ab6`, `matlab/demoSeq/writeSpiral.m`, 293 lines, MIT. Read without
reference to PR23, then compared.

> Agreement between PR23 and a new implementation is not independent evidence. Agreement with
> this is.

---

## 1. What it does

**Path.** A classical Archimedean spiral, `K(t) = φ(t)/(2π·FOV)` — radius proportional to azimuth,
built in polar coordinates with **no scanner limits involved**:

```text
deltak   = 1/fov
kRadius  = Nx/2                                  turns
kSamples = round(2π·kRadius)·adcOversampling     samples on the outermost turn
r = deltak·c·interleaves/kSamples/tos            a = c·2π/kSamples/tos
```

**Traversal.** A separate pass. `traj2grad` differentiates the path, then the time axis is
re-parameterised so neither limit is exceeded:

```text
dt_gabs = |g|/(maxGrad·0.97)·du
dt_sabs = sqrt(|s|/(maxSlew·0.97))·du
dt_opt  = max(dt_gabs, dt_sabs)          then interpolate the path onto the gradient raster
```

**Start condition.** A `slowStarting` reparameterisation, `c − K·log(1+c/K)`, plus a half-raster
safety factor on the first step — because "the extreme non-smoothness of the trajectory start"
would otherwise exceed slew.

**End condition.** `'first', 0, 'last', spiral_grad_shape(:,end)` — it **starts at zero and ends
at full gradient**. The ramp-down lives in an extended-trapezoid spoiler that begins at the
readout's final amplitude.

**ADC.** `mr.calcAdcSeg(adcSamplesDesired, adcDwell, sys)` against `adcSamplesLimit = 8192`, and
`MaxAdcSegmentLength` is written into the file as a definition so the scanner runs it without
tweaking. The ADC is advanced by `adcDwell/2` so **the first sample is k ≈ 0**.

**Interleaves.** `interleaves` enters the **path** (`r ∝ c·interleaves/kSamples`). Rotation is
applied afterwards by `mr.TransformFOV` over a block range, plus a per-block `mr.rotate('z', phi)`
for the base orientation.

**Spoiling.** z spoiler of `4·Nx·deltak`; x and y ramp from the readout's last amplitude to zero.

## 2. What it does **not** do

Searched: zero occurrences of spiral-in, in-out, out-in, reverse, flyback or multi-echo.

```text
variable density     NO -- FOV is one scalar
return to origin     NO -- there is no rewinder; k is left at the edge
spiral-in            NO
in-out / out-in      NO
multi-echo           NO
flybacks / joins     NO
```

## 3. The classification asked for

| concern | verdict |
|---|---|
| path and traversal as **separate passes** | **reference-supported physical requirement.** Independently arrived at, and for the same reason |
| Archimedean path, `r ∝ φ` | reference-supported |
| slew- and amplitude-limited time reparameterisation | reference-supported |
| `interleaves` changes the **path** | reference-supported — confirms the coupling analysis in `findings.md` §1 |
| starting at zero gradient | reference-supported (by a different mechanism: a slow-start on the *path*, not a boundary condition on the traversal) |
| ADC segmentation against a sample limit | **reference-supported**, and Pulseq ships `calcAdcSeg` plus a `MaxAdcSegmentLength` file convention |
| k = 0 is the first sample of a spiral-out | reference-supported |
| **ending at zero gradient** | **PR23-only.** The reference ends at full gradient |
| **rewinder / return to origin** | **PR23-only** |
| **variable density** | **PR23-only** (its lineage is `vds.m`, a third source) |
| **in / in-out / out-in** | **PR23-only** |
| **multi-echo, flybacks, joins** | **PR23-only** |
| rotation as a **per-call** argument | SeqCraft architectural choice. The reference rotates at sequence level; `RadialReadout` already settled per-call here, on its own evidence |
| join correction on the assembled waveform | PR23-only, and a consequence of its own raster-centre representation |

## 4. What this changes

**PR23's failure A is confirmed by the reference, by exhibiting it.** `writeSpiral.m` ends at full
gradient, and it is spiral-out only. Those two facts are the same fact: an arm that ends at full
gradient cannot be time-reversed, so there is no spiral-in to have. PR23's `v = 0 at both ends` is
a **design decision beyond the reference**, and it is precisely what buys the variants.

**Roughly half of PR23's design has no independent witness**, and it is the half that drives the
open `variant`/`echoes` question. An acceptance claim measured against `writeSpiral.m` can cover
the path, the traversal, the ADC segmentation and a spiral-out trajectory — and **nothing else**.

**`density` has no witness here at all.** Two independent sources would be needed before treating
sampled-FOV-multipliers as anything but one implementation's convention, and PR23's own failure B
shows the convention question is where the silent errors live.
