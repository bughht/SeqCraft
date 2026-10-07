# Examples

| | What it covers |
|---|---|
| [`01_getting_started.ipynb`](01_getting_started.ipynb) | Blocks, `Opts` and `compile`; the overlap rules, provenance, checking against physics, writing a file — and `sc.Module` at the end, once there is a reason for one. **Uses no modules**, which is the point: it is the demonstration that the compile path stands alone. |
| [`gre_2d/`](gre_2d/) | A complete spoiled 2D gradient echo, built three ways and then simulated and reconstructed. The best place to start for what a module is and what it composes. |
| [`mprage_2d/`](mprage_2d/) | An inversion-prepared segmented GRE, and where an inversion time gets placed to the microsecond. Introduces `IRPrep` and `GRE2D.time_to_center_line`; defines `MPRAGE2D` in its own notebook. |
| [`mp2rage_2d/`](mp2rage_2d/) | One inversion, two trains, and the `SET` label that separates them — plus the ratio that cancels the receive field. Defines `MP2RAGE2D` in its own notebook. |
| [`se_2d/`](se_2d/) | A spin echo, and the one rule that makes it one: between refocusing centres, area before the echo equals area after it. Introduces `Refocusing` and `CartesianLine(prephase=False)`; defines `SE2D` in its own notebook. |
| [`fse_2d/`](fse_2d/) | The same composition at 1, 16 and 72 echoes — TSE and HASTE — and what an echo train costs in blurring, ghosting and signal. Defines `FSE2D` in its own notebook. |
| [`gre_epi_2d/`](gre_epi_2d/) | The whole of k-space in one shot: the centred sampling window, the blip on the zero crossing, ramp sampling and the operator that undoes it, off-resonance and the N/2 ghost, and GRAPPA. Introduces `EPI2D`; defines `GREEPI2D` in its own notebook. |
| [`se_epi_2d/`](se_epi_2d/) | The same readout after a refocusing pulse — and the measurement that a spin echo **does not** fix EPI distortion. Defines `SEEPI2D` in its own notebook. |
| [`megre_2d/`](megre_2d/) | Eight echoes off one excitation, monopolar and bipolar, fitted for T2\* and ΔB0 against a phantom that carries the ground truth for both — the one directory here whose output is a *map* rather than an image. Plus what it takes to flow compensate every echo of a train, and why the two readouts pay very differently for it. |
| [`gre_radial_2d/`](gre_radial_2d/) | One spoke and an angle schedule you write yourself — equal increments and a golden angle out of the same `RadialReadout`. Its `02` is the introduction to non-Cartesian reconstruction. |
| [`fat_sat/`](fat_sat/) | A spectrally selective saturation and the spoiler that destroys what it made — the chemical-shift sign traced from ppm to the emitted `freq_offset`, and the same protocol across four field strengths. Introduces `SaturationPrep`. |
| [`gre_spiral_2d/`](gre_spiral_2d/) | A spoiled gradient echo on a spiral: one arm, eight interleaves, and the **TE that this layer owns** — measured from the excitation's effective centre to the crossing the readout reports. Then what a 21.7 ms readout pays in off-resonance. |
| [`se_spiral_2d/`](se_spiral_2d/) | A spin echo in front of a spiral, and the alignment it needs: the trajectory's origin crossing has to land on the physical echo. Three placements that all compile and only one of which is right, and then what refocusing is actually worth. |
| [`dwi_se_epi_2d/`](dwi_se_epi_2d/) | Diffusion-weighted imaging: a Stejskal–Tanner gradient pair around a 180°, a single-shot EPI readout, and an ADC map recovered from a phantom whose diffusion coefficient is known. Introduces `DiffusionSEPrep` and `sc.b_value`. |
| [`t2prep_gre_2d/`](t2prep_gre_2d/) | A T2 preparation in front of a spoiled gradient echo: adding controlled T2 weighting before a readout that would otherwise follow T2\*, and the measurement that the weighting really is T2 rather than T2\*. Introduces `T2Prep`. |
| [`pc_gre_2d/`](pc_gre_2d/) | Phase contrast: a bipolar pair that leaves still spins alone and gives moving ones a phase proportional to their velocity, and the two acquisitions whose difference is a velocity map. Introduces `VelocityEncode`. |

## `gre_2d/`

| | |
|---|---|
| [`01_build.ipynb`](gre_2d/01_build.ipynb) | The sequence, three times over: raw pypulseq events, the four leaf modules composed inline, and the same composition written as a module. Every arithmetic check, three sampling patterns, three `.seq` files. **Needs nothing but `seqcraft`.** |
| [`02_simulate_and_reconstruct.ipynb`](gre_2d/02_simulate_and_reconstruct.ipynb) | Those three files against a BrainWeb phantom — PD, T1, T2, T2′, D and a synthesised B0, all six simulated — with an eight-element receive ring, then one reconstruction across three samplings. **Needs `seqcraft[sim,recon]`** and a one-off ~19 MB phantom download. |
| [`03_flow_comp.ipynb`](gre_2d/03_flow_comp.ipynb) | Why a gradient echo's first moment at the echo gives moving spins a velocity-dependent phase, and the reshaped pre-echo gradients that remove it — on `x`, `y` and `z`, with the zeroth-moment job each axis already had. Why the reference instant matters on the phase-encode axis. Two `.seq` files, the compensated one on all three axes. **Needs nothing but `seqcraft`.** |

**On where flow compensation lives.** `03` is a later notebook in this family rather than a
directory of its own, because flow compensation is a reusable physical capability of a gradient
echo, not a different acquisition. A new top-level example directory is for a distinct sequence
family; a capability added to one belongs beside it.

`01` deliberately does **not** import `GRE2DTR` or `GRE2D`. It writes them, because the third pass
is where a reader sees a working composition become something reusable, and a pass that imported
the answer would teach nothing.
[`tests/modules/test_notebook_matches_the_package.py`](../tests/modules/test_notebook_matches_the_package.py)
asserts that what the notebook writes and what the package ships produce identical events, so the
tutorial cannot drift from the library without CI noticing.

## `mprage_2d/` and `mp2rage_2d/`

| | |
|---|---|
| [`mprage_2d/01_build.ipynb`](mprage_2d/01_build.ipynb) | The timing, which is where the risk is: TI runs from `IRPrep.time_to_center()` to `GRE2D.time_to_center_line()`, and both ends are module methods because a tree of events cannot know either. Segmentation, dummy shots, the two minima that stop a train overlapping its own inversion, then `MPRAGE2D` written out. Two `.seq` files, one per line ordering. **Needs nothing but `seqcraft`.** |
| [`mprage_2d/02_simulate_and_reconstruct.ipynb`](mprage_2d/02_simulate_and_reconstruct.ipynb) | The null point at `TI = T1·ln2`, measured through a real shot — the one check arithmetic cannot make — and what the line ordering costs. **Needs `seqcraft[sim,recon]`.** |
| [`mp2rage_2d/01_build.ipynb`](mp2rage_2d/01_build.ipynb) | Two trains at two inversion times, and the stateful-label trap: `SET` has to be emitted before *every* train, because setting it once leaves every later shot's first train wearing the previous shot's value. The broken version is built, and the compiler refuses it. **Needs nothing but `seqcraft`.** |
| [`mp2rage_2d/02_simulate_and_reconstruct.ipynb`](mp2rage_2d/02_simulate_and_reconstruct.ipynb) | Split by `SET` read back out of the file, then `UNI = Re(S₁·conj(S₂))/(|S₁|²+|S₂|²)` — simulated against two different receive arrays to show that the magnitudes move and `UNI` does not. **Needs `seqcraft[sim]`.** |

`MPRAGE2D` and `MP2RAGE2D` are **defined in those notebooks and do not ship**: one consumer each,
and a module with one consumer belongs where that consumer is.
[`tests/modules/test_mprage_notebooks.py`](../tests/modules/test_mprage_notebooks.py) runs the two
build notebooks and asserts against what they defined, so a tutorial that stays a tutorial still
cannot drift silently.

## `gre_3d/`

| | |
|---|---|
| [`gre_3d/01_build.ipynb`](gre_3d/01_build.ipynb) | A complete 3D acquisition out of `sc.modules.GRE3DTR` and **two `for` loops** — non-selective and slab-selective, the signed `kz` coupling printed partition by partition, and TE measured constant across the volume. **Needs nothing but `seqcraft`.** |
| [`gre_3d/02_simulate_and_reconstruct.ipynb`](gre_3d/02_simulate_and_reconstruct.ipynb) | The volume, reconstructed with a 3D FFT and shown in three planes — which is where a reversed `kz` or a `ky`/`kz` swap stops being a number. Also compares the kernel's **combined** z winder against a **sequential** slab-rephase-then-encode arrangement: same image to 0.0014 of the peak, 100 µs less TE per repetition. **Needs `seqcraft[sim]`**; runs in about a minute. |

A complete 3D acquisition is `GRE3DTR` plus two loops, and `01` writes the loops out so that the
ordering stays the caller's.

The four `.seq` files under `gre_3d/seq/` are **two pairs answering two questions**, not four
sequence variants:

| | question | files |
|---|---|---|
| `01` | which **excitation mode**? | `gre_3d_nonselective.seq`, `gre_3d_slab.seq` |
| `02` | within the slab-selective mode, how is the **z moment realised**? | `gre_3d_slab_combined.seq`, `gre_3d_slab_sequential.seq` |

Both files in the second pair are slab-selective, and `gre_3d_slab_combined.seq` is the same
physics as `01`'s `gre_3d_slab.seq` — the production combined-z path. It is written separately
only so the comparison has both halves at one common TR.

## `bssfp_2d/`

| | |
|---|---|
| [`bssfp_2d/01_build.ipynb`](bssfp_2d/01_build.ipynb) | A complete 2D balanced-SSFP acquisition out of `sc.modules.bSSFP2DTR`: zero net gradient area on every axis between consecutive RF effective centres, symmetric timing at `TE ≈ TR/2`, and the **symmetric gradient realisation** behind it — `M1z` zero over the interval, `M1x` constant across the table, `M1y` line-dependent, which is the structure the bSSFP literature assumes and which `M0 = 0` alone does not give. Then a 20-step start-up ramp and an interrupted four-segment acquisition as ordinary loops. **Needs nothing but `seqcraft`.** |
| [`bssfp_2d/02_simulate_and_reconstruct.ipynb`](bssfp_2d/02_simulate_and_reconstruct.ipynb) | Where a balanced *waveform* stops and a *steady state* begins. Single voxels with every repetition sampled, comparing no preparation against `α/2 – TR/2` and a 20-step flip-angle ramp, aligned at the first full-flip repetition; the off-resonance response and its nulls near ±1/(2 TR); then the continuous and interrupted-segmented acquisitions reconstructed and compared line by line on a single-tissue object, with the measured transient compared against Scheffler's analytic decay rate. **Needs `seqcraft[sim]`**; about 4 minutes. |
| [`bssfp_2d/03_flow_and_motion.ipynb`](bssfp_2d/03_flow_and_motion.ipynb) | **Why conventional balanced SSFP often tolerates ordinary flow, and where that breaks down.** Built on one progression: a single spin tells you the *phase* a gradient moment produces, a voxel of spins tells you whether that phase costs *magnitude*, and a train of lines tells you whether the image is *consistent*. A common velocity can rotate voxel phase without reducing magnitude, while intravoxel phase dispersion can cause cancellation. Then Bieri's per-TR increment and the view-order dependence, the failure modes including outflow, what this means for bSSFP angiography — bright blood from `T2/T1` rather than inflow — and only then explicit compensation and which endpoint condition each kind satisfies. **Needs `seqcraft[sim]`**; about half a minute. |

`01` teaches the repetition: what balance is, where the timing origin is, and how start-up and
segmentation are composed around it. `02` is about the magnetisation the train drives, and `03`
about what a *moving* spin sees — three questions that are easier apart than together.

`bssfp_2d/seq/` holds the two acquisitions `01` writes. The probes `02` and `03` build — single
lines read hundreds of times, a preparation pulse sampled on its own, one view order against
another — go in `bssfp_2d/seq/diagnostics/` and are named for what they measure. **None of them
is a flow-compensated bSSFP sequence**: the only compensated variant among them nulls the first
moment at the echo, which is a different condition from the between-excitations one `03` Part C
describes and which is not implemented.

## `bssfp_3d/`

| | |
|---|---|
| [`bssfp_3d/01_build.ipynb`](bssfp_3d/01_build.ipynb) | A 3D balanced-SSFP acquisition out of `sc.modules.bSSFP3DTR`, where `z` selects the slab, balances that selection **and** encodes `kz`. The physics point is *ownership*: each repetition discharges its own partition moment, so balance holds across a view order whose neighbouring repetitions change both indices. The acquisition point is everything above the repetition — a 24-step flip ramp run **at the first acquired encoding**, the `0/pi` phase cycle, and a bidirectional sequential view order whose smoothness is measured over the **whole RF train**, not just the acquired views. Segmentation is defined by a fixed **53-view shot** -- the echo train length Morita et al. (2008) report -- and **comes with its catalyzation**: the same ramp runs before *every* segment, including the first, because an interruption leaves the magnetisation off the periodic state and re-establishing it is part of what the acquisition is rather than a repair applied to it. The continuous acquisition is the same policy code with one segment. Writes a representative `192 x 192 x 16` acquisition and a reduced surrogate at the same TR and shot length, each continuous and segmented, then re-reads all four and asserts the artifact contract `02` depends on. **Needs nothing but `seqcraft`.** |
| [`bssfp_3d/02_simulate_and_reconstruct.ipynb`](bssfp_3d/02_simulate_and_reconstruct.ipynb) | **What a waveform check cannot establish.** Opens with the quantity that decides how long anything takes: the transient decays at an eigenvalue of the repetition's transfer matrix (Hargreaves 2001), computed rather than guessed — so CSF is then excluded from settling claims, its eigenvalue implying ~2500 repetitions. Chooses the preparation by measuring `alpha/2 - TR/2` against a flip ramp **across the off-resonance band**, where the non-selective half-angle pulse loses, and under `+-10%` B1, where it does not; the ramp length comes from a sweep kept in the notebook. Reconstructs the reduced volume in three planes for geometry — against a phantom that covers the encoded `z` extent, with FOV, slab, extent, matrix and voxel printed and asserted against each other. Then segmentation, on the **representative** files, as two deliberately separate measurements: the **per-view echo-magnitude departure** (a 0.1 mm point at the origin has flat k-space by construction, so `|echo|` per repetition is the acquisition's own modulation -- one sampled echo's amplitude, explicitly *not* the magnetisation vector) and the reconstructed point response, which depends additionally on where in the table the modulated views fall and on phase. Prices the catalyzation: 58 shots, 57 breaks, 45% of the acquired views spent re-establishing the periodic state -- time inside a fixed gating window. Keeps one sensitivity demonstration, 53 against 64 views per shot, because a shot length that divides the 192-line table repeats its modulation at the same few `ky`; it is a reason to measure a segmented protocol, not a knob to turn. Ends with a differential off-resonance sweep, a tissue comparison that sweeps the coherence-state budget before reading itself -- white matter converges at 200 states, CSF is still moving at 1600 and is excluded rather than quoted -- and the set of quantities none of these numbers transfer across. **Needs `seqcraft[sim]`**; about 17 minutes, or 11 with `BRAIN_IMAGES = False`.

## `se_2d/` and `fse_2d/`

| | |
|---|---|
| [`se_2d/01_build.ipynb`](se_2d/01_build.ipynb) | One interval placed by hand out of the shipped modules, then the whole shot, then `SE2D`. Every check is `sc.kspace` at the echo *sample*, signed, because `|k|` is symmetric and a k-space extent check passes on a mirrored image. **Needs nothing but `seqcraft`.** |
| [`se_2d/02_simulate_and_reconstruct.ipynb`](se_2d/02_simulate_and_reconstruct.ipynb) | The one claim arithmetic cannot make: the spin echo recovers the phantom's **T2** and the same sweep without the 180 recovers **T2\***. Then the image, and which weighting the protocol landed on. **Needs `seqcraft[sim,recon]`.** |
| [`fse_2d/01_build.ipynb`](fse_2d/01_build.ipynb) | `FSE2D`, turbo 1 → 16 → 72 on one instance, three orderings as data, the echo-band warning, and HASTE with `partial_fourier`. Four `.seq` files. **Needs nothing but `seqcraft`.** |
| [`fse_2d/02_simulate_and_reconstruct.ipynb`](fse_2d/02_simulate_and_reconstruct.ipynb) | The echo envelope measured, the point-spread width per ordering, contrast and blurring as two separate knobs, and **the ghost a scattered table makes** — a periodic modulation of `ky`, which is a replica of the object rather than a blur. Then HASTE with POCS. **Everything is measured on one spin**: a ghost is a modulation of `ky`, and a point object's k-space *is* that modulation — 1.3 s per sequence instead of 173 s, and no reconstruction in between. **Needs `seqcraft[sim,recon]`**; runs in 16 s. It closes with one representative image built from `sc.modules.FSE2D`, so the **package path** is shown to reach a valid end-to-end acquisition — the study itself is not repeated, because event identity against `01` is stronger evidence that the two agree. |

| [`fse_2d/03_module_api.ipynb`](fse_2d/03_module_api.ipynb) | **How to build an FSE today**: `sc.modules.TSEShot` for one shot, `sc.modules.FSE2D` for the scan, three orderings as data, and HASTE as a configuration. Start here if you want to write one rather than understand one. **Needs nothing but `seqcraft`.** |

**`sc.modules.TSEShot` and `sc.modules.FSE2D` ship**, and `03_module_api.ipynb` is the
recommended way to use them. `SE2D` is still **defined in its notebook and does not ship**, for
the reason `MPRAGE2D` and `MP2RAGE2D` do not: one consumer each.

`fse_2d/01_build.ipynb` keeps its **own** `FSE2D` on purpose and should not be replaced with an
import. It is the working implementation the package classes were extracted from, and
[`tests/modules/test_fse_notebook_matches_the_package.py`](../tests/modules/test_fse_notebook_matches_the_package.py)
compares the two event for event at turbo 1, 8, 16 and HASTE. Editing it: see
[*Preserved reference notebooks*](../docs/writing_examples.md).

[`tests/modules/test_se_notebooks.py`](../tests/modules/test_se_notebooks.py) runs both build
notebooks, asserts k at every echo of every train length, and pins `FSE2D(echoes=1)` against what
`se_2d/01` writes — event for event, which is what makes two example directories safe.

## `gre_epi_2d/` and `se_epi_2d/`

| | |
|---|---|
| [`gre_epi_2d/01_build.ipynb`](gre_epi_2d/01_build.ipynb) | The two rules an EPI is: the sampling window **exactly centred** in its lobe, and the blip **centred on the seam** with `sc.barrier()` forced there. Both are built the wrong way as well as the right way, because both failures compile, pass every k-space check and simulate correctly. Then what ramp sampling, partial echo, oversampling and `blip_lines` each cost, `GREEPI2D`, and parallel imaging — including why the calibration data is a **Cartesian reference** rather than more EPI shots, and what a *segmented* EPI needs that a single-shot one does not — with **one echo time, one flip and one TR across every imaging file**, checked and printed, so `02`'s pictures compare acquisitions rather than contrasts. Nine `.seq` files. **Needs nothing but `seqcraft`.** |
| [`gre_epi_2d/02_simulate_and_reconstruct.ipynb`](gre_epi_2d/02_simulate_and_reconstruct.ipynb) | The regridding operator and, before it, the **sampling** question no interpolator can answer; the two reconstruction failures; off-resonance measured against a `B0 = 0` control; the navigator correction and its Fourier dual that does nothing; interleaved against blocked at the same lines and the same shots; a brain under the phantom's own field map; GRAPPA at R = 3 and 4, with the separately acquired calibration band priced against a self-calibration no scan can have; and the transient a four-shot EPI is acquired on, decomposed over three files that differ only in the magnetisation each shot sees — as first built it put a 9.6× replica at `Ny/4` that looks exactly like an under-sampling artefact. **Needs `seqcraft[sim,recon]`** and `pygrappa`. |
| [`se_epi_2d/01_build.ipynb`](se_epi_2d/01_build.ipynb) | `EPI2D` unchanged, after a refocusing pulse: `prephase=False` and both dephasers **before** the 180, with the sign the conjugation will flip. TE as **two numbers** — the gradient echo on the ADC raster, the spin echo on the gradient raster — and why the pulse is placed with `Raster.nearest`. Four `.seq` files, **all four at one echo time**, checked and printed: a single shot, a partial echo, and the same 128 lines in two and in four interleaved shots. **Needs nothing but `seqcraft`.** |
| [`se_epi_2d/02_simulate_and_reconstruct.ipynb`](se_epi_2d/02_simulate_and_reconstruct.ipynb) | That it is a **T2** echo, and then the section the notebook exists for: **a spin echo does not fix EPI distortion.** Same shift, same prediction, same ghost as the gradient echo — the pulse refocuses at one instant and the distortion accumulates *between* echoes. What it does fix, measured against a flat-envelope floor rather than assumed. Then partial echo with POCS, and **segmentation**: the same lines in one, two and four interleaved shots, where the train falls 89.6 → 46.1 → 23.7 ms and the residual against the shimmed image falls 0.127 → 0.073 → 0.047 with it — the same short train an `R = 4` acceleration buys, with no calibration band, no kernel and no unfolding residual. **Needs `seqcraft[sim,recon]`.** |

`EPI2D` is the one module that ships out of this pair; `GREEPI2D` and `SEEPI2D` are **defined in
those notebooks and do not ship**, one consumer each.
[`tests/modules/test_epi_notebooks.py`](../tests/modules/test_epi_notebooks.py) runs both build
notebooks and asserts the things whose failure is silent: every readout gradient still a `trap`,
signed `k` at the `k = 0` sample of every echo, the reverse echoes landing on the forward grid, and
the spin echo where the kernel says it is.

Both `01`s are in [`tools/run_notebook_smoke.py`](../tools/run_notebook_smoke.py); the `02`s are
not, for the reason the other simulation notebooks are not.

## `megre_2d/`

| | |
|---|---|
| [`01_build.ipynb`](megre_2d/01_build.ipynb) | The four numbers that are wrong *quietly*: the reverse lobe's grid (**0.2051 Δk** on the lobe that is right for a single line, built the wrong way first and plotted so a reader sees that the two trajectories *look* identical); the fly-back's area, with all three plausible wrong answers built and their trajectories drawn; the dwell quantum and why it is the RF raster; and `te_s` against `echo_spacing_s`, whose consequence — **0.75 % low everywhere** in a two-point field map — is computed from `te_s` alone, with no simulator. Plus the barrier table, which reports a prediction that turned out **false**. Three `.seq` files, **one flip, one TR, one matrix and one slice across all of them**, checked and printed. **Needs nothing but `seqcraft`.** |
| [`02_simulate_and_reconstruct.ipynb`](megre_2d/02_simulate_and_reconstruct.ipynb) | **One reconstruction for all three files** — ESPIRiT coil maps estimated once (a coil sensitivity belongs to the coil, and using one set is what makes echo *e*'s phase comparable with echo 0's), POCS for the partial echo, then SENSE, because RSS throws away the phase that is half of this notebook. All eight echoes shown in magnitude *and* phase before anything is fitted, then T2\* and ΔB0 fitted for **all three readouts** against the phantom's own `T2dash_map` and `field_hz()` — maps, residuals and scatter side by side, not just the best one. **The phantom is 192 mm across and the protocol images 220 mm**, so the oracle is resampled onto the reconstruction's pixel grid before anything is subtracted, by nearest neighbour so it stays a set of values the phantom contains; the alignment is checked with a Dice overlap rather than assumed. Finally what partial echo costs and what POCS gets back, measured against the full-echo image. **Needs `seqcraft[sim,recon]`.** |
| [`03_flow_comp.ipynb`](megre_2d/03_flow_comp.ipynb) | What "flow compensated" means when there are four echoes rather than one: **every** acquired echo, measured at the arithmetic floor on all three axes of both polarities. `y` and `z` come free — they are silent across the train, and a first moment about a fixed origin stops accruing where the gradient is zero. `x` does not, and the inter-echo recurrence $M_1^{(i+1)} = M_1^{(i)} + \mu_i$ says what it needs: every interval nulling $\Delta M_0$ **and** $\Delta M_1$, a condition on the interval as a whole rather than on the fly-back alone. Monopolar reshapes a fly-back it already had and pays **3 %** of its echo spacing; bipolar has to create one between abutting lobes and pays **59 %**, which **reverses** the uncompensated ranking of the two readouts. Minimum spacing derived from the waveform, so it moves with the amplifier across three hardware profiles. Two `.seq` files. **Needs nothing but `seqcraft`.** |

**This is the one example directory that defines no class**, and that is the result rather than an
omission. `MPRAGE2D`, `SE2D`, `FSE2D`, `GREEPI2D` and `SEEPI2D` each exist because their sequence
is a composition the package does not have. A multi-echo gradient echo is not: it is a spoiled GRE
whose readout was told to read the line more than once, so the whole protocol is

```python
sc.modules.GRE2D(opts=opts, fov_mm=220.0, matrix=(128, 128), thickness_mm=3.0, flip_deg=15.0,
                 bandwidth_hz_px=500.0, tr_s=40e-3, echoes=8, polarity='monopolar')
```

Inventing a `MEGRE2D` to say that would be exactly the wrapper the rule above warns about — one
whose extraction shortens no notebook. `CartesianLine` gained three arguments and `GRE2DTR` and
`GRE2D` gained three pass-throughs and **no arithmetic**; that last part is the design's own
falsification test, and `tests/modules/test_gre_2d_tr.py` is where it would have failed.

[`tests/modules/test_megre_notebooks.py`](../tests/modules/test_megre_notebooks.py) runs `01` and
asserts what fails silently: every readout gradient still a `trap`, signed `k` at the `k = 0`
sample of every echo, the reverse echoes landing on the forward grid, `ECO` and `REV` read off the
*compiled* sequence, and the written `te_s` matching the module's to a nanosecond — which is the
one that matters most, because `02` fits against those numbers and a wrong echo time is a scaled
map with nothing to look wrong. `01` is in
[`tools/run_notebook_smoke.py`](../tools/run_notebook_smoke.py); `02` is not, for the reason the
other simulation notebooks are not.

## `gre_radial_2d/`

| | |
|---|---|
| [`01_build.ipynb`](gre_radial_2d/01_build.ipynb) | `sc.modules.RadialReadout` is **one spoke**, already oriented, and an acquisition is that spoke plus a list of angles — so the repetition is assembled in the notebook out of `Excitation`, the spoke, a spoiler and TR fill rather than by a class. Equal increments and a golden angle from one readout instance, the trajectory drawn for both, the centre sample measured **exactly on `k = 0`** and the spoke angles measured against the ones asked for, then `partial_fourier` walked from a full spoke to centre-out. Two `.seq` files. **Needs nothing but `seqcraft`.** |

**This directory defines no class, and that is the finding.** `MPRAGE2D`, `SE2D`, `FSE2D`,
`GREEPI2D` and `SEEPI2D` each exist because their sequence is a composition the package does not
otherwise have. A radial GRE is not: strip the angle list out of it and what remains is a spoiled
gradient echo whose readout happens to be a spoke. A `RadialGRE` class would own the schedule —
golden angle or equal increment, how many spokes, in what order — and a schedule is a sequence
choice, which is exactly the thing `sc.modules` does not take from the caller. So the notebook
writes the loop, and it is four lines.

| [`02_simulate_and_reconstruct.ipynb`](gre_radial_2d/02_simulate_and_reconstruct.ipynb) | **The introduction to the non-Cartesian path.** Why `ifft2` is not the reconstruction for these samples — 4096 samples for a 4096-point grid with a third of the grid empty — where the physical trajectory enters, what density compensation buys and what it is assumed to be and is not, why the spokes are one inverse problem, and a point at four known positions that fixes the image's axis convention. **Needs `seqcraft[recon]`, `MRzeroCore` and a phantom download.** |

**The reconstruction lives in [`noncartesian_recon.py`](noncartesian_recon.py)**, beside
[`phantom.py`](phantom.py) and for the same reason: it is example infrastructure, not API.
SeqCraft builds sequences, MRzeroCore simulates them, SigPy reconstructs them, and a
reconstruction dependency has no business on the compile path. Its adapter — physical k times the
field of view — is checked against an explicit dense DFT that calls no SigPy in
[`tests/examples/test_noncartesian_recon.py`](../tests/examples/test_noncartesian_recon.py).

Two results from `02` that nothing else in this repository states. A non-Cartesian reconstruction
comes out indexed **`[x, y]`**, not the `[y, x]` every Cartesian `02` produces — a property of the
adapter's axis order, established by putting a point somewhere known. And **preconditioning
changes the path and not the answer**: preconditioned and plain conjugate gradient agree to 0.0003
at convergence, so density compensation is a speed choice there rather than a better estimator.

`01` is in [`tools/run_notebook_smoke.py`](../tools/run_notebook_smoke.py); `02` is lab-tier.

## `fat_sat/`

| | |
|---|---|
| [`01_build.ipynb`](fat_sat/01_build.ipynb) | `sc.modules.SaturationPrep`: the chemical shift is **signed and relative to water**, converted once, and the notebook traces it from `shift_ppm` through `offset_hz` to the `freq_offset` on the compiled file. The band-to-water margin printed across 1.5 / 2.89 / 3 / 7 T — 8 Hz of clearance at 1.5 T against 816 Hz at 7 T — and the refusal when a band would reach water. Then the preparation in composition with `GRE2DTR`, and the same GRE without it, at one protocol. Two `.seq` files. **Needs nothing but `seqcraft`.** |

**The check the notebook exists to make is one line of output:**

```text
fat_sat  block 1 rf: use='saturation'   freq_offset=   -424.5 Hz
plain    block 1 rf: use='excitation'   freq_offset=     +0.0 Hz
```

A sign error in the chemical shift produces legal Pulseq, legal timing, legal gradients and a
perfectly ordinary waveform, and saturates **water**. Nothing downstream notices — the sequence
compiles, every k-space check passes, and the image comes back with the wrong tissue suppressed.
`write_epi_se_rs.py` and `writeEpiSpinEchoRS.m` reach the same number by different routes, one carrying the sign in
the ppm constant and the other applying it at the point of use, so an implementation that mixed
the conventions would be exactly this wrong. That is why the module owns the conversion and why
the notebook reads the answer off the *compiled sequence* rather than off the constructor.

CEST saturation trains, spatial saturation slabs and water excitation are out of scope: each is a
different physical contract rather than a mode of this one. `01` is in
[`tools/run_notebook_smoke.py`](../tools/run_notebook_smoke.py).

## `gre_spiral_2d/`

| | |
|---|---|
| [`01_build.ipynb`](gre_spiral_2d/01_build.ipynb) | `sc.modules.SpiralReadout` in a complete acquisition: `Excitation` + one arm + an angle schedule + spoiling and TR fill, assembled in the notebook. Eight interleaves from one readout instance, the trajectory measured off the compiled file, and the alignment check the notebook exists for. One `.seq`. **Needs nothing but `seqcraft`.** |
| [`03_flow_comp.ipynb`](gre_spiral_2d/03_flow_comp.ipynb) | Flow compensation on a spiral-in gradient echo, where the first moment the prephaser must cancel **rotates with the interleaf** on both in-plane axes at once. One designed pre-echo interval replaces the prephaser and the slice rephaser, compensating `x`, `y` and `z` with one echo time across all eight interleaves. Introduces `sc.PhysicalDesignScope`. One `.seq`. **Needs nothing but `seqcraft`.** |

**The join is the point.** `SpiralReadout` says where the trajectory crosses the origin; it does
not say where the *echo* is, and it must not — a readout that claimed to own TE would hide the
failure where a legal sequence has its echo time wrong by about one arm. So this layer computes

```text
TE = readout start  +  crossing offset  -  RF effective centre
```

from two numbers each module reports about itself, and then checks the result against the
**compiled** sequence: `t_excitation` and `t_adc` both come from `sc.kspace`, which shares no code
with either module. TE comes out identical across all eight interleaves **to 0.000 ns**, and
matches the reported value exactly.

A spiral GRE is one arm, an excitation and a list of angles, so `01` writes those three lines out
rather than wrapping them.

[`02_simulate_and_reconstruct.ipynb`](gre_spiral_2d/02_simulate_and_reconstruct.ipynb) is **the
off-resonance notebook**, and it reuses `noncartesian_recon.py` with no additions. `01` writes a
single-shot file as well as the interleaved one, because the claim needs a long readout: at 100 Hz
the 2.73 ms interleaved acquisition holds one voxel in one pixel and the 21.7 ms single-shot one
spreads it over 108. Eight time segments put it back — and **four segments is worse than none**,
because a four-point phase screen aliases 2.17 cycles of phase, which is the wrong lesson a reader
would otherwise draw from a real measurement.

`01` is in [`tools/run_notebook_smoke.py`](../tools/run_notebook_smoke.py); `02` is lab-tier.

## `se_spiral_2d/`

| | |
|---|---|
| [`01_build.ipynb`](se_spiral_2d/01_build.ipynb) | **Semantic timing ownership.** A 90, a `Refocusing`, and an `'in-out'` arm placed so that its origin crossing lands on the spin echo — measured against the compiled sequence, with pypulseq's own refocusing conjugation applied, at 2.0 µs residual and 0.02 Δk. Then **three wrong placements that all compile**. Two `.seq` files. **Needs nothing but `seqcraft`.** |
| [`02_simulate_and_reconstruct.ipynb`](se_spiral_2d/02_simulate_and_reconstruct.ipynb) | What the 180 is worth: a $T_2'$ sweep in which the spin echo does not move and the gradient echo follows $e^{-\mathrm{TE}/T_2'}$ to 1 %, off-resonance refocused to 0.0000 cycles, and the brain. **Needs `seqcraft[recon]`, `MRzeroCore` and a phantom download.** |

**A crossing through `k = 0` is not a spin echo**, and the two files exist to say so with one
variable. Their in-plane trajectories are bit-identical and both sample the origin at 14.002 ms;
one has a refocusing pulse and the other does not, and at that instant one carries zero
off-resonance phase and the other carries 0.7 of a cycle at 50 Hz.

```text
SpiralReadout   owns   where the trajectory crosses k = 0, and reports it
this layer      owns   where the readout starts, so that a crossing lands on the spin echo
```

Nothing enforces the join, which is why `01` builds the wrong versions too. Starting the readout
at the echo refocuses the magnetisation at $k_{\max}$, 32 Δk out. "The echo is in the middle of
the readout" is correct for `'in-out'` and wrong for the other three variants — worst for
`'out-in'`, whose midpoint falls *between* its two crossings.

The composition is the placement arithmetic and an angle list, both written out in `01`.

`SpiralReadout.time_to_echo` is measured from the start of the block `build()` returns, not from
the start of the arm — the prephaser sits between them, and for this readout that is 300 µs.
Measuring from the wrong one is the first of `01`'s three placements.

`01` is in [`tools/run_notebook_smoke.py`](../tools/run_notebook_smoke.py); `02` is lab-tier.

## `dwi_se_epi_2d/`

| | |
|---|---|
| [`01_build.ipynb`](dwi_se_epi_2d/01_build.ipynb) | What a b-value is, the same-polarity gradient pair a spin echo needs, and the trapezoid expression that turns a target `b` into a waveform — including the finite-ramp terms. Then the weighting the sequence actually delivers, the background weighting the imaging gradients contribute, an EPI readout with its centre-line echo on the spin echo, and why every b-value is acquired at one echo time. Four `.seq` files and a nominal `.npz`. **Needs nothing but `seqcraft`.** |
| [`02_simulate_and_reconstruct.ipynb`](dwi_se_epi_2d/02_simulate_and_reconstruct.ipynb) | The mono-exponential model, a known diffusion coefficient recovered to 0.06 %, what a mismatched echo time costs, and a brain ADC map against the phantom's own diffusion map — with an honest account of why one shot per b-value is not yet quantitative. **Needs `MRzeroCore`, `torch` and a phantom download.** |

`DiffusionSEPrep` owns the excitation, the refocusing pulse and both diffusion lobes together,
because the lobe separation Δ is measured **across** the 180°: neither lobe can be designed
without knowing where that pulse sits and how long it lasts.

```text
DiffusionSEPrep   owns   the lobes, the echo time, and where the 180 goes
the acquisition   owns   the echo time across the b-value list
```

The second line is the part no module can decide. A diffusion coefficient comes from a ratio, so
every b-value has to be acquired at the same echo time or the ratio carries `T2` as well as `D` —
`01` §6 derives the bias and `02` §3 measures it.

`sc.b_value` measures the weighting a whole tree delivers, by integrating the emitted gradients.
It is what shows that a nominally unweighted acquisition is not unweighted, which no check on the
diffusion lobes alone could see.

`01` is in [`tools/run_notebook_smoke.py`](../tools/run_notebook_smoke.py); `02` is lab-tier.

## `t2prep_gre_2d/`

| | |
|---|---|
| [`01_build.ipynb`](t2prep_gre_2d/01_build.ipynb) | What a T2 preparation is made of — fourteen hard pulses at one B1, an MLEV-4 composite refocusing train, a single `-90` tip-up — the two instants `prep_time_s` is measured between, and the simple and composite tip-up choices with their `prep_time_s` endpoints and their RF cost. Then the emitted block checked for interval symmetry and for the absence of any gradient while the magnetisation is transverse, and a composition with a shipped `GRE2D`. `.seq` files with both tip-ups. **Needs nothing but `seqcraft`.** |
| [`02_simulate_and_reconstruct.ipynb`](t2prep_gre_2d/02_simulate_and_reconstruct.ipynb) | A known T2 recovered from the slope of `ln S` against `prep_time_s`, within 1 % from 40 to 250 ms — and the control that matters: T2\* swept by a factor of five moves the answer by 3 %. Then where the residual bias comes from, and a B1 / off-resonance sweep across both tip-ups that finds no consistent advantage for either. **Needs `MRzeroCore` and `torch`.** |

A conventional gradient echo does not refocus static dephasing, so it follows T2\* rather than
refocused T2 decay — and at the short echo times a fast scan needs there is hardly any decay at
all. A **preparation** adds controlled T2 weighting before the imaging train, as a factor on the
longitudinal magnetisation:

```text
90x  ->  refocusing train  ->  -90x  ->  spoiler  ->  any readout at all
     |<------ prep_time_s ---->|
```

```text
T2Prep           owns  the pulses, the train, the timing solve and the spoiler
the acquisition  owns  the line ordering, and it matters
```

The second line is the part worth reading twice. What the preparation stores is *longitudinal*
magnetisation, which recovers at T1 — so every millisecond between the tip-up and the acquisition
of the centre of k-space erodes the contrast. A linear ordering reaches k = 0 after 315 ms of a
618 ms train; a **centric** one reaches it after 6. `01` measures both.

`01` is in [`tools/run_notebook_smoke.py`](../tools/run_notebook_smoke.py); `02` is lab-tier.

## `pc_gre_2d/`

| | |
|---|---|
| [`01_build.ipynb`](pc_gre_2d/01_build.ipynb) | What a velocity encoding is made of — two lobes of equal area and opposite polarity, solved in closed form from the gradient limits — the factor of two between `delta_m1` and one toggle's first moment, and what `polarity` means. Then the emitted waveform checked for net area zero, for the first moment `venc_m_s` implies, and for independence of where the pair is placed; then a phase-contrast repetition composed by hand from leaves. Two `.seq` files, one per polarity. **Needs nothing but `seqcraft`.** |
| [`02_simulate_and_reconstruct.ipynb`](pc_gre_2d/02_simulate_and_reconstruct.ipynb) | A moving spin simulated through both toggles: the phase difference is `pi * v / venc` to five decimal places. Then a static control that shows the leftover sign is the simulator's phase convention rather than the encoding's, what happens above VENC, and a velocity map of three voxels moving at three speeds. **Needs `MRzeroCore`, `torch` and `seqcraft`.** |

A gradient echo measures how much signal a voxel has, not how fast it is moving. A bipolar pair
adds the second:

```text
-A          +A            net area zero      -> stationary spins untouched
____        ####          first moment != 0  -> moving spins phase-shifted
    \      /
     \____/
```

One acquisition is not enough, because a voxel's phase also holds coil phase, off-resonance and
susceptibility. So the pair is played twice with opposite sign and the two are subtracted:

```text
VelocityEncode   owns  the waveform, the closed-form design, and both toggles
the acquisition  owns  which polarity is played when, and the phase subtraction
```

`venc_m_s` is the velocity that gives exactly `pi` of phase *difference*. The quantity behind it
is `delta_m1 = 1 / (2 * venc)`, the **change** in first moment across the toggles — each toggle
carries half, and a caller who confuses the two is out by a factor of two. `polarity` is defined
as the sign of the emitted first moment, which is something `sc.moments` can read off the events;
which sign of flow a reconstruction calls forward is its own choice.

`01` is in [`tools/run_notebook_smoke.py`](../tools/run_notebook_smoke.py); `02` is lab-tier.

## Requirements

Building needs only `seqcraft`. Simulating and reconstructing need `MRzeroCore`, `torch` and
`sigpy` — `pip install "seqcraft[sim,recon]"`.

Build notebooks run in seconds. Simulation notebooks vary from about 15 s to a few minutes; the
expensive parts sit behind a named switch that is off by default, and each one prints its budget
in the first cell. Writing a new one: see
[*What a simulation notebook is allowed to cost*](../docs/writing_examples.md).

## `phantom.py`

**The BrainWeb slab every 2D simulation runs against, prepared in one place.** Load, resample in
plane, cut the slab, move it to isocentre, hang an eight-element receive ring on it — twenty lines
that five notebooks had each grown their own copy of, and that had drifted: two different download
URLs and two spellings of the same slice range. The *parameters* stay in the notebooks, because
those are real choices with reasons attached — [`gre_2d/02`](gre_2d/02_simulate_and_reconstruct.ipynb)
wants four slices, [`fse_2d/02`](fse_2d/02_simulate_and_reconstruct.ipynb) wants two because it runs
three acquisitions — but the preparation does not.

Two things in it are worth knowing about:

- **`field_hz`** returns the **phantom's own B0 map** — the one `slab()` leaves in place and every
  example already simulates against. It is worth knowing what that map is, because the file name
  invites the wrong assumption: `subject05.npz` carries PD, T1, T2, T2′ and D and **no `B0_map`**,
  so MRzero's loader falls back to its own `generate_B0_B1`, whose comment says *"generate a
  somewhat plausible B0 and B1 map, visually fitted"*. Two Lorentzians, demeaned over the proton
  density.

  So it arrives centred on zero — nothing here shims it — and it spans about −11 to +45 Hz, which
  at 128 lines and a 700 µs echo spacing is roughly −1 to +4 pixels of displacement. Real, mild,
  and the same map the other five notebooks run against, which is what makes their numbers and the
  EPI ones comparable.
- **`same_frame`** is the one transpose, and the orientation paragraph above it is the part worth
  reading. The phantom is indexed `[x, y]` with **x left-right and y anterior-posterior, anterior
  at high y** — read off the anatomy rather than assumed, because a first draft had those two
  swapped, and a picture drawn on that assumption is silently sideways. A reconstruction is `image[ky, kx]`, so it
  needs no transform at all and a *map* needs one; and every example draws with `origin='lower'`,
  without which the frontal lobe comes out at the bottom.

## `data/`

Phantoms, shared by every example and downloaded on first use. `gre_2d/02` fetches BrainWeb
subject 05 (~19 MB) from the MRzero-Core repository into `examples/data/`; a second example wanting
a phantom should look there before fetching its own copy.

Not tracked by git — it is somebody else's data, it is the same for everyone, and one download is
cheaper than carrying it in the history. Outputs under any `seq/` are build products and are not
tracked either.

## What used to be here

Two DTI acquisitions — a spiral and an EPI — sat in `_parked/` as the acceptance test for whatever
module set was written next, along with the `lib/` helpers they needed. That job is discharged:
[`docs/adr/003`](../docs/adr/003-scanner-and-module-reform.md) asked for a module set chosen by
making real scans work, and `gre_2d/` is the one that was written against.

They are in git history, and the physics they depended on — the b-value solve, the
variable-density spiral, the EPI ramp-sampling moment integral — is in [`salvage/`](../salvage/) as
plain functions. Two of the three have since been written as modules, in
[`gre_spiral_2d/`](gre_spiral_2d/) and [`dwi_se_epi_2d/`](dwi_se_epi_2d/).

`mr0_bridge.py` went with them: every `02` uses `mr0.Sequence.import_file` on the written `.seq`,
which is a stronger check than converting the tree, because it tests the file a scanner would
actually play.

## Contributing

Writing or editing an example? See [`docs/writing_examples.md`](../docs/writing_examples.md).
