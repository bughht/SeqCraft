# Examples

| | What it covers |
|---|---|
| [`01_getting_started.ipynb`](01_getting_started.ipynb) | Blocks, `Opts` and `compile`; the overlap rules, provenance, checking against physics, writing a file — and `sc.Module` at the end, once there is a reason for one. **Uses no modules**, which is the point: it is the demonstration that the compile path stands alone. |
| [`gre_2d/`](gre_2d/) | A complete spoiled 2D gradient echo, built three ways and then simulated and reconstructed. This is where most of `sc.modules` came from. |
| [`mprage_2d/`](mprage_2d/) | An inversion-prepared segmented GRE: where `IRPrep` and `GRE2D.time_to_center_line` came from, and where an inversion time gets placed to the microsecond. Defines `MPRAGE2D` in its own notebook. |
| [`mp2rage_2d/`](mp2rage_2d/) | One inversion, two trains, and the `SET` label that separates them — plus the ratio that cancels the receive field. Defines `MP2RAGE2D` in its own notebook. |
| [`se_2d/`](se_2d/) | A spin echo, and the one rule that makes it one: between refocusing centres, area before the echo equals area after it. Where `Refocusing` and `CartesianLine(prephase=False)` came from. Defines `SE2D` in its own notebook. |
| [`fse_2d/`](fse_2d/) | The same composition at 1, 16 and 72 echoes — TSE and HASTE — and what an echo train costs in blurring, ghosting and signal. Defines `FSE2D` in its own notebook. |
| [`gre_epi_2d/`](gre_epi_2d/) | The whole of k-space in one shot: the centred sampling window, the blip on the zero crossing, ramp sampling and the operator that undoes it, off-resonance and the N/2 ghost, and GRAPPA. Where `EPI2D` came from. Defines `GREEPI2D` in its own notebook. |
| [`se_epi_2d/`](se_epi_2d/) | The same readout after a refocusing pulse — and the measurement that a spin echo **does not** fix EPI distortion. Defines `SEEPI2D` in its own notebook. |
| [`megre_2d/`](megre_2d/) | Eight echoes off one excitation, monopolar and bipolar, fitted for T2\* and ΔB0 against a phantom that carries the ground truth for both. **Defines no class**, which no other directory here can say. |
| [`gre_spiral_2d/`](gre_spiral_2d/) | The three things a spiral is — the path, the traversal, and what the amplifier actually does — with three wrong versions built too, because all three compile. Variable density as data, the ADC split the sample limit forces, and a reconstruction that puts the field map *inside* the forward model. Where `Spiral2D` came from. Defines `GRESpiral2D` in its own notebook. |
| [`se_spiral_2d/`](se_spiral_2d/) | Spiral in/out after a refocusing pulse, where **the seam, the `k = 0` sample and the spin echo are one instant** — and the same sequence built with the readout placed against its first sample instead, which compiles. Defines `SESpiral2D` in its own notebook. |
| [`matlab/`](matlab/) | The same small GRE 2D built two ways with official MATLAB Pulseq events: direct `LogicBlock` assembly and an example-only single-repetition `GRE2DTR` Module. Both compile through Python and return `mr.Sequence`. |

## `gre_2d/`

| | |
|---|---|
| [`01_build.ipynb`](gre_2d/01_build.ipynb) | The sequence, three times over: raw pypulseq events, the four leaf modules composed inline, and the same composition written as a module. Every arithmetic check, three sampling patterns, three `.seq` files. **Needs nothing but `seqcraft`.** |
| [`02_simulate_and_reconstruct.ipynb`](gre_2d/02_simulate_and_reconstruct.ipynb) | Those three files against a BrainWeb phantom — PD, T1, T2, T2′, D and a synthesised B0, all six simulated — with an eight-element receive ring, then one reconstruction across three samplings. **Needs `seqcraft[sim,recon]`** and a one-off ~19 MB phantom download. |

`01` deliberately does **not** import `GRE2DTR` or `GRE2D`. It writes them, because the third pass
is where a reader sees a working composition become something reusable, and a pass that imported
the answer would teach nothing.
[`tests/modules/test_notebook_matches_the_package.py`](../tests/modules/test_notebook_matches_the_package.py)
asserts that what the notebook writes and what the package ships produce identical events, so the
tutorial cannot drift from the library without CI noticing.

Every module in `sc.modules` was extracted from this pair. That is the rule the library is built
on: a module that cannot be extracted without altering the sequence is not a module, and one whose
extraction does not shorten the notebook is a wrapper.

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

## `se_2d/` and `fse_2d/`

| | |
|---|---|
| [`se_2d/01_build.ipynb`](se_2d/01_build.ipynb) | One interval placed by hand out of the shipped modules, then the whole shot, then `SE2D`. Every check is `sc.kspace` at the echo *sample*, signed, because `|k|` is symmetric and a k-space extent check passes on a mirrored image. **Needs nothing but `seqcraft`.** |
| [`se_2d/02_simulate_and_reconstruct.ipynb`](se_2d/02_simulate_and_reconstruct.ipynb) | The one claim arithmetic cannot make: the spin echo recovers the phantom's **T2** and the same sweep without the 180 recovers **T2\***. Then the image, and which weighting the protocol landed on. **Needs `seqcraft[sim,recon]`.** |
| [`fse_2d/01_build.ipynb`](fse_2d/01_build.ipynb) | `FSE2D`, turbo 1 → 16 → 72 on one instance, three orderings as data, the echo-band warning, and HASTE with `partial_fourier`. Four `.seq` files. **Needs nothing but `seqcraft`.** |
| [`fse_2d/02_simulate_and_reconstruct.ipynb`](fse_2d/02_simulate_and_reconstruct.ipynb) | The echo envelope measured, the point-spread width per ordering, contrast and blurring as two separate knobs, and **the ghost a scattered table makes** — a periodic modulation of `ky`, which is a replica of the object rather than a blur. Then HASTE with POCS. **Everything is measured on one spin**: a ghost is a modulation of `ky`, and a point object's k-space *is* that modulation — 1.3 s per sequence instead of 173 s, and no reconstruction in between. **Needs `seqcraft[sim,recon]`**; runs in 16 s. |

`SE2D` and `FSE2D` are **defined in those notebooks and do not ship**, for the same reason
`MPRAGE2D` and `MP2RAGE2D` do not: one consumer each.
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

## `gre_spiral_2d/` and `se_spiral_2d/`

| | |
|---|---|
| [`gre_spiral_2d/01_build.ipynb`](gre_spiral_2d/01_build.ipynb) | **The path and the traversal are two questions**, and `v = 0` at both ends of the arm is the design decision the module rests on — it is what makes one arm reversible and therefore four variants one class. The slew-limited regime measured rather than assumed, with the `max_grad` crossover *found*; seven densities side by side, with the turn spacing measured off the built trajectory; the four variants and `echoes`; the ADC split at one shot and the 6.6 `Δk` hole each seam costs. Then **three wrong versions**, because all three compile: a rephaser overlapping the arm, a join measured in isolation (0.024 `Δk` of residual `k`), and a naive per-axis prephaser at 45° (250 T/m/s on a 180 T/m/s amplifier). `GRESpiral2D`, six imaging `.seq` files at **one flip, one TE, one TR and one matrix**, and the dual-echo Cartesian reference that calibrates them. **Needs nothing but `seqcraft`.** |
| [`gre_spiral_2d/02_simulate_and_reconstruct.ipynb`](gre_spiral_2d/02_simulate_and_reconstruct.ipynb) | The calibration first — B0 map, ESPIRiT maps, and the wrap check as an assertion. Then the blur **on one spin**, against a `B0 = 0` control and against readout length, because a point object's point spread *is* the answer. The brain uncorrected, corrected, and **corrected with the sign flipped** — the failure that reads as "correction does not help". The exact operator against the segmented one at 64², four interleaves solved jointly against four solved separately, the variable-density file at `edge_undersampling = 2.0`, T2\* from the three-echo file, and the gradient-delay sweep, which is the one calibration neither sidecar nor field map supplies. **Needs `seqcraft[sim,recon]`.** |
| [`se_spiral_2d/01_build.ipynb`](se_spiral_2d/01_build.ipynb) | `Spiral2D(variant='in-out', prephase=False)` after a refocusing pulse, with the 2D dephaser **before** the 180 carrying the sign the conjugation will flip — one oblique trapezoid, because two independent ones reach `sqrt(2)·max_slew` and compile cleanly. TE as **two numbers**, the gradient echo and the spin echo, brought within half a dwell of each other. **Built wrong first**: the readout placed against its first sample rather than its seam, which compiles, whose every `k` is right, and which acquires the centre of k-space a whole arm from the spin echo. Four `.seq` files at one echo time. **Needs nothing but `seqcraft`.** |
| [`se_spiral_2d/02_simulate_and_reconstruct.ipynb`](se_spiral_2d/02_simulate_and_reconstruct.ipynb) | That it is a **T2** echo and not a T2\* one, swept over `T2'`. Then the section the directory exists for: **an in-out spiral's point spread is nearly real and an out spiral's is not**, measured on a point object at the same echo time and the same arm — and what that does to a brain. Correction on both and the residual after it, then the honest negative: correction fixes the blur and does nothing for the through-plane dephasing a 3 mm slice in the same field is also suffering. **Needs `seqcraft[sim,recon]`.** |

`Spiral2D` is the one module that ships out of this pair; `GRESpiral2D` and `SESpiral2D` are
**defined in those notebooks and do not ship**, one consumer each. Both `01`s are in
[`tools/run_notebook_smoke.py`](../tools/run_notebook_smoke.py).

[`noncartesian_recon.py`](noncartesian_recon.py) sits beside [`phantom.py`](phantom.py) for the
same reason that one does: it is shared example code, and it is **not** in `src/`, because the
package builds sequences and a reconstruction module there would put an ESPIRiT dependency on
the compile path. Its operator is composed from `sigpy.linop` primitives rather than
hand-written, so the adjoint is *derived* — a conjugate gradient on a nearly adjoint pair does
not raise, it converges to the wrong image and looks under-regularised.

## `megre_2d/`

| | |
|---|---|
| [`01_build.ipynb`](megre_2d/01_build.ipynb) | The four numbers that are wrong *quietly*: the reverse lobe's grid (**0.2051 Δk** on the lobe that is right for a single line, built the wrong way first and plotted so a reader sees that the two trajectories *look* identical); the fly-back's area, with all three plausible wrong answers built and their trajectories drawn; the dwell quantum and why it is the RF raster; and `te_s` against `echo_spacing_s`, whose consequence — **0.75 % low everywhere** in a two-point field map — is computed from `te_s` alone, with no simulator. Plus the barrier table, which reports a prediction that turned out **false**. Three `.seq` files, **one flip, one TR, one matrix and one slice across all of them**, checked and printed. **Needs nothing but `seqcraft`.** |
| [`02_simulate_and_reconstruct.ipynb`](megre_2d/02_simulate_and_reconstruct.ipynb) | **One reconstruction for all three files** — ESPIRiT coil maps estimated once (a coil sensitivity belongs to the coil, and using one set is what makes echo *e*'s phase comparable with echo 0's), POCS for the partial echo, then SENSE, because RSS throws away the phase that is half of this notebook. All eight echoes shown in magnitude *and* phase before anything is fitted, then T2\* and ΔB0 fitted for **all three readouts** against the phantom's own `T2dash_map` and `field_hz()` — maps, residuals and scatter side by side, not just the best one. **The phantom is 192 mm across and the protocol images 220 mm**, so the oracle is resampled onto the reconstruction's pixel grid before anything is subtracted, by nearest neighbour so it stays a set of values the phantom contains; the alignment is checked with a Dice overlap rather than assumed. Finally what partial echo costs and what POCS gets back, measured against the full-echo image. **Needs `seqcraft[sim,recon]`.** |

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

## Requirements

Building needs only `seqcraft`. Simulating and reconstructing need `MRzeroCore`, `torch` and
`sigpy` — `pip install "seqcraft[sim,recon]"`.

### What a simulation notebook is allowed to cost

MRzero's inner loop is dense complex linear algebra and it takes every core it is offered — sixteen
threads by default on a 32-core machine. A draft of [`fse_2d/02`](fse_2d/02_simulate_and_reconstruct.ipynb)
ran several 16-echo, 9-shot acquisitions in a single cell, minutes of saturated CPU with no output
until it finished, and took a workstation down.

So the spin-echo and EPI `02` notebooks each set `SIM_THREADS` **before importing torch**, print
their budget in the first cell, and keep anything expensive behind a named switch that is off by
default. `se_2d/02` is 38 s. `fse_2d/02` is 15 s of measurement plus 4 minutes of images, and the
images are the part behind the switch — every *number* in it is measured on a single spin, because a
ghost is a modulation of `ky` and a point object's k-space is that modulation.

`megre_2d/02` is three acquisitions and 24 SENSE reconstructions, and it has no switch because it
does not need one: it **picks up a GPU when there is one**, and the MRzero simulation — which is
essentially all of its cost — is a torch computation that moves there wholesale. Both `DEVICE` and
sigpy's `RECON_DEVICE` are chosen by *probing* rather than by asking, because
`torch.cuda.is_available()` reports the driver and comes back true on a machine whose devices are
all hidden, and because cupy needs a toolkit to compile against where torch needs only a driver to
talk to — so the two can disagree, and a machine can simulate on its GPU with no working cupy at
all. Either way the notebook is the same notebook, and every step prints what it cost.

The two EPI `02`s take the same bargain at a different scale: about a minute of one-voxel
measurement each, then `BRAIN_IMAGES` for the slab. `gre_epi_2d/02`'s switch also covers §7, which
is the one section that **cannot** be measured on a single spin — parallel imaging is a statement
about receive sensitivity, so it needs coils and an object or it needs nothing.

If you add a simulation to an example, price it first: the cost is voxels × repetitions × states, and
a 16-echo train has all three. Two of the three are not negotiable — the repetitions are the
sequence, and dropping `max_state_count` from 200 to 32 is 5 % wrong because a CPMG train's
stimulated pathways *are* the physics. The slab thickness was the one that turned out to be free, and
only because it was measured rather than assumed.

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
plain functions, which is where a future spiral or EPI module should start from. `mr0_bridge.py`
went with them: `02` uses `mr0.Sequence.import_file` on the written `.seq`, which is a stronger
check than converting the tree, because it tests the file a scanner would actually play.
