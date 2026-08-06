# Changelog

## Unreleased — compiler revision, phases A and B

Plan and reproductions: [PLAN_COMPILER_V2.md](PLAN_COMPILER_V2.md). Every item below has a test that
fails without it.

### Fixed — silent wrong output

- **Unknown event types were dropped, not rejected.** `rot3D`, `soft_delay` and `rf_shim` matched none
  of the compiler's positive branches, so they were placed with a zero-width reservation, collected by
  nothing, and vanished. A tree carrying a rotation extension compiled clean, reported no issues, and
  played **unrotated**. Now a whitelist: every pypulseq type is either emitted or rejected by name,
  with the way round. An introspective test reads pypulseq's own source, so a new event type in a
  future release fails loudly here instead of disappearing.
- **A label between two ADCs could address the wrong one.** A label is a running register, so it may
  live in any block after the previous ADC and at or before its target — but the compiler attached it
  by containment, and a boundary pushed later by a gradient put it in the *previous* readout's block,
  overwriting that k-space address. Measured: `LIN = [7, 7]` where `[1, 7]` was meant. Labels now
  attach to the block holding the first ADC at or after their own time, which is independent of where
  boundaries land. A barrier is never crossed, so it remains a way to order labels explicitly.

### Fixed — legitimate sequences rejected

- **A trigger crossing a block boundary failed to compile.** Triggers held their block open while
  contributing nothing to boundary selection, so a physio trigger overlapping a readout raised
  `block 0 spans 660.0 us but its events need 2500.0 us` and blamed the module. Reservations now
  separate two properties that were conflated: *indivisible* (no cut inside — RF, ADC, trigger) and
  *exclusive* (at most one per block — RF and ADC only). Triggers are the first; several may share a
  block, as pulseq allows.

### Added

- **An exact waveform oracle** (`tests/compiler/fidelity.py`). Every pulseq gradient is piecewise
  linear, so two are equal everywhere iff they agree at the union of their knots — a proof rather than
  a sampled approximation. Blocks are *concatenated*, not summed, since summing double-counts every
  seam by the amplitude there. Includes a seam-continuity check and a self-test that a 1 % corruption
  is detected, because an oracle that cannot fail proves nothing.
- **Two new compile-time invariants.** Per-axis **m1** referenced to the sequence start, which catches
  a gradient playing at the wrong time — m0 is exactly what a time shift preserves. Computed in closed
  form on both sides: `g(t)·t` is quadratic between knots, so a trapz-based m1 disagreed with itself by
  0.1 % on nothing more than a split. And **label addresses** against the fold of the tree's labels,
  which catches an addressing shifted by one readout — the duplicate check only sees collisions.
- **Errors where there was silence.** A barrier inside an indivisible span; a mandatory gap nothing can
  be cut in; order-dependent labels aimed at one readout (pypulseq sorts a block's extensions by
  library id, so `add_block(set, inc)` and `add_block(inc, set)` are the same block — verified — and
  only order-independent groups are expressible). A label with no following ADC is reported as a
  warning.

### Performance

- **Boundary selection was quadratic on EPI**, the case its own fallback exists for: a continuous
  readout gradient makes both natural candidates unacceptable, so every echo needed a gap midpoint and
  the mark set was re-sorted each time. One monotone scan instead — byte-identical output, 32× faster
  at 6400 echoes (1250 ms → 39 ms), ratio per doubling 3.7 → 2.

### Fixed — precision (phase C)

- **Splitting or merging a gradient no longer resamples it.** Both operations are exact on a
  piecewise-linear representation, and every pulseq gradient *is* piecewise linear — so the sum of
  several is evaluated at the union of their knots, where it can bend, instead of on a uniform
  raster grid. Splitting a spiral by a barrier used to round 2.5 % off its peak; it is now exact to
  1e-14 %, and the peak is preserved digit for digit. Sampling was exact only when every knot
  happened to land on the grid, and an arbitrary gradient's samples sit at raster *centres*, so
  they never did.
- **A split arbitrary gradient stays an arbitrary gradient.** A block boundary is on the raster, so
  cutting there leaves each piece's samples at the centres of its own raster intervals, with the
  seam amplitude becoming one piece's `last` and the other's `first`. Recognising that pattern is
  what keeps a spiral a spiral.
- **The one case pulseq cannot represent is now reported, not silent.** A trapezoid's corners are on
  raster edges and an arbitrary gradient's samples at raster centres; their sum bends at both and no
  pulseq gradient event has room for it. That resample now raises a `grad_resample` warning carrying
  a *measured* bound on how far the waveform moved — exact, because two piecewise-linear functions
  differ most at a knot of one of them.
- **A gradient started off the gradient raster is an error.** It used to be snapped silently by up
  to half a raster: a gradient asked for at 5 µs played at 10 µs. There is no correct snap to make,
  and which way to round is the caller's decision.
- **m0 was blind by construction.** Both sides of the invariant were integrated by the same
  approximation, so their errors cancelled — the 2.5 % peak loss above left m0 agreeing to 1e-14.
  Moments are now computed from exact knots on both sides, by Gauss-Legendre over each segment,
  which is exact for every order the API offers.

### Changed — one implementation each

- `events.knots_of` and `events.pwl_moment` are the single exact-gradient primitives; the compiler's
  private copies are gone. `moment_of` is now exact and its `raster` argument is ignored (kept so
  existing calls work). `waveform_of` remains for plotting, which genuinely wants uniform samples.
- `CompiledSequence.moments()` integrates from knots rather than raster samples, so it is now correct
  for arbitrary waveforms at every order.
- Definition merging had two implementations, one of them dead and better. The compiler now uses
  `validate.merge_definitions`, so a conflict names *which two sources* claimed the key rather than
  reporting "already" and "also".
- `_sample` and `_reduce_corners` deleted. The second existed only to undo the damage the first did:
  a uniform resample turned a merged trapezoid into hundreds of collinear points that then had to be
  reduced back to corners. The exact union of knots is minimal by construction.

## 0.3.0 — the logic-block rewrite

A complete redesign. The previous architecture is gone, not deprecated.

### The new model

Three concepts replace the previous five layers.

- **`LogicBlock`** — a tree of pulseq events and nested blocks, each with a start time. Two
  attributes (`tag`, `nodes`) and one method (`add`). `duration` is measured from the nodes, never
  declared. Overlap is legal everywhere; making it legal for pulseq is the compiler's job.
- **`Module`** — `__init__` designs, `build()` returns a logic block, and timing a caller needs in
  order to place it is a plain property. No class variables, no categories to register, no hooks.
  `core/module.py` went from about 1000 lines to under 200.
- **`sc.compile`** — finds block boundaries, sums same-axis gradients, and validates the compiled
  waveform. About 700 lines, and the only complicated thing in the package.

### Removed

`BuildResult`, `BlockSpec`, `SegmentRef`, `Concurrent`, `Sequential`, `Sequence`, `Prewinder`,
`Rewinder`, `MomentSpec`, `CATEGORY`, `ANCHORS`, `EMIT_PARAMS`, `SEGMENTS`, `HAS_ADC`, `REQUIRED`,
`_build`, `_emit`, `_validate`, `with_()`, the frozen-parameter machinery, the design cache and its
fingerprint, and all of `timing.py`'s gap solving. TE is now arithmetic you can read:
`t0 + te - ro.time_to_echo`.

There is no `Sequence` class. A sequence *is* a logic block; compiling it produces a
`CompiledSequence` holding the pypulseq object, the report and per-block provenance.

### Compiler behaviour

- Gradients on **different axes** overlapping: silent. It is the normal way to build a sequence, and
  warning about it would only teach people to ignore warnings.
- Gradients on the **same axis**: warned, then summed, with both sources named.
- Two **RF** or two **ADC** events overlapping, or an RF and an ADC: an error naming both tag paths
  and the overlap in microseconds. Dead times and ringdown are inside the compared spans, so an RF
  too close to a readout is caught by name rather than by pypulseq 40 000 blocks later.
- **Amplitude and slew are measured on the compiled waveform**, because two individually legal
  gradients on one axis can sum to an illegal one — an area-100 plus an area-200 trapezoid on a
  40 mT/m, 150 T/m/s system reaches 189 % of the slew limit. Vector-norm violations are warnings,
  since real amplifiers permit root-2 times the per-axis slew on two axes.
- An event's own `delay` is **preserved and reserved**, not folded away: the reservation begins at
  the node time, which already contains the RF dead time or the ADC dead time. In-block delays are
  quantised onto each event's own raster and computed in integer picoseconds.

### New modules

`SlabExcitation`, `HardExcitation`, `HardRefocusing`, `GaussSaturation`, `AdiabaticInversion`,
`SpiralVDS`, `NoiseAcquisition`, `Prephaser`, `MonopolarDiffusion`, `BipolarDiffusion`,
`ArbitraryDiffusion`, `FatSat`, `InversionRecovery`, `Barrier`, `SLRExcitation`,
`SLRRefocusing`, plus `dti_directions` and `direction_condition_number`.

There are **no recipes**. A recipe is somebody else's sequence choices in library code: to change
your scan you would edit a package. The notebooks assemble sequences from modules instead.

### Physics fixed along the way

Each of these was found by measurement, and each has a regression test.

- **Diffusion b-value was 2-4 % high.** The published ramp correction
  `+- eps^3/30 -+ delta*eps^2/6` appears with several incompatible conventions for whether `delta`
  includes the ramps. Replaced with an exact piecewise integral, now verified against numerical
  integration of the built waveform to 0.5 %. A 3 % b-value error biases every diffusivity by 3 %,
  and nothing in the reconstruction would reveal it.
- **The spiral generator ignored the `theta''` term**, so the realised slew ran several times over
  the limit right at the origin — where a diffusion measurement's best samples are. Rate-limiting in
  both directions fixed the violation and introduced a subtler problem in its place, described next.
- **The spiral was not slew-limited at all — it ran at constant angular velocity.** Rate-limiting the
  *descent* by the leftover tangential budget looks symmetric and is not: at the centripetal ceiling
  the whole slew budget is already going into turning, so the leftover is zero, so `theta'` could
  never change again. It froze at whatever rate it reached near the origin and coasted. The realised
  slew ramped linearly from 5 % of the limit at the centre to 89 % only at the very edge, averaging
  **45 %**, and the readout ran a third longer than the amplifier could deliver — with nothing in the
  sequence looking wrong, because every limit check passed. Peak is now 100 % and the **mean is
  99.8 %**, and a 96-matrix single-shot readout went from 28.95 ms to 19.36 ms for free.
- **Normalising the finished waveform onto `k_max` silently coarsened the sampling by 16 %.** The
  gradient used to be differentiated analytically, so the cumulative sum of the *sampled* waveform
  fell short of the analytic radius — near the origin the spiral turns an eighth of a radian per
  raster and the polygon cuts the corners. Scaling the result to land on `k_max` scaled the **turn
  spacing** with it, so a spiral asked for Nyquist delivered 16 % coarser and the aliasing looked like
  undersampling nobody requested. The waveform is now the exact finite difference of the analytic
  trajectory, so `cumsum(g) * raster` reproduces `k` to the last bit and no scaling is needed. The
  step is chosen against the **discrete** constraints the hardware actually applies rather than a
  continuous-time model of them, capped at the sustainable rate — without the cap a greedy pass
  accelerates past what the next turn can hold and, since shedding speed also costs slew, stalls with
  no legal move at all.
- **The trajectory sidecar was written on the gradient raster** while the ADC has ten times as many
  samples, leaving the interpolation to whoever wrote the reconstruction. It is now written at ADC
  sample times, with `t_adc_s` alongside, and checked against pypulseq's own `calculate_kspacePP`
  to 0.02 % of `k_max`.
- **The refocusing pulse had no CPMG phase.** A 180 needs a quarter turn of carrier phase from the
  excitation to be first-order insensitive to its own flip-angle error. Bloch-simulated: at 20 % low
  B1 the echo keeps 98 % of its amplitude with the `pi/2` phase and 83 % without. Now the default on
  every refocusing class, composed with any phase the caller sets rather than overwriting it.
- **`_flip_and_phase` in the MRzero bridge reported an SLR refocusing pulse as 0 degrees.** It
  integrated `2*pi * int B1 dt`, which is the small-tip approximation — and an SLR `'se'` pulse is
  designed so its *refocusing* profile is flat, leaving its envelope with almost zero net area. The
  whole sequence simulated as producing no signal, with nothing resembling a pulse problem to point
  at. Bloch-simulated over intervals of the pulse's own `rf.t` grid now.
- **The spiral's rewinder was designed for the unrotated end point.** Rotating an interleaf moves
  area between axes, so an axis needing `kx_end` unrotated can need `hypot(kx_end, ky_end)` at some
  other angle. Designed for the worst case and scaled down now.
- **The spiral's block was too short for its ADC's trailing dead time**, so the compiler had nowhere
  to put a boundary and summed the rewinder into the readout — a step discontinuity in the middle of
  the acquisition.
- **`BipolarDiffusion` did not null m1.** Inverting the second lobe made the effective waveform
  `(+ -) gap (+ -)`, which is not even about the encoding centre. Both lobes are built identically
  now; the 180 does the inverting, giving `(+ -) gap (- +)`. Measured, not asserted.
- **`dti_directions` left a 14-degree pair at n=30.** A golden-angle spiral is uniform in area, which
  is not the same as well separated. Relaxed under antipodal repulsion now; the design-matrix
  condition number is 1.58 for every n, against a floor of about 1.58.
- **The spin-echo recipe's winders had the wrong sign.** A 180 inverts accumulated phase, so a
  prephaser placed ahead of it must prephase the other way. With the conventional sign the readout
  ended at `3*k_max` and the phase encode was mirrored — the latter invisible in a k-space extent
  check, because `|k|` is symmetric.
- **The DTI recipe was missing the excitation's slice rewinder.** A spiral starts at k=0 in x and y,
  which makes it tempting to think nothing needs rephasing on z. But the excitation's slice-select
  gradient leaves through-slice dephasing equal to its own tail, and the refocusing pulse does not
  undo it: that pulse's gradient is symmetric about its centre, so its lead and tail cancel each
  other. `k_z` at the echo was -525 1/m — 2.1 cycles across a 4 mm slice, about 95 % of the signal —
  with nothing else looking wrong.
- **The DTI recipe cycled the spiral interleaf across volumes** instead of acquiring all interleaves
  of each volume, leaving every volume undersampled by the interleaf factor.
- **TR was applied per shot rather than per slice.** TR is the time between exciting the *same*
  slice, so every slice belongs inside one TR period. The bug did not produce a wrong image; it
  produced a correct one that took `n_slices` times longer — 5.5 hours instead of 16.5 minutes for
  30 directions and 20 slices.
- **The compiler's duration invariant used a fixed 1 ns tolerance**, which is below float64
  resolution on a sequence over an hour, so it flagged every long acquisition and printed two
  identical-looking numbers. Scales with the total now.
- **Gradient endpoints from resampling landed at about 1e-6 Hz/m instead of exactly zero.** pypulseq
  tests `first != 0` exactly and then demands the previous block continue it, so a rounding artifact
  became a continuity error hundreds of blocks away.
- **`FatSat` defaulted to a time-bandwidth product of 4**, giving a 500 Hz pulse against a 434 Hz
  fat-water separation at 3 T — it would have saturated water too, which looks like poor SNR rather
  than a bug. The default is 1.6 now, and a too-wide pulse is refused.

### Notes for anyone reading the sequences

- **A spiral's k-space is not in the `.seq` file** in a form a reconstruction can use, so
  notebook 2 writes a `.traj.npz` alongside it — at ADC sample times, not on the gradient raster.
- **`CartesianLine`'s prephaser ignored the readout gradient's ramp.** It was minus ``k_max``, which
  is the tempting value and short by ``amplitude * rise_time / 2`` -- the ramp carries area too. Every
  line was therefore displaced in k-space by that much: one ``dk`` on a wide readout and **3.95 dk**
  on the field map's short low-bandwidth one. A k-space offset is a linear phase ramp across the
  image, so a magnitude image looks entirely normal and anything reading the phase is wrong.

  The prephaser and ``time_to_echo`` now both come from one moment integral, so they cannot disagree,
  and ``time_to_echo`` points at k=0 rather than at the middle of the ADC window. The test that
  asserted ``pre.area == -k_max`` was asserting the bug; it now checks where k=0 actually lands, in
  pypulseq's own trajectory, and is parametrised over partial echo.
- **`CartesianLine` supports partial echo.** ``partial_echo`` is the fraction of the *leading* half of
  k-space acquired: 0.75 starts at ``-0.75 k_max``, still ends at ``+k_max``, and brings TE forward by
  the samples dropped. Resolution is unchanged, since ``dk`` per sample is fixed by the FOV -- what is
  given up is the conjugate-symmetric part, which the reconstruction has to fill in. The sample count
  snaps to the ADC divisor and the attribute reports the fraction achieved rather than the one asked
  for. It is also the case where ``time_to_echo`` stops being the window centre: k=0 arrives earlier,
  which is the entire point, and taking the centre would put TE late by half the dropped samples.
- **A spiral readout needs a field map, so notebook 2 writes one.** A dual-echo low-resolution GRE at
  the same FOV and slice, 3.1 s, with both echoes in one TR. Five things it has to get right, each of
  which produced a wrong map first:

  - **Echo spacing** is chosen as `1 / (fat-water offset)`, which puts fat back in phase with water so
    it cannot corrupt the phase difference, and sets the unambiguous range to `±½ΔTE⁻¹` so nothing
    needs unwrapping. The *readout duration* is picked to land the spacing there, not the reverse.
  - **A flyback, not a reversed gradient.** Reversing mirrors the sample positions and the readout's
    sampling is not exactly symmetric about k=0, so the echoes land on k grids a fraction of a sample
    apart. A k offset between echoes is a linear phase ramp, which over `2π·ΔTE` becomes an apparent
    field gradient of `1/ΔTE` across the FOV — ±210 Hz here, the whole unambiguous range, and it looks
    exactly like a badly shimmed magnet.
  - **The flyback area is the gradient's own area**, read off the event — not `2·k_max`. The ADC samples
    the flat top only, so the ramps add area the sampled span never shows: 272.7 against 266.7, and the
    6 1/m left over is 1.4 sample spacings.
  - **The receiver phase must equal the transmit phase.** RF spoiling changes it every shot;
    `CartesianLine.build(rf_phase_rad=...)` now sets the ADC's. Invisible in a magnitude image, fatal
    in a phase map — and invisible in simulation, because the bridge aligns the two on your behalf.
  - **Dummy excitations.** Longitudinal magnetisation settles to 43 % of equilibrium at this flip and
    TR, so without them the first phase-encode line carries 2.3 times the signal of the last — a filter
    along k-space, and the dominant artefact.

  With all five, the map recovers to 2.0 Hz RMS against the truth. `angle(echo0 · conj(echo1))` gives
  `+Δf`, because the signal accumulates `exp(-2πi·Δf·t)`; the other way round returns a map of the right
  shape and range and the wrong sign, which makes the correction add the blur it should remove.
- **A dark band through a simulated spiral image was the fat sat, not off-resonance.** MRzero's pulses
  are instantaneous ideal rotations, so they cannot be frequency selective: a −419 Hz fat-sat pulse
  simulates as an on-resonance 90° that saturates *water*. Worse, its spoiler lies on the slice axis,
  which a one-voxel-thick phantom forces the bridge to drop — so the leftover transverse magnetisation
  survives into the readout and interferes with the echo, cancelling it wherever the field takes the
  wrong value. The giveaway was that a *uniform* off-resonance offset, which cannot dephase anything
  spatially, quadrupled the total signal energy. Notebook 3 starts the simulation at the excitation;
  on the scanner the pulse is selective and its spoiler works, so the sequence is right and only the
  simulation had to change.

  The general lesson, which cost two bugs: **a z-only spoiler is invisible to a 2D simulation.** If a
  sequence relies on slice-axis spoiling, the simulation sees no spoiling at all.
- **Two gradient moments were being silently discarded by the MRzero bridge.** `to_mr0` resampled the
  gradient *at* event edges and applied the trapezoidal rule -- and a trapezoid starts and ends at
  zero, so any block shorter than about twice the event length integrated to nearly nothing, and one
  spanning a whole prephaser integrated to **exactly** zero. Cartesian k-space then started at the
  origin instead of its corner. Areas are now taken by differencing a cumulative integral, which is
  exact at any subdivision.
- **The reconstructed image came out transposed, and the off-resonance sign was inverted.** MRzero
  indexes its phantom `[x, y]` while every reconstruction convention treats rows as y, and
  `SpiralOffresonance` computed `exp(-2i pi (k.r - df t))` where the standard convention -- and what a
  scanner delivers -- has both terms negative. With the sign inverted the correction *adds* the phase
  it should remove. Both are fixed at the source: nothing at the call site transposes or negates
  anything. Together they made every number in an ADC table wrong without making any absurd.
- **The exact reconstruction operator does not scale.** 67 000 samples at 128 x 128 is over a billion
  complex exponentials per matrix-vector product. `SegmentedOffresonance` interpolates the
  off-resonance term between a few instants and leaves the rest to a NUFFT -- seconds instead of hours,
  agreeing with the exact operator to 0.011 where both can be run. `segments_for` derives the segment
  count from an explicit error bound (`dtheta^2 / 8`) rather than a round number of residual cycles,
  and the field's mean is removed exactly, since a uniform off-resonance is a function of time alone.
- **A ripple in a corrected image is more likely in the operator's inputs than in the operator.**
  Measured on this sequence: segmentation contributed 0.0004 of RMSE, nearest-neighbour resampling of
  the B0 map 0.0072, and conjugate-gradient iteration count 0.069 -- and the last with the *opposite*
  sign to intuition, because correcting worsens the conditioning and CG then fits the model's own
  mismatch. Which term dominates depends on the other two, so vary them together; one at a time finds
  whichever was looked at first.
- **`Readout` and `reconstruct_shot` are the whole reconstruction interface.** Four arrays via
  `Readout.from_sidecar`, and one call taking `(n_coils, n_samples)` however it was obtained, so
  simulated and twix data take the same path. `to_mr0` reads a written `.seq` back given its `System`,
  and `first_rep` selects one shot -- MRzero materialises a samples x voxels array, so two long shots
  against a fine phantom can ask for gigabytes.
- **Undersampling is paid for in coil sensitivity, not in algorithm.** One uniform channel at 4x has
  0.78 samples per unknown: RMSE 0.50 against 0.052 with twelve channels, on identical data.
- **ADC dwell does not shorten a spiral readout.** The gradient waveform sets the duration and the
  spiral is slew-limited, so dwell only decides how finely it is sampled -- 28.95 ms at 1, 2, 4 and
  8 us alike. The constraint on dwell is Nyquist *along the arm*: `gamma*|G|*dwell` must stay under
  `1/FOV`, which 8 us violates.
- **The six signed axes cannot fit a diffusion tensor.** Each direction contributes the row
  `[dx^2, dy^2, dz^2, 2dxdy, 2dxdz, 2dydz]`, and along an axis the last three are zero, so the
  off-diagonal elements are undetermined and the design matrix is singular. Three mutually
  perpendicular directions look ideally separated -- 90 degrees apart -- which is the case where the
  minimum pairwise angle is not merely a weak figure of merit but an actively misleading one.
  Notebook 2 acquires the axes *plus* the 12 face diagonals, at both polarities: condition number
  `sqrt(2)`, better than a repulsion-optimised set of the same size. Direction schemes live in the
  notebook, not the package -- `MonopolarDiffusion` takes a unit vector and produces gradients.
- **`synthetic_hardware()` is not a stand-in for a PNS verdict.** On the DTI sequence it reports
  2.44x the stimulation limit while the real Cima.X descriptor reports **0.95x** -- the difference
  between "not runnable" and "passes". Acting on the synthetic number would mean derating the
  diffusion lobes and lengthening TE to fix a problem the scanner does not have. Notebook 2 now loads
  the vendor `.asc` when `SEQCRAFT_ASC_DIR` points at one and labels the synthetic result as unfit for
  judging runnability. The file is never copied into the repository; only its name and sha256 are
  recorded.
- **PNS attribution is not additive, and the two models disagree about which term dominates.** With
  everything at slew 0.65 the peak sits in the spiral; derate the readout alone and it *moves* onto
  the refocusing crushers, so the limiting gradient changes as you tune. Against the real descriptor,
  spoiling on three axes instead of one takes the peak from 0.95 to 1.50 -- from passing to over --
  while sampling density changes it by 0.03 across a 4x difference in readout length. The synthetic
  model ranked those two the other way round. Both facts are measurements, and neither is guessable.
- **The vector-norm warning is not a proxy for peripheral nerve stimulation.** In the DTI example it
  pointed at the three-axis spoiler while the actual PNS was dominated by the spiral: dropping the
  spoiler to one axis changed PNS by 3 %, while derating the slew from 0.65 to 0.15 was what cleared
  it. Measure PNS with `CompiledSequence.pns()` against your own hardware model.
- Simulation and reconstruction helpers live in `examples/lib/`, **not** in the package.

### Tests

544 tests and doctests, over `tests/logic`, `tests/compiler`, `tests/modules` and
`tests/integration`. `tests/compiler` is the heart of it: one case per overlap rule, plus the
adversarial ones — a gradient straddling an RF, a boundary that would fall inside an ADC window, a
split mid-ramp, two RFs whose dead times overlap but whose waveforms do not.

---

## 0.2.1 and earlier

Superseded by 0.3.0. The `SeqModule` / `Sequence` / segments architecture and its 547 tests were
replaced wholesale; nothing from them is importable.
