# Changelog

## Unreleased — `pp.Opts` is the scanner, and there is no module library

Three concepts are removed, and none could go alone: they held each other up. Full rationale in
[ADR-003](docs/adr/003-scanner-and-module-reform.md).

### Removed — `System`, `Limits`, presets and named regimes (`core/system.py`, 652 lines)

The scanner is described by the **official `pypulseq.Opts`** and nothing else. The compiler reads
eight fields of it and touched `System` at three call sites; everything `System` stored was already
an `Opts` field.

```python
out = sc.compile(tree, system, regime='epi')     # was
out = sc.compile(tree, opts)                     # is
```

- `CompiledSequence` stores `opts` instead of `system` + `regime`.
- `pns(hardware)` now **requires** the model: PNS prediction is analysis, not compilation, and a
  response model has nothing to do with a limit set.
- The provenance sidecar records `vars(opts)` under `opts`, replacing `system` and `regime`.
- Named regimes become a second `Opts`: `sc.opts.derate(opts, grad=0.85)`. It copies every field it
  was not asked to change, which is what the multi-regime consistency check used to guarantee.
  **Do not hand-write `pp.Opts(max_grad=...)` to derate** — `Opts` fills every omitted argument from
  the *process-global* default, so the hand-written version silently returns your dead times to zero.
- **No compatibility shim.** A `System` forwarding to `Opts` would preserve exactly the concept being
  removed.

`sc.opts` holds the two operations `Opts` makes awkward or unsafe — `derate`, and `from_scanner`,
which looks a scanner up in [PulseqSystems](https://github.com/nimpulseq/PulseqSystems) (optional
extra `seqcraft[systems]`). There is deliberately **no wrapper around the `Opts` constructor**:
build one the ordinary way. `from_scanner` takes `rf_dead_time`, `rf_ringdown_time`, `adc_dead_time`
and `max_b1` as *required* keyword arguments, because a vendor database cannot supply them and
pypulseq defaults the first three to **zero** — a sequence built on those compiles cleanly,
validates cleanly, and is refused or silently mangled at the console.

`load_hardware` / `synthetic_hardware` moved to `sc.hardware`, out of `core`. `load_hardware` now
returns just the model, with its provenance string on `.source`; the acoustic-resonance bands it
used to return had no consumer anywhere.

### Removed — the module library (`seqcraft/modules/`, 27 classes, 5 762 lines)

**seqcraft ships no concrete modules.** `sc.modules.SincExcitation`, `EPIReadout`, `SpiralVDS`,
`MonopolarDiffusion` and the rest are gone, along with the flat `sc.*` re-exports.

The physics worth keeping was **lifted out as plain functions** into `salvage/` before the deletion,
with no scanner object and no base class in their signatures: the exact b-value integral and its two
solvers, the variable-density spiral trajectory, the EPI ramp-sampling moment integral, and the DTI
direction tables.

### Added — `sc.Module`, the contract, at `src/seqcraft/module.py`

Parameters in, one `LogicBlock` out. Four members, and each earns its place:

```python
class PhaseEncode(sc.Module):
    def __init__(self, *, opts, fov_mm, matrix, tag=None):
        super().__init__(opts=opts, tag=tag)
        ...
    def build(self, *, line=0) -> sc.LogicBlock:
        ...

pe = PhaseEncode(opts=opts, fov_mm=250, matrix=64)
tr.add(t, pe(line=17))                       # __call__ is the interface; build is what you write
```

- `opts` is **required**, because pypulseq's fallback is the process-global `Opts.default`.
- `build` is **abstract** — a subclass that does not produce a block is not a module. The old base
  declared no abstract method and so guaranteed nothing.
- A `build` returning a non-block raises `TypeError` **naming your class**, rather than failing
  hundreds of lines later inside `add`.
- It lives beside `core`, not inside it. A test asserts `seqcraft.core` never imports it.

### Changed — a module declares no duration, so a build argument may change one

The rule is **inverted**. There is no `duration` property, so there is no second source of truth for
a build argument to invalidate. Build the block, then place by it:

```python
exc_block = exc()
tr.add(0.0, exc_block)
tr.add(exc_block.duration + gap, readout(line=k))
```

Building is cheap — the design happened in `__init__` — so nothing is lost. `sc.testing`'s
`assert_duration_is_honest` and `assert_timing_properties_in_range` are deleted with the property
they policed.

### Changed — `sc.testing` and `sc.display` take an `Opts`

`assert_raster`, `assert_within_limits`, `assert_compiles` and `assert_output` take `opts` and lose
every `regime=` parameter. `assert_all(module, **args)` calls `module(**args)` and reads
`module.opts`. `plot_block(root, opts)` likewise.

`plot_trajectory` now takes **an iterable of `(kx, ky)` array pairs** instead of an object with
`.trajectory()`, `.n_interleaves` and `.k_max_per_m` — a duck type only one class in one library ever
satisfied, which quietly made a plotting helper depend on that library.

### Removed — `seqcraft.ordering`, and five unused validators

Four of `ordering`'s six functions had never had a caller, and ordering tables are
sequence-programming choices rather than physics. The module moves to `salvage/`; the two the tests
used are three lines each and are written out where they are needed.

`validate.py` loses `check_units`, `require_positive`, `require_int_in`, `require_divides` and
`suggest_field` — the deleted base was their only consumer. `merge_definitions` (the compile path)
and the plausibility bands `geometry` uses remain.

### Changed — every test fixture is raw pypulseq

`tests/compiler/test_fidelity.py` built its realistic trees out of library classes, which made
compiler coverage depend on whatever the library contained. All fixtures are now `pp.make_*` calls,
and that is a standing rule rather than a state of the tree. The integration tier keeps a gradient
echo and a spin echo, both raw; the DTI and EPI-DWI tiers went with the library they were built on.

The compiler baseline was re-captured, since the four module-built recipes it froze no longer exist.

### Changed — examples

`01_getting_started.ipynb` is rewritten and **uses no modules at all** — raw pypulseq events into a
`LogicBlock`, then `sc.compile`. It doubles as the proof that the core stands alone. `sc.Module`
appears at the end, once there is a reason for one.

`dti_spiral/` and `dti_epi/` are **parked** under `examples/_parked/`. They do not run against this
version. They are kept because they are the acceptance test for whatever module set is written next
— rewriting them is what should *drive* that library rather than follow it.

`examples/lib/` survives: `to_mr0` takes `opts=` instead of `system=`, and the reconstruction helpers
already took explicit arguments rather than reading attributes off module objects.

## Earlier unreleased work — a batch form for `LogicBlock.add`

### Added — `add` accepts a table of `[time, *items]` rows

The chained form is fine when the times are literals. It is the wrong shape when the schedule is
*computed*, because the rows already exist as data and the only way to hand them over was a loop
that called `add` once per row. The batch form makes the schedule a first-class value:

```python
seq.add([[t0, rf, gz],
         [t1, gzr]])                       # == seq.add(t0, rf, gz).add(t1, gzr)

plan = [[t, *events] for t, events in schedule]
seq.add(plan)                              # a computed score, in one call
```

- **Dispatch is on `list`/`tuple`, and nothing else.** Not "any iterable": a `LogicBlock` is itself
  iterable and must never be mistaken for a table of rows.
- **A bare row needs no outer brackets.** `add([t0, rf, gz])` — a first element that is a number
  means one row. `numbers.Real`, so numpy scalars are times without a cast, and `nodes` still stores
  plain floats.
- **Rows are never sorted by time.** They are appended in the order given and items within a row keep
  theirs, because insertion order is what `flatten` and the compiler's tie-breaking rely on.
- **`[t0]` with no items is a legal no-op**, consistent with `add(t0)`. So is an empty table.
- **Errors name the row**, so a thirty-row table points at the offending line rather than at `add`.
  A mixed call (`add([[...]], rf)`) is rejected: the batch form takes exactly one argument.
- `@overload` stubs declare both call shapes, so a type checker sees them and rejects the mixed call
  statically as well.

`nodes` keeps exactly the same shape, so the compiler, `flatten`, `display` and `testing` are
untouched. This is input-side sugar over an unchanged data model, and `LogicBlock` still has two
attributes and one method.

### Fixed — `add(block)` was a silent no-op

Forgetting the start time — `seq.add(exc.build())` rather than `seq.add(0.0, exc.build())` — fell
into the `*items` loop with nothing to iterate and returned `self` having added nothing, so the
event went missing with no error anywhere. An item as the first argument is now a
`ConfigurationError` naming the type and the two spellings that work.

## Unreleased — `LogicBlock` is the interface

`Module` was described as one of seqcraft's three concepts and lived in `core` beside the compiler,
which read as a requirement: to write a sequence component you inherited a base class and implemented
one abstract `build(**args) -> LogicBlock`. It never was a requirement — the compiler's input has
always been a `LogicBlock` and it has never asked what produced one — but the architecture said
otherwise, and the contract it imposed is wrong for components with more than one output. A diffusion
encoding's two lobes became `build(part='pre')`, a keyword standing in for two methods that should
have been named `pre` and `post`.

**The concepts are now two: `LogicBlock` and `compile`.** Modules are a provided library on top of
them.

### Changed — `Module` moved to `seqcraft.modules.base` and stopped being a contract

- **`core/module.py` is gone.** The class is `seqcraft.modules.base.Module`, re-exported as
  `seqcraft.modules.Module`. `core` now holds only what is on the path from a block to a validated
  `.seq`, which was always the stated membership rule.
- **No `abc.ABC`, no abstract `build`.** The base declares no abstract method, so nothing is required
  of a subclass. It keeps what a reusable module actually needs and nothing more: `system`, `regime`,
  the resolved `opts`, the unit check that runs when a subclass's `__init__` returns, `params()`,
  `submodules()` and `repr`. Existing `build()` methods are unchanged — the *requirement* is what was
  removed, not the method.
- **`sc.Module` and `sc.core.Module` are gone, with no shim.** Use `sc.modules.Module` or
  `from seqcraft.modules import Module`. A compatibility import would keep asserting the membership
  this release removes, and at 0.3.0 alpha the project has broken cleanly before rather than
  accumulate shims.
- `RFPulse` now declares `abc.ABC` itself. Its abstract `_design` is a real requirement — a pulse
  shape is the whole of what its subclasses add — where the `build` requirement never was.

A component may now be anything that returns a `LogicBlock`: a function, or a class with as many
domain-shaped methods as it likes (`readout.readout()` and `readout.prephaser()`,
`diffusion.pre()` and `diffusion.post()`). `tests/modules/test_module_base.py` holds that line,
including a test that `core` never imports `seqcraft.modules`.

### Changed — `seqcraft.testing` splits into block-level and convention-level

The assertions took a module and called `.build()`, so testing a component that had neither meant
not testing it.

- **`assert_output(make, system, *, regime=)`** is the new general entry point: it takes any callable
  returning a block and runs the whole block-level suite on it. **`assert_block`**, **`assert_raster`**,
  **`assert_within_limits`** and **`assert_compiles`** now take `(block, system)` rather than a module,
  and **`assert_deterministic`** takes the callable.
- **`assert_pure(component, make)`** and **`assert_duration_is_honest(component, make)`** take the call
  under test; **`assert_timing_properties_in_range(component)`** needs no call at all.
- **`assert_all(module, **build_args)`** is unchanged in signature and is now a wrapper over those. It
  reads `build`, `system`, `regime` and `duration` off the object and never checks its type, so
  inheritance is not required.
- **`all_modules()` → `module_subclasses()`**, renamed because the old name claimed a universe it
  never described. It enumerates what inherits the optional base — the right question for
  parametrising the library's contract suite, the wrong one for deciding what seqcraft accepts.
  Private bases (`_AreaTrapezoid`) are skipped by name now that they are no longer abstract, and the
  library's own coverage assertion filters by package so a user's subclass cannot contaminate it.

### Fixed — `assert_pure` could not see the bug it was written for

It hashed the stored events, called the builder **twice**, and compared. The canonical mutation it
exists to catch — the reference implementation's `self.gx.amplitude = -self.gx.amplitude` inside a
readout loop — is an involution, so two calls left every hash where it started and the check passed.
Now checked after each call.

### Changed — `assert_deterministic` walks the whole tree

It compared a block's direct children, so a component that nests — `FatSat`, `EPIReadout` — had its
actual events skipped and passed vacuously. The same fix `assert_within_limits` and `assert_raster`
got in the EPI release, applied to the one that was missed.

## Unreleased — EPI

Plan: [PLAN_EPI_V1.md](PLAN_EPI_V1.md), which is phase 1 of
[PLAN_MODULES_V1.md](PLAN_MODULES_V1.md). Every item below has a test that fails without it.

### Added — `EPIReadout`

An echo-planar train with its own prephasers, blips and labels. Ramp sampling, partial Fourier along
the phase encode, partial echo along the readout, interleaved segmentation, blip-up/blip-down, and
`LIN`/`SEG`/`REV`. `modules/readout/epi.py`; nothing in `core` changed to accommodate it.

**The train is one gradient event per axis, not one per echo.** Per-echo trapezoids put the tail of
echo *n* and the head of *n+1* in the same compiled block, and the compiler warns whenever a block
holds two gradients on an axis — **96 warnings per shot**, 1824 for a 19-volume acquisition, every one
describing correct output. One `make_extended_trapezoid` per axis gives each block one piece to split,
and the compiled waveform is identical.

**Echo spacing is slew-limited, and the amplitude limit never binds.** For a lobe of fixed area the
minimum-time trapezoid is a triangle at `G = √(A·S)` — 48.9 mT/m against a 170 mT/m limit at 240 mm on
a 128 matrix. So there is no flat top to sample, ramp sampling is not optional, and the only levers on
echo spacing are slew, `k_max` and partial echo. `flat_time_us` is an override for deliberate
derating, not the parameter. A test asserts that doubling `max_slew` shortens the echo by √2 while
doubling `max_grad` changes nothing.

**What must cover k-space is the *sampled* extent, not the lobe's area.** The blip is centred on the
junction where the readout gradient crosses zero and the ADC skips `blip/2` at each end, so `ky` is
constant while a line is read — measured drift **0.000000 1/m**. The lobe therefore carries 540.53 1/m
to deliver 533.33 1/m of sampled extent, and the prephaser cancels the difference exactly. Minus
`k_max` is the tempting value and would displace every shot by 1.7 `dk`; a k-space offset is a linear
phase ramp across the image, so the magnitude image looks perfect and every phase-derived quantity is
wrong. This is `CartesianLine`'s documented ramp-area trap in a second form.

**`time_to_echo` is neither the midpoint nor the first sample.** With partial Fourier 0.75 the `ky = 0`
echo is number 32 of 96 — 17.26 ms into a 49.92 ms train. Taking the midpoint puts TE 7.7 ms late,
silently; taking the first sample, as a spiral would, is 17 ms out. It is found by inverting the
lobe's moment, and it is `p // n_shots`, which is **shot-independent** — required, since a build
argument may not change a timing property, and also the physics: segmented EPI plays every shot at the
same TE so they carry the same T2 weighting.

**Three rasters constrain the ADC, and each was found by being caught.** The blip is an even number of
gradient rasters, so `blip/2` lands on one — a 50 µs blip failed with *"The last time point must be on
a gradient raster"*. The ADC's node goes at the lobe start so it lands on the **block** raster: a node
off it is snapped with a `raster` warning, which moved sampling against the gradient by up to half a
raster and put `kx` at ±268.66 where ±266.67 was asked for. And the sampling offset lives in the ADC
event's **own delay**, on the **RF** raster, which is what pypulseq's `check_timing` requires of an ADC
delay — 33.2 µs produced one error per echo. Carrying the offset in the node time instead made two
delays *add*, because `pp.make_adc` silently raises a delay below `adc_dead_time` up to it while
seqcraft preserves an event's own delay: sampling began 40 µs into the lobe rather than 30, leaving
`kx` 4.88 1/m asymmetric about the echo with every limit check passing.

**The sampled window must be centred on the lobe, exactly.** Otherwise the leading and trailing
unsampled corners differ, they stop cancelling between a forward echo and the reversed one after it,
and the reversed echoes' `kx` drifts. Measured at **800 ns of asymmetry on a two-shot design: 0.28 1/m,
6.8 % of `dk`** — while a single shot was exact, because its leftover happened to split evenly. An
alternating k-space offset is the worst kind: a phase ramp of alternating sign is a Nyquist ghost, and
it looks exactly like an uncalibrated gradient delay rather than an arithmetic mistake. The sample
count is now chosen so the leftover splits onto the RF raster, and the trajectory agrees with
`calculate_kspacePP` to **0.0000 1/m** for 1, 2 and 4 shots.

**The dwell is the longest that satisfies Nyquist at the apex**, which is where `k` moves fastest —
2000 ns and 228 samples per echo here, against a 2000 ns bound. Maximising the *sampled span* instead,
so the gradient is as low as possible, is measurably the wrong trade: it picks 1000 ns and 460 samples
to buy 4 µs more span and **0.2 %** less gradient, doubling the data the reconstruction and the
simulator both scale with.

`EPIBlip` and `ReversedPolarityPair`, which PLAN_MODULES_V1 listed alongside, are **not** included. A
blip is one trapezoid of a chosen area, which is `PhaseEncode` — and the train cannot delegate it
anyway, since the blip's duration sets the unsampled corner, which sets `G`, which sets the echo
spacing. A blip-up/blip-down pair is one train played twice with the phase encode reversed, which
changes no duration and is therefore `build(pe_polarity=-1)`. Flyback EPI is out of scope for v1, and
the docstring says why.

### Fixed — `check_limits` reported 2441 % where the truth was 94 %

Three independent errors in `events.check_limits`, all reachable only through a gradient whose knots
are neither uniform nor on the raster — which is what an EPI train has to be:

- **The vector-norm slew divided `diff(norm)` by the raster** rather than by `diff(t)`. An EPI echo
  junction's knots are 260 µs apart, so dividing by 10 µs overstated the slew by 26×, and 4882.9 T/m/s
  is exactly the peak amplitude over one raster. The per-axis path divided correctly, which is why this
  surfaced only as a warning — 97 of them on a 96-echo train, which is the state in which a warning
  stops meaning anything.
- **The axes were stacked by sample index**, so `gx` at 260 µs was combined with `gy` at 490 µs, and
  the shorter axis was zero-padded rather than held at its `last`. That hides a real violation as
  easily as it invents one, and the test demonstrates the hiding direction: two axes peaking at the
  same instant reach 106 % of the amplitude limit and were reported as 84 %.
- **A second event on an axis overwrote the first** in the norm's dictionary, so only the last was ever
  combined.

Everything is now measured on the union of the events' knots, where a piecewise-linear function's
slope is piecewise constant — the same fact the compiler's `_superpose` rests on. `check_limits` gained
an optional `starts=`, needed when two events on one axis are passed without their node times;
`testing.assert_within_limits` passes them, so two lobes of a bipolar pair no longer land on top of
each other and sum to zero.

**`waveform_of` promised a uniform raster and did not deliver one**, which is the root cause: for
`type='grad'` it returns `event.tt` verbatim and ignores `raster`. True for `make_arbitrary_grad`
(raster centres), false for `make_extended_trapezoid`. Now documented as what it is — a curve to plot
or interpolate, carrying its own time axis — with `knots_of` and `pwl_moment` named as the exact path.

### Changed — the contract suite walks the whole tree

`testing.assert_within_limits` and `assert_raster` iterated a built block's **direct children**, so a
module that nests — `FatSat`, and now `EPIReadout` — had its actual gradients skipped and passed
vacuously. Both use `flatten` now, and `assert_raster` checks gradients only, since an RF or ADC
carries its dead time in its own delay and answers to the RF raster.

### Changed — `examples/` is one folder per scan

```
examples/  01_getting_started.ipynb   lib/   dti_spiral/   dti_epi/
```

Each folder holds `01_build.ipynb`, `02_simulate_and_reconstruct.ipynb` and its own `seq/`, and is
runnable from inside itself with nothing above it. `spiral_recon.py` → **`noncartesian_recon.py`** and
`SpiralOffresonance` → **`NoncartesianOffresonance`**, wholesale and with no shim, because the module
serves both readouts: a spiral's blur and an EPI's phase-encode shift are the same `exp(-2πi·df·t)`
term, and correcting either needs no new operator. The field-map section is duplicated between the two
build notebooks; `examples/README.md` records that as a decision and what it buys.

### Fixed — three things `examples/lib` had inferred from a spiral

Each was correct for a spiral for one reason — a spiral's first sample is its echo — and wrong for any
readout whose echo falls mid-train.

- **`Readout` rebased the sample times on `min(t_adc_s)`.** For an EPI train that puts
  `2π·df·17 ms` of phase into the operator, and since `df` varies with position it is an image artefact
  rather than a global constant — up to 2.9 cycles over a 172 Hz field range. The echo reference is now
  an explicit `t_echo_s`, read from the sidecar, and a sidecar whose times look *absolute* is refused
  outright rather than reconstructed.
- **`density_compensation` weights by `|k|`.** A ramp-sampled EPI line crosses `k = 0` at full
  gradient, which is where its samples are **sparsest**, so `|k|` weights this trajectory backwards.
  `speed_compensation` computes `|dk/dt|`, the correct Jacobian for both readouts;
  `reconstruct(dcf=...)` selects, and the default is unchanged so the spiral notebook does not move.

  Measured, the cost is not where one would guess, and the notebook shows both cases rather than
  asserting one: for a **single density-compensated back-projection** — what gridding does — the
  Jacobian has to be right, because nothing downstream repairs it. For a **converged
  conjugate-gradient solve** the weights are a preconditioner on a problem whose solution they do not
  change, so a wrong one costs iterations rather than accuracy. That distinction is why "the density
  compensation barely matters" and "the density compensation is critical" are both repeated in the
  literature, each true of one case.
- **`to_mr0` rebased each repetition's k on its own first ADC sample.** An EPI readout's prephaser sits
  inside the same repetition as the readout, so that subtracts the prephaser away and reports a
  trajectory starting at the origin instead of at the corner of k-space. `k_reference='excitation'`
  accumulates from the excitation and flips the sign at each refocusing pulse — the physical
  convention, and what `calculate_kspacePP` computes. No rebasing is then needed for any readout, and
  no constant offset could have substituted for the missing sign flip: a b=1000 encoding's lobes are
  many times `k_max`, and with the same sign they swamp the readout instead of cancelling.

### Fixed — the MRzero bridge placed every ADC sample half a dwell late

**This one produced a Nyquist ghost, and it was found by looking at the reconstructed image.** Three
faults, compounding, all of them invisible on a spiral:

- **Each sampled event ended on a dwell *boundary* rather than at its sample's centre.** MRzero applies
  an event's gradient moment and *then* records, so where the event ends is where the sample sits in
  k-space — half a dwell late, here **2.14 1/m, 0.514 `dk`**. A spiral never reverses its readout, so
  this is a small radial offset that merely blurs. An EPI train reverses every echo, so the offset
  **alternates sign with echo parity** — a phase ramp of alternating sign, which is a Nyquist ghost at
  half the FOV, and which looks exactly like an uncalibrated gradient delay rather than an arithmetic
  mistake. Events now end at sample centres.
- **The lead and tail areas of an ADC block were handed to the simulator but left out of the reported
  k.** On a spiral the lead is 10 µs of almost nothing. On an EPI train the lead and tail are *exactly*
  where the phase-encode blips sit, so the reported `ky` never advanced: **95 `dk` short over a
  96-echo train**, while the simulated signal was correct throughout.
- **`_areas` interpolated a cumulative integral linearly**, having built it by the rectangle rule on
  raster-centre samples. Between raster points that integral is a quadratic, so it under-read on a
  ramp — and an ADC's sample centres fall precisely between raster points. Another 0.024 `dk` of
  alternating error. It now integrates the gradient's **exact piecewise-linear knots**, the same
  primitive the compiler uses.

Together: the simulator's trajectory and the module's now agree to **0.0001 `dk`**, from 0.514.
`tests/examples/test_bridge_areas.py` covers the area arithmetic and the sample placement against
`calculate_kspacePP`, including an explicit assertion that the *alternating* component is zero — a
constant k offset is a linear phase and nothing worse, and only the alternating one makes stripes.
The area tests need no simulator, so they run in the ordinary suite.

### Corrected — what a wrong echo reference actually costs

The plan claimed a first-sample time reference puts a spatially varying phase into the operator and
therefore an image artefact. The first half is right and the conclusion was not, and the notebook now
does the algebra: shifting every sample time by `Δ` multiplies the operator by
`diag(exp(-2πi·df(r)·Δ))`, so the solution comes back multiplied by `exp(+2πi·df(r)·Δ)` — a spatially
varying phase on the **image**, with the magnitude untouched exactly. Measured: magnitude RMSE 0.000000,
phase error peak-to-peak as predicted by `2π·df·Δ`. It costs a magnitude image and an ADC map nothing,
and it is fatal for a field map, for flow, and for the phase navigators multi-shot DWI needs.

### Added — `examples/dti_epi/`

The same diffusion encoding as the spiral example — 3 axes plus 12 face diagonals, condition number
√2 — through single-shot and two-shot ramp-sampled EPI at 1.88 mm with partial Fourier 0.75. Measured
on a Cima.X at 240 mm:

| | single-shot | two-shot |
|---|---|---|
| echo spacing | 520 µs | 520 µs |
| echoes, train | 96, 49.92 ms | 48, 24.96 ms |
| `ky = 0` at | echo 32, 17.26 ms in | echo 16, 8.94 ms in |
| TE | 59.40 ms | 42.80 ms |
| PE bandwidth per pixel | **20.03 Hz** → 5.0 px per 100 Hz | 40.06 Hz → 2.5 px |
| samples per shot | 21 888 | 10 560 |
| PNS, synthetic model | 276 % | — |

**Derating the readout does not fix PNS on its own.** From 1.00× to 0.30× readout slew the synthetic
figure falls 276 % → 175 % while the train grows 49.9 → 90.2 ms: the readout dominates, but the
diffusion lobes on a face diagonal put the whole gradient vector on two axes and set a floor of their
own. The notebook prints the whole trade and refuses to call any row runnable — on the spiral sequence
the synthetic model said 2.44× where the real descriptor said 0.95×, and the error is neither small nor
in a known direction.

The reconstruction notebook plays all three files against a phantom with a 217 Hz field — chosen so the
displacement is 7 pixels rather than the spiral example's 3, which is not enough to see — and goes
through to an ADC map. Against a true `D` of 1.50 / 0.38 / 3.75 e-3 mm²/s it recovers **1.50 / 0.37 /
3.74** single-shot and 1.50 / 0.37 / 3.77 two-shot. **No unwarping step appears anywhere**: the field
map goes into the operator and the solve removes the distortion, which is the same operator and the
same code that deblurs the spiral.

Three measurements in it are worth keeping, because each contradicts a reasonable guess:

- **Off-resonance correction barely moves the ADC where signal is plentiful.** An ADC is a ratio of two
  images with the same readout, so the distortion largely divides out: corrected and uncorrected agree
  to 0.01 e-3 for the background. The exception is free water — 3.74 against 3.46 uncorrected — where
  `exp(-b·D)` leaves 2.4 % of the signal and a residual is no longer small against it. So the
  correction buys least where there is most signal, and what it always buys is geometric fidelity,
  which no table shows.
- **A radial density compensation is worse than none at all on this trajectory.** Relative RMSE for a
  single back-projection: `|dk/dt|` **0.052**, no weighting 0.093, `|k|` **0.906** — seventeen times
  the true Jacobian's error, and still four times it after 30 iterations. An EPI line is close to
  uniformly spaced in k except at its two ends, so the radial argument is not merely suboptimal.
- **A wrong echo reference is a pure image-phase error.** Predicted 3.66 cycles peak-to-peak from
  `2π·df·Δ`; observed residual against that prediction **0.0000 cycles RMS**, magnitude RMSE 0.000030.

### Performance — the EPI reconstruction notebook, 433 s to 301 s

Profiled rather than guessed at. `simulate` is 3.6 s a shot and `to_mr0` 0.3 s, so the simulation was
never the cost: **a corrected solve is 37 s**, because `segments_for` returns **88** here against 13
for the spiral. The segmented operator does `segments × n_coils` oversampled FFTs per application and
the segment count scales with `field range × readout span` — EPI loses on both, at 49.9 ms of train
against 17.1 ms and a field deliberately widened to 217 Hz so the displacement is 7 pixels rather than
3. That much is inherent to the demonstration.

What was not inherent, and is fixed:

- **The density-compensation comparison ran through the corrected operator** — five solves including a
  150-iteration reference, ~300 s of the total. A density compensation is a property of the
  *trajectory*; off-resonance has nothing to do with it. At `b0 = 0` the operator needs the 8-segment
  floor, which is 8× faster **and** the correct isolation: with the field map in, a reader could
  attribute one effect to the other. The reference is now 200 iterations because it became cheap.
- **A duplicate simulation** — `acquire` was called a second time purely to fetch trajectory metadata
  that the first call had already returned.
- **A duplicate solve** — the `speed` variant recomputed what section 4 had already produced.

The notebook now states the measured costs and the order the levers are worth pulling in: do not
correct what you are not studying, do not solve twice, and only then reduce iterations — which is the
only one of the three that trades away accuracy.

Its committed size also came down from 608 kB to 277 kB: the figure was rendering 128² images at
100 dpi, upsampling them three and a half times, and the repository's own `check-added-large-files`
backstop is 512 kB.

### Changed — the spiral example, which the bridge fix improved

`t_echo_s` is written into the `.traj.npz` rather than left to be inferred, which for a spiral is
`adc.delay` and not `min(t_adc_s)` — a half-dwell difference, and the reason the inference existed at
all.

Re-executed, the spiral's **field map now recovers to 0.95 Hz RMS against the truth where 0.3.0
recorded 2.0 Hz**, and its ADC reads 1.50 / 0.38 / 3.66 against a true 1.50 / 0.38 / 3.75. Nothing in
the sequence changed: the improvement is the sample-placement fix above, which a spiral suffered from
too — just as a blur rather than as a ghost, which is why it had gone unnoticed. The inline note in the
build notebook that quoted the old figures for one- against three-axis spoiling has been rewritten to
point at what the reconstruction notebook prints, rather than carrying numbers that no longer hold.

### Tests

768 tests and doctests, up from 671. `tests/modules/test_epi.py`, `tests/core/test_events.py` and
`tests/examples/test_bridge_areas.py` are new; `tests/integration` gained an SE-EPI DWI built from the
same helpers as the spiral DTI, including the assertion that the two agree on
`DiffusionLobeDuration`, `DiffusionBigDelta` and `AchievedbValues` while differing on TE.

Two of the EPI tests are parametrised over shot count **because a single shot passed while the design
was wrong**: the asymmetric-window bug left `kx` exact for one shot and 6.8 % of `dk` out for two.
Where a bug is invisible in the default configuration, the parametrisation is the test.

## Unreleased — core revision (v3): less is more

Plan: [PLAN_CORE_V3.md](PLAN_CORE_V3.md). The compiled output does not move — the integration suite
asserts physics and byte-identical rewrites, and it is unchanged.

### Changed — `core/units.py` rewritten as one function

- **`convert(value, from_unit, to_unit=None, *, gamma=, f0=)`** replaces fifteen one-argument helpers
  (`mm`, `us`, `deg`, `mT_per_m`, `s_per_mm2`, …). Each of those encoded one conversion in one
  direction from one unit, so µs → ms, Hz/m → G/cm and ppm → Hz were simply not expressible; nothing
  outside `core/` called any of them. The new signature is the one `pypulseq.convert` uses,
  generalised to eleven dimensions — time, length, angle, frequency/field, gradient, slew, k-space
  (and its rate and area), b-value, ratio — in both directions, over 51 units.
- **Field and frequency are one dimension**, because in pulseq they are: B1 is carried in hertz. So
  `convert(12, 'uT', 'Hz')`, `convert(3.0, 'T', 'MHz')` (Larmor) and `convert(-434, 'Hz', 'ppm')` are
  all the same call.
- **`ppm ↔ Hz` needs a Larmor frequency** and says so, rather than guessing. `System.convert` fills in
  the scanner's own `gamma` and `f0 = gamma·B0`, so a chemical shift or a B1 limit cannot silently pick
  up the proton value on a system configured for another nucleus.
- **Conversions are exact where decimals are exact.** Each unit's scale is a `Fraction`, so the ratio
  between two units is computed exactly and collapsed to a float once — and where the ratio is a
  reciprocal integer the single operation is a division by an exactly-representable integer. `4200 µs`
  is `0.0042 s`, not `0.004200000000000001 s`, which is the difference between a value that compares
  equal to its raster and one a `ceil` pushes up 10 µs.
- **The hand-rolled conversions in the module layer are gone.** `/ gamma * 1e3` appeared nine times
  across `diffusion.py`, `spiral.py` and `pulses.py`, `* 1e-6 * gamma * b0_T` twice, and `/ 1e3`,
  `/ 1e6`, `* 1e6` at every `duration_us` and `fov_mm` site. All of them are now `convert`.
- **One vocabulary.** `validate.DEFAULT_RANGES` names a unit per field-name suffix; a test asserts
  every one of those names — and every alias name — is a unit `convert` knows, so an error message
  cannot quote a unit the reader is unable to pass back in.

### Changed — `core/raster.py` → `core/timing.py`, and the raster is an object

- **`Raster`** carries `dt`, a name, and the seven operations: `holds`, `ceil`, `floor`, `nearest`,
  `count`, `at`, `require`. Six free functions taking a bare `raster: float` are gone, and with them
  the call-site shape `ceil_to(self.duration_ms / 1e3, self.system.block_raster_s)` — two conversions,
  one of them a magic number, and a raster passed positionally. It is now
  `self.system.block_raster.ceil(convert(self.duration_ms, 'ms', 's'))`.
- **`System.grad_raster / rf_raster / adc_raster / block_raster`** return `Raster` objects and replace
  the `*_raster_s` floats. `.dt` is the float when array arithmetic needs one. One spelling per raster.
- **Nothing assumes 10 µs.** Rasters are values the scanner supplies; the doctests and tests exercise
  GE's 4 µs, Philips' 6.4 µs and two invented values, and a 4 µs/2 µs system is compiled end to end.
  The only floor is one tick (1 ps), 10⁵ finer than pulseq's finest raster, and a finer raster raises
  with that sentence in the message.
- **`picoseconds()` / `seconds()` → `to_ticks()` / `from_ticks()` / `TICKS_PER_SECOND`.** They read as
  unit conversions and were not; ticks are the exact-integer time domain, and the docstring now says so
  and points at `units.convert` for the actual conversion. `sum_exact`/`sub_exact` are `exact_sum`/
  `exact_diff`.
- **Fixed: `nearest` rounded negative times away from zero.** `-1.4` rasters became `-2`, not `-1`; the
  sign branch it came from is gone. Unreachable from the compiler, which rejects negative times before
  rounding, but wrong.

### Removed — the module registry

- **`core/registry.py`, 34 `@register()` decorators and `register`/`registered`/`lookup` are deleted.**
  The registry's only consumer was the contract suite's `parametrize`; `lookup` had no callers outside
  a docstring. It accelerated nothing and could not: it is a dict filled at import time, and no code
  path is shortened by it. A registry earns its place when a *string* must become a *class* at run
  time — a YAML front end, plugin entry points, a `--readout=` flag — and seqcraft has none of those on
  purpose. Meanwhile it created the failure it existed to prevent: a new module that forgot the
  decorator silently lost the entire contract suite.
- **`seqcraft.testing.all_modules()`** replaces it, walking `Module.__subclasses__()`. Subclassing *is*
  the registration, so it cannot be forgotten. `RFPulse` and `RefocusingPulse` now declare `_design`
  abstract, which is both honest and how discovery skips them.

### Changed — `core` holds only what compiles a sequence

`core` goes from 15 modules to 11. Public import paths via `seqcraft` are unchanged (`sc.ordering`,
`sc.plot_block`, `sc.testing`); `from seqcraft.core.X import …` changes for the three that moved.

- `core/ordering.py` → **`seqcraft/ordering.py`** — view orders, golden angle and RF-spoil phase are
  sequence-programming vocabulary, not infrastructure; the compiler and the data model never
  reference them.
- `core/provenance.py` → **`seqcraft/provenance.py`** and `core/display.py` → **`seqcraft/display.py`**
  — output tooling, neither on the path from a logic block to a `.seq`. `sc.provenance` is now a lazy
  attribute alongside `sc.display`. Two latent import bugs in `provenance.py` are fixed along the way:
  a `TYPE_CHECKING` import of a module that does not exist, and a missing `Mapping` import.

### Added

- `tests/core/test_units.py` — known equivalences, exhaustive round-trips over every pair in every
  dimension, exactness assertions, gamma cancellation within the tesla family, and the error messages.
- `tests/core/test_timing.py` — the float traps as tests, across eight vendor and invented rasters.
  Includes the measured one: accumulating an exactly-legal 1.5 ms TR leaves the 10 µs raster after
  9813 repetitions, 14.7 s in, after which more than half of all later start times are illegal.

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
