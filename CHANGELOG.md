# Changelog

## Unreleased — flow compensation and velocity encoding as physical intent

Two new public names, `sc.FlowCompensation` and `sc.VelocityEncoding`, and one new keyword each on
`GRE2DTR` and `GRE3DTR`:

```python
tr = sc.modules.GRE3DTR(
    ..., flow_comp=sc.FlowCompensation(axis='y'),
    velocity_encode=sc.VelocityEncoding(venc_m_s=1.5, axis='z'),
)
```

They are **values, not modules**: no `opts`, nothing built, and constructing one commits to
nothing. A caller says what physics they want; the repetition answers whether and how.

**Which axes work is the repetition's answer.** `GRE2DTR` carries a first-moment requirement on
`y`, its own phase-encode winder, and on `x`, where `CartesianLine` solves it locally with the
analytic two-lobe prephaser. Its `z` gradient is the slice rephaser, which `Excitation` owns, so a
claim there is refused before anything is designed — naming what that axis carries rather than
saying "unsupported". `GRE3DTR` owns `y` and `z`, because its z winder already carries the
partition encode. The same `FlowCompensation(axis='z')` value is honoured by one and refused by the
other, which is why that check cannot live on the value.

The two routes are invisible from outside. `flow_comp=sc.FlowCompensation(axis=('x', 'y'))` sends
one axis to a leaf's closed-form solve and the other to the repetition designer, and the emitted
repetition has both moments nulled.

**Flow compensation is the common mode, not a per-state zero.** Alone, the common mode is the one
acquired state and it reads as `m1 = 0`. Beside a velocity encoding it fixes the *mean* of the two
states while the encoding fixes their difference, so the pair comes out at `±Δm1/2` — derived, not
specified. Nothing computes that half; it falls out of two independent constraints on the same
moment, which is what lets the two compose without a precedence rule.

**Velocity encoding separates what a state means from when it is acquired.** The intent owns the
VENC relation — read from `VelocityEncode.delta_m1_for`, still the single source — and the identity
of the two states; the caller owns whether both are acquired and in what order; the repetition
realises the one it is handed:

```python
for state in venc.states:
    scan.add(at, tr(line=line, partition=partition, encoding_state=state))
```

A state is required exactly when velocity encoding is present, and refused when it is not. An
ignored `encoding_state` would hand a caller two identical repetitions and show up as a
subtraction that comes back zero, so both directions raise.

`sc.modules.VelocityEncode` is unchanged and stays public: it is the standalone bipolar, and
constructing one is not a prerequisite for the joint path.

Internally: the claim vocabulary — `CommonModeClaim`, `DifferenceClaim`, `JointProblem`,
`Schedule`, the realisation families — stays private, and public intent is translated to it in one
place. The value-domain checks moved to `design/validation.py`, below `modules`, which is what
keeps the intents from closing an `augmentation -> modules -> kernel -> augmentation` import cycle.

## Unreleased — the readout nulls its own first moment

`CartesianLine` can now null the **first** gradient moment at the echo on its own axis, as well
as the zeroth. The zeroth is what a prephaser already cancels, so a stationary spin picks up no
phase from that axis; the first is not zero, and a spin moving along the readout arrives with
`2 pi m1 v`. The capability reshapes the prephaser into two lobes of opposite sign so that both
moments are zero at the echo.

**No new public parameter.** The capability is reached through the repetition-level augmentation
API, which is where a caller says what physics they want and the repetition decides which of its
components owns the problem. A per-leaf keyword would be the first instance of a
leaf-by-augmentation surface the design deliberately avoids, so the switch that selects it is
private.

**Not a new module.** The fine scan recorded flow compensation as `EXTEND_EXISTING`, because
gradient moment nulling is a modification of a waveform that already belongs to someone — and
here that is `CartesianLine`, which owns both the winder and the readout lobe and can integrate
their moments. No moment-requirement interface, no solver, no compiler or `LogicBlock` change, and
nothing outside the readout is read or rewritten.

The solve is closed form. Taking the echo as the origin, the readout's own contribution up to it
is a property of the lobe alone, computed once from its knots. Two winders of equal duration `D`
played back to back have centroids at known offsets, and a symmetric trapezoid's first moment
about an external origin is its area times its centroid exactly, so both conditions are **linear**
in the two areas:

```text
A1 + A2                            = -m0_ro
A1 (t_rs - 3D/2) + A2 (t_rs - D/2) = -m1_ro
```

`D` is found by bisection on the gradient raster, and the monotonicity that licenses bisection is
derived rather than assumed: before the first echo the readout has one sign, so
`C = m0_ro e - M` is the integral of `g(t) t` over an interval where both factors are
non-negative, which makes `A1 = C/D + |S|/2 > 0` and `A2 = -(3|S|/2 + C/D) < 0` with both
magnitudes falling as `D` grows, against a trapezoid capacity that rises. The predicate is the
lobes pypulseq actually builds.

**The duration claim is scoped.** This construction gives the two winders the *same* duration and
the search returns the shortest such pair on the raster. Whether unequal durations could be
shorter is not established, no two-duration search was added, and no minimum-TE property is
claimed anywhere.

`prephaser_area_per_m` keeps its identity with `-area_to_echo_per_m` and is documented as the
**total** across the lobes. `prephaser_lobes` is new and exposes the events themselves, which is
what the record says a leaf should offer a composing module instead of a cached moment; it also
replaced two reaches into a private attribute, one of them from a notebook.

Validated in two layers. **Layer 1**, 45 tests, every moment integrated from the emitted events'
own knots, truncated at the echo and summed across every event on the axis — which is the part no
per-event check reaches, since the two winders and the readout's ramp-up only cancel together.
`m1` comes out at the arithmetic floor where an ordinary readout carries 0.03 to 0.37 s/m, over
full-echo, partial-Fourier, monopolar-train and bipolar-train geometries. Tolerances are
dimensioned per moment order. Seven of nine deliberate mutations fail the suite; the two survivors
are the gradient and slew checks inside the feasibility predicate, which are unreachable one
raster below the minimum because pypulseq refuses to build the lobe at all. **Layer 2**,
`examples/flowcomp_gre_2d/01_build.ipynb`, both readouts in a complete gradient echo: 0.680 ms
added to the echo time at the reference protocol.

**No Layer 3.** Given `m0 = 0`, the removal of the constant-velocity phase term follows from
`m1 = 0` by `phi = 2 pi (m0 x0 + m1 v)`, and a moving-spin simulation already exercised that
relation for `VelocityEncode`. The claim a simulation would add is image-level artefact
reduction, which the record explicitly does not make.

Readout axis only, and the first echo of a train only — later echoes accumulate their own first
moment from the lobes between them, and a test measures that rather than leaving it to the prose.
Still deferred: acceleration and higher orders, which need a fourth lobe; the slice-selection
case, which needs the excitation's moment from the RF isodelay point and is kernel work; and the
merged minimum-TE designs.

## Unreleased — velocity, written into the phase

`VelocityEncode` (`encoding/`): a bipolar pair on one logical axis. Two lobes of equal area and
opposite polarity, so the net area is zero and a stationary spin comes out with the phase it went
in with, while a moving one is somewhere else for the second lobe than it was for the first and
comes out ahead. Beside `PhaseEncode`, whose folder contract — gradients, no ADC, imposing a phase
you intend to sample — describes it exactly.

**`venc_m_s` in, and never a moment.** The public quantity is the velocity that accumulates `pi`
of phase *difference between the toggles*, and the quantity behind it is the **change** in first
moment, `delta_m1 = 1 / (2 * venc)`. Each toggle carries half of that. A caller who sets one
toggle's first moment to `delta_m1` is out by a factor of two and gets a plausible velocity map at
half the VENC they asked for, so the module takes `venc_m_s` and the mistake has nowhere to
happen.

No gamma appears in that relation. SeqCraft's gradients are in Hz/m, so a moving spin accumulates
`phi = 2 pi (m0 x0 + m1 v)` and `delta_m1 = 1 / (2 venc)` holds for any nucleus. The Handbook's
`pi / (gamma VENC)` is the same number in T/m units.

**`polarity` is the sign of the emitted first moment** — not of a velocity, and not of a
reconstructed phase. That is the candidate record's one open API question, answered: one module
holds the design and both toggles, and `build(polarity=...)` emits one, which is `PhaseEncode`'s
design-once/build-per-acquisition idiom. Defining the sign on the waveform is what makes it
checkable, because `sc.moments` reads it straight back off the events and nothing has to agree a
convention with a reconstruction to assert it.

The **continuous** hardware-limited design is closed form rather than a search. With lobe area
`A`, rise `r` and flat `f`, the two lobe centres are `2r + f` apart and `A (2r + f) = delta_m1 / 2`:
a cubic in the amplitude when the lobes are triangular, a quadratic in the flat time when the
amplitude pins at `max_grad`. Both times are then rounded **up** onto the gradient raster and the
amplitude is solved again against the realised ones, which keeps `delta_m1` exact and keeps the
gradient and slew inside their limits — rounding up can only lower the amplitude the target needs.

**That is the whole of the guarantee, and it is worth being exact about.** Rounding a continuous
optimum up is not searching the discrete lattice, so the emitted pair is legal and exact but
**not** the shortest legal pair the raster admits — on the nominal system at `venc = 1.0 m/s` it is
1.100 ms where a lattice search finds a legal 1.080 ms. No search was added, no `min_duration_s` is
published, and a test asserts the gap so the claim cannot become true by accident.

Once `m0 = 0` the first moment is independent of the time origin, so the pair can be placed
anywhere in a repetition and still deliver what it promised. That is what makes this a leaf: it
never has to be told where the echo is.

Validated in three layers. **Layer 1**, 64 tests, every moment integrated off the emitted events
by `sc.moments` rather than asked of the module: `m0 = 0`, the change in first moment across the
toggles, the half each one carries, `polarity` as the sign of `m1`, origin independence, the
hardware and raster behaviour, and the refusals. Seven deliberate mutations of the module — the
factor of two, the lobe order, the amplitude solve, the raster direction — were each checked to
fail the suite. **Layer 2**, `examples/pc_gre_2d/01_build.ipynb`, the pair inside a phase-contrast
gradient-echo repetition composed by hand from leaves, which is where its cost in echo time is
visible. **Layer 3**, `02_simulate_and_reconstruct.ipynb`, against a phantom whose spins move:

```text
|dphi| / pi against v / venc, 0 to 0.9 m/s     agrees to 1e-5 or better
a velocity map of three voxels at three speeds recovered to 2e-4 m/s
above VENC                                     wraps, as an angle must
```

Layer 3 also settled a sign. The measured phase difference is *negative* where the textbook
relation is positive, which could have been the module emitting the wrong sign or the simulator
reporting the opposite phase convention. A **static control** — a stationary spin at a known
offset under a deliberately un-nulled zeroth moment — shows the same `-1` ratio, so it is the
simulator's convention throughout and not something the encoding did. The notebook establishes
that before it divides by anything.

That check was worth running for a second reason: `sc.moments(..., 1)` is the instrument Layer 1
measures everything with, and its first-moment path had little independent exercise before this.
A moving-spin simulation is a different route to the same number, so the agreement checks the
analyser as much as the module.

Scope is the approved **appended bipolar** form and nothing else. The merged designs, where the
encoding lobe and an already-compensated imaging lobe become one shorter waveform, reshape
gradients this module does not own and are the next stage's problem; acceleration and higher-order
encoding, concomitant-field phase, and any claim about achievable minimum TE are all still
deferred. No compiler or `LogicBlock` change, and no moment-requirement IR or numerical optimiser
— the canonical case is algebra, which is the finding the candidate record rests on.

## Unreleased — T2 contrast for a sequence that has none of its own

`T2Prep` (`preparation/`): a T2-weighting preparation. Tip down, refocus while the magnetisation
decays at **T2 rather than T2\***, tip back up, spoil what is left. What the imaging train then
reads is longitudinal magnetisation already scaled by `exp(-prep_time_s / T2)`, and the readout
can be anything at all, because the contrast is no longer its job.

Fourteen hard pulses at **one B1 amplitude**: a 90 tip-down, four `90x 180y 90x` composites at
1/8, 3/8, 5/8 and 7/8 of the preparation in the MLEV-4 pattern `+ + − −`, and a single `-90`
tip-up. At a fixed amplitude a flip angle is a duration, which is what lets three pulses be one
composite — and it makes the peak-B1 contract one number for the whole block. The **base 90** is
what is rounded onto the RF raster, not each pulse separately, so the 1 : 2 : 3 : 4 ratio and the
single amplitude are exact for any `refocus_duration_s` rather than only for convenient ones;
`refocus_duration_s` reports the 180 that was realised and `requested_refocus_duration_s` keeps
the request.

**The tip-up is a realisation, not the contract.** What a T2 preparation *is* comes to four
things — tip into the transverse plane, refocus while it decays there, tip back, spoil — and how
the tip-up is built is not one of them. `tip_up='simple'`, a single `-90`, is the default and the
canonical form. `tip_up='composite_270_360'`, `270x` then `360(-x)`, is an established
realisation offered as a named alternative. It is deliberately not a boolean and deliberately not
named `robust` or `b1_robust`: it is a choice of realisation, and a future adiabatic family would
change what `prep_time_s` *means* rather than how one pulse is built, so it is not a third value
of this parameter either.

`min_prep_time_s` is the true floor, not a convenient one. It is solved from the real distances
between a composite's refocusing instant -- the centre of its 180 -- and each end of its group,
which are not equal: a group is a run of raster-ceiled slots and every pulse carries a dead time
in front and a ringdown behind, so that instant is the group's midpoint only by coincidence. Half
the group quoted a floor 0.3 ms long on the default scanner and 1.1 ms long at a 300 us dead time
with a 20 us ringdown, and a test at those lopsided times now holds the binding seam tight.

**`prep_time_s` is the transverse period**, and both of its ends are plain RF centres — nothing
weighted and nothing derived. It runs from the centre of the `+90` to the centre of the `-90`
under `'simple'`, and to the centre of the **270** under `'composite_270_360'`, which is the
convention the implementation that uses that realisation books its own preparation time by.
Implementations differ at the first end too: booking edge to edge instead differs by one
tip-pulse duration, a few per cent of the quantity being measured. Every RF group starts on the
gradient raster — the seam between two pulses is where a block boundary has to be cut — so the
placement is quantised first and `prep_time_s` is **measured off it**.

Validated in three layers. **Layer 1**, 60 tests, every claim read back out of the emitted events
or the compiled blocks rather than from an attribute, and **both realisations asserted** rather
than the default tested and the other trusted: the pulse counts, MLEV-4 phases, the net `-90` of
each tip-up summed off the emitted waveform, the gaps either side of every refocusing pulse, the
single B1 under awkward `refocus_duration_s` values, no gradient at all while the magnetisation is
transverse, the spoiler outside the interval, and the short-preparation and unknown-`tip_up`
refusals. **Layer 2**, `examples/t2prep_gre_2d/01_build.ipynb`, in front of a shipped `GRE2D`.
**Layer 3**, `02_simulate_and_reconstruct.ipynb`:

```text
a known T2 recovered from the slope of ln S against prep_time_s   within 1 %, 40 to 250 ms
T2* swept 80 -> 16 ms, a factor of five, moves the answer by        3 %
the residual bias is T1     +2 % at a 500 ms T1, 0.1 % at 20 s
```

The second line is the one the module exists for: the weighting follows the tissue's T2 and not
the magnet's contribution to T2*. The third is a property of the sequence rather than of the
implementation, and it is why this is a contrast preparation and not a T2 mapping method.

Layer 3 also sweeps B1 and off-resonance across both realisations, as a measurement rather than a
claim, and finds **no consistent advantage for the composite**: both become unreliable away from
nominal B1, their errors are not identical, and at B1 = 0.8 the composite is the worse of the two.
That is enough for the only conclusion drawn from it — the composite realisation is not a
robustness feature — and it establishes nothing about which part of the preparation dominates B1
sensitivity, which shipping the module does not need. That sweep needed a different instrument. MRzero's phase-distribution graph estimates which coherence pathways matter
from the *nominal* flip angles, and a 360 is the identity only at nominal B1 — so the sweep uses
`mr0.isochromat_sim`, and the notebook shows the two instruments disagreeing before it does.

**The line ordering is the caller's, and it matters.** What the preparation stores is longitudinal
magnetisation, which recovers at T1, so every millisecond before k = 0 is acquired erodes the
contrast. Linear ordering reaches k = 0 after 315 ms of a 618 ms train; centric reaches it after
6, and the example uses centric.

Shipped as **`mlev4-hard` only**. Adiabatic BIR-4/AHP tip pairs and the B1-robustness claim that
lives with them, a two-composite train, an arbitrary refocusing count, the B1-robust variant whose
rephasing lobe lands inside the *imaging* sequence, and diffusion preparation are all still
deferred, and none was added to make the first mode fit. **No quantitative B0 or B1 robustness
claim is made by either realisation**, and choosing the composite tip-up does not add one.

The plan records were reconciled before any code was written: one acceptance bullet still carried
the older edge-to-edge timing wording, which contradicts the semantic definition by one tip-pulse
duration. The durable lesson from the realisation/contract split is recorded in the fine-scan
playbook's *Realisation is not contract*, not as a rule local to this module, and the source
registry gained an `educational-synthesis` role for material that is a canonical teaching account
rather than evidence that can overrule a domain reference. Direction for the phase after this is
`tools/module_mining/plans/current/2026-09-23_post_v0_implementation_to_authoring_roadmap.md`.

## Unreleased — SeqCraft does not emit RF over `max_b1`

**`pp.Opts()` defaults `max_b1` to 851.52 Hz — 20 µT — so the limit was always live.** Nothing
enforced it except `Refocusing`, which meant every other RF path could emit a pulse the scanner's
transmit chain forbids. An ordinary **1 ms 90° sinc peaks at 130 %** of that default.

The contract now has two layers:

| | |
|---|---|
| **each RF module checks its own pulse** | at construction, where it still knows what to suggest changing. `ConfigurationError`, with the measured peak in both Hz and µT |
| **the compiler checks the emitted sequence** | `check_rf_amplitude`, beside `check_event_sizes` and called from the same place. `HardwareLimitError`, naming the worst block and its origin path |

The backstop is what makes it a contract rather than four habits: it catches a raw
`pp.make_sinc_pulse` added straight to a `LogicBlock`, and any future RF module that forgets its
own check.

The invariant is literal — `peak_b1_hz(rf) <= opts.max_b1` — with no sentinel, so **`max_b1 = +inf`
is how a caller designs without a transmit limit, and only `+inf` is.** Zero, a negative, `NaN` and
an **absent** limit are each reported as unusable rather than treated as unlimited, at **both**
layers. `NaN` matters most — every comparison against it is False, so a check written as
`worst > limit` would let any pulse through silently — and `None` is the subtle one, because
`pp.Opts(max_b1=None)` falls back to the 20 µT default rather than to no limit. An absent limit
needs a guard before the pulse is built, since pypulseq's shaped factories compare against it while
designing one.

Zero is not pypulseq's convention here either: unlike `adc_samples_limit`, whose `0` is documented
as "no limit", `make_sinc_pulse` compares `rf_amplitude > system.max_b1` unconditionally, so a zero
limit makes pypulseq warn at `inf %`.

**The measurement is shared and the remedy is not.** `_support.check_peak_b1` measures
`max(abs(rf.signal))` against the limit; each module supplies its own fixes, because the right
advice depends on the pulse family — and a quoted number is only called a *floor* where rebuilding
at it is proved by a test:

| family | how the peak responds | what the refusal says |
|---|---|---|
| sinc, or gauss pinned by a **time–bandwidth product** | envelope stretches, so exactly `1 / duration` | the duration that fits, as a **floor** |
| gauss pinned by an explicit **`bandwidth`** | held to that bandwidth: 1000.0 Hz at 1 ms *and* at 2 ms | a **starting point** |
| SLR | the filter is recomputed | the same number, as a **starting point** to re-check |
| `SaturationPrep` | `bandwidth_hz` is fixed, so TBW moves with duration | likewise a starting point |
| `hypsec` | set by the sweep — **`bandwidth` does nothing**, 563.7 Hz at 40, 10 and 2 kHz alike, and duration does nothing either | `beta`, `mu` or `adiabaticity` |
| `wurst` | `sqrt(bandwidth / duration)`, so duration helps only as `1/sqrt` | `bandwidth` or `adiabaticity` |

Generalising `Refocusing`'s `1 / duration` estimate to every family would have produced advice
that does not work. Lowering `adiabaticity` is not free either — it is what B1 robustness is
bought with — and the messages say so.

Peak amplitudes are reported in microtesla alongside hertz, converted with **`opts.gamma`**, so
the numbers stay right on a non-proton system.

Two paths were emitting physically impossible pulses **with no warning from anywhere**:

- `IRPrep(pulse='wurst')` at pypulseq's default 40 kHz sweep asks 187 % of 20 µT at 10 ms and
  419 % at 2 ms. `make_adiabatic_pulse` does not warn.
- `Excitation(pulse='slr', pulse_opts={'filter_type': 'min'})` returns a waveform peaking at
  **33 492 Hz ≈ 790 µT** — 3933 % of the limit — and delivering 84° rather than the requested 90°.
  `make_slr_pulse`'s min/max-phase branch does not warn either. Both are upstream characteristics;
  this release is the first thing that refuses them.

Every example notebook still compiles unchanged, so no working protocol was in the way. Eleven
existing tests used an over-limit pulse as a **stand-in** for something else — a 1 ms 90 to give
the compiler a block to schedule around, a minimum-phase SLR to give the rephaser an asymmetric
envelope — and now run with `max_b1 = inf` through a new `unbounded_b1` fixture that says why.
One, `IRPrep`'s WURST case, narrows the sweep instead, because that is the remedy the family has.

## Unreleased — the examples teach, the docstrings document, and the history lives elsewhere

Editorial only: no behaviour, no API shape, no test behaviour. 27 of 29 example notebooks touched,
26 of them markdown-only with code cells and outputs byte-identical; the twenty-seventh re-executed
to identical numbers.

**Sequence first, subtlety second.** The rule is written down in `examples/README.md` →
*Writing example notebooks*, and it is what most of this diff is. A reader who has seen only the
title and the first two paragraphs should know what sequence this is, what it is for, and what the
notebook is about — before meeting the clever part. Eleven notebooks opened on the clever part
instead, so their titles moved from conclusions to sequences:

| was | is |
|---|---|
| A crossing through k = 0 is not a spin echo | What the refocusing pulse changes in spiral imaging |
| A spiral pays for its long readout in off-resonance | Off-resonance in spiral imaging |
| Which instant is the echo aligned to | Spin-echo spiral imaging |
| A 3D gradient echo, out of one kernel and two `for` loops | 3D Cartesian gradient echo |
| Three planes, which is where a reversed `kz` stops being a number | Reconstructing a 3D gradient-echo volume |
| The samples are not on a grid, and that is the whole problem | Reconstructing a non-Cartesian acquisition |
| Does it produce the image we expect? | Simulating and reconstructing a 2D gradient echo |
| Does the magnetisation do what the timeline says? | Simulating MPRAGE: the inversion recovery, and what the ordering costs |
| Two contrasts, one ratio, and no bias field | Reconstructing the MP2RAGE UNI image |
| Saturating one spectral component, and throwing it away | Fat saturation |
| The whole of k-space in one shot | Single-shot echo-planar imaging |

None of the catchy sentences was deleted. Each moved to where a reader can already decode it — *a
crossing through the origin is not, by itself, a spin echo* now closes the section that measures
0.7 of a cycle of unrefocused phase, instead of standing unexplained at the top.

Also removed from notebooks: openings that defended an abstraction (`gre_3d/01`,
`gre_radial_2d/01`, `se_spiral_2d/01` all began on why some class does not exist), "the check this
notebook exists for", "a prediction that turned out to be wrong", reference-repository
comparisons, and `What this notebook established`. Kept: the physics, and the intentionally wrong
cases that teach MRI — three echo placements that all compile, an uncentred readout lobe, and the
demonstration that a spin echo does **not** fix EPI distortion.

**Public docstrings.** `SaturationPrep` was a design record with physics in it and is now API
documentation with the same physics. `CartesianLine`, `SpiralReadout`, `TSEShot`, `Refocusing`,
`GRE3DTR`, `IRPrep` and `Excitation` lost their abstraction-selection arguments and kept their
equations. `sc.b_value`'s docstring now states the contract — integration over the whole tree, the
excitation as phase origin, conjugation at each refocusing centre, and what `end_s` chooses — and
its validation history moved to the module-mining findings.

**Corrections found by rewriting**, each verified against the implementation:

| | |
|---|---|
| `gre_spiral_2d/02` said time segmentation treats Δf as constant within each interval | it does not: the code applies phase screens at representative times and **interpolates linearly between them**, so the approximation is set by temporal resolution. The "four segments aliases" claim became the measured 0.72 of a cycle between adjacent screens |
| `analysis.py` said "four answers" | there are five; `b_value` is now in the exactness table, as numerical rather than exact |
| `b_value` implied universal correctness | it models one excitation and zero or more refocusing pulses. Stimulated echoes and general coherence pathways are now stated as out of scope |
| `b_value`'s `end_s` said integrating past the echo "reports a number no experiment measures" | later samples do accumulate further weighting; `end_s` chooses *which* instant's `b` |
| `DiffusionSEPrep` called `axis=('x', 'y')` an oblique direction | it applies **equal** per-axis weighting, so that is the equal-component diagonal. Corrected, and the test renamed |
| `DiffusionSEPrep` put tensor encoding above this layer | tensor-valued encoding changes the waveform. Now: not implemented by this class |
| `_SLEW_FRACTION` implied 0.7 was source-backed | the source motivates staying below maximum slew and prescribes no fraction |
| `dwi_se_epi_2d/01` imported `seqcraft.modules.kernel.diffusion_se` | an internal path in a teaching notebook; the equation is now written out in the notebook |
| eight constructor parameters of `DiffusionSEPrep` were undocumented | documented |

A final pass took the remaining openings that still leaned on a preceding example:
`gre_2d/01` → **2D spoiled gradient echo**, `se_2d/01` → **2D spin-echo imaging**, `mprage_2d/01`
→ **MPRAGE — Magnetization-Prepared Rapid Gradient Echo**, `mp2rage_2d/01` → **MP2RAGE —
Magnetization-Prepared 2 Rapid Acquisition Gradient Echoes**, `se_epi_2d/01` → **Spin-echo
echo-planar imaging**. Each now states its own sequence structure before any cross-reference, and
the guideline gained the corollary: *do not assume the reader has completed a preceding example.*
Special-purpose documents — API guides, preserved references, migration guides — may still lead
with their document purpose, which is why `fse_2d/03_module_api` is unchanged.

Two stale factual claims fixed: `01_getting_started` said "seqcraft ships no concrete modules",
which has not been true for some time, and linked to an `examples/_parked/` that no longer exists;
`examples/README.md` said every module was extracted from the `gre_2d/` pair, which is no longer
true of Diffusion, Spiral or Saturation.

`AGENTS.md` added as a four-line pointer, so `examples/README.md` is the single source of the
notebook rules and `CLAUDE.md`, `AGENTS.md` and the module-mining skill only point at it.

## Unreleased — a b-value in, and the analyser that had to be corrected to prove it

`DiffusionSEPrep` (`kernel/`): a diffusion-weighted spin echo. `b_s_per_mm2` in; the lobe width,
the gradient amplitude and the echo time out.

**The first kernel whose reason for existing is a pulse it does not own.** The two encoding lobes
have the *same* polarity and sit either side of a 180 that conjugates what came before them, so
neither lobe can be designed without knowing where that pulse is and how long it lasts — and the
pulse belongs to `Refocusing`. One owner holds all three.

**It needed no new interface to do that.** Every fact the design requires was already a
construction-time property of the two RF leaves — `time_to_center`, `time_to_rephaser`,
`rephaser_duration_s`, `time_to_crusher`, `crush_duration_s`. The kernel constructs its leaves,
reads them, completes the coupled design, and only then calls `build()`. Nothing inspects emitted
events and nothing rewrites another module's output. `LogicBlock`, the compiler and
`GradientEvent` are untouched.

**The design is closed form end to end.** Substituting `Δ = C + δ + ε` into the Handbook's
trapezoid expression makes it a cubic in the lobe width — the source's own Example 9.2 — and both
timing windows grow as TE/2, so the shortest echo time is a maximum of two linear expressions. No
optimizer, no search. `salvage/bvalue.py` says "a closed form does not exist because `Delta`
depends on `delta`" and steps up the raster; it stopped one substitution short. The two agree to
1e-9, which is a test.

**The ramp terms are in it.** The correction is negative, so the rectangular formula over-states
`b`. The MIT reference implementation's own helper is documented "for trapezoid gradients: TODO".

New: `sc.b_value(tree, opts, end_s=...)` — the b-value per axis from the emitted gradients, which
shares no code with the module it checks. That mattered twice. It caught a lobe-centre separation
one ramp short, where every lobe was legal and the amplitude solve had quietly compensated by
delivering 4 % too much `b`.

**And then a Bloch simulation caught two defects in `b_value` itself.** It reported nine times the
weighting the magnetisation experienced on an unweighted acquisition, and a *different* number for
two files whose encoding is identical and whose only difference is 26 ms of dead time — which is
the diagnostic, because `k` is zero while nothing is playing.

| | |
|---|---|
| the phase origin was the start of the tree | a slice-select rephaser undoes the lobe area *after* the RF centre, not half the lobe, so half a lobe was left as a constant `k` |
| refocusing pulses were located without `event.delay` | `pp.calc_rf_center` measures from the start of the waveform; the conjugation was applied 700 µs early |

Neither was visible on the diffusion axis, where `k` is zero before the encoding and flat across
the refocusing block. Both are fixed and pinned by regression tests. **An independent validator is
only independent of what it was pointed at.**

New examples in [`examples/dwi_se_epi_2d/`](examples/dwi_se_epi_2d/): `01` builds a complete
single-shot diffusion-weighted spin-echo EPI at b = 0, 500 and 1000 s/mm², **all at one echo
time**, because a diffusion coefficient comes from a ratio and a mismatch adds
`(TE_b − TE_0)/(T2·b)` to every reported ADC. `02` plays them through MRzeroCore and recovers a
known diffusion coefficient to 0.06 %.

Evidence and the architecture argument:
[`tools/module_mining/plans/diffusion/`](tools/module_mining/plans/diffusion/).

## Unreleased — one reversible spiral arm, and the four ways to traverse it

`SpiralReadout` (`readout/`): the Nyquist path, its slew- and amplitude-limited traversal, the
raster waveform, the segmented ADC the interpreter's sample limit forces, and the trajectory
measured off what is emitted.

**One module, not four.** Everything rests on one decision — an arm begins and ends at rest, so it
can be played backwards and two arms join with no connector. `out`, `in`, `in-out` and `out-in`
are then a choice of how many arms and in which order. The independent reference,
`pulseq/pulseq`'s `writeSpiral.m`, ends its spiral-out at **full gradient** and has no spiral-in:
those are the same fact.

That policy costs **0.20–0.78 %** of readout duration, measured over seven protocols. It is a
braking distance, so it is fixed in absolute terms and shrinks as readouts lengthen; k-space extent
is identical and peak gradient is unchanged or lower.

**`out-in` crosses the origin twice.** So `origin_crossing_samples` is a sequence from the first
line of code — a scalar would work for three variants and change meaning on the fourth. An origin
crossing is a fact about the trajectory; whether the physical echo lands there belongs to the
kernel above.

**The reported trajectory is integrated from the emitted knots**, not the design pass. It agrees
with `sc.kspace`, which shares no code with it, to **~1e-12 1/m**.

Three things the implementation found, each a silent failure:

| | |
|---|---|
| the rewinder sized against the *design* end-point | left **0.015 `dk`** behind; sized against the *emitted* end it leaves zero |
| limits measured on raster-spaced differences | the emitted first and last intervals are **half** a raster, so the real slew is twice what that reports |
| the acquisition outlasting its gradient | the compiler holds the last block open and pads the waveform with area nobody designed |

A fourth, found by minimizing to a single failing protocol: each ADC event is rounded up to the
gradient raster **as well as** spending a lead delay and a trailing dead time, so sizing the sample
budget by subtracting only the dead times under-estimates the overhead by one raster per segment.
The acquisition then outlasts its gradient by 10 µs, the compiler holds the last block open, and
the waveform is padded with 0.045 1/m nobody designed. The budget is now sized against the
placement rather than an estimate of it. 96 protocols — four matrices by three shot counts by four
variants — all compile.

A fifth, in the same function, found by building the reconstruction against it: that budget only
ever searched **downward** from its conservative estimate, so it discarded up to five samples per
segment. Harmless for `out`, whose origin is at the first sample — and not harmless for `in`,
which brakes to rest at k = 0 and therefore keeps its last `dk` in its last 20 µs. The discarded
tail *was* the echo: at `matrix=128` the final sample landed at 1.004 half-steps from the origin
and `origin_crossing_samples` came back empty for a variant whose contract promises one crossing.
The budget now grows as well as shrinks. Six of 24 swept protocols were affected;
`test_origin_crossing_count` is parametrized over six protocols rather than asserted at one,
because the assertion was already there and the single fixture sat on the lucky side of a knife
edge.

A sixth, found by composing a spin echo against the module rather than by testing the module:
`time_to_echo` and `sample_times_s` were documented as measured from the readout block's start and
were measured from the **arm's**. Those differ by the prephaser, so `'in'` and `'in-out'` reported
every sample time and every origin crossing 300 µs early — in a sequence that compiles, whose
trajectory is correct to 1e-12, and whose only symptom is that a kernel composing TE from
`time_to_echo` puts the echo 300 µs off. `'out'` has no prephaser, which is why the one shipped
consumer never saw it, and `RadialReadout.time_to_center` had been block-relative all along.
Reporting is now on the block's clock and `prephaser_duration_s` is the bridge to the arm's; no
emitted waveform changes.

`examples/gre_spiral_2d/01_build.ipynb` is the complete acquisition: Excitation, one arm, eight
interleaves and spoiling, assembled in the notebook rather than by a class. It exists for the
join — **TE is owned there, not by the readout** — and checks it against the compiled sequence,
where it comes out identical across all eight interleaves to 0.000 ns.

`echoes > 1` refuses, naming the contract; `out-in` already exercises the plural machinery at
`echoes=1`. No 3D, no anisotropic FOV, no density presets, no `GRESpiral2D`/`SESpiral2D`, and the
three PR23 helpers stay private: one consumer each.

## Unreleased — a pulse whose purpose is to destroy what it makes

`SaturationPrep` (`preparation/`): a spectrally selective pulse, offset from water by a chemical
shift, followed by a spoiler and no rephasing. The sibling `IRPrep` was always going to have —
`modules/__init__.py` defines the folder as `rf.use in {inversion, saturation, preparation}` and
only the first seat was taken.

**It is not `Excitation` with a different label.** The two operations are opposites: an excitation
makes transverse magnetisation to be read and its selective form owns rephasing so the signal
survives to the echo; a saturation makes it to be thrown away. `Excitation(use='saturation')` was
considered and rejected — a semantic-purpose flag turns a module with one physical contract into a
generic RF pulse whose label carries the meaning. Reuse of waveform-design machinery is not reuse
of the public physical abstraction.

**The chemical-shift sign is the failure that looks fine.** `shift_ppm` is signed and relative to
water; fat is below it. Get the sign wrong and you get legal Pulseq, legal timing, legal gradients,
a correct-looking waveform — and water saturated instead of fat, with nothing downstream noticing.
The two published references reach the same number by different routes, one carrying the sign in
the ppm constant and the other applying it at the point of use, so an implementation that mixed
the conventions would be exactly this wrong. The module converts once, reports `offset_hz`, and
the tests trace `shift_ppm -> offset_hz -> emitted rf.freq_offset` as far as the **compiled
sequence**.

`freq_ppm` was considered for the emitted event and rejected: letting the interpreter resolve ppm
against its own B0 is more portable and moves the one number we are most likely to get wrong off
the file we can inspect.

**`flip_deg`, `duration_s` and `bandwidth_hz` have no defaults.** Two independent references
disagree on all three — 110° against 90°, 8 ms against 12, a bandwidth derived from the offset
against a fixed 200 Hz — and both work. That is what a protocol parameter looks like; picking one
lab's number would encode a protocol as physics. What the module does own is the relationship: it
refuses a pulse whose excited band reaches water, and reports `band_edge_hz` so the margin is
visible. At 1.5 T that margin is 8 Hz, against 816 Hz at 7 T.

No selection gradient, therefore no rephaser and no `thickness_mm`. Spatial saturation slabs and
CEST trains are deliberately not folded in: what they share with this is "prepare, then spoil",
which is a waveform silhouette rather than a physical solve.

`examples/fat_sat/01_build.ipynb` shows it in composition with `GRE2DTR`, against the same GRE
without it, and reads `use` and `freq_offset` off both compiled files.

## Unreleased — the m1 tolerance gets a floor, because the relative term sat in the noise

`verify_against_tree`'s first-moment check refused legal sequences, and which ones was not
predictable from the physics. An arbitrary gradient the compiler has to split and rebuild does not
round-trip exactly, and **the residual is proportional to peak `|g|` and to nothing else** — not
to knot count, not to duration, not to `max_grad`. The relative term scales as traversed area
times the horizon, which grows faster, so the two cross.

Measured on raw `make_arbitrary_grad` waveforms with no module involved: **six of ten legal
sign-changing gradients were refused**. `n=2000` at 90 % passed with a residual of 7e-15;
`n=4000` at 90 % failed with 1.2e-7. Same physics, opposite verdicts.

The tolerance is now `max(1e-9 · traversed · horizon, 1e-11 · peak)`. Both bounds are properties
of the representation and of the failure being caught, not of any one module:

| | |
|---|---|
| measured residual | `8.4e-14 · peak`, constant to 5 % across 2000–8000 knots, 20–90 % amplitude, two `max_grad` |
| the floor | clears it by **~120×** |
| one raster of displacement — what m1 exists to catch | `traversed · raster`, **500–50 000×** above the floor |

`tests/compiler/test_m1_floor.py` asserts the measurement the tolerance is derived from, both
margins, and that a trapezoid moved by one raster is still refused.

**m0 is deliberately unchanged.** Its `magnitude` accumulates `|∫g|` under a comment describing
`∫|g|`, and those differ by six orders of magnitude for one sign-changing event — but no
configuration made the check actually fail, so the mismatch is recorded and not patched. A
correctness alignment needs its own reproducer.

## Unreleased — a 3D repetition, and the axis where two moments become one gradient

`GRE3DTR` (`kernel/`), a **sibling** of `GRE2DTR` rather than a wrapper around it: a kernel that
contained one would have to reach back inside decisions that kernel has already made — its winder,
its TE, its rewinding — to change the z axis. Nor is a slice-selective 2D acquisition the centre
partition of a 3D slab. Both compose the same leaves, independently.

On x and y a 3D repetition is a 2D one. On z it is not, and `slab_thickness_mm` picks the case —
a physical quantity rather than a `slab_selective=True` flag, mirroring `Excitation`:

| | |
|---|---|
| `None` | non-selective, as every official Pulseq 3D reference is. The z axis carries a partition encode and nothing else |
| a thickness | slab-selective. The rephasing the slab implies and the partition encoding are two moments on **one axis in one window**, solved as `A_z(p) = A_slab + A_partition(p)` and realised as a single gradient |

Played in sequence those two cost two windows of echo time; `fmrifrey/lps` does exactly that and
its own `te_min` pays for both. Played as one they cost one. The non-selective case is the same
expression with `A_slab = 0`, not a second code path.

**Both terms are signed, so the limiting partition is a result and never an index.** With
`A_slab = -120` and partitions `-200 | 0 | +200` the low edge limits at `-320`; reverse the slab
gradient and the high edge limits at `+320`. Every partition is enumerated *after* the signed
combination — `_limiting_index` is a pure function so that a test can hand it a slab term of
either sign — and nothing infers sign from array order, which the Pulseq 3D demo would punish
anyway: it reverses its z lattice relative to y for a vendor reconstruction.

**One winder duration serves every partition.** Letting each take its own shortest would make TE a
function of `kz` — a contrast gradient across the volume that no k-space check shows and no
reconstruction expects.

**When the worst partition needs longer than x and y do, the window lengthens.** That is design,
not legalization: the kernel holds `opts`, so how long a moment needs is its question, while the
compiler may never rescue an infeasible module by stretching a gradient or moving a readout. With
`te_s=None` the result is the shortest legal design; an explicit request below it raises, naming
the partition responsible and its combined moment.

### The excitation mode chooses the pulse

`slab_thickness_mm=None` now defaults to a **short hard block pulse**, 0.2 ms, rather than
inheriting `Excitation`'s 3 ms shaped sinc. A non-selective excitation with a shaped pulse is the
worst of both — it spends a soft pulse's duration and selects nothing — and it was putting 3 ms
into every echo time for it. Measured at the reference geometry, the minimum TE drops from 4.377
to **2.977 ms**. Every official Pulseq 3D reference uses a block pulse for the same reason.

A slab still gets a shaped sinc, and `rf_pulse=` overrides either — a shaped but spatially
non-selective excitation is a real thing to want. Asking for a time-bandwidth product without a
shape is refused here, naming the fix, because this module chose the block pulse and the caller
should not have to discover that from a message about a pulse they never asked for.

### `Excitation` states its rephasing requirement, and can decline to realise it

`rephaser_area_per_m` is the signed moment that compensates a selective pulse, integrated from the
RF's effective centre — the same thing `CartesianLine.area_to_echo_per_m` is for a readout, and it
agrees with the rephaser the factory designs to 6e-14. `build(rephase=False)` omits that gradient
so a caller can take the realisation over; it does **not** mean the pulse needs no rephasing, and
the docstring says so.

Together they are what stops the slab term being applied twice. **Default behaviour is unchanged**
— `rephase=True` — and a test pins the default events by content hash.

Validated against `pulseq/pulseq`'s shipped `matlab/demoSeq/gre3d.seq`, which needs no MATLAB run:
Δk agrees on all three axes at 5.000 / 5.000 / 6.250 1/m, and the encoded index of every probed
line and partition lands where asked. Compared on lattice and semantic centre rather than sample
position, because the reference prephases with `-gx.area/2` and has no sample at the centre of
k-space where `CartesianLine` puts one.

Spoiling defaults to `('x',)` where `GRE2DTR` uses `('x', 'z')`: the z axis already carries the
partition rewinder in the tail, and two gradients each designed at the full slew limit do not sum
to a legal one. `writeGradientEcho3D.m` spoils on x for the same reason.

## Unreleased — an ADC sample count the receiver can actually digitise

`CartesianLine` now refuses a sampling geometry whose ADC sample count is not a multiple of
`opts.adc_samples_divisor`, at construction, naming the two arguments that produced it and a
combination that works.

**The rule is on the count, not on the parity of the matrix.** `matrix=65` at
`partial_fourier=1.0` gives 65 samples; an even matrix gives an illegal count just as easily, and
0.6 of 128 is 77 — which is why `examples/fse_2d` runs HASTE at 0.625, a fact the notebook stated
in prose and nothing enforced. Checked here because this is the layer that *computes* the count.
The compiler still refuses an illegal sequence, but it does so after the caller has built one and
names a block index rather than `matrix` and `partial_fourier`.

**Nothing is rounded.** Moving 65 samples to 64 would quietly change the partial Fourier fraction,
the k-space extent, the semantic centre sample and the sampling density.

`RadialReadout` inherits the refusal by composing the line rather than repeating the rule.

## Unreleased — a radial spoke, built out of the line it already is

`RadialReadout` (`readout/`): a prephaser, a readout gradient and an ADC, oriented in the x-y
plane. A leaf, because it determines every event it emits from its own parameters — given the
field of view, the matrix, the dwell and the angle, nothing above it is consulted.

**It is built out of `CartesianLine`.** A 2D spoke *is* a Cartesian line pointing somewhere other
than along an axis, and the prephaser, the dwell arithmetic, the ADC and — the hard one — which
sample carries `k = 0` are the same problem, already solved and already tested. Designing them
again would be a second source of truth for one piece of arithmetic. What the new module adds is
the rotation, the trajectory geometry a caller needs, and a contract stated as a spoke.

**`partial_fourier` spans full spoke to centre-out**, reusing the existing parameter rather than
adding `readout_asymmetry` beside it. At `matrix=64`: `1.0` gives 64 samples with the centre at
32, `0.75` gives 48 with the centre at 16, `0.5` gives 32 with the centre at 0 — a centre-out
spoke — and `Δk = 1/FOV` never moves. The range is closed at both ends, unlike `require_range`'s
half-open convention, because 0.5 is not degenerate here: it is a sequence family.

**The centre sample is not the ADC midpoint.** A full spoke with an even matrix has its centre one
sample past the middle, `-32Δk … +31Δk`, which is `PhaseEncode`'s `matrix // 2` convention and the
official reference's; a centre-out spoke has it at sample zero. `center_sample`,
`time_to_center()`, `dk_per_m`, `k_first_per_m`, `k_last_per_m` and `k_max_per_m` are the module's
answers. They exist because of what the alternative looks like: OpenMRF's radial readout compiles
a probe sequence and reads its trajectory back to discover where its own centre sample landed.

**The block `build` returns is already oriented**, so the module's semantic properties and its
emitted events have one owner rather than two. One canonical spoke is designed along x and fresh
rotated copies are derived per call — with the stored `area` scaled too, which both references
that extend a readout's flat time warn in a comment is otherwise left wrong. A spoke along an axis
emits one gradient, not one plus a 1e-11 Hz/m ghost.

No kernel and no trajectory layer. A radial GRE is currently `Excitation` → `RadialReadout` →
spoiler → delay, which is composition rather than a shared physical solve that no leaf can do
alone — the thing `GRE2DTR` and `TSEShot` exist for. Angle schedules, spoke counts and ordering
stay with the caller, as OpenMRF's own `phi_mode` does by building its angle list before any
waveform exists.

Extracted the same way as the last one: the references were measured first, and the
rotation-equivariance test — that the spoke at φ is the spoke at 0 rotated — was written against
the official PyPulseq reference and passed there **before this module existed**, so it cannot be
encoding this module's behaviour. It now passes against the package at ~1e-14 /m, across eight
sweep cases spanning two matrices, two fields of view, three dwell times and the full
partial-Fourier range.

[`examples/gre_radial_2d/01_build.ipynb`](examples/gre_radial_2d/01_build.ipynb) shows it, and shows the
paragraph above rather than asserting it: the repetition is assembled in the notebook out of
`Excitation`, the spoke, a spoiler and TR fill, and equal-increment and golden-angle schedules are
two list comprehensions over one readout instance. It is **build and trajectory visualisation
only**. A radial image needs a non-Cartesian reconstruction and the example suite has none to
reuse, so writing one to complete a `01`/`02` pair would be new reconstruction infrastructure
justified by a directory listing; the geometric claims are measured on the compiled trajectory
instead, and the agreement with the official reference above is the stronger evidence anyway.

## Unreleased — a turbo spin echo, split where the information is

`TSEShot` (`kernel/`) and `FSE2D` (`imaging/`), extracted from
[`examples/fse_2d/01_build.ipynb`](examples/fse_2d/01_build.ipynb) with **every event held fixed**:
`tests/modules/test_fse_notebook_matches_the_package.py` compares the two by absolute time and
content hash at turbo 1, 8, 16 and at a 72-echo single-shot HASTE, and they are identical. The
notebook still writes its own class, for the reason the GRE notebook does.

`echoes=1` is a conventional spin echo and one long train with `partial_fourier` below 1 is HASTE.
Neither gets a class, and the fine scan behind this change found no physics in either that the
composition does not already express.

**The split is by which layer has enough information to determine the value**, which is a sharper
rule than "avoid duplication" and occasionally disagrees with it:

| | owns |
|---|---|
| `CartesianLine` | how asymmetric it is about its own echo — new: `area_after_echo_per_m` and `echo_moment_imbalance_per_m` beside the existing `area_to_echo_per_m` |
| `Refocusing` | its crusher pair, balanced about its own effective centre, on the selection axis. **Unchanged** |
| `TSEShot` | the crusher window three axes share, the echo spacing it implies, the readout-axis lobes, and the shift that keeps every echo at the midpoint |
| `FSE2D` | how many shots, which lines each acquires, how many dummies, and therefore the effective TE |

The readout-axis lobe pair straddles a refocusing pulse, so it looks like `Refocusing`'s work. It
is not: its size depends on the *readout's* geometry as well as the pulse's, and a refocusing pulse
that knew about readouts would be wrong for SE-EPI, for a diffusion spin echo and for
spectroscopy. The readout states the geometry, the kernel decides the compensation — and
`(gx.area - 2*area_to_echo)/2`, written by hand in two notebooks, is now a property of the readout
that wrote it.

The window solve is two small pure functions rather than a shared abstraction. `GRE2DTR` performs
the same coupling under the name "winder coupling", and two examples is not enough to define the
contract for a third.

**Three invariants a train has that a single refocused echo does not**, each a legal block
structure when violated and each measured in `tests/modules/test_tse_shot.py`: the refocusing
pulses are uniformly spaced, every echo sits at the midpoint between its two pulses — primary and
stimulated echoes coincide only there — and `ky` returns to zero at every pulse, which is why the
blip goes *after* the refocusing pulse and the rewinder *before the next one*. A single-echo spin
echo may legally put the blip first; a train that does makes every later pulse flip the encode.

The extraction was argued from measurements rather than from reading, against the official PyPulseq
`write_tse.py` and the Pulseq MATLAB `writeTSE.m`. At a protocol both can express, this package and
the PyPulseq demo are the same sequence — identical event inventory, echo times agreeing to 5e-14 s
— in 185 pulseq blocks against 350. Two differences are worth knowing about rather than
reconciling: the demo has **no DC sample**, its 64 samples straddling `k = 0` by half a dwell, and
its construction is exact only while `rf_dead_time == rf_ringdown_time`, which it sets them to be.
`Refocusing` has never assumed that, and the notebook this module came from runs at a 30 us
ringdown, where the demo cannot be built at all.

## Unreleased — the MATLAB frontend is withdrawn

Out of the package and parked whole in [`salvage/matlab-frontend/`](salvage/matlab-frontend/):
`matlab/`, `examples/matlab/`, `src/seqcraft/exchange/`, `src/seqcraft/cli.py`,
`src/seqcraft/__main__.py`, `tests/exchange/`, `tests/fixtures/lbtx/` and
`docs/matlab_interoperability.md`. Out of the build with them: the `jsonschema` dependency, the
`seqcraft` console script, and the wheel's schema force-include. Full argument in
[ADR-006](docs/adr/006-the-matlab-frontend-is-withdrawn.md); [ADR-005](docs/adr/005-logicblock-tree-exchange.md)
is kept and marked superseded, because what it records is what the argument rests on.

**`sc.compile` does not move.** The compiler never imported the exchange package, so the compile
path, the `Module` contract, the warning policy and the return type are untouched, and every
remaining test passes unedited. That is the one good consequence of the boundary having been thin.

Three reasons, and the first is the one that matters:

| | |
|---|---|
| **A module is a callable, and MATLAB cannot say that** | `readout(line=17)` is the interface; `__call__` is the single place the block is checked and named. MATLAB's only equivalent is overloading `subsref`, which breaks object arrays, so the port split one lifecycle into `build` and a protected `buildImplicit`. A second name for one idea, in the class whose entire purpose is to have one |
| **Nothing MATLAB ran in CI, and nothing could** | The MATLAB tests are `crossval` -- public runners carry no licence. The only guard in the normal tier was a **handwritten** fixture, so it asserted what we believed MATLAB emits rather than what it emits. `writeLBTX` hardcoded `"pulseq": "1.5"` as well, so the compatibility version was a claim and never a measurement |
| **The example did not demonstrate what it claimed** | `seqcraft_examples.GRE2DTR` was named after the Python kernel and documented as its counterpart, but designed every event inline instead of composing leaves, spoiled two cycles per voxel where Python defaults to four, had no `te_s`, and was compared only against its sibling script -- never against `sc.modules.GRE2DTR` |

A MATLAB user reaches Pulseq as they did before: official MATLAB Pulseq, or a `.seq` written by
Python and read with `mr.Sequence.read`. The `.seq` file is already the portable artefact between
the two languages, and it is the one both toolchains validate.

`docs/serialization.md` returns to saying nothing there is implemented, which is true again, and
the deferred list -- Python writer, MATLAB reader, round-trip, semantic hash, provenance,
diagnostics -- is void rather than pending. The implementation is parked rather than deleted because
the cost of a revival was never the code: it is how a module gets called in MATLAB without a second
lifecycle, how the wire format is tested on every pull request, and who maintains it.
[`salvage/matlab-frontend/README.md`](salvage/matlab-frontend/README.md) says where every path came
from and what to restore in the build.
## Unreleased — the front page shows a sequence diagram coloured by LogicBlock

[`docs/figures/gre-tr.svg`](docs/figures/gre-tr.svg): one GRE repetition drawn the way a pulse
sequence always is — RF, Gz, Gy, Gx, ADC, Label — and then coloured and boxed by **the module that
wrote each waveform** rather than by event type. That one substitution is the argument, and it
needs no sentences beside it:

| | |
|---|---|
| `Excitation` is one box across **two lanes** | a slice-selective pulse is one idea, and the RF and the slice gradient are not two components |
| three boxes **overlap** across one window on three axes | the slice rephaser, the phase-encode blip and the readout prephaser, written by three modules that know nothing of each other |
| the outer box is `GRE2DTR` | a module made of modules, emitting the `LIN` label itself |
| the strip along the bottom is the **only** part nobody wrote | `sc.compile` put those four boundaries there, straight through the boxes above |

[`tools/draw_sequence_figure.py`](tools/draw_sequence_figure.py) draws it from the events the
modules build and the blocks the compiler returns — the boxes are the tree's own nodes, so no
clustering heuristic decides what belongs together — and
[`tests/test_readme_figure.py`](tests/test_readme_figure.py) **regenerates it and compares bytes**.
A stale figure fails in the commit that made it stale, with a message naming the command that fixes
it. It is the argument `tools/check_api_reference.py` already makes for `docs/api_reference.md`,
applied to the page everyone reads first.

Every gradient lane is scaled to its own axis rather than to one shared maximum, because a spoiler
is sixteen times a phase-encode blip and one scale draws the blip as a flat line. Lanes nothing
writes to are not drawn at all, so a spoiled GRE gets no empty `Trigger` row — and a sequence that
fires one gets the lane back without anybody editing a list.

### Added — `rf_duration_s` on `GRE2DTR` and `GRE2D`

Forwarded unchanged to `Excitation`, defaulting to `None`, which means **the leaf's own default**
rather than a copy of it here. The kernel could set the flip angle and the slice but not the pulse
length, and that is usually the largest term in `min_te_s`: a 3 ms pulse puts 1.5 ms into TE before
a single gradient plays, which no readout can shorten. It adds no arithmetic — the winder and TE
measure the excitation they were handed rather than predicting it, which is what made forwarding
one argument the whole change.

## Unreleased — a Cartesian line, read more than once

A multi-echo gradient echo, as **three arguments on one module and three helper promotions.**
No new module, no new class — including in the example.

```python
mono = sc.modules.CartesianLine(opts=opts, fov_mm=220.0, matrix=128, bandwidth_hz_px=500.0,
                                echoes=8, polarity='monopolar')

assert len(mono.te_s) == 8                                         # k = 0, once per echo
assert abs(mono.flyback_area_per_m + float(mono.gx.area)) < 1e-9   # the WHOLE lobe
assert [mono.polarity_of(n) for n in range(4)] == [1, 1, 1, 1]     # 1, -1, 1, -1 under 'bipolar'

megre = sc.modules.GRE2D(opts=opts, fov_mm=220.0, matrix=(128, 128), thickness_mm=3.0,
                         flip_deg=15.0, bandwidth_hz_px=500.0, tr_s=40e-3,
                         echoes=8, polarity='monopolar')           # and that is the whole sequence
```

| | |
|---|---|
| `CartesianLine` | `echoes`, `polarity` and `echo_spacing_s`; and `te_s`, `min_echo_spacing_s`, `polarity_of`, `echo_sample`, `flyback_area_per_m`, `flyback_duration_s`. **`build` does not change** — not one argument, not one default — because the echo count is a property of the *design*, as `matrix` and `partial_fourier` are |
| `halve_onto`, `dwell_quantum`, `require_count` | promoted out of `readout/epi_2d.py` into `modules/_support.py`. Not on the usual different-folders argument — both callers are in `readout/` — but on the **direction**: `CartesianLine` is what `GRE2DTR`, `MPRAGE2D` and three notebooks already stand on, and `from .epi_2d import _halve` would point that dependency backwards through a private name |
| `GRE2DTR`, `GRE2D` | the same three arguments, forwarded. **No arithmetic in either**, which was the design's own falsification test and is what says a multi-echo GRE needs no class |
| `MEGRE2D` | **does not exist.** `examples/megre_2d/` is the first example directory that defines no class, and that is the result rather than an omission |
| everything else | unchanged. `design/`, `compiler/`, `analysis`, `EPI2D`, `Excitation`, `Refocusing`, `PhaseEncode`, `IRPrep` and `spoiler` are untouched, and `tests/modules/test_epi_2d.py` passes **unedited**, which is the acceptance test for the promotion |

### `echoes=1` is byte-identical, and it is pinned

Same dwell, same lobe, same prephaser, same block, same duration, same `time_to_echo()` — asserted
against `content_hash` digests captured from the module **before** the train existed, at three
bandwidths and three `partial_fourier` values. A change there is a claim about every existing GRE,
MPRAGE and SE example rather than about this feature, so it is the first test in the file. `ECO`
and `REV` are emitted only above one echo for the same reason: two extra events are not identical.

### Four numbers that are wrong quietly

**The reverse lobe's grid.** A trapezoid is symmetric about its own midpoint, so a reverse lobe
samples the forward lobe's grid only if the window is *exactly* centred — `T == 2*guard +
N*dwell`, with equality. The flat-top rounding that is right for a single line, and right for a
monopolar train because every lobe in it is that same lobe, puts the two polarities **0.9324 1/m =
0.2051 Δk** apart. Every k-space extent check passes on it; `max|k|` is identical to four decimals,
and only the *signed* difference separates them. It is a linear phase ramp between odd and even
echoes, which is exactly what a field map measures — so it is never an artefact, only a wrong
number. `'bipolar'` redesigns the lobe and the same measurement is **1.7 × 10⁻¹³ 1/m**.

**The fly-back's area.** Minus the lobe's *total* area, ramps included. The dangerous wrong answer
is `-area_to_echo_per_m`, because it is a number the module already computes and already exposes,
and it is **294.638695 against 585.664336 1/m** — half, near enough to look right in a plot, and
each echo would then start half a k-space further along than the last.

**The dwell.** The guard is the ADC's `delay`, and pypulseq's `_check_timing_block` divides those
by `rf_raster_time`, not by the ten-times-finer `adc_raster_time`. So `num_samples * dwell_s` must
be an even number of RF rasters, and under `'bipolar'` the dwell is snapped **up** onto that
quantum — 500 ns at 128 samples. 500 Hz/px asked for is **500.8 achieved monopolar and 488.3
bipolar**, which is a fact about the file and belongs in the protocol table.

**The echo times.** `k = 0` is a *sample*, not an instant: index `pre_echo_samples` forward and
`num_samples - 1 - pre_echo_samples` reverse, which are different distances into identical lobes.
So a bipolar train's ΔTE alternates by two dwells — measured **32.00 µs** at 500 Hz/px and
**1953.00 µs** at `partial_fourier=0.75` and 250 Hz/px, both matching `2·|N-1-2·p|·τ` to the
nanosecond. Nothing checks this. A two-point field map divided by `echo_spacing_s` instead of by
`te_s[1] - te_s[0]` is **0.75 % low in every voxel**; an eight-echo least-squares fit is 0.0358 %
low, because a symmetric alternation largely cancels. `te_s` is the only sanctioned source of echo
times, and `echo_spacing_s` is documented as the seam-to-seam period and nothing else.

### One prediction that turned out to be false

The design document said the barriers were load-bearing: that without them the compiler would take
the midpoint of the gap between consecutive ADCs, land inside a lobe, and split it into arbitrary
waveforms. **Measured, that is wrong for this design.** Two barriers, one and none compile to the
identical block count — 16 monopolar, 9 bipolar at eight echoes — with every gradient still a
`trap` and no `merge` warning, at 250, 500 and 1000 Hz/px. `find_boundaries` already offers every
*gradient edge* as an opportunistic candidate and accepts one that falls strictly inside no
gradient, and every seam here is exactly that: the lobe ends where the fly-back starts and the
fly-back ends where the next lobe starts, with nothing overlapping either join. That is the
difference from `EPI2D`, whose blip is **centred on** the seam and therefore covers it.

The barriers stay, because that rule is documented as *a preference, not a requirement* and one
barrier per gradient edge makes the block structure a property of the design rather than of a
preference. They cost nothing — the emitted file is byte-identical with and without — and the test
asserts the **measured** fact rather than the claim it replaced.

### `examples/megre_2d/`

Two notebooks, and **neither defines a class**. `01` builds the wrong lobe before the right one
and plots both, derives the fly-back out loud with all three wrong areas' trajectories drawn
beside it, and computes what assuming a uniform ΔTE costs **from `te_s` alone, with no
simulator**. It writes three files — monopolar, bipolar, and monopolar at `partial_fourier=0.75`
— each carrying `echoes`, `polarity`, the echo samples, the echo polarities and the `te_s` array
in its `[DEFINITIONS]`, so `02` reads the echo times *out of the file* rather than rebuilding the
module.

`02` puts all three through **one reconstruction** — ESPIRiT coil maps estimated once, POCS for
the partial echo, then SENSE — and fits both maps against oracles that are in the **phantom**
rather than in the fit: `subject05.npz` carries `T2dash_map`, so `1/(1/T2 + 1/T2′)` is the true
T2\* voxelwise, and `phantom.field_hz()` is the map the simulation ran against. All three readouts
land on a field-map slope of 1.001 and about half a hertz of scatter, and what is left of the T2\*
spread is partial volume at tissue boundaries, named rather than explained away.

One set of coil maps for every echo is not economy: it is what makes echo *e*'s phase comparable
with echo 0's, which is the entire ΔB0 measurement. And the oracle is resampled onto the
reconstruction's pixel grid before anything is subtracted — the phantom is 192 mm across and the
protocol images 220 mm, so an image pixel is not the phantom voxel with the same index; scored in
the phantom's own frame the same fits give a slope of 0.963 and seven times the scatter, which
reads as a calibration problem and is a geometry one.

### Out of scope, named

No ramp sampling — `area_until` would make it come out right and nothing here needs it, which is
also why there is no `k_read_per_m`: the ADC never opens on a ramp, so sample *i* of echo *e* is at
`(i - echo_sample(e)) * dk_per_m * polarity_of(e)`. No through-zero seam, worth one rise time per
echo — 40 µs of a 2130 µs period — and costing every lobe its `trap`; `EPI2D` made the same call.
No Dixon: what a water–fat separation needs from this design is the exact `te_s`, and that ships.
And **`polarity='bipolar'` with `prephase=False` raises**, because balancing reversed lobes across
a refocusing pulse is a different question and nobody has measured it.

## Unreleased — `EPI2D`, and the two rules that stop an EPI ghosting on its own

The whole of k-space in one shot, as **one new module and one promoted helper.**

```python
epi = sc.modules.EPI2D(opts=opts, fov_mm=220.0, matrix=(128, 128), dwell_s=2.5e-6)

assert 2 * epi.guard_s + epi.num_samples * epi.dwell_s == epi.echo_spacing_s   # exactly
assert epi.k_read_per_m[epi.echo_sample(0)] == 0.0        # k = 0 is a sample, not a moment
assert epi.polarity(1) == -1                              # so echo 1 is k_read_per_m[::-1]
```

| | |
|---|---|
| `EPI2D` | prephasers, alternating readout lobes, blips on the zero crossings, one ADC per echo, the barrier that keeps the lobes whole, and the `LIN`/`REV`/`NAV`/`SEG`/`REF`/`IMA` labels a reconstruction reads back. `readout/`, because it contains ADCs and no RF |
| `require_pair` | `_pair`, promoted out of `kernel/gre_2d_tr.py` into `modules/_support.py`. It now has two callers in **different folders**, which is the argument `shift_slice` already makes there. `GRE2DTR` imports it and its output is byte-identical |
| everything else | unchanged. `design/`, `compiler/`, `Excitation`, `Refocusing`, `PhaseEncode`, `CartesianLine`, `spoiler`, `IRPrep`, `GRE2DTR` and `GRE2D` are untouched |

### The two rules

**Rule 1 — the sampling window is exactly centred in its lobe**: `T = 2*guard + N*dwell`, with
equality. A trapezoid is symmetric about its own midpoint, so only a centred window makes sample
`i` of a forward lobe and sample `N-1-i` of the reverse lobe land on the *same* k. Write the lobe
the way `CartesianLine` correctly writes a **line** — flat top rounded up onto the gradient raster,
slack after the last sample — and the two polarities sample grids up to half a raster apart. That
is an N/2 ghost the sequence made itself, and no build-notebook check would catch it.

**Rule 2 — the blip is centred on the seam, and `sc.barrier()` is forced there.** The seam is where
the readout gradient passes through zero: the least the two axes ask of the amplifier together, and
the only placement needing half a blip of guard rather than a whole one — 700 µs of echo spacing
against 760. But a centred blip *straddles* the boundary the compiler has to cut between two ADCs,
and no instant in that gap cuts nothing. Left to choose, the compiler takes the midpoint and splits
**every readout lobe in the train** into arbitrary gradients: 127 of 129 `gx` events, on a sequence
that still compiles, still passes every k-space check and still simulates correctly. The only
report of it is a `merge` warning naming the readout axis. `EPI2D` is the first module in the
library to emit a barrier, and it is deliberate.

What the barrier costs is one `merge` warning on the **blip** axis per compile — the two halves of
each split blip summed back together, reported instead of hidden. `FSE2D` already produces 71 of
the same kind.

### Three things that are arguments, and one that is not

`oversampling` defaults to **2** and costs no echo spacing at all: twice the samples at half the
dwell is the same sampling duration, so the lobe, the guard, the peak gradient and the bandwidth
per reconstructed pixel do not move. What it buys is the sample *spacing* — at `oversampling=1` a
ramp-sampled flat top leaves gaps above one reconstruction `dk`, so the readout does not sample the
grid it is reconstructed onto and no interpolator can recover what was never measured.

`ramp_sampling` is not a bandwidth trick: it is 880 µs of echo spacing against 700 at the same
nominal rate, bought with a higher peak gradient and paid for with a non-uniform `k_read_per_m`.

`blip_lines` **prices an ordering before it is built** — 700 µs at one line and half as much again
at the 126 a centric single shot steps — which is information no table of integers carries on its
own.

And parallel imaging adds **nothing** to the readout. `blip_lines=R` is the acceleration and `lines`
is the table, both of which already exist for their own reasons. The single addition is
`reference=`, which emits the `REF`/`IMA` pair, because a reconstruction has to tell a calibration
shot from an imaging one and counting shots is how it comes to guess wrong.

What the calibration band is **not** is more EPI shots. A table of length one is a Cartesian
gradient echo on the EPI’s own readout lobe -- same gradient, same ADC, same ramp sampling, so one
regridding operator serves both, and no train, no blips, no alternating polarity and no
`ky`-versus-time ramp for GRAPPA to mistake for receive sensitivity. Because it *is* a Cartesian
gradient echo it is written as one: a small flip, the shortest echo time the lobe allows, and RF
spoiling. `gre_epi_2d/02` §7 prices it against a self-calibration taken from the fully sampled
data's own centre -- the ceiling no acquisition can reach -- and it costs 1.18x that at R = 3.
The term that does *not* appear is the one the argument predicted: acquiring the band at 60 Hz
costs 1.00x, because displacement follows the time from excitation and one lobe is 0.04 px where
the R = 3 train is 1.9.

### `EPI2D.build` takes `phase_deg`, which it should have had from the start

`CartesianLine.build` has always taken a receiver phase, and its docstring has always said it is
not optional in a spoiled sequence. `EPI2D.build` did not take one at all, so a spoiled EPI could
not be written: the transmitter would run a schedule the receiver knew nothing about, and that
schedule's quadratic phase would land in `ky` rather than in the carrier. One number reaches every
ADC of the train -- navigators included -- because one train is one excitation.

Found the way the other two below were: by building a **segmented** acquisition and simulating it.
Four interleaved shots put consecutive `ky` lines in different shots, so anything that differs
between shots is a period-4 modulation of `ky`, which is a replica of the object at `Ny/4` -- and
it reads as an under-sampling artefact. On a **shimmed** brain, with the ordering, the off-resonance
and the echo time all held fixed, `gre_epi_2d/02` §8 measures that replica across three files:

| file | TR | run-in | RF spoiling | replica at `Ny/4` |
|---|---|---|---|---|
| `gre_epi_4shot_naive` | 200 ms | 1 shot | none | 0.1499 |
| `gre_epi_4shot_tr1s` | 1 s | 1 shot | none | 0.0377 |
| `gre_epi_4shot` | 1 s | 2 shots | 117° | 0.0156 |

The mechanism is one number. A perfectly spoiled 90° excitation reaches its steady state in one TR
whatever the TR is, so the transient is entirely what the gradient spoiler fails to destroy, and
that is `exp(-TR/T2)`: **8.2%** at 200 ms against 4e-6 at 1 s. Hence the ordering of the table --
the TR is worth 4.0x and the run-in the rest, for 9.6x in total.

RF spoiling is therefore worth nothing measurable *at a 1 s TR*, and it is in the sequence because
it is correct and because at the 200 ms TR of a real segmented protocol it is not optional. Only
that third row is the module's; the other two are the notebook's, and this one could not have been
fixed in the notebook at all.

### Two refusals that were not refusals

Both found by using the module rather than by reading it, which is the argument for the notebooks
being part of the change rather than documentation of it.

**An over-limit readout amplitude escaped as pypulseq's exception.** `_check_hardware` ran *after*
the lobe was built, and `pp.make_trapezoid` refuses an over-limit amplitude first — with
`ValueError: Amplitude violation (133%)` and nothing else: no bandwidth, no field of view, no note
that ramp sampling puts the peak above the nominal `dk/dwell`, and no way to tell which of the
module's several trapezoids raised it. The check now runs on the amplitude *before* each build, so
the refusal is the module's own and carries its remedies. Surfaced by moving the examples to a
128 matrix, where the amplitude binds before the ramps do.

**A navigator echo could not be addressed at all.** `polarity`, `echo_sample` and `time_to_echo`
refused every negative index, but a reconstruction of a navigated shot needs exactly those: it
regrids the navigator readouts, which needs their polarity. `echo` is now an offset from the first
imaging echo, so `-navigator_echoes … -1` are the navigators and a reconstruction can write
`range(-navigator_echoes, len(lines))` to get the file's readouts in order. `time_to_echo` needed
more than a relaxed bound — a navigator is **not** `echo * echo_spacing_s` before imaging echo
zero, because the blip-axis prephaser plays between them, and that offset is the reason it is a
method and not a formula. Measured at 0 ps across the whole train.

### The four notebooks, and the measurement each one is for

`GREEPI2D` and `SEEPI2D` are written in `examples/gre_epi_2d/01_build.ipynb` and
`examples/se_epi_2d/01_build.ipynb` **and stay there** — one consumer each, the same rule that kept
`SE2D` and `FSE2D` out.

`gre_epi_2d/01` builds both rules the **wrong** way as well as the right way, because both failures
compile: dropping the barrier turns 63 of 65 `gx` events into arbitrary gradients, and letting the
blip end at the seam costs 570 µs of echo spacing against 510. It also prints `sc.pns` against
`sc.hardware.synthetic_hardware()` at 130, 150 and 180 T/m/s — the slew rate this family runs at is
what makes ramp sampling and partial echo build *together*, and it costs 19 % more predicted
stimulation.

`gre_epi_2d/02` measures the regridding operator against **exact** k-space and finds the answer is
`scipy.interpolate.CubicSpline` materialised as a matrix — it reproduces scipy to 3.1e-15 and
applies 18× faster because it does not depend on the data. Before that it measures the question
that comes first: at `oversampling=1` every interpolator fails, including on the band-limited object
the least-squares sinc is supposed to be optimal for. Then off-resonance against a `B0 = 0` control
(2 px at 60 Hz, predicted 1.96), the navigator correction and the Fourier **dual** of it that does
nothing at all, and interleaved against blocked at the same lines and the same shots — **0.059
against 0.152**, from the ordering alone.

`se_epi_2d/02` is the one worth reading for the result that is not what the name suggests: **a spin
echo does not fix EPI distortion.** Same shift, same prediction, same ghost as the gradient echo at
the same echo spacing, because the pulse refocuses at one instant and the distortion accumulates
between echoes. What it *does* fix is the envelope: 4.3× the signal on the lines that carry the
contrast, and half the excess point-spread width — measured against a **flat-envelope floor**,
because a finite window blurs on its own and without that row neither claim could be sized.

And since only a shorter train touches the distortion, §5 shortens one: the same 128 lines at the
same echo time, in one, two and four **interleaved shots**. The train falls 89.6 → 46.1 → 23.7 ms
and the residual against the shimmed image falls 0.127 → 0.073 → 0.047 with it. That is the train
an `R = 4` acceleration would buy, without a calibration band, a kernel or an unfolding residual --
affordable here precisely because a spin-echo TR is long for reasons that have nothing to do with
EPI. What segmentation *does* cost is that the shots have to agree, and the replica at `Ny/S` that
measures whether they do is 0.021 at two shots and 0.022 at four, against the run-in every file in
the notebook now carries.

**All four files are written at one echo time**, and both EPI notebooks now check that invariant and
print the residual spread -- 7.50 us across nine files in `gre_epi_2d`, 0.00 across four in
`se_epi_2d`. Every figure that had grown a per-family grey scale is back to one scale throughout.
At one T2 of echo time a contrast difference between two panels dwarfs anything the acquisition
does, and a figure that shares a ceiling across panels that do not share a contrast is one that
lies.

### The examples share one phantom, and one field map that is not a formula

`examples/phantom.py` is new and is not a module of the library: it is the BrainWeb slab every 2D
example simulates against, prepared once. Five notebooks had grown five copies of the same twenty
lines and they had drifted — two download URLs, two spellings of one slice range. The parameters
stay in the notebooks, because a slice count with a reason attached is a choice; the preparation
does not.

`field_hz` hands back **the phantom's own B0 map**, and the interesting part is what that map
turns out to be. `subject05.npz` carries PD, T1, T2, T2' and D and **no `B0_map`**; MRzero's loader
notices and calls its own `generate_B0_B1`, whose comment reads *"generate a somewhat plausible B0
and B1 map, visually fitted"*. Two Lorentzians, demeaned over the proton density. So it arrives
centred on zero -- nothing needs to shim it -- and it is analytic rather than measured, which the
EPI notebook says out loud before showing a distortion figure, because the file name invites the
opposite assumption.

It spans about -11 to +45 Hz, which at 128 lines and a 700 us echo spacing is -1 to +4 pixels of
displacement: real, mild, and the same field the other five notebooks run against, which is what
makes their numbers and the EPI ones comparable.

A first version of this module synthesised its own map instead -- a susceptibility volume with a
frontal sinus in it, through the dipole kernel, shimmed to second order, plus a skull shell and a
smoothing term to cover for BrainWeb being a brain in a vacuum. It produced a stronger and more
interesting picture and it was **~150 lines of physics nobody asked for**, sitting beside a map the
phantom already had and every other example already used. It is gone; `field_hz` is four lines.

**Orientation is the other thing in it, and it took three corrections to get right.** The phantom
is indexed `[x, y]` with **x left-right and y anterior-posterior, anterior at high y** — read off
the anatomy (the cerebellum and the splenium sit at low `y`; the inferior slices sit ten voxels
posterior of the brain's centroid) rather than assumed. A reconstruction is `image[ky, kx]`, so it
needs no transform and a *map* needs one, which is `same_frame`. And every example draws with
`origin='lower'`, without which the frontal lobe comes out at the bottom.

Getting the first of those backwards is not a cosmetic error: it placed the frontal sinus against
one ear, which showed up as a field map that looked asymmetric rather than as a mistake in a
coordinate. `tests/examples/test_phantom.py` now pins the anterior direction, the midline and the
transpose, because nothing else does — notebooks are not run by the fast test tier, and a phantom
that turns sideways makes every number in them wrong without making any of them fail.

## Unreleased — `Refocusing`, and the one rule that makes a spin echo one

Spin echo, turbo spin echo and HASTE, and the smallest thing they needed: **one new module and one
boolean.**

```python
refoc = sc.modules.Refocusing(opts=opts, thickness_mm=6.25, crush_voxel_mm=5.0)
line = sc.modules.CartesianLine(opts=opts, fov_mm=256.0, matrix=128,
                                bandwidth_hz_px=200.0, prephase=False)

assert refoc.area_to_center_per_m == refoc.area_from_center_per_m   # 600.000000 both
assert refoc.time_to_center() == refoc().duration / 2               # exactly, to 0 ns
```

| | |
|---|---|
| `Refocusing` | a 180 and its crusher pair, as **one continuous waveform** on the selection axis, symmetric about the RF's effective centre in time *and* in area. `rf/`'s second membership value — `use == 'refocusing'` — which until now had no member |
| `CartesianLine(prephase=False)` | drops the prephaser event and **nothing else**. In a CPMG train the readout's dephasing is half of a pair straddling a *refocusing pulse*, so there is no event of its own to cancel it |
| `CartesianLine.area_to_echo_per_m` | new, and always available: what the readout accumulates by the echo. Minus the prephaser when there is one, and the number a crusher pair is balanced around when there is not |
| `Refocusing.time_to_crusher()` | where the trailing crusher window begins, so the readout block starts **inside** the refocusing block |
| everything else | unchanged. `Excitation`, `PhaseEncode`, `spoiler`, `GRE2DTR`, `GRE2D`, `IRPrep`, `design/` and `compiler/` are untouched |

### The one rule

A refocusing pulse conjugates k, so between consecutive refocusing centres **each axis's gradient
area before the echo equals its area after it**. Everything in the design is one consequence of
that, and the two worth writing down are the two the references get wrong.

**"Equal-area crushers" means measured to the RF's *effective centre*, not the same trapezoid
twice.** The selection plateau's own halves are unequal whenever the transmit dead time and the
ringdown differ, or the pulse is asymmetric. On a 100 µs / 30 µs system that residual is 11.2 1/m —
and since $k_n = -k_{n-1} + \delta$ it **alternates sign echo to echo**: the odd/even modulation an
FSE is famous for, which reads as a hardware fault. `writeTSE.m` and `write_tse.py` avoid it only by
setting their dead time and ringdown to the same 100 µs.

**Two fixes, because they fix different halves.** The plateau is *symmetrised* about the effective
centre (`time_to_center() == duration / 2`, which is what puts the echo at the midpoint between two
pulses and what the stimulated-echo pathways depend on), and the crusher areas are then *solved by
integration* rather than by formula. With a symmetrised plateau the solve returns two equal
amplitudes — which is the point: it is then the check that says the symmetrisation worked, on a
pulse shape nobody has tried.

### Two sequences that deliberately did not become modules

`SE2D` and `FSE2D` are written in `examples/se_2d/01_build.ipynb` and
`examples/fse_2d/01_build.ipynb` **and stay there** — one consumer each, the same rule that kept
`MPRAGE2D` and `MP2RAGE2D` out. `Refocusing` ships because two notebook pairs use it, and it is the
obvious base for SE-EPI, diffusion and T2 preparation.
`tests/modules/test_se_notebooks.py` runs both build notebooks, asserts k at every echo for
`echoes` in `(1, 8, 16, 72)`, and pins `FSE2D(echoes=1)` against what `se_2d/01` writes — event for
event, which is what makes two example directories safe.

### Five things the extraction found

**The excitation needs no change at all.** A first draft added `Excitation.build(rephase=False)` so
the slice rephasing could be folded into the first crusher "as the references do". That was a
consequence of measuring the crusher area from the *block edge* rather than from the *RF centre*:
once `area_to_center_per_m` includes the selection plateau's own half — which it must, because that
is where the conjugation happens — the first-interval condition on `z` becomes *net zero from the
excitation's centre*, which is the ordinary rephaser.

**The readout block starts at the refocusing block's trailing crusher, not after it.** The slice
crusher, the phase-encode blip and the readout's own lobe all live in **one window on three axes** —
the same coupling `GRE2DTR` resolves for its winder, moved one sequence over. Laying the two blocks
end to end costs a measured **1.4 ms per echo, 13 %** of the minimum echo spacing, and 30 extra
pulseq blocks per eight echoes. `time_to_crusher()` exists for that and nothing else.

**The crusher belongs on the slice axis, and only there.** The two readout-axis lobes that bracket
a readout are a *balance*, not a crusher: the conjugation cancels whatever they share, so their
common area has no k-space consequence at all and is set to **zero** — an imperfect 180's FID is
already dephased by three cycles across the slice. What is *not* optional is their **difference**,
which is one half-dwell's worth of readout area (0.87 of one `dk` on the reference protocol), split
`+skew` and `-skew`. `writeTSE.m` and `write_tse.py` do add a readout-axis crusher (`rSpoilFac`) and
also drop the half-dwell term, which makes their two lobes identical at the price of a half-pixel
image shift. `examples/se_2d/01` derives the difference from the readout's own knots, so partial
Fourier and ramp sampling produce it too.

**There are three minima on the echo spacing and one of them is not about timing.** The train, the
first interval, and **`z` clearance under the ADC**: the next refocusing block's leading crusher may
not start before the readout's trailing one, or the slice gradient is on while the ADC is open. That
is signal loss rather than an illegal block, so it shows up in none of the k-space checks and the
composite has to bound it explicitly. The train and the first interval cross between 250 and
300 Hz/px, and above that the first interval binds.

**`sc.compile` does not check RF amplitude at all**, and a 180 is where that first bites: peak B1
scales with the flip angle, so a 2 ms sinc 180 at TBW 4 asks 130 % of a 20 µT `max_b1`.
`Refocusing` refuses it with the duration that fixes it. `Excitation` has the same hole and is
deliberately left alone — widening a new module's change to patch an old one is how a one-module
change becomes five.

### What the simulator could and could not settle

`examples/se_2d/02` measures the headline: the spin echo recovers the phantom's **T2** — 80.0 and
100.0 ms fitted against 80 and 100 — and the same sweep without the 180 recovers **T2\*** at 21.8
and 28.6 ms. `examples/fse_2d/02` measures the point-spread width per ordering and separates
contrast from blurring.

**And it found a ghost that every arithmetic check in `01` passes through.** Shot 0 starts from
thermal equilibrium and every later shot starts from whatever the last TR recovered to. With
interleaved segmentation shot 0 owns *every eighth k-space line*, so that single difference is a
periodic modulation of `ky` — and a periodic modulation is a **replica of the object**, at FOV/8.
Nothing in `01` can see it: the k-space *positions* are exact to 10⁻⁷ 1/m, and it is the
*amplitudes* that differ. It is absent from HASTE, because one shot cannot be inconsistent with
anything, and that asymmetry is what identifies it.

Measured on a **one-voxel** phantom, because a ghost is a modulation of `ky` and the k-space of a
point object *is* that modulation: shot 0's lines come back **11.5 % brighter** than the other seven
shots', and the FOV/8 replica sits at 0.0149 of the peak. One dummy shot takes those to **−1.8 %**
and **0.0057**. `examples/fse_2d/01` therefore builds every multi-shot file with `DUMMY_SHOTS = 1`,
and `02` reads the comb back off `ky` to confirm it. It is the same effect
`examples/mprage_2d/02` found for a segmented GRE, one sequence over.

It was also confirmed on the object -- **0.2313 against 0.0301**, a factor 7.7. The one-spin number
is what a notebook can afford to *measure*, because it costs 1.3 s
against 66 s and has none of the failure modes an image-space metric has: three separate attempts to
read this ghost off a picture were wrong before one was right, first because ESPIRiT maps are ~0
exactly where a replica lands, then because the object's own edge dominated the region, then because
thermal noise at SNR 60 is the same size as the ghost.

**And the measurement replaced a simulation that had to go.** An earlier draft of
`examples/fse_2d/02` reconstructed five 16-echo, 9-shot acquisitions against a 128 x 128 x 4 phantom
in a single cell -- 173 s of saturated CPU each, sixteen threads by default, eight minutes with no
output until it finished. It took a 32-core workstation down three times. Peak memory was 1.2 GB, so
it was load and not RAM; and cutting `max_state_count` was not available, because 200 -> 32 saves
half the time and is 5 % wrong -- a CPMG train's stimulated pathways are the sequence.

What made the images affordable was the slab: `gre_2d/02` uses four slices because a phantom one
voxel thick returns exactly zero signal under a slice-select gradient, and **two** turns out to
reproduce the four-slice replica to 4 % (0.0301 against 0.0289, ratio 7.7 against 7.8) for half the
voxels. So `fse_2d/02` draws all three orderings in **4 minutes** at
`SIM_THREADS = 4`, against eight at sixteen threads for fewer pictures -- 5.8x less work at a
quarter the instantaneous load. Both `02` notebooks now set their thread count **before importing
torch**, print their budget in the first cell, free each phase graph after use, and keep the
expensive part behind a named switch.

Two things are **named rather than worked around**. The `z` **crusher balance** is not simulable
here: a slab four voxels thick has nothing for a slice gradient to dephase across within a voxel, so
an unbalanced crusher simulates as perfectly fine and is wrong on a scanner — which is why it is
asserted arithmetically, with the alternating failure mode as its own test. And MRzero's phase-graph
model is **insensitive to the CPMG phase relation**: sweeping the refocusing phase from 0° to 90° at
`B1 = 0.8` changes the echo magnitudes by 6 × 10⁻⁸. Both notebooks say so where a reader would
otherwise assume coverage.

### Measured

| | |
|---|---|
| worst error in k at the echo, all three axes | **2.0 × 10⁻⁷ 1/m** over 8 and over 16 echoes, **6.0 × 10⁻⁸** over a 72-echo HASTE shot |
| `Refocusing` time symmetry | **0 ns**; area to centre against area from centre, **600.000000** both |
| refocusing pulse spacing uniformity | **0 ps**; echo displacement from the midpoint **+1.95 µs**, within one gradient raster |
| minimum echo spacing, reference protocol at 200 Hz/px | **10.70 ms**, train-bound (train / first interval / ADC clearance = 10.684 / 9.320 / 10.616 ms) |
| not sharing one crusher window across three axes | **+1.4 ms per echo, 13 %**, and 49 pulseq blocks instead of 19 |
| plateau padding to symmetrise about the effective centre | **70 µs**, exactly the difference between the dead time and the ringdown |
| pulseq blocks | **5** for a single-echo spin echo, **296** for 8 shots × 16 echoes, **293** for a 72-echo HASTE shot |
| scan time, 128 lines | **4.3 min** at one echo per TR, **18 s** at sixteen with a run-in, **0.8 s** for HASTE |
| shot-to-shot comb depth, interleaved | **+11.5 %** with no dummy shots, **-1.8 %** with one |
| the FOV/8 replica it puts in the point spread | **0.0149** against **0.0057** |
| `examples/fse_2d/02`, end to end | **4 min** at 4 threads with every image drawn, or **15 s** with `BRAIN_IMAGES = False`; the draft that crashed took ~8 min at 16 threads for less |
| one-spin measurement against one brain image | **1.3 s** against **66 s**, same conclusion |

---

## Unreleased — `IRPrep`, and an inversion time placed to the microsecond

MPRAGE and MP2RAGE, and the smallest thing they needed: **one new module**, one new method, and one
argument that moved.

```python
inv = sc.modules.IRPrep(opts=opts, thickness_mm=None, spoil_voxel_mm=5.0)
gre = sc.modules.GRE2D(opts=opts, fov_mm=220.0, matrix=(128, 128), thickness_mm=5.0, flip_deg=8.0)

t_train = inv.time_to_center() + ti_s - gre.time_to_center_line(lines=segment)
```

| | |
|---|---|
| `IRPrep` | an inversion pulse and the crusher after it, in `modules/preparation/` — a new folder, because `rf.use == 'inversion'` is not `rf/`'s membership rule |
| `GRE2D.time_to_center_line(lines=…, dummies=…)` | block start to the readout of `center_line`. Depends on the ordering, so it is a method taking the same arguments `build` does |
| `GRE2D.build(dummies=…)` | **moved from `__init__`.** `dummies` describes the acquisition, not a waveform, and the timing query above needs it beside `lines` |
| `GRE2DTR(spoil_axis=…)` | spoiling is `('x', 'z')` by default rather than `'z'` alone |
| `shift_slice` | promoted from `rf/excitation.py` to `modules/_support.py`, now that a second folder needs it |

### Two sequences that deliberately did not become modules

`MPRAGE2D` and `MP2RAGE2D` are written in `examples/mprage_2d/01_build.ipynb` and
`examples/mp2rage_2d/01_build.ipynb` **and stay there.** Each has exactly one consumer — its own
notebook — and a module with one consumer belongs where that consumer is. `IRPrep` shipped because
two notebooks use it; `GRE2D` already shipped and now has a second consumer.
`tests/modules/test_mprage_notebooks.py` runs both build notebooks and asserts against the classes
they defined, so neither can drift without CI noticing.

### Four things the extraction found

**An inversion time is measured from the pulse's effective centre, and that is 5.1 ms into a 10 ms
hyperbolic secant.** `IRPrep.time_to_center()` is load-bearing rather than bookkeeping: referencing
TI to the block start instead is a 5 ms error in the one quantity the sequence exists to control,
and the crusher after the pulse makes "half the block" wrong by a different amount again.

**Where k = 0 sits inside a train is not a constant.** Centric ordering acquires it first, linear
mid-train, and the same TI places those two trains most of a train-length apart —
`time_to_center_line` is a method for that reason, and its `dummies` term is the one that gets
dropped: the block starts at the first dummy, so omitting it under-reports by `dummies · TR`
silently.

**Spoiling on one axis is not spoiling.** A gradient dephases only along its own direction, so
winding cycles on `z` alone leaves the residual coherent in `x`. The readout axis is the one most
often missed, because after the echo only half the readout's area is left on `kx` — half a cycle
across one voxel. `spoil_axis` now defaults to `('x', 'z')`, with cycles counted per axis across
that axis's own voxel, which differ by a factor of three between the slice and the in-plane voxel.

**Interleaved segmentation turns a shot-to-shot difference into a replica of the object.** Each shot
owns a comb in k-space, so an amplitude difference between shots is a periodic modulation and
therefore a ghost at FOV/n_shots. Measured on the example scan: 296 % spread across shots with no
dummy shots, 48 % with one, 8.7 % with three. That is what `dummy_shots` is for, and it is *not* a
spoiling problem — turning the extra spoiler axes off moves the image background by 0.01 %.

### What the simulator could and could not settle

`examples/mprage_2d/02` measures the null point at `TI = T1·ln2` through a real shot, to a few
milliseconds, for two tissues — which validates the inversion, `time_to_center` and the TI
placement in one measurement. It also found that **MRzero applies a pulse as an instantaneous
rotation by its integrated envelope**, so an adiabatic inversion does not invert there at all: a
10 ms hyperbolic secant arrives as 289°. `IRPrep` keeps `'hypsec'` as its default because that is
right on a scanner; the examples pass `'block'` and say why, and `time_to_center` absorbs the
4.5 ms difference without any other number in the timeline moving.

### One thing the compiler caught on its own

Pulseq labels are stateful, so emitting `SET=0` once at the start of an MP2RAGE and `SET=1` before
each second train leaves every later shot's first train wearing the previous shot's `SET=1` — half
the data in the wrong contrast, in a file that is perfectly legal. `sc.compile` refuses it: 96
readouts share one `(LIN, SET)` address, found by the same check that catches a repeated
phase-encode line. `examples/mp2rage_2d/01_build.ipynb` builds the broken version to show it.

## Unreleased — `sc.modules`, extracted from a 2D GRE

seqcraft ships concrete modules again. Six names, all of them extracted from one working sequence
rather than designed:

```python
gre = sc.modules.GRE2D(opts=opts, fov_mm=220.0, matrix=(64, 64), thickness_mm=5.0)
sc.compile(gre(lines=range(64)), opts, name='gre_2d').write('gre_2d.seq')
```

| | |
|---|---|
| `Excitation` | an RF pulse and, when selective, its selection gradient and rephaser |
| `PhaseEncode` | one Cartesian phase-encode blip, designed once and scaled per line |
| `CartesianLine` | prephaser, readout gradient and ADC as one design |
| `spoiler` | *n* turns of phase across a voxel — a function, not a class |
| `GRE2DTR` | one repetition of a spoiled 2D gradient echo |
| `GRE2D` | the complete scan |

The previous library — 27 classes, 5 762 lines — was deleted rather than migrated
([ADR-003](docs/adr/003-scanner-and-module-reform.md)). What is here now came back under a rule the
old one had no way to satisfy: **write the sequence in raw pypulseq first, simulate it until the
image is right, and only then extract the module with the compiled output held fixed.**
`examples/gre_2d/01_build.ipynb` still builds the same sequence out of raw events beside them, and
`tests/modules/test_notebook_matches_the_package.py` asserts that the class the notebook writes and
the class the package ships produce identical events.

### Four things the extraction found

**The prephaser is the integrated moment up to the echo.** `-gx.flat_area / 2` drops the ramp-up
entirely; `-gx.area / 2` includes it and is still half a Δk short, because k = 0 belongs at sample
index `matrix // 2`, whose centre time is half a dwell after the midpoint of the ADC window.
`CartesianLine` integrates the gradient's own knots instead, which produces both terms for free and
generalises to ramp sampling and partial echo with no second derivation.

**The receiver is phase-locked to the transmitter.** `GRE2DTR` gives the same `phase_deg` to the
excitation and to the readout. Leaving the ADC at zero while the RF runs an RF-spoiling schedule
writes that schedule's quadratic phase into `ky`: measured on a single off-centre voxel, it
scattered a point source across the whole phase-encode direction and moved its peak by thirteen
pixels, with the readout direction perfectly correct beside it.

**The winder is three axes wide.** The slice rephaser on `z`, the phase-encode blip on `y` and the
readout prephaser on `x` all play at once. Waiting for the rephaser before starting the other two
adds its whole duration to TE — 520 µs on the example protocol, out of 3.9 ms — and buys nothing.
Each leaf reports its minimum and accepts an override; the kernel takes the maximum and passes it
down.

**pypulseq's slice rephaser is already correct**, including for a minimum-phase SLR pulse whose
effective centre is at the end of the waveform. `Excitation` was written expecting to compute its
own and does not; the invariant is asserted in `tests/modules/test_excitation.py` against an
asymmetric pulse, where "half the total area" is wrong by a factor of sixty.

### The layout, and what it forbids

Five role folders (`rf/`, `encoding/`, `readout/`, `kernel/`, `imaging/`) plus the top level, each
with a membership rule, and a **flat** re-export so no import path names a folder.
`tests/modules/test_layout.py` asserts the rules, the flat re-export, and that `modules/` imports
nothing downstream of `design` — no compiler, no `analysis`, no `display`, no `scanner`. The
compiler's own tests still use no modules at all.

### Also

- `pp.scale_grad` replaces the hand-written `derive(amplitude=…, area=…, flat_area=…)` in
  `design/module.py`, `docs/writing_a_module.md` and `tests/module/test_module.py`. A trapezoid
  carries three numbers that must move together and an arbitrary gradient four; writing them out is
  that many chances to update all but one, and the result still compiles.
- `make_slr_pulse`, `make_gauss_pulse` and `make_block_pulse` join `_compat._REQUIRED_TOPLEVEL`.
- `examples/_parked/` and `examples/lib/` are removed. The two DTI notebooks were kept as the
  acceptance test for whatever module set came next; that job is discharged, they are in git
  history, and the physics they depended on is still in `salvage/`. `mr0_bridge.py` went with them:
  `examples/gre_2d/02` simulates the written `.seq` through `mr0.Sequence.import_file`, which tests
  the file a scanner would play rather than the tree that produced it.

### One thing to know about the simulation

MRzero's `.seq` importer parses each event's `freq` field and then ignores it. A slice offset or a
readout FOV shift is therefore **invisible** to a simulation through a written file — both are
frequency offsets. `examples/gre_2d/02` is a single slice at isocentre for that reason, and
`CartesianLine` and `Excitation` say so where they implement the offsets.

---

## Unreleased — compile returns a `pypulseq.Sequence`

`sc.compile(tree, opts)` returns a `pypulseq.Sequence`. Not a wrapper, not a pair of a sequence and
a report. **If the tree cannot become a legal sequence it raises; if the compiler had to change a
waveform to make it legal it warns.** There is nothing to unpack and nothing to remember to check.

```python
seq = sc.compile(tree, opts)          # was: out = sc.compile(...); out.check().raise_if_failed()
seq.write('gre.seq')                  # was: out.write(path) -> WriteResult
```

The old shape had one failure mode and it was the same one `pSeq_Base.get_report()` had, one
indirection later: findings on an object nobody has to look at. `get_report()` printed and returned
`None`; `CompiledSequence.check()` returned a `Report` that a caller could simply not call. Either
way the sequence is written, and the console refuses it an hour later.

**Nothing about the emitted bytes changed.** `build_gre` and `build_se` write `.seq` files with the
same sha256 as before the revision began, and every structural field of
`tests/baselines/compiler_phase0.json` — block durations, event content, provenance digest, m0, m1,
m2 — matches. That is checked at three gates and asserted on every test run.

### Every legality failure raises

| Was | Is |
|---|---|
| a `grad_limit` / `slew_limit` `Issue` | `HardwareLimitError`, on the **summed** waveform |
| an `adc_samples_limit` / `rf_samples_limit` `Issue` | `HardwareLimitError` |
| a `label` `Issue` (duplicate k-space address) | `CompileError` |
| a `timing` `Issue` from `check_timing` | `CompileError`, less the `TotalDuration` artifact |
| a `duration` / `moment` / `address` `Issue` | `CompilerContractError` |

Each message names the offending number, the tag path it came from, the time it happens, and two
remedies with the values already computed:

```text
HardwareLimitError: slew 189% of the 150 T/m/s limit on axis x.
  from   :  tr.readout.prephaser
  at     :  2.340 ms (block 117)
  reached:  283.5 T/m/s
  fix
    lengthen the lobe, or lower the readout bandwidth
    or design that part against sc.opts.derate(opts, slew=0.52)
```

`HardwareLimitError` was declared and exported but raised nowhere; it is live now.

### Everything else is a `SeqCraftWarning`

A `UserWarning` subclass, through the standard `warnings` machinery — so `simplefilter('error')`
makes any of them fatal and `pytest.warns` asserts on them. Five categories: `merge`, `resample`,
`snap`, `orphan_label`, `norm`.

**One aggregated warning per category**, naming the count and the sites. Python's default filter
shows a warning once per unique `(message, category, module, lineno)`, so one `warn` per merge
would print the first and silently swallow the other eleven — strictly worse than the count the
report gave. A 64-line acquisition merging the same two modules 64 times is one fact, not
sixty-four, so identical entries are counted rather than repeated.

```text
SeqCraftWarning: 8 same-axis gradient merges: tr.rewinder+tr.prephaser (axis x) x8
```

### New — `seqcraft.analysis`

Four functions, one entry shape: **give it a tree, get numbers back.** You should not have to know
that PNS prediction needs a compiled sequence while a moment does not, so none of them asks for
one; `kspace` and `pns` compile internally.

| | |
|---|---|
| `sc.sample(tree, opts)` | the tree on a uniform grid. **Lossy**, on purpose — for looking |
| `sc.moments(tree, order=0)` | m0/m1/m2 per axis, from exact knots. Takes no `Opts` and does no compile |
| `sc.kspace(tree, opts)` | the trajectory at true ADC sample times, with a named return |
| `sc.pns(tree, opts, hw)` | the full SAFE curve — `ok`, `peak`, `norm`, `components`, `t` |

`moments` is not built on `sample`, and the reason is subtler than "sampling is lossy". Linear
interpolation errs *antisymmetrically* about each knot, so the halves cancel under the integral and
m0 comes back bit-identical while the peak is visibly rounded off. A `moments` built on `sample`
would look correct on every test anyone would think to write.
`tests/analysis/test_analysis.py` asserts both halves of that.

`sc.pns` now returns the whole curve rather than a verdict: when `ok` is `False`, `peak` says how
much but `norm` and `t` say *where*, which is what you need to fix it.

### Removed

| Removed | Use instead |
|---|---|
| `sc.CompiledSequence`, `sc.WriteResult` | the returned `pypulseq.Sequence` |
| `out.seq` | the returned object itself |
| `out.n_blocks` | `len(seq.block_events)` |
| `out.duration_s` | `seq.duration()[0]` |
| `out.definitions` | `seq.definitions` — set by the compile, not by `write()` |
| `out.check()`, `out.report`, `sc.Issue`, `sc.Report`, `sc.ReportFailed` | the compile raised, or warned |
| `out.write(path)` | `seq.write(str(path))` |
| `out.origin(i)` | gone; provenance paths still name the source in every error message |
| `out.moments(1)` | `sc.moments(tree, 1)` |
| `out.kspace()` | `sc.kspace(tree, opts)` |
| `out.pns(hw)` | `sc.pns(tree, opts, hw)` |
| `sc.plot_sequence(out)` | `seq.plot()` — pypulseq's own, better maintained |
| `sc.plot_trajectory(shots)` | `sc.kspace(...)` plus three lines of matplotlib |
| `sc.testing.*` | the compiler raises; for purity see [`docs/writing_a_module.md`](docs/writing_a_module.md) |
| `sc.provenance`, the JSON sidecar | **nothing yet** — see below |
| `sc.UnitSanityError` | raised nowhere; it guarded `Geometry`'s ranges, and `salvage/geometry.py` has its own |
| `seqcraft.compiler.definitions` | two sources need a collision check, not a merge algorithm |

> **The provenance sidecar is a real gap.** Nothing currently records the git commit, dirty flag or
> package versions that produced a `.seq`. It went with the result wrapper deliberately, so that
> what replaces it starts from a blank slate rather than from code shaped around a type that no
> longer exists.

### Moved — an exception lives with the code that raises it

Only what more than one package raises stays in `seqcraft.errors`. **Every one is still
re-exported at the root**, and `tests/test_layering.py` asserts the *identity* rather than the
presence, so `except sc.CompileError` is unchanged.

| Exception | Now defined in |
|---|---|
| `CompileError`, `HardwareLimitError`, `DefinitionConflict` | `seqcraft.compiler.errors` |
| `RasterError` | `seqcraft.design.timing` |
| `UnknownFieldError` | `seqcraft.scanner.opts` |
| `SeqCraftError`, `ConfigurationError`, `MissingExtraError`, `SeqCraftWarning` | `seqcraft.errors` |

### The layers

```text
errors  ─►  design  ─►  compiler  ─►  analysis  ─►  display
   └──────────►  scanner  (independent of all five)
```

`result/`, `report.py`, `testing.py` and `design/sampling.py` are gone. The edge the layering test
now forbids in particular is `compiler → analysis`: the compiler must not reach for the
measurements taken *of* what it produced. It would be an easy one to add — the self-check wants a
moment per axis — and it would be wrong, because that moment walks the tree while the self-check
needs the compiled side.

### Also

- **The synthetic PNS model says what it is.** `synthetic_hardware()` carries `is_synthetic=True`
  and a `repr` reading `NOT a real scanner; never use it to clear a human scan`. The caveat used to
  live in a docstring that this revision deleted, and it is the one safety note in the package.
- **`tools/capture_compiler_baseline.py` works again.** It monkeypatched `compiler._place` and
  `._axis_gradient`, neither of which has existed since the compiler was split into stages, so it
  had been raising `AttributeError` rather than arbitrating anything.
- **A real defect surfaced.** `make_arbitrary_grad` without `first=`/`last=` extrapolates from the
  end samples, so a "lone spiral" test tree started at −2680 Hz/m — a step from zero the amplifier
  cannot play. It compiled silently before; `check_timing` inside the compile now stops it.
- **An off-by-one in an error message.** The event-size check indexed a 0-based `origins` list with
  pypulseq's 1-based block id, so it named the next block along and could raise `IndexError` on the
  last one.

---

## Unreleased — four packages, and 600 fewer lines

The package was one 4 700-line `core/` directory whose membership rule — *"what is required to get
a logic block to a validated `.seq`"* — is a property of the package rather than of any file in it.
It admitted the scheduler, the unit table, the geometry, the report type and the exception hierarchy
alike, and gave no reason to keep anything out. It is now four packages named for four questions,
ordered by the one direction the dependencies run.

```text
scanner/     what you build against     compiler/    the transform
design/      what you build             result/      what compile returns

errors, report  ─►  design  ─►  result  ─►  compiler
```

`errors.py` and `report.py` stay at the root as a pair — hard failures raise, soft findings report
— because both are cross-cutting and neither is *produced* by any one layer. `Report` in particular
is the vocabulary the compiler writes findings in: five of its stage modules build `Issue` objects
long before a `CompiledSequence` exists, so filing it under "what compile returns" would have made
`compiler/` import `result/` for a type.

`tests/test_layering.py` asserts that order per file, read from the source rather than from
`sys.modules`, so a `TYPE_CHECKING`-only import still counts and one satisfied by import order does
not. Nothing in this section changes what the compiler produces: the frozen baseline
(`tools/capture_compiler_baseline.py`) matches block for block.

### Moved — every import path

| was | is |
|---|---|
| `seqcraft.core.logic` | `seqcraft.design.logic` |
| `seqcraft.core.events` | `seqcraft.design.events` |
| `seqcraft.core.timing` / `.units` / `.geometry` | `seqcraft.design.timing` / `.units` / `.geometry` |
| `seqcraft.module` | `seqcraft.design.module` |
| `seqcraft.core.compiler` (+ `core._compiler`) | `seqcraft.compiler` |
| `seqcraft.core.report` | `seqcraft.report` |
| `seqcraft.provenance` | `seqcraft.result.provenance` |
| `seqcraft.core.errors` | `seqcraft.errors` |
| `seqcraft.core.validate` | split: `seqcraft.compiler.definitions` + `seqcraft.design.geometry` |
| `CompiledSequence`, `WriteResult` | `seqcraft.result` |

**The top-level names are otherwise unchanged.** `sc.compile`, `sc.LogicBlock`, `sc.Module`,
`sc.Raster`, `sc.convert`, `sc.opts`, `sc.hardware`, `sc.display`, `sc.provenance`, `sc.testing`,
`sc.events`, `sc.timing`, `sc.units` all still spell the same way, and `seqcraft` remains the only
global re-export layer. Only code that reached past it needs editing — plus `sc.Geometry`, below.

- **`sc.sample(tree, opts)` is new** — a tree as arrays, gradients per axis plus the RF/ADC spans.
  It was `display._sample`; turning a tree into numbers is useful without drawing it, and moving it
  out means nothing numeric is left behind the matplotlib import.
- **`sc.validate` is gone.** `merge_definitions` is compiler-internal; the plausibility bands went
  with `Geometry`, below.

### Removed — `sc.Geometry`, and `compile(geometry=)` with it

The compiler now takes pulseq's own definition keys and nothing else:

```python
out = sc.compile(tree, opts, geometry=geometry)                       # was
out = sc.compile(tree, opts, definitions=geometry.definitions())      # is
```

`geometry=` existed to call `definitions()` on a `Geometry` and merge the eight keys that came
back — which is what `definitions=` already did, for any source, with no dataclass in between. So
~450 lines of dataclass plus range framework sat inside the package to produce one dict. FOV, matrix
and slice order are decisions about the *scan* you are running; the compiler turns a tree into legal
pulseq blocks and is indifferent to why the tree looks the way it does, and keeping `Geometry` in
the package made it look like a required input it never was.

- `CompiledSequence.geometry` is gone. `CompiledSequence.definitions` is unchanged and is now the
  only place scan metadata lives, which is the point: what the file says and what it plays come from
  one mapping.
- **The class is preserved whole** in `salvage/geometry.py`, standalone — no seqcraft imports, its
  two error types plain `ValueError` subclasses — so it can be copied into a module library or into
  user code as it stands. Unlike the rest of `salvage/` it is *expected back*: a geometry of that
  shape is wanted the moment the module library gets its infrastructure, and this is the design to
  start from.
- The two tests asserting the range bands' unit names are ones `convert` knows went with it. They
  guarded against one package growing two spellings of a unit; with the bands outside the package
  there is no shared vocabulary left to keep in step. `salvage/geometry.py` records the check so
  whoever adopts it can re-establish it.

### The compiler is seven files instead of one 1 724-line module

`compile_sequence` is now the pass that orchestrates them, and each stage is a module with a name:
`placement` → `boundaries` → `legalization` → `emission`, with `verification`, `model` and
`definitions` beside them. The `_Placed` / `_GRAD` / `_place` alias shim is gone; the real names are
at the ~180 use sites.

**`CompiledSequence._verify` moved to `compiler.verification.verify_against_tree`.** It took
`Sequence[PlacedEvent]` — the compiler's *private IR* — while living on the result type, and that
one method was the only reason `result/` depended on `compiler/`. It is a compile stage that happens
to run last, not an accessor, and moving it removed ~150 lines from the result type and the cycle at
once.

**`in_block_delay` lives in `compiler/model.py`**, not in `emission.py` as first planned:
legalization needs it to delay a gradient that passes through untouched and emission needs it for
every RF, ADC and label, so putting it in either stage would make the two import each other.

### Removed — 543 lines with no reachable consumer, and 59 quarantined

Each of these was found by grepping `src/`, `tests/`, `examples/`, `tools/` and `salvage/` for a
caller and finding none.

- `events`: `sanitise`, `duration_of`, `kinds_of`, `channels_of`, `moment_of` — four of nine public
  functions had no caller and no doctest, and a fifth had only a doctest.
- `report`: `merge`, `combine`, `to_dict`. `to_dict`'s docstring claimed it was "for the provenance
  sidecar"; `provenance.py` never called it.
- `errors`: `PurityError` — documented as "raised by the test helpers", which never raised it.
- `_compat`: `probe`, `has`, `rotate_3d`, `_rotate_3d_available`, `_OPTIONAL`, `supported_rf_uses` —
  the entire optional-capability chain. `require()` does its own `hasattr` checks and never called
  `probe()`. The module now does exactly one job: fail once, at import, with a complete list.
- `geometry`: `round_half_up`, `dk_per_m`, `res_mm`, `slice_position_m`, `check_module`, `params`,
  `describe`.
- `testing`: `assert_block`, `assert_raster`, `assert_compiles`, `module_subclasses`. The first two
  check what the compiler already refuses with a better message — `LogicBlock.add()` rejects a
  malformed node at construction, and an off-raster gradient start raises naming the nearest raster
  point above *and* below plus two fixes. `assert_compiles` was one line, now inlined into
  `assert_output`. Nine public assertions become five.
- `display`: `plot_kspace` — it was `plot_trajectory([compiled.kspace()['k_adc'][:2]])` plus
  decimation and a title. `plot_trajectory` absorbed both and gained `max_points=`.

**Quarantined to `salvage/geometry_pe.py`:** the phase-encode table (`pe_first_index`, `pe_lines`,
`n_pe_acquired`, `pe_lines_for_shot`, `par_first_index`, `par_lines`). `geometry` was defended as
holding the *one* index computation `kspace_center_line` and the `LIN` label would share — but the
sharing never happened, because the module library that would have consumed it was deleted. What is
worth keeping is one non-obvious line, `skip += (c - skip) % r`, without which an even skip and an
odd centre miss k = 0 entirely; it is preserved with the reasoning. `kspace_center_line` stays,
because `definitions()` writes it.

### Fixed — three stale references

- `pyproject.toml` declared `[project.scripts] seqcraft = 'seqcraft.cli:main'`. There is no
  `cli.py`, so `pip install` shipped a console script that raised `ModuleNotFoundError`.
- Two places named `sc.opts.from_specs()`, which has never existed; the name is `from_scanner`.
- The `viz` extra required `seqeyes-python>=0.2.9`, which nothing imports. Dropped; it belongs there
  the day a module in `src/` actually imports it.

### Corrected — a documented invariant that was not true

`import seqcraft` was documented as not importing matplotlib. It does, and always did:
`pypulseq/Sequence/calc_grad_spectrum.py` imports matplotlib at module level, so *any* import of
pypulseq pulls it in. What seqcraft can promise, and what is now asserted, is that `display.py` is
the only file here that imports it and that `display` and `testing` stay behind `__getattr__`.

---

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
