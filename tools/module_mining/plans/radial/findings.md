# RadialReadout Fine Scan — Code Archaeology and Semantic Decomposition

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** `RadialReadout`
**Date:** 2026-09-18
**Passes:** 1 — code archaeology; 2 — semantic decomposition; 3 — the boundary questions
**Status:** evidence only. **No class has been written.**

References: [`reference_inventory.md`](reference_inventory.md). Measurable contract:
[`invariant_table.md`](invariant_table.md). Measurements reproduced by
`tools/module_mining/candidates/radial/run_archaeology.py`.

---

## 1. What the user configures

| | R1 / R2 | R2b fast | R3 OpenMRF | R4 UTE |
|---|---|---|---|---|
| geometry | `fov`, `n_x` | `fov`, `Nx`, `ro_os` | `FOV.fov_xy`, `Nxy`, `ro_os` | `fov`, `n_x`, oversampling |
| readout | flat time fixed in body | `ro_dur` | `adcTime` | `readout_duration` |
| spokes | `n_spokes` | `Nr` | `NR` | `n_spokes` |
| angle | `pi/n_spokes` in the loop | same | **`phi_mode`**: equal, random, golden | `pi/n` in the loop |
| asymmetry | — | — | — | **`readout_asymmetry`** |
| spoiling | area constants | `ro_spoil`, `sl_spoil` | `spoil_ro`, `spoil_sl`, `mode_rewind` | fixed area |
| contrast | `te`, `tr`, `flip` | derived minimum | config | `tr`, derived TE |

Two of these are the whole fine scan: **`readout_asymmetry`** (§5) and **`phi_mode`** (§6).

## 2. Which calculations encode MRI physics

1. **Sample spacing from the FOV.** `deltak = 1/fov`, readout flat area `n_x * deltak`.
   Universal across all four.
2. **The prephaser, which is where the references differ** — three formulations of one
   requirement, that a sample land exactly on `k = 0`:

   | | |
   |---|---|
   | R1 / R2 | `-gx.area/2 - deltak/2` |
   | R2b | `-gx.amplitude * (dwell*(n_samples/2 - 0.5) + 0.5*gx.riseTime)`, commented *"0.5 is necessary to account for the Siemens sampling in the center of the dwell periods"* |
   | R3 | `-Nxy*deltak/2 - (gx.area - gx.flatArea)/2 - kx_corr`, where `kx_corr` is **measured** by compiling a probe sequence and reading `calculateKspacePP` at the centre sample |

   All three are after the same number. R1's `-deltak/2` and R2b's `-0.5` dwell term are the same
   half-sample shift written two ways; R3 declines to write it at all and measures the residual
   instead. **Three independent derivations agreeing on one invariant is the strongest evidence
   this scan has**, and it is why the invariant is stated as a measurement (R1 in the table)
   rather than as a formula.
3. **Rotation.** R1, R2, R2b and R3 use Pulseq's block-level `rotate` helper. R4 does not: it
   copies the designed trapezoid and scales its amplitude by `cos(phi)` and `sin(phi)` onto `x`
   and `y`. Since a trapezoid's ramps are unchanged by that scaling, the two are equivalent for
   trapezoids — and **R4's is the only one a SeqCraft module can use**, because `rotate` is a
   pulseq block helper and the compiler is out of scope.
4. **Readout spoiling by extending the flat top** (R2b, R3): the readout gradient simply runs on
   past `k_max`. Both implementations warn in a comment that the event's stored `area` is then
   wrong — a hazard a module should not reproduce.

## 3. Which code exists only because Pulseq is flat

Less than in the TSE references, and of a different kind:

- the per-block `rotate(...)` calls, which exist because a block is the unit a rotation can be
  applied to;
- R2b's `splitGradientAt` / `addGradients` fusing of the slice spoiler into the selection
  gradient — the same hand-stitching `Refocusing` already removes on the selection axis;
- R3's probe-compile-and-measure loop, which is a workaround for not being able to ask the
  design what k it will produce. **A module that exposes its own semantic `k = 0` makes that
  loop unnecessary**, which is an argument for question Q3 below.

## 4. What is acquisition policy

- the angle schedule — equal-increment, golden-angle, randomised (R3 offers six);
- the number of spokes, and dummy count;
- RF spoiling increment and its linear/quadratic form;
- view sharing and interleaving, where a protocol has them;
- which spokes are acquired versus played as dummies.

## 5. Semantic decomposition

```text
Excitation (+ slice rephase)
      |
      |  prephaser along the spoke  (x,y scaled by cos/sin)   } same block
      |  slice rephaser on z                                  }
      v
   TE delay
      |
      v
 [ readout gradient along the spoke + ADC ]    <- the candidate
      |
      v
 [ spoiler along the spoke + z spoiler + TR fill ]
```

| component | class |
|---|---|
| excitation and slice rephaser | `EXISTING_SEQCRAFT` (`Excitation`) |
| prephaser along the spoke | **`CANDIDATE_OWNED`** — its area is the readout's own geometry |
| readout gradient + ADC, oriented | **`CANDIDATE_OWNED`** |
| where `k = 0` falls among the samples | **`CANDIDATE_OWNED`** (Q3) |
| TE / TR timing | kernel, later — a `RadialGRETR` if one is justified |
| spoiler along the spoke | kernel: its **direction** is the readout's angle but its **area** is a spoiling choice, and only the layer holding both can place it |
| z spoiler, RF spoiling schedule | caller / kernel as today |
| angle schedule, spoke count, dummies | `CALLER_POLICY` |
| `rotate(...)` blocks, `splitGradientAt` fusion, probe-compile correction | `REFERENCE_ARTIFACT` |

## 6. The four boundary questions, answered from measurement

### Q1 — full spoke vs centre-out: one abstraction or two?

**One.** `write_ute.py` already spans the continuum with a single argument and a single formula,
and the trajectory measurements confirm it does so cleanly:

| `readout_asymmetry` | centre sample | signed k range (1/m) | Δk (1/m) |
|---|---|---|---|
| 0.00 | 64 / 128 | −128.0 … +126.0 | 2.000 |
| 0.25 | 48 / 128 | −76.8 … +126.4 | 1.600 |
| 0.50 | 32 / 128 | −42.7 … +126.7 | 1.333 |
| 0.75 | 16 / 128 | −18.3 … +126.9 | 1.143 |
| 1.00 | 0 / 128 | −0.0 … +127.0 | 1.000 |

The centre sample walks linearly from the middle to the first sample, `k_max` stays put, and Δk
scales to hold the resolution. Nothing changes kind at either end.

**And SeqCraft already has this parameter.** `CartesianLine.partial_fourier` does exactly this to
a Cartesian line: it moves the echo sample towards the start of the readout. The mapping is
`partial_fourier = 1 - asymmetry/2`, so a full spoke is `1.0` and centre-out is `0.5` — the same
convention, the same meaning, on a different trajectory.

*Recommendation:* one `RadialReadout` with `partial_fourier` in SeqCraft's existing sense. A
`centre_out=True` flag or a sibling `CenterOutReadout` would split a continuum that the evidence
shows is not split.

### Q2 — does the leaf own angle/orientation?

**Yes, as a per-call `build` argument.** The test is whether the component can determine the
value from its own contract: given the angle, the FOV, the matrix and the bandwidth, the readout
determines every event it emits, and nothing above it is consulted. That is exactly the relation
`PhaseEncode` has to `line` — designed once, scaled per call.

The *schedule* of angles is not the leaf's (Q4).

### Q3 — does the leaf own the semantic `k = 0`?

**Yes, emphatically, and R3 is the evidence.** OpenMRF has to compile a probe sequence and read
its trajectory back to discover where its own centre sample lands, because nothing in its design
will tell it. That is a layer reverse-engineering k-space from waveforms — the failure mode the
ownership rule exists to prevent.

`CartesianLine` already answers this for a line (`echo_sample`, `time_to_echo`,
`area_to_echo_per_m`), and the TSE pilot showed what it buys: the kernel asked the leaf for its
geometry instead of re-deriving it. A radial readout should answer the same questions.

### Q4 — spoke ordering and golden-angle schedules?

**Policy, above the leaf.** R3 makes the separation visible: `phi_mode` builds a list of angles
before any waveform exists, and `RAD_add` consumes one per repetition. The readout never needs
the list. This is the same conclusion the TSE pilot reached about `segments`, reached the same
way.

## 7. What the measurements say about the reference

`write_radial_gre.py`, 8 spokes of 64 samples at FOV 260 mm — all seven radial invariants pass,
most at the 1e-12 level. Full numbers in [`invariant_table.md`](invariant_table.md); the two
worth stating here:

- **a sample lands exactly on `k = 0`** (2.2e-12 1/m), which is the opposite of what the TSE
  fine scan found for the Cartesian demo, where the centre of k-space fell *between* two samples.
  In radial it must: every spoke shares that one point, and it is the most densely sampled point
  in the acquisition;
- **every spoke is spoke 0 rotated** (3.3e-12 1/m) — the metamorphic test, passing on the
  reference before any candidate exists, which is what makes it usable as an acceptance test
  afterwards.

The full spoke is **asymmetric by one sample**: with 64 samples the centre is index 32, so k runs
−32Δk … +31Δk. Any reconstruction, and any test, has to know that.

## 8. Open questions carried forward

1. Does `RadialReadout` own its prephaser, or should it be optional as `CartesianLine`'s is?
   A radial spoke has no meaning without one, which argues for owning it — but a caller wanting
   the prephaser fused with a slice rephaser needs its duration, exactly as
   `CartesianLine.prephaser_duration_s` provides.
2. Does the candidate emit the rotated events, or emit an unrotated design plus an angle that
   something else applies? R4 shows scaling works; the question is whose job it is.
3. Is a `RadialGRETR` kernel justified, or is radial GRE a composition? The spoiler's direction
   couples to the readout angle, which is real coupling — but it may be the only one.
4. Does the readout own `ro_spoil` (extending its own flat top past `k_max`)? R2b and R3 both do
   it, and both corrupt their own area bookkeeping doing so.
5. 3D radial and `rot3D` are out of scope here and should stay out until 2D is settled.
