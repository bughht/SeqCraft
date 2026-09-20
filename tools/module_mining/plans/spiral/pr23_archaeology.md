# PR23 archaeology — design decisions and failures, classified

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-20, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**PR #23** `feat/spiral-2d`, +7365/−15 over 17 files, written against `d27217d` (pre-PR #22 merge)
and now **36 commits behind `main`**. Read as a historical experiment, not a specification.

Classification per the brief: `PHYSICAL_REQUIREMENT` · `ARCHITECTURAL_REQUIREMENT` ·
`IMPLEMENTATION_CHOICE` · `WORKAROUND` · `FAILURE_EVIDENCE` · `UNPROVEN_CLAIM`.

---

## 1. Design-decision inventory

| # | decision | class | note |
|---|---|---|---|
| D1 | Path (Nyquist geometry) and traversal (slew/amplitude-limited time) are two passes | **PHYSICAL_REQUIREMENT** | see §3; the reversibility argument is independent of the API |
| D2 | Traversal has `v = 0` at **both** ends | **PHYSICAL_REQUIREMENT** | this, not the two-pass split, is what makes an arm reversible |
| D3 | Neither design pass imports pypulseq | ARCHITECTURAL_REQUIREMENT | keeps geometry testable without a scanner model |
| D4 | `density` is sampled FOV multipliers, not polynomial coefficients | **PHYSICAL_REQUIREMENT** | §2-B; the alternative expresses the opposite trend silently |
| D5 | `out` / `in` / `in-out` / `out-in` in one class via `variant` | IMPLEMENTATION_CHOICE | follows from D2 but is not forced by it; see §4 |
| D6 | `echoes` with no `polarity` | IMPLEMENTATION_CHOICE | mirrors `CartesianLine`'s monopolar/bipolar trade |
| D7 | `k_per_m` measured off **built events**, design `k` discarded | **PHYSICAL_REQUIREMENT** | §2-C; the strongest idea in the PR |
| D8 | ADC split into segments for the interpreter's sample limit | **PHYSICAL_REQUIREMENT** | hardware limit, already has `check_event_sizes` on main |
| D9 | An acquisition region is not an arm — one ADC spans joined arms | **PHYSICAL_REQUIREMENT** | §2-D; what keeps `k=0` inside a window |
| D10 | `adc_segments=1` refuses rather than splitting, naming the smallest `shots` | IMPLEMENTATION_CHOICE | good refusal design; not forced |
| D11 | `SET` label carries which ADC event of a split region | ARCHITECTURAL_REQUIREMENT | `check_label_addresses` requires it |
| D12 | No `REV` label emitted | ARCHITECTURAL_REQUIREMENT | "an emitted label with no reader is a promise nobody keeps" |
| D13 | Joins built **inside** the arm's arbitrary gradient as triangles | WORKAROUND | compensates for raster-centre vs raster-edge lattices; correct, but it is a workaround for a representation mismatch |
| D14 | Join area corrected on the **assembled** waveform by one rescale | **PHYSICAL_REQUIREMENT** | §2-E |
| D15 | Three helpers promoted to `_support` | UNPROVEN_CLAIM | §5 — one consumer each |
| D16 | `magnitude` = ∫\|g\| rather than \|∫g\| | UNPROVEN_CLAIM *as a defect* | `compiler_claims.md` §3 |
| D17 | m1 coefficient 1e-9 → 1e-7 | **CONFIRMED defect**, form disputed | `compiler_claims.md` §2 |
| D18 | Reconstruction in `examples/`, not `src/` | ARCHITECTURAL_REQUIREMENT | "a reconstruction module there would put an ESPIRiT dependency on the compile path" |
| D19 | One `Spiral2D` owns 13 concerns | **open** | §4 |
| D20 | No `Spiral3D`, no anisotropic FOV, no density presets | ARCHITECTURAL_REQUIREMENT | scope refusals, all well argued |

## 2. Failure inventory

The failures are the most valuable part of the PR, because each is a **silent** failure: the
sequence compiles and looks right.

| | what looked reasonable | why it was wrong | observable that detects it | still exposed on `main`? |
|---|---|---|---|---|
| **A** path/traversal conflated | integrate one ODE in time, as `vds.m` does | trajectory ends at full gradient, so its time-reverse *begins* at full gradient — a waveform no block can start with | first/last waveform sample ≠ 0 on a reversed arm | **N/A** — no spiral on main; would recur in any new implementation |
| **B** density as polynomial coefficients | `(1.0, 0.5)` reads like "half the FOV at the edge" | under `vds.m`'s convention it means FOV *grows* by half: oversampled outer ring, 60 % longer readout, nothing raises | measured turn spacing vs `1/density` | **N/A**, and it is an API-shape trap any new design can re-enter |
| **C** trusting design-space `k` | the design pass computed `k`, so report it | the NUFFT is *told* positions; `k` one raster from where it played is a blur, a shading and a wrong field map | design `k` vs `sc.kspace` on the compiled sequence | **live** — nothing on main forbids reporting a design trajectory |
| **D** splitting an acquisition wrongly | split per arm | `k=0` lands in the ADC dead-time guard between two events instead of inside a window | `\|k\|` at the nearest sample to the region's centre | **live** for any future segmented readout |
| **E** joins designed in isolation | each join exact on its own | the assembled piecewise-linear waveform carries area the isolated design does not count — measured **0.108 1/m**, 0.024 `dk`, invisible in every plot | end-of-block `k` vs origin | **live** |
| **F** spin-echo spiral placed on the wrong semantic sample | place the readout by its start or its midpoint | legal sequence, TE wrong by ≈ one arm | compiled RF-centre-to-`k=0` time vs reported TE | **live** — same class as the GRE3D TE error already found |
| **G** m0 scaled by net moment | "area" in the comment | sign-changing arbitrary gradients collapse the tolerance | `traversed / net` ratio on one event | **live**, but the failure **did not reproduce** — see `compiler_claims.md` §3 |
| **H** m1 tolerance at the round-trip floor | 1e-9 looks strict and safe | it sits *inside* the representation noise | sweep legal waveforms and count refusals | **live and CONFIRMED** — 6 of 10 refused |
| **I** arbitrary-gradient joins causing resampling | superpose a trapezoid on the arm | resample moves the trajectory silently while the sequence stays legal | `resample` warning; m0 shift of 0.046 1/m | **live** |

**C, E, F and H are the four that should become acceptance tests** for any Spiral extraction, and
F and H are not spiral-specific at all.
