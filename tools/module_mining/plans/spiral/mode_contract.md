# SpiralReadout — the variant mode contract

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-20, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Purpose:** to establish that `out`, `in`, `in-out` and `out-in` are four physical modes of one
readout, not four implementations behind a string flag. Written **before** implementation, and
required by rule A because each mode moves the semantic echo.

Reference coverage sets the **evidence tier** of a claim. It does not set the abstraction
boundary.

---

## 1. The arm, which is what makes a family possible

One **arm** is a continuous traversal of the density-modulated path between the origin and
`k_max`, designed once and **at rest at both ends** (`g = 0`).

That endpoint condition is the whole design. An arm at rest at both ends can be played backwards,
so `in` is `out` reversed; two arms can be laid end to end with no connector, because both are at
`g = 0` where they meet. Every mode below is a choice of *how many arms and in which order*.

## 2. The mode table

| | `out` | `in` | `in-out` | `out-in` |
|---|---|---|---|---|
| **path order** | 0 → k_max | k_max → 0 | k_max → 0 → k_max | 0 → k_max → 0 |
| **arms** | 1 | 1 (reversed) | 2 | 2 |
| **semantic k = 0** | **first** sample | **last** sample | the **seam**, mid-region | **both** ends — two instants |
| **k = 0 count** | 1 | 1 | 1 | **2** |
| **gradient start** | `g = 0` | `g = 0` | `g = 0` | `g = 0` |
| **gradient end** | `g = 0` | `g = 0` | `g = 0` | `g = 0` |
| **starts at k =** | origin | k_max | k_max | origin |
| **ends at k =** | k_max | origin | k_max | origin |
| **prephase** | not required | **required** | **required** | not required |
| **rewind** | **required** | not required | **required** | not required |
| **ADC topology** | one region | one region | **one region spanning both arms** | one region spanning both arms |
| **internal join** | none | none | at **k = 0** | at **k_max** |
| **join continuity** | n/a | n/a | both arms at rest → continuous, no connector | both arms at rest → continuous, no connector |
| **TE semantic point** | region start | region end | the seam | **two** — region start and end |
| **multi-echo transition** | flyback k_max → 0 | flyback 0 → k_max | **none** — consecutive arms share an endpoint | **none** |

## 3. Same invariant / mode-specific consequence / implementation detail

**Invariant across all four modes** — these are the contract:

```text
every arm begins and ends at g = 0
the path is the Nyquist condition and nothing else; the traversal is the only thing that
    knows the scanner
the reported trajectory is measured off the BUILT events, never the design pass
every instant at which the path crosses the origin is reported, and each falls INSIDE a
    sampling window, never in an ADC dead-time guard
turn spacing matches 1/density everywhere
the arm at angle phi is the arm at 0, rotated
```

**Mode-specific physical consequence** — these differ, and the table above is their statement:
the number and location of k = 0 instants; which endpoint needs a prephaser and which needs a
rewinder; where the internal join falls; whether a multi-echo transition costs a flyback.

**Implementation detail** — not part of the contract: whether the join is realised inside the
arm's arbitrary gradient or as a separate event; how ADC segmentation distributes samples across
segments; which pypulseq factory builds what.

## 4. The finding that constrains the first implementation

**`out-in` has two k = 0 instants.** Every other mode has one. A module whose semantic echo is a
scalar works for three modes and breaks on the fourth — and it would break by *changing meaning*,
which is exactly the temporary model that must not be built.

> **The durable contract requires the k = 0 instants to be a sequence from the first line of
> code**, even while `out` is the only implemented mode and its sequence always has length one.

Concretely: `echo_sample(i)`, `te_s` as a tuple, `time_to_center(echo=0)`. Not `echo_sample`,
not a scalar `te_s`. This is the single most important thing the mode table produced, and it is
free if done first and expensive later.

Second, smaller: **prephase and rewind are two faces of one rule.** Each mode needs exactly one of
them, decided by which endpoint is away from the origin. They follow the precedent this codebase
already has twice — `CartesianLine(prephase=False)` and `Excitation.build(rephase=False)` — so the
module **states the requirement and realises it by default**, and a caller who wants to own the
realisation says so. That is one rule covering all four modes, not a per-mode switch.

## 5. The endpoint policy, decided

```text
DURABLE POLICY: every arm begins and ends at g = 0, in every mode, including `out`.
```

`writeSpiral.m` ends its spiral-out at **full gradient** and ramps down inside the spoiler. That
is a coherent choice for a mode family of one, and it is **incompatible with reversibility**: an
arm ending at full gradient time-reverses into one starting at full gradient, which no block can
begin with. The reference has no spiral-in for exactly this reason.

So SeqCraft's `out` **deliberately differs from the reference at its tail**, and the consequence
for acceptance must be stated rather than discovered:

- **comparable against the reference:** path geometry, turn spacing, k-space extent, the
  limit-respecting property of the traversal, ADC segmentation against the sample limit;
- **not comparable:** the final samples, total arm duration, and any event-level similarity near
  the end.

Copying the reference endpoint to obtain event-level similarity would buy a stronger-looking
comparison for `out` and make the family internally inconsistent. It is refused.

## 6. Evidence tier, per mode

| mode | independent reference | analytic | metamorphic | emitted-trajectory | simulation |
|---|---|---|---|---|---|
| `out` | **yes** — `writeSpiral.m`, for path, extent, traversal, ADC segmentation | endpoint policy | rotation equivariance | required | deferred |
| `in` | **no** | reversal identity | **`in` is `out` reversed** — the primary check | required | deferred |
| `in-out` | **no** | seam continuity; k = 0 inside a window | reversal of the first arm | required | deferred |
| `out-in` | **no** | outer turnaround; two k = 0 instants | reversal of the second arm | required | deferred |

**Nothing but `out` may be described as validated against `writeSpiral.m`.** The other three rest
on analytic and metamorphic evidence, which is the tier `GRE3DTR`'s slab-selective path used and
is sufficient — a reversed piecewise-linear waveform's integral is computable without any external
implementation.

The metamorphic checks are the strong ones here, and they are definable **before** the module
exists: reverse the emitted `out` waveform, integrate it, and require the trajectory to be the
`out` trajectory reversed. That is a property of the design, testable against the reference's own
spiral-out arm.

## 7. Does it still fit one `SpiralReadout`? — yes

Against the flag-explosion criteria:

| | |
|---|---|
| one path/traversal model | **yes** — identical for all four; the mode never re-enters the solver |
| one hardware-limit solve | **yes** — the arm is designed once |
| one trajectory representation | **yes** — measured off built events, k = 0 instants as a sequence |
| one ADC ownership model | **yes** — one region per contiguous non-zero stretch, segmented against the sample limit |
| one set of correctness concepts | **yes** — §3's invariant list is mode-independent |
| do unrelated parameters become conditional? | **no** — `density`, `shots`, `dwell_s`, `fov_mm`, `matrix` are untouched by `variant` |
| different ownership rules per mode? | **no** — prephase/rewind is one rule evaluated per mode |
| effectively separate machines? | **no** — one arm, played in a different order |

**One public leaf, `SpiralReadout`, owning all four variants.** The internal boundary stays
`path → traversal → realization`, and remains internal.

---

## 8. Proposed implementation order

Sequence, not scope. The stage target is the four-variant family.

```text
1. the arm, and the durable shapes
       path -> traversal -> realization, internal
       echo instants as a SEQUENCE           <- must be right from the first line
       echo_spacing_s DERIVED, not stored
       prephase / rewind as stated requirements with default realisation
   validated by: the analytic invariants, and rotation equivariance against the reference

2. variant='out'
   validated by: writeSpiral.m for path, turn spacing, extent, traversal limits, ADC
                 segmentation.  NOT for the tail -- the endpoint policy differs on purpose.

3. variant='in'
   validated by: the metamorphic check that carries this whole family --
                 reverse the emitted `out` waveform, integrate it, and require the
                 trajectory to be `out`'s reversed.  Definable BEFORE the module exists.

4. variant='in-out', 'out-in'
   validated by: seam continuity; k = 0 inside a sampling window for in-out; the two
                 origin crossings reported for out-in.

5. echoes > 1
   validated by: echoes_contract.md -- transitions free for two-arm variants and a flyback
                 for one-arm ones, one ADC region per contiguous non-zero stretch.
```

Steps 3 and 4 need no new external evidence, which is why the absence of a reference for them is
not a blocker — it is a claim-scope statement.

**What must not happen at step 2:** a scalar echo time, a stored `echo_spacing_s`, or an endpoint
policy copied from the reference. Each would make step 3 a change of meaning rather than an
addition, which is the failure mode the incremental order exists to avoid.
