# Repetition-level joint design — spike findings, before any framework code

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-24, so that the pull request carries the direction behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../../README.md) for which one to edit.

**Status:** findings for review. No framework code written, no API proposed.
**Date:** 2026-09-24
**Reads:** [`2026-09-24_repetition_level_joint_waveform_design.md`](2026-09-24_repetition_level_joint_waveform_design.md)

Measured on `spike/repetition-joint-design`, branched from PR #39's head. Every number below
comes from building real events and integrating the emitted waveform.

---

## 1. The shape already exists, for one order

`GRE3DTR` solves a per-state gradient-moment requirement jointly and shares one window:

```text
combined_z_area_per_m(p)   signed requirement per state       (slab + partition)
_limiting_index(...)       enumerate states, find the worst
one winder duration        shared, so TE does not vary with kz
_z_winder(p)               realise this state inside that window
```

That is the design document's loop, restricted to **`m0` on one axis**. What is missing is not the
architecture; it is that a requirement is currently a scalar area rather than a value at an order
and a semantic instant, and that the realisation is one lobe scaled per state.

**The scaling is the part that does not generalise.** Scaling one waveform by `s` scales every
moment by `s`, so a single scalable lobe can only reach targets on a ray through the origin. An
`(m0, m1)` pair needs two degrees of freedom.

## 2. Two lobes in a fixed window reach any `(m0, m1)`

With two lobes of fixed duration in a fixed window, areas `A1, A2` free, the map to
`(m0, m1)` at a semantic instant is linear and invertible — the centroids differ, so the
2x2 system is non-singular. Measured, in one window of 2 x 700 us on `y`, echo at 3.0 ms:

```text
requirement                     m0 error     m1 error    peak
flow-comp, ky = 0                0.0e+00      0.0e+00      0 %
flow-comp, ky = +16              1.4e-14      2.8e-17     66 %
flow-comp, ky = -32              8.5e-14      1.1e-16     81 %
venc 1.5 m/s, +, ky = 0          2.8e-14      5.6e-17     70 %
venc 1.5 m/s, -, ky = 0          2.8e-14      5.6e-17     70 %
venc 1.5 m/s, + AND ky = +16     7.1e-14      2.8e-17     97 %
venc 0.5 m/s, +, ky = 0          2.3e-13      1.1e-16     93 %
```

Every row is the **same solve** with a different target vector, inside the **same** window, sized
once for the limiting state. Phase encoding, flow compensation and velocity encoding are not three
mechanisms here; they are three values of `(m0, m1)`.

## 3. What that says about the three questions

**WHAT must be true** is a small vector per axis per state: `(m0, m1)` at a named instant. The C3
scan already said flow compensation and velocity encoding "share one requirement form, different
target" — this is that sentence made arithmetic.

**WHO owns it** is the repetition kernel, because only it knows the fixed contributions, the free
window and the state. Neither leaf can see both.

**HOW** is a linear solve of `order + 1` lobes in a fixed window, plus a monotone search for the
window length. **No optimizer is required by any case examined.**

## 4. What a kernel would have to expose — and what it would not

Everything the probe needed is a number the kernel already computes:

```text
per axis, per state:
    the semantic instant                     time_to_echo(), already public
    the fixed moments up to it                integrated by the kernel from its own events
    the free window (start, max duration)     the kernel's own timing solve
    the base requirement for this state       ky/kz encode area, already computed
    hardware limits                           opts
```

**No event handles cross the boundary.** The kernel hands over moments and a window; the designer
hands back areas; the kernel places the lobes. An augmentation never sees a gradient object, never
inspects a built `LogicBlock`, and never names a kernel class.

## 5. PR #39's construction, restated

`CartesianLine(null_moment_order=1)` is exactly the `order = 1` case of section 2, with the target
fixed to `(-m0_ro, -m1_ro)` and the window fixed to the prephaser. Its solve, its monotone duration
search and its tests are all reusable as the realisation primitive.

What is **not** general is its public surface: a per-leaf `null_moment_order` is one axis, one
target, one augmentation. Repeating that shape for velocity encoding, and again for the next
augmentation, on each gradient-producing leaf, is the `N x M` growth the design document exists to
avoid.

## 6. The origin correction — and what it invalidated

**A moment target needs an origin as well as an endpoint.** Under a shift of origin
`m1' = m1 - dt * m0`, so "m1 = 0 at the echo" is only unambiguous when `m0 = 0` there.

| | |
|---|---|
| readout axis, PR #39 | `m0 = 0` at the echo, so `m1` is **origin-free**. Measured about the echo, the block start and -5 ms: `1.03e-01` for all three uncompensated, `~1e-16` for all three compensated. **PR #39 is unaffected.** |
| phase-encode / partition axes | `m0 != 0` by design, so the origin decides the number |

The size of the ambiguity, at `ky = +16` with a 2.5 ms excitation-to-echo interval, is
`dt * m0 = 0.182 s/m` — larger than a whole VENC target of `0.167 s/m` at 1.5 m/s. So **the
ky/kz rows of section 2 were measuring something that had no physical name**, and are withdrawn.
Re-solved with `origin = excitation centre`, `endpoint = echo`, the same two-lobe basis hits every
one of them exactly, and the window it needs is *different* — 370 us rather than 440 at
`ky = +16` — so the origin changes the design, not only the bookkeeping.

The model therefore carries both instants:

```text
MomentTarget(axis, order, origin, endpoint, value)
```

`order = 0` does not read the origin; that is not a reason to leave it out of the semantics.

## 7. Composition: the two augmentations constrain orthogonal components

The C3 record is explicit that velocity encoding's public quantity is the **change** between
states -- "the module emits the two toggles whose first moments differ by delta_m1" -- and that
the `+-delta_m1/2` split comes from the standalone realisation making the toggles exact waveform
negations. So:

```text
VelocityEncode   m1(+) - m1(-)        = delta_m1        the DIFFERENCE
FlowComp         (m1(+) + m1(-)) / 2  = 0               the COMMON MODE
```

Composed, those give `m1(state) = state * delta_m1 / 2`, so the symmetric pair is **derived**
rather than assumed. And because the two augmentations write **different components** of the
state-indexed `m1` vector, they cannot issue contradictory absolute targets and no precedence
rule is needed. The rule that *is* needed is narrow: each `(axis, order, component)` may be
claimed once, and a second claim on the same component is a refusal.

Measured on one shared window of 2 x 630 us, over ten `(ky, polarity)` states:

```text
difference   0.333333 against 0.333333 wanted, error <= 4.4e-16
common mode  <= 1.4e-16
m0           exactly the encode area, every line
```

The limiting state is `(ky = -32, polarity = +1)` -- neither an extreme of `ky` nor a fixed
polarity, which is `GRE3DTR`'s own "the limiting partition is a result of the coupling, never a
property of the index" showing up again in a bigger state space.

## 8. The real kernels already own the window decision

```text
GRE2DTR   winder_s = ceil_raster(max(ro.prephaser_duration_s,
                                     pe.min_duration_s,
                                     exc.rephaser_duration_s))
          then handed to each leaf as prephaser_duration_s / duration_s
```

That is one shared window sized from the limiting demand, across three axes, already. `GRE3DTR`
does the same across partitions. **Neither needs a new mechanism; both need the requirement they
size against to carry an order and an origin.**

Measured on the reference protocols:

| | GRE2DTR | GRE3DTR |
|---|---|---|
| window today | 520 us | 320 us |
| origin (`exc.time_to_center()`) | 1.600 ms | 0.200 ms |
| echo | 6.218 ms | 3.198 ms |
| window for `m1` too | 1440 us on `y` (venc 1.5 + flow comp) | 580 us on `z` (flow comp) |
| TE cost | +0.920 ms | +0.260 ms |

The shared-window-from-the-limiting-state property survives in both.

## 9. Slice selection needs no new leaf API

With the origin at the excitation centre, the excitation's own contribution is just another fixed
moment. Measured on `GRE2DTR`, from the excitation centre to the echo:

```text
z   m0  +0.0000 1/m      the slice rephaser already nulls it
z   m1  -0.425956 s/m    the fine scan's `ms`
```

`-0.426 s/m` is 2.5x a VENC-1.5 target, so it is not negligible. It is computed by the **kernel**
from events it already places, which is the same way the readout's fixed moments are obtained.
No leaf publishes a moment; the fine scan's open question about what a leaf must expose is
answered by not needing anything new.

## 10. The basis formulation

What the probes actually did *is* the basis formulation:

```text
basis column j        a unit shape placed in the window
M[order][j]           that column's order-th moment about the origin
solve  M a = target - fixed
```

Two fixed-duration lobes are two columns and a 2x2 matrix whose coefficients are the areas. It
subsumes what exists: `GRE3DTR`'s single reference lobe scaled per partition is the **1x1** case,
with `scale_grad` as the coefficient. It also places the slice-selection term correctly -- as part
of `fixed`, not as a column. Multi-echo is where it stops being square: more rows than columns,
which needs more lobes rather than a different formulation.

## 11. Open, and not answered by this spike

- **`wave-gre-flow-comp` has not been read**, because it is not in this workspace and is not a
  registered corpus in `sources.yaml`. Nothing here is evidence about it, and the question it was
  reserved for -- whether the boundary describes a real hard-coded case without sequence-specific
  logic in the shared designer -- is **unanswered**. It needs a checkout.
- Multi-echo targets, which are more rows than a two-lobe basis has columns.
- Whether the `y` window can be lengthened inside `GRE2DTR` without disturbing its TE/TR refusals
  in some protocol; the cost was measured, the refusal paths were not exercised.
- Whether `exc.time_to_center()` is the right origin for a non-selective or adiabatic excitation,
  where the effective centre is not the pulse centre.
