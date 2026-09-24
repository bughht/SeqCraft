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

## 6. Open, and not answered by this spike

- Whether the free window on `y` and `z` in `GRE2DTR`/`GRE3DTR` can be re-tasked without
  disturbing their TE/TR solves. The probe used a standalone window, not those kernels' own.
- Multi-echo targets, which are more constraints than a two-lobe window has freedom for.
- `wave-gre-flow-comp` as a stress case; not yet read.
- Whether the slice-selection case, which needs the excitation's moment from the RF isodelay
  point, fits the same boundary or needs a different one.
