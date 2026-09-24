# Repetition-level joint design — current spike findings

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-24, so that the pull request carries the direction behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../../README.md) for which one to edit.

**Status:** current state of `spike/repetition-joint-design`. Internal and provisional; no public
API is proposed here.
**Date:** 2026-09-24
**Governed by:** [`2026-09-24_repetition_physical_design_architecture.md`](2026-09-24_repetition_physical_design_architecture.md)

Every number below is measured by building real events and integrating the emitted waveform.

---

## 1. What exists

`src/seqcraft/modules/_joint.py` — typed stages composed as functions, no solver hierarchy:

```text
claims + states      resolve_claims  -> one absolute target per state
targets + fixed      JointProblem    -> what must be true, and what is already there
a candidate          Schedule        -> semantic origin, endpoint, window start and length
a family             realise_*       -> events, or None, with utilisation and what limited it
the orchestration    attempt/design  -> the first schedule every claimed axis can meet
diagnostics          utilisation     -> how far over a shape is when it does not fit
```

**`GRE2DTR` and `GRE3DTR` are rewired.** Each supplies a semantic origin, a window, its own
per-state base requirement and the fixed contribution it already places; neither names an
augmentation. With no claims the local path is untouched — same winder, same TE, same events.

27 tests in `tests/modules/test_joint_moment_design.py`; 1474 in the suite.

## 2. The model

```text
MomentTarget(axis, order, origin, endpoint, value)
```

Both instants are required: under a shift of origin `m1' = m1 - dt m0`, so a first moment is not a
number until the origin is named. `order = 0` does not read the origin, which is why omitting it
survived the first pass. The origin is **supplied by the repetition**, never read off a particular
class.

Scope of the current backend: **orders 0 and 1**. For a symmetric lobe `m0` is its area and `m1`
its area times its centroid, exactly; that shortcut does not extend to `m2`, where the lobe's own
internal shape contributes. A second-order family must compute its own exact response.

`target = None` means **unconstrained**; `target = 0.0` means **constrained to zero**. They are
different requirements.

## 3. Composition: augmentations resolve state, the realisation sees absolutes

```text
VelocityEncode   m1(+) - m1(-)       = delta_m1      the difference
FlowComp         (m1(+) + m1(-)) / 2 = 0             the common mode
```

Different components of the same state-indexed vector, so they compose without a precedence rule
and the symmetric pair `+-delta/2` is **derived** rather than assumed. A second claim on the same
`(axis, order, component)` refuses. `VelocityEncode.delta_m1_for(venc)` is the single source of
the relation, available without building a waveform.

The realisation layer receives one absolute target for one state and never learns where it came
from.

## 4. One schedule, several axes

A candidate is a window **and** the fill an explicit TE would need at that window. Every claimed
axis is solved against the schedule that results, and a candidate is accepted only when all of
them fit. One `JointProblem` per axis — never one multidimensional solver.

Measured, `GRE3DTR` with velocity encoding on `z` and flow compensation on `y`, selective slab:

```text
one shared schedule, two designs      y: two-lobe        z: two-lobe/0.513
whole-repetition residuals            m0 2.4e-13         m1 2.3e-15
limiting states                       y (0, +1)          z (2, +1)
```

## 5. Timing

```text
kernel                winder us   min TE ms   min TR ms   family
GRE2DTR plain               520       4.618       9.650   -
GRE2DTR joint              1390       5.488      10.950   two-lobe/0.532
GRE3DTR plain               320       2.998       6.630   -
GRE3DTR joint              1020       3.698       7.390   two-lobe/0.520
GRE3DTR slab joint         1260       5.338      10.670   two-lobe/0.524
```

`te_s=None` recomputes the minimum after the enabled physics. An explicit time below the new
minimum refuses, naming `te_s`, and is never silently rounded up.

## 6. Realisation families

```text
single lobe                m0 only.  Declines an m1 it cannot produce
adjacent two lobes         2 DOF
base + zero-area bipolar   2 DOF, the decoupled coordinates the stress case uses
raster-aligned split       the first family with real shape freedom
```

**Representation and shape family are different things.** The adjacent and decoupled bases span
the *same* `(m0, m1)` set and report identical utilisation — two coordinate systems for one
two-dimensional space, so the decoupled one is better conditioned, not more capable. The split
search is 8 to 13 per cent shorter at the same target.

## 7. The faithful stress harness

Built from the admitted stress case at `0e1ec51`: **fixed echo spacing**, the inter-echo
first-moment recurrence `m1 = -m0 * ESP`, and the mixed design profiles the reference actually
declares.

```text
physical                        80 mT/m    200 T/m/s
sequence envelope  x0.90/x0.70      72         140
lowPNS (most)      x0.90/x0.41      72          82
lowPNS2 (FC)       x0.60/x0.35      48          70
```

At `ESP = 3.0 ms`, leaving a 1400 us inter-echo window, worst-case `ky` needing `m1 = -0.436 s/m`:

```text
design profile              u = max(G_req/G_lim, S_req/S_lim)     fits
sequence envelope                    0.85                          yes
lowPNS (most parts)                  1.46                          no
lowPNS2 (FC module)                  1.71                          no
```

**The same waveform, in the same window, fits the sequence envelope with 15 per cent headroom and
fails its assigned design profile by 71 per cent.**

Sweeping the FC profile's slew at that fixed window gives `u` = 0.85, 1.00, 1.20, 1.46, 1.71, 1.99,
2.39, 2.99, 3.98, 5.98 — **smooth, with the crossing between 0.60 and 0.50 of physical slew.** No
jump anywhere.

And AUTO timing buys the difference cheaply. The shortest window that fits, per profile, and the
echo spacing it implies:

```text
slew fraction   window us   implied ESP   vs fixed ESP
0.70                 1330      2.930 ms         0.98x
0.41                 1590      3.190 ms         1.06x
0.35                 1680      3.280 ms         1.09x
0.15                 2220      3.820 ms         1.27x
0.10                 2550      4.150 ms         1.38x
```

**At the FC module's own profile a 9 per cent longer echo spacing restores feasibility.**

### Which of the four mechanisms

```text
1  genuine physical frontier            NOT the explanation -- it fits the envelope
2  waveform-family artificial cliff     NOT observed -- u degrades smoothly everywhere
3  fixed-schedule artifact              YES, and mild: +9 % ESP removes it
4  profile-assignment artifact          YES, and dominant: 0.85 under the envelope, 1.71 under
                                        the FC profile, for the same waveform
```

**The large-duration pathology is not reproduced.** No GrOpt comparison was run, because no cliff
was found to compare at.

### What this harness does not model

The recurrence target and the mixed profiles are faithful; the wave gradients' own fixed
contributions and the real readout geometry are not — `fixed` is zero here. So it is faithful in
structure rather than in every term, and a full reproduction would need the wave waveform itself.

## 8. Corrections discovered during the spike

Each was found by measuring the emitted waveform, and none would have shown in a per-lobe check.

| | what was wrong | how it showed |
|---|---|---|
| **origin omitted** | `(m0, m1)` at an instant, with no origin | at `ky = +16` the ambiguity is `dt m0 = 0.182 s/m`, larger than a whole VENC target of 0.167 |
| **slab double-counted** | `combined_z_area_per_m` already contains the slab rephasing, and was also used as the target | m0 residual 17 1/m on a selective slab |
| **wrong excitation integrated** | `fixed` used `exc()` while the repetition emits `exc(rephase=False)` | m1 residual 0.026 s/m |
| **zero target vs unconstrained** | the single-lobe family asked *whether* an m1 was requested, not what it would emit | flow compensation reported success while emitting `m1 = -0.066` |
| **explicit-TE translation** | designed at the minimum schedule, then fill shifted the waveform | m1 wrong by **0.264 s/m** against a 0.167 target at TE + 2 ms |
| **bad cliff diagnostic** | utilisation designed against a permissive system, so `make_trapezoid` chose near-rectangular shapes | slew inflated 5x, manufacturing a cliff that is not there |

A seventh is a hazard rather than a failure: solving one axis at its own minimum and letting
another axis widen the winder afterwards is the same staleness as the explicit-TE bug. **I could
not construct a GRE protocol where it bit** — the joint axis wanted the longest window in every
case tried — so the invariant is pinned directly instead: the design's schedule is the schedule
the kernel emits.

## 9. Open

- The public augmentation API. `joint_claims` and `encoding_states` are the spike's internal
  representation and must not ship.
- `VelocityEncode` designs its standalone bipolar eagerly; a joint realisation may be feasible
  where that standalone one is not.
- A whole-waveform safety evaluator: seam and ownership settled, nothing implemented.
- Pathway-aware moments, deferred: pathway-aware `m0` first, `m1` later, neither now.
