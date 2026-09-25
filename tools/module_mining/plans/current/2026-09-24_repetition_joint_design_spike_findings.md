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

## 7. Stress evidence

Two separate things, and the second is the gate.

### 7a. The adversarial stress matrix

Not a handful of mild derating points: a broad sweep over the current architecture, searching
deliberately for a region where a small input change produces a disproportionate jump in the
minimum feasible window. Behaviour like `1.8 ms, 1.9 ms, 2.0 ms -> 18 ms` is what it was looking
for.

```text
gradient limit        72  48  32  20  12  8 mT/m        (and 3 to 72 in the venc pass)
slew limit           140 100  70  45  28 18 11 T/m/s    (and 6 to 200 in the venc pass)
m0 target              0  40  145  400 1/m
m1 target              0  0.05  0.167  0.44  1.2 s/m
fixed contribution    (0,0)  (0,-0.43)  (120,-0.03)  (-200,0.25)   -- slab and wave-like tails
origin-to-window lead  0  1.5  5.0 ms                   -- what an explicit TE fill inserts
window duration        every gradient raster step
```

```text
10 080   coarse matrix cells, 9 897 feasible under 6 ms
 1 920   whole-kernel GRE2DTR builds over venc, matrix, bandwidth, echo count
   972   multi-echo builds over 1..8 echoes, polarity, AUTO and explicit ESP
   648   standalone VelocityEncode systems, 3-72 mT/m x 6-200 T/m/s x venc 5.0-0.01
 3 072   cells for the hole census, 100 raster steps above each first feasible window
37 440   points comparing the continuous diagnostic against the emitter
   243   two-state joint designs with a schedule-dependent wave-like fixed contribution
```

plus fine 1-D refinements at 2 % steps in the hardware limits, and 20 us steps in the lead.

**Minimum feasible window is monotone in the hardware limits, and the steps are proportionate.**
Refined locally at the five largest coarse transitions:

```text
refined step                          steepest adjacent jump
slew, 1 T/m/s steps, 30 -> 16              x1.140
max_grad, 0.5 mT/m steps, 13 -> 7          x1.175
fixed m0, 20 1/m steps, +130 -> -210       x1.531
m0 target, 15 1/m steps, 130 -> 420        x1.745
origin lead, 20 us steps, 1.40 -> 1.80 ms  x1.154
```

The three above `x1.1` are all the same thing and it is not a cliff: they sit where the **net**
requirement passes near zero, so the ratio is large because the denominator is small. At
`m0 = 145` against a fixed `+120` the net is 25 1/m and the window is 470 us; one step later the
net is 40 and the window is 820. The window tracks the net requirement, which is what it should
do.

End to end through `GRE2DTR`, the same statement in TE:

```text
adjacent step        min TE ratio  median   p99     max
slew limit                          1.059   1.251   1.258
gradient limit                      1.049   1.192   1.214
venc                                1.076   1.310   1.341
matrix 128 -> 256                   1.134   1.301   1.303
```

**No jump above `x1.35` in minimum TE anywhere in the kernel sweep.**

### 7b. What the search did find

**A turning point in the lead, not a cliff.** Pushing the semantic origin further before the
window helps until it does not. At `G = 8 mT/m`, `m0 = 400`, `m1 = 0.44`, fixed `(0, -0.43)`:

```text
lead ms    1.40  1.46  1.50  1.54  1.56  1.58  1.60  1.70  1.80
window us  1390  1340  1300  1260  1360  1570  1690  2070  2310
```

Longer lead gives each unit of area more `m1` leverage, until the `m0` target pins the area and
its own `m1 = area x centroid` overshoots, so a growing opposing correction is needed. The curve
is a smooth V with its minimum at about 1.54 ms; at 20 us resolution the steepest step is
`x1.154`. **"A longer TE fill always helps" is false**, and a schedule search that assumes it
will walk the wrong way past the optimum.

**Holes in the feasible set.** Over 3 072 cells, scanning 100 raster steps above each cell's
first feasible window:

```text
cells with at least one hole      48  (1.6 %)
hole width                        median 3, p95 11, max 11 raster steps
maximum OBSERVED width            110 us, in this sampled census
```

A window can refuse where a shorter one succeeded. The two-lobe split must land on the gradient
raster, and for some window lengths no raster-aligned split solves the 2x2 system inside the
limits while both neighbours do. This is mechanism 2, an artificial realisation-family effect.

Two things this is **not**. It is not a bound: 110 us is the widest hole seen over the grid above,
and a target outside that grid could do worse. And it is not something a numerical backend would
necessarily meet as well -- the gradient raster is fundamental, but *this particular* missing
two-lobe split is not. A richer raster-aligned family, with a third lobe or an asymmetric one,
could plausibly fill these holes without leaving the raster.

The kernel walks windows upward and takes the first that fits, so the observed cost is up to
~110 us of TE.

**The continuous diagnostic was measuring the wrong lattice.** `utilisation` sampled 48 fixed
fractions of the window while `realise_split_search` emits on the gradient raster. `u(T)` is not
smooth -- where `m1 / m0` equals a lobe centre the two-lobe solve degenerates to one comfortable
lobe, a narrow minimum -- and a coarse fixed grid steps over those minima:

```text
window us   u(48 fractions)   u(raster lattice)
   710            0.919             0.694
   720            1.171             0.653
   730            0.805             0.636
   800            0.953             0.525
```

Overstated by up to **x1.82**, reporting `u = 1.171` where the lattice reaches 0.653. Against the
emitter over 37 440 points, `u <= 1` agreed 99.82 % of the time before and 99.93 % after; the
aggregate rate understates it, because what matters is the **size** of the error away from the
boundary, and that is what manufactures a cliff. Fixed, and pinned by two tests.

This is the second time a diagnostic rather than the design produced a false cliff. A validator is
only independent of what it saw: the metric and the emitter must search the same lattice.

**On the tested GRE phase-encode axis, adding readout echoes does not lengthen the compensation
window.** The narrow statement, because it is a statement about one axis of one kernel family:

```text
for the GRE phase-encode axis, with no additional gradient on that axis after the
winder, and with M1 defined about the fixed excitation origin, adding readout echoes
does not increase the required compensation window
```

972 kernel builds over 1 to 8 echoes, bipolar and monopolar, AUTO and two explicit echo spacings,
three gradient and four slew limits, flow compensation and two vencs: **all 972 built, none
refused**, and the winder ratio across every echo-count step was exactly `1.000` over 756 pairs.

**This does not generalise to "multi-echo costs nothing."** It holds because of a condition that
an axis with inter-echo gradients does not satisfy: the phase encode is silent after the winder,
and `m1` about a **fixed origin** stops accruing when the gradient stops, so once it is nulled from
the excitation it stays nulled at every echo. An axis that plays between echoes — a readout, or a
wave gradient — still needs the checkpoint/propagation model, and nothing here tests that.

```text
m1 on y about the excitation centre      echo 1   0.000000   echo 2   0.000000   echo 3   0.000000
m1 on y about each echo                 -0.739636          -1.171636          -1.603636
                                                   drift    -0.432000  -0.432000  per echo
                                        -m0 * ESP = -163.6364 * 2.640 ms = -0.432000
```

The recurrence `m1 += -m0 * ESP` is real and reproduces to six digits — **in the echo-referenced
frame**. It is not a correction the phase-encode axis owes, because the velocity-dependent phase a
spin accumulates is `m1` about the **excitation**, which is the one already nulled. This is the
concrete reason a moment target carries an origin as well as an endpoint: the same waveform is
compensated in one frame and off by 0.74 s/m in the other, and only naming the origin distinguishes
them.

### 7c. The narrow conclusion, and how the search was done

**Extensive stress testing did not reveal a large artificial feasibility cliff.** It did find a
family artefact whose widest observed instance was 110 us, a smooth turning point in the lead, and
a diagnostic that overstated utilisation by up to x1.82 until it was corrected.

This is not a claim that cliffs are impossible. The families here are low-dimensional by design,
and a richer requirement -- a second-order moment, a third coupled axis, a fixed contribution with
structure the two-lobe solve cannot oppose -- could produce one. The claim is about **this sampled
and refined parameter space**.

The method matters to how far the claim reaches. The coarse matrix finds each cell's minimum
window by bisecting to a feasible window and then walking down until 60 consecutive windows
refuse. Bisection alone would be unsound, because this same search establishes that the feasible
set has holes; the walk-down bounds the error at 60 raster steps, which is five times the widest
hole observed. The fine 1-D refinements scan every raster step and assume nothing.

So the numbers above are **observed maxima over a sampled grid**, not bounds over all targets. The
harness is committed and re-runnable:

```sh
python tools/module_mining/stress_repetition_design.py --quick
python tools/module_mining/stress_repetition_design.py --full --json evidence.json
```

It does not run in CI. The conclusions are held in place by the regression tests in
`tests/modules/test_joint_moment_design.py`, not by the harness.

### 7d. Which of the four mechanisms

```text
1  genuine physical frontier         observed and correct -- the window tracks the net
                                     requirement, proportionately, in every refinement
2  waveform-family artificial cliff  PRESENT and small: raster-aligned split holes, 1.6 % of
                                     cells, widest OBSERVED 110 us.  A sampled maximum, and
                                     a property of this family rather than of the raster
3  fixed-schedule artifact           present and mild: AUTO timing removes it for +9 % ESP
4  profile-assignment artifact       dominant in the historical case: 0.85 under the envelope
                                     and 1.71 under the FC profile, for the same waveform
```

Mechanism 4 is a policy question, not a physics one, which is why it belongs to the
envelope/profile split rather than to the realisation layer.

### 7e. The historical fixed-ESP harness

Kept because it is where mechanism 4 was measured, **not as a gate** -- reproducing
`wave-gre-flow-comp` is not a success criterion. Built from the admitted stress case at `0e1ec51`:
fixed echo spacing, the inter-echo recurrence `m1 = -m0 * ESP`, and the mixed design profiles that
implementation declares.

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
0.70                 1340      2.940 ms         0.98x
0.41                 1600      3.200 ms         1.07x
0.35                 1680      3.280 ms         1.09x
0.25                 1880      3.480 ms         1.16x
0.15                 2220      3.820 ms         1.27x
0.10                 2550      4.150 ms         1.38x
```

**At the FC module's own profile a 9 per cent longer echo spacing restores feasibility.**

These numbers are re-derived against the corrected diagnostic; two windows moved by one raster
step and nothing else changed, because this target has `m0 = 0` and so never reaches the
degenerate branch the old grid was missing.

The recurrence target and the mixed profiles are faithful; the wave gradients' own fixed
contributions and the real readout geometry are not — `fixed` is zero here. So it is faithful in
structure rather than in every term.

### 7f. No GrOpt comparison was run

Per the governing decision, a numerical-backend comparison is warranted only if the **current**
stress matrix exposes a suspicious gap — not because the historical sequence once showed one. It
did not, so none was run.

The holes are the one candidate, and they are not being dismissed on the grounds that an optimiser
could not help. **The reason is that the observed artefact is small enough not to justify a
numerical-backend investigation**, not that an optimiser has been shown incapable of improving it.
A richer family, numerical or analytic, might well close them; nobody has tried, and this record
should not be read as saying otherwise.

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
| **diagnostic on the wrong lattice** | utilisation sampled 48 fixed fractions of the window; the emitter searches the gradient raster | `u` overstated by up to **x1.82** -- 1.171 reported where the lattice reaches 0.653 |

A seventh is a hazard rather than a failure: solving one axis at its own minimum and letting
another axis widen the winder afterwards is the same staleness as the explicit-TE bug. **I could
not construct a GRE protocol where it bit** — the joint axis wanted the longest window in every
case tried — so the invariant is pinned directly instead: the design's schedule is the schedule
the kernel emits.

## 9. Open

- The public augmentation API. `joint_claims` and `encoding_states` are the spike's internal
  representation and must not ship.
- `VelocityEncode` designs its standalone bipolar eagerly. That is a **cost, not a capability
  gap**: over 648 extreme systems (3-72 mT/m x 6-200 T/m/s x venc 5.0-0.01 m/s) the standalone
  pair never refused -- its closed-form solve grows the lobe until it fits -- so there is no case
  in which the joint path reaches a venc the standalone one cannot. An earlier record here said
  the opposite; it was a guess, and the search for a counterexample found none.
- A whole-waveform safety evaluator: seam and ownership settled, nothing implemented.
- Pathway-aware moments, deferred: pathway-aware `m0` first, `m1` later, neither now.

Closed since the stress pass:

- **Capability and refusal.** Implemented at the repetition adapter rather than deferred to the
  public API: each kernel declares the axes it materialises and refuses the rest before designing
  anything. `GRE2DTR` owns `('y',)`; `GRE3DTR` owns `('y', 'z')`.
- **Search exhaustion.** Running out of candidate windows now reports a design search limit rather
  than infeasibility.
- **Reproducibility.** `tools/module_mining/stress_repetition_design.py`.
