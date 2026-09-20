# What reversibility costs

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-20, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

The durable family rule is `g = 0` at **both** ends of every arm. `writeSpiral.m` ends its
spiral-out at full gradient and ramps down inside the spoiler — and has no spiral-in, which is the
same fact. This measures the difference so the policy is adopted with a number beside it.

Reproducer: `tools/module_mining/candidates/spiral/run_endpoint_cost.py`. Same path, same scanner,
same solver, same dwell; **only the terminal boundary condition differs.** The traversal is a
standard time-optimal forward–backward sweep over arclength under the gradient and slew limits.

---

## 1. The cost

| protocol | free / ms | at rest / ms | penalty | | samples + |
|---|---|---|---|---|---|
| reference-like 256/4 | 34.125 | 34.227 | 0.102 ms | **0.30 %** | 40 |
| modest 128/4 | 14.357 | 14.469 | 0.112 ms | **0.78 %** | 45 |
| single shot 128/1 | 57.245 | 57.356 | 0.112 ms | **0.20 %** | 45 |
| high res 256/8 | 20.735 | 20.871 | 0.137 ms | **0.66 %** | 55 |
| small FOV 128/4 | 19.618 | 19.755 | 0.137 ms | **0.70 %** | 55 |
| weak gradients 128/4 | 20.310 | 20.438 | 0.128 ms | **0.63 %** | 51 |
| strong gradients 128/4 | 12.434 | 12.530 | 0.097 ms | **0.78 %** | 39 |

**Worst case 0.78 %, on the shortest readout.** Nothing is pathological in any regime tried —
weak gradients, strong gradients, single-shot, high resolution, small FOV.

## 2. The shape of the cost, which is the part that generalises

**The penalty is a fixed ~0.1 ms, not a fraction of the readout.** It ranges 0.097–0.137 ms
across readouts spanning 12 to 57 ms. It is a braking distance: the time to decelerate from the
arm's terminal speed to rest under the slew limit, and it does not care how long the arm was.

So the relative penalty **shrinks as readouts lengthen** — 0.78 % at 12 ms, 0.20 % at 57 ms — and
the worst case is bounded by the shortest readout anyone would design.

Split as the brief asks:

```text
path cost         ZERO, by construction.  k-space extent is identical between policies, and the
                  path is the same object -- the policy is a traversal boundary condition and
                  touches no geometry.

traversal cost    0.097 .. 0.137 ms, fixed.  This is the whole cost.

ADC consequence   39 .. 55 extra samples, which is exactly the traversal cost divided by the
                  dwell.  Nothing independent.
```

## 3. Two things the study found that were not asked for

**The peak gradient does not rise, and sometimes falls.** Unchanged in four of seven protocols,
**lower** in three. The rest policy brakes where a spiral is fastest — the outer edge — so it can
only reduce the peak, never raise it. Ending at rest is not a hardware cost.

**The entire penalty lands in the outer tenth of k-space.**

| protocol | outer samples, free | outer samples, at rest |
|---|---|---|
| reference-like 256/4 | 2243 | 2283 |
| modest 128/4 | 837 | 881 |
| single shot 128/1 | 3345 | 3389 |
| high res 256/8 | 1305 | 1360 |

Every extra sample is spent at high `|k|`, because that is where the braking happens. That is not
a loss: the outer ring is where a spiral's SNR per sample is worst, so the extra dwell there is
mildly useful. It also means the policy does **not** perturb the centre of k-space at all, which
is where contrast and low-order phase live.

## 4. Conclusion

```text
cost              0.2 - 0.8 %, fixed in absolute terms, worst on the shortest readout
pathology         none found across seven protocols and two hardware regimes
peak gradient     unchanged or reduced
outer-k           mildly oversampled, where oversampling is cheapest
k-space extent    identical
```

**The endpoint policy is adopted.** Under 1 % of readout duration buys `in`, `in-out` and `out-in`
as members of one family rather than three separate implementations, plus joins that need no
connector. That is a good trade, and it is now a measured one.

**It does not become a reason to weaken the reference comparison.** `out`'s tail still differs
from `writeSpiral.m` deliberately, and `reference_review.md` §3 states where. A measured
justification for diverging is not permission to claim similarity.
