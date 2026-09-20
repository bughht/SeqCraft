# Spiral — the two compiler claims, tested without any Spiral module

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-20, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Method.** Every number below comes from a raw `pp.make_arbitrary_grad` in a bare `LogicBlock`.
No Spiral module, no PR23 code, no cherry-pick. Reproducers in
`tools/module_mining/candidates/spiral/`.

The structural property both claims depend on is not "spiral" — it is **one arbitrary gradient
whose amplitude changes sign inside itself**. A spiral is merely the first module here that has
one.

---

## 1. What did not reproduce

A sign-changing arbitrary gradient, alone in a block, compiles **clean** on current `main`:

```text
n = 4000 knots, net |∫g| = 0.006 1/m, traversed ∫|g| = 39031 1/m
m0 error 0        m1 error 0        compile: OK
```

So the bare event is not enough. `sc.compile` does not compress or resample a single unsplit
event, and the compiled side matches the tree exactly.

**The trigger is a split.** Adding `sc.barrier()` mid-event makes the compiler cut the gradient
across blocks, the pieces are rebuilt, and a real discrepancy appears. That matters for the
architecture assessment too: a spiral's ADC *must* be segmented (sample-count limits), so its
readout necessarily spans blocks. The compiler behaviour and the module's hardest constraint meet
at the same place.

## 2. m1 — **CONFIRMED**, and worse than PR23 described

One arbitrary gradient, split once, swept over knot count and amplitude:

| knots | amp | m1 error (s/m) | current tol `1e-9·net·T` | verdict |
|---|---|---|---|---|
| 1000 | 30 % | 4.02e-08 | 4.06e-09 | **FAIL** |
| 1000 | 90 % | 0 | 1.22e-08 | pass |
| 2000 | 30 % | 4.02e-08 | 8.13e-09 | **FAIL** |
| 2000 | 90 % | 7.11e-15 | 2.44e-08 | pass |
| 4000 | 30 % | 4.02e-08 | 1.63e-08 | **FAIL** |
| 4000 | 90 % | 1.20e-07 | 4.88e-08 | **FAIL** |
| 8000 | 30 % | 4.01e-08 | 3.25e-08 | **FAIL** |
| 8000 | 90 % | 1.20e-07 | 9.76e-08 | **FAIL** |
| 16000 | 30 % | 4.01e-08 | 6.50e-08 | pass |
| 16000 | 90 % | 1.20e-07 | 1.95e-07 | pass |

**Six of ten legal waveforms are refused, and which ones is not predictable from the physics.**
`n=2000` at 90 % passes with an error of 7e-15; `n=4000` at 90 % fails with 1.2e-7. That is
PR23's "for one density and not the next", reproduced on a waveform that has nothing to do with
spirals.

**The separation is enormous and the tolerance sits on the wrong side of it.**

```text
representation floor      ~4e-8 .. 1.2e-7 s/m      (what the round-trip costs)
current tolerance          4e-9 .. 2e-7 s/m        (straddles the floor -- the bug)
one raster of displacement 0.016 .. 0.78 s/m       (what the check exists to catch)
```

Five to six orders of magnitude between floor and signal, and the tolerance is placed inside the
noise.

### Where I disagree with PR23

**The floor is absolute, not relative.** It sits at ~1.2e-7 s/m almost regardless of knot count
and amplitude — it does not scale with the waveform. PR23 fixes this by raising a *relative*
coefficient to `1e-7`, which clears the floor by 40×–10373× for spiral-scale waveforms. But a
relative coefficient against an absolute floor is structurally fragile: a short, weak arbitrary
gradient with `traversed × duration ≈ 1` gets a tolerance of 1e-7 against a 4e-8 floor — **2.5×**,
which is the same class of mistake at a smaller scale.

**A cleaner fix is a floor term**, e.g. `max(1e-9 · scale, 5e-7)`, which is claim-driven rather
than calibrated on one family. I am **not proposing the change here** — §8 of the brief puts that
behind its own review, and this document exists to establish the defect, not to fix it.

**And the two PR23 changes are not independent.** m1's scale derives from `magnitude`, the same
quantity the m0 change redefines. I tested whether the m0 fix alone would rescue m1 at the
existing `1e-9`: it does not. `n=1000` at 30 % still fails (`1e-9 × traversed × T` = 1.62e-8
against a 4.02e-8 error), and the passing margins elsewhere are 2×–104×, still knot-dependent.
**So the coefficient — or better, a floor — is needed on its own merits, and my first hypothesis
that the traversed fix subsumed it was wrong.**

## 3. m0 — **structurally sound, failure NOT reproduced**

The code/comment mismatch is real and checkable by inspection. `verification.py` accumulates
`abs(area)` — `|∫g|` — under a variable named `magnitude` whose own comment says *"the total area
traversed, not the net"*. For a trapezoid the two agree; for one sign-changing event they differ
by up to **1.6e6×** in my reproducers.

**But I could not make it fail.** Across every configuration tried, including the near-zero-net
case where the tolerance collapses to its `max(..., 1.0)` floor:

```text
near-zero net, split:   m0 error 1.8e-12   tolerance 1e-6      pass
large net, split:       m0 error 6.0e-06   tolerance 1.2e-3    pass
```

PR23 reports an m0 error of 1.1e-5 1/m from shape compression against a collapsed 1e-6 tolerance.
I did not observe compression of that size; my split errors are 6e-6 at worst and the tolerance
was never the collapsing one at the same time. Their failure plausibly needs **superposition** —
their own changelog mentions a gradient "laid across the spiral's own axes" forcing a resample
that moves m0 by 0.046 1/m — which my reproducers do not exercise.

**Assessment:** `UNPROVEN_CLAIM` as a defect; `PHYSICAL_REQUIREMENT` as a correctness statement.
The change makes the code do what its comment says and is strictly more permissive only where an
event changes sign. It should not be merged as a *bug fix* on this evidence, because the bug did
not reproduce. It could be merged as a *correctness alignment* with a test asserting
`traversed ≠ |net|` for a sign-changing event — but that is a separate, smaller argument than
PR23 makes, and it should not travel inside a Spiral extraction.

## 4. Conclusion

```text
m1 defect       CONFIRMED independently, no Spiral module involved
                current tolerance straddles the representation floor
                6 of 10 legal waveforms refused, unpredictably
                -> a compiler fix is justified ON ITS OWN, before and separate from any Spiral work
                -> but PR23's relative 1e-7 is not the best form; an absolute floor is
                   claim-driven where a bigger coefficient is calibrated to one family

m0 defect       NOT REPRODUCED
                the comment/code mismatch is real and worth correcting
                the failure it was introduced to fix did not occur in a candidate-free reproducer
                -> do not port as a bug fix; re-examine if superposition reproduces it
```

Both belong in their own PR with their own reproducer tests, **not** in a Spiral extraction. A
compiler fix that arrives with a module is a compiler fix nobody can review independently.
