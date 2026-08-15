# Parked examples

**These notebooks do not run against the current `seqcraft`.** They are kept, unmodified, because
they are the specification for what comes next — not because they are expected to work today.

## Why they are here

Both folders were built on `seqcraft.modules`: `SincExcitation`, `MonopolarDiffusion`, `SpiralVDS`,
`EPIReadout` and the rest. That library was deleted rather than migrated
([ADR-003](../../docs/adr/003-scanner-and-module-reform.md)) — it was written against a base class that no longer
exists, a `System` that has been deleted, and a declared `duration` that has been removed, and a
compatibility shim would have kept asserting all three.

## Why they were not deleted with it

They are the **acceptance test for whether the next primitive set is right**. A module library
chosen in the abstract acquires classes nobody needs and misses the ones everybody writes by hand;
one chosen by making these two scans work again does not. So rewriting these is what should *drive*
the new library, not follow it — which means they have to survive until then.

Read them for what a real DTI acquisition needs end to end:

| | What it covers |
|---|---|
| `dti_spiral/` | Single-shot spin-echo **spiral** DTI at 1.88 mm, and the two-echo field map its 67 ms readout cannot do without |
| `dti_epi/` | The **same diffusion encoding** through a ramp-sampled **EPI** train — single-shot and two-shot, partial Fourier 0.75 |

The physics inside the deleted classes they depend on — the b-value solve, the variable-density
spiral trajectory, the EPI ramp-sampling moment integral — was lifted out as plain functions before
the deletion and is in [`salvage/`](../../salvage/).

## What still works

`examples/lib/` is unaffected. `mr0_bridge.py` and `noncartesian_recon.py` take explicit arguments
rather than reading attributes off module objects, so they carry over to whatever builds the
sequences next.
