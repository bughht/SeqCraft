# `salvage/` — physics lifted out of the deleted module library

**This directory is not part of the `seqcraft` package.** It is not installed, not imported, not
covered by CI, and nothing in `src/` refers to it. It exists so that a small amount of genuinely
hard-won physics survived the deletion of `src/seqcraft/modules/`.

## Why it exists

The old module library was 5 762 lines across 27 classes, written against a base class that no
longer exists, a `System` that has been deleted, and a declared `duration` that has been removed.
It was deleted rather than migrated ([ADR-003](../docs/adr/003-scanner-and-module-reform.md)), because a shim
layer would have kept asserting all three.

But the classes were scaffolding around a few pieces of real work — closed-form integrals and
solvers that took measurement and argument to get right, and that are wrong in subtle,
image-degrading ways when re-derived carelessly. Those are lifted here **as plain functions**,
with no scanner object, no base class and no `Opts` in their signatures. They take numbers and
return numbers.

## What is here

| File | Lifted from | What it is |
|---|---|---|
| `bvalue.py` | `modules/encoding/diffusion.py` | Exact b-value of a trapezoidal monopolar or bipolar pair, and the two solvers (shortest lobe / required amplitude) |
| `epi_moment.py` | `modules/readout/epi.py` | The unit-amplitude trapezoid moment integral and its closed-form inverse — what makes ramp sampling and `time_to_echo` agree |
| `spiral_vds.py` | `modules/readout/spiral.py` | Variable-density spiral trajectory design, integrated under slew and amplitude limits |
| `directions.py` | `modules/encoding/diffusion.py` | Electrostatic-repulsion DTI direction tables and their condition number |
| `ordering.py` | `seqcraft/ordering.py` | k-space ordering and RF-spoiling tables — moved whole, not lifted |

`ordering.py` is here for a different reason from the rest. It was not buried inside a class and it
still worked; it was removed because four of its six functions had **never had a caller**, and
because [ADR-003](../docs/adr/003-scanner-and-module-reform.md) assigns ordering tables to the opinionated
library — they are sequence-programming choices, not physics. The two the tests did use,
`interleaved_slice_order` and `rf_spoil_phase`, are three lines each and are now written out in
`tests/integration/conftest.py` where they are needed.

## What is deliberately *not* here

The class scaffolding: `__init__` signatures, unit checks, `build()` methods, `duration`
properties, `params()`. That was the part being deleted, and carrying it forward is what a shim
would have done.

## Where it is going

These belong in the separate opinionated module library described in
[ADR-003](../docs/adr/003-scanner-and-module-reform.md) — diffusion, EPI and spiral are all on the
"not in the core" side of that split. When that repository exists, these move into it and this
directory goes away.

Until then: read them, copy from them, port them. Do not import them from `seqcraft`.
