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
| `geometry_pe.py` | `core/geometry.py` | The partial-Fourier / accelerated / multi-shot phase-encode table, and the residue nudge that keeps k = 0 sampled |
| `geometry.py` | `core/geometry.py` + `core/validate.py` | `Geometry` — FOV, matrix, slice layout, the `[DEFINITIONS]` derived from them, and the unit-plausibility bands |

`ordering.py` is here for a different reason from the rest. It was not buried inside a class and it
still worked; it was removed because four of its six functions had **never had a caller**, and
because [ADR-003](../docs/adr/003-scanner-and-module-reform.md) assigns ordering tables to the opinionated
library — they are sequence-programming choices, not physics. The two the tests did use,
`interleaved_slice_order` and `rf_spoil_phase`, are three lines each and are now written out in
`tests/integration/conftest.py` where they are needed.

`geometry_pe.py` is the same case one layer down. `Geometry` was defended as holding *one*
phase-encode index computation so that `kspace_center_line` and the `LIN` label could not disagree
— but the sharing never happened, because the module library that would have consumed the table was
deleted. The fully-sampled `range(matrix)` the integration recipes actually need is written out at
its call site for the same reason.

`geometry.py` — the rest of the class — followed shortly after, for a different reason again: not
that it had no consumer, but that its **only** consumer was `compile(geometry=)`, which called
`definitions()` on it and merged the eight keys that came back. `compile(definitions=...)` already
does that for any source, so ~450 lines of dataclass and range framework were sitting inside the
package to produce one dict. FOV, matrix and slice order are decisions about the *scan*; the
compiler turns a tree into legal pulseq blocks and is indifferent to why the tree looks the way it
does.

Unlike the rest of this directory, `geometry.py` is **expected back**. A geometry of that shape is
wanted the moment the module library gets its basic infrastructure — a readout and a phase encoder
both need FOV and matrix to size a gradient, and the whole argument for deriving the definitions
from the same fields is that the two cannot then disagree. It is the design to start from, and it is
standalone (no seqcraft imports, plain `ValueError` subclasses) so it can be copied as it stands.

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
