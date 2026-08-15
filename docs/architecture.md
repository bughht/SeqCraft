# Architecture

Two concepts, and no more.

```
   any Python you like                                  seqcraft does this
  ─────────────────────────                          ──────────────────────────
  a function                ─┐
  a class of your own shape  ├─►  LogicBlock  ──►    sc.compile()
  sc.modules.SincExcitation ─┘    lb.add(t, ...)       finds block boundaries
                                                       sums same-axis gradients
                                                       checks the amplifier
                                                       ──► pypulseq.Sequence
```

**`LogicBlock` is the interface.** It is where seqcraft imposes structure, and the only place. What
produces one — a function, a class, a notebook cell, a module from the library — is not seqcraft's
business, and nothing on the path from a tree to a `.seq` inspects the producer's type.

---

## `LogicBlock` — two attributes

```python
class LogicBlock:
    tag: str                 # optional; for plots and error messages
    nodes: list[Node]        # Node(start, item); item is an event or another LogicBlock
```

A tree of children, each with a start time relative to this block's own start. That is the whole
model. `duration` is a read-only property measured from the nodes, so a block cannot advertise a
length that disagrees with what it plays.

**Overlap is legal everywhere.** Two gradients on one axis, a gradient across an RF, anything. Making
it legal *for pulseq* is the compiler's job, and keeping those two concerns apart is what lets a
module say what it means without also knowing where the block boundaries will land.

CRUD is Python's, not seqcraft's: create with `add`, read with `lb.nodes[i]` or iteration, update with
`lb.nodes[i].start += dt`, delete with `del lb.nodes[i]`. One method of our own.

### What it deliberately lacks

| Absent | Where the job went |
|---|---|
| `marks` / anchors | **Whatever produced the block.** A block cannot know where it sits, so it cannot know when its echo occurs — and pinning the answer in would stop the same block being reused. `readout.time_to_echo` is a property of the readout. |
| `labels` field | A pulseq label is an event, so it is a node: `add(t, pp.make_label('LIN', 'SET', k))`. |
| a declared `duration` | Measured. To make a block longer, add a delay event — which is how you lengthen a block in pulseq anyway. |
| `then`, `over`, `chain`, `scaled` | `add` is shorter than any algebra and needs nothing learned. |
| `walk`, `axes`, `moment` | The compiler. Reading a tree is not part of being one. |

---

## `sc.compile` — the scheduler

The only deliberately complicated subsystem in the project: one public façade, explicit internal
stages, and tests at each responsibility boundary. The Phase 0 implementation still concentrates
placement, boundary selection, gradient legalization, emission, and verification in
`core/compiler.py`; the refactor extracts those responsibilities without changing the public API.
Internally this is treated as a deterministic constrained-scheduling algorithm: tree placement,
interval legality, boundary selection, piecewise-linear gradient transformation, and mechanical
emission.
Boundaries come from where RF and ADC events fall and nowhere else under the current policy;
same-axis gradients are summed with a warning; different axes are silent; limits are checked on the
*compiled* waveform because that is the only place a merge's effect is visible. Full detail in
[`compiler.md`](compiler.md).

Returns a `CompiledSequence` holding the `pypulseq.Sequence`, the compile report, and per-block
provenance. There is **no `Sequence` class** in seqcraft — a sequence *is* a logic block, and
compiling it produces the artifact.

Its signature is `compile(root: LogicBlock, system: System, ...)`, and that is the whole of what it
knows about the world upstream of it. A test asserts that `core` never imports `seqcraft.modules`,
because the direction of that arrow is the design.

The compatible implementation boundary is:

```text
src/seqcraft/core/
├── compiler.py              public façade and orchestration
└── _compiler/               private stage contracts and implementations
    ├── model.py
    ├── placement.py
    ├── legalization.py
    ├── emission.py
    └── verification.py
```

The private package is introduced incrementally. It must not expose a second compile path or force
users to import stage types.

---

## Modules — a library, not an interface

`sc.modules` is a **provided library of reusable design components**: excitations, readouts, phase
encodes, spoilers, diffusion encodings. Everything in it is written on `sc.modules.Module`, an
**optional** base:

```python
class Module:
    def __init__(self, system, *, regime='default') -> None
    @property
    def opts(self) -> Opts       # the resolved pypulseq limits for `regime`
    def params(self) -> dict     # for the provenance sidecar and repr
    def submodules(self) -> dict[str, Module]
```

It declares no abstract method. There is no `build()` requirement, no class variable to declare, no
category to register and no hook to override — inheriting it is a way of not rewriting the same
bookkeeping, and nothing more. What it does behind your back is unit validation: `check_units` runs
when a subclass's `__init__` returns, so `fov_mm=0.22` fails at construction with a hint rather than
producing a sequence with a 22 cm error in it. `check_units` takes any object, so a component that
inherits nothing can call it itself.

The library's own modules follow one convention that the base does not enforce: `__init__` designs
once, `build(**args)` is cheap and returns a block, and its arguments select the *variant* — which
diffusion lobe, which line, which slice, what phase — without changing the block's duration. That
convention suits a module with one output. A component with two — `diffusion.pre()` and
`diffusion.post()`, `readout.readout()` and `readout.prephaser()` — should say so in its method
names, and seqcraft has no opinion about it either way. See
[`writing_a_module.md`](writing_a_module.md).

**Why it is not in `core`.** `core` is what gets a logic block to a validated `.seq`. A convenience
base for writing components is not on that path, and while it lived there it read as a requirement:
"three concepts" implied a sequence had to be expressed as modules. It never did.

---

## The rest of `core`

`core` holds what is required to get from a logic block to a legal, validated `.seq` — and nothing
else. That is the whole membership rule, and it is what keeps the layer small enough to read.

Membership in `core` does not require every responsibility to remain in one flat file. Conversely,
moving a file is not evidence that its responsibility improved. Phases 1–6 only organize compiler
internals under `core/_compiler`; they do not relocate unrelated core modules. Phase 7 records a
dependency- and cohesion-based boundary audit, and any broader move requires a separate decision and
compatibility plan. The audit scope and evidence rules are defined in
[`refactor/core_package_boundary_audit.md`](refactor/core_package_boundary_audit.md).

| Module | What it is for |
|---|---|
| `system` | `System` holding **named `Opts` regimes**, so a diffusion encoding can run at full amplitude while the readout is derated. Asserts every regime agrees on the rasters, gamma and B0 — nothing upstream checks that, and a disagreement produces an unplayable file. Also the source of every `Raster` and of `system.convert`. |
| `timing` | `Raster` — the raster as an object, with `ceil / floor / nearest / count / at / holds / require`. Arithmetic in integer ticks, because `1.5e-3 / 1e-5` is `149.99999999999997` and `250 * 1e-7` is not `2.5e-5`; both errors propagate into ADC dwells and off-raster block durations. Nothing here assumes 10 µs: rasters are values the scanner supplies. |
| `units` | One function — `convert(value, from_unit, to_unit, gamma=, f0=)` — over eleven dimensions, in both directions, the shape `pypulseq.convert` uses. Scales are exact `Fraction`s, so `4200 us` is `0.0042 s` and not `0.004200000000000001 s`. |
| `events` | `derive()` — the one sanctioned way to copy a pypulseq event. Also `knots_of` and `pwl_moment` (the exact-gradient primitives), `waveform_of` (a curve to plot, **not** a uniform raster — a `grad` event carries its own sample times), `moment_of`, `content_hash`, `check_limits`. |
| `geometry` | FOV, matrix, slices, and one authoritative phase-encode index computation shared by `kspace_center_line` and the `LIN` label values, so the two cannot disagree. |
| `validate` | Unit plausibility bands inferred from a field's name suffix, plus explicit `require_*` helpers. Its unit names are the ones `convert` knows, asserted by a test — one vocabulary, not two. |
| `report` | `Issue` and `Report`. Immutable; `ok` is true when there are no errors, so warnings inform without failing. |

## What sits *outside* `core`, and why

None of these is on the path from a tree to a file, so none of them is core. They are ordinary
top-level modules, and `sc.ordering`, `sc.plot_block`, `sc.testing` are unchanged as import paths.

| Module | What it is for |
|---|---|
| `modules` | The library of reusable design components, and `Module`, the optional base they share. |
| `ordering` | `interleaved_slice_order`, `centric_order`, `golden_angle`, `rf_spoil_phase` — sequence-programming vocabulary, the closed forms that otherwise get copy-pasted per notebook. Becomes a package when the trajectory work adds a second file beside it. |
| `provenance` | The JSON sidecar: versions, git commit and dirty flag, definitions, sha256. Output tooling; the compiler imports it lazily at `write()` time. Takes a mapping, so a component that reports its own parameters however it likes is not shut out. |
| `display` | **The only module allowed to import matplotlib**, and lazily at that, so `import seqcraft` stays cheap. Every function returns a figure and never calls `show()`. |
| `testing` | Two tiers. `assert_output(make, system)` and the block-level assertions it composes take a callable and a block, so they ask nothing about ancestry; `assert_all(module, **build_args)` adds the checks that only mean something for the `build()` convention. `module_subclasses()` enumerates what inherits the base, which is what the library's own contract suite parametrises over. |

**There is no registry.** A registry earns its place when something must turn a *string* into a
*class* at run time: a YAML front end, plugin entry points, a `--readout=spiral_vds` flag. seqcraft
has none of those, on purpose. Its one consumer was the contract suite's `parametrize`, and
`Module.__subclasses__()` serves that with no decorator anywhere — and cannot be forgotten, because
subclassing *is* the registration. A `@register()` that a new module omitted silently lost the whole
contract suite, which is the failure a registry was supposed to prevent.

`module_subclasses()` is scoped accordingly: it answers *what inherits the base*, which is the right
question for parametrising the library's contract suite and the wrong one for deciding what seqcraft
accepts. Nothing outside the tests consults it.

---

## What is not in the package

`examples/lib/` holds the MRzero simulation bridge and the sigpy off-resonance reconstruction. They
are **not** part of seqcraft, on purpose: the package builds sequences, and simulating or
reconstructing them are downstream jobs with heavy dependencies of their own. Folding them in would
make every user pay for them.

**There are no recipes either.** A recipe is somebody else's sequence choices in library code: to
change your scan you would edit a package, and the package would accumulate everyone's variants. The
notebooks assemble sequences from modules instead. With module properties the timing is four lines --

```python
first_half  = (exc.duration - exc.isodelay) + delta + refoc.isodelay
second_half = (refoc.duration - refoc.isodelay) + delta + readout.time_to_echo
te = 2 * max(first_half, second_half)
```

-- so a `solve_timing` function was never earning its place.

Also deliberately absent: any solver for gap timing (`t0 + te - ro.time_to_echo` is arithmetic you
can read), a YAML or GUI front end, and a `.seq` importer.

---

## Layout

The current public layout remains:

```
src/seqcraft/
  core/          logic  compiler  system  geometry  events
                 timing  units  validate  errors  report
  modules/       base.py (the optional Module)
                 rf/  readout/ (cartesian, epi, spiral)  encoding/  prep/  control/
  ordering.py    view orders, golden angle, RF-spoil phase
  provenance.py  the JSON sidecar
  display.py     the only matplotlib importer
  testing.py     assertions you can point at any component of your own

tests/        core/  logic/  compiler/  modules/  integration/
examples/     01_getting_started
              dti_spiral/  dti_epi/   one scan each: 01_build, 02_simulate_and_reconstruct, seq/
              lib/                    sim + recon helpers, not the package
docs/         architecture  compiler  writing_a_module
```

During the compiler refactor, `core/compiler.py` remains import-compatible while private stage files
are added beneath `core/_compiler/`. The current flat compiler implementation is removed only after
the differential and public-import gates pass.
