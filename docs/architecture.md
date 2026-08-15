# Architecture

Two concepts, and no more — plus one object that is pypulseq's.

```
   any Python you like                                  seqcraft does this
  ─────────────────────────                          ──────────────────────────
  a function                ─┐
  a class of your own shape  ├─►  LogicBlock  ──►    sc.compile(tree, opts)
  an sc.Module subclass     ─┘    lb.add(t, ...)       finds block boundaries
                                                       sums same-axis gradients
                                                       checks the amplifier
                                                       ──► pypulseq.Sequence
```

**`LogicBlock` is the interface.** It is where seqcraft imposes structure, and the only place. What
produces one — a function, a class, a notebook cell — is not seqcraft's business, and nothing on the
path from a tree to a `.seq` inspects the producer's type.

**The scanner is a `pypulseq.Opts`.** Not a seqcraft class wrapping one. The compiler reads eight
fields of it; the same object configures `pp.make_trapezoid`. A part designed against derated limits
carries a second `Opts` (`sc.opts.derate`), which is one object more and one concept fewer than the
named "regimes" it replaced.

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

`add` has two call shapes and one meaning. `add(t, *items)` is one instant; `add(rows)` takes a table
of `[time, *items]` rows, for when the schedule is computed rather than written out. Rows are appended
in the order given and never sorted by time — insertion order is what `flatten` and the compiler's
tie-breaking rely on. `nodes` ends up identical either way, so nothing downstream knows which was used.

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

The only deliberately complicated subsystem in the project: one public façade, one file per stage,
and tests at each responsibility boundary. Internally this is a deterministic constrained-scheduling
algorithm — tree placement, interval legality, boundary selection, piecewise-linear gradient
transformation, and mechanical emission — and each of those is now a module with a name.
Boundaries come from where RF and ADC events fall and nowhere else under the current policy;
same-axis gradients are summed with a warning; different axes are silent; limits are checked on the
*compiled* waveform because that is the only place a merge's effect is visible. Full detail in
[`compiler.md`](compiler.md).

Returns a `CompiledSequence` holding the `pypulseq.Sequence`, the compile report, and per-block
provenance. There is **no `Sequence` class** in seqcraft — a sequence *is* a logic block, and
compiling it produces the artifact.

Its signature is `compile(root: LogicBlock, opts: pp.Opts, ...)`, and that is the whole of what it
knows about the world upstream of it. A test asserts that `compiler/` never imports
`seqcraft.design.module`, because the direction of that arrow is the design.

The implementation is one file per stage, in the order a tree passes through them:

```text
src/seqcraft/compiler/
├── __init__.py       compile_sequence — the pass that orchestrates the rest
├── placement.py      tree -> absolute-time events
├── boundaries.py     where the blocks are cut, and which readout each label addresses
├── legalization.py   superpose, represent, resample, measure limits
├── emission.py       ready blocks -> pypulseq.Sequence
├── verification.py   the IR contracts, and the compile checked against the tree
├── model.py          the IR itself, plus block-format policy and time policy
└── definitions.py    merge_definitions
```

None of the stage modules is re-exported. The public surface is `compile_sequence` and nothing
else, so a user never imports a stage type.

---

## `Module` — a contract, and no library behind it

**seqcraft ships no concrete modules.** There is no `SincExcitation`, no `EPIReadout`, no
`MonopolarDiffusion`. What it ships is the shape a reusable component takes when you write one:

```python
class Module(ABC):
    def __init__(self, *, opts: pp.Opts, tag: str | None = None) -> None
    def __call__(self, *args, **kwargs) -> LogicBlock     # the interface
    @abstractmethod
    def build(self, *args, **kwargs) -> LogicBlock        # what you write
```

Four members, and each earns its place. `opts` is required rather than defaulted, because
pypulseq's fallback is the *process-global* `Opts.default` and a module designing against that makes
a sequence depend on import order. `tag` is identity, and becomes the provenance path. `__call__` is
the single place the framework gets to check and name what came back. `build` is abstract, because a
subclass that does not produce a block is not a module.

**Three things it deliberately does not do.** It declares no `duration` — the block measures itself,
so there is no second source of truth, and a build argument may therefore change the duration
freely. It checks no units. It scrapes no `__dict__` for provenance. Each of those was a real
feature; none is required by *every* module, and a base class carrying what only some subclasses
need is how a small idea becomes a framework.

**It is also not a requirement.** The compiler takes a `LogicBlock` and never asks what produced it,
so a plain function is a component. A class with two outputs — `diffusion.pre()` and
`diffusion.post()` — names them for what they are rather than passing `build(part=…)`, and seqcraft
has no opinion about it. See [`writing_a_module.md`](writing_a_module.md).

**Why it is in `design/` and not on the compile path.** `Module` is a shape for *writing* the thing
that produces a tree, so it belongs with the tree. It is emphatically not part of the transform:
while it lived beside the compiler it read as a requirement, and "three concepts" implied a sequence
had to be expressed as modules. It never did.

**Why there is no library.** The previous one — 27 classes, 5 762 lines — grew from whichever
sequences happened to get built, and the compiler's own tests came to depend on it. It was deleted
rather than migrated. What replaces it should be *chosen*: each primitive written only when a real
sequence needs it, with the raw-pypulseq path kept beside it so the module has to earn its place.

---

## The four packages

Each is named for a question, and they are ordered by the one direction the dependencies run:

```text
errors, report  ─►  design  ─►  result  ─►  compiler
        └──────────►  scanner  (independent)      display ─► design, result
```

`tests/test_layering.py` asserts it, per file, from the source rather than from `sys.modules` —
so an import that only fires under `TYPE_CHECKING` still counts, and one that happens to be
satisfied by import order does not. Two things it forbids in particular: `result/` importing
`compiler/`, and anything under `compiler/` reaching `display`, `provenance`, `testing`,
`module` or `scanner`.

### `design/` — what you build

| Module | What it is for |
|---|---|
| `logic` | `LogicBlock`, `Node`, `flatten`, `span`, `barrier`. A tree of events with relative start times, and nothing else. |
| `module` | `sc.Module`: the standard shape for a reusable component. Four members, and not a gate. |
| `events` | `derive()` — the one sanctioned way to copy a pypulseq event. Also `knots_of` and `pwl_moment` (the exact-gradient primitives), `waveform_of` (a curve to plot, **not** a uniform raster — a `grad` event carries its own sample times), `content_hash`, `check_limits`, and the kind vocabulary every stage classifies by. |
| `timing` | `Raster` — the raster as an object, with `ceil / floor / nearest / count / at / holds / require`. Arithmetic in integer ticks, because `1.5e-3 / 1e-5` is `149.99999999999997` and `250 * 1e-7` is not `2.5e-5`; both errors propagate into ADC dwells and off-raster block durations. Nothing here assumes 10 µs: rasters are read from the `Opts` at each call site. |
| `units` | One function — `convert(value, from_unit, to_unit, gamma=, f0=)` — over eleven dimensions, in both directions, the shape `pypulseq.convert` uses. Scales are exact `Fraction`s, so `4200 us` is `0.0042 s` and not `0.004200000000000001 s`. |
| `sampling` | `sample(root, opts)` — a tree as arrays. Was private inside `display`; it is useful without drawing. |

**There is no `Geometry` class.** `compile()` used to take `geometry=` and call `definitions()` on
it, which is what `definitions=` already does for any source — so ~450 lines of dataclass and range
framework sat in the package to produce one eight-key dict. FOV, matrix and slice order are
decisions about the *scan*; the compiler turns a tree into legal pulseq blocks and is indifferent to
why the tree looks the way it does. Hand it pulseq's own keys:

```python
sc.compile(tree, opts, definitions={'FOV': [0.25, 0.25, 0.005], 'kSpaceCenterLine': 32, ...})
```

The dataclass is preserved standalone in [`salvage/geometry.py`](../salvage/geometry.py), including
the plausibility bands, because a geometry of that shape is wanted again when the module library
gets its infrastructure — a readout and a phase encoder both need FOV and matrix, and deriving the
definitions from the same fields is what stops the file disagreeing with what it plays.

### `result/` — what compile returns

| Module | What it is for |
|---|---|
| `__init__` | `CompiledSequence` and `WriteResult`. The `pypulseq.Sequence`, the report, per-block provenance, and the questions you ask afterwards: `moments`, `check`, `kspace`, `pns`, `write`. |
| `provenance` | The JSON sidecar: versions, git commit and dirty flag, definitions, the `Opts`, sha256. Written by default; `write(sidecar=False)` suppresses it. Takes a mapping, so a component that reports its own parameters however it likes is not shut out. |

`Issue` and `Report` are **not** in `result/`, though a compile returns one — five compiler stage
modules build `Issue`s long before a `CompiledSequence` exists, so keeping them here made
`compiler/` import `result/` for a type. They are at the root beside `errors`.

`CompiledSequence._verify` used to live here and took `Sequence[PlacedEvent]` — the compiler's
private IR. That single method was the only reason `result/` depended on `compiler/`. It is now
`compiler.verification.verify_against_tree`, a compile stage that happens to run last, and the
result types import nothing from the transform that made them.

### Beside all four

| Module | What it is for |
|---|---|
| `errors` | The exception hierarchy and `format_error`. At the root, because every package raises from it. |
| `report` | `Issue` and `Report` — findings as data. The other half of the pair: hard failures raise, soft findings report. Immutable; `ok` is true when there are no errors, so a 98 %-of-limit warning informs without failing. |
| `scanner` | `sc.opts` — `derate` and `from_scanner`, the two operations `pp.Opts` makes awkward or unsafe. And `sc.hardware` — the PNS response model, which is not a limit and which only post-compile analysis reads. |
| `display` | **The only module allowed to import matplotlib**, and lazily at that. Every function returns a figure and never calls `show()`. It stays at the root because it spans two stages by signature: `plot_block(LogicBlock)` against `plot_sequence(CompiledSequence)`. |
| `testing` | `assert_output(make, opts)` takes a callable and a block, so it asks nothing about ancestry; `assert_all(module, **args)` adds the checks that only mean something for the `Module` convention. |

**`system.py` is gone.** 652 lines of `System`, `Limits`, presets, regimes and hardware loading,
replaced by `pp.Opts`. Everything it held is already a field of `Opts`, the compiler read exactly
eight of them, and it touched `System` at three call sites. What is genuinely lost is a named
scanner catalogue — which moves to
[PulseqSystems](https://github.com/nimpulseq/PulseqSystems), reached through
`sc.opts.from_scanner` — and the guarantee that several regimes agree, which went with the regimes.

**There is no registry.** A registry earns its place when something must turn a *string* into a
*class* at run time: a YAML front end, plugin entry points, a `--readout=spiral_vds` flag. seqcraft
has none of those, on purpose. `module_subclasses()` was written to parametrise a contract suite
over them; no such suite was ever written and it had no callers, so it is gone too.

**There is no ordering module, and no phase-encode table.** `interleaved_slice_order`,
`centric_order`, `golden_angle` and the rest were removed: four of the six had never had a caller,
and ordering tables are sequence-programming choices rather than physics. `Geometry.pe_lines`
followed them for the same reason — it was defended as the one computation `kspace_center_line`
and the `LIN` label would share, but nothing consumed it once the module library was deleted. Both
live in [`salvage/`](../salvage/) until the opinionated library exists to hold them.

---

## What is not in the package

`examples/lib/` holds the MRzero simulation bridge and the sigpy off-resonance reconstruction. They
are **not** part of seqcraft, on purpose: the package builds sequences, and simulating or
reconstructing them are downstream jobs with heavy dependencies of their own. Folding them in would
make every user pay for them.

**There are no recipes either.** A recipe is somebody else's sequence choices in library code: to
change your scan you would edit a package, and the package would accumulate everyone's variants. A
sequence is assembled where it is used, and the timing is arithmetic you can read --

```python
first_half  = (exc_block.duration - isodelay) + delta + isodelay_180
second_half = (refoc_duration - isodelay_180) + delta + time_to_echo
te = 2 * max(first_half, second_half)
```

-- so a `solve_timing` function was never earning its place.

Also deliberately absent: any solver for gap timing (`t0 + te - time_to_echo` is arithmetic you can
read), a YAML or GUI front end, and a `.seq` importer.

---

## Layout

```
src/seqcraft/
  __init__.py    the façade -- the only global re-export layer
  errors.py      the exception hierarchy      } how a problem
  report.py      Issue and Report             } is communicated
  display.py     the only matplotlib importer
  testing.py     assertions you can point at any component of your own

  scanner/       what you build against
                 opts.py      derate, from_scanner
                 hardware.py  PNS response models

  design/        what you build
                 logic.py  module.py  events.py  sampling.py
                 timing.py  units.py

  compiler/      the transform
                 __init__.py  compile_sequence
                 placement.py  boundaries.py  legalization.py  emission.py
                 verification.py  model.py  definitions.py

  result/        what compile returns
                 __init__.py  CompiledSequence, WriteResult
                 provenance.py

tests/        design/  logic/  compiler/  module/  opts/  integration/
              test_layering.py           the layout, asserted
examples/     01_getting_started.ipynb   uses no modules, on purpose
              _parked/                   two DTI scans, kept as the spec for the next library
              lib/                       sim + recon helpers, not the package
salvage/      physics lifted out of the deleted library; not packaged, not imported
docs/         architecture  compiler  writing_a_module  testing
```

The four packages are the four questions in the order a sequence passes through them, and the
arrow between them only ever points forward. That is the whole of the membership rule, and unlike
"is this on the compile path?" — which put the unit table, the geometry and the report in one
directory with the scheduler — it can be checked mechanically, which is what `test_layering.py`
does.
