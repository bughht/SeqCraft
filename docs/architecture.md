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

Its signature is `compile(root: LogicBlock, opts: pp.Opts, ...)`, and that is the whole of what it
knows about the world upstream of it. A test asserts that `core` never imports `seqcraft.module`,
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

**Why it is not in `core`.** `core` is what gets a logic block to a validated `.seq`. A base class
for writing components is not on that path, and while it lived there it read as a requirement:
"three concepts" implied a sequence had to be expressed as modules. It never did.

**Why there is no library.** The previous one — 27 classes, 5 762 lines — grew from whichever
sequences happened to get built, and the compiler's own tests came to depend on it. It was deleted
rather than migrated. What replaces it should be *chosen*: each primitive written only when a real
sequence needs it, with the raw-pypulseq path kept beside it so the module has to earn its place.

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

Phase 1 establishes the private `PlacedEvent` and `PulseqReadyBlock` contracts. Phase 2 moves tree
traversal and absolute-time resolution into `core/_compiler/placement.py`, returning an ordered
immutable tuple without constructing a PyPulseq sequence. The remaining `compile_sequence` control
flow stays authoritative. Time arithmetic and shallow immutability policies are recorded in
[ADR-001](adr/001-compiler-internal-time-policy.md) and [ADR-002](adr/002-compiler-ir-contracts.md).

| Module | What it is for |
|---|---|
| `timing` | `Raster` — the raster as an object, with `ceil / floor / nearest / count / at / holds / require`. Arithmetic in integer ticks, because `1.5e-3 / 1e-5` is `149.99999999999997` and `250 * 1e-7` is not `2.5e-5`; both errors propagate into ADC dwells and off-raster block durations. Nothing here assumes 10 µs: rasters are read from the `Opts` at each call site. |
| `units` | One function — `convert(value, from_unit, to_unit, gamma=, f0=)` — over eleven dimensions, in both directions, the shape `pypulseq.convert` uses. Scales are exact `Fraction`s, so `4200 us` is `0.0042 s` and not `0.004200000000000001 s`. |
| `events` | `derive()` — the one sanctioned way to copy a pypulseq event. Also `knots_of` and `pwl_moment` (the exact-gradient primitives), `waveform_of` (a curve to plot, **not** a uniform raster — a `grad` event carries its own sample times), `moment_of`, `content_hash`, `check_limits`. |
| `geometry` | FOV, matrix, slices, and one authoritative phase-encode index computation shared by `kspace_center_line` and the `LIN` label values, so the two cannot disagree. Reached only through the optional `geometry=` argument, which contributes `[DEFINITIONS]`; it is a candidate to leave `core`. |
| `validate` | `merge_definitions`, which is on the compile path, plus the unit-plausibility bands `geometry` uses. Its unit names are the ones `convert` knows, asserted by a test — one vocabulary, not two. |
| `report` | `Issue` and `Report`. Immutable; `ok` is true when there are no errors, so warnings inform without failing. |

**`system.py` is gone.** 652 lines of `System`, `Limits`, presets, regimes and hardware loading,
replaced by `pp.Opts`. Everything it held is already a field of `Opts`, the compiler read exactly
eight of them, and it touched `System` at three call sites. What is genuinely lost is a named
scanner catalogue — which moves to
[PulseqSystems](https://github.com/nimpulseq/PulseqSystems), reached through
`sc.opts.from_scanner` — and the guarantee that several regimes agree, which went with the regimes.

## What sits *outside* `core`, and why

None of these is on the path from a tree to a file, so none of them is core.

| Module | What it is for |
|---|---|
| `module` | `sc.Module`: the standard shape for a reusable component. Four members. Beside `core`, not inside it, and not inside any library. |
| `scanner` | `sc.opts` — `derate` and `from_scanner`, the two operations `pp.Opts` makes awkward or unsafe. And `sc.hardware` — the PNS response model, which is not a limit and which only post-compile analysis reads. |
| `provenance` | The JSON sidecar: versions, git commit and dirty flag, definitions, the `Opts`, sha256. Output tooling; the compiler imports it lazily at `write()` time. Takes a mapping, so a component that reports its own parameters however it likes is not shut out. |
| `display` | **The only module allowed to import matplotlib**, and lazily at that, so `import seqcraft` stays cheap. Every function returns a figure and never calls `show()`. Takes arrays and an `Opts`, never a module object. |
| `testing` | Two tiers. `assert_output(make, opts)` and the block-level assertions it composes take a callable and a block, so they ask nothing about ancestry; `assert_all(module, **args)` adds the checks that only mean something for the `Module` convention. |

**There is no registry.** A registry earns its place when something must turn a *string* into a
*class* at run time: a YAML front end, plugin entry points, a `--readout=spiral_vds` flag. seqcraft
has none of those, on purpose. `module_subclasses()` serves the one real consumer — parametrising a
contract suite — with no decorator anywhere, and cannot be forgotten, because subclassing *is* the
registration. A `@register()` that a new module omitted silently lost the whole contract suite,
which is the failure a registry was supposed to prevent.

**There is no ordering module either.** `interleaved_slice_order`, `centric_order`, `golden_angle`
and the rest were removed: four of the six had never had a caller, and ordering tables are
sequence-programming choices rather than physics. They live in [`salvage/`](../salvage/) until the
opinionated library exists to hold them.

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
  core/          the compile path, and nothing else
                 logic  compiler  _compiler/  events
                 timing  units  validate  errors  report  geometry
  module.py      sc.Module -- the component contract
  scanner/       opts.py     derate, from_scanner
                 hardware.py PNS response models
  provenance.py  the JSON sidecar
  display.py     the only matplotlib importer
  testing.py     assertions you can point at any component of your own

tests/        core/  logic/  compiler/  module/  opts/  integration/
examples/     01_getting_started.ipynb   uses no modules, on purpose
              _parked/                   two DTI scans, kept as the spec for the next library
              lib/                       sim + recon helpers, not the package
salvage/      physics lifted out of the deleted library; not packaged, not imported
docs/         architecture  compiler  writing_a_module  testing
```

During the compiler refactor, `core/compiler.py` remains import-compatible while private stage files
are added beneath `core/_compiler/`. Phase 1 adds `model.py` and the verification skeleton; Phase 2
adds the authoritative `placement.py`. Legalization and emission remain in the façade until their
numbered extraction phases. The remaining flat compiler implementation is removed only after the
differential and public-import gates pass.
