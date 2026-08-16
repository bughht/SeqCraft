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

**Returns a `pypulseq.Sequence`, and nothing else.** Not a wrapper, not a pair of a sequence and
a report. `seq.write(...)`, `seq.plot()`, `seq.block_events` and `seq.definitions` are pypulseq's
own and are what the rest of the ecosystem already reads; the definitions you passed are set on it
during the compile, so nothing has to survive until write time to put them there. There is **no
`Sequence` class** in seqcraft — a sequence *is* a logic block, and compiling it produces the
artifact.

That is only tenable because **every legality failure raises**. A returned object carrying a report
is a way of not noticing: it writes a `.seq` the console refuses an hour later, and the explanation
is on an object nobody looked at. What the compile *did* rather than refused — summed two gradients
on an axis, resampled one onto the raster — is a `SeqCraftWarning`, one aggregated line per
category, so the standard `warnings` machinery decides what happens to it.

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
├── verification.py   the IR contracts, the finished-sequence checks, the compile vs the tree
├── errors.py         CompileError, HardwareLimitError, DefinitionConflict
└── model.py          the IR itself, plus block-format policy and time policy
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

## The layers

Each is named for what it answers, and they are ordered by the one direction the dependencies run:

```text
errors  ─►  design  ─►  compiler  ─►  analysis  ─►  display
   └──────────►  scanner  (independent of all five)
```

`tests/test_layering.py` asserts it, per file, from the source rather than from `sys.modules` —
so an import that only fires under `TYPE_CHECKING` still counts, and one that happens to be
satisfied by import order does not. The edge it forbids in particular is `compiler → analysis`:
the compiler must not reach for the measurements taken *of* what it produced. It would be an easy
one to add — the self-check wants a moment per axis, and `analysis.moments` computes one — and it
would be wrong, because that moment walks the tree while the self-check needs the compiled side. A
self-check that shares code with the thing it checks compares a number with itself.

### `design/` — what you build

| Module | What it is for |
|---|---|
| `logic` | `LogicBlock`, `Node`, `flatten`, `span`, `barrier`. A tree of events with relative start times, and nothing else. |
| `module` | `sc.Module`: the standard shape for a reusable component. Four members, and not a gate. |
| `events` | `derive()` — the one sanctioned way to copy a pypulseq event. Also `knots_of` and `pwl_moment` (the exact-gradient primitives), `waveform_of` (a curve to plot, **not** a uniform raster — a `grad` event carries its own sample times), `content_hash`, `check_limits`, and the kind vocabulary every stage classifies by. |
| `timing` | `Raster` — the raster as an object, with `ceil / floor / nearest / count / at / holds / require`. Arithmetic in integer ticks, because `1.5e-3 / 1e-5` is `149.99999999999997` and `250 * 1e-7` is not `2.5e-5`; both errors propagate into ADC dwells and off-raster block durations. Nothing here assumes 10 µs: rasters are read from the `Opts` at each call site. |
| `units` | One function — `convert(value, from_unit, to_unit, gamma=, f0=)` — over eleven dimensions, in both directions, the shape `pypulseq.convert` uses. Scales are exact `Fraction`s, so `4200 us` is `0.0042 s` and not `0.004200000000000001 s`. |

`sample` used to be a sixth module here. It is in `analysis` now, beside the other three ways of
measuring a tree.

`RasterError` lives in `timing.py`, with the only code that raises it. An exception stays with its
raiser unless more than one package needs it; every one is re-exported as `sc.<Name>`, and
`test_layering.py` asserts the *identity*, not just the presence, so the move is invisible to a
caller's `except`.

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

### `analysis` — measuring a tree

Four functions, **one entry shape: give it a tree, get numbers back.** You should not have to know
that PNS prediction needs a compiled sequence while a moment does not, so none of them asks for
one; `kspace` and `pns` compile internally.

| Function | Basis | Exact? |
|---|---|---|
| `sample(tree, opts)` | uniform raster grid, **interpolated** | **No** — for looking, and for approximate numeric work |
| `moments(tree, order=0)` | `knots_of` + `pwl_moment` over `flatten(tree)` | **Yes** — never routed through `sample` |
| `kspace(tree, opts)` | compiled, then `calculate_kspacePP()` | Yes, at true ADC sample times |
| `pns(tree, opts, hw)` | compiled, then `calculate_pns()` | pypulseq's validated SAFE model |

`moments` looks like it could be built on `sample` now that they sit in one file. It must not be,
and the reason is subtler than "sampling is lossy". Linear interpolation errs *antisymmetrically*
about each knot, so the halves cancel under the integral and m0 comes back bit-identical — while
the peak is visibly rounded off. A `moments` built on `sample` would therefore look correct on
every test anyone would think to write, and the compiler's own self-check would be comparing two
differently-wrong numbers. `tests/analysis/test_analysis.py` asserts both halves of that.

### `result/` is gone

It held `CompiledSequence`, `WriteResult` and the JSON provenance sidecar. `compile` returns the
`pypulseq.Sequence` directly now, so there is nothing left to wrap: `n_blocks` is
`len(seq.block_events)`, `duration_s` is `seq.duration()[0]`, `check()` is the compile raising, and
`write()` is pypulseq's.

**The sidecar leaves a real gap.** Nothing currently records the git commit, dirty flag or package
versions that produced a `.seq`. That is deferred rather than solved — see
[`serialization.md`](serialization.md) — and it starts from a blank slate rather than from sidecar
code shaped around a type that no longer exists.

### Beside all four

| Module | What it is for |
|---|---|
| `errors` | `SeqCraftError`, `ConfigurationError`, `MissingExtraError`, `SeqCraftWarning` and `format_error` — only what more than one package needs. Everything else lives with the code that raises it, and is re-exported here. |
| `scanner` | `sc.opts` — `derate` and `from_scanner`, the two operations `pp.Opts` makes awkward or unsafe. And `sc.hardware` — the PNS response model, which is not a limit and which only `analysis.pns` reads. |
| `display` | **The only module allowed to import matplotlib**, and lazily at that. One public function, `plot_block`, and it returns a figure rather than calling `show()`. It draws the *tree*; for a compiled sequence `seq.plot()` is pypulseq's own, is better maintained, and is the picture everyone else in the ecosystem reads. |

**`report.py` is gone.** `Issue` and `Report` were the soft half of "hard failures raise, soft
findings report", and the soft half turned out not to exist: every finding they carried was either
a legality failure — which now raises, because there is no legal sequence to hand back — or
something the compiler *did*, which is now a `SeqCraftWarning` through the standard machinery.

**`testing.py` is gone.** Of the five assertions it shipped, three are things the compiler now
checks with a better message and on the *summed* waveform. The two it cannot make — determinism and
purity — are in `tests/conftest.py`, because a package that ships assertions has to keep them
working and these are forty lines that only this repository's own module tests use. What they are
*for* is documented in [`writing_a_module.md`](writing_a_module.md), which is better than naming a
function.

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
  errors.py      what more than one package raises, and format_error
  analysis.py    sample, moments, kspace, pns -- measuring a tree
  display.py     plot_block; the only matplotlib importer

  scanner/       what you build against
                 opts.py      derate, from_scanner, UnknownFieldError
                 hardware.py  PNS response models

  design/        what you build
                 logic.py  module.py  events.py
                 timing.py (Raster, RasterError)  units.py

  compiler/      the transform
                 __init__.py  compile_sequence
                 placement.py  boundaries.py  legalization.py  emission.py
                 verification.py  model.py  errors.py

tests/        analysis/  design/  logic/  compiler/  module/  opts/  integration/
              conftest.py                the two checks the compiler cannot make
              test_layering.py           the layout, asserted
examples/     01_getting_started.ipynb   uses no modules, on purpose
              _parked/                   two DTI scans, kept as the spec for the next library
              lib/                       sim + recon helpers, not the package
salvage/      physics lifted out of the deleted library; not packaged, not imported
docs/         architecture  compiler  writing_a_module  testing
```

The layers are the questions in the order a sequence passes through them, and the arrow between
them only ever points forward. That is the whole of the membership rule, and unlike "is this on the
compile path?" — which put the unit table, the geometry and the report in one directory with the
scheduler — it can be checked mechanically, which is what `test_layering.py` does.
