# Architecture

Three concepts, and no more.

```
      you write this                                     seqcraft does this
  ─────────────────────────                          ──────────────────────────
  Module.__init__   designs the waveforms
  Module.build()    returns a LogicBlock
  lb.add(t, ...)    places blocks in a tree     ──►   sc.compile()
                                                        finds block boundaries
                                                        sums same-axis gradients
                                                        checks the amplifier
                                                        ──► pypulseq.Sequence
```

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
| `marks` / anchors | The **module**. A block cannot know where it sits, so it cannot know when its echo occurs — and pinning the answer in would stop the same block being reused. `readout.time_to_echo` is a property. |
| `labels` field | A pulseq label is an event, so it is a node: `add(t, pp.make_label('LIN', 'SET', k))`. |
| a declared `duration` | Measured. To make a block longer, add a delay event — which is how you lengthen a block in pulseq anyway. |
| `then`, `over`, `chain`, `scaled` | `add` is shorter than any algebra and needs nothing learned. |
| `walk`, `axes`, `moment` | The compiler. Reading a tree is not part of being one. |

---

## `Module` — where the logic lives

```python
class Module(abc.ABC):
    def __init__(self, system, *, regime='default') -> None
    @property
    def opts(self) -> Opts                  # the resolved pypulseq limits
    @abc.abstractmethod
    def build(self, **args) -> LogicBlock   # the one method you write
    def params(self) -> dict                # for the provenance sidecar
```

Works like `torch.nn.Module`: `__init__` designs once, `build` is cheap and returns a value. Its
arguments select the *variant* — which diffusion lobe, which line, which slice, what phase — and must
not change the block's duration.

There are no class variables to declare, no categories to register, no hooks to override. What a
module *is*, is a thing that can build a logic block. See
[`writing_a_module.md`](writing_a_module.md).

The one thing that happens behind your back is unit validation: `check_units` runs when a subclass's
`__init__` returns, so `fov_mm=0.22` fails at construction with a hint rather than producing a
sequence with a 22 cm error in it.

---

## `sc.compile` — the scheduler

The only complicated thing in the project, deliberately: one hard thing, in one file, tested hard.
Boundaries come from where RF and ADC events fall and nowhere else; same-axis gradients are summed
with a warning; different axes are silent; limits are checked on the *compiled* waveform because that
is the only place a merge's effect is visible. Full detail in [`compiler.md`](compiler.md).

Returns a `CompiledSequence` holding the `pypulseq.Sequence`, the compile report, and per-block
provenance. There is **no `Sequence` class** in seqcraft — a sequence *is* a logic block, and
compiling it produces the artifact.

---

## The rest of `core`

| Module | What it is for |
|---|---|
| `system` | `System` holding **named `Opts` regimes**, so a diffusion encoding can run at full amplitude while the readout is derated. Asserts every regime agrees on the rasters, gamma and B0 — nothing upstream checks that, and a disagreement produces an unplayable file. |
| `raster` | Integer-picosecond time arithmetic. `10 * 1e-6 != 1e-05` in float64, and `250 * 1e-7` is not `2.5e-5`; those errors propagate into ADC dwells and off-raster block durations. |
| `events` | `derive()` — the one sanctioned way to copy a pypulseq event. Also `waveform_of`, `moment_of`, `content_hash`, `check_limits`. |
| `geometry` | FOV, matrix, slices, and one authoritative phase-encode index computation shared by `kspace_center_line` and the `LIN` label values, so the two cannot disagree. |
| `validate` | Unit plausibility bands inferred from a field's name suffix, plus explicit `require_*` helpers. |
| `report` | `Issue` and `Report`. Immutable; `ok` is true when there are no errors, so warnings inform without failing. |
| `ordering` | `interleaved_slice_order`, `centric_order`, `golden_angle`, `rf_spoil_phase` — the closed forms that otherwise get copy-pasted per notebook. |
| `provenance` | The JSON sidecar: versions, git commit and dirty flag, definitions, sha256. |
| `display` | **The only module allowed to import matplotlib**, and lazily at that, so `import seqcraft` stays cheap. Every function returns a figure and never calls `show()`. |

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

```
src/seqcraft/
  core/       logic  compiler  module  system  geometry  events  raster  units
              validate  errors  report  provenance  ordering  registry  display
  modules/    rf/  readout/  encoding/  prep/  control/
  testing.py  assertions you can point at your own modules

tests/        logic/  compiler/  modules/  integration/
examples/     01 getting started   02 DTI spiral (builds + writes it)   03 simulate + reconstruct
              seq/ (what notebook 2 writes)   lib/ (sim + recon helpers, not the package)
docs/         architecture  compiler  writing_a_module
```
