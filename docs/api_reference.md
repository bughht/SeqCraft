# seqcraft API reference

Every public name in `src/seqcraft`, by layer. Signatures are the ones the code has; every code
block below is executed by `tools/check_api_reference.py`, so this document cannot drift from the
package without CI noticing.

Read §0 first — the rest is reference.

- [0. The whole package in one page](#0-the-whole-package-in-one-page)
- [1. `design/` — what you build](#1-design--what-you-build)
- [2. `scanner/` — what you build against](#2-scanner--what-you-build-against)
- [3. `compiler/` — the transform](#3-compiler--the-transform)
- [4. `analysis` — measuring a tree](#4-analysis--measuring-a-tree)
- [5. `display` — looking at a tree](#5-display--looking-at-a-tree)
- [6. `errors` — the exception hierarchy](#6-errors--the-exception-hierarchy)
- [7. Compiler internals](#7-compiler-internals)
- [8. Migration](#8-migration)
- [9. Index](#9-index)

---

## 0. The whole package in one page

Three things you write, one function that transforms them, one toolbox that measures them.

```python
import numpy as np
import pypulseq as pp
import seqcraft as sc

opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
               rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)

gx = pp.make_trapezoid('x', flat_area=64 * 4.0, flat_time=3.2e-3, system=opts)
adc = pp.make_adc(num_samples=64, duration=3.2e-3, delay=gx.rise_time, system=opts)

tr = sc.LogicBlock('tr')                                   # 1. a tree of events
tr.add(0.0, gx, adc)                                       #    anything may overlap anything

seq = sc.compile(tr, opts)                                 # 2. legal pulseq blocks, or an exception
assert type(seq).__module__.startswith('pypulseq')         #    it is a pypulseq.Sequence

m = sc.moments(tr)                                         # 3. measure the tree
assert round(m['x']) == 258      # 256 of flat area, plus what the two ramps contribute
```

| Layer | Module | You use it to |
|---|---|---|
| **design** | `logic` `module` `events` `timing` `units` | say what you mean |
| **scanner** | `opts` `hardware` | say what the machine can do |
| **compiler** | `compile()` | turn meaning into legal pulseq blocks |
| **analysis** | `sample` `moments` `kspace` `pns` | measure a tree |
| **display** | `plot_block` | look at a tree |

**The central contract.** `sc.compile(tree, opts)` returns a `pypulseq.Sequence` and nothing else.
If the tree cannot become a legal sequence, it **raises**. If the compiler had to change a waveform
to make it legal, it **warns**. There is no report object to inspect and no result wrapper to unpack.

### The complete façade

Everything reachable as `sc.<name>`. This is the whole of what `__init__.py` re-exports.

| Name | Kind | Section |
|---|---|---|
| `compile`, `compile_sequence` | function | [3.1](#31-sccompileroot-opts--name-definitions---pypulseqsequence) |
| `LogicBlock`, `Node`, `Item` | class, class, type alias | [1.1](#11-logicblock) |
| `flatten`, `span`, `barrier` | function | [1.2](#12-flatten-span-barrier) |
| `Module` | class | [1.3](#13-module) |
| `events` | module | [1.4](#14-designevents--pypulseq-event-arithmetic) |
| `timing`, `Raster` | module, class | [1.5](#15-designtiming--exact-time-arithmetic) |
| `units`, `convert` | module, function | [1.6](#16-designunits--one-conversion-function) |
| `scanner`, `opts`, `hardware` | module | [2](#2-scanner--what-you-build-against) |
| `analysis`, `sample`, `moments`, `kspace`, `pns` | module, function | [4](#4-analysis--measuring-a-tree) |
| `display`, `plot_block` | module, function | [5](#5-display--looking-at-a-tree) |
| `SeqCraftError`, `ConfigurationError`, `MissingExtraError` | exception | [6](#6-errors--the-exception-hierarchy) |
| `CompileError`, `HardwareLimitError`, `DefinitionConflict` | exception | [6](#6-errors--the-exception-hierarchy) |
| `RasterError`, `UnknownFieldError`, `CompilerContractError` | exception | [6](#6-errors--the-exception-hierarchy) |
| `SeqCraftWarning` | warning | [3.3](#33-warnings) |
| `__version__` | str | — |

`display` is resolved on first access, so `import seqcraft` does not pull matplotlib in through
seqcraft. Everything else is eager.

---

# 1. `design/` — what you build

Nothing here knows what a pulseq block is. A logic block lets anything overlap anything, a raster
quantises a time, a unit converts a number — and none of them has an opinion about where the
boundaries will fall.

## 1.1 `LogicBlock`

A tree of pulseq events with relative start times. Two attributes, one method.

```python
class LogicBlock:
    tag: str                 # optional label; becomes the provenance path and appears in errors
    nodes: list[Node]        # Node(start: float, item: Event | LogicBlock)
```

`Node.start` is measured **from this block's own start**, so a block does not know where it sits.
That is what lets the same block be reused at three different times in one sequence.

### `add(start, *items)` / `add(rows) -> Self`

Two call shapes, one meaning. Returns `self`, so calls chain.

```python
lb = sc.LogicBlock('excitation')
lb.add(0.0, gx)                            # one item at one instant
lb.add(4.0e-3, gx, adc)                    # two items at one instant

lb2 = sc.LogicBlock('table')               # the table form, for a computed schedule
lb2.add([[0.0, gx],
         [5.0e-3, gx, adc]])
assert len(lb2) == 3             # three nodes: every item gets its own, whatever the row shape
```

Rows are **appended in the order given and never sorted by time.** Insertion order is what
`flatten` and the compiler's tie-breaking rely on, so two events at the same instant keep the order
you wrote.

`add` validates: anything that is neither a pypulseq event nor a `LogicBlock` raises
`ConfigurationError` naming the type it got. You therefore cannot construct a malformed tree.

```python
try:
    sc.LogicBlock('bad').add(0.0, 42)
except sc.ConfigurationError as err:
    assert 'int' in str(err)
```

### `duration` — a read-only property

**Measured** from the nodes: `max(start + span(item))`. A block cannot advertise a length that
disagrees with what it plays.

```python
spoil = sc.LogicBlock('spoil').add(0.0, pp.make_trapezoid('z', area=500.0, system=opts))
assert round(spoil.duration * 1e3, 2) == 0.57
```

To make a block *longer*, add a delay — which is how you lengthen a block in pulseq anyway:

```python
spoil.add(0.0, pp.make_delay(10e-3))
assert round(spoil.duration * 1e3, 3) == 10.0
```

### The rest is Python's

| Call | Does |
|---|---|
| `lb.nodes[2].start += 40e-6` | move the third child later |
| `del lb.nodes[0]` | remove the first |
| `for node in lb:` | `__iter__` yields `Node`s |
| `len(lb)` | number of **direct** children |
| `lb.copy()` | a new block with a **new node list, sharing the same items** — see below |
| `lb.describe(indent=0)` | indented text view of the whole tree |
| `repr(lb)` | `LogicBlock(tag, n nodes, d ms)` |

```python
tree = sc.LogicBlock('tr').add(0.0, sc.LogicBlock('inner').add(0.0, gx))
assert len(tree) == 1                               # one direct child, which is a block
assert 'inner' in tree.describe()

variant = tree.copy()
variant.nodes[0].start = 1e-3                       # the node list is new: this is safe
assert tree.nodes[0].start == 0.0
assert variant.nodes[0].item is tree.nodes[0].item  # but the *items* are shared
```

**`copy()` is one level deep.** It gives you a new node list over the same items, which is what you
want when adding one block object at many times: 1700 references to one spoiler rather than 1700
copies. Reordering, retiming or deleting nodes on the copy leaves the original alone; reaching
through a node into a nested block does not. For an independent sub-tree, copy that block too.

### `Node`

```python
class Node:
    start: float             # seconds, from the enclosing block's own start
    item: Item               # a pypulseq event, or a nested LogicBlock
```

A frozen-in-spirit pair; `start` is writable, which is how you move a child.

`Item` is the type alias `SimpleNamespace | LogicBlock`.

### What it deliberately lacks

| Absent | Where the job went |
|---|---|
| `marks` / anchors | Whatever produced the block. A block cannot know where it sits, so it cannot know when its echo occurs — `readout.time_to_echo` is a property of the readout. |
| a `labels` field | A pulseq label *is* an event: `lb.add(t, pp.make_label('LIN', 'SET', k))`. |
| a declared `duration` | Measured, so there is no second source of truth. |
| `then` / `over` / `chain` | `add` is shorter than any algebra and needs nothing learned. |

## 1.2 `flatten`, `span`, `barrier`

### `flatten(root, t0=0.0, path=()) -> Iterator[tuple[float, Event, tuple[str, ...]]]`

Walks the whole tree and yields every **event** (never a block) with its absolute time from the root
and its provenance path — the chain of non-empty tags enclosing it.

```python
inner = sc.LogicBlock('spoiler').add(0.0, gx)
outer = sc.LogicBlock('tr').add(2e-3, inner)
assert [(round(t, 5), p) for t, _, p in sc.flatten(outer)] == [(0.002, ('tr', 'spoiler'))]
```

This is the function every measurement is built on. Use it, not `lb.nodes`, whenever you need the
events themselves — a component that nests would otherwise have its actual events skipped and your
check would pass vacuously.

`t0` and `path` are the recursion's own accumulators; callers pass neither.

### `span(item) -> float`

The duration of an event or a block, in seconds. For an event this is pypulseq's own
`calc_duration` logic; for a block it is `block.duration`.

```python
assert sc.span(gx) == pp.calc_duration(gx)
assert sc.span(spoil) == spoil.duration
```

### `barrier(tag='barrier') -> SimpleNamespace`

A zero-duration marker event. Adding one forces the compiler to place a block boundary at that
instant:

```python
long_g = pp.make_trapezoid('x', area=2000.0, duration=4e-3, system=opts)
split = sc.compile(sc.LogicBlock('t').add(0.0, long_g).add(2e-3, sc.barrier('mid')), opts)
assert len(split.block_events) == 2
```

Use it when a boundary must fall somewhere the compiler would not otherwise choose — for instance
so that a trigger lands at the start of its own block, or so a later reconstruction step can find a
known seam.

`seqcraft.design.logic.BARRIER` is the `type` string these carry (`'seqcraft_barrier'`). It is in
`HANDLED_KINDS` and is not a pypulseq type.

## 1.3 `Module`

The standard shape for a reusable component. **A convention, not a gate:** the compiler takes a
`LogicBlock` and never asks what produced it, so a plain function is a component too.

```python
class Module(ABC):
    def __init__(self, *, opts: pp.Opts, tag: str | None = None) -> None
    def __call__(self, *args, **kwargs) -> LogicBlock      # what callers use
    @abstractmethod
    def build(self, *args, **kwargs) -> LogicBlock         # what you write
```

Four members, each earning its place:

- **`opts` is required, not defaulted.** pypulseq falls back to the *process-global* `Opts.default`,
  so a module designing against that makes a sequence depend on import order.
- **`tag` is identity** and becomes the provenance path. An untagged block returned by `build` is
  named after the class.
- **`__call__` is the single place** the framework checks and names what came back. It rejects a
  `build` that returned the wrong type, with *your* class in the message.
- **`build` is abstract**, because a subclass that does not produce a block is not a module.

```python
class PhaseEncode(sc.Module):
    def __init__(self, *, opts, fov_m, matrix, duration_s, tag=None):
        super().__init__(opts=opts, tag=tag)
        self.matrix = matrix
        # Design once, in __init__: the *shape* is fixed, only the amplitude varies per line.
        self.blip = pp.make_trapezoid('y', area=1.0 / fov_m, duration=duration_s, system=opts)

    def build(self, line: int) -> sc.LogicBlock:
        k = line - self.matrix // 2                    # signed index about k = 0
        scaled = sc.events.derive(self.blip, amplitude=self.blip.amplitude * k)
        return sc.LogicBlock().add(0.0, scaled)

pe = PhaseEncode(opts=opts, fov_m=0.22, matrix=64, duration_s=600e-6)
assert pe(line=17).tag == 'PhaseEncode'                        # named after the class
assert PhaseEncode(opts=opts, fov_m=0.22, matrix=64,
                   duration_s=600e-6, tag='ky')(line=17).tag == 'ky'
```

**The rule that matters:** `build` must not mutate `self`. Note `sc.events.derive` above rather than
`self.blip.amplitude = ...` — see §1.4. The compiler cannot check this: it validates a tree and
never sees the second call, so a module that mutates itself hands it a perfectly legal tree every
time. The three-line check is in [`writing_a_module.md`](writing_a_module.md).

Three things `Module` deliberately does not do: declare a `duration` (the block measures itself),
check units, or scrape `__dict__` for provenance. Each was a real feature; none is needed by *every*
module, and a base class carrying what only some subclasses need is how a small idea becomes a
framework.

`_finalize(block)` is the private hook `__call__` runs: it type-checks and names. Tests call it
directly; nothing else should.

## 1.4 `design/events` — pypulseq event arithmetic

seqcraft does **not** wrap pypulseq events. They are `SimpleNamespace` objects and stay that way, so
interop with the rest of the pulseq ecosystem is free. `Event` is an alias for `SimpleNamespace`,
not a class. What this module adds is the arithmetic the compiler needs and one safe way to copy.

### `derive(event, **changes) -> Event`

**The only sanctioned way to modify a pypulseq event.**

```python
flipped = sc.events.derive(gx, amplitude=-gx.amplitude)
assert gx.amplitude == -flipped.amplitude
```

Two properties of pypulseq 1.5 make ad-hoc copying dangerous:

1. `Sequence.set_block` writes a registration cache **onto the event object**
   (`event._pypulseq_sequence_event_cache`), keyed by `id(sequence)`. A `copy.deepcopy` carries that
   stale key along, and CPython reuses addresses — so a copied event can be silently mistaken for a
   registered one in a *different* sequence.
2. Assigning to a shared event mutates it for every tree that holds it.

`derive` strips `id`, `shape_IDs` and `_pypulseq_sequence_event_cache`, then applies your changes to
a fresh namespace. **Never write `event.amplitude = x`.**

It raises `AttributeError` for a name the event does not already have, so a typo is caught rather
than silently adding a field:

```python
try:
    sc.events.derive(gx, amplitide=1.0)          # note the typo
except AttributeError as err:
    assert 'amplitide' in str(err)
```

Waveform arrays are shared by reference on purpose — they are never mutated in place, and sharing is
what lets pypulseq's per-event registration cache hit for repeated events.

### `knots_of(event, t0=0.0) -> (times, amps)`

The exact-gradient primitive. A gradient is piecewise linear; this returns its corner points in
absolute time.

For an **arbitrary** gradient the samples sit at raster *centres*, so `first` and `last` appear as
two extra knots half a raster outside them. For an **extended trapezoid** `tt` already reaches both
edges, so they coincide with existing knots and are not duplicated. Getting that wrong is what makes
an arbitrary gradient appear to lose amplitude.

### `pwl_moment(times, amps, order=0) -> float`

Integrates the polyline analytically. `order` 0 gives area (1/m), 1 the first moment (s/m), 2 the
second (s²/m).

```python
t, a = sc.events.knots_of(gx)
assert round(sc.events.pwl_moment(t, a, 0)) == 258          # area, 1/m: flat top plus ramps
```

**Use these, never raster sampling, for any arithmetic.** An arbitrary gradient's samples sit at
raster *centres* and an extended trapezoid's knots are not uniformly spaced, so interpolating onto a
uniform grid costs percent-level amplitude error on a spiral. The compiler uses knots throughout; so
should you. See §4.5.

Moments are **linear in the waveform**, which is why the compiler can sum per-piece moments whether
the pieces overlap (in a tree, where they superpose) or abut (in compiled blocks, where they
concatenate). No union of knots needs building.

The closed form matters: `g(t) * t` is quadratic between knots, so the trapezoidal rule is *wrong*
for `order >= 1`.

### `waveform_of(event, raster) -> (times, amps)`

A curve suitable for plotting or interpolating — *not* a uniform raster. A `trap` is sampled onto
`raster`; a `grad` event carries its own sample times and they come back unchanged.

**Never assume the returned times are uniform**; use the `t` that comes back. `check_limits` once did
not, and divided by `raster` where it should have divided by `diff(t)` — on an EPI train's extended
trapezoid that overstated the slew by 26×.

```python
t, g = sc.events.waveform_of(gx, opts.grad_raster_time)
assert bool(np.allclose(np.diff(t), opts.grad_raster_time))         # a trapezoid: uniform
```

### `content_hash(event) -> str`

A stable sha256 of an event's numeric content, excluding registration state. This is what makes
mutation detectable:

```python
before = sc.events.content_hash(pe.blip)
pe(line=17)
assert sc.events.content_hash(pe.blip) == before, 'build() mutated a stored event'
```

### `check_limits(events, opts, raster=0.0, *, starts=None) -> list[tuple[str, str, float, float]]`

Amplitude and slew against `opts`, per axis and as a vector norm. Each entry is
`(kind, where, value, limit)`, where `kind` is one of `grad`, `slew`, `grad_norm`, `slew_norm`.

```python
strong = {'amplitude': 0.9 * opts.max_grad, 'duration': 1e-3, 'system': opts}
sx, sy = pp.make_trapezoid('x', **strong), pp.make_trapezoid('y', **strong)
assert [k for k, *_ in sc.events.check_limits([sx], opts)] == []          # legal alone
assert 'slew_norm' in [k for k, *_ in sc.events.check_limits([sx, sy], opts)]
```

**Pass `starts=` whenever the events do not all begin at zero.** Two lobes of a bipolar pair are both
on one axis; without their start times they are taken to play simultaneously and sum to zero, which
hides a real violation.

The `_norm` entries are informational for a tree: the vector norm across simultaneous axes routinely
exceeds the per-axis limit and is legal on real amplifiers. The compiler treats per-axis violations
as **errors** and reports norm violations as **warnings**.

### The event vocabulary

Classifications of pulseq's own `type` strings, kept here because they describe what an event *is*:

| Constant | Members |
|---|---|
| `GRADIENT_KINDS` | `trap`, `grad` |
| `POINT_KINDS` | `labelset`, `labelinc`, `trigger`, `output` |
| `LABEL_KINDS` | `labelset`, `labelinc` |
| `HANDLED_KINDS` | every type the compiler can place: the above plus `rf`, `adc`, `delay`, `seqcraft_barrier` |
| `AXES` | `('x', 'y', 'z')` |
| `ADDRESS_KEYS` | `('SLC', 'LIN', 'PAR', 'AVG', 'REP', 'SEG', 'ECO', 'SET')` — together they address a k-space location |

The two genuinely *compiler* constants — what may not be cut, what may not share a block — live in
`compiler/model.py` instead, because they are block-format policy rather than event identity.

> `events.trapz` is a numpy-compat alias (`np.trapezoid` on numpy ≥ 2, `np.trapz` below), kept so the
> declared `numpy>=1.24` floor is real. Deliberately absent from `__all__`: it is numpy's function
> under a compatibility name, not seqcraft API.

## 1.5 `design/timing` — exact time arithmetic

Pulseq quantises time on four rasters, every one supplied by the scanner. Floating point is not good
enough here: `1.5e-3 / 1e-5` is `149.99999999999997` and `250 * 1e-7` is not `2.5e-5`. Both errors
propagate into ADC dwells and off-raster block durations.

### `Raster(dt, name='')`

The raster as an object. All arithmetic goes through integer ticks (10¹² per second) internally.

```python
r = sc.Raster(opts.grad_raster_time, 'gradient')
assert r.holds(1.5e-3)                       # is this a multiple?
assert r.ceil(1.5049e-3) == 0.00151          # round up to the raster
assert r.floor(1.5049e-3) == 0.0015
assert r.nearest(1.5049e-3) == 0.0015
assert r.count(1.5e-3) == 150                # how many rasters is this?
assert r.at(150) == 0.0015                   # the time of raster 150
assert repr(r) == 'Raster(gradient, 10 us)'
```

`require(t, *, what='duration')` raises `RasterError` naming the value and the raster:

```python
try:
    r.require(1.505e-3)
except sc.RasterError as err:
    assert str(err).startswith('duration 1.505000 ms is not a multiple of the 10.0 us gradient raster')
```

**Nothing here assumes 10 µs.** Siemens' 10 µs, GE's 4 µs and Philips' 6.4 µs are all just numbers
read from the `Opts` at each call site.

### `to_ticks` / `from_ticks` / `exact_sum` / `exact_diff`

Integer-tick arithmetic, for when a sum of times must land exactly on a raster:

```python
t = 0.0                                                 # what a TR loop does
for _ in range(1000):
    t += 1.997e-3
assert t != 1.997                                       # and it has drifted by 23 fs
assert sc.timing.exact_sum([1.997e-3] * 1000) == 1.997  # ticks do not

terms = [1e-5] * 121                                    # all exact multiples of 10 us
assert sum(terms) != 0.00121
assert sc.timing.exact_sum(terms) == 0.00121
assert sc.timing.exact_diff(3e-3, 1e-3) == 0.002
assert sc.timing.to_ticks(1.5e-3) == 1_500_000_000
assert sc.timing.from_ticks(1_500_000_000) == 0.0015
```

The accumulating loop is the honest demonstration, and it is the shape the real bug had: a plain
subtraction of absolute times produced an RF delay of 129.9999999986 µs at 39 s into a sequence,
which pypulseq rejects.

> Do not demonstrate this with a short `sum()`. **CPython 3.12 gave `sum()` Neumaier compensated
> summation**, so `sum([1.3e-4, 5.0e-4, 3.7e-4, 2.1e-4])` is exactly `0.00121` there and drifts on
> 3.11 — a claim that flips between interpreters. `+=` gets no such compensation, which is both why
> the loop above is stable and why the accumulation in real code is the thing to worry about.

| Constant | Value | Meaning |
|---|---|---|
| `EPS` | `1e-9` | the compiler's time tolerance — one nanosecond, well below any raster |
| `TICKS_PER_SECOND` | `10**12` | one tick is a picosecond; float64 holds integers exactly to 2.5 hours of them |

`RasterError` is defined here, with the only code that raises it, and re-exported as
`sc.RasterError`.

## 1.6 `design/units` — one conversion function

```python
convert(value, from_unit, to_unit=None, *, gamma=GAMMA_1H, f0=None) -> float
```

Eleven dimensions, both directions, the shape `pypulseq.convert` uses. Scales are exact `Fraction`s,
so `convert(4200, 'us', 's')` is `0.0042` and not `0.004200000000000001`.

```python
assert sc.convert(220, 'mm', 'm') == 0.22
assert sc.convert(40, 'mT/m', 'Hz/m') == 1703040.0          # at the proton gamma
assert round(sc.convert(20, 'uT', 'Hz'), 2) == 851.52       # B1, carried in hertz
assert sc.convert(220, 'mm') == 0.22                        # to_unit=None -> the SI base unit
```

Use it so a module never spells out a factor.

| Helper | Returns |
|---|---|
| `dimension_of(unit)` | the dimension a unit belongs to, e.g. `'gradient'` for `'mT/m'` |
| `known_units(dimension=None)` | every unit, or every unit of one dimension |
| `dimensions()` | `('angle', 'bvalue', 'frequency', 'gradient', 'kspace', 'kspace_area', 'kspace_rate', 'length', 'ratio', 'slew', 'time')` |
| `GAMMA_1H` | `42576000.0` Hz/T |

```python
assert sc.units.dimension_of('mT/m') == 'gradient'
assert 'ms' in sc.units.known_units('time')
assert len(sc.units.dimensions()) == 11
```

An unknown unit raises `ConfigurationError` listing what is available.

---

# 2. `scanner/` — what you build against

**The scanner is a `pypulseq.Opts`.** Not a seqcraft class wrapping one. The compiler reads eight
fields of it, and the same object configures `pp.make_trapezoid`. Build it the ordinary way:

```python
opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s', B0=3.0,
               rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
```

> **Set the dead times.** pypulseq defaults `rf_dead_time`, `rf_ringdown_time` and `adc_dead_time` to
> **zero**, which is wrong on every real scanner: the sequence compiles cleanly, validates cleanly,
> and is refused or silently mangled at the console. They are properties of your *installation*, so
> no preset and no vendor database can supply them.

The eight fields the compiler reads:

```text
opts.max_grad          opts.rf_ringdown_time    opts.gamma
opts.max_slew          opts.rf_raster_time      opts.adc_dead_time
opts.grad_raster_time  opts.adc_raster_time     (and opts.block_duration_raster)
```

## 2.1 `sc.opts.derate(opts, *, grad=1.0, slew=1.0, **overrides) -> Opts`

Because deriving one `Opts` from another by hand loses the rest of it — `Opts` fills every argument
you omit from the *process-global* `Opts.default`, so the hand-written version silently returns to
zero dead times, pypulseq's rasters and a foreign gamma.

```python
soft = sc.opts.derate(opts, grad=0.85, slew=0.70)       # design a diffusion lobe against this
assert soft.max_grad == 0.85 * opts.max_grad
assert soft.rf_dead_time == opts.rf_dead_time           # the hand-written version loses this
assert pp.Opts(max_grad=opts.max_grad * 0.85).rf_dead_time == 0.0
```

`derate` copies all eighteen fields and changes only what it was asked to. `**overrides` sets any
other `Opts` field; an unknown one raises `UnknownFieldError` with a spelling suggestion.

The finished sequence is still compiled against the full `opts` — the derated one is a design
constraint, not a validation ceiling.

## 2.2 `sc.opts.from_scanner(...) -> Opts`

```python
from_scanner(manufacturer, model, gradient=None, *,
             rf_dead_time, rf_ringdown_time, adc_dead_time, max_b1, **overrides) -> Opts
```

Looks `max_grad`, `max_slew` and `B0` up in
[PulseqSystems](https://github.com/nimpulseq/PulseqSystems) so they need not be copied off a spec
sheet. The four site constants are **required keyword arguments** because no vendor database has
them — they belong to the installation.

```python
opts = sc.opts.from_scanner(
    'Siemens', 'Prisma', 'XR',
    rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
    max_b1=sc.convert(20, 'uT', 'Hz'),
)
```

Raises `UnknownFieldError` for an override that is not an `Opts` field, and `MissingExtraError` if
the `systems` extra is not installed.

`UnknownFieldError` is defined here and re-exported as `sc.UnknownFieldError`.

## 2.3 `sc.hardware` — PNS response models

**A hardware model is not a limit.** `Opts` says how hard the amplifier may be driven; this
describes how the *body* responds — three exponential time constants per axis, a stimulation
threshold, and the forbidden acoustic-resonance bands. Nothing on the compile path reads it;
`sc.pns(tree, opts, hw)` takes it as a third argument.

### `synthetic_hardware(name='synthetic_generic') -> SimpleNamespace`

```python
hw = sc.hardware.synthetic_hardware()
assert hw.is_synthetic is True
assert 'NOT a real scanner' in repr(hw)
```

> **`synthetic_hardware()` is not a real scanner.** It is a conservative stand-in so PNS checks can
> run without a vendor file, shaped like `pypulseq.utils.siemens.asc_to_hw`'s output but carrying
> the illustrative coefficients from pypulseq's own `safe_example_hw()`. It **must never be used to
> clear a human scan.** On one DTI example it reported 2.44 where the site's own descriptor reported
> 0.95 — the difference between "not runnable" and "passes with a thin margin".
>
> The caveat travels with the object: `is_synthetic=True` and a `repr` that says so.

### `load_hardware(filename, *, cardiac_model=False, directory=None) -> SimpleNamespace`

**No vendor hardware file is ever read from inside this repository.** Siemens `.asc` descriptors
carry proprietary response coefficients, so this resolves them through the `SEQCRAFT_ASC_DIR`
environment variable only (`ASC_ENV_VAR` is that name). `filename` must be a bare name; a path
raises `ConfigurationError`. The returned model carries `.source`, a provenance string of the form
`'<filename> sha256:<12 hex>'` — the file *name* and hash only, never the contents.

```python
hw = sc.hardware.load_hardware('CimaX.asc')             # needs $SEQCRAFT_ASC_DIR
```

---

# 3. `compiler/` — the transform

## 3.1 `sc.compile(root, opts, *, name='', definitions=None) -> pypulseq.Sequence`

Also spelled `sc.compile_sequence`, for callers who would rather not shadow the builtin.

```python
seq = sc.compile(tr, opts, name='gre', definitions={'TE': 8e-3, 'TR': 20e-3})
assert seq.definitions['TE'] == 8e-3
assert seq.definitions['Name'] == 'gre'
assert 'TotalDuration' in seq.definitions
```

**Returns a `pypulseq.Sequence`.** Yours to use directly — `seq.write(path)`, `seq.plot()`,
`seq.block_events`, `seq.duration()[0]` are all pypulseq's own. The definitions you passed are
already set on it, along with `Name` and `TotalDuration`, so nothing has to survive until write time.

### What it does

A logic block lets anything overlap anything. A pulseq block holds **at most one RF, one ADC and one
gradient per axis**, must last a whole number of block rasters, must leave the RF dead time before a
pulse and the ringdown after it, and must have gradients that join continuously across boundaries.
The compiler finds boundaries satisfying all of that and combines whatever lands inside each one.

**Boundaries come from two independent properties**, which the reservation model keeps apart:

- **Indivisible** — no boundary may fall strictly inside it. An RF's dead time and ringdown, an ADC's
  window and trailing dead time, a trigger's pulse: each is one hardware action, so splitting it is
  meaningless.
- **Exclusive** — a block may hold at most one, so a boundary is *required* between two. RF and ADC
  only. Triggers are `TRIGGERS` extensions and pulseq accepts several per block.

Putting a boundary somewhere in the gap between each pair of consecutive exclusive reservations
guarantees at most one RF and one ADC per block **by construction**, with nothing left to check.

### Three rules on overlap

| Case | Behaviour | Why |
|---|---|---|
| Gradients on **different** axes | silent | A phase-encode blip beside a slice rewinder beside a readout prephaser is the normal way to build a sequence. Warning would teach people to ignore warnings. |
| Gradients on the **same** axis | **warned, then summed** | Summing is almost always what was meant — a rewinder sharing time with a prephaser. But it is the one place the compiler changes a waveform, so it says so. |
| Two RF or two ADC overlapping | **raises** | You cannot transmit twice at once, sample twice at once, or transmit and receive at once. |

### Limits are checked *after* summing

Two individually legal gradients on one axis can sum to an illegal one, and no module can see that in
isolation: adding an area-100 and an area-200 trapezoid on a 40 mT/m, 150 T/m/s system reaches 93 % of
the amplitude limit but **189 %** of the slew limit. So amplitude and slew are measured on the
*compiled* waveform, which is the only place the truth is visible.

### `definitions=`

Pulseq's own vocabulary, exactly as it will be written — `FOV`, `SliceThickness`, `kSpaceCenterLine`,
`TE`, `TR`, and anything else your acquisition wants recorded. Merged with the sequence name under a
collision check: passing `Name` with a value different from `name`/`root.tag` raises
`DefinitionConflict` rather than one silently winning.

```python
try:
    sc.compile(tr, opts, name='gre', definitions={'Name': 'something_else'})
except sc.DefinitionConflict as err:
    assert 'gre' in str(err) and 'something_else' in str(err)
```

**There is no geometry argument.** FOV, matrix and slice order are decisions about the scan you are
running; the compiler turns a tree into legal blocks and is indifferent to why the tree looks the way
it does. `salvage/geometry.py` holds a `Geometry` dataclass with a `definitions()` method if you want
one.

## 3.2 The exception contract

**Every legality failure raises.** There is no report to forget to check.

| Exception | Raised when | Fix |
|---|---|---|
| `CompileError` | Two RF or two ADC overlap; an absolute start is negative; a gradient starts off the gradient raster; no boundary can be cut in the gap between two exclusive events; a block boundary would fall inside a gradient an ADC is sampling; two ADCs write the same k-space address; an unsupported or unknown event type; `check_timing` fails | Fix the tree — the message names the event, its provenance path, and usually two concrete remedies |
| `HardwareLimitError` | The *summed* waveform exceeds `max_grad` or `max_slew` on an axis; an ADC or RF event exceeds the interpreter's per-event sample limit | Lengthen the lobe, derate the design, or split the readout into several ADCs |
| `DefinitionConflict` | `name=` and `definitions['Name']` disagree | Pass one or the other |
| `CompilerContractError` | The compiled sequence does not match the tree: total duration, m0, m1 or a label address drifted; or a stage broke an IR contract | **A compiler bug.** Report it with the tree that produced it |
| `ConfigurationError` | `add()` got something that is not an event or a block; a unit is unknown; an `Opts` is unusable | Fix the call |
| `RasterError` | A time does not land on a raster it must | `sc.Raster(opts.grad_raster_time).ceil(t)` where the time is computed |

Every message follows one shape, so it is scannable and testable:

```text
HardwareLimitError: slew 189% of the 150 T/m/s limit on axis x.
  from   :  tr.readout.prephaser
  at     :  2.340 ms (block 117)
  reached:  283.5 T/m/s
  fix
    lengthen the lobe, or lower the readout bandwidth
    or design that part against sc.opts.derate(opts, slew=0.52)
```

**Catching them:**

```python
big = pp.make_trapezoid('x', amplitude=0.6 * opts.max_grad, duration=2e-3, system=opts)
try:
    sc.compile(sc.LogicBlock('tr').add(0.0, big, big), opts)
except sc.HardwareLimitError as err:
    assert str(err).startswith('gradient 120% of the 40 mT/m limit on axis x.')
```

All exception classes are re-exported at the package root — `sc.CompileError`,
`sc.HardwareLimitError`, and so on — regardless of which submodule defines them.

## 3.3 Warnings

Things the compiler **did** rather than refused. These are `SeqCraftWarning`, a `UserWarning`
subclass, emitted through the standard `warnings` machinery.

| Category | Meaning |
|---|---|
| `merge` | Two or more gradients on one axis were summed |
| `resample` | A gradient was resampled onto the raster to join a boundary continuously |
| `snap` | An RF or ADC reservation was snapped to the block raster |
| `orphan_label` | A label has no ADC after it, so its placement is undefined |
| `norm` | The vector norm across simultaneous axes exceeds the per-axis limit — normal for multi-axis gradients, and legal on real amplifiers |

**One aggregated warning per category, at the end of the compile.** Python's default filter shows a
warning once per source line, so one warning per merge would show the first and swallow the rest.
Instead you get a single line naming the count and the sites, with identical sites counted rather
than repeated:

```text
SeqCraftWarning: 8 same-axis gradient merges: tr.rewinder+tr.prephaser (axis x) x8
```

Handle them the standard way:

```python
import warnings

gentle = {'duration': 2e-3, 'rise_time': 400e-6, 'system': opts}
merged = sc.LogicBlock('m').add(0.0, pp.make_trapezoid('x', area=100.0, **gentle),
                                pp.make_trapezoid('x', area=200.0, **gentle))

with warnings.catch_warnings(record=True) as caught:      # inspect them
    warnings.simplefilter('always')
    sc.compile(merged, opts)
assert any('same-axis gradient merge' in str(w.message) for w in caught)

try:                                                       # or make them fatal
    with warnings.catch_warnings():
        warnings.simplefilter('error', sc.SeqCraftWarning)
        sc.compile(merged, opts)
    raise AssertionError('should have raised')
except sc.SeqCraftWarning:
    pass
```

In a test, `pytest.warns(sc.SeqCraftWarning, match='merge')` is the positive form; for the negative —
*no* warning of a kind — use `catch_warnings(record=True)`, because `pytest.warns` has no clean
negative.

---

# 4. `analysis` — measuring a tree

Four functions, one entry shape: **give it a tree, get numbers back.** You should not have to know
that PNS needs a compiled sequence while moments do not.

## 4.1 `sample(tree, opts) -> (grid, grads, marks)`

The tree on a uniform gradient raster.

```python
grid, grads, marks = sc.sample(spoil, opts)
# grid  : ndarray of times, seconds, spaced by opts.grad_raster_time
# grads : {'x': ndarray, ...} in Hz/m, summed per axis, only axes in use
# marks : [(kind, start, end, label), ...] for each RF, ADC and barrier
assert sorted(grads) == ['z'] and marks == []
assert bool(np.allclose(np.diff(grid), opts.grad_raster_time))
```

`marks` entries are `('rf' | 'adc' | 'barrier', start, end, label)`, in tree order; a barrier has
`start == end`.

**It samples the tree, not the compiled sequence** — deliberately. It shows what you *meant*, before
the compiler chose block boundaries, which is what you want when the question is "did I place this
correctly". Use `seq.plot()` to see what the compiler made of it.

Typical uses beyond plotting:

```python
slew = {ax: np.diff(g) / opts.grad_raster_time for ax, g in grads.items()}   # slew per axis
gaps = grid[np.flatnonzero(np.abs(grads['z']) < 1e-9)]                       # find dead time
```

> **`sample` is lossy on purpose.** The grid is uniform; an arbitrary gradient's samples sit at
> raster *centres* and an extended trapezoid's knots are not uniformly spaced at all, so both are
> **interpolated** onto the grid. For a moment, a split or a sum use `moments`, or `knots_of` +
> `pwl_moment`, which are exact. See §4.5.

## 4.2 `moments(tree, order=0) -> dict[str, float]`

Whole-tree gradient moment per axis, integrated from **exact knots**. Takes no `opts` and does no
compile — a moment is a property of the waveform, not of the scanner.

```python
assert round(sc.moments(spoil)['z'], 6) == 500.0            # m0, in 1/m
assert 'z' in sc.moments(spoil, 1)                          # m1, in s/m
```

| `order` | Quantity | Unit | Reads on |
|---:|---|---|---|
| 0 | area — k-space displacement | 1/m | Is my readout balanced? Did the spoiler survive? |
| 1 | first moment | s/m | Flow/motion sensitivity; whether a bipolar pair is truly nulled |
| 2 | second moment | s²/m | Acceleration sensitivity |

A worked check — is this diffusion pair nulled at the echo?

```python
m1 = sc.moments(tr, 1)
assert all(isinstance(v, float) for v in m1.values())
```

Because it does no compile, it works on a tree the compiler would refuse:

```python
off_raster = sc.LogicBlock('t').add(5e-6, gx)               # 5 us is off the gradient raster
assert round(sc.moments(off_raster)['x']) == 258            # still measurable
```

## 4.3 `kspace(tree, opts) -> dict[str, np.ndarray]`

The k-space trajectory, in 1/m. Compiles internally, then uses pypulseq's own calculation, so the
sample times are the **true ADC sample times** rather than a raster approximation.

```python
k = sc.kspace(tr, opts)
assert set(k) == {'k_adc', 't_adc', 'k', 't_k', 't_excitation', 't_refocusing'}
assert k['k_adc'].shape[0] == 3                             # (3, n_samples)
assert k['k_adc'].shape[1] == k['t_adc'].size
```

| Key | Meaning |
|---|---|
| `k_adc` | `(3, n_samples)` at the ADC sample times |
| `t_adc` | their times |
| `k`, `t_k` | the dense trajectory and its timebase |
| `t_excitation`, `t_refocusing` | RF centre times |

> **Why this exists rather than calling pypulseq directly:** `calculate_kspacePP` returns its tuple
> in a **different order** from `calculate_kspace`, and getting that wrong silently swaps the
> trajectory for its timebase — a wrong answer with no error. The named return makes that
> unmistakable.

## 4.4 `pns(tree, opts, hardware) -> dict[str, Any]`

Peripheral-nerve-stimulation prediction. Compiles internally, then delegates to pypulseq's SAFE
model implementation.

```python
r = sc.pns(spoil, opts, sc.hardware.synthetic_hardware())
assert set(r) == {'ok', 'peak', 'norm', 'components', 't'}
assert isinstance(r['ok'], bool) and r['peak'] >= 0.0
```

| Key | Meaning |
|---|---|
| `ok` | `True` if the model predicts the sequence stays under the stimulation threshold |
| `peak` | Peak stimulation as a fraction of the limit — `0.63` is 63 % |
| `norm` | The full time-resolved stimulation curve |
| `components` | Per-axis contributions |
| `t` | Timebase for `norm` and `components` |

**The full return matters.** When `ok` is `False`, `peak` tells you *how much* but `norm` and `t`
tell you *where*, which is what you need to fix it:

```python
if not r['ok']:
    worst_at = r['t'][int(np.argmax(r['norm']))]
    print(f'peak stimulation at {worst_at * 1e3:.1f} ms')
```

It **delegates** rather than reimplementing. `dG/dt` convolved with three exponentials is about
sixty lines and looks reachable from `sample`, but it is a *safety* calculation, pypulseq's
implementation is validated against vendor behaviour, and a second one that can silently drift is
the wrong thing to own.

> **Never used to clear a human scan with a synthetic model.** See §2.3.

## 4.5 Which function is exact, and which is not

The single most important thing to get right in this module:

| Function | Basis | Exact? |
|---|---|:--:|
| `sample` | uniform raster grid, **interpolated** | **No** — for looking, and for approximate numeric work |
| `moments` | `knots_of` + `pwl_moment` over `flatten(tree)` | **Yes** — never routed through `sample` |
| `kspace` | compiled, then `calculate_kspacePP()` | **Yes**, at true ADC sample times |
| `pns` | compiled, then `calculate_pns()` | pypulseq's validated SAFE model |

`moments` looks like it could be built on `sample` now that they sit together. It must not be — and
the reason is subtler than "sampling is lossy". Linear interpolation errs **antisymmetrically** about
each knot, so the halves cancel under the integral and **m0 comes back bit-identical** while the peak
is visibly rounded off:

```python
w = np.sin(np.pi * np.linspace(0.0, 1.0, 13)) * 0.1 * opts.max_grad
lobe = pp.make_arbitrary_grad('x', waveform=w, first=0.0, last=0.0, system=opts)
lobe_tree = sc.LogicBlock('lobe').add(0.0, lobe)

grid, grads, _ = sc.sample(lobe_tree, opts)
peak_exact = float(np.max(np.abs(sc.events.knots_of(lobe, 0.0)[1])))
peak_sampled = float(np.max(np.abs(grads['x'])))

assert (peak_exact - peak_sampled) / peak_exact > 0.01      # the peak is visibly lost
assert abs(float(sc.events.trapz(grads['x'], grid))
           - sc.moments(lobe_tree)['x']) < 1e-6             # the area is not
```

So a `moments` built on `sample` would look correct on every test anyone would think to write, and
the compiler's own self-check would be comparing two differently-wrong numbers.
`tests/analysis/test_analysis.py` asserts both halves.

---

# 4b. `modules` — the concrete building blocks

Everything above is the tree, the compiler and the contract, and none of it knows what a slice is.
`sc.modules` is the other half: MR physics with its arithmetic attached, one class per thing.

```python
gre = sc.modules.GRE2D(opts=opts, fov_mm=220.0, matrix=(64, 64), thickness_mm=5.0)
seq = sc.compile(gre(lines=range(64)), opts, name='gre_2d')
```

Re-exported **flat**, so no import path names a folder:

| Name | What |
|---|---|
| `Excitation` | an RF pulse and, when selective, its selection gradient and rephaser |
| `PhaseEncode` | one Cartesian phase-encode blip, designed once and scaled per line |
| `CartesianLine` | prephaser, readout gradient and ADC as one design |
| `spoiler` | a gradient winding *n* turns of phase across a voxel — a **function**, not a class |
| `GRE2DTR` | one repetition of a spoiled 2D gradient echo |
| `GRE2D` | the complete scan |

The folders behind that table are taxonomy, and each has a rule: `rf/` is `rf.use`, `encoding/` is
gradients with no ADC, `readout/` contains an ADC, `kernel/` composes more than one leaf folder,
`imaging/` composes kernels, and the top level holds what is not an `sc.Module` subclass.
`tests/modules/test_layout.py` asserts all of it.

**Nothing here was designed in the abstract.** Every module was extracted from
[`examples/gre_2d/01_build.ipynb`](../examples/gre_2d/README.md), which builds the same sequence
out of raw pypulseq beside it — so a module that cannot be extracted without altering the sequence
is not a module, and one whose extraction does not shorten the notebook is a wrapper.

## What each one knows that a block cannot

`LogicBlock.duration` is measured, so a module never declares one. What a module *does* answer is
the questions a tree of events cannot:

| | |
|---|---|
| `Excitation.time_to_center()` | to the RF's effective centre — for a minimum-phase pulse, nowhere near the midpoint |
| `Excitation.time_to_rephaser()` | where the slice rephaser begins, so another axis can start there |
| `CartesianLine.time_to_echo()` | to k = 0, which is not the middle of the ADC window |
| `GRE2DTR.min_te_s` / `min_tr_s` | feasibility, known at design time; a shorter request raises |

## The two couplings

**The winder** is three gradients on three axes playing at once — the slice rephaser on `z`, the
phase-encode blip on `y`, the readout prephaser on `x`. Each leaf reports its minimum and accepts
an override; `GRE2DTR` takes the maximum and passes it down. Stretching the short ones keeps TE at
its minimum, which inserting a delay would not.

**The receiver is phase-locked to the transmitter.** `GRE2DTR` gives the same `phase_deg` to the
excitation and to the readout. Moving one and not the other writes the RF-spoiling schedule's
quadratic phase into `ky`, which scatters a point source across the phase-encode direction while
the readout direction stays perfectly correct.

## The sampling pattern is a build argument

```python
gre(lines=range(ny))                               # fully sampled
gre(lines=sorted(acs | set(range(0, ny, 3))))      # undersampled with a calibration block
gre(lines=reversed(range(ny)))                     # ordering, with no new argument
```

One argument replaces acceleration, phase-encode partial Fourier, multi-shot and ordering, because
each of those is a different list — and **no generators ship**, because which lines to acquire is a
sequence-programming choice ([ADR-003](adr/003-scanner-and-module-reform.md)). Out-of-range or
duplicate indices raise; a pattern omitting the centre of k-space warns.

---

# 5. `display` — looking at a tree

## `plot_block(root, opts, *, title='', figsize=(10.0, 4.5))`

The **one** picture pypulseq cannot draw: the tree, before the compiler chose boundaries.

```python
fig = sc.plot_block(tr, opts, title='readout')       # needs matplotlib
fig.savefig('readout.png')
```

Gradients per axis, with RF and ADC windows shaded and barriers marked. Returns the figure and never
calls `show()`, so the caller decides — a notebook displays it, a script saves it, a test discards it.

**This is the only module allowed to import matplotlib,** and it does so lazily inside each function,
so `import seqcraft` stays cheap. Without matplotlib installed you get `MissingExtraError` naming the
extra: `pip install "seqcraft[viz]"`. A test asserts that no other file in the package imports it.

For a **compiled** sequence, use pypulseq's own plotter — it is better maintained and it is the same
picture everyone else in the ecosystem reads:

```python
seq = sc.compile(tr, opts)
seq.plot()                     # block boundaries and all
```

The difference between the two pictures is exactly the block structure the compiler chose.

---

# 6. `errors` — the exception hierarchy

```text
SeqCraftError                  # the base; catch this to catch everything
├── ConfigurationError         # a call is wrong: bad type, unknown unit, unusable Opts
│   ├── RasterError            #   a time is off a raster            (design/timing.py)
│   └── UnknownFieldError      #   an Opts override that is not a field (scanner/opts.py)
├── CompileError               # this tree cannot become a legal sequence (compiler/errors.py)
├── HardwareLimitError         # the machine cannot play it          (compiler/errors.py)
├── DefinitionConflict         # two sources claim one [DEFINITIONS] key (compiler/errors.py)
└── MissingExtraError          # an optional dependency is needed; names the extra

CompilerContractError(RuntimeError)   # a compiler bug (compiler/verification.py)
SeqCraftWarning(UserWarning)          # the compiler changed something
```

**Exceptions live with the code that raises them**, and only those raised from more than one package
stay at the root. Every one is re-exported as `sc.<Name>`, so the spelling never depends on the
layout — and the layering test asserts the *identity*, not just the presence:

```python
import seqcraft.compiler.errors, seqcraft.design.timing, seqcraft.scanner.opts
assert sc.CompileError is seqcraft.compiler.errors.CompileError
assert sc.RasterError is seqcraft.design.timing.RasterError
assert sc.UnknownFieldError is seqcraft.scanner.opts.UnknownFieldError
assert issubclass(sc.HardwareLimitError, sc.SeqCraftError)
```

## `format_error(headline, fields=None, fixes=(), *, sections=()) -> str`

The shared formatter. Use it if you raise from your own component and want the same shape.

```python
from seqcraft.errors import format_error
print(format_error('TE = 4.00 ms is too short.', {'minimum': '6.34 ms'}, ["te_s='min'"]))
```

```text
TE = 4.00 ms is too short.
  minimum:  6.34 ms
  fix
    te_s='min'
```

`sections` adds extra named blocks after `fields`, for showing one ledger per constraint chain.

---

# 7. Compiler internals

You do not need these to use seqcraft. They are documented because the package asserts a layering
over them, because a `CompilerContractError` points here, and because the error messages name them.

**None of these is re-exported.** `seqcraft.compiler.__all__` is `['compile_sequence']` and nothing
else; reach a stage module by its full path.

| Stage | File | Responsibility |
|---|---|---|
| **place** | `placement.py` | Walk the tree, resolve every event to absolute active and reservation intervals |
| **bound** | `boundaries.py` | Choose block edges: exclusivity, indivisibility, label retiming, gap feasibility |
| **legalize** | `legalization.py` | Schedule events, transform gradients, check limits, and return validated ready blocks |
| **emit** | `emission.py` | Mechanically hand validated ready blocks to `Sequence.add_block` |
| **verify** | `verification.py` | IR contracts, finished-sequence checks, and the against-the-tree self-check |

## 7.1 `compiler/model.py` — the IR and the policy

```python
@dataclass(frozen=True)
class PlacedEvent:
    node_t: float            # the node's own time, before the event's delay
    start: float; end: float             # the active interval
    res_start: float; res_end: float     # the reservation, including dead times
    event: Any                           # a *reference* to the pypulseq event, never a copy
    path: tuple[str, ...]                # provenance
```

Properties: `kind` (the event's `type`), `where` (the dotted path), `source_path`, `duration`,
`reservation_duration`, and `summary()`. It holds a *reference* rather than a copy: seqcraft treats
pypulseq events as read-only during a compile.

```python
@dataclass(frozen=True)
class PulseqReadyBlock:
    index: int
    start: float; end: float; duration: float
    events: tuple[Any, ...]
    source_paths: tuple[tuple[str, ...], ...]
    origin: tuple[str, ...]
```

Properties: `kinds`, and `summary()`.

```python
@dataclass(frozen=True)
class LegalizationResult:
    blocks: tuple[PulseqReadyBlock, ...]
    notes: tuple[tuple[str, tuple[str, ...]], ...]
```

`notes` is explicit immutable stage output: emission neither discovers nor mutates warnings.

| Name | What it is |
|---|---|
| `EXCLUSIVE_KINDS` | `{'rf', 'adc'}` — a block holds at most one, so a boundary is *required* between two |
| `INDIVISIBLE_KINDS` | `{'rf', 'adc', 'trigger', 'output'}` — no boundary may fall strictly inside |
| `in_block_delay(p, block_start, opts)` | the event's delay within its block, snapped to its own raster in integer ticks |
| `time_equal`, `time_before`, `time_at_or_before`, `time_strictly_between` | comparisons at `EPS` tolerance |
| `interval_duration(start, end)` | exact tick arithmetic, not `end - start` |

Comparisons use tolerance; interval arithmetic uses integer ticks. Keeping those apart is
[ADR-001](adr/001-compiler-internal-time-policy.md).

## 7.2 `compiler/placement.py`

- `place_events(root, opts) -> tuple[PlacedEvent, ...]` — the tree to absolute-time IR.
- `UNSUPPORTED_KINDS` — `{type: (description, hints)}` for pypulseq events seqcraft rejects **by
  name** rather than dropping: `rot3D`, `soft_delay`, `rf_shim`. A tree carrying a rotation extension
  once compiled clean, reported nothing, and played unrotated.

## 7.3 `compiler/boundaries.py`

| Function | Does |
|---|---|
| `check_exclusive(placed) -> None` | Raise if two RF or ADC reservations overlap — the one conflict no choice of boundaries can fix |
| `find_boundaries(placed, total, raster, max_block) -> list[float]` | The sorted block edges |
| `label_targets(placed) -> dict[int, float]` | Per label, the time it should be *assigned* at — its target ADC's reservation start |
| `label_order_conflict(placed, targets) -> CompileError \| None` | Reject label groups pulseq cannot express, because it sorts a block's extensions by library id |
| `orphan_label_notes(placed, targets) -> list[str]` | One note per label with no ADC after it |

Labels are assigned by **target**, not by containment: a boundary pushed past an ADC's reservation
end would otherwise put a label into the previous readout's block and silently overwrite its k-space
address.

## 7.4 `compiler/legalization.py`

| Function | Does |
|---|---|
| `superpose(pieces, a, b) -> (rel_ticks, amps)` | The exact PWL superposition over `[a, b]`, in integer ticks |
| `axis_gradient(axis, pieces, a, b, opts, notes, block_index)` | The single gradient event for one axis over one interval, or `None` |
| `check_limits(events, opts, block_index, origin, notes, *, start=0.0) -> None` | Raise `HardwareLimitError` on a per-axis violation; append norm findings to `notes` |
| `required_duration(events, opts) -> float` | The shortest block that can hold `events`, by PyPulseq's rules |
| `common_path(paths) -> tuple[str, ...]` | The longest tag path common to everything in a block |
| `legalize_blocks(edges, placed, targets, opts, raster) -> LegalizationResult` | Produce and validate the complete immutable ready-block sequence and its notes |

The sum of PWL functions is PWL with knots at the union of theirs, so both operations — split at a
boundary, sum on an axis — are exact. `axis_gradient` picks whichever pulseq representation holds the
knots untouched (extended trapezoid for on-raster knots, arbitrary gradient for the raster-centre
pattern). Only a waveform bending both on and off the raster fits neither, and that one is resampled
with a **measured** bound on how far it moved.

## 7.5 `compiler/emission.py`

| Function | Does |
|---|---|
| `emit_blocks(seq, blocks) -> None` | Add validated ready blocks without scheduling, transforming, or checking limits |

The ready-block origin is used to add source context when PyPulseq rejects a block. Emission has no
access to placed events, boundaries, label targets, the scanner limits, or transformation notes.

## 7.6 `compiler/verification.py`

| Function | Does |
|---|---|
| `verify_placed_events(events) -> tuple[ContractViolation, ...]` | Structural contract on the first IR |
| `verify_ready_blocks(blocks, *, expected_first_index=0, expected_start=None)` | Structural contract on the second |
| `require_valid_contract(name, violations) -> None` | Raise `CompilerContractError` if any |
| `check_event_sizes(seq, opts, origins=()) -> None` | Raise `HardwareLimitError` above the interpreter's per-event sample limit |
| `check_label_addresses(seq) -> None` | Raise `CompileError` if two imaging ADCs write the same k-space address |
| `expected_addresses(placed, targets) -> list[dict[str, int]]` | Fold the tree's labels the way the interpreter will |
| `verify_against_tree(placed, targets, *, duration_s, tree_duration_s, moments, label_states) -> None` | The four semantic invariants |

**The self-check** asserts four invariants after emission, each catching a class of compiler bug no
individual test case would:

1. **Total duration** matches the tree's — guards a boundary merge that dropped time.
2. **m0 per axis** matches the sum over flattened events — a split that lost a tail changes it.
3. **m1 per axis**, referenced to the sequence start — m0 survives a *time shift*, so it cannot see a
   gradient that plays at the wrong moment; m1 can, because area `A` displaced by `dt` changes it by
   `A·dt`.
4. **Label addresses** match the fold of the tree's labels — a duplicate-address check only fires on
   a *collision*, but an addressing shifted by one and still unique passes it.

All four are measured before anything is raised, because one mistimed stage usually breaks several
axes and the pattern is the diagnosis.

Tolerances scale with the total, because they must: float64 resolves about 4 ns at 20 000 s, and the
moment tolerance is scaled by total area *traversed* rather than the net, since a readout and its
prephaser nearly cancel.

`_sequence_moments(seq, order)` is the compiled-side integrator. It is private and deliberately
**not** shared with `analysis.moments`, which walks the tree: a self-check that calls the same code
as the thing it checks compares a number with itself.

`ContractViolation(location, message)` is the finding type; `CompilerContractError(name, violations)`
is the exception.

## 7.7 `_compat`

`require()` probes the pinned pypulseq fork for the functions seqcraft needs and fails once with a
complete list, rather than letting the first caller that needs a missing one fail with an opaque
`AttributeError` halfway through a build. `PYPULSEQ_VERSION` is the expected version string. It runs
at import.

---

# 8. Migration

| Removed | Replacement |
|---|---|
| `out = sc.compile(...)` returning `CompiledSequence` | `seq = sc.compile(...)` returning `pypulseq.Sequence` |
| `out.seq` | `seq` |
| `out.n_blocks` | `len(seq.block_events)` |
| `out.duration_s` | `seq.duration()[0]` |
| `out.definitions` | `seq.definitions` — set by the compiler |
| `out.check()`, `out.report`, `Issue`, `Report`, `ReportFailed` | the compile raised, or warned |
| `out.report.of_kind('grad_resample')` | `pytest.warns(sc.SeqCraftWarning, match='resample')` |
| `out.check().raise_if_failed()` | nothing — failures already raised |
| `out.write(path)`, `WriteResult` | `seq.write(str(path))` |
| `out.write(path, sidecar=True)`, `sc.provenance` | nothing yet — see [`serialization.md`](serialization.md) |
| `out.origin(i)` | gone; provenance paths still name the source in every error message |
| `out.moments(1)` | `sc.moments(tree, 1)` |
| `out.kspace()` | `sc.kspace(tree, opts)` |
| `out.pns(hw)` | `sc.pns(tree, opts, hw)` — now returns the full curve |
| `sc.plot_sequence(out)` | `seq.plot()` |
| `sc.plot_trajectory(interleaves)` | `sc.kspace(...)` + three lines of matplotlib |
| `sc.testing.assert_*` | the compiler raises; for purity see [`writing_a_module.md`](writing_a_module.md) |
| `sc.Geometry`, `compile(geometry=)` | `salvage/geometry.py`, passed via `definitions=` |
| `sc.validate.*`, `merge_definitions` | gone; two sources need a collision check, not a merge |
| `sc.UnitSanityError` | raised nowhere; `salvage/geometry.py` has its own |
| `seqcraft.design.sampling` | `seqcraft.analysis` |
| `seqcraft.result`, `seqcraft.report`, `seqcraft.testing` | deleted |

[ADR-004](adr/004-compile-returns-a-sequence.md) records why.

---

# 9. Index

| Name | Where | Kind |
|---|---|---|
| `ADDRESS_KEYS` | `design.events` | constant |
| `ASC_ENV_VAR` | `scanner.hardware` | constant |
| `AXES` | `design.events` | constant |
| `BARRIER` | `design.logic` | constant |
| `CartesianLine` | `modules` | class |
| `CompileError` | `compiler.errors` | exception |
| `CompilerContractError` | `compiler.verification` | exception |
| `ConfigurationError` | `errors` | exception |
| `ContractViolation` | `compiler.verification` | class |
| `DefinitionConflict` | `compiler.errors` | exception |
| `EPS` | `design.timing` | constant |
| `EXCLUSIVE_KINDS` | `compiler.model` | constant |
| `Event` | `design.events` | type alias |
| `Excitation` | `modules` | class |
| `GAMMA_1H` | `design.units` | constant |
| `GRADIENT_KINDS` | `design.events` | constant |
| `GRE2D` | `modules` | class |
| `GRE2DTR` | `modules` | class |
| `HANDLED_KINDS` | `design.events` | constant |
| `HardwareLimitError` | `compiler.errors` | exception |
| `INDIVISIBLE_KINDS` | `compiler.model` | constant |
| `Item` | `design.logic` | type alias |
| `LABEL_KINDS` | `design.events` | constant |
| `LegalizationResult` | `compiler.model` | class |
| `LogicBlock` | `design.logic` | class |
| `MissingExtraError` | `errors` | exception |
| `Module` | `design.module` | class |
| `Node` | `design.logic` | class |
| `POINT_KINDS` | `design.events` | constant |
| `PYPULSEQ_VERSION` | `_compat` | constant |
| `PhaseEncode` | `modules` | class |
| `PlacedEvent` | `compiler.model` | class |
| `PulseqReadyBlock` | `compiler.model` | class |
| `Raster` | `design.timing` | class |
| `RasterError` | `design.timing` | exception |
| `SeqCraftError` | `errors` | exception |
| `SeqCraftWarning` | `errors` | warning |
| `TICKS_PER_SECOND` | `design.timing` | constant |
| `UNSUPPORTED_KINDS` | `compiler.placement` | constant |
| `UnknownFieldError` | `scanner.opts` | exception |
| `axis_gradient` | `compiler.legalization` | function |
| `barrier` | `design.logic` | function |
| `check_event_sizes` | `compiler.verification` | function |
| `check_exclusive` | `compiler.boundaries` | function |
| `check_label_addresses` | `compiler.verification` | function |
| `check_limits` | `design.events` | function |
| `check_limits` | `compiler.legalization` | function |
| `common_path` | `compiler.legalization` | function |
| `compile_sequence` | `compiler` | function |
| `content_hash` | `design.events` | function |
| `convert` | `design.units` | function |
| `derate` | `scanner.opts` | function |
| `derive` | `design.events` | function |
| `dimension_of` | `design.units` | function |
| `dimensions` | `design.units` | function |
| `emit_blocks` | `compiler.emission` | function |
| `exact_diff` | `design.timing` | function |
| `exact_sum` | `design.timing` | function |
| `expected_addresses` | `compiler.verification` | function |
| `find_boundaries` | `compiler.boundaries` | function |
| `flatten` | `design.logic` | function |
| `format_error` | `errors` | function |
| `from_scanner` | `scanner.opts` | function |
| `from_ticks` | `design.timing` | function |
| `in_block_delay` | `compiler.model` | function |
| `interval_duration` | `compiler.model` | function |
| `knots_of` | `design.events` | function |
| `known_units` | `design.units` | function |
| `kspace` | `analysis` | function |
| `label_order_conflict` | `compiler.boundaries` | function |
| `label_targets` | `compiler.boundaries` | function |
| `legalize_blocks` | `compiler.legalization` | function |
| `load_hardware` | `scanner.hardware` | function |
| `moments` | `analysis` | function |
| `orphan_label_notes` | `compiler.boundaries` | function |
| `place_events` | `compiler.placement` | function |
| `plot_block` | `display` | function |
| `pns` | `analysis` | function |
| `pwl_moment` | `design.events` | function |
| `require` | `_compat` | function |
| `require_valid_contract` | `compiler.verification` | function |
| `required_duration` | `compiler.legalization` | function |
| `sample` | `analysis` | function |
| `spoiler` | `modules` | function |
| `span` | `design.logic` | function |
| `superpose` | `compiler.legalization` | function |
| `synthetic_hardware` | `scanner.hardware` | function |
| `time_at_or_before` | `compiler.model` | function |
| `time_before` | `compiler.model` | function |
| `time_equal` | `compiler.model` | function |
| `time_strictly_between` | `compiler.model` | function |
| `to_ticks` | `design.timing` | function |
| `verify_against_tree` | `compiler.verification` | function |
| `verify_placed_events` | `compiler.verification` | function |
| `verify_ready_blocks` | `compiler.verification` | function |
| `waveform_of` | `design.events` | function |
