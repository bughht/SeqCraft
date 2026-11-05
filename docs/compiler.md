# The compiler

A logic block lets anything overlap anything. A pulseq block holds **at most one RF, one ADC and
one gradient per axis**, must last a whole number of block rasters, must leave the RF dead time
before a pulse and the ringdown after it, and must have gradients that join continuously across its
boundaries.

The compiler's job is to find block boundaries satisfying all of that, and to combine whatever lands
inside each one. It is the only complicated subsystem in seqcraft, which is deliberate: one hard
thing, one file per stage, tested at every responsibility boundary.

The current stage architecture is stable. For ownership and review steps when adding an event or
changing boundary policy, use [`compiler_maintenance.md`](compiler_maintenance.md); changes to the
architecture itself require the evidence described in the
[freeze record](refactor/compiler_architecture_freeze.md).

---

## How an event's own `delay` is handled

This is the first question people ask, and it is worth answering precisely, because pypulseq puts
meaningful physics in that field.

A slice-selective RF pulse is the clearest case. `make_sinc_pulse` returns an `rf` whose `delay` is
`max(rf_dead_time, gz.rise_time)` — 130 µs for a 5 mm slice on a 40 mT/m system — and a `gz` whose
delay is zero. The delay is not padding. It is there so the slice-select gradient can **reach its
plateau before the pulse starts**, and so the transmitter's dead time is respected.

seqcraft treats the delay as **part of the event**, never folds it away, and reserves it:

```
node time (what you passed to add)          t
event's own delay                           + rf.delay
------------------------------------------------------
the pulse actually starts at                t + rf.delay
```

and the *reservation* — the interval no block boundary may fall inside — starts at **`t`**, not at
the pulse. So the dead time is inside the protected span:

```python
lb = sc.LogicBlock('exc').add(5e-3, rf, gz)      # placed at 5 ms
```

| | node time | active from | active to | reserves |
|---|---|---|---|---|
| `rf` | 5000 µs | 5130 µs | 6130 µs | 5000 → 6160 µs |
| `gz` | 5000 µs | 5000 µs | 6260 µs | 5000 → 6260 µs |

The RF's reservation runs from 5000 µs (its node time, which already contains the dead time) to
6160 µs (through the 30 µs ringdown). The compiled block then comes out as:

```
block 1:  5000.0 us   (empty -- the gap before the pulse)
block 2:  1260.0 us   rf.delay=130.0   gz.delay=0.0
```

**`rf.delay` is still exactly 130 µs inside the block.** It was preserved, not recomputed, because
the boundary landed at or before the node time — which is guaranteed by the reservation.

Three consequences worth knowing:

- **You place the pulse where you want it to *be*, not where its waveform starts.** `add(5e-3, rf)`
  means "this pulse belongs at 5 ms"; the 130 µs of dead time and gradient ramp is the pulse's own
  business.
- **A conflict inside a dead time is caught.** An RF whose waveform ends 5 µs before an ADC starts
  looks fine and is not: the ringdown and the ADC dead time overlap. Because the comparison is
  between *reservations*, that is an error naming both, rather than something pypulseq rejects
  40 000 blocks later.
- **In-block delays are quantised onto the event's own raster** — 1 µs for RF, 100 ns for ADC, 10 µs
  the rasters `Opts` carries — and computed in integer ticks
  (`design.timing`). A plain subtraction of absolute times drifts:
  at 39 s into a sequence it produced an RF delay of `129.9999999986 µs`, which pypulseq rejects.

An ADC's `delay` works the same way: `make_adc` sets it to at least `adc_dead_time`, the reservation
starts at the node time, and the in-block delay comes out unchanged.

A `delay` event is the one exception. pypulseq stores its length in `delay` and it has no waveform,
so it occupies `[t, t + delay]` rather than starting after its own delay. That is what makes
`lb.add(0.0, pp.make_delay(4.2e-3))` the way to hold a block open — a b=0 diffusion volume that must
fill the same slot as an encoded one.

---

## Two properties decide where a boundary may and must go

Keeping these apart is what makes triggers work, and it is the difference between "cannot be
split" and "cannot share a block".

**Indivisible** — no boundary may fall strictly inside the span. An RF's dead time and ringdown, an
ADC's window and trailing dead time, a trigger's pulse. Splitting any of them is meaningless: each
is one hardware action, not a waveform.

**Exclusive** — a block may hold at most one, so a boundary is *required* between two. RF and ADC
only, because pulseq stores one of each per block.

| Event | Indivisible | Exclusive |
|---|---|---|
| `rf`, `adc` | ✅ | ✅ |
| `trigger`, `output` | ✅ | ❌ — pulseq accepts several per block |
| `trap`, `grad` | ❌ — split exactly, then summed | ❌ — summed exactly |
| `labelset`, `labelinc` | — | — retimed onto their target ADC's block |

Treating a trigger as exclusive would invent a constraint the hardware does not have. Treating it as
divisible — which is what happened before this distinction existed — let a boundary land inside one,
and the block then came out shorter than the trigger it had to hold.

---

## How boundaries are chosen

A boundary is only placed where it **cuts nothing**. Three rules produce them.

**Mandatory.** One boundary somewhere in the gap between each pair of consecutive *exclusive*
reservations. That alone guarantees at most one RF and one ADC per block — pulseq's rule, satisfied
by construction rather than by checking. Zero and the total duration are boundaries too, as are
explicit `sc.barrier()` markers. If nothing in a gap can legally be cut — a trigger stretched right
across it — that is genuinely impossible in pulseq, and the error names the trigger.

**Opportunistic.** Every gradient edge and reservation edge is a *candidate*, accepted only if it
falls strictly inside neither a gradient nor an indivisible span. This is what keeps a trapezoid a
trapezoid: a slice rephaser on z overlapping a phase blip on y stays two clean events in one block,
rather than being cut where the blip happens to end. It is also why a readout gradient survives
whole — its own edges are candidates, and nothing else's edge lands inside it.

**Forced.** Intervals longer than pulseq's fixed-width duration field are subdivided.

A `sc.barrier()` inside an indivisible span is an error rather than a silent later failure: it is an
explicit request for a boundary somewhere one cannot exist.

Snapping is directional: **starts round down, ends round up**, so a block only ever grows to fit its
contents. Nearest-rounding an end is a real bug — an ADC whose window plus trailing dead time lands
4 µs past a raster edge would get a block 6 µs too short, which pypulseq reports as an unaligned
duration rather than as the missing microseconds it is.

---

## Three rules on overlap

| What overlaps | What happens | Why |
|---|---|---|
| Gradients on **different axes** | Nothing at all | A phase blip beside a rephaser beside a prephaser is the normal way to build a sequence. Warning about it would only teach you to ignore warnings. |
| Gradients on the **same axis** | **Warning**, then summed | Summing is almost always what was meant. But it is the one place the compiler changes a waveform, so it says so, naming both sources. |
| Two **RF**, two **ADC**, or RF and ADC | **Error** | You cannot transmit twice at once, sample twice at once, or transmit and receive at once. No choice of boundaries fixes that. |

---

## Labels are placed by target, not by containment

A pulseq label is a **running register**, not an instant: the interpreter applies a block's labels on
reaching that block, and an ADC in the same block then samples the state. So a label may live in any
block strictly after the previous ADC's and at or before its target ADC's — all equivalent.

The compiler therefore has to *choose*, and choosing by containment was wrong. A boundary pushed
later than an ADC's reservation end put a label placed comfortably after one readout into that
readout's own block, silently overwriting its k-space address. Each label is instead attached to the
block holding the **first ADC at or after its own time**, which makes the result independent of where
boundaries land.

Two consequences worth knowing:

- **A barrier is never crossed.** A barrier is an explicit request for a seam, so a label is not
  retimed over one. This is also what makes a barrier a working way to order two labels.
- **Intra-block label order carries no meaning.** pypulseq sorts a block's extensions by library id,
  so `add_block(set, inc)` and `add_block(inc, set)` build the identical block. Only
  order-independent groups are expressible — labels on different keys, or several `labelinc` on one
  key, which commute. Anything else is an error naming the barrier remedy, because guessing would
  silently pick one of two different addressings.

A label with no ADC after it addresses nothing, so it keeps containment placement and is **reported
as a warning** — it can still land in the preceding readout's block.

---

## Limits are checked after summing

Two individually legal gradients on one axis can sum to an illegal one, and no module can see that
in isolation. On a 40 mT/m, 150 T/m/s system, adding an area-100 and an area-200 trapezoid reaches
93 % of the amplitude limit and **189 % of the slew limit** — each one perfectly legal alone.

So amplitude and slew are measured on the *compiled* waveform, which is the only place the truth is
visible. That is also why the per-module `validate()` of an earlier design was checking the wrong
thing.

The **vector norm** across simultaneous axes is a warning, not an error: two axes ramping together
reach √2 times the per-axis slew in vector magnitude and three reach √3, which real amplifiers
permit. A three-axis spoiler routinely reports 139 % of the amplitude limit on the norm and is fine.

**The vector-norm warning is not a proxy for peripheral nerve stimulation.** They answer different
questions — one about an instant, the other about a history, and a waveform can sit inside every
instantaneous limit while accumulating a PNS response several times over threshold. Measure it with
`sc.pns(tree, opts, hw)` against a real hardware model rather than reading it off the warnings.

**And `synthetic_hardware()` cannot give you a verdict.** It exists so CI always has *a* model, and
it is an *illustrative* one — pypulseq's public `safe_example_hw()` — not a bound and not
representative of any particular scanner. An illustrative response model may differ from site
hardware in the absolute peak, in which axis dominates, and in how it ranks one waveform against
another, so a number from it is a smoke test that the calculation ran, not a verdict on whether a
sequence is runnable. Load the vendor `.asc` with `sc.hardware.load_hardware()` and
the path to it for anything scanner-specific; the file stays outside the repository, and only its
name and sha256 are recorded.

Which term dominates is not guessable and it moves as the sequence changes, so it cannot be
reasoned about, only measured — and measured against the model you actually care about. Two things
that catch people out, and both are properties of the response rather than of any one model:

**PNS depends on the recent slew history, not on elapsed time.** Changing how long a sequence
takes, or changing one waveform component, therefore does not predictably move the peak — the
realised waveform has to be re-evaluated rather than reasoned about.

**Attribution is not additive.** Changing one gradient can move the peak onto a different one, so
the gradient that limits you changes as you tune. `origin(block_index)` at the peak time is how you
find out which one it is.

One further caveat: a `.asc` also declares forbidden acoustic resonance bands, and nothing in seqcraft
checks them. A spiral is a narrowband gradient waveform and can sit squarely inside one.

---

## Four invariants that run on every compile

**Total duration** must equal the tree's. The tolerance scales with the total, because it has to:
float64 resolves about 4 ns at 20 000 s, so demanding nanosecond agreement on a long acquisition
would flag every sequence over an hour.

**m0 and m1 per compiled block and axis**, each about the block's own start. m0 is what catches a
split that lost a tail, a merge that dropped a piece, or area moved from one block into the next —
a whole class of compiler bug that no individual test case would have to think of. m1 catches a
gradient at the wrong moment, which m0 cannot see: m0 is exactly what a time shift preserves, and a
piece of area `A` displaced by `dt` moves m1 by `A·dt`.

**Why per block.** Referenced to the start of the sequence, m1 moved with placement. Registration
moves a stored shape's area by `δm0`, so the same legal waveform placed at `t` moved the absolute
m1 by `t·δm0`, and a compressible shape accepted at 0 s was refused at 1 s. A tolerance wide enough
for that at 10 s would hide a 4.5 1/m blip played one raster late. About each block's own start the
lever arm is at most a block long, and a correct waveform compares the same wherever it sits. The
sequence's whole moments are sums of the block ones — `M0 = Σ m0_k` and `M1 = Σ (m1_k + a_k·m0_k)`,
with `a_k` the block start — so a whole-sequence comparison that tolerated every block's error could
never fire, and a narrower one would bring `a_k·δm0_k` back. There is none. The decision is
[ADR-007](adr/007-against-tree-moments-per-block.md).

Both are compared in **two stages**, each against its own input:

| stage | compares | its tolerance covers |
|---|---|---|
| compile | the tree's placed events, restricted exactly to each block, with the ready blocks emission hands PyPulseq | SeqCraft's own arithmetic: float summation; the float of superposing pieces; the endpoint snap; and where a knot's time lands in ticks |
| registration | the ready blocks, measured before emission, with what PyPulseq stored, read back block by block | PyPulseq's storage precision: compression rounding, and the six-digit key it deduplicates events by |

A refusal names its stage, block, time, difference and tolerance, so it says whether SeqCraft or
PyPulseq moved the waveform. Stored against tree is bounded by the two tolerances together.

**Resampled blocks.** Where legalization resampled a block's gradient onto the raster, it says so
as data, not only in the warning's text. The compile stage then compares that block and axis with
the tree resampled the same way, onto the gradient raster from the block's start. Neither of the
alternatives is sound: a resample moves area between the blocks of a stretch, so the tree as given
refuses it; and merging a stretch's inner seams is exact only where both its ends are flat, since a
sloped end leaves a boundary term of about `r²·g′/8` in the stretch. The other axes of a resampled
block are compared as usual, and registration needs no special case: what it compares is the
resampled waveform itself.

**Each tolerance, and where it holds.** Every term is a stated bound, not a constant tuned to a
case.

- *Float*, both stages. A sum of `N` segment terms is within `(N + 10)·u·Σ|term|` of exact, with
  `u = 2⁻⁵³`, for any `N` and any summation order.
- *Superposition*, compile. Each superposed value is within `(k + 8)·u` of the summed end values
  of the segments its `k` pieces are interpolated from. `superpose` interpolates in ticks from the
  block's start, so this does not depend on where the block sits. A raster resample interpolates
  in float seconds from the block's start, which adds `4·u·L` times the summed peak slope for a
  block of length `L`.
- *Endpoint snap*, compile. `axis_gradient` sets a superposed block's first or last value to zero
  when it is below `1e-6·max_grad`, removing a triangle of that height out to the next knot. This is
  modelled from the tree's own value at each block edge, wherever the compile's value could fall
  under the threshold. The snap stays because PyPulseq's timing check refuses an end above 1e-9 Hz/m
  that the next block does not continue: it is what keeps a summed end that cancels to nearly zero
  from being refused.
- *Time representation*, compile. The raster check puts a gradient's start within half a tick of
  its tick, and superposing turns each knot's absolute float time into ticks with a few more
  roundings of about `u·t`. A knot can therefore land a tick (1 ps) from where the other side puts
  it only when its piece's start lies within those roundings of half a tick, or, for every piece in
  a block, once `u·t` itself approaches a tick, past about 1100 s. The term is that shift times
  each such piece's variation and end values, and zero for every other piece.
- *Compression*, registration. PyPulseq stores a shape as differences of its peak-normalised
  samples, each rounded to `1e-7`. A stored sample is within half a quantum, except when the first
  sample lies exactly on a half quantum: rounding half away from zero then carries one more quantum
  into every later sample. The bound is the edge case, `1.5e-7` of the peak, over the shape's whole
  span. It needs no condition on duty, which is what closed the low-duty gap of the earlier
  traversed-area tolerance.
- *Deduplication*, registration only. Every PyPulseq library keys an entry by its values printed
  with `%.6g`, and stores an event that prints like an earlier one as that one. A value then moves
  by under `1e-5` of itself; a shape difference of at least a tenth of the peak, which prints more
  coarsely than a quantum, moves by under `1e-5` of itself, and that carries into every later sample.
  Two legal trapezoids at 100000.0 and 100000.4 Hz/m are therefore accepted — and the second still
  plays 100000.0. This allowance never enters the compile stage.
- *Domain.* Event time parameters and shape lengths below `1e6` gradient rasters, which `%.6g`
  prints exactly. The tree's gradients lie inside the sequence and each block's gradient inside its
  block; content outside either is reported in its own right.

**What these checks cannot detect.** m0 and m1 are two linear functionals per block. They are not
waveform fidelity.

- A change that keeps every block's area and centroid: a distortion orthogonal to `1` and `t − a`,
  a zero-net piece moved inside its block, or two equal pieces swapped.
- Anything below its tolerance. In the registration stage that includes any change within
  PyPulseq's six-digit precision, about `1e-5` of a block's `|g|` integral.
- A resample the compiler reported: it is accepted as the raster interpolant of the tree, which is
  what the `resample` warning states.
- RF and ADC timing, amplitude and slew limits, and seam continuity. Those are checked elsewhere.
- The written file. Serialisation rounds again, to six significant digits, and nothing here reads
  it back.

Both sides integrate **exact piecewise-linear segments** in closed form, on integer ticks, rather
than raster samples. The trapezoidal rule on samples is exact for m0 and only for m0: a trapz-based
m1 disagreed with itself by 0.1 % on nothing more than a split. And integrating both sides the same
approximate way makes their errors cancel: a resample that rounded 2.5 % off a spiral's peak left a
sampled m0 agreeing to 1e-14.

**Label addresses** must equal the fold of the tree's labels. The duplicate-address check in
`check()` only fires when two addresses *collide*; an addressing shifted by one readout but still
unique passes it, which is exactly how mis-retimed labels escaped notice. Skipped when a label has no
ADC after it, since its placement is then reported as a warning rather than defined.

---

## Escape hatches

`sc.barrier(tag)` as a node forces a boundary where the compiler would not have put one — a trigger
that must sit alone, or a gradient you want split at a known instant so a later reconstruction step
can find the seam. It occupies no time, so it changes block structure and never timing.

The returned `pypulseq.Sequence` is the object itself, not a wrapper around one — modify it,
`plot()` it, `write()` it. Nothing forces you to wait for a module to be written.

---

## What a compile tells you

**It raises, or it warns. There is no report.**

```python
seq = sc.compile(tree, opts)          # opts is a pypulseq.Opts, not a seqcraft wrapper
```

That line either returns a legal `pypulseq.Sequence` or throws. There is nothing to check
afterwards, which is the point: an object carrying findings is an object whose findings can go
unread, and the failure mode is a `.seq` where the problem surfaces when someone tries to run the file rather than when it was built.

### It raises

| Exception | When |
|---|---|
| `CompileError` | Two RF or ADC events overlap; an absolute start is negative; a gradient starts off the gradient raster; no boundary can be cut in the gap between two exclusive events; a boundary would fall inside a gradient an ADC is sampling; two imaging ADCs write the same k-space address; `check_timing` fails on anything but the `TotalDuration` float-equality artifact, which pypulseq emits even on pulseq's own approved files. |
| `HardwareLimitError` | The **summed** waveform exceeds `max_grad` or `max_slew` on an axis, or one ADC or RF event carries more samples than the interpreter accepts. |
| `DefinitionConflict` | `name=` and `definitions['Name']` disagree. |
| `CompilerContractError` | An IR contract broke, or the compiled sequence does not match the tree in duration, a block's m0 or m1, or a label address. A compiler bug — or PyPulseq storing a gradient beyond its own precision, which the message names as registration; report it with the tree. |

Every message names the event, its **provenance path** and the time, and offers two concrete
remedies with the numbers filled in:

```text
HardwareLimitError: slew 189% of the 150 T/m/s limit on axis x.
  from   :  tr.readout.prephaser
  at     :  2.340 ms (block 117)
  reached:  283.5 T/m/s
  fix
    lengthen the lobe, or lower the readout bandwidth
    or design that part against sc.opts.derate(opts, slew=0.52)
```

Where several components share a block the path is their **common ancestor** rather than one
picked arbitrarily, which is the honest answer to "where did this come from".

### It warns

Things it **did** rather than refused, as `SeqCraftWarning` — a `UserWarning` subclass, through the
standard `warnings` machinery.

| Category | Means |
|---|---|
| `merge` | Two or more gradients shared an axis in one block and were summed. Names every source. The sum itself is exact. |
| `resample` | A trapezoid was summed with a raster-centre waveform. Their sum bends both on and off the gradient raster and no pulseq gradient event can hold that, so it was resampled. The message carries a **bound on how far the waveform moved**, measured rather than estimated. This is the only place the compiler is knowingly inexact. |
| `snap` | An RF or ADC reservation had to move onto the block raster. |
| `orphan_label` | A label has no ADC after it, so its placement is undefined. |
| `norm` | The vector norm across axes exceeds the per-axis limit. Normal for multi-axis gradients, and legal on real amplifiers — which is why it is not an error. |

**One aggregated warning per category, at the end of the compile**, naming the count and the
sites. That is not cosmetic: Python's default filter shows a warning once per unique
`(message, category, module, lineno)`, so one `warn` per merge would print the first and silently
swallow the other eleven — strictly worse than a count.

```python
import warnings

with warnings.catch_warnings():
    warnings.simplefilter('error', sc.SeqCraftWarning)      # treat any as fatal
    seq = sc.compile(tree, opts)

# or, in a test
with pytest.warns(sc.SeqCraftWarning, match='merge'):
    sc.compile(tree, opts)
```

---

## What the compiler reads from the scanner

Eight `Opts` fields, and nothing else:

```text
opts.max_grad          opts.rf_ringdown_time    opts.gamma
opts.max_slew          opts.rf_raster_time      opts.adc_dead_time
opts.grad_raster_time  opts.adc_raster_time
```

plus `adc_samples_limit` / `rf_samples_limit` in the post-compile size check, and
`block_duration_raster` for every boundary and duration.

That measured coupling is why there is no scanner class. `compile_sequence(root, opts, ...)` takes
the same `pp.Opts` you pass to `pp.make_trapezoid`, so a block and the compile that schedules it
cannot be built against different limits by accident. A part designed against derated limits carries
its own `Opts` from `sc.opts.derate`; the compiler validates the *combined* waveform against the one
it was given, which is what it always did.
