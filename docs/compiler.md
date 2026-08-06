# The compiler

A logic block lets anything overlap anything. A pulseq block holds **at most one RF, one ADC and
one gradient per axis**, must last a whole number of block rasters, must leave the RF dead time
before a pulse and the ringdown after it, and must have gradients that join continuously across its
boundaries.

The compiler's job is to find block boundaries satisfying all of that, and to combine whatever lands
inside each one. It is the only complicated thing in seqcraft, which is deliberate: one hard thing,
in one file, tested hard.

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
  for gradients on Siemens, whatever the scanner reports elsewhere — and computed in integer ticks
  (`core.timing`). A plain subtraction of absolute times drifts:
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
`CompiledSequence.pns()` against a real hardware model rather than reading it off the warnings.

**And `synthetic_hardware()` cannot give you a verdict.** It exists so CI always has *a* model and it
is deliberately conservative. On the DTI example it reports **2.44** while the site's own Cima.X
descriptor reports **0.95** — the difference between "not runnable" and "passes with a thin margin".
Acting on the synthetic number would mean derating the diffusion lobes and lengthening TE to fix a
problem the scanner does not have. Load the vendor `.asc` with `load_hardware()` and
`$SEQCRAFT_ASC_DIR`; the file stays outside the repository, and only its name and sha256 are recorded.

Which term dominates is not guessable, it moves as the sequence changes, and **the two models rank
the terms differently** — so this cannot be reasoned about, only measured. On the DTI example against
the real descriptor:

| change | peak PNS |
|---|---|
| slew 0.65, spoil on z, single shot | 0.95 |
| slew 0.80 | 1.12 |
| spoil on x, y **and** z | 1.50 |
| density 0.25 instead of 1.0 (4× shorter readout) | 0.92 |

Sampling density moves it by 0.03 across a fourfold change in readout length, because PNS responds to
how hard the gradients slew over the recent past rather than to how long the readout lasts. Spoiling
on three axes instead of one is what pushes it over. The synthetic model ranked those two the other
way round.

Attribution is also not additive. With everything at 0.65 the peak sits inside the spiral; derate the
readout alone and the peak *moves* onto the refocusing crushers, so the gradient that limits you
changes as you tune. `origin(block_index)` at the peak time is how you find out which one it is.

One further caveat: a `.asc` also declares forbidden acoustic resonance bands, and nothing in seqcraft
checks them. A spiral is a narrowband gradient waveform and can sit squarely inside one.

---

## Four invariants that run on every compile

**Total duration** must equal the tree's. The tolerance scales with the total, because it has to:
float64 resolves about 4 ns at 20 000 s, so demanding nanosecond agreement on a long acquisition
would flag every sequence over an hour.

**Per-axis m0** must equal the sum over the flattened events. This is the one that catches a split
which lost a tail or a merge which dropped a piece — a whole class of compiler bug that no
individual test case would have to think of. Its tolerance scales with the total area *traversed*,
not the net: a readout and its prephaser very nearly cancel, and a relative tolerance on the net
would demand exactness that float summation over thousands of pieces cannot deliver.

m0 alone is **not** evidence of fidelity, and it is worth being precise about why. Area is exactly
what linear resampling preserves, so when both sides were integrated the same approximate way their
errors cancelled: a resample that rounded 2.5 % off a spiral's peak left m0 agreeing to 1e-14. Both
sides are now integrated from exact knots, which is what makes the agreement mean something.

**Per-axis m1**, referenced to the start of the sequence. m0 is exactly the quantity that survives a
*time shift* — a lobe playing a whole raster early leaves it untouched — so m0 alone cannot see a
gradient at the wrong moment. m1 can: a piece of area `A` displaced by `dt` changes it by `A·dt`.

Both sides are computed in **closed form**, not by quadrature. `g(t)·t` is quadratic between two
knots, so the trapezoidal rule is wrong there by an amount depending on knot spacing — enough that a
trapz-based m1 disagreed with itself by 0.1 % on nothing more than a split, which would have made it
useless as an invariant.

**Label addresses** must equal the fold of the tree's labels. The duplicate-address check in
`check()` only fires when two addresses *collide*; an addressing shifted by one readout but still
unique passes it, which is exactly how mis-retimed labels escaped notice. Skipped when a label has no
ADC after it, since its placement is then reported as a warning rather than defined.

---

## Escape hatches

`sc.barrier(tag)` as a node forces a boundary where the compiler would not have put one — a trigger
that must sit alone, or a gradient you want split at a known instant so a later reconstruction step
can find the seam. It occupies no time, so it changes block structure and never timing.

`RawEvents` wraps arbitrary pypulseq events as a module, and `CompiledSequence.seq` is the pypulseq
object itself. Nothing forces you to wait for a module to be written.

---

## Reading a compile report

```python
out = sc.compile(tree, system)
report = out.check()
report.raise_if_failed()

for issue in report.issues:
    print(issue.severity, issue.kind, issue.where, issue.message)
```

| `kind` | Severity | Means |
|---|---|---|
| `grad_merge` | warning | Two or more gradients shared an axis in one block and were summed. Names every source. The sum itself is exact. |
| `grad_resample` | warning | A trapezoid was summed with a raster-centre waveform. Their sum bends both on and off the gradient raster and no pulseq gradient event can hold that, so it was resampled. The message carries a **bound on how far the waveform moved**, measured rather than estimated. This is the only place the compiler is knowingly inexact. |
| `grad_limit`, `slew_limit` | **error** | The merged waveform exceeds the amplifier per axis. |
| `grad_norm_limit`, `slew_norm_limit` | warning | The vector norm across axes exceeds the per-axis limit. Normal for multi-axis gradients. |
| `raster` | warning | An RF or ADC start had to move onto the block raster. |
| `duration`, `moment` | **error** | An invariant failed. This is a compiler bug; please report it with the tree. |
| `timing` | error, or info | From `Sequence.check_timing`. The `TotalDuration` float-equality artifact is downgraded to info, because pypulseq emits it even on pulseq's own approved files. |
| `label` | **error** | Two imaging ADCs write the same k-space address. |

`out.origin(block_index)` gives the tag path that produced a block. Where several modules share a
block it reports their **common ancestor** rather than picking one arbitrarily, which is the honest
answer to "where did this come from".
