# SaturationPrep — Boundary Analysis

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Stage:** prospective supervised run. Outcome open throughout.

---

## 1. What the references agree on

Across two independent designs (`reference_inventory.md` §1):

```text
a spectrally selective RF pulse, played with NO gradient
its centre frequency offset from water by the fat chemical shift
its bandwidth narrow enough not to touch water
NO rephasing afterwards -- the transverse magnetisation is unwanted
a strong spoiler immediately after, at a spatial scale far finer than a voxel
rf.use = 'saturation'
```

That last line is in both designs and is not decorative: `use` is what tells a reader, a scanner
and any analysis tool what the pulse is *for*.

## 2. What they disagree on

| | R1 / R2 (one witness) | R3 (independent) |
|---|---|---|
| flip angle | **110°** | **90°** |
| duration | **8 ms** | **12 ms** |
| bandwidth | `abs(sat_freq)` ≈ 424 Hz, *derived from the offset* | **200 Hz fixed**, via `timeBwProduct = BW × duration` |
| apodization | factory default | **0.42** |
| chemical shift | `-3.45 ppm`, sign carried in the constant | `3.5 ppm`, sign applied at use: `freqOffset = -fatOffres` |
| offset phase | R2 compensates, **R1 does not** | not compensated |

Two of these matter architecturally and the rest are protocol:

**The sign path.** One reference carries the sign in the ppm constant, the other applies it at the
point of use. Both arrive at a negative offset, by different routes, and a third implementation
combining the two conventions would get a **positive** offset and saturate the wrong side of the
spectrum — a failure that compiles, passes every timing and k-space check, and produces an image
with fat *enhanced*. This is a value a caller can get wrong.

**The offset-induced phase.** R2 compensates `φ = -2π·Δf·t_center`, R1 does not. For a saturation
pulse whose transverse magnetisation is immediately spoiled, the compensation is arguably
irrelevant — but **the references disagree and neither says why**, which is `REFERENCES_DISAGREE`
territory and must be recorded rather than resolved by preference.

Flip angle, duration and bandwidth are **protocol parameters**: two independent designs choosing
differently, both working, is what a caller-owned parameter looks like. An abstraction that picked
one would be encoding one lab's protocol.

## 3. What SeqCraft can already express, tested

Not assumed — built and inspected on this branch.

```python
Excitation(opts=opts, flip_deg=110.0, thickness_mm=None, pulse='gauss', duration_s=8e-3,
           pulse_opts={'bandwidth': 424.0, 'freq_ppm': -3.45})
```

This **works today** and produces an 8 ms gauss with `freq_ppm = -3.45`, no gradient and no
rephaser. `Excitation`'s own docstring already names the case: *"a sinc or a gauss played with no
gradient is how water excitation works."* `spoiler(opts, cycles_per_voxel=…, voxel_mm=0.1)` gives
the spoiler.

So the pulse and the spoiler are both covered. Which raises the wrapper question directly: **what
would a module add?**

## 4. The gap, which is smaller and sharper than expected

```text
>>> rf.use
'excitation'
```

`Excitation` hard-codes `use='excitation'` and lists `use` and `freq_offset` in its `_RESERVED`
set, so a caller **cannot** label the pulse. A fat saturation built the way §3 shows emits a
`.seq` in which the saturation pulse is declared an **excitation**.

That is precisely what rule B exists to catch: *if the emitted file were inspected without reading
the Python, would it tell the same physical story as the API?* Here it would not. And pypulseq
supports the right value — `get_supported_rf_uses()` returns
`('excitation', 'refocusing', 'inversion', 'saturation', 'preparation', 'other', 'undefined')`.

**The taxonomy already reserves the slot.** `modules/__init__.py` defines `preparation/` as
`rf.use in {inversion, saturation, preparation}`, and only `ir_prep.py` exists. The folder was
designed with this candidate's seat empty.

## 5. Ownership

> If a component cannot determine the correct value using only its own physical contract and
> parameters, it should not own that event.

| value | who can determine it |
|---|---|
| `rf.use = 'saturation'` | the preparation module — it is *what the module is* |
| chemical-shift offset in ppm, and its sign | the module, from a named species; the caller cannot be expected to re-derive the sign convention each time |
| bandwidth not overlapping water | needs the offset **and** the pulse duration — both the module's |
| spoiler area | `spoiler()` already, given a spatial scale |
| flip angle, duration, apodization | **the caller** — the references disagree, which is what makes them protocol |
| when the preparation is played, and how often | the caller |

Nothing here needs to reach across a leaf boundary, so this is a **leaf**, in `preparation/`,
sibling of `IRPrep` — not a kernel. The same shape `IRPrep` already has: a prep pulse, its `use`,
and a crusher.

## 6. The wrapper objection, taken seriously

`RED - WRAPPER_ONLY` asks: does extraction shorten the notebook, and is there a value the wrapper
determines that its caller could not?

**Against the candidate:** a saturation is three lines today, and two of the three already work.
The module's owned content is thin — a label, a ppm constant with its sign, and a bandwidth check.

**For it:** `spoiler()` is the precedent that thin is not disqualifying — it is a *function*, not a
Module, and it exists because one piece of arithmetic should not be got right twice. `IRPrep` is
the precedent for the Module form. And the label is not cosmetic: it is the difference between a
file that describes the experiment and one that misdescribes it.

**This is the live question for review**, and it has a third answer neither `NEW_LEAF` nor `RED`
covers: the gap in §4 could be closed by letting `Excitation` state what its pulse is for, with no
new Module at all. That would be `EXTEND_EXISTING`, and it would also fix the mislabelling for any
other spectrally selective use. It is listed in `acceptance_options.md` §2 as a real alternative,
not a strawman.

## 7. Out of scope, stated so the boundary is explicit

- **CEST saturation trains** (`pulseq-cest`): long trains, duty cycle and B1rms constraints — a
  different physical problem that should not be folded in to make one module look more general.
- **Regional / spatial saturation slabs (REST):** these *do* play a gradient and select a region;
  the shared part with spectral saturation is "prepare then spoil", which is the waveform-silhouette
  trap. Not assumed to be the same candidate.
- **Water excitation / binomial pulses:** related spectral physics, different intent.
