# SaturationPrep — Acceptance-Claim Options and Proposed Classification

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**For human review. Nothing has been implemented.**

---

## 1. Acceptance-claim options

### Option A — event identity against a SeqCraft notebook

Not available. No SeqCraft notebook builds a fat saturation; there is nothing to promote.

### Option B — waveform and spectral agreement against R1 at matched protocol *(recommended)*

R1 is executable Python, MIT, and produces the pulse directly. The claim would be:

```text
at R1's protocol -- 2.89 T, -3.45 ppm, 110 deg, 8 ms, bandwidth = |sat_freq| --
the candidate emits an RF waveform whose envelope, duration, centre frequency and
resulting spectral profile match R1's, and a spoiler of the same moment.
```

Checkable without new tooling: RF envelope comparison on the RF raster, the offset in Hz after ppm
resolution, and gradient moment integration — all existing comparator primitives.

**Claim scope, both halves.** It would establish that the candidate reproduces R1's *pulse*. It
would **not** establish that 110° / 8 ms is correct in any absolute sense (R3 disagrees), nor that
fat is actually suppressed in an image, nor anything about R2's phase compensation.

### Option C — analytic spectral invariants, independent of any implementation

Complementary to B and cheap:

```text
the pulse centre frequency equals ppm * 1e-6 * B0 * gamma, with the sign convention stated
the excited band does not overlap water:  |offset| - bandwidth/2 > 0
no gradient is played during the pulse
no rephasing gradient follows it
the spoiler moment dephases at least N cycles across the stated spatial scale
rf.use == 'saturation' in the emitted file
```

The last line is a **rule-B check promoted to an invariant**, and it is the one the current gap
would fail.

### Option D — simulation

Available (`MRzeroCore`) and would show fat suppression directly. Not proposed as an *acceptance*
path for a first pass: it would test a contrast outcome rather than the contract, and Layers 1–2
carry this candidate adequately. Recorded as available if review wants it.

### The metamorphic check worth defining before any code exists

```text
changing B0 changes the Hz offset proportionally and changes nothing else
changing the species' ppm moves the offset and nothing else
```

Definable against **R1 today**, before the candidate exists — the property that made Radial's
rotation test credible.

---

## 2. The three outcomes on the table

This candidate's interest is that a reasonable reviewer could pick any of these.

| | what it says | the case for | the case against |
|---|---|---|---|
| **`NEW_LEAF`** — `preparation/saturation_prep.py` | a sibling of `IRPrep` owning the offset, the label and the spoil | the taxonomy reserves the slot (`rf.use in {inversion, saturation, preparation}`); `IRPrep` is the exact structural precedent; the sign convention is a real trap | thin: most of the content already works via `Excitation` |
| **`EXTEND_EXISTING`** — let `Excitation` state what its pulse is for | closes the §4 gap with no new class, and fixes mislabelling for every spectrally selective use | smaller change, no new surface, and `use` is arguably `Excitation`'s to expose anyway | leaves the ppm/sign arithmetic and the spoil pairing uncomposed; a caller still assembles three things |
| **`RED - WRAPPER_ONLY`** | the module owns a label and a constant | a saturation is three lines today | `spoiler()` and `IRPrep` both set the precedent that thin arithmetic earns a home when getting it wrong is silent |

## 3. Proposed classification

```text
status:        NEW_LEAF
traffic_light: YELLOW
reason_code:   REFERENCES_DISAGREE
```

**Why `NEW_LEAF` and not the alternatives.** The ownership analysis puts every value cleanly in one
layer, the taxonomy has a seat reserved, and `IRPrep` is the same shape. `EXTEND_EXISTING` is a
genuine contender and is the recommendation I would expect review to weigh hardest — it is cheaper
and fixes the mislabelling more broadly — but on its own it leaves the caller assembling pulse,
offset arithmetic and spoiler each time, which is what `IRPrep` exists to avoid for inversion.

**Why `YELLOW` and not a straight proceed.** Two independent designs disagree on flip angle,
duration and bandwidth derivation, and the two implementations of *one* design disagree on offset
phase compensation. My reading is that the first three are protocol — caller parameters — and that
the phase compensation is irrelevant behind a strong spoiler. **Both readings are judgement, not
measurement**, and getting them wrong shapes the API: if bandwidth is derived rather than passed,
the module owns a value one reference treats as free.

So the YELLOW is narrow and named: **confirm that flip / duration / bandwidth are caller
parameters and that offset-phase compensation may be omitted.** With those two confirmed, this
proceeds to extraction under Option B + C.

**What would make it `RED`.** If extraction shows the module owns nothing but `use='saturation'`
once `Excitation` is allowed to set it — that is `WRAPPER_ONLY`, and the right answer is the
`EXTEND_EXISTING` row instead.

## 4. Calibration note

`SaturationPrep` was chosen to test a different capability from the previous run: whether the
workflow can extract a comparatively simple concept **without over-designing**. The relevant
pressure is visible already — the temptation to widen the candidate until it covers CEST trains
and regional saturation slabs, which `findings.md` §7 refuses explicitly.

No terminal outcome was encoded before the evidence was collected, and the proposal above was
reached from it. The `EXTEND_EXISTING` alternative emerged from testing what the package can
already do, not from the coarse scan.
