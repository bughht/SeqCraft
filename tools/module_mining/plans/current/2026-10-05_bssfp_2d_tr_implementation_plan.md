# `bSSFP2DTR` — implementation plan

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-10-05, so that the pull request carries the direction behind the module it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../../README.md) for which one to edit.

**Scope:** one pull request, adding one public class. The 3D sibling is a separate pull
request on a separate branch and is out of scope here, including any abstraction built in
anticipation of it.

**Record this implements:** [`../bssfp/findings.md`](../bssfp/findings.md) and
[`../bssfp/candidate.yaml`](../bssfp/candidate.yaml), as corrected on 2026-10-05. Hypothesis A
is the kernel, layering D is selected, and the public API is one class per dimension.

---

## 1. The contract

`bSSFP2DTR` builds **one balanced repetition** of a 2D Cartesian bSSFP acquisition, measured
RF centre to RF centre.

### 1.1 What makes it balanced

```text
definition          the net gradient area on each axis is zero over the RF-to-RF interval
realisation         TE = TR/2 is the canonical symmetric case, and the default
variants            an asymmetric TE is still balanced if the definition still holds
```

The defining condition is the first line. `TE = TR/2` is a default, not a definition, and the
tests must encode it that way: there has to be a test that sets an explicit asymmetric TE and
then measures that the net area per axis is still zero. The bSSFP records already separate
these three levels — §5 of `findings.md` and trap (2) in `candidate.yaml` — so no record
needs reconciling before the tests are written. That reconciliation was anticipated; it turned
out not to be needed.

### 1.2 Two constraint sets, not one

This is the structural difference from `GRE2DTR`, and it drives everything below.

```text
echo constraints     over [RF centre -> echo]:  k_x = 0, k_y = requested, k_z = 0
balance constraints  over [RF centre -> next RF centre]:  M0 = 0 on x, y and z
```

A spoiled gradient echo has only the first set. The second set is what `bSSFP2DTR` adds, and
it reaches **before** the echo as well as after it: the half of the slice-selection lobe that
plays before the RF centre contributes to the balance condition even though it contributes
nothing to the echo condition, because the transverse magnetisation it would act on does not
exist yet. That is why a balanced repetition carries a slice winder on the leading side as
well as the trailing side, and it is the cleanest available demonstration that the two
constraint sets are genuinely different things.

`PhysicalDesignScope` expresses the first set (it designs a window and states `k_at_echo`).
It does not express the second. **Do not extend it in this pull request.** Realise the balance
condition inside the module and record what the module had to do for itself; whether that
belongs in the public scope is a question for after two implementations exist, not before one
does.

### 1.3 Balanced is not the same as in steady state

A repetition whose waveform satisfies the balance condition can be emitted, compiled, measured
and shown correct with the magnetisation nowhere near steady state. Steady state is a property
of a train of repetitions and of the tissue; balance is a property of one repetition's
waveform. Every test, docstring and notebook cell in this pull request must keep the two
apart, and the Layer 3 notebook must say which of the two it is showing.

---

## 2. Ownership

### 2.1 Owned by `bSSFP2DTR`

- slice-selective excitation, and its flip angle, duration and phase
- the z axis: selection, rephasing, **and** the leading winder that closes the balance
- the readout axis: prephasing, the readout itself, and the area that closes the balance
- the phase-encode axis: the requested `ky` at the echo, and the rewind that closes the balance
- TE placement, with `TR/2` as the default
- TR timing, `min_te_s` and `min_tr_s`
- the RF phase progression across repetitions
- the matching receiver phase, set from the same value in the same place
- shortest legal timing, and a refusal that names the achievable minimum when asked for less

### 2.2 Not owned by `bSSFP2DTR`

- steady-state establishment of any kind
- the start-up or catalyzation method, and how many start-up repetitions there are
- segmented k-space, the segment count, and what happens at a segment boundary
- `ky` ordering and acquisition ordering
- interruption and restart policy
- reconstruction
- 3D partition encoding

Each of these is a property of a whole acquisition, and `GRE2DTR` already sets the precedent:
its dummy repetitions live in the example's loop, not in the kernel.

### 2.3 Segmentation goes in the notebook

A segmented acquisition is another loop around the repetitions. It gets no class. The example
shows the continuous case first, then segmentation as ordinary authoring code:

```python
for segment in segments:
    startup(...)
    for ky in segment:
        scan.add(at, bssfp_2d_tr(line=ky, phase_deg=phase(n)))
```

If that loop later turns out to need a correctness condition of its own that a caller keeps
getting wrong, that is the revisit trigger. A second loop is not one.

---

## 3. The axis-by-axis design

Each axis has to satisfy its echo constraint and its balance constraint at once. The symmetric
realisation satisfies the balance constraint by construction — the trailing side is the time
mirror of the leading side — and that is the first thing to build, because it makes the
balance condition true for a structural reason that can be stated in one sentence.

| axis | leading side | at the echo | trailing side |
|---|---|---|---|
| `z` | winder, closing the balance for the selection lobe's pre-origin half | `k_z = 0` | rephaser, closing the balance for the post-origin half |
| `x` | prephaser, cancelling the readout's pre-echo area | `k_x = 0` | area equal to `area_after_echo_per_m`, with the opposite sign |
| `y` | blip to the requested `k_y` | `k_y` as requested | rewind of the same blip |

`CartesianLine` already reports `area_to_echo_per_m`, `area_after_echo_per_m` and
`prephaser_area_per_m`; `Excitation` already reports `rephaser_area_per_m` and
`rephaser_duration_s`. The leading z winder is the one quantity no existing leaf reports,
because no shipped module has needed it: for a symmetric pulse it is the same area as the
rephaser, but **derive it from the selection lobe rather than assuming the pulse is
symmetric**, and say in the docstring which it is.

The asymmetric-TE case is where "mirror the head" stops being enough and the trailing side has
to be solved rather than copied. Build the symmetric case first and the asymmetric case
second, so the test that distinguishes definition from realisation has something to fail
against.

---

## 4. The boundary question, settled from emitted waveforms

`candidate.yaml` defers implementation pending a decision about the cost of forcing gradients
to terminate at a module boundary — recorded there as 220 µs per axis per repetition for that
protocol. This pull request settles it, from what the compiler emits rather than from the
recorded wording.

### 4.1 The two readings to separate

```text
A  the efficient waveform shape depends on solving the bSSFP timing problem jointly
   -> joint physical realisation, and it belongs in the kernel / Module layer

B  two already-decided physical waveforms can be fused with no change to moments, to RF or
   ADC timing, to declared requirements or to semantics
   -> a candidate for a general, semantics-preserving compiler optimisation
```

Nothing in the compiler is to be taught anything called "bSSFP" under either reading.

### 4.2 What is already measured

Two abutting same-polarity trapezoids, added to one `LogicBlock` on one axis and compiled
(40 mT/m, 150 T/m/s, 1000 1/m each):

```text
emitted knot amplitudes (normalised)   0, 1, 1, 0, 1, 1, 0
```

The compiled waveform ramps to zero between them and back up. Replacing that pair with one
trapezoid of twice the area, on the same limits, gives:

```text
M0        2000      -> 2000       unchanged to 2e-13
M1        1.72      -> 1.45       changed by 16%
duration  1720 µs   -> 1450 µs    270 µs shorter
```

So on this evidence the fusion is not a substitution that leaves semantics alone: it changes
the first moment about the origin, and it changes the duration of the region it replaces. It
would qualify as reading B only on an axis with no declared first-moment requirement, and only
where the shortened duration is absorbed somewhere deliberate. That is a much narrower
optimisation than "fuse adjacent same-polarity lobes", and it is narrow for a reason that has
nothing to do with bSSFP.

The 270 µs here and the 220 µs in the record are the same effect at different areas and
limits; neither number is a property of the module.

### 4.3 What this pull request must do

Measure the same thing on the real composition rather than on a two-lobe probe: build the
repetition, compile it, and report per axis how much of the RF-to-RF interval is spent ramping
to and from zero at the places where the balance structure puts two same-polarity lobes next
to each other. Then classify, in the pull request and in the record:

- if the shape that removes the cost can only be found by solving the bSSFP timing problem as
  a whole, it is reading A, it stays in this module, and the compiler is not touched;
- if it is a local rewrite that provably preserves every declared requirement, it is reading
  B, and it becomes a **separately scoped** proposal with its own evidence — not part of this
  pull request.

Ship `bSSFP2DTR` with the honest cost either way. A correct balanced repetition that pays for
its ramps is shippable; an uncosted claim about where the optimisation belongs is not.

---

## 5. Reuse policy

**Build nothing shared in this pull request.** No `_BSSFPTRBase`, no dimension-switching
engine, no `mode=` parameter, no parameter that exists only so that the 3D case can set it
later. `bSSFP2DTR` is written as if it were the only one.

The 3D pull request begins by re-stating the 3D physical contract from the domain evidence,
not by copying this file and adding `partition=`. Its diff against this module is then the
evidence for what, if anything, is genuinely duplicated, and the smallest private machinery
that removes the duplication is extracted at that point. The z axis is where the two are
expected to differ most — selection plus rephasing plus balance here, against slab role plus
partition encoding plus balance there — so a shared base written now would be written around
the one axis least likely to be shareable.

---

## 6. Validation

Every measurement is made on the **emitted composition** — compile the repetition and
integrate the waveforms — not on a moment the module reports about itself. A module that
reports its own correctness is checking its arithmetic, not its output.

| # | claim | how it is measured |
|---|---|---|
| 1 | balance | `M0x`, `M0y`, `M0z` over the full RF-to-RF interval are zero |
| 2 | encoding | `k_y` at the ADC sample at the echo is the requested value |
| 3 | echo placement | the echo falls at the declared TE from the RF effective centre |
| 4 | default | with no `te_s`, the declared TE equals `TR/2` |
| 5 | definition, not realisation | an explicit asymmetric `te_s` still satisfies claim 1 |
| 6 | RF phase | the carrier phase across repetitions follows the declared progression |
| 7 | receiver phase | the ADC phase offset matches the RF phase, repetition by repetition |
| 8 | minimum timing | `min_te_s` and `min_tr_s` are achievable, and a shorter request is refused with the achievable value named |
| 9 | hardware | amplitudes and slews are within `opts`, and every event lands on its raster |
| 10 | balance is not steady state | claim 1 holds on a single repetition built in isolation, with no train around it |

Claim 10 is not redundant with claim 1. It is the one that makes the distinction in §1.3
checkable rather than merely stated.

### Layer 3

The simulation notebook shows a complete 2D acquisition and must distinguish, in its prose and
in what it plots, between *the waveform is balanced* (claims 1–10, true of each repetition from
the first) and *steady state has been established* (true only after enough repetitions, and a
property of the train and the tissue). Showing the approach to steady state is the honest way
to make the second visible; asserting it from a balanced waveform is the error this notebook
exists to prevent.

---

## 7. The teaching path

| file | what it carries |
|---|---|
| `src/seqcraft/modules/kernel/bssfp_2d_tr.py` | the module |
| `src/seqcraft/modules/__init__.py` | export, `__all__`, and the folder-map comment |
| `tests/modules/test_bssfp_2d_tr.py` | claims 1–10 |
| `docs/api_reference.md` | the class entry, and the rows in the shipped-module tables |
| `examples/bssfp_2d/01_build.ipynb` | continuous acquisition first, segmentation second |
| `examples/bssfp_2d/02_simulate_and_reconstruct.ipynb` | Layer 3 |
| `tests/modules/test_bssfp_notebooks.py` | the notebooks match the package |
| `examples/README.md`, `README.md` | the example is listed where the others are |
| `tools/run_notebook_smoke.py` | the notebooks run in CI |
| `CHANGELOG.md` | what shipped |
| `tools/module_mining/plans/bssfp/candidate.yaml` | the deferral closed, and the §4 classification recorded |

Notebooks follow `docs/writing_examples.md`; the module and its docstring follow
`docs/writing_a_module.md`, including the attribution rule. Each notebook is written in the
JSON encoding convention its own file already uses.

---

## 8. Ordered work list

1. Name the class and the file, and get an empty module importing and exported.
2. Build the symmetric case: z winder, selection, rephaser; x prephaser, readout, balance
   lobe; y blip and rewind. `TE = TR/2`.
3. Write claims 1, 2, 3 and 4 against the compiled output. Nothing else proceeds until these
   pass.
4. Add explicit `te_s`, solve the asymmetric trailing side, and write claims 5 and 10.
5. Add the RF and receiver phase progression, and write claims 6 and 7.
6. Add `min_te_s`, `min_tr_s` and the refusals, and write claims 8 and 9.
7. Run the §4.3 measurement on the finished composition, classify A or B, and write both into
   the pull request and the record.
8. Write `01_build.ipynb`: continuous first, segmented second, with the §2.3 loop.
9. Write `02_simulate_and_reconstruct.ipynb` with the §6 Layer 3 distinction.
10. Register the example, update `docs/api_reference.md` and `CHANGELOG.md`, and close the
    deferral in `candidate.yaml`.

Step 7 is a gate, not a postscript: if it comes out as reading B, the optimisation is still
not in this pull request, but the pull request has to say so and say why.

---

## 9. Open questions this plan does not decide

- **Whether any private machinery is shared with `bSSFP3DTR`.** Left to that pull request's
  diff, by §5.
- **Whether `PhysicalDesignScope` should be able to state a balance constraint.** The module
  realises it privately here. Two implementations and a composition that wants it are the
  evidence that would justify a public seam; one module is not.
- **Whether `GRE2DTR` and `bSSFP2DTR` share a repetition skeleton.** Already an open question
  in `candidate.yaml`, and refactoring a shipped kernel on the strength of this one is not
  warranted.
