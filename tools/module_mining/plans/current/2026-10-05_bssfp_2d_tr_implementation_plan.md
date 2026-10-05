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

### 1.1a What the default argument actually does

`TE = TR/2` is not only the conceptual canonical case. It is the API behaviour, and it is
**deliberately unlike the GRE-style reading of `None`**:

```text
te_s=None, tr_s=None        choose the shortest legal SYMMETRIC repetition, then TE = TR/2
te_s=None, tr_s=explicit    TE = TR/2
te_s=explicit               the caller is asking for an asymmetric-TE realisation;
                            solve around that TE and keep RF-to-RF balance
```

In `GRE2DTR`, `te_s=None` means *the earliest echo this protocol can reach*. Here it does not.
`tr_s=None` still asks for the shortest legal repetition, but the thing being shortened is the
**symmetric** one, so the minimum `bSSFP2DTR` reports is the shortest TR whose half is a
reachable TE — not the shortest TR that can place an echo anywhere. Writing the default as
"minimum TE" would silently ship an asymmetric acquisition to every caller who passed nothing,
and asymmetry is a choice with a contrast consequence.

The exception is scoped to the trajectory, not to encoding in general:

```text
canonical symmetric readout / trajectory      TE = TR/2 by default
intentional asymmetric acquisition            TE may differ from TR/2
```

An intentionally asymmetric acquisition is one whose *readout or trajectory* is asymmetric:
readout-direction partial Fourier (partial echo), an explicitly placed asymmetric echo, or a
UTE-like trajectory. `CartesianLine.partial_fourier` is readout-direction partial echo, so it
is one of these. **Phase-encode-direction partial Fourier is not**: it changes which `ky` lines
are acquired and leaves the echo where it was, so it does not by itself call for an asymmetric
TE. Nothing in the module should couple the two.

### 1.2 Two constraint sets, not one

This is the structural difference from `GRE2DTR`, and it drives everything below.

```text
echo constraints     over [RF centre -> echo]:  k_x = 0, k_y = requested, k_z = 0
balance constraints  over [RF centre -> next RF centre]:  M0 = 0 on x, y and z
```

A spoiled gradient echo has only the first set. The second set is what `bSSFP2DTR` adds. The
two intervals begin at the same instant and end at different ones: both start at RF centre `n`,
the echo interval ends at the echo, and the balance interval runs on to RF centre `n+1`. So
what the balance interval contains is:

```text
from   RF centre n        the post-centre half of THIS repetition's selection gradient
then                      every gradient played during the repetition
to     RF centre n+1      the pre-centre half of the NEXT repetition's selection gradient
```

The pre-centre half of *this* repetition's selection gradient is outside the interval: it
belongs to the interval that ended at RF centre `n`, which is repetition `n-1`'s.

That is why a balanced repetition carries a slice winder on the leading side as well as the
trailing side — and it is why the z balance is settled by **four** terms contributed by **two**
repetitions:

```text
A_post(n)    + R_tail(n)    +    W_lead(n+1)  + A_pre(n+1)   = 0
\______________________/         \______________________/
  contributed by n                 contributed by n+1
```

```text
A_post(n)     area of repetition n's selection gradient after RF centre n
R_tail(n)     the trailing z lobe repetition n emits
W_lead(n+1)   the leading z winder repetition n+1 emits
A_pre(n+1)    area of repetition n+1's selection gradient before RF centre n+1
```

For identical symmetric selection gradients `A_post = A_pre = A`, and the natural symmetric
realisation is:

```text
R_tail = -A        W_lead = -A        A - A - A + A = 0
```

**Each balancing lobe cancels one selection half, not two.** Setting
`R_tail = -(A_post + A_pre)` while the leading winder is still non-zero double-counts the
neighbour's half, and the interval then carries a net `-A`. The leading winder contributes
nothing to this repetition's echo condition either — the transverse magnetisation it would act
on does not exist until the RF centre — so it is a lobe that serves only a balance constraint
shared with a neighbour, which is the clearest demonstration available that the two constraint
sets are different things.

**Consequence for ownership: the contract is compositional.** No single emitted repetition is
self-contained. What `bSSFP2DTR` guarantees is:

> A `bSSFP2DTR` configuration defines a repetition contract such that two **compatible**
> consecutive repetitions built from that configuration satisfy zero net gradient area on each
> axis over the RF-centre-to-RF-centre interval between them.

**Compatible** has to be stated, not assumed. Adjacent repetitions must preserve the
selection-gradient geometry that meets at the boundary:

```text
required to match        slice thickness
                         RF duration / pulse design
                         selection-gradient waveform and timing

free to differ           flip angle, when it changes RF amplitude and leaves Gz unchanged
                         ky, which is a y-axis encoding; the rewind of n and the blip of n+1
                         still have to satisfy the y balance across the interval
```

A neighbour built by something else, or with a different slice thickness or pulse duration,
is outside what this contract covers. The docstring says that rather than implying a guarantee
the module cannot hold.

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
- the z axis: selection, the trailing lobe cancelling its own post-centre half, and the
  leading winder cancelling its own pre-centre half
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
| `z` | `W_lead = -A_pre`, cancelling **this** lobe's pre-centre half, inside the interval that ended at this RF centre | `k_z = 0` | `R_tail = -A_post`, cancelling **this** lobe's post-centre half, inside the interval that starts at this RF centre |
| `x` | prephaser, cancelling the readout's pre-echo area | `k_x = 0` | area equal to `area_after_echo_per_m`, with the opposite sign |
| `y` | blip to the requested `k_y` | `k_y` as requested | rewind of the same blip |

Read the z row as three pieces an emitted repetition plays, each settling **its own** selection
half. Neither balancing lobe reaches across the RF centre to cancel a neighbour's half: the
leading winder closes the interval that ended at this RF centre and the trailing lobe opens the
one that starts there, and each carries `-A`, not `-2A`. For a symmetric pulse both lobes have
the same area, which is also what a spoiled gradient echo's rephaser carries — the extra
gradient in a balanced repetition is a second lobe, not a larger one.

`CartesianLine` already reports `area_to_echo_per_m`, `area_after_echo_per_m` and
`prephaser_area_per_m`; `Excitation` already reports `rephaser_area_per_m` and
`rephaser_duration_s`. The pre-centre half of the selection lobe is the one quantity no
existing leaf reports, because no shipped module has needed it: for a symmetric pulse it is
the same area as the rephaser, but **derive it from the selection lobe rather than assuming
the pulse is symmetric**, and say in the docstring which it is.

The asymmetric-TE case is where "mirror the head" stops being enough and the trailing side has
to be solved rather than copied. Build the symmetric case first and the asymmetric case
second, so the test that distinguishes definition from realisation has something to fail
against.

---

## 4. The boundary cost, and which layer may pay it

`candidate.yaml` defers implementation pending a decision about the cost of forcing gradients
to terminate at a module boundary — recorded there as 220 µs per axis per repetition for that
protocol. This pull request settles it, from what the compiler emits rather than from the
recorded wording.

### 4.1 Where the line between the layers is

```text
compiler representation optimisation
    the emitted physical G(t) is UNCHANGED
    RF and ADC timing are UNCHANGED
    only the pulseq representation and the block partition change

physical redesign / joint realisation
    anything else -- gradient shape, duration, M1, any higher moment, any other
    physical property of the waveform
```

The test is the emitted waveform, not the declared requirements. "This rewrite preserves every
requirement the current module happened to declare" is **not** a sufficient criterion, because
it licenses the compiler to reason as follows:

```text
M1 is undeclared here, therefore changing M1 is harmless
```

That is an MRI-physics judgement about which moments matter for a pathway, and the compiler is
the layer that does not make those. A rewrite that changes `G(t)` is physical redesign whether
or not anyone wrote down a requirement it would have violated.

Nothing in the compiler is to be taught anything called "bSSFP" under any outcome below.

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

`G(t)` changes, the duration changes and `M1` changes. By §4.1 that settles it: **this is
physical redesign, not representation optimisation.** The minimal reproducer has rejected a
general gradient-fusion pass in the compiler, and this pull request is not to revisit that.

The 270 µs here and the 220 µs in the record are the same effect at different areas and
limits; neither number is a property of the module.

### 4.3 What this pull request must do

The remaining question is not *which layer* — §4.2 answered that — but what this module does
about a cost that is now known to live in its own layer. Measure it on the real composition
rather than on a two-lobe probe: build the repetition, compile it, and report per axis how much
of the RF-to-RF interval is spent ramping to and from zero where the balance structure puts two
same-polarity lobes next to each other. Then choose, in the pull request and in the record:

```text
A  jointly realise the more efficient balanced waveform inside the Module / kernel layer,
   where changing G(t) is a decision this layer is allowed to make

C  keep the composable zero-at-the-boundary realisation, and quantify and accept the TR cost
```

**Both are valid outcomes, and C is not a failure.** The composable realisation is the one that
lets every leaf keep its own correctness condition, and paying measured microseconds for that
may well be worth it. Do not delay shipping a correct kernel
because a hand-written waveform would be shorter; ship the measured cost alongside it, and let
a protocol that cannot afford it be the evidence for A.

Under neither outcome does the compiler learn anything candidate-specific.

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
| 1 | balance | build `N + 1` **compatible** consecutive repetitions, compile them, and integrate `Gx`, `Gy`, `Gz` over each of the first `N` RF-centre-to-RF-centre intervals: `M0x = M0y = M0z = 0` |
| 2 | encoding | `k_y` at the ADC sample at the echo is the requested value |
| 3 | echo placement | the echo falls at the declared TE from the RF effective centre |
| 4 | default | with no `te_s`, the declared TE equals `TR/2` |
| 5 | definition, not realisation | an explicit asymmetric `te_s` still satisfies claim 1 |
| 6 | RF phase | the carrier phase across repetitions follows the declared progression |
| 7 | receiver phase | the ADC phase offset matches the RF phase, repetition by repetition |
| 8 | minimum timing | `min_te_s` and `min_tr_s` are achievable, and a shorter request is refused with the achievable value named |
| 9 | hardware | amplitudes and slews are within `opts`, and every event lands on its raster |
| 10 | balance is not steady state | claim 1 holds on the **first** RF-to-RF interval of a two-repetition sequence, where no magnetisation is anywhere near steady state |
| 11 | balance across unlike neighbours | claim 1 holds for adjacent repetitions with **different `ky`**, and for adjacent repetitions with a **different flip angle and an unchanged selection gradient** |

Claim 10 is not redundant with claim 1. It is the one that makes the distinction in §1.3
checkable rather than merely stated: balance is a property the waveform has from the first
interval, before any train has had time to establish anything.

Claim 11 is what keeps claim 1 from being satisfied by an accident of symmetry. Measuring a
pair of identical repetitions cannot distinguish a module that balances the physical RF-to-RF
interval from one that merely zeroes each `LogicBlock` edge to edge, because for identical
neighbours the two agree. They stop agreeing when the neighbours differ, so the pair must
differ:

```text
different ky                  the y rewind of n and the blip of n+1 still have to satisfy the
                              y balance across the interval; this separates a rewind that
                              belongs to the repetition that encoded it from one that merely
                              returns a block edge to zero

different flip angle,         the z balance has to survive a neighbour with a different RF
unchanged Gz                  amplitude, because the flip angle is carried by the RF and not
                              by the selection gradient -- which is what makes a start-up ramp
                              compatible in the sense of §1.2
```

Both pairs are **compatible** by §1.2: they differ only in ways the contract allows. A pair
that differed in slice thickness or pulse duration would be outside the contract, so a failure
there would say nothing about the module.

The acceptance procedure, stated once:

```text
build N+1 compatible consecutive repetitions
measure the first N RF-centre-to-RF-centre intervals
verify M0x = M0y = M0z = 0
```

The last repetition has no successor and therefore no interval of its own. **The test measures
the physical RF-to-RF definition, not a module-container boundary that happens to resemble
it.** A test written against block edges would pass on a module that is wrong in exactly the
way §1.2 describes — and, because each balancing lobe cancels one selection half rather than
two, it would also pass on a module whose trailing lobe double-counts.

### Layer 3

The simulation notebook shows a complete 2D acquisition and must distinguish, in its prose and
in what it plots, between *the waveform is balanced* (claims 1–11, true of each interval from
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
| `tests/modules/test_bssfp_2d_tr.py` | claims 1–11 |
| `docs/api_reference.md` | the class entry, and the rows in the shipped-module tables |
| `examples/bssfp_2d/01_build.ipynb` | continuous acquisition first, segmentation second |
| `examples/bssfp_2d/02_simulate_and_reconstruct.ipynb` | Layer 3 |
| `tests/modules/test_bssfp_notebooks.py` | the notebooks match the package |
| `examples/README.md`, `README.md` | the example is listed where the others are |
| `tools/run_notebook_smoke.py` | the notebooks run in CI |
| `CHANGELOG.md` | what shipped |
| `tools/module_mining/plans/bssfp/candidate.yaml` | the deferral closed, and the §4.3 outcome recorded |

Notebooks follow `docs/writing_examples.md`; the module and its docstring follow
`docs/writing_a_module.md`, including the attribution rule. Each notebook is written in the
JSON encoding convention its own file already uses.

---

## 8. Ordered work list

1. Name the class and the file, and get an empty module importing and exported.
2. Build the symmetric case: z winder, selection, rephaser; x prephaser, readout, balance
   lobe; y blip and rewind. `TE = TR/2`.
3. Write claims 1, 2, 3 and 4 against `N + 1` **compatible consecutive** compiled
   repetitions, measuring the first `N` intervals RF centre to RF centre. Nothing else
   proceeds until these pass, and nothing in the test may use a `LogicBlock` edge as an
   interval endpoint.
4. Write claims 10 and 11 — the two-repetition sequence, the differing `ky`, the differing
   flip angle with unchanged `Gz` — before adding any further parameter. These are the claims
   that catch both §1.2 errors, the wrong interval and the double-counted half, so they are
   worth more early than late.
5. Resolve `te_s=None` / `tr_s=None` to the shortest legal **symmetric** repetition, add
   explicit `te_s`, solve the asymmetric trailing side, and write claim 5.
6. Add the RF and receiver phase progression, and write claims 6 and 7.
7. Add `min_te_s`, `min_tr_s` and the refusals, and write claims 8 and 9. The minimum TR is
   the shortest symmetric one, per §1.1a.
8. Run the §4.3 measurement on the finished composition, choose A or C, and write the number
   and the choice into the pull request and the record.
9. Write `01_build.ipynb`: continuous first, segmented second, with the §2.3 loop.
10. Write `02_simulate_and_reconstruct.ipynb` with the §6 Layer 3 distinction.
11. Register the example, update `docs/api_reference.md` and `CHANGELOG.md`, and close the
    deferral in `candidate.yaml`.

Step 8 is a gate on the *record*, not on shipping: outcome C ships the same module with the
cost written down. What may not happen is shipping with the cost unmeasured.

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
- **Whether a protocol exists that cannot afford outcome C.** That protocol, not a shorter
  hand-written waveform, is the evidence that would justify outcome A. If §4.3 comes out as C,
  this becomes the revisit trigger.
