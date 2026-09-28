# PNS-aware physical design — Phase A, public synthetic hardware only

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-28, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../../README.md) for which one to edit.

**Status:** evidence spike complete, then **corrected and re-run** — see §0. No public API is
recommended yet. The code added is a re-runnable harness under `tools/module_mining/spikes/`,
explicitly marked experimental.
**Date:** 2026-09-28, revised the same day
**Measured against:** `main` at `5bf1953` (PR #42 merged)
**Hardware model:** `sc.hardware.synthetic_hardware()` only — pypulseq's own illustrative
`safe_example_hw()` coefficients, public provenance, **Level 1**. Nothing here is a
scanner-specific claim and nothing here clears anything for human scanning.

---

## 0. Model provenance correction and re-run

**The first run used SeqCraft's earlier hand-maintained `synthetic_hardware()` implementation.**
That implementation did not match the pypulseq `safe_example_hw()` example model named in its
documentation: `safe_example_hw()` gives the three axes distinct response parameters, and the
SeqCraft reconstruction applied **`x`'s** `tau` and `a` to `y` and `z` and used different
`g_scale` values on all three — `y` most of all. Nothing compared the two objects, so nothing
reported it. Git history shows the two differed from the first commit that introduced the
function, so this was a mismatch from the start rather than a later divergence.

`synthetic_hardware()` now derives its response parameters from `safe_example_hw()`, eliminating
SeqCraft's second coefficient table, and a regression compares the returned objects directly.
Every table below was therefore **repeated using the current SeqCraft model derived from that
public example model**; first-run numbers are quoted only where the difference is the point.

Note what this is *not*. The two response models involved — pypulseq's example model and any site
descriptor — are different models and are expected to give different estimates; neither is
established as wrong by disagreeing with the other. What was wrong was SeqCraft describing one
model and supplying another.

What that invalidated, and what it did not:

```text
SUPERSEDED   the phase-encode winder dominating the LOCAL peak, and the size of that effect
             the local/global ratio of 1.44
             the particular slew rate at which local-pass / global-fail appeared
             the ranking of the representative waveform families

STANDS       PNS must evaluate a realised waveform, not a schedule
             widening an owned region and adding TE fill are different interventions
             an owner may or may not have a lever that moves the failing peak (RETRY vs STOP)
             a leaf sees only part of its own repetition
             feedback returns through design, never by stretching emitted events
```

The waveform-specific conclusions were the model-dependent ones and they did not survive. The
architectural conclusions did not depend on the coefficients and did.

## 1. The question

`sc.pns(tree, opts, hardware)` exists and runs after the whole sequence is built. What does not
exist is **PNS-constrained minimum timing**. `design.scope.design_scope` already carries a private
`admissible` seam, documented as *"the seam a local PNS criterion will attach to, and nothing uses
it yet"*, which is handed a fully realised candidate.

So: is that seam the right place, and is it sufficient?

## 2. Every family sits near the illustrative limit at ordinary limits

At 40 mT/m / 150 T/m/s, against the current `safe_example_hw()`-derived model, every family sits
close to the example model's stimulation limit and several cross it. That is a statement about this illustrative model
and this protocol, not about any scanner.

| family | local | global | g/l | peak at | dominant axis |
|---|---|---|---|---|---|
| GRE 2D | 0.978 | 0.984 | 1.01 | 6.00 ms | x |
| MEGRE monopolar plain | 0.977 | 0.982 | 1.01 | 13.50 ms | x |
| MEGRE monopolar all-echo FC | 0.976 | **1.083** | **1.11** | 14.78 ms | x |
| MEGRE bipolar plain | 1.010 | 1.015 | 1.01 | 12.26 ms | x |
| MEGRE bipolar all-echo FC | 1.008 | **1.084** | 1.08 | 16.96 ms | x |
| EPI 2D single shot | 0.802 | 0.802 | 1.00 | 1.92 ms | y |
| Diffusion SE-EPI b=1000 | 0.993 | 0.993 | 1.00 | 51.52 ms | x |

**The dominant axis is `x` almost everywhere — the spoiler, at the end of TR** — not the
phase-encode winder the first run reported.

**Cost:** one repetition takes 7–50 ms to evaluate; a 64-line scan takes 0.22–0.42 s. A local
check is affordable inside a design search. A global check is not.

## 3. The peak is not where the sequence work was

Decomposing by which axis was compensated, on a four-echo monopolar train:

```text
flow_comp axes     local  global   peak at      x      y      z   ESP/us
None               0.977   0.982    13.50ms  0.718  0.059  0.660   2500.0
x                  0.976   0.977    14.22ms  0.718  0.027  0.660   2580.0
y                  0.976   1.316    14.15ms  0.718  0.024  0.660   2500.0
z                  0.976   0.976    14.53ms  0.718  0.017  0.660   2500.0
('x', 'y')         0.976   1.333    14.40ms  0.718  0.023  0.660   2580.0
('x', 'y', 'z')    0.976   1.083    14.78ms  0.718  0.017  0.660   2580.0
```

Three things follow:

- **Compensating `x` lowers the global peak**, 0.982 → 0.977. The inter-echo transitions PR #42
  built — the whole cost of the all-echo contract — are not PNS-binding on this model. Spreading
  the fly-back's area over a longer, gentler waveform is easier on the nerve than the fly-back was.
  This is the one waveform-specific result that also held in the first run.
- **Compensating `y` costs nothing locally and 34 % globally**, 0.982 → 1.316, with the *local*
  peak unmoved at 0.976. So its cost is **accumulation across repetitions**, not anything visible
  in one TR — which is a different mechanism from the one the first run reported, and a more
  interesting one.
- **Adding `z` to `('x','y')` reduces the global peak**, 1.333 → 1.083, because the longer echo
  time it forces lets every axis be realised more gently.

The baseline peak is at the **end of TR on x** — the spoiler, not the imaging gradients.

This remains the concrete justification for the brief's first principle. No function of nominal
timing parameters predicts that compensating one axis lowers the peak while compensating another
raises it only in the assembled scan. The evaluator has to see the realised waveform.

## 4. PNS does change AUTO timing, materially

Holding `max_grad` at 40 mT/m and lowering the slew rate until each becomes admissible:

```text
 max_slew  plain local  plain glob   FC local   FC glob  FC ESP/us  FC TR/ms
      150       0.977*      0.982*     0.976*    1.083      2580.0     16.06
      100       0.769*      0.771*     0.769*    0.829*     2670.0     16.78
       70       0.585*      0.586*     0.585*    0.630*     2780.0     17.70
       50       0.437*      0.438*     0.437*    0.472*     2940.0     18.94
       35       0.319*      0.320*     0.319*    0.346*     3480.0     21.55
       25       0.239*      0.240*     0.239*    0.262*     3620.0     23.18
```

Everything here is admissible from 100 T/m/s downward; the only failure is the assembled
compensated scan at the default 150. Under this model PNS binds in a narrow band near the default
rather than across a wide range — which is a weaker statement than the first run supported, and
the one the corrected evidence actually carries.

## 5. Local and global are not the same question

Ratios run from 1.00 (single-shot EPI and diffusion, where local *is* global) to **1.11** for the
all-echo compensated trains. The regime where that matters is present at the default slew rate:

> At `max_slew = 150 T/m/s`, the all-echo compensated repetition passes locally (0.976) and the
> assembled 64-line scan fails (1.083). Compensating `y` alone is a larger case: local 0.976,
> global 1.316.

So a locally PNS-admissible design family is **not** by itself a globally admissible scan. The
narrow statement the evidence supports is: *local PNS feedback is useful during design, but it is
not a proof about the assembled sequence unless composability is separately established.* How
often the two disagree in practice is not something one illustrative model and one protocol can
say.

## 6. The seam: the lever converges, the policy does not

The `admissible` predicate is a closure supplied by the caller, so it can capture `scope.before`
and `scope.after` and rebuild `before + designed + after` exactly as `PhysicalDesign.build` does.
**Visibility is not the problem.** The problem is what rejecting a candidate actually does.

Widening the owned window on a spiral scope, measuring the designed region alone and the whole
candidate:

```text
   flow compensation
   window/us  region alone  whole candidate  peak in arm
        1450         0.718            0.736         True
        2000         0.554            0.735         True
        3000         0.436            0.734         True
        5000         0.283            0.734         True
        9000         0.169            0.734         True
       15000         0.101            0.734         True
```

- **The region's own contribution falls monotonically**, 0.718 → 0.101, with flow compensation and
  without it (0.225 → 0.026). Widening is a genuinely convergent lever *on what the designer owns*.
- **The whole candidate does not move at all**, 0.736 → 0.734 over a tenfold window increase,
  because the peak is inside the spiral arm — which is `after`, and not the designer's to change.

This is the **STOP** case, and it is the clearest one in the study: the owner's contribution
collapses by a factor of seven and the number being tested does not move. A **RETRY** case is the
complement — an owner whose region *is* the dominant contribution, where the same widening would
move the whole candidate. Both shapes remain visible after the correction.

So the AUTO search walks to its 4000-window ceiling and refuses with *"the requirement may still be
realisable beyond 40.0 ms; if that is plausible here, raise the search ceiling"*. That advice is
wrong: no window would ever have worked.

**Widening the owned region and spending the same time as TE fill are different interventions**,
and the current model makes the point more sharply than the first run did. An explicit `te_s`
keeps the solved window and adds fill in front of it, and what that does to the global peak
depends on what is being compensated:

| requested TE / ms | `flow_comp='y'` global | `flow_comp=('x','y','z')` global |
|---|---|---|
| shortest | 1.316 | 1.083 |
| 5.504 | 1.327 | 1.038 |
| 9.004 | 1.362 | 0.994 |
| 16.004 | **1.380** | **1.002** |

Fill makes the `y`-only case steadily **worse** and the three-axis case **better**, and nothing
about the timing predicts which. Widening the owned region, by contrast, reduced the owner's
contribution monotonically in every case tested.

So the durable statement is not *"delaying is harmful"* — that came from the first run, against
the earlier SeqCraft implementation. It is: **extra elapsed time alone is not a reliable response
to a PNS failure, and widening the region the owner actually controls is.**

## 7. A leaf cannot use the same evaluator on the same terms

`CartesianLine`'s transition search has the same shape as `design_scope`'s — walk the raster, take
the first feasible — so a predicate would drop in structurally. But a leaf sees only its own block:

```text
case                          readout alone  repetition  leaf / rep
MEGRE mono plain                      0.688       0.977        0.70
MEGRE mono all-echo FC                0.559       0.976        0.57
MEGRE bipolar all-echo FC             0.846       1.008        0.84
```

A leaf-local check sees between 57 % and 84 % of what its repetition plays, so it would pass a
readout whose repetition fails. The excitation, the winder and the spoiler are all invisible to
it — and the spoiler is where the peak actually is. **The same *conceptual* evaluator can serve both, but only if the leaf is
handed context it does not currently have**, which is a larger change than this spike justifies.

## 8. Answers to the ten questions

1. **Smallest object a local evaluator needs:** a realised candidate *in the context it will
   play* — `before + designed + after`, compiled. The closure already supplies that; no new object
   is required.
2. **Is `admissible` sufficient?** Mechanism yes, policy no. It can reject, and rejection drives a
   convergent lever on the owned region, but it cannot say *"my region is already negligible and
   the peak is not mine"*, so an unreachable target exhausts the search and refuses misleadingly.
3. **Can one evaluator serve `PhysicalDesignScope` and leaf-local designers?** Conceptually yes,
   practically not yet: a leaf sees 57–84 % of its repetition. Not this spike's problem to solve.
4. **How should failure request redesign without post-hoc mutation?** It already does — rejecting
   a candidate returns to the search, which realises again at a longer window. Nothing is
   stretched. The gap is a *stopping* rule, not a redesign mechanism.
5. **Deterministic and raster-stable?** Yes. The predicate filters candidates the existing raster
   walk produces; it adds no continuous search and no new quantisation.
6. **Which becomes PNS-limited first?** On this model, the **assembled scan** of an all-echo
   compensated train — local 0.976, global 1.083 — and `flow_comp='y'` alone is worse globally
   (1.316) while costing nothing locally. The baseline peak is the spoiler on `x`. Single-shot EPI
   is the *lowest* of the seven families here (0.802), not the highest.
7. **How often does local-pass fail global?** The ratio runs 1.00–1.11 across the families, and
   the disagreement is real at the default slew rate. How often it matters in practice is not
   something one illustrative model and one protocol can establish, and this record does not
   claim it is common.
8. **What provenance is needed?** The peak's **time** and **per-axis components**, which `sc.pns`
   already returns, plus the owning region — which nothing currently records. The three peaks here
   are the spoiler, the compensated winder and the spiral arm; all three were identified by hand
   from `t` and `components`. That identification is what a global failure would need automated.
9. **Is new public API required yet?** **No.**
10. **What remains before this can drive production timing?** See §10 for what is undecided
    and §11 for validating against real hardware.

## 9. Recommendation: small internal generalisation, no public API

**The first thing this spike produced was not an architecture change but a correction**, §0: the
public example model had to be made to match the public source it named before any of the rest
was worth reading. That is done, and a regression holds it.

After that: keep `admissible` where it is and what it is. The one change the evidence justifies is that the
predicate should be able to return **why**, not just a bool — enough for `design_scope` to stop
when the owned region's contribution is already negligible and to refuse naming the region that
actually drove the peak, instead of advising a longer search.

That is an internal signature change to a private seam. It needs no public surface, no new
abstraction, and no compiler involvement.

**Explicitly not recommended yet:** any public `pns_model=` spelling, a global optimiser, a
redesign loop driven by whole-sequence failure, or leaf-local PNS. The evidence does not yet say
what those should look like.

## 9a. Where the numbers in this record come from

Every table is `tools/module_mining/spikes/pns_design_spike.py`, run against
`synthetic_hardware()` **after** it was corrected to mirror `safe_example_hw()`. The claims are
scoped to that illustrative public model, this protocol, and the realisation family this readout
implements. They are not scanner truth, and the model is explicitly not an upper bound.

## 10. What remains undecided

```text
the stopping rule            what "my region is already negligible" should mean numerically
the global loop              whether local-admissible + a global recheck is a usable workflow, or
                             whether a local/global gap of this size leaves local admissibility
                             worth having at all
leaf context                 what a leaf would have to be handed, and whether that is worth it
the public spelling          nothing decided; no model may ever be implicit
whether PNS belongs in AUTO  at all, or only as an opt-in that a caller asks for explicitly
```

## 11. Validating against real hardware

Everything above is one **example** response model, and §0 is the standing warning: a locally
maintained reconstruction whose provenance is never checked can produce confident,
model-dependent conclusions. If SeqCraft claims to expose an external example model, it should
compare the returned fields directly rather than maintain a second parameter table. Before any of
this drives production timing it should also be checked against a site descriptor — and a
disagreement there is two models disagreeing, not proof that either is wrong.

That work does not belong in this repository, and its results do not either. The rules it runs
under:

```text
the .asc stays outside the repository, read-only, and its path is passed explicitly
no raw coefficients, no mirror of the file, no perturbed "generic scanner" default
no field name relied on by name unless it is identifiable from public pypulseq or SAFE source
scanner-specific findings stay wherever the study was done, not here
```

What such a study can usefully return to a public record is **generic**: whether a conclusion is
model-dependent or structural. The architectural findings in §6 and §7 are structural — they are
about which gradients an owner controls — and should survive any response model. The
waveform-specific ones in §2 and §3 are not, and §0 is what happens when that distinction is not
made.
