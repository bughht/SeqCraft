# PNS-aware physical design — Phase A, public synthetic hardware only

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-28, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../../README.md) for which one to edit.

**Status:** evidence spike complete. **No production change is proposed and no public API is
recommended yet.** The only code added is a re-runnable harness under
`tools/module_mining/spikes/`, explicitly marked experimental.
**Date:** 2026-09-28
**Measured against:** `main` at `5bf1953` (PR #42 merged)
**Hardware model:** `sc.hardware.synthetic_hardware()` only — pypulseq's own illustrative
`safe_example_hw()` coefficients, public provenance, **Level 1**. Nothing here is a
scanner-specific claim and nothing here clears anything for human scanning.

---

## 1. The question

`sc.pns(tree, opts, hardware)` exists and runs after the whole sequence is built. What does not
exist is **PNS-constrained minimum timing**. `design.scope.design_scope` already carries a private
`admissible` seam, documented as *"the seam a local PNS criterion will attach to, and nothing uses
it yet"*, which is handed a fully realised candidate.

So: is that seam the right place, and is it sufficient?

## 2. Everything fails the synthetic model at ordinary limits

At 40 mT/m / 150 T/m/s, every representative waveform exceeds the synthetic model's stimulation
limit. That is a property of the model being conservative, not of the sequences being wrong.

| family | local | global | g/l | peak at | dominant axis |
|---|---|---|---|---|---|
| GRE 2D | 1.131 | 1.211 | 1.07 | 6.00 ms | x |
| MEGRE monopolar plain | 1.132 | 1.211 | 1.07 | 13.50 ms | x |
| MEGRE monopolar all-echo FC | 1.792 | 2.572 | **1.44** | 3.97 ms | **y** |
| MEGRE bipolar plain | 1.168 | 1.213 | 1.04 | 13.50 ms | x |
| MEGRE bipolar all-echo FC | 1.793 | 2.572 | **1.44** | 3.97 ms | **y** |
| EPI 2D single shot | 1.968 | 1.968 | 1.00 | — | — |
| Diffusion SE-EPI b=1000 | 2.085 | 2.085 | 1.00 | — | — |

**Cost:** one repetition takes 7–50 ms to evaluate; a 64-line scan takes 0.22–0.42 s. A local
check is affordable inside a design search. A global check is not.

## 3. The peak is not where the sequence work was

The two all-echo compensated trains give **identical** peaks to three decimals despite completely
different inter-echo waveforms. That is the finding, not a coincidence: the peak is at 3.97 ms,
**before the first echo**, and it is on `y`.

Decomposing by which axis was compensated, on a four-echo monopolar train:

```text
flow_comp axes     local  global   peak at      x      y      z   ESP/us
None               1.132   1.211    13.50ms  0.821  0.130  0.768   2500.0
x                  1.126   1.131    14.22ms  0.821  0.060  0.768   2580.0
y                  2.203   3.482     3.84ms  0.009  2.197  0.154   2500.0
z                  1.125   1.127    14.53ms  0.821  0.038  0.768   2500.0
('x', 'y')         2.296   3.490     3.84ms  0.648  2.197  0.154   2580.0
('x', 'y', 'z')    1.792   2.572     3.97ms  0.399  1.646  0.587   2580.0
```

Three things follow, and none of them is what the sequence work of PR #42 would have predicted:

- **Compensating `x` LOWERS the global peak**, 1.211 → 1.131. The inter-echo transitions that PR
  #42 built — the whole cost of the all-echo contract — are *not* PNS-binding. Spreading the
  fly-back's area over a longer, gentler waveform is easier on the nerve than the fly-back was.
- **Compensating `y` is the entire PNS cost**, 1.211 → 3.482. The compensated phase-encode winder
  is a two-lobe reversal that has to carry a first moment in the interval before the echo.
- **Adding `z` to `('x','y')` reduces the peak**, 3.490 → 2.572, because the longer echo time it
  forces lets every axis be realised more gently.

The uncompensated baseline's peak is at the **end of TR on x and z** — the spoiler, not the
imaging gradients.

This is the concrete justification for the brief's first principle. No function of nominal timing
parameters predicts that compensating one axis lowers PNS while compensating another triples it.
The evaluator has to see the realised waveform.

## 4. PNS does change AUTO timing, materially

Holding `max_grad` at 40 mT/m and lowering the slew rate until each becomes admissible:

```text
 max_slew  plain local  plain glob   FC local   FC glob  FC ESP/us  FC TR/ms
      150       1.132       1.211      1.792     2.572      2580.0     16.06
      100       0.890*      0.928*     1.364     1.914      2670.0     16.78
       70       0.676*      0.737*     1.017     1.431      2780.0     17.70
       50       0.504*      0.594*     0.782*    1.066      2940.0     18.94
       35       0.366*      0.469*     0.573*    0.781*     3480.0     21.55
       25       0.274*      0.370*     0.445*    0.584*     3620.0     23.18
```

A plain train is admissible from 100 T/m/s; an all-echo compensated one needs 35. Getting there
moves the compensated echo spacing from 2580 to 3480 µs — **+35 % over the amplifier-limited
minimum**. PNS is the binding constraint over a wide and realistic range, not a rare edge case.

## 5. Local does not predict global

Ratios run from 1.00 (single-shot EPI and diffusion, where local *is* global) to **1.44** for the
all-echo compensated trains. And there is a regime where the distinction bites:

> At `max_slew = 50 T/m/s`, the all-echo compensated repetition passes locally (0.782) and the
> assembled 64-line scan fails (1.066).

So a locally PNS-admissible design family is **not** a globally admissible scan, exactly as the
brief anticipated. Local admissibility is useful for finding a candidate family; it cannot be the
acceptance criterion.

## 6. The seam: the lever converges, the policy does not

The `admissible` predicate is a closure supplied by the caller, so it can capture `scope.before`
and `scope.after` and rebuild `before + designed + after` exactly as `PhysicalDesign.build` does.
**Visibility is not the problem.** The problem is what rejecting a candidate actually does.

Widening the owned window on a spiral scope, measuring the designed region alone and the whole
candidate:

```text
   flow compensation
   window/us  region alone  whole candidate  peak in arm
        1450         0.862            1.901         True
        2000         0.665            1.896         True
        3000         0.528            1.896         True
        5000         0.350            1.896         True
        9000         0.223            1.897         True
       15000         0.135            1.897         True
```

- **The region's own contribution falls monotonically**, 0.862 → 0.135, with flow compensation and
  without it. Widening is a genuinely convergent lever *on what the designer owns*.
- **The whole candidate does not move at all**, 1.901 → 1.897 over a twentyfold window increase,
  because the peak is inside the spiral arm — which is `after`, and not the designer's to change.

So the AUTO search walks to its 4000-window ceiling and refuses with *"the requirement may still be
realisable beyond 40.0 ms; if that is plausible here, raise the search ceiling"*. That advice is
wrong: no window would ever have worked.

**A distinction worth recording, because it nearly went the other way.** Requesting a longer echo
time through `te_s=` is *not* the same operation, and it moves PNS the other way:

| TE / ms | winder total \|area\| 1/m | winder peak \|G\| kHz/m | local PNS |
|---|---|---|---|
| 3.764 | 534 | 595 | 2.203 |
| 5.504 | 764 | 652 | 2.408 |
| 9.004 | 1108 | 711 | 2.519 |
| 16.004 | 1615 | 907 | 2.640 |

An explicit TE keeps the window and spends the extra time as fill, while the first moment the
winder must cancel grows with the lever arm to the echo — so the same waveform duration has to
carry more area, and PNS rises. **Widening the window and delaying the echo are opposite
interventions**, and only the first is a convergent response to a PNS rejection.

## 7. A leaf cannot use the same evaluator on the same terms

`CartesianLine`'s transition search has the same shape as `design_scope`'s — walk the raster, take
the first feasible — so a predicate would drop in structurally. But a leaf sees only its own block:

```text
case                          readout alone  repetition  leaf / rep
MEGRE mono plain                      0.786       1.132        0.69
MEGRE mono all-echo FC                0.639       1.792        0.36
MEGRE bipolar all-echo FC             0.967       1.793        0.54
```

A leaf-local check would pass a readout whose repetition fails, by a factor of up to three. The
excitation, the winder and the spoiler are all invisible to it — and the spoiler is the peak in
the uncompensated case. **The same *conceptual* evaluator can serve both, but only if the leaf is
handed context it does not currently have**, which is a larger change than this spike justifies.

## 8. Answers to the ten questions

1. **Smallest object a local evaluator needs:** a realised candidate *in the context it will
   play* — `before + designed + after`, compiled. The closure already supplies that; no new object
   is required.
2. **Is `admissible` sufficient?** Mechanism yes, policy no. It can reject, and rejection drives a
   convergent lever on the owned region, but it cannot say *"my region is already negligible and
   the peak is not mine"*, so an unreachable target exhausts the search and refuses misleadingly.
3. **Can one evaluator serve `PhysicalDesignScope` and leaf-local designers?** Conceptually yes,
   practically not yet: a leaf sees 36–69 % of its repetition. Not this spike's problem to solve.
4. **How should failure request redesign without post-hoc mutation?** It already does — rejecting
   a candidate returns to the search, which realises again at a longer window. Nothing is
   stretched. The gap is a *stopping* rule, not a redesign mechanism.
5. **Deterministic and raster-stable?** Yes. The predicate filters candidates the existing raster
   walk produces; it adds no continuous search and no new quantisation.
6. **Which becomes PNS-limited first?** Not a sequence family — an **axis**. The compensated
   phase-encode winder, on any train that asks for `y` compensation. EPI and diffusion have the
   highest absolute peaks (1.97, 2.09) but no repetition structure to accumulate.
7. **How often does local-pass fail global?** Across the tested range, whenever the ratio exceeds
   the headroom: 1.00 for single-shot, 1.07 for plain trains, **1.44** for compensated trains, with
   a demonstrated regime (50 T/m/s) where local passes and global does not.
8. **What provenance is needed?** The peak's **time** and **per-axis components**, which `sc.pns`
   already returns, plus the owning region — which nothing currently records. The three peaks here
   are the spoiler, the compensated winder and the spiral arm; all three were identified by hand
   from `t` and `components`. That identification is what a global failure would need automated.
9. **Is new public API required yet?** **No.**
10. **What should Phase B measure?** See §10.

## 9. Recommendation: small internal generalisation, no public API

Keep `admissible` where it is and what it is. The one change the evidence justifies is that the
predicate should be able to return **why**, not just a bool — enough for `design_scope` to stop
when the owned region's contribution is already negligible and to refuse naming the region that
actually drove the peak, instead of advising a longer search.

That is an internal signature change to a private seam. It needs no public surface, no new
abstraction, and no compiler involvement.

**Explicitly not recommended yet:** any public `pns_model=` spelling, a global optimiser, a
redesign loop driven by whole-sequence failure, or leaf-local PNS. The evidence does not yet say
what those should look like.

## 10. What remains undecided

```text
the stopping rule            what "my region is already negligible" should mean numerically
the global loop              whether local-admissible + a global recheck is a usable workflow, or
                             whether the 1.44 ratio makes local admissibility not worth having
leaf context                 what a leaf would have to be handed, and whether that is worth it
the public spelling          nothing decided; no model may ever be implicit
whether PNS belongs in AUTO  at all, or only as an opt-in that a caller asks for explicitly
```

## 11. Phase B, when a private `.asc` is supplied

Not started, and not to be started until this is reviewed. When it is: the file stays outside the
repository, read-only, and the study reports **which SAFE fields are consumed, field-to-model
mapping, per-axis versus global behaviour, normalised sensitivity and qualitative ranking** —
never raw coefficients, never a perturbed "generic scanner", and never a field name whose public
provenance is not established from pypulseq or SAFE source.

The specific question Phase B should answer, given §3: **does a real model agree that the
compensated phase-encode winder, rather than the readout train, is the PNS-dominant element?** If
the synthetic and real models disagree on that ranking, the synthetic model is not usable even for
algorithm development, and §9's recommendation would have to be revisited.
