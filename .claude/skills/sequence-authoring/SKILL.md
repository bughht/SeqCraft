---
name: sequence-authoring
description: Write a complete, defensible MRI acquisition with SeqCraft and PyPulseq from a natural-language request — choose a representative protocol, build it, validate the emitted .seq against the physics it claims, and report what was and was not established. Use when asked to write, extend, or revise an MRI sequence. Not for editing the SeqCraft package itself.
---

# Sequence authoring

A compilable `.seq` is not a correct acquisition. PyPulseq will happily emit a legal file for a
sequence that never reaches steady state, encodes the wrong k-space, or measures nothing. The
governing rule:

> **Build the acquisition the user actually needs, then prove with the emitted file that it is
> that acquisition.**

A correct repetition is not a sequence. Most of what makes an acquisition right — preparation,
phase cycling, view ordering, segmentation, recovery, shot structure — lives *above* the
repetition, and none of it is checked by the compiler.

## What this covers

Writing a complete acquisition for a user: interpreting the request, choosing a representative
protocol, implementing it, validating the artifact, and delivering it with honest limits.

**Not covered.** Do not edit, patch, or monkey-patch the installed SeqCraft package. Do not run or
invoke the `module-mining` workflow — it is a maintainer workflow for deciding what belongs in the
package, and nothing here depends on it. Do not propose a compiler change to make a sequence work.

If a capability genuinely seems missing, say so in one sentence, give the local fallback you used,
and continue. A missing public abstraction is an observation, not a blocker.

---

## 1. Read the request

**A minimal request is a valid request.** "Write an EPI sequence" is enough to start. Do not demand
a protocol form, scanner specification, or validation plan first.

Resolve routine unspecified parameters yourself, with defensible choices you state and the user can
revise. Keep four categories visibly distinct when you report:

```text
requested       the user said it
literature      a named source specifies it for this acquisition
adapted         derived from a source but changed, with the reason
default         your choice, because the request did not constrain it
```

Never invent provenance. If you cannot reach a source, say which parts are grounded and which are
your own judgement. Write a citation from the source in hand, not from memory — author lists,
titles and page ranges are exactly what recall gets confidently wrong.

**Ask a question only when the answer changes the physics.** A focused question is warranted when
information is missing or contradictory *and* it materially changes the acquisition, cannot be
settled by a documented default, or is needed for safe deployment on real hardware. Offer the
meaningful alternatives when you ask. Questions may arise mid-build; you do not have to front-load
them.

Everything else — field strength conventions, a sensible matrix, a plausible TE — you decide and
disclose.

## 2. Choose the acquisition

Before substantial implementation, consult textbooks, reviews, original methods papers and mature
reference implementations. Then commit to a **specific representative acquisition** and say what it
is:

```text
excitation mode            selective or not, and the RF family that implies
contrast / preparation     what creates the contrast, and what it costs
encoding and readout       trajectory, dimensionality, sampling
TE / TR regime             and what constrains the minimum
RF and ADC phase           the cycling convention, and whether it survives interruptions
view ordering              and where the centre of k-space falls in it
segmentation / recovery    including the catalyzation policy, if interrupted
```

A published special case is evidence that a choice exists — not that it is the default. Do not pick
an unrepresentative protocol because it is cheaper to simulate; if cost forces a reduced version,
that is a *validation surrogate*, and the representative acquisition is still what you author.

**Design check.** Before writing much code, you should be able to state: what complete acquisition
this is, its key physical constraints, why the acquisition-level choices suit the goal, which
choices were sourced versus defaulted, and which unresolved inputs actually matter. If you cannot,
research more or ask. This is your own check, not a user approval gate.

## 3. Build it

**Separate the two kinds of correctness.** Kernel physics is within and across the repetition —
gradients, RF, ADC, timing, phase. Acquisition policy is everything above it — preparation,
cycling, ordering, segmentation, recovery. Mixing them is how a correct kernel ends up in a wrong
sequence.

State the physical requirement before you build the waveform, and anchor it to a meaningful
instant: the RF effective centre, the echo, a true ADC sample. Not a waveform midpoint that happens
to be convenient.

**Implementation preference**, following physical ownership rather than tidiness:

```text
1. a shipped SeqCraft Module
2. a sequence-local Module or function, in the user's own project
3. direct LogicBlock composition or raw PyPulseq events
```

A good acquisition often creates no new abstraction at all. Do not manufacture a class to make
generated code look modular.

**Use the shipped facilities within their actual scope.** Read the contract before relying on one;
none of these is a general solver, and none is a validation oracle:

| Facility | What it is |
|---|---|
| `sc.PhysicalDesignScope` / `PhysicalDesign` | Opt-in, **region-scoped** gradient design. You declare a window, what plays either side, which axes it may use, and which states the family comes in; it puts `k` at the echo where you asked and meets the declared intents. It is not a layer and never rewrites a finished block. |
| `sc.FlowCompensation`, `sc.VelocityEncoding` | The moment-constrained intents a scope or a kernel accepts. |
| `sc.timing` — `Raster`, `to_ticks`, `exact_sum` | Exact raster arithmetic without drift. Not an optimizer. |
| Module `min_te_s`, `min_tr_s`, `time_to_echo()` | Each kernel computes its own feasible minima and reference instants. Gap timing is arithmetic you write and can read. |
| `sc.analysis` — `moments`, `kspace`, `b_value`, `sample` | Measurement of a compiled tree. This is your Layer 1 evidence. |
| `sc.pns`, `sc.hardware.*` | A PNS **prediction** over a compiled tree, given a response model. Not a constraint the designer honours, and `synthetic_hardware()` carries illustrative public coefficients only. |

The package deliberately contains no recipes, no gap-timing solver, no `.seq` importer, and no
simulation or reconstruction — see [`docs/architecture.md`](../../../docs/architecture.md). Inherit
that; do not reimplement one locally and present it as SeqCraft capability.

Keep the public architecture intact: `LogicBlock.add(t, ...)` places events and nests blocks,
modules own physical design, you own protocol policy, and `sc.compile(tree, opts)` schedules and
checks legality.

**Use explicit `pypulseq.Opts`.** If the user gave no scanner limits, choose conservative
illustrative ones and state the values and where they came from. Never call them scanner-safe.

## 4. Validate the emitted file

Write the real `.seq`, read it back, and inspect what was actually emitted. A module's
self-reported `te_s` is not evidence that the file contains that TE.

The three layers are defined in [`docs/testing.md`](../../../docs/testing.md). Match the evidence to
the claim:

- **Layer 1** — measure the emitted events. TE and TR at the true reference instants, RF and ADC
  phase, gradient moments, k-space at real ADC samples, encoding signs and labels.
- **Layer 2** — the whole acquisition builds, writes, rereads and inspects: view count, train
  structure, the policy you said you implemented, the definitions you declared.
- **Layer 3** — what the spins do, when the claim needs it. Preparation and transient, contrast,
  off-resonance, flow, point response, reconstruction. Simulation reads the **written file**
  through `mr0.Sequence.import_file`, which is the stronger check because it tests what a scanner
  would play.

Rules that catch the common failures:

- **Legal Pulseq is necessary, not sufficient.** A balanced TR does not make a steady state. A
  segmented acquisition is not a continuous one with gaps cut in; its catalyzation policy is part
  of the acquisition's definition.
- **Simulate the file you shipped.** A reduced protocol may supplement it, but only if it preserves
  the quantity under study — not merely TR. Say what you changed and why it still answers the
  question.
- **Different evidence is not interchangeable.** Single-spin phase, voxel-ensemble magnitude,
  per-view echo amplitude, and reconstructed image behaviour are four different measurements. Equal
  echo modulation can still produce different image structure, because where the modulation lands
  in k-space matters.
- **A number that moves with the numerical budget is not a result.** Sweep the budget, say what
  converged, and exclude what did not rather than quoting it.
- **Name the quantity you measured.** A summed sidelobe magnitude is not energy and not resolution.

When evidence contradicts the design, find out which stage is wrong — protocol choice, physical
design, event implementation, emission, simulation, or reconstruction — fix it, and revalidate the
artifact rather than the belief.

## 5. Deliver

Give the user the authoring source, the generated `.seq`, the measurements you ran, the settings
and where each came from, and the parameters they can change.

**Evidence check.** Does each physical statement rest on evidence you actually produced? Report in
plain prose:

```text
validated             the claims made are supported by the evidence run
partially validated   some claims rest on evidence not yet run -- say which
blocked               a physical or technical obstacle prevents the claim
```

Say which layers ran. "Partially validated" is an honest and useful result when simulation is
unavailable; **"compiles" is never validation**, and nothing here is scanner-ready. Gradient and
slew legality under chosen `Opts`, or a PNS prediction against an illustrative model, is a design
check — not RF/SAR, acoustic, vendor-interpreter or site approval. Say so when you hand it over.

Put generated artifacts under `scratch/sequence_authoring/<case>/`, which is git-ignored. Do not
commit `.seq` files. Only a reviewed teaching acquisition belongs in `examples/`, under
[`docs/writing_examples.md`](../../../docs/writing_examples.md).

For follow-up requests, continue from the established design. Do not reopen settled literature or
architecture decisions unless the new requirement actually disturbs them.

## When to stop and ask

```text
the physics genuinely forks and no default is defensible
the request is self-contradictory
real hardware deployment is requested and site limits are unknown
a requirement cannot be realised or validated with available capabilities
```

In the last case, report the obstacle. Do not change the package, and do not quietly deliver
something weaker under the original name.
