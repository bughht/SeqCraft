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
phase cycling, view ordering, segmentation, recovery — lives *above* the repetition, and none of it
is checked by the compiler.

## What this covers

Writing a complete acquisition for a user: interpreting the request, choosing a representative
protocol, implementing it, validating the artifact, and delivering it with honest limits.

**Not covered.** Do not edit, patch, or monkey-patch the installed SeqCraft or PyPulseq packages.
Do not run or invoke the `module-mining` workflow — it is a maintainer workflow for deciding what
belongs in the package, and nothing here depends on it. Do not propose a compiler change to make a
sequence work.

If a capability genuinely seems missing, say so in one sentence, give the local fallback you used,
and continue. A missing public abstraction is an observation, not a blocker.

---

## 0. Establish the environment first

**Two environments are in scope, and they are not the same job.**

```text
a SeqCraft source checkout     the repository is present: its docs, examples and the
                               module source are readable, and the installed package may
                               be an editable install of exactly that tree

an ordinary Python project     SeqCraft and PyPulseq are installed, e.g. by pip, and there
                               is no SeqCraft repository anywhere
```

Before designing anything, find out which one you are in and what is actually installed:

```python
import importlib.metadata as md
import seqcraft as sc, pypulseq as pp
print(sc.__version__, md.version('pypulseq'))
print(sorted(n for n in dir(sc.modules) if not n.startswith('_')))
```

Then check the specific surfaces you intend to rely on — a module's constructor signature, whether
a facility exists at all — by introspecting the installed package rather than by recalling what the
development branch has. **An installed release does not contain everything on `main`.** If
something you expected is absent, say so and compose what you need locally instead.

This file is self-contained. Where it points at repository documents, those are *supplementary*
and exist only in a source checkout; nothing here requires them at runtime.

## 1. Read the request

**A minimal request is a valid request.** "Write an EPI sequence" is enough to start. Do not demand
a protocol form, scanner specification, or validation plan first.

**Say what you are about to build, in one or two sentences, then start.** When the request does not
set the scope or the depth of checking, name both and keep going — the user can redirect you while
you work:

> I'll build a compact, complete EPI example and validate the emitted sequence, with focused
> physical and image-level checks where they matter. Say so if you'd rather have a quick
> unvalidated draft, or a different scope.

That is a courtesy, not a gate. Do not wait for an answer, do not turn it into a questionnaire or
a template, and do not repeat it on every follow-up turn.

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

### How big should it be?

**Scope follows the request.** For a short or open-ended one, build the **smallest complete,
physically defensible, modifiable acquisition** that answers it — and *complete* is the operative
word. A readout waveform on its own is not an acquisition; whatever preparation, encoding and
train-level policy the requested contrast physically requires is part of the smallest version, not
an extra.

```text
"write an EPI sequence"            a compact but complete EPI, sensible stated defaults
"write a BOLD fMRI EPI protocol"   BOLD-specific, multislice if that is what it takes
"...and validate it thoroughly"    same compact demo, deeper checks
```

**Optional features need a reason, not a precedent.** Fat saturation, many slices or volumes,
navigators, reference or calibration scans, elaborate preparation — none of these belongs in the
build merely because a paper or a reference implementation has it. Each adds RF, timing,
assumptions and validation you then owe. If the user's goal or the physics calls for one, add it
and say why. If not, leave it out and mention that you did.

The failure to avoid is answering a narrow question with a clinical protocol: asked for one
continuous EPI readout, do not deliver thirty slices and fat saturation whose mechanism nothing in
your evidence then tests.

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
1. a shipped SeqCraft Module          check it exists in the installed version
2. a sequence-local Module or function, in the user's own project
3. direct LogicBlock composition or raw PyPulseq events
```

A good acquisition often creates no new abstraction at all. Do not manufacture a class to make
generated code look modular. Sequence-local code belongs to the user's project, never to the
installed package.

**Use the shipped facilities within their actual scope.** None is a general solver; confirm each is
present in the installed version before relying on it:

| Facility | What it is |
|---|---|
| `sc.LogicBlock.add(t, ...)` | Places an event or nests a block at a time you compute. The composition primitive. |
| `sc.compile(tree, opts)` | Schedules and checks legality, returning a `pypulseq.Sequence`. Refusal is a hardware or raster violation, not a physics review. |
| `sc.PhysicalDesignScope` / `PhysicalDesign` | Opt-in, **region-scoped** gradient design. You declare a window, what plays either side, which axes it may use, and which states the family comes in; it puts `k` at the echo where you asked and meets the declared intents. It is not a layer and never rewrites a finished block. |
| `sc.FlowCompensation`, `sc.VelocityEncoding` | The moment-constrained intents a scope or a kernel accepts. |
| `sc.timing` — `Raster`, `to_ticks`, `exact_sum` | Exact raster arithmetic without drift. Not an optimizer. |
| Module `min_te_s`, `min_tr_s`, `time_to_echo()` | Each kernel computes its own feasible minima and reference instants. Gap timing is arithmetic you write and can read. |

SeqCraft deliberately ships no recipes, no gap-timing solver, no `.seq` importer, and no simulation
or reconstruction. Inherit that; do not reimplement one locally and then present it as SeqCraft
capability.

**Use explicit `pypulseq.Opts`.** If the user gave no scanner limits, choose conservative
illustrative ones and state the values and where they came from. Never call them scanner-safe.

## 4. Validate — and know what each measurement looks at

**Validation is on by default, and its depth is a separate question from the acquisition's size.**
A one-slice demo can deserve careful checking; an elaborate protocol can honestly ship partly
validated. Do not let one dimension set the other.

By default: emit the `.seq`, read it back, independently check the invariants that matter, and add
targeted simulation or reconstruction where a signal-level or image-level claim needs it. That is
the floor, not a licence to run every sweep — a convergence study, a sensitivity analysis or a
hardware investigation is warranted by a specific risk or an unresolved discrepancy, not by
thoroughness for its own sake.

**If the user asks for a quick draft without validation, give them that.** Hand over what exists,
mark it **unvalidated** where they will see it, list the checks you did not run, and claim nothing
about physical accuracy. Never imply an unperformed check passed. What you must not do, in either
mode, is call anything scanner-ready.

**Unvalidated means independent validation was not performed** — the read-back, the checks on the
emitted file, any simulation or reconstruction. It does not switch off the checks that building,
compiling or writing the file run by themselves. Never bypass those to get a file out, and report
truthfully that they ran and what they found. Passing them is not independent validation of the
acquisition, and does not make it scanner-ready.

Three different things get measured, and conflating them is the most common way an authoring
session convinces itself it is finished.

```text
the tree you built        sc.moments   exact, over flatten(tree)
                          sc.sample    interpolated onto a uniform raster grid -- not exact
                          sc.b_value   numerical, inherits sample's grid

compiled internally       sc.kspace    compiles, then calculate_kspacePP(); exact at true ADC times
                          sc.pns       compiles, then calculate_pns(); a model-based prediction

the emitted file          read the written .seq back with pypulseq and measure it, or import it
                          into a simulator.  Nothing in sc.analysis does this for you.
```

**`sc.analysis` measures your tree, not your file.** That matters, because the file is what a
scanner plays. Write the `.seq`, read it back, and derive your expectations from the protocol and
from physics — `k_y = (line - n/2) / FOV_y`, TE from the RF effective centre to the echo, TR
between successive excitations. A module's self-reported `te_s` is not evidence that the file
contains that TE.

Evidence comes in three layers. In a source checkout, `docs/testing.md` defines them and each
workflow's acceptance policy lives with that workflow; the working summary:

- **Layer 1 — analytic.** Invariants measured on the emitted events: TE and TR at true reference
  instants, RF and ADC phase, gradient moments, k-space at real ADC samples, encoding signs and
  labels.
- **Layer 2 — build and readback.** The complete acquisition builds, writes a real `.seq`, and is
  read back and inspected: view count and order, train structure, the policy you said you
  implemented, the definitions you declared.
- **Layer 3 — simulation and reconstruction.** What the spins do, when the claim needs it.
  Preparation and transient, contrast, off-resonance, flow, point response, reconstruction.

Layer 3 is not owed by every primitive. It *is* owed when the claim is about train-level or
image-level behaviour, because nothing below it can test that claim.

If you simulate, import the **written file** (for example `mr0.Sequence.import_file`) rather than
the tree, since that tests what a scanner would play. Importing cleanly is evidence that the file
is well formed and that its physics simulates as expected — it is **not** evidence that a scanner
will execute it.

Rules that catch the common failures:

- **Legal Pulseq is necessary, not sufficient.** A balanced TR does not make a steady state. A
  segmented acquisition is not a continuous one with gaps cut in; its catalyzation policy is part
  of the acquisition's definition.
- **Coverage is not ordering.** Checking that every view was acquired says nothing about the order
  they were acquired in, and view order decides which lines are collected while the magnetisation
  is still settling.
- **A check that imports its expectation from your own authoring code agrees with itself.** Reading
  the protocol back out of the module or the JSON you just wrote measures emission, which is worth
  measuring — but it cannot catch a wrong protocol choice or a mistake common to both sides. Derive
  at least the central invariants from the physics instead: `k_y = (line - n/2) / FOV_y`, TE from
  the RF centre to the echo.
- **Simulate the file you shipped, and say what the surrogate dropped.** A reduced protocol may
  supplement it if it preserves the quantity under study — not merely TR. A one-slice extract with
  a preparation pulse replaced by an equal delay is fine evidence about encoding and useless
  evidence about that preparation; name which it is.
- **Different evidence is not interchangeable.** Single-spin phase, voxel-ensemble magnitude,
  per-view echo amplitude, and reconstructed image behaviour are four different measurements. Equal
  echo modulation can still produce different image structure, because where the modulation lands
  in k-space matters.
- **A number that moves with a numerical control is not a result — and the budget is not the only
  control.** Raising a state count while leaving a pruning or truncation threshold alone can leave
  a biased answer looking stable. Vary every control the claim depends on, and where a mixture of
  tissues is pruned against averaged relaxation, compute per tissue instead. Say what converged and
  exclude what did not rather than quoting it.
- **A measured gradient moment is not an impossibility theorem.** A large crusher dephases the
  pathways you analysed under the timing you analysed. "No pathway can ever refocus" is a different
  and much stronger claim; it needs an ensemble that includes the RF imperfections you are implying
  it survives.
- **Legal under your `Opts` is not legal on a scanner.** Per-axis limits say nothing about the
  vector norm when axes ramp together, and an oblique prescription is rotated onto physical axes
  before it plays — recheck there, and do not treat a `1/√2` derating as a guarantee. ADC
  segmentation that satisfies PyPulseq's arithmetic is still subject to a particular interpreter's
  limits.
- **Steady state is a declared target, not a universal requirement.** A finite dummy train does not
  have to bring every tissue within some tolerance to be a sound acquisition; long-`T1`/`T2` species
  may stay visibly transient, and over-catalyzing to fix that can cost more time than the artifact
  is worth. State the criterion, what you measured against it, what it does to the image, and what
  it cost — then let that be the answer.
- **Suspect your own reproduction before you blame a library.** If an upstream call behaves
  strangely, reduce it to a minimal case, record the input and the measured output, and report that.
  Do not describe it as a PyPulseq or SeqCraft bug on one surprising observation, and do not silently
  work around it without saying you did.
- **Name the quantity you measured.** A summed sidelobe magnitude is not energy and not resolution.

When evidence contradicts the design, find out which stage is wrong — protocol choice, physical
design, event implementation, emission, simulation, or reconstruction — fix it, and revalidate the
artifact rather than the belief.

### When to stop adding evidence

Stop when the requested acquisition is complete and physically defensible, every claim you are
making has evidence behind it, and what you have not established is stated. Not when you run out
of things that could be measured — there is always another sweep.

How it was built does not enter into it. Shipped Modules, sequence-local Modules or functions, and
LogicBlocks composed directly are all finished work, and one acquisition will often use all three.
Reuse what exists where it fits, but never create an abstraction in order to be done.

```text
look at the shipped Modules first, before researching a design you may not need to make
read enough literature to justify the choices you are about to state, then stop reading
run the highest-value emitted-file checks first
add a specialist analysis for a specific risk or an unresolved discrepancy, not for completeness
```

A contradiction or a non-converging number is the one thing worth spending more time on.
Investigate it, or report it unresolved — never quietly drop it. And none of this licenses skipping
a physical check to finish sooner: the aim is effort proportionate to the request, with the rigour
intact.

## 5. Deliver

Give the user the authoring source, the generated `.seq`, the measurements you ran, the settings
and where each came from, and the parameters they can change.

**Where output goes.** Use the directory the user asked for. If they did not name one, create a
clearly named directory inside their current project — not inside the installed package, and not
inside a SeqCraft checkout's source tree — and tell them where you put it. Generated code, `.seq`
files, simulation artifacts and figures are the user's project's, never the package's.

**Evidence check.** Does each physical statement rest on evidence you actually produced? Report in
plain prose:

```text
validated             the claims made are supported by the evidence run
partially validated   some claims rest on evidence not yet run -- say which
unvalidated           the user declined independent validation, and none was performed --
                      built-in checks may still have run while the file was built
blocked               a physical or technical obstacle prevents the claim
```

Say which layers ran. "Partially validated" is an honest and useful result when simulation is
unavailable; **"compiles" is never validation**, and nothing here is scanner-ready. Gradient and
slew legality under chosen `Opts`, or a PNS prediction against an illustrative model, is a design
check — not RF/SAR, acoustic, vendor-interpreter or site approval. Say so when you hand it over.

An **unvalidated** hand-off still describes what it is: the acquisition you intended, the
independent checks that would normally run, the fact that none of them did, and any built-in checks
that ran while the file was built and written. A draft the user asked for is a legitimate
deliverable. A draft presented as though it had been checked is not.

For follow-up requests, continue from the established design. Do not reopen settled literature or
architecture decisions unless the new requirement actually disturbs them.

## When to stop and ask

```text
the physics genuinely forks and no default is defensible
the request is self-contradictory
real hardware deployment is requested and site limits are unknown
a requirement cannot be realised or validated with available capabilities
an expected API is absent from the installed version and no local fallback fits
```

In the last two cases, report the obstacle. Do not change the installed package, and do not quietly
deliver something weaker under the original name.
