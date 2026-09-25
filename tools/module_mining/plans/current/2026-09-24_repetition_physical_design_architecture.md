# SeqCraft Repetition Physical Design Architecture

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-24, so that the pull request carries the direction behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../../README.md) for which one to edit.

**Status:** architecture decision and implementation reference  
**Project:** SeqCraft  
**Date:** 2026-09-24  
**Scope:** repetition-level joint waveform/timing design for reusable augmentations such as velocity encoding and flow compensation

---

## 0. Project scope

> SeqCraft is a high-level MRI sequence authoring framework that converts physical acquisition
> intent into hardware-aware, fully timed Pulseq sequences. It owns reusable MRI physical design
> needed to make that translation practical — including coupled waveform realization, minimum
> timing, and whole-sequence validation — while scanner-specific safety models remain optional
> evaluators rather than core sequence semantics.
>
> SeqCraft aims to produce realistic scanner-ready research sequences, but does not replace vendor
> safety validation or certify scan safety.

```text
scanner-runnable / physically realistic     a goal
scanner-certified / guaranteed safe         NOT a SeqCraft claim
```

**Hide complexity; do not deny complexity.** MRI sequence design is intrinsically complicated. A
framework that refuses to handle timing coupling, hardware feasibility or flow compensation does
not remove that work — it returns it to every user's script. A framework that exposes all of it
becomes harder to use than Pulseq.

So the internal physical-design machinery may become sophisticated, provided ordinary authoring
stays close to

```python
GRE3D(..., flow_comp=..., velocity_encode=..., te_s=None, tr_s=None)
```

and never requires `DifferenceClaim`, `CommonModeClaim`, `JointProblem`, `Schedule`, waveform
families, PNS internals or solver selection.

---

## 1. Decision summary

SeqCraft should separate three responsibilities that all involve time but are not the same problem:

```text
MRI physical / repetition design
    decides what waveform and timing are physically feasible
    decides minimum required windows, TE, echo spacing and TR
        ↓
LogicBlock
    records the already-decided event placement
        ↓
compiler
    lowers the timed tree into legal Pulseq blocks
    split / superpose / emit / validate legality
```

The new joint machinery is **not a new timing layer between LogicBlock and the compiler**. It lives in the existing Module/kernel layer and is used only when a physical problem is coupled enough that local module design is no longer sufficient.

The governing rule is:

```text
self-contained physical problem
    → local module design

cross-component / cross-state coupled physical problem
    → repetition-level shared design

both materialize final events only after the physical solve
    → LogicBlock
    → compiler
```

The shared design composes reusable physical requirements **before** waveform realization. It must not patch already-built LogicBlocks or rewrite private events owned by another module.

---

## 2. Why this architecture is needed

The original problem is the `N sequences × M augmentations` growth.

A bad long-term architecture would require each combination to be implemented independently:

```text
GRE2DTR × VelocityEncode
GRE2DTR × FlowComp
GRE3DTR × VelocityEncode
GRE3DTR × FlowComp
EPI     × VelocityEncode
EPI     × FlowComp
...
```

The desired maintenance shape is approximately:

```text
N repetition/base-sequence adapters
+
M reusable augmentation requirement contributors
+
shared realization machinery
```

A compatible new sequence should not rederive VENC or flow-compensation physics. A new augmentation should not require editing every existing compatible sequence.

---

## 3. Existing SeqCraft already contains the first version of this idea

This work does not introduce physical timing design from nothing.

`GRE2DTR` already computes:

```text
readout prephaser minimum
phase-encode minimum
excitation rephaser minimum
        ↓
common winder duration
        ↓
minimum TE
        ↓
requested TE fill
        ↓
minimum TR
```

`GRE3DTR` already goes further:

```text
slab rephase requirement
+
partition-dependent kz requirement
        ↓
enumerate partitions
        ↓
find limiting state
        ↓
one common z-winder duration
        ↓
constant TE across partitions
```

That is already a small repetition physical-design problem. The shared joint machinery generalizes this pattern only where coupling justifies it.

---

## 4. Timing ownership

Three meanings of timing must stay distinct.

### 4.1 MRI physical timing design

This layer answers:

```text
What is the minimum legal winder duration?
What is the minimum TE?
What is the minimum echo spacing?
What is the minimum TR?
If FlowComp or VENC is enabled, how much must those minima move?
Which state is limiting?
```

This belongs to the Module/kernel layer because the answer depends on MRI physics.

### 4.2 LogicBlock scheduling

This layer records concrete placement:

```text
RF begins at t = ...
gradient begins at t = ...
these blocks overlap
ADC centre is at t = ...
```

LogicBlock **carries a timing decision**. It does not decide why a TE is feasible or minimal.

### 4.3 Compiler timing/lowering

This layer answers:

```text
Where are Pulseq block boundaries?
How are overlapping gradients split and superposed?
How is the already-timed tree emitted legally?
```

The compiler may validate legality, but it must not redesign MRI protocol timing.

---

## 5. Core rule: take over coupled degrees of freedom before materialization

The shared designer must not work like this:

```text
module builds waveform
    ↓
shared solver inspects it
    ↓
shared solver overrides / patches it
```

Instead:

```text
module exposes physical requirement / fixed fact
        ↓
repetition owner detects a coupled problem
        ↓
shared designer owns the coupled degrees of freedom
        ↓
final waveform is materialized once
```

Existing precedent: selective excitation inside `GRE3DTR`. The excitation can normally realize its own slice rephaser, but `GRE3DTR` combines slab rephasing with partition encoding before materialization rather than patching a finished excitation block.

---

## 6. Local path versus joint path

### 6.1 Local path

Simple modules keep local design when the physical problem is self-contained.

Examples:

```text
PhaseEncode
    designs its own trapezoid

CartesianLine
    designs its ordinary prephaser

Excitation
    designs its ordinary rephaser

VelocityEncode
    can still emit its standalone canonical bipolar
```

Do not rewrite simple closed-form modules into a generic constraint engine merely for architectural uniformity.

### 6.2 Joint path

When several requirements constrain the same waveform/timing freedom:

```text
PhaseEncode M0 requirement
+
VelocityEncode ΔM1 requirement
+
FlowComp common-mode M1 requirement
+
fixed readout / excitation / wave contributions
+
hardware limits
        ↓
repetition-level joint physical design
```

The relevant local realizations are not built first. The shared designer realizes the coupled waveform, and the repetition kernel incorporates the required window into its own TE/TR policy.

---

## 7. Reusable requirement semantics

The current evidence supports a minimal semantic model built around a named interval and a physical target.

Conceptually:

```text
MomentTarget(
    axis,
    order,
    origin,
    endpoint,
    value/component,
)
```

Both origin and endpoint matter. For a shift of time origin:

```text
m1' = m1 - Δt * m0
```

Therefore `m1 = 0 at the echo` is only origin-independent when `m0 = 0`.

For ordinary GRE single-path analysis, the repetition supplies the effective excitation instant as the origin and the target echo as the endpoint. The joint designer must not hard-code `Excitation.time_to_center()`; a future RF family may supply a different physically meaningful semantic instant.

Pathway-aware moments remain deferred.

---

## 8. Velocity encoding and flow compensation compose at the requirement layer

For a two-state velocity acquisition:

```text
VelocityEncode:
    m1(+) - m1(-) = Δm1

FlowComp:
    (m1(+) + m1(-)) / 2 = 0
```

with:

```text
Δm1 = 1 / (2 * venc)
```

Composing them derives:

```text
m1(+) = +Δm1/2
m1(-) = -Δm1/2
```

The VENC relation remains single-source in `VelocityEncode`.

The waveform realization layer receives only the resolved absolute target for one concrete repetition state. It must not know whether that target came from PE, VENC, FlowComp or another acquisition policy.

A second claim on the same semantic component should refuse rather than silently override another claim.

---

## 9. “Generic joint solver” means a generic boundary, not one universal algorithm

The generic part is the **physical design boundary and orchestration**, not one solver that eventually handles every MRI design problem.

Use this rule:

```text
new need still has:
    linear gradient moments
    + basis coefficients
        → extend/use the moment realization family

new need changes the mathematical problem:
    b-value / b-tensor
    eddy-current convolution
    PNS constraints
    concomitant-field / Maxwell constraints
    arbitrary fixed endpoints
    genuinely nonlinear objective
        → use a different realization backend

new need only contributes another physical target
        → add/compose a requirement contributor
        → do not modify the waveform solver
```

Do not build a universal `BaseSolver` inheritance hierarchy now.

---

## 10. Preferred implementation style: function composition, typed contracts

Prefer compiler-like functional stages over an inheritance tree.

Conceptually:

```text
requirements/state
    ↓
resolve claims
    ↓
resolved physical targets
    ↓
construct repetition physical problem
    ↓
choose schedule candidate
    ↓
try realization backend/family
    ↓
feasibility + utilization + events
    ↓
final repetition timing
```

Possible internal pieces:

```text
resolve_claims(...)
JointProblem / ResolvedTargets
ScheduleCandidate
RealizationAttempt
JointDesignResult

realize_moment_basis(...)
realize_hb_m01(...)
optional realize_with_gropt(...)
```

Only introduce a small `Protocol` for interchangeable realization objects if two genuinely different backends later require runtime interchange.

---

## 11. Current moment realization scope

The current spike proves a useful M0/M1 realization but should not imply arbitrary moment-order support.

For a symmetric lobe:

```text
M0 = area
M1 = area × centroid
```

is exact. The same shortcut is not sufficient for `M2+`, where the waveform's internal higher moments matter.

Therefore the first internal moment backend should state its actual scope: **M0/M1**.

If M2 is added later, compute the exact response of the chosen basis and provide enough independent degrees of freedom.

---

## 12. Realization families

The joint layer should not decide one universal waveform shape.

A realization family owns:

```text
what waveform degrees of freedom exist
how their exact moments are computed
how hardware feasibility is evaluated
how events are materialized
```

Useful current families include:

```text
single-lobe M0 realization
two-degree-of-freedom M0/M1 realization
base-lobe + zero-area bipolar M1 correction
small raster-aligned shape family
```

`wave-gre-flow-comp` shows why this boundary matters: its base component carries M0, a zero-area bipolar carries pure M1, and an internal raster-aligned shape parameter may be searched to improve gradient/slew utilization.

That shape search is a realization-family detail, not augmentation semantics.

---

## 13. Do not hard-code one duration-search algorithm into the generic layer

Neither generic bisection nor generic upward raster search is universally justified.

Bisection needs a valid monotonicity argument for the specific feasibility predicate.

Simple upward raster search is robust and understandable, but a narrow waveform family may show artificial feasibility cliffs.

Instead separate:

```text
schedule / duration candidate
        ↓
realization family at that candidate
        ↓
best feasible waveform + utilization
```

A realization family may use closed form, small parameter enumeration, raster search, proof-backed bisection, or a small numerical optimization depending on its mathematics.

---

## 14. Two-level design: schedule outside, waveform realization inside

The important architectural conclusion is that the problem is not only waveform search.

Product-like sequence behavior requires the timing schedule itself to move when a feature makes the old schedule infeasible.

```text
OUTER: repetition schedule design
    TE
    echo spacing / ΔTE
    TR
    adjustable window durations

        ↓ for one candidate schedule

INNER: waveform realization
    fixed contributions
    physical targets
    analytic realization family
    richer low-dimensional family
    optional numerical backend

        ↓

    feasible / infeasible
    gradient utilization
    slew utilization
    diagnostics
```

The outer design searches for the earliest feasible repetition schedule allowed by the user's timing policy.

This directly addresses the failure mode seen in fixed-TE `wave-gre-flow-comp`: a slightly longer TE or echo spacing may restore a short, healthy waveform instead of forcing a pathological tens-of-milliseconds realization.

---

## 15. Product-like TE/TR semantics

Preserve SeqCraft's useful distinction between automatic timing and explicit timing.

Conceptually:

```text
te_s=None
    → TE is AUTO / minimum feasible after all enabled physics

te_s=<number>
    → hard requested TE
    → refuse if enabled physics makes it infeasible
```

Likewise for TR and future echo-spacing controls.

Enabling FlowComp/VENC may therefore change:

```text
min_te_s
min_echo_spacing_s
min_tr_s
```

when those values are automatic.

It must not silently change an explicitly requested timing.

Example:

```text
feature off:
    min TE = 10.00 ms

feature on:
    min TE = 10.37 ms
```

Then:

```text
te_s=None
    → achieved TE becomes 10.37 ms

te_s=10.00 ms
    → ConfigurationError
      requested TE is shorter than new minimum
```

---

## 16. Repetition timing policy remains owned by the kernel

The shared realization machinery answers:

```text
For this coupled physical problem,
what waveform can be realized,
and how large must this local window be?
```

The repetition kernel still answers:

```text
Given all x/y/z windows and RF/readout constraints,
what are min TE, echo spacing and TR?
```

Ownership therefore remains:

```text
realization family
    owns local waveform feasibility

shared joint design
    owns coupled physical requirements

kernel
    owns repetition timing policy

LogicBlock
    owns concrete placement representation

compiler
    owns Pulseq lowering/legality
```

---

## 17. No feedback loop through LogicBlock/compiler

Reject this architecture:

```text
Module
    ↓
LogicBlock
    ↓
compiler discovers waveform does not fit
    ↑
asks Module to increase TE
    ↑
rebuild
```

MRI physical feasibility must be settled before final LogicBlock materialization.

The compiler may catch illegal output as a backstop, but it must not optimize protocol timing.

---

## 18. Multi-echo

`wave-gre-flow-comp` suggests that multi-echo does not necessarily require one train-wide optimizer.

A useful structure is checkpoint-to-checkpoint propagation:

```text
echo n target satisfied
    ↓
propagate moment state over ESP
    ↓
add fixed waveform contribution in this interval
    ↓
solve adjustable inter-echo window
    ↓
echo n+1 target satisfied
```

For a maintained non-zero M0:

```text
M1 changes by -M0 * ESP
```

This can remain a sequence of small local solves. A global optimizer remains possible if a future problem truly couples all windows simultaneously, but multi-echo alone does not establish that need.

---

## 19. Optional numerical backend / GrOpt

SeqCraft currently has **no implemented GrOpt adapter or optimizer hook**.

Previous work only established:

```text
GrOpt is a legitimate possible future realization backend.
It must not define SeqCraft's public physical semantics.
```

The architecture should preserve that possibility.

A future realization cascade may be:

```text
candidate schedule
    ↓
cheap canonical analytic family
    ↓ if unavailable / infeasible
richer low-dimensional family
    ↓ if still infeasible and justified
optional GrOpt adapter
```

GrOpt remains an optional external backend, not a core dependency. Its GPL licensing/distribution implications must be reviewed before an adapter is shipped.

Do not implement a universal optimizer interface now merely to reserve a seat for GrOpt.

---

## 20. Waveform feasibility cliffs

When small hardware-limit changes produce a jump from a few milliseconds to tens of milliseconds, distinguish:

```text
true physical feasibility cliff
```

from:

```text
artificial realization-family cliff
```

A narrow preselected waveform family may fail over a broad duration interval even though a slightly richer shape could succeed close to the original schedule.

Therefore fixed-duration evaluation should preferably report:

```text
feasible?
peak gradient utilization
peak slew utilization
limiting constraint
chosen realization family/shape
```

rather than only a boolean.

### 20.1 What the adversarial stress matrix found

A deliberate search for a cliff over the current architecture — gradient and slew limits, M0 and
M1 targets, fixed contributions, the origin-to-window lead, window duration, explicit and AUTO TE,
fixed and AUTO echo spacing, acquisition state, axis and echo count — with local refinement at
every large coarse transition. The evidence is in the spike findings, §7.

**Extensive stress testing did not reveal a large artificial feasibility cliff.** This is a
statement about a sampled and refined parameter space, not a proof that cliffs cannot occur; the
families here are low-dimensional by design, and a richer requirement could produce one.

Reproduce it with `tools/module_mining/stress_repetition_design.py` (`--quick` or `--full`), which
does not run in CI.

Three things it did find, all worth carrying into the design:

```text
a family artefact               a raster-aligned two-lobe split does not exist for every
                                window, so the feasible set has holes -- 1.6 % of cells,
                                widest OBSERVED 110 us.  A sampled maximum, not a bound,
                                and a property of this family: the raster is fundamental,
                                this particular missing split is not, and a richer
                                raster-aligned family could plausibly fill it.

a turning point in the lead     pushing the origin further before the window helps until the
                                M0 target pins the area and its own M1 overshoots.  A smooth V.
                                "A longer TE fill always helps" is false.

a diagnostic on the wrong       `utilisation` sampled 48 fixed fractions while the emitter
lattice                         searches the gradient raster, overstating u by up to x1.82.
                                Corrected.  The metric and the emitter must search the same
                                lattice, or the diagnostic manufactures the cliff it exists
                                to rule out.
```

The third is the second time a measurement rather than a design produced a false cliff, which is
the argument for reporting the five fields above rather than a boolean: a boolean cannot be caught
being wrong by 82 per cent.

---

## 21. Baseline sequences do not need immediate migration

Do not rewrite ordinary FLASH/GRE/MPRAGE timing merely for consistency.

Use the principle:

```text
simple physical problem
    → local analytic implementation is good

coupled reusable problem
    → shared repetition design

genuinely complex constrained problem
    → numerical backend when evidence requires it
```

A baseline kernel may later migrate to shared machinery if that clearly simplifies the implementation or improves correctness, but migration is not an architectural requirement.

---

## 22. PR #39 disposition

PR #39 remains useful evidence and should remain open during the architecture integration.

It establishes:

```text
CartesianLine
    readout axis
    first echo
    M0 = 0
    M1 = 0
    local equal-duration two-winder realization
```

Its realization and tests should be preserved.

**Decided 2026-09-25: the selector is private.** `_null_moment_order` is internal, and the public
way to ask is at the repetition:

```python
GRE2DTR(..., flow_comp=sc.FlowCompensation(axis='x'))
```

which routes to this solve. A public per-leaf keyword would have been the first instance of a
leaf-by-augmentation surface with no principled place to stop.

The evidence increasingly suggests the solve is useful as an internal realization primitive, while the public per-leaf option may represent the N×M API shape the joint architecture is intended to avoid.

Do not remove or preserve the public API until the real-kernel integration demonstrates the final user-facing path.

---

## 22b. Sequence envelope versus constituent design profile

One `Opts` conflates two different things, and the stress case separates them:

```text
physical scanner capability
        ↓  operational margin
sequence admissibility envelope      the final summed waveform must fit this
        ↓
constituent design profiles          what a given component is designed under
```

The reference implementation declares a broad sequence envelope and then designs different
constituents under different, tighter policies — gentle readouts may use a larger fraction of the
envelope, while short aggressive transitions (blips, ramps, prephasers, spoilers, flow-compensation
lobes) are held well below it to limit PNS and gradient stress.

Measured from the admitted stress case at `0e1ec51`. **These numbers are test input, not a
proposed API** — they are one author's derating choices for one sequence, recorded so the stress
harness can reproduce the conditions that motivated this section:

```text
physical                       80 mT/m     200 T/m/s
sequence envelope   x0.90/x0.70   72         140
lowPNS  (most parts) x0.90/x0.41  72          82
lowPNS2 (FC module)  x0.60/x0.35  48          70
```

**The durable concepts are the envelope and the profile.** Names like `sys`, `sys_lowPNS`,
`sys_lowPNS2` are that implementation's spelling, not an abstraction to reuse. SeqCraft ships no
named profile set, no derating factors and no table like the one above; a profile is whatever
`Opts` a component is designed under, and the envelope is whatever the compiler validates
against.

Ownership:

```text
compiler                validates the final summed waveform against the sequence envelope
Module / repetition     chooses and satisfies the design profile for what it realizes
```

For Stage C this separation stays **internal**. No public scanner/profile framework yet.

## 22c. Safety evaluation is a separate whole-waveform concern

```text
Design envelopes are policies.
Safety evaluation is a separate whole-waveform concern.
```

A design profile is a static heuristic. A real PNS model has **temporal memory**, so it is not an
instantaneous slew check and it cannot be satisfied by a per-lobe rule. The long-term seam:

```text
initial repetition design
        ↓ assemble the complete gradient waveform
optional whole-waveform PNS evaluation
        ↓ acceptable?  yes -> accept
          no -> identify the dominant regions / axes
                redesign those degrees of freedom
                (lower slew, wider window, different shape family,
                 or a longer AUTO TE / ESP / TR)
                and re-evaluate the complete waveform
```

Iterative feedback, not two fixed passes, and **not implemented in Stage C** — only the seam and
the ownership are being settled here.

### The parameter interface, and what it is not

A SAFE-style model needs a compact per-axis parameter set, not a confidential vendor file:

```text
per axis    tau1, tau2, tau3    a1, a2, a3    stim_limit    stim_thresh    g_scale
```

The core model uses `tau1..3`, `a1..3`, `stim_limit` and `g_scale`; `stim_thresh` rides along
because a vendor file carries it and a caller supplying parameters by hand will have it, not
because the model needs it.

so an eventual interface can support either a **local importer** (parse the vendor file privately,
keep only the compact parameters, never commit the file) or **user-supplied parameters directly**.
Derived parameters are not automatically non-confidential; users remain responsible for their own
vendor restrictions.

Three modes, each labelled honestly:

```text
no model              static conservative profiles; no scanner-specific number reported
synthetic reference   development, CI and waveform comparison; NOT scanner-specific
user parameters       a scanner-specific estimate; still NOT a certification
```

**Any such integration is a prediction, never a guarantee.**

## 22d. RF scope: peak B1 is a hard limit, and is not SAR

```text
peak B1               a hard instantaneous RF/hardware constraint -- already enforced,
                      per RF path, with a compiler backstop
RF energy proxy       a possible future reporting metric, e.g. integral |B1(t)|^2 dt,
                      for comparing candidate designs.  A burden proxy, not SAR
SAR                   no claim, unless and until a real SAR model exists
```

Stage C does not extend into SAR modelling.

## 23. Pathway-aware moments

Explicitly deferred.

This architecture uses the intended single sequence pathway and semantic intervals.

Future work order:

```text
pathway-aware M0 / k-space bookkeeping first
pathway-aware M1 later
```

Do not introduce EPG/coherence-pathway machinery as part of Stage C.

---

## 24. Implementation plan from the current spike

### Step A — tighten the current prototype boundary

Refactor provisional joint code so that requirement resolution/state semantics are separate from M0/M1 realization and schedule search.

Do not imply arbitrary moment-order support.

Prefer typed data plus function composition over solver inheritance.

### Step B — make fixed contributions state-aware

The physical problem must support fixed contributions that vary with state and semantic interval.

This is required by real sequence cases such as state-/echo-dependent wave-gradient contributions.

### Step C — separate schedule candidates from waveform realization

Do not keep the old echo endpoint fixed while making the winder longer.

For every candidate schedule:

```text
derive actual semantic times
re-evaluate fixed moments at those times
solve the waveform
check hardware
```

Timing and waveform feasibility are coupled.

### Step D — compare realization families near the feasibility boundary

At fixed candidate timing, compare at least:

```text
simple adjacent two-lobe M0/M1
decoupled base + zero-area bipolar M1
small raster-aligned shape search
```

Measure:

```text
feasibility
gradient utilization
slew utilization
```

and perturb hardware limits around the boundary.

Use this to choose the first robust internal realization cascade.

### Step E — keep optional optimizer fallback conceptual unless needed

If all real cases remain well served by analytic/low-dimensional families, do not add GrOpt yet.

If a concrete case demonstrates a material failure that GrOpt solves, write the smallest concrete adapter.

### Step F — rewire real kernels

Integrate the shared design into actual:

```text
GRE2DTR
GRE3DTR
```

without augmentation-specific physical derivations inside the kernels.

The augmentation/acquisition-state layer supplies reusable claims.

The kernel supplies:

```text
base state
semantic origin / endpoint
fixed contributions
available timing windows
repetition timing policy
hardware
```

### Step G — verify product-like timing behavior

Exercise:

```text
AUTO TE / TR
explicit TE at new minimum
explicit TE between old and new minimum → refuse
explicit TR at new minimum
previously legal TR made too short by feature → refuse
representative and limiting repetition states
```

Validate moments on the **complete emitted repetition**.

### Step H — decide PR #39 public surface

**Done 2026-09-25: internal.** The solve stays, the selector became private, and the repetition-level
intent is the public way to reach it. Re-exposing it later, if a standalone-readout need is ever
demonstrated, is a one-line non-breaking change; retiring a released keyword would not have been.

### Step I — finish Stage C before bSSFPTR

Do not start bSSFPTR until this repetition-level boundary is stable.

---

## 25. Success criteria

The architecture is successful if:

```text
VelocityEncode physics is defined once.
FlowComp semantics are defined once.

A new compatible sequence does not reimplement existing augmentations.
A new augmentation does not require editing every compatible sequence.

Local simple modules remain simple.

Coupled degrees of freedom are designed before LogicBlock materialization.

AUTO timing recomputes minimum feasible TE/ESP/TR.
Explicit requested timing refuses when infeasible.

LogicBlock remains a thin timed-event tree.
Compiler remains responsible for Pulseq lowering and legality only.

The final emitted repetition is independently measured for correctness.

A numerical optimizer remains optional until a concrete case genuinely needs it.
```

---

## 26. Final mental model

```text
Protocol / acquisition intent
        ↓
Module / kernel physical design
        │
        ├─ local analytic design
        │       for self-contained cases
        │
        └─ shared repetition design
                for coupled cases
                │
                ├─ resolve physical requirements
                ├─ choose/search repetition timing
                ├─ try waveform realization families
                └─ determine minimum feasible windows
        ↓
Final timed events
        ↓
LogicBlock tree
        ↓
Compiler
    expand
    split / merge / superpose
    Pulseq legality
        ↓
pypulseq.Sequence
```

The shortest summary is:

> **The Module/kernel layer decides what MRI waveform and timing are physically feasible; LogicBlock records that decision; the compiler lowers it. Simple modules keep local timing design, while cross-module or cross-state coupled degrees of freedom are handed to a shared repetition designer before materialization.**
