# The physical-design scope

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-26, so that the pull request carries the direction behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../../README.md) for which one to edit.

**Status:** implemented. A **thin public declaration surface** sits above an internal substrate;
`seqcraft.design.scope` and `seqcraft.design.joint` stay private.
**Date:** 2026-09-26
**Supersedes:** the spike record `2026-09-26_repetition_family_scope_spike.md`

---

## 1. What it is

The region of a sequence whose physical degrees of freedom a designer may own.

```text
                         packaged kernel
                               |
                               | adapter
                               v
FlowComp / VENC ------> PhysicalDesignScope ------> shared physical designer
                               ^                              |
                               | adapter                      v
                               |                      timed LogicBlock
                      user composition                        |
                                                              v
                                                          compiler
```

`ScopeGeometry` + `AxisRequirement` in `design/scope.py`, over `design/joint.py`. The middle never
learns which side produced it.

### The two layers

```text
PUBLIC     seqcraft/physical_design.py
           PhysicalDesignScope     where a designer may work, in a composition you wrote
           design_repetition       design the family
           RepetitionDesign        .repetition(state), .at(te_s=...), .te_s, .window_s

INTERNAL   seqcraft/design/scope.py    ScopeGeometry, AxisRequirement, design_scope
           seqcraft/design/joint.py    Schedule, JointProblem, claims, realisation families
```

The public layer is a **translation**, not a second designer. It derives the two callables the
internal layer needs rather than asking for them: `fixed` by measuring the user's own blocks
where the schedule puts them, and `target` from the intent plus one default — *the designed
region brings `k` back to the origin at the echo on every axis it owns*.

That default is also the public surface's sharpest limit. A composition wanting a **non-zero** `k`
at the echo on a designed axis, as a Cartesian phase encode does, is not expressible; that is what
the packaged kernels are for. Adding it later is one optional field, not a redesign.

Nothing public names a `Schedule`, a `JointProblem`, a claim, a realisation family, an ownership
table or a raw `(m0, m1)` tuple. `design_scope(admissible=...)` stays internal: it is the future
local-PNS seam (§9), not a user constraint hook.

## 2. It is optional, and that is a hard principle

```text
events / Modules / functions  ->  LogicBlock  ->  compile
```

must stay valid with no scope anywhere, and it does. A scope says one thing:

> *for this part of my sequence, design the physics and the timing together before anything is
> materialised.*

It is **not** a required Module, a level of the sequence hierarchy, necessarily a class,
necessarily one literal TR, necessarily the whole sequence, or something the compiler knows about.
It is not a `LogicBlock` that already exists and is then patched: the adapter materialises once,
from a finished design.

**Not every useful composition has to become a packaged kernel to get physical design.**
`examples/gre_spiral_2d/03_flow_comp.ipynb` is the proof: a repetition with no kernel class reaches
the same designer `GRE2DTR` uses, through a function returning typed data.

## 3. Who the physical owner is

> The **physical owner** of a requirement is the smallest pre-materialization scope that sees all
> of the coupled requirements and all of the adjustable freedom they compete for.

Smaller than that and it cannot solve the coupling. Larger and it takes over waveforms it has no
reason to touch. For a Cartesian gradient echo the answer differs per axis, and a caller writing
`flow_comp=sc.FlowCompensation(axis=...)` does not need to know which:

```text
x   CartesianLine's own prephaser.  One component sees the whole problem, so it stays local
y   the phase-encode winder, which also carries the encode area -- two requirements, one
    waveform, so the repetition owns it
z   the slice rephasing.  Excitation realises it by default; when a scope claims z it asks for
    rephase=False and designs rephasing and compensation together, because on that axis they are
    one problem
```

That last line is the general mechanism, not a special case. `Excitation.build(rephase=False)`,
`CartesianLine(prephase=False)` and `SpiralReadout(prephase=False)` are the same seam: a component
**states** the moment it needs and hands the realisation to whoever can see more of the problem.

## 4. The adapter invariant

> Semantic instants, fixed contributions and adjustable regions supplied by an adapter must all
> refer to the same physical decomposition.

If a scope takes ownership of a component's prephaser, winder or rephaser, timing metadata that
component reports **for its own configuration** — one that still contains the region now owned
elsewhere — cannot be reused unchanged.

Found the hard way. `SpiralReadout.time_to_echo` is measured from a block including its own
prephaser; an adapter that sets `prephase=False` and designs that prephaser itself must subtract
`prephaser_duration_s`. Not doing so put the endpoint past the end of the arm and inside the
rewinder, where the `fixed` callable and the emitted block disagreed about what was playing. It
does not raise — it shows up as a residual, `7.3e-2` where `1e-13` was expected, which is why
`test_the_owned_region_must_be_taken_out_of_the_borrowed_timing` pins both numbers.

A second instance of the same class, found in this pass and present on `main` before it: fixed
contributions must be integrated **from the semantic origin**, not from the start of the block. On `x` and `y` those agree. On `z`
half the selection lobe plays before the RF centre and dephases nothing, because there is no
transverse magnetisation yet; counting it nulls the wrong interval and leaves `k_z` off zero while
reporting `m1 = 0`.

`GRE3DTR` had the same bug for a **selective slab** and was not migrated in this pass -- it was
fixed in place, because it is one line and the alternative was leaving known-wrong physics on
`main` to preserve scope discipline. Non-selective was always right, which is why nothing caught
it: with no slab lobe the two integration starts agree.

## 5. `design_states` is a proposal

```text
design_states  ->  propose a candidate schedule
               ->  verify across the supported family
               ->  widen if necessary
```

An adapter may choose it by analytic bound, signed extrema, corners of a multidimensional state
space, a representative set, or full enumeration. **Correctness never depends on the choice being
right.** A poor proposal costs search time; verification is what keeps it from shipping an
infeasible state, and it is needed because the feasible set has small raster holes — "needs the
most area" does not imply "the window that serves it serves everything".

Two justified proposals in tree:

```text
Cartesian PE   the two signed ENDS.  The requirement is linear and signed in the index, so the
               ends bound every index between them -- and both ends rather than the largest |k|,
               because a signed fixed contribution makes one end harder than the other

spiral         the axis-aligned angles.  The arm is one trajectory rotated, so on x the
               requirement is a state-independent magnitude times cos(angle); the bound is
               attained at angle 0 on x and a quarter turn on y
```

Family-level timing versus state-level realisation is the contract: the search runs over the
proposal, and every state is realised at the one schedule that comes out. Avoiding repeated
*timing search* is the point, not avoiding per-state realisation.

## 6. Re-realisation, not stretching

`realise_scope_at` realises a whole family at a schedule someone else chose. A protocol holding
two families takes the longer echo time and asks both to design again at it — no acquisition-wide
timing manager. Going back through realisation is the point: every moment is measured to the echo,
and the echo moved.

## 7. Current representational limits

Real, and refused rather than worked around:

```text
one adjustable window    a family needing two independently adjustable regions -- a prephaser and
                         a rewinder designed together -- does not fit, because two scopes cannot
                         share one schedule.  Such a family should stay outside this substrate
                         until evidence justifies the extension

a constant tail          what plays between the window and the echo may not depend on the state.
                         A rotated spiral arm is fine; an arm whose length varies per state is not

m0 and m1 only           design/joint.py's scope.  m2 needs a family that computes its own exact
                         response and carries enough independent freedom
```

These are limits of the current representation, not definitions of physical design.

## 8. Edge boundaries

```text
external gradient / RF history   a scope guarantees physics only from the contributions it is
                                 given.  The two principled futures are to ENLARGE the scope, or
                                 to pass the contribution in explicitly as part of `fixed`.
                                 Nothing inspects the surrounding sequence

pathway-aware moments            deferred until a real sequence requires it.  No pathway model

MRF and other irregular          may simply bypass this abstraction.  A sequence with no useful
sequences                        fixed repetition structure is not a counterexample to it

multi-echo                       a scope names ONE endpoint.  An axis with inter-echo gradients
                                 needs the checkpoint model, which has no place here yet
```

## 9. The PNS seam

PNS is **different in kind** from `FlowCompensation` and `VelocityEncoding`, and must not become
another intent of the same type:

```text
FlowComp / VENC   define WHAT physical relation must hold
PNS               constrains WHICH realised waveforms and timings are admissible, and has
                  temporal history across repetition boundaries
```

**Level 1, local, the seam that exists today.** `design_scope(admissible=...)` — documented,
unused, no implementation. It receives a **fully realised candidate**, not a schedule, because
admissibility depends on the realised gradient:

```text
candidate schedule -> candidate realisation -> optional local admissibility -> accept / reject
```

One realisation family may fail at a schedule where another passes, and rejecting a candidate
leaves longer schedules free to be tried. Two local contexts are worth keeping reachable — zero
history, and periodic steady state for a homogeneous repeated family — and both are computable
from a `ScopeDesign` alone. Arbitrary incoming filter state is not representable, and should stay
that way.

**Level 2, global, later and authoritative.** The ordered whole-sequence waveform, evaluated
continuously.

```text
local design -> assemble -> global check
    if over limit: return overage + provenance / a conservative derating estimate
                   -> feed a stricter design profile back to the affected scope
                   -> redesign ONCE -> rebuild -> re-check
```

Deliberately not an iterative optimizer. If the corrected design passes, no claim of globally
minimal TE/TR. If it still fails, report and let the caller choose a stronger policy.

**Never time-stretch an already-realised waveform by `1 / overage`.** Feedback goes back through
realisation so M0/M1, the VENC relation, echo alignment and timing are solved together again. A
PNS-derived derating factor is an estimate entering as a design profile, not a multiplier on a
finished shape. A future evaluator should be able to name the dominant scope, axis or interval and
redesign only that, rather than derating the whole acquisition.

## 10. The compiler boundary

Same-axis gradient overlap is **not** intrinsically illegal. The designer reasons about the
composite waveform inside its declared scope — summed `G(t)`, summed `dG/dt`, the M0/M1 targets,
the hardware limits — and the compiler stays authoritative for the whole tree.

```text
physical designer   waveform and timing feasibility inside the declared scope
LogicBlock          final placement representation
compiler            whole-tree lowering, split / superpose / emit, final legality backstop
```

Local design passing while the compiler later refuses is **expected** when a caller overlays
something outside the declared scope. The answer is diagnostics that make the provenance
readable, not moving global compilation into the designer. Both notebooks emit same-axis merge
warnings for exactly this reason: the designed winder merging with the leaf beside it is the
compiler's business to lower, and it does.
