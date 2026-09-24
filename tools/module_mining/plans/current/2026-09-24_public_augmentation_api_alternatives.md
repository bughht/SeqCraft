# Public augmentation API — alternatives for review

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-24, so that the pull request carries the direction behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../../README.md) for which one to edit.

**Status:** proposal for review. **Nothing here is implemented.**
**Date:** 2026-09-24
**Governed by:** [`2026-09-24_repetition_physical_design_architecture.md`](2026-09-24_repetition_physical_design_architecture.md) §0

`joint_claims` and `encoding_states` are the spike's internal representation. They must not ship.
What must not appear in any public signature:

```text
CommonModeClaim   DifferenceClaim   JointProblem   Schedule
composite state keys   waveform families   solver selection
```

The target, from the scope statement:

```python
GRE3D(..., flow_comp=..., velocity_encode=..., te_s=None, tr_s=None)
```

---

## A — keyword intent on the kernel

```python
GRE2D(
    opts=opts, fov_mm=220.0, matrix=(64, 64), thickness_mm=5.0, flip_deg=15.0,
    flow_comp='y',                 # or ('y', 'z'), or True meaning "every axis I own"
    velocity_encode=(1.5, 'z'),    # venc in m/s, and the axis
)
```

The acquisition layer then yields two encoding states per line, and reconstruction subtracts them.

```text
+ smallest surface; reads like a protocol card
+ no new public type at all
- two features, two spellings, and a third would want a third keyword
- `flow_comp=True` has to mean something, and "every axis I own" is a policy the user cannot see
- awkward the moment an augmentation needs more than one number
```

## B — a small intent object per augmentation

```python
GRE2D(..., augmentations=[
    sc.FlowCompensated(axis='y'),
    sc.VelocityEncoded(venc_m_s=1.5, axis='z'),
])
```

Each object carries physical intent only — no waveform, no claim vocabulary. Internally it
contributes the claims the resolver already understands.

```text
+ one extensible spelling; a third augmentation adds a type, not a keyword
+ the object is where per-augmentation documentation and refusals live
+ `VelocityEncoded` is intent, so nothing is realised until the kernel asks
- a new public name per augmentation
- "augmentations=[...]" invites the plugin bus the architecture rules out; the list must stay
  a closed, documented set rather than an extension point
```

## C — intent on the acquisition, not the repetition

```python
GRE2D(..., flow_comp='y')
scan = sc.modules.PhaseContrast2D(gre, venc_m_s=1.5, axis='z')
```

Velocity encoding is, unlike flow compensation, also an **acquisition** decision: it doubles the
scan and defines the subtraction. This puts the state-generating augmentation where the states are
scheduled and leaves the purely-per-repetition one on the repetition.

```text
+ puts the thing that changes the acquisition at the acquisition layer, which is where
  ordering, reconstruction and the number of shots already live
+ the repetition keeps a single spelling for augmentations that do not add states
- two places to look
- needs the acquisition to pass state down into build(), which today it does not
```

**Recommendation: B for the repetition-level surface**, with C's observation kept in mind — a
velocity acquisition is the thing that knows there are two states, so whatever generates them
belongs beside the line ordering rather than inside one repetition.

The stress pass supplied an argument that was not available when these three were written. The
spike accepts a claim on an axis the repetition does not own and either silently does nothing or
refuses with the wrong reason (see **Capability and refusal** below). Whatever ships has to
validate the axis, and the three alternatives differ in whether there is anywhere to put that:

```text
A  flow_comp='q'                      a bare string; the only place to check it is inside the
                                      kernel, mixed in with the physics
B  sc.FlowCompensated(axis='q')       the object validates its own axis against the repetition
                                      it is handed to, and owns the message that explains it
C  PhaseContrast2D(gre, axis='q')     same as B for the state-generating case, and nothing for
                                      the per-repetition one
```

B is the only one of the three with an obvious home for the refusal, which is the deciding
argument rather than a tie-breaker: this is a class of bug the stress matrix caught once already.

**Still open for the reviewer:** whether `augmentations=[...]` is worth the plugin-bus risk noted
above, or whether the closed set should be named keywords taking intent objects
(`flow_comp=sc.FlowCompensated('y')`), which keeps the refusal home and removes the list.

---

## `VelocityEncode`: eager realisation, and whether it matters

The constructor designs its standalone bipolar eagerly. Swept over 648 systems — 3 to 72 mT/m,
6 to 200 T/m/s, venc 5.0 down to 0.01 m/s — it **never refused**:

```text
builds      648 / 648
pair duration   0.600 ms  (72 mT/m, 200 T/m/s, venc 5.0)
             to 28.5 ms   (3 mT/m, 6 T/m/s, venc 0.01)
```

`_solve_times` solves the continuous hardware-limited problem in closed form and rounds up onto
the raster, so the lobe grows until it fits rather than failing. Eager realisation is therefore a
**cost**, not a correctness blocker: it designs a waveform the joint path may never emit. It is
also not a capability gap in the other direction — there is no venc the joint path can reach that
the standalone pair refuses, because the standalone pair does not refuse.

```text
A  make the standalone realisation lazy      works, but a property that builds on first
                                             access is a surprise in a module whose other
                                             attributes are all cheap
B  separate specification from realisation   already true: `VelocityEncode.delta_m1_for(venc)`
                                             is the relation with no waveform and no instance
```

**So no refactor is needed.** The public intent should carry `venc_m_s` and an axis, and the
repetition should read the relation through `delta_m1_for`. Constructing a `VelocityEncode` is
then what a caller does when they want the *standalone bipolar* — which remains a real and useful
thing — not a prerequisite for the joint path. The VENC relation is not duplicated either way.

---

## Capability and refusal

An augmentation must not be assumed to work on every axis of every sequence. The repetition owns
the answer, and it is a capability question rather than a class-name question:

```text
does this repetition own an adjustable window on that axis,
between the semantic origin and the endpoint?
    yes  -> it can carry a moment claim there
    no   -> refuse, naming the axis and what it would need
```

`GRE2DTR` owns `y` (the phase-encode winder) and, if it took its readout prephaser over, `x`.
`GRE3DTR` owns `y` and `z`. A slice-selective axis whose rephasing the excitation still realises
locally is **not** owned until the kernel folds it in, which is exactly what `GRE3DTR` already
does for a slab.

Unsupported combinations should refuse with the axis, the reason, and what would make it
supported — never silently do nothing.

**The spike does not do this yet, and the stress pass caught it.** Nothing validates the claimed
axis against the axes the repetition owns, so `GRE2DTR` accepts a claim on an axis that does not
exist, at 40 mT/m and 150 T/m/s:

```text
joint_claims                     min TE      winder
(none)                           4.662 ms     520 us
CommonModeClaim('y', 1)          5.732 ms    1590 us
CommonModeClaim('x', 1)          5.652 ms    1510 us
CommonModeClaim('z', 1)          6.542 ms    2400 us
CommonModeClaim('q', 1)          4.662 ms     520 us    <- accepted, identical to no claim
CommonModeClaim('zzz', 1)        4.662 ms     520 us    <- accepted, identical to no claim
```

Two failure shapes, both from the same missing check:

```text
target resolves to zero      reports success having emitted nothing -- a silent no-op that
  (common mode alone)        claims the repetition is flow-compensated on an axis it has not
                             got, which is the worst of the two

target is nonzero            refuses, but with the physics message: "allow a longer window",
  (a difference claim)       "relax the target", "widen the gradient limits".  The real cause
                             is that pypulseq cannot make a trapezoid on channel 'q', which
                             `_lobe` reads as "did not fit" and the cascade reads as infeasible
```

The fix is a capability check where the claims arrive, phrased physically — the axis must be one
the repetition owns an adjustable window on — not an `isinstance` test and not a hard-coded axis
letter list. It belongs with the public API, so it is recorded here rather than patched into the
spike.

---

## Scanner and profile configuration

Out of scope for this proposal, but the ownership it assumes:

```text
sequence source          expresses MRI intent
scanner / profile        hardware envelope, design policies, optional safety parameters
```

The same sequence source should be redesignable against a different declared scanner without
rewriting its MRI logic. No API is proposed here.
