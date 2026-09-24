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

**Leaning:** B for the repetition-level surface, with C's observation kept in mind — a velocity
acquisition is the thing that knows there are two states, so whatever generates them belongs
beside the line ordering rather than inside one repetition.

---

## `VelocityEncode`: eager realisation, and whether it matters

The constructor designs its standalone bipolar eagerly. Measured on a deliberately weak system
(8 mT/m, 25 T/m/s):

```text
venc 0.05 m/s   standalone builds, 8.00 ms pair
venc 0.02 m/s   standalone builds, 12.46 ms pair
```

It does not refuse — the duration search grows the window until it fits — so eager realisation is
a **cost**, not a correctness blocker: it designs a waveform the joint path may never emit.

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

---

## Scanner and profile configuration

Out of scope for this proposal, but the ownership it assumes:

```text
sequence source          expresses MRI intent
scanner / profile        hardware envelope, design policies, optional safety parameters
```

The same sequence source should be redesignable against a different declared scanner without
rewriting its MRI logic. No API is proposed here.
