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

**The reason is representation, not validation.** A small intent object gives one documented
spelling of physical intent and keeps the claim algebra — `DifferenceClaim`, `CommonModeClaim`,
`JointProblem`, `Schedule` — off the public surface entirely. That is the whole of the argument.

It is explicitly **not** that the intent object is where capability is validated. Capability
belongs to the repetition adapter and is now implemented there (see **Capability and refusal**
below). An intent object may check that an axis name is syntactically well formed; it must not
know which axes a particular repetition owns, because it cannot — `z` is owned by `GRE3DTR` and
not by `GRE2DTR`, and the object is the same object in both cases.

### The two remaining B spellings

```python
# B1 -- a list
GRE2D(..., augmentations=[FlowComp(axis='y'), VelocityEncoding(venc_m_s=1.5, axis='y')])

# B2 -- named keywords carrying intent objects
GRE2D(..., flow_comp=FlowComp(axis='y'), velocity_encode=VelocityEncoding(venc_m_s=1.5))
```

```text
                          B1 list                        B2 named keywords
a third augmentation      add a type                     add a type and a keyword
reads like                a set of things applied        a protocol card
two of the same kind      possible, and meaningless      impossible by construction
                          -- has to be refused           -- the signature refuses it
discoverability           one name to find               each shows in the signature
                          the set behind it              and in help()
risk                      invites "register your own"    none of that shape
```

**Leaning B2**, narrowly. The closed set is the point, and a keyword per augmentation makes the
set closed *in the signature* rather than by documentation and a runtime check. B1's advantage —
a third augmentation costs one type instead of a type and a keyword — is small when the set is
expected to stay at three or four.

**If B1 is chosen instead**, then `augmentations=[...]` must be defined as a **closed,
SeqCraft-owned set of intent types**: a fixed union the kernel knows by name, not a registration
point, not a plugin API, and not an extension hook. An unknown type is a refusal, not a
pass-through. The architecture rules out a plugin bus and the spelling must not smuggle one back.

### Where a state-generating intent is consumed

Velocity encoding differs from flow compensation in a way the spelling has to respect: it
**creates states**, and a repetition realises one state but owns neither the ordering nor the
existence of the pair.

```text
acquisition layer     knows there are two states, orders the shots, doubles the scan,
                      and defines the subtraction reconstruction will do
        |  passes one state down into build()
        v
repetition            realises THAT state: the claim resolves to one absolute target per
                      state, and the joint design serves the state it was handed
        |
        v
realisation           sees a number.  It never learns a state existed
```

So `VelocityEncoding(venc_m_s=...)` on the repetition declares *what a state means* — the venc
relation, through `VelocityEncode.delta_m1_for` — while the acquisition decides *that there are
two of them* and in what order. Today `build(..., encoding_state=...)` is the seam and the
acquisition does not yet pass anything down; that gap is the concrete work the public API needs,
and it is larger than the keyword question.

---

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

**Implemented at the repetition adapter, not deferred to the public API.** An augmentation must
not be assumed to work on every axis of every sequence, and the question is physical:

```text
does this repetition own an adjustable window on that axis, between the semantic
origin and the endpoint, that it will actually MATERIALISE?
    yes  -> it can carry a moment claim there
    no   -> refuse, before any waveform is designed, naming the axis
```

The last clause is the one the stress pass had to teach. "Owns gradients on that axis" is not the
same question as "materialises a joint design on that axis", and the difference was expensive:

```text
GRE2DTR, at 40 mT/m and 150 T/m/s, before the check existed

joint_claims                winder     min TE      emitted m1 on the claimed axis
(none)                       520 us   4.662 ms     --
CommonModeClaim('y', 1)     1590 us   5.732 ms     -0.000000   designed and emitted
CommonModeClaim('x', 1)     1510 us   5.652 ms      0.597804   designed, costed, DROPPED
CommonModeClaim('z', 1)     2400 us   6.542 ms     -0.736067   designed, costed, DROPPED
CommonModeClaim('q', 1)      520 us   4.662 ms      --         silent no-op
```

`x` and `z` are legal axes `GRE2DTR` plays gradients on. It does not own their winders: the
readout prephaser belongs to `CartesianLine` and the slice rephaser to `Excitation`. So a claim on
either was designed, grew the winder, pushed TE out by up to **1.9 ms that the caller paid for**,
and was then dropped unemitted — leaving the axis reported as compensated and carrying
`m1 = -0.74 s/m`. A nonexistent axis was worse in a different way: accepted, identical to no claim,
a silent no-op claiming compensation on an axis that is not there.

### The owned axes, verified by what is emitted

```text
GRE2DTR     ('y',)         the phase-encode winder, materialised through `_encode`
GRE3DTR     ('y', 'z')     the phase-encode blip and the z winder carrying the partition
                           encode; `_z_winder` folds the slab rephasing in, which is what
                           makes z genuinely owned here and not in the 2D kernel
```

Each kernel declares this as `joint_axes` and calls `_joint.require_owned_axes` **before**
designing anything. The check is on the axis the repetition emits, not on the claim's type and not
on whether the letter is one of `x`, `y`, `z`; `require_owned_axes` never inspects the claim
beyond reading its `axis`, so the same claim object is accepted by `GRE3DTR` and refused by
`GRE2DTR`, which is the behaviour a class-name test could not produce.

Three tests hold it: every advertised axis is designed, **emitted**, and met on the complete
repetition; every unowned axis refuses before realisation; and the same claim is accepted or
refused depending only on the repetition it is handed to.

### Search exhaustion is not infeasibility

A related honesty problem in the same refusal path. The schedule search stops after
`_joint.SEARCH_LIMIT_WINDOWS` candidate windows, and reaching that ceiling used to raise "no
candidate schedule realises these moments" — a claim the search never established, since the next
window along was simply never tried. The two now read differently, because they need different
fixes:

```text
windows ran out       "design search limit reached ... not a proof of infeasibility:
                       windows longer than 40.0 ms were never tried"
a schedule failed     "no candidate schedule realises these moments on y"
```

The ceiling stays a fixed guard rather than becoming a growth policy. One consequence is worth
knowing: because the split search tries every raster-aligned split at every candidate window, it
is quadratic in window length, and walking all 4000 windows to reach the ceiling takes the better
part of a minute. An absurd target therefore refuses slowly. That is a cost of the exact search,
and making it dynamic would trade it for the risk of stepping over the true minimum — which this
same stress pass showed is not always where a monotone search would expect.

---

## Scanner and profile configuration

Out of scope for this proposal, but the ownership it assumes:

```text
sequence source          expresses MRI intent
scanner / profile        hardware envelope, design policies, optional safety parameters
```

The same sequence source should be redesignable against a different declared scanner without
rewriting its MRI logic. No API is proposed here.
