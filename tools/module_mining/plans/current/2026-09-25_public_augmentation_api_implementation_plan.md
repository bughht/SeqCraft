# Public augmentation API — implementation plan

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-25, so that the pull request carries the direction behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../../README.md) for which one to edit.

**Status:** **implemented.** Approved in direction, revised after the correction pass, and built
on `feat/augmentation-api`. One adjustment was accepted during implementation and is recorded in
§1.1b; everything else shipped as planned.
**Date:** 2026-09-25, revised and implemented 2026-09-25
**Governed by:** [`2026-09-24_repetition_physical_design_architecture.md`](2026-09-24_repetition_physical_design_architecture.md)
**Evidence:** [`2026-09-24_repetition_joint_design_spike_findings.md`](2026-09-24_repetition_joint_design_spike_findings.md),
frozen at tag `stage-c-evidence` (`67c3e3c`)

Stage C is closed. This plan is written so that implementation can proceed mechanically without
reopening anything settled there. Where it names a decision the reviewer still owns, it says so.

---

## 0. Why this is the smallest sufficient API

Two public types, two keywords, one attribute, and nothing else.

```text
sc.FlowCompensation(axis=...)              what physics is wanted; one axis or several
sc.VelocityEncoding(venc_m_s=..., axis=...)   one axis -- two is a different acquisition

GRE2DTR(..., flow_comp=..., velocity_encode=...)     where it is wanted
tr(..., encoding_state=state)                        which state is being realised
venc.states                                          what states exist
```

Everything else the architecture needs already exists and stays internal: `CommonModeClaim`,
`DifferenceClaim`, `JointProblem`, `Schedule`, the realisation families, the schedule search, the
capability check.

It is sufficient because the four questions a caller actually has to answer are exactly the four
above: *which physics, on which axis, of which sequence, in which state.* It is smallest because
removing any one of them makes a real case unauthorable:

```text
drop the intent types        the keyword has to carry a tuple or a string, and `venc` needs
                             two numbers -- back to alternative A, which fails as soon as an
                             augmentation needs more than one value

drop the named keywords      a list arrives, and with it the question of what an unknown entry
                             means; that question is a plugin registry with the serial numbers
                             filed off

drop `encoding_state`        the repetition cannot realise a two-state acquisition at all

drop `venc.states`           the caller writes `for state in (+1, -1)` and the meaning of the
                             states leaves the type that owns the physics
```

And there is no room to make it smaller by *merging* the two intent types, because they differ in
kind: flow compensation constrains one repetition, velocity encoding defines a pair of them.

What it deliberately does **not** buy: a plugin point, a capability DSL, a scanner-profile surface,
an acquisition framework. Section 12 lists those as non-goals.

---

## 1. Public surface

### 1.1 The intent types

```python
@dataclass(frozen=True)
class FlowCompensation:
    """Ask that the COMMON-MODE first gradient moment be zero at the echo."""
    axis: str | tuple[str, ...]          # one axis, or several
    @property
    def axes(self) -> tuple[str, ...]: ...

@dataclass(frozen=True)
class VelocityEncoding:
    """Ask that the first moment differ by a known amount between two acquired states."""
    venc_m_s: float
    axis: str

    @property
    def delta_m1_s_per_m(self) -> float: ...
    @property
    def states(self) -> tuple[int, int]: ...      # (+1, -1)
    @property
    def axes(self) -> tuple[str, ...]: ...        # always one; see 1.1b
```

### 1.1b One axis or several — the accepted adjustment

**Accepted during implementation, 2026-09-25.** `FlowCompensation.axis` takes one axis or several;
`VelocityEncoding.axis` stays single.

The draft gave both a single `axis`, and implementation showed that one `flow_comp` keyword then
could not express flow compensation on two axes at once -- something the spike could do and a real
protocol wants, since compensating the phase-encode and slice axes together is ordinary. That is a
capability regression rather than a simplification.

```python
FlowCompensation(axis='y')                # one
FlowCompensation(axis=('y', 'z'))         # two independent requirements
```

**This does not reopen B2.** It is still one keyword carrying one intent object per augmentation,
with the set closed in the signature; only the intent's own field is richer. The spelling is the
one `spoil_axis` already uses on these kernels, so it is a convention the codebase has rather than
a new one.

The asymmetry with `VelocityEncoding` is principled, not an oversight:

```text
flow compensation, two axes    two independent requirements on two waveforms.  Asking for both
                               means the same as asking for each, and they do not interact

velocity encoding, two axes    a different acquisition.  Encoding two directions needs more than
                               two states -- four, or a three-directional scheme -- so it is not
                               this intent applied twice, and spelling it as though it were would
                               hide that the scan just got longer
```

Frozen dataclasses: they are values, not objects with behaviour. Two of them comparing equal means
the same physics was asked for, which is what a caller reading a protocol expects.

They carry physical intent and nothing else. Explicitly absent, and to stay absent:

```text
solver selection      realization family      timing policy       PNS policy
hardware profile      internal claims         Schedule            JointProblem
opts                  duration                any waveform at all
```

`VelocityEncoding` has no `opts`. This is the point of the type: **constructing it costs nothing
and designs nothing.** `delta_m1_s_per_m` delegates to `VelocityEncode.delta_m1_for(venc_m_s)`,
which stays the single source of the VENC relation (§4).

### 1.1a What FlowCompensation actually means

**Its reusable contract is common-mode first-moment nulling**, not `m1 = 0` for every state. The
distinction is invisible in the common case and load-bearing in the composed one:

```text
alone, single-state acquisition        the common mode IS the one state, so it reduces to
                                       m1 = 0 -- which is what a reader expects "flow
                                       compensated" to mean, and it is true here

with VelocityEncoding, same axis       it constrains the MEAN of the two states to zero while
                                       the difference claim constrains their separation, so
                                       the pair comes out at +delta/2 and -delta/2.
                                       Neither state has m1 = 0, and neither should
```

The symmetric plus-and-minus half-delta is therefore **derived, not specified**. Nothing in the
intent type, the translation layer or the kernel computes a half-delta; it falls out of resolving a
common-mode claim and a difference claim on the same axis and order. That is the property which
lets a third augmentation compose without a precedence rule, and the docstring must say *common
mode* rather than *zero* or the contract will be re-broken the first time someone reads it.

A consequence for the tests (§10): the whole-repetition assertion for flow compensation alone is
`m1 == 0`; combined with velocity encoding it is `mean(m1 over states) == 0`, and asserting
`m1 == 0` per state there would be asserting the wrong physics.

### 1.2 Where they live

New module `src/seqcraft/augmentation.py`, exported at **top level**:

```python
import seqcraft as sc
sc.FlowCompensation(axis='y')
sc.VelocityEncoding(venc_m_s=1.5, axis='z')
```

Not under `sc.modules`. `seqcraft.modules` is documented as reusable MR building blocks —
*parameters in, one `LogicBlock` out* — and every name in it is a `Module` subclass or a factory
returning one. An intent type takes no `opts` and emits nothing; putting it there would make that
contract false for two of its members, and the first question a reader would ask is what
`FlowCompensation(...).build()` does.

Top level already holds the vocabulary that is not itself a building block — `Raster`, `span`,
`barrier`, the error types — so a physical-intent value belongs there.

### 1.2a Prerequisite: move the shared validators to a neutral home

**A separate, mechanical commit, landing before any augmentation code.**

`require_axis` and `require_positive` live in `seqcraft/modules/_support.py`. `augmentation.py`
needs them and must not import them from there, for a reason stronger than taste:

```text
augmentation.py  imports  modules._support     (for the validators)
modules/kernel/* imports  augmentation         (for the keyword types and the routing)
```

That is an import cycle, and it is the kind that works until someone changes an import order. The
fix is to put the pure value validators below both:

```text
new         src/seqcraft/design/validation.py
moves       require_axis, require_positive, require_count, require_pair, require_range
re-exported from modules/_support.py so existing call sites keep working unchanged
```

Why `design/validation.py` and not `design/units.py`: `units.py` is the unit-conversion table and
`convert`, a different job. Why the whole `require_*` family and not only the two needed: they are
one coherent thing — check a value's domain, raise `ConfigurationError` — and splitting a family
across two modules to save two moves is how the next reader ends up looking in both places.

The move is safe: these five touch nothing but `ConfigurationError`, `format_error` and `AXES`,
all of which already live in `errors` and `design.events`, below `modules`. Call sites number 12
for `require_axis` and 17 for `require_positive`; keeping the re-export from `_support` means none
of them change in this commit, and they can be repointed later or never.

`design/validation.py` stays **internal**. It is not added to `docs/api_reference.md` §9 and not
exported at top level; `_support` re-exporting it preserves exactly today's public surface, which
is what makes this commit mechanical and separately reviewable.

### 1.3 The kernel keywords

```python
class GRE2DTR(Module):
    def __init__(self, *, opts, fov_mm, matrix, thickness_mm, ...,
                 flow_comp: FlowCompensation | None = None,
                 velocity_encode: VelocityEncoding | None = None,
                 tag=None): ...
```

Same two keywords on `GRE3DTR`. One keyword per augmentation, each accepting exactly its own type
or `None`. The closed set is closed **in the signature**: there is no place to put a third thing
without editing the signature, which is the property a list does not have.

Passing the wrong type is a `ConfigurationError`, not a `TypeError` with a traceback:

```text
flow_comp expects a FlowCompensation, and got VelocityEncoding.
  fix
    velocity encoding goes in velocity_encode=
```

### 1.4 How they are documented without becoming a plugin system

```text
docs/api_reference.md §9     both names added to the index, under `augmentation`
                             (a new entry in tools/check_api_reference.py `_INDEXED`)

docs/writing_a_module.md     unchanged.  These are not modules and the module-writing rules
                             do not apply to them

class docstrings             the physics each one asks for, the units, and what it does NOT
                             promise -- specifically that whether it can be honoured is the
                             repetition's answer, not the intent's
```

The thing that would turn this into a plugin system is a documented way for a caller to add a
third intent type. There is none, and the docstrings should not imply one: no "augmentations
implement the following protocol", no abstract base class, no registry. The two types share no
supertype beyond `object`. If a third is ever justified, it is a SeqCraft-owned type added to
`augmentation.py` with its own keyword, by the same process that added these two.

---

## 2. Capability belongs to the owner

Stage C implemented this and it does not change. Restated so the plan is self-contained:

```text
intent      says WHAT is wanted:  "first moment zero at the echo, on z"
owner       says WHETHER and HOW: "this repetition owns no adjustable z interval"
```

An intent validates only what is knowable without a sequence:

```text
FlowCompensation      axis is a well-formed axis name        via require_axis
VelocityEncoding      axis is well formed; venc_m_s > 0      via require_axis, require_positive
```

Those validators live in `seqcraft/modules/_support.py` today, and **that is not where
`augmentation.py` may reach for them.** See §1.2a: they move first, in their own commit.

It must not know, and has no way to learn, that `GRE2DTR` owns `y` jointly, that `GRE3DTR` also
owns `z`, that `CartesianLine` owns the readout problem locally, or that `Excitation` owns slice
rephasing. The same `FlowCompensation(axis='z')` value is honoured by `GRE3DTR` and refused by
`GRE2DTR`; a check living on the value could not produce that.

The implementation already present:

```text
kernel.joint_axes                  the axes the kernel MATERIALISES a joint design on
_joint.require_owned_axes(...)     refuses the rest, before any waveform is designed
```

§3 extends this by one step, from *which axes are jointly owned* to *which owner handles which
axis*.

---

## 3. Local-versus-joint routing

### 3.1 The rule, unchanged

```text
self-contained physical problem       -> the local owner solves it
cross-component or cross-state        -> repetition-level joint design
```

A readout-axis first moment is self-contained: the prephaser and the readout lobe are one
component's problem, and `CartesianLine` already solves it analytically (PR #39). A phase-encode
first moment is not: it shares a window with the encode area, the slice rephaser and, when
velocity encoding is present, a second state.

### 3.2 The mechanism — one private table per kernel

```python
# in GRE2DTR.__init__, before anything is designed
self._moment_owners: dict[str, str] = {'y': 'joint', 'x': 'readout'}

# in GRE3DTR.__init__
self._moment_owners: dict[str, str] = {'y': 'joint', 'z': 'joint', 'x': 'readout'}
```

A dict literal in the constructor. Not a registry, not a DSL, not public, not extensible at
runtime, and not something an intent can read. `joint_axes` becomes a derived view of it:

```python
@property
def _joint_axes(self) -> tuple[str, ...]:
    return tuple(a for a, owner in self._moment_owners.items() if owner == 'joint')
```

**Ownership metadata is internal implementation state, not API.** Stage C left `joint_axes` as a
public attribute; the API work renames it to `_joint_axes` and it stays private. A caller never
needs to ask a kernel which axes it owns — they state an intent and get either a sequence or a
refusal that explains itself (§8.1). Publishing the table would invite branching on it, and the
capability answer would then live in two places.

Two in-tree readers to update with the rename: the capability tests in
`tests/modules/test_joint_moment_design.py`, and `probe_capability` in
`tools/module_mining/stress_repetition_design.py`. Both are ours.

Routing is then a lookup, not a branch on the augmentation's class:

```python
def _route(self, intent) -> str:
    """Which owner handles this intent's axis, or a refusal naming the axis."""
    owner = self._moment_owners.get(intent.axis)
    if owner is None:
        _refuse_no_owner(type(self).__name__, intent.axis, self._moment_owners)
    return owner
```

`_route` reads `intent.axis` and nothing else about the intent. Adding a third augmentation
changes no routing code.

### 3.3 What each route does

```text
FlowCompensation(axis='x')  on GRE2DTR / GRE3DTR
    _route -> 'readout'
    the kernel already constructs its CartesianLine; it passes the readout's first-moment
    nulling through to it.  PR #39's analytic solve is reused unchanged.

FlowCompensation(axis='y')  on GRE2DTR / GRE3DTR
    _route -> 'joint'
    translated to CommonModeClaim('y', 1) and handed to the existing joint design

FlowCompensation(axis='z')  on GRE3DTR   -> 'joint'   (the z winder carries the partition
                                                       encode and the slab rephasing)
                            on GRE2DTR   -> no entry -> refusal (§8)

VelocityEncoding(axis=...)  -> 'joint'.  On an axis whose owner is 'readout' it is refused
                               for v1 (§3.3a), with a reason that says why (§8.2).
```

### 3.3a Why readout-axis velocity encoding is unsupported in v1

The reason is the **current local realisation's scope, not an architectural rule**:

```text
CartesianLine's local solve implements first-moment NULLING -- one fixed target, m1 = 0.
It does not implement arbitrary state-resolved M1 targets, which is what a difference
between two encoding states requires.
```

So `VelocityEncoding(axis='x')` is refused because *this* local owner cannot express that target,
and the refusal should say so rather than implying impossibility.

**Do not turn this into "cross-state intent can never be realised locally."** It is not true, and
writing it down would foreclose a cheap future: extending `CartesianLine` to accept a target `m1`
rather than only zero would make readout-axis velocity encoding locally realisable, and the
routing table above would need no change to pick it up. Whether that extension is worth doing is a
separate question with its own evidence, and nothing here should prejudge it.

The check itself is one condition in the translation layer — *this intent needs a state-resolved
target; does the routed owner offer one?* — not a per-augmentation branch in the router.

### 3.4 Combining routes

`flow_comp=FlowCompensation(axis='x')` and `velocity_encode=VelocityEncoding(axis='y')` in the
same repetition route to different owners and do not interact: different axes, one physical
problem per axis, which is the architecture's own rule. Both on the *same* axis is a conflict and
is refused (§8) — except for the one case the requirement layer composes by design, flow
compensation and velocity encoding on the same axis, which resolves to a common-mode claim and a
difference claim on the same axis and order and is exactly what `resolve_claims` exists for.

---

## 4. VelocityEncoding: state meaning versus state scheduling

### 4.1 Ownership, unchanged

```text
VelocityEncoding intent    the VENC relation, and the identity of the two states
caller / acquisition       whether both are acquired, and in what order
repetition                 realises the ONE state it is handed
joint designer             sees a resolved absolute target and no state semantics at all
```

### 4.2 The smallest user-facing pattern

No new acquisition framework. The seam already exists — `tr(..., encoding_state=state)` — and what
is missing is only that nothing tells the caller which states there are. `VelocityEncoding.states`
is that, and it is one property:

```python
venc = sc.VelocityEncoding(venc_m_s=1.5, axis='z')

tr = sc.modules.GRE3DTR(opts=opts, fov_mm=FOV_MM, matrix=MATRIX, velocity_encode=venc)

scan = sc.LogicBlock()
for partition in range(nz):
    for line in range(ny):
        for state in venc.states:
            scan.add(at, tr(line=line, partition=partition, encoding_state=state))
```

The ordering above — state innermost — is the caller's choice and the library does not impose it.
That is the point of leaving scheduling with the caller: state-innermost minimises the time
between the two acquisitions that will be subtracted, state-outermost minimises eddy-current
history differences, and which is right is a protocol decision.

`states` returns `(+1, -1)`. The representation is reviewable (an enum, or opaque sentinels, would
also work); the ownership is not.

### 4.3 What must not be required

Constructing `sc.modules.VelocityEncode` must not be a prerequisite for joint velocity encoding.
The stress evidence settled why this is only a cost and not a correctness issue — the standalone
pair never refuses — but it remains a waveform designed and thrown away. `VelocityEncoding` holds
no `opts` and therefore cannot build one even by accident.

`sc.modules.VelocityEncode` stays public and unchanged: it is the standalone bipolar, which is a
real thing to want on its own, and `delta_m1_for` on it stays the single source of the relation
that both paths read.

---

## 5. Internal claim translation

### 5.1 Where it happens

One new private module, `src/seqcraft/modules/_augment.py`:

```python
def claims_and_states(flow_comp, velocity_encode) -> tuple[tuple[object, ...], tuple[object, ...]]:
    """Public intent -> internal claims and the states they are stated over.

    The ONLY place an intent type is matched to a claim type.  Everything downstream --
    resolve_claims, JointProblem, the realisation families -- sees claims and numbers.
    """
    claims, states = [], ('only',)
    if velocity_encode is not None:
        states = velocity_encode.states
        claims.append(DifferenceClaim(velocity_encode.axis, 1,
                                      velocity_encode.delta_m1_s_per_m, states))
    if flow_comp is not None:
        # CommonModeClaim, not a per-state zero: this is the whole of FlowCompensation's
        # contract (1.1a).  With a difference claim beside it the two resolve to +-delta/2,
        # and nothing here computes that half.
        claims.append(CommonModeClaim(flow_comp.axis, 1))
    return tuple(claims), states
```

The translation is a direct restatement of §1.1a: *common mode* maps to `CommonModeClaim`, and the
`m1 = 0` a single-state caller sees is `resolve_claims` evaluating a mean over one state. There is
no branch on how many states there are, and there must not be one — that branch is exactly how a
per-state zero would creep back in and break the composed case.

### 5.2 The layering it preserves

```text
public intent        FlowCompensation / VelocityEncoding
       |   _augment.claims_and_states   <- the one dispatch point, outside the numerics
       v
internal claims      CommonModeClaim / DifferenceClaim
       |   _joint.resolve_claims
       v
absolute targets     one (m0, m1) per state
       |
       v
realisation          sees numbers.  Never learns an augmentation existed
```

Two properties to hold, both testable:

- **The kernel does not re-derive VENC physics.** `delta_m1_s_per_m` is read from the intent,
  which reads `VelocityEncode.delta_m1_for`. Grep-able: `1.0 / (2.0 *` appears once in the source.
- **No `isinstance` on an augmentation below the translation layer.** `_joint.py` and the kernels'
  design methods never import `augmentation`. The existing test
  `test_the_kernel_adapter_does_not_branch_on_the_augmentation` already asserts the shape of this
  and should be extended to the new entry point.

---

## 6. Migration off the spike parameters

### 6.1 What disappears

```text
GRE2DTR(joint_claims=..., encoding_states=...)      both removed
GRE3DTR(joint_claims=..., encoding_states=...)      both removed
```

They were never released — they exist only on the spike branch — so there is no deprecation cycle
and no compatibility shim. They are deleted, not aliased.

### 6.2 What replaces them

```text
joint_claims        -> flow_comp= / velocity_encode=, translated by _augment.claims_and_states
encoding_states     -> derived from velocity_encode.states, or ('only',) when there is none
encoding_state=     -> STAYS on build().  This is the seam, not a spike artefact
_moment_owners      -> new, private, per kernel (§3.2)
joint_axes          -> becomes _joint_axes, private; derived from _moment_owners
```

### 6.3 Tests

17 call sites across `tests/modules/test_joint_moment_design.py`, in three groups:

```text
requirement-layer tests        keep using CommonModeClaim / DifferenceClaim directly.  They
(resolve_claims, duplicate     test the internal layer and should go on testing it in its own
 claims, missing states)       vocabulary.  No change.

realisation-family tests       never touched the kernels.  No change.
(realise_*, utilisation)

whole-repetition and timing    rewrite to the public keywords.  These are the tests that
tests, capability tests        should read the way a user writes, and after this they do.
```

**No private injection path is proposed.** The requirement-layer tests already construct claims
directly and call `_joint` functions directly; nothing needs to push a raw claim *through a
kernel*. If implementation finds a case that genuinely does, add `_claims=` as a clearly private
keyword with a docstring saying it is a test seam — but the expectation is that it is not needed,
and adding it speculatively would leave the spike API public under a new spelling.

---

## 7. PR #39: the public-surface decision

### 7.1 What is uncontested

Keep, unchanged: the readout-axis M0/M1 solve, the analytic realisation, the
`prephaser_lobes` accessor, `_compensated_areas`, the validation and the 45 tests. This is the
substrate `FlowCompensation(axis='x')` routes to, and it is good work independent of how it is
spelled.

### 7.2 The decision — **made: internal**

The reviewer has decided: `null_moment_order` becomes internal before PR #39 merges.

The asymmetry that decided it:

```text
ship it public, retire it later     a deprecation cycle on a released name, and a period
                                    where the docs must explain two ways to ask for the
                                    same physics and when to use which

keep it private, expose it later    a one-line change, non-breaking, whenever a real
                                    standalone need is demonstrated
```

PR #39 is unmerged, so choosing privacy now costs nothing and choosing it later costs a
deprecation. The architecture §22 N×M argument also applies: a public `null_moment_order` on
`CartesianLine` would be the first instance of *per-leaf, per-augmentation* keywords, with no
principled place to stop.

### 7.3 What was weighed against it

*(Paths in this section are as they were written. That notebook is now
`examples/gre_2d/03_flow_comp.ipynb`.)*

Recorded because it was a real argument, not to reopen the decision. `CartesianLine` is a public
leaf usable without a kernel, and the `flowcomp_gre_2d` notebook built one directly; a caller
assembling a custom readout has no repetition to ask. The decision accepts that cost on the
grounds that re-exposing a private parameter later is cheap, and that whether an external caller
wants a standalone first-moment-nulled readout is still unknown.

### 7.4 Concretely, either way

The finalization commit on `feat/flow-compensation`:

```text
src/seqcraft/modules/readout/cartesian_line.py
    null_moment_order          -> _null_moment_order, constructor parameter made private
    the docstring section      -> marked internal; the physics prose is kept, the "you can
                                  pass this" framing removed
    prephaser_lobes            -> KEEP.  It is a readable accessor for a two-lobe prephaser
                                  and does not advertise the option

docs/api_reference.md          -> drop the public mention of null_moment_order

tests/modules/test_flow_compensation.py, test_cartesian_line.py
                               -> KEEP all 45.  They test the solve, not the spelling;
                                  update the parameter name only

examples/gre_2d/03_flow_comp.ipynb      -> MOVE to the integration PR, where it can be rewritten
                                  against flow_comp=sc.FlowCompensation(axis='x').  It
                                  cannot be written against the public API before that API
                                  exists, and it must not ship demonstrating a private
                                  parameter

CHANGELOG.md                   -> reword to describe the capability rather than the keyword
```

The notebook is the one thing that moves rather than lands. That is the right trade: the analytic
solve and its 45 tests are the valuable part of PR #39 and should merge early, and a demonstration
notebook is exactly the kind of thing that belongs with the surface it demonstrates.

---

## 8. Refusal and debug behaviour

Every message explains the physical situation. None of them names a claim type, a schedule, a
family or a solver. Proposed wording:

### 8.1 No physical owner for the axis

```text
GRE2DTR cannot apply flow compensation on 'z'.
  axis            :  z
  axes it can use :  ('y', 'x')
  fix
    this repetition does not own an adjustable z interval between the excitation and the
    echo -- its z gradient is the slice rephaser, which Excitation realises for itself
    GRE3DTR does own z, because its z winder carries the partition encode
```

### 8.2 Augmentation unsupported on this route

```text
GRE2DTR cannot velocity encode on 'x'.
  axis  :  x
  fix
    the readout axis is designed by CartesianLine, which solves one repetition at a time
    velocity encoding is a difference between two acquisitions, so it needs an axis this
    repetition designs jointly: ('y',)
```

### 8.3 Missing encoding state

```text
this repetition velocity encodes, so build() needs to be told which state it is realising.
  states  :  (1, -1)
  fix
    pass encoding_state=+1 or encoding_state=-1
    acquire both, and subtract their phase to get velocity
```

### 8.4 Unknown encoding state

```text
encoding_state=0 is not one of this repetition's states.
  given   :  0
  states  :  (1, -1)
  fix
    iterate the states the intent defines:  for state in venc.states:
```

### 8.4a A state supplied where none is meant

The symmetric case, and the one easiest to leave out:

```text
this repetition does not velocity encode, so encoding_state=+1 has no meaning here.
  given  :  1
  fix
    drop encoding_state -- this repetition has one state and build() already realises it
    or pass velocity_encode=... if the acquisition is meant to be phase contrast
```

Both directions refuse, so **state semantics are never silently ignored**. Accepting an
`encoding_state` that nothing acts on is the worse of the two failures: the caller writes a
two-state loop, gets two identical repetitions, and finds out when the subtraction is zero.

The check is a pair of conditions in `build()`, not a per-augmentation branch:

```text
states != ('only',)  and  encoding_state is None       -> 8.3
states == ('only',)  and  encoding_state is not None    -> 8.4a
```

### 8.5 Duplicate or conflicting intent

```text
flow compensation and velocity encoding both ask for the first moment on 'y', and they
disagree about what it should be.
```

— only when they genuinely conflict. The composable case (common mode plus difference on the same
axis and order) is **not** an error and must not be refused; it is the case the requirement layer
was built for.

### 8.6 Explicit TE made infeasible by the augmentation

Existing behaviour, unchanged, and already correct: the refusal is the kernel's and names `te_s`,
not the designer's naming an axis.

```text
te_s = 4.000 ms is shorter than this repetition can achieve.
  te_s     :  0.004
  min_te_s :  0.00573
  fix
    pass te_s >= 0.00573
    flow compensation on y lengthened the winder, and with it the shortest echo time
```

The last fix line is new and worth adding: it connects the number to the cause.

### 8.7 Design search limit

Implemented in Stage C, unchanged. Says "design search limit reached ... windows longer than
40.0 ms were never tried", never "infeasible".

---

## 9. Examples, as acceptance criteria

The test of this API is whether the authoring code reads like a protocol. Each of these is
proposed verbatim as what the user writes.

### 9.1 GRE2D + flow compensation on y

```python
tr = sc.modules.GRE2DTR(
    opts=opts, fov_mm=220.0, matrix=(128, 128), thickness_mm=5.0, flip_deg=15.0,
    flow_comp=sc.FlowCompensation(axis='y'),
)
scan = sc.LogicBlock()
for line in range(128):
    scan.add(line * tr.tr_s, tr(line=line))
```

### 9.2 GRE3D + flow compensation on z

```python
tr = sc.modules.GRE3DTR(
    opts=opts, fov_mm=(220.0, 220.0, 120.0), matrix=(128, 128, 16),
    flow_comp=sc.FlowCompensation(axis='z'),
)
```

### 9.3 GRE3D + velocity encoding

```python
venc = sc.VelocityEncoding(venc_m_s=1.5, axis='z')
tr = sc.modules.GRE3DTR(..., velocity_encode=venc)

for partition in range(nz):
    for line in range(ny):
        for state in venc.states:
            scan.add(at, tr(line=line, partition=partition, encoding_state=state))
```

### 9.4 GRE3D + both

```python
venc = sc.VelocityEncoding(venc_m_s=1.5, axis='y')
tr = sc.modules.GRE3DTR(
    ...,
    flow_comp=sc.FlowCompensation(axis='y'),      # the common mode: no background flow phase
    velocity_encode=venc,                          # the difference: the velocity signal
)
```

Same axis, deliberately, and the two compose rather than conflict. Flow compensation fixes the
**mean** of the two states at zero; velocity encoding fixes their **difference**; the symmetric
plus-and-minus half-delta falls out. Neither state has `m1 = 0` here, and neither should — that is
what the velocity signal *is*. A caller does not need to know any of that sentence: they need the
code above to work, and they need "flow compensation" to keep meaning the same thing whether or not
velocity encoding is present beside it.

### 9.5 Flow compensation on x, routed to the local realisation

```python
tr = sc.modules.GRE2DTR(..., flow_comp=sc.FlowCompensation(axis='x'))
```

Identical spelling to 9.1. That it reaches `CartesianLine`'s analytic solve rather than the
repetition designer is invisible, which is the routing requirement.

### 9.6 A refused combination

```python
tr = sc.modules.GRE2DTR(..., flow_comp=sc.FlowCompensation(axis='z'))
# ConfigurationError, §8.1
```

### 9.7 The acceptance criterion

A sequence author writing any of the above never needs to know: M0/M1 claim decomposition, common
versus difference components, the schedule search, realisation families, or whether their axis was
served locally or jointly. If review finds a case where they do, the API is wrong, not the caller.

---

## 10. Validation plan

Stage C evidence is reused, not re-derived. Implementation tests:

```text
translation        FlowCompensation      -> CommonModeClaim(axis, 1)
                   VelocityEncoding      -> DifferenceClaim(axis, 1, delta, states)
                   delta matches VelocityEncode.delta_m1_for exactly
                   no isinstance on an augmentation below _augment

routing            x  -> CartesianLine local solve, emitted readout m1 nulled at the echo
                   x  -> VelocityEncoding refused, with the 3.3a reason and not an
                         impossibility claim
                   y  -> joint design on both kernels
                   z  -> joint on GRE3DTR, refused on GRE2DTR
                   the route is invisible in the emitted result: same physics either way

physics            whole-repetition M0 and M1 residuals at the achieved echo, per state,
                   off-centre lines and partitions -- the Stage C thresholds
                   (m0 < 1e-6 relative, m1 < 1e-11) carried over unchanged

                   flow compensation ALONE      m1 == 0             (one state, so the
                                                                     common mode is it)
                   velocity encoding ALONE      m1(+) - m1(-) == delta
                                                and the common mode is UNCONSTRAINED --
                                                assert the difference, not the values
                   both, same axis              mean(m1 over states) == 0
                                                AND m1(+) - m1(-) == delta
                                                the +-delta/2 symmetry is then DERIVED, so
                                                assert it as a consequence rather than
                                                testing it as the specification

timing             AUTO TE/TR moves when an augmentation is enabled and not otherwise
                   explicit TE below the new minimum refuses, naming te_s
                   explicit TE above it does not stale the first moment (the Stage C bug)

refusals           each of §8.1-8.5 by message content, not by exception type alone
                   both state refusals: state supplied with no velocity encoding, and
                   velocity encoding with no state (§8.3, §8.4a)

no-op              with neither keyword, every emitted repetition is byte-identical to
                   today's.  This is the regression that matters most and is cheap:
                   compare the flattened event streams

integration        sc.compile accepts a jointly designed repetition on both kernels
                   docs/api_reference.md examples execute; §9 index matches
                   notebook smoke passes
```

**Not** an implementation gate: the adversarial stress matrix. It is the architecture regression
tool and is re-run only if the implementation changes the solver itself — which this plan does not.
`tools/module_mining/stress_repetition_design.py --quick` is a reasonable sanity check before the
integration PR, not a required one.

---

## 11. Branch and PR execution

### 11.1 Where things stand

```text
origin/main                        982b308   (PR #38 merged)
origin/feat/flow-compensation      9967ffc   PR #39, 4 commits, OPEN, untouched
spike/repetition-joint-design      67c3e3c   8 further commits, FROZEN
                                             tagged stage-c-evidence
feat/augmentation-api              67c3e3c   new, branched from the frozen spike
```

The spike was branched from **PR #39's tip**, not from main, so it already contains PR #39's four
commits. That is the fact the rebase has to respect.

### 11.2 The one thing that makes the rebase clean

The eight spike commits touch **no file that PR #39 touches**:

```text
PR #39               cartesian_line.py, test_cartesian_line.py, test_flow_compensation.py,
                     flowcomp_gre_2d/01_build.ipynb, megre_2d/01_build.ipynb,
                     CHANGELOG.md, docs/api_reference.md, examples/README.md,
                     flow_compensation.yaml, run_notebook_smoke.py

spike after 9967ffc  _joint.py, velocity_encode.py, gre_2d_tr.py, gre_3d_tr.py,
                     test_joint_moment_design.py, stress_repetition_design.py,
                     the plans, sources.yaml, the module_mining READMEs
```

Disjoint. The integration work will add `CHANGELOG.md` and `docs/api_reference.md` entries and so
will overlap there; those are the only conflicts to expect, and they are additive.

### 11.3 The sequence

The PR #39 surface decision is made, so the integration branch rebases onto the **updated** PR #39
tip *before* implementation, not onto main afterwards. Implementing against the final local
substrate is the point: `FlowCompensation(axis='x')` routes into `CartesianLine`, and it should be
written against the private parameter it will actually call.

```sh
# ---- done already -----------------------------------------------------------------
#   tag stage-c-evidence at 67c3e3c            spike frozen
#   branch feat/augmentation-api from it       this plan committed there
#   PR #39 untouched at 9967ffc

# ---- step 1: finalize PR #39 (its own work, on its own branch) ---------------------
git switch feat/flow-compensation
#   null_moment_order -> _null_moment_order, docs and tests per 7.4,
#   notebook moved out to the integration PR
git commit ...                                 # -> NEW PR #39 TIP, call it <pr39-final>

# ---- step 2: rebase the integration branch onto that tip ---------------------------
git rebase --onto feat/flow-compensation 9967ffc feat/augmentation-api
#   replays the 8 spike commits + this plan onto <pr39-final>.
#   9967ffc is the ORIGINAL fork point and stays the rebase base in every step below:
#   it is what "everything after this" means, regardless of what PR #39's tip becomes.

# ---- step 3: implement the API on feat/augmentation-api ---------------------------
#   1.2a validator move first, then augmentation.py, routing, translation, migration

# ---- step 4: merge PR #39 to main -------------------------------------------------
#   the repo merges with merge commits (PRs #36-#38)

# ---- step 5: rebase onto main, dropping everything through the updated tip ---------
git fetch origin
git rebase --onto origin/main <pr39-final> feat/augmentation-api

# ---- step 6: verify nothing was reintroduced --------------------------------------
git log --oneline origin/main..feat/augmentation-api    # spike + plan + API commits only
git diff --stat origin/main...feat/augmentation-api     # no cartesian_line.py unless the
                                                        # API work deliberately touched it

# ---- step 7: open the integration PR against main ---------------------------------
```

Two different bases, and using the wrong one is the failure this section exists to prevent:

```text
step 2   --onto feat/flow-compensation  9967ffc        the ORIGINAL fork point.  Everything
                                                       after it is spike + plan, none of which
                                                       PR #39 has.

step 5   --onto origin/main            <pr39-final>    the UPDATED tip.  Everything through it
                                                       is now in main via the merge, so
                                                       replaying it would duplicate PR #39.
```

`--onto` with an explicit base in both, never a bare `git rebase`. A bare one relies on git
skipping already-applied commits by patch-id, which a merge commit makes *usually* true and which
fails outright if PR #39 is squashed — and PR #39 is being modified in step 1, so its patches will
not match anyway.

### 11.4 Current tips, for the record

```text
origin/main                     982b308     unchanged
feat/flow-compensation          9967ffc     PR #39, OPEN, UNCHANGED -- step 1 not yet done
spike/repetition-joint-design   67c3e3c     frozen, tagged stage-c-evidence
feat/augmentation-api           <this plan> branched from 67c3e3c

integration base after step 2   <pr39-final>, which does not exist yet
rebase base for steps 2         9967ffc      (fixed, the original fork point)
rebase base for step 5          <pr39-final>
```

**PR #39 has not been modified.** This pass revised the plan only, as instructed; step 1 is the
next action and needs one word to start.

## 12. Non-goals for this phase

Not in this work, and not to be added opportunistically:

```text
PNS implementation        SAR model              RF-energy reporting
GrOpt                     generic optimizer      scanner-profile public API
acquisition framework     plugin registry        pathway-aware moments
M2 and above              bSSFPTR                public capability DSL
public requirement objects
```

---

## 13. Decisions, and what is left open

Recorded as settled. Implementation follows these without reopening them.

```text
sc.FlowCompensation           APPROVED    the intent type and the name
sc.VelocityEncoding           APPROVED    the intent type and the name; the near-collision
                                          with sc.modules.VelocityEncode is accepted, and
                                          the two are distinguished in the docs by what they
                                          are -- an intent, and the standalone bipolar Module
states = (+1, -1)             APPROVED    plain ints, consistent with
                                          VelocityEncode.build(polarity=...)
PR #39 null_moment_order      INTERNAL    before merge, per 7.2

sc.modules.VelocityEncode     UNCHANGED   stays public, stays the standalone bipolar Module,
                                          stays the single source of delta_m1_for
joint_axes                    INTERNAL    ownership metadata is implementation state
```

Nothing is open that blocks implementation. Two things are deliberately deferred and named here so
they are not mistaken for oversights:

```text
readout-axis velocity encoding      unsupported in v1 because CartesianLine's local solve
                                    implements nulling, not state-resolved targets (3.3a).
                                    An implementation limit with a known cheap fix, not an
                                    architectural rule, and not scheduled.

a third augmentation                no mechanism exists and none is planned.  Adding one is
                                    a SeqCraft-owned type, a keyword and a routing entry, by
                                    the same process that added these two.
```
