# C3 gap note — what v0 could not represent, and what is actually missing

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-22, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way.

Companion to [`velocity_encoding.yaml`](velocity_encoding.yaml) and
[`flow_compensation.yaml`](flow_compensation.yaml). Those hold what v0 can honestly express; this
holds what it could not, and why. **Phase E evidence. No Skill or schema change was made.**

---

## 1. The schema holds one status, and the candidate was two

C3 was one scan of one question — *are velocity encoding and flow compensation one reusable
contract with different target values?* The answer is **one requirement model, two owners**:

```text
velocity encoding    NEW_LEAF        an encoding leaf; needs nothing from outside
flow compensation    EXTEND_EXISTING not a module; it reshapes a waveform a leaf already owns
```

A record has one `status` and one `traffic_light`, so either would have been a distortion. Rather
than pick one, the scan produced **two records with a shared findings document**. That is not a
workaround — it is the honest shape, and it worked. What v0 lacks is only a way to say *"these two
records came from one scan and one question"*; the folder and this note say it in prose.

**Not a proposed schema change.** A `scan_id` linking records would be a field invented from one
case, which is the rule the Skill sets for itself.

## 2. The capability that is genuinely missing

Phase B saw a real 4D-flow implementation state a residual moment requirement and hand it to an
external constrained optimiser, and inferred that SeqCraft needs a way to express requirements
before realisation. **The domain source corrects that.**

The Handbook reduces the hard canonical case — velocity-compensating a slice-selection waveform,
which must account for the first moment the excitation already accumulated from the RF isodelay
point — to *"a quadratic equation for the unknown lobe area"* with an explicit physically
significant root (§10.4, p. 343). Target, neighbour's contribution, hardware limits, in; lobe
areas and widths out. **Algebra, not search.**

So the correction to Phase B is precise, and it is easy to over-state:

```text
WRONG      Open4DFlow uses GrOpt, therefore SeqCraft needs a general solver layer
CORRECT    the canonical cases do not justify introducing a general solver layer
ALSO WRONG coupled gradient design never needs an optimizer
```

The Handbook gives closed forms for canonical constructions and **does not establish them for the
general problem.** Several moment orders at once, several semantic timing windows, a minimum-TE
search, fixed waveform segments, non-zero start or end gradients, arbitrary pre-existing
waveforms, eddy-current, concomitant-field or PNS constraints — for those a numerical optimizer
remains a legitimate *realization backend*. What the evidence rules out is founding SeqCraft's
public semantics on one.

Three questions that must be kept apart, because answering one does not answer the others:

```text
WHAT must be true?   the physical requirement -- a target moment at a semantic instant
WHO owns it?         which abstraction holds the coupled design
HOW is it realised?  a closed-form construction, or a numerical backend
```

The same owner can serve a closed-form case and an optimizer-backed one.

What *is* missing is narrower and sharper:

> **The most efficient physical waveform may span boundaries that SeqCraft currently assigns to
> separate leaf/module instances.**

Stated that way rather than as "one abstraction cannot own a waveform several leaves emit",
because in the preferred resolution those leaves do not emit the contested pieces at all — a joint
owner holds their coupled realization **before any events are materialised.**

Two independent instances, from two different fine scans, neither looking for it:

| | the efficient form | what it needs | measured cost of not having it |
|---|---|---|---|
| **C2 bSSFP** | the rewinder of repetition *n* and the prephaser of *n+1* are one continuous lobe | a waveform continuing across a repetition boundary | **220 µs per axis per repetition**, ~5.5 % of a 4 ms TR |
| **C3 flow** | the bipolar encoding lobe merged with a compensated imaging lobe — *"five lobes can be combined into three"*, which *"always reduces the minimum TE"* (§9.2.3, p. 290) | one abstraction reshaping lobes that `Excitation` or `CartesianLine` emit | two lobes of ramp time; at 40 mT/m and 180 T/m/s the ramp scale is 222 µs, so of order 0.4–0.9 ms depending on the design |

Same cause, different sequences: **the physical object spans what SeqCraft has made into separate
modules.** That is the ownership question Phase E should look at, and it is independent of how any
particular waveform is computed.

**And it is not the only way to read C2's case.** Where two already-decided events could be fused
with *no change to physical semantics* — same integrated waveform, same RF and ADC timing, same
requirements, merely a more efficient continuous representation — that is a candidate for a
semantics-preserving optimisation in the compiler, not a physical redesign. C3's merged waveform
is emphatically not of that kind: its final shape depends on VENC and on a moment target, so it
must be designed with the physics in view. Phase E should keep the two apart:

```text
physics-aware joint realization      -> Module / kernel layer
semantics-preserving waveform fusion -> possible future compiler optimisation
```

**bSSFP physics must not end up in the compiler** as a way of recovering the 220 µs.

## 3. Where v0's ownership model stopped being expressive

The ownership test asks whether a component can determine a value from its own contract and
parameters, and the layer table says the kernel owns "coupling between leaves that no leaf can
see". Both held up here — flow compensation's cross-leaf case *is* kernel work by that definition,
and the test said so cleanly.

What the model has no vocabulary for is the distinction between two kinds of kernel:

```text
a PLACEMENT kernel                 GRE2DTR, TSEShot -- it decides offsets and hands each
                                   leaf's events through unchanged

a JOINT-REALIZATION kernel         what a merged flow waveform or a cross-repetition
                                   balanced lobe would need -- it takes physical facts and
                                   requirements from several leaves and jointly designs the
                                   replacement waveform
```

Every shipped kernel is the first kind. Nothing in v0 says the second kind is allowed, forbidden,
or a different layer. That is the gap, stated once, from two cases.

**The second kind does its work before materialisation.** It does not build child `LogicBlock`s,
take them apart and rewrite their events; it owns the coupled physical design and then emits the
final events once. That keeps the existing architecture intact — `Module → build → LogicBlock →
tree → compiler` — and puts the new capability *inside* the Module layer rather than adding a
layer to it.

## 4. Where the validation model stopped being expressive

One concrete thing. The acceptance claim that matters for flow compensation is:

> the n-th moment of **everything emitted on the axis**, from a semantic origin to the echo,
> equals the target

That is a claim about a **sum across events from several modules**, not about any module's own
output. Rules E and F are both per-module — measure on the lattice *you* emit, validate occupancy
after *your* events are placed. Neither reaches a summed cross-module quantity, and `sc.moments`
takes a tree, so the measurement is available while the rule that would require it is not.

Not a proposed rule. One case is one case, and C2's equivalent claim is about a waveform's
continuity rather than a sum, so the two are not yet the same shape.

One thing the shape does settle: **the validator must integrate the emitted gradients itself.**
Asking a module for its own `M1` and checking that against its own target tests arithmetic, not
physics — and a moment is not an intrinsic property of an event anyway, since it depends on the
placement time, the semantic origin, the endpoint and the order. The invariant is a property of
the *composition*, so the measurement must be too.

## 5. What the records deliberately do not contain

No `RequirementIR`, no `AugmentationPlugin`, no `ConstraintSolver`, no new schema field, no new
vocabulary, and no change to `.claude/skills/module-mining/`. The hand-added `evidence_class:
domain-reference` appears for the third time and is still hand-added.
