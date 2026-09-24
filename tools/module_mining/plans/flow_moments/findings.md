# Velocity encoding + flow compensation — Phase C third fine scan

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-22, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way.

**Date:** 2026-09-22 · **Skill:** v0, frozen and unmodified
**`VelocityEncode` implemented** 2026-09-24; `flow_compensation.yaml` is still unbuilt and is the
next stage. The scan record below is left as it was written and
[§8](#8-what-implementing-velocityencode-added-2026-09-24) records what it added.
**Result:** one requirement model, **two owners** — so two records:

```text
velocity_encoding.yaml     NEW_LEAF        / APPROVED_FOR_IMPLEMENTATION
flow_compensation.yaml     EXTEND_EXISTING / APPROVED_FOR_IMPLEMENTATION
```

Nothing implemented. Domain card:
[`../../domain_evidence/gradient_moments.md`](../../domain_evidence/gradient_moments.md).
What v0 could not express: [`gap_note.md`](gap_note.md).

## 1. The question, and the answer

*Are these one reusable contract with different target values, or two visually similar gradient
families?* **Neither.** They share one requirement form:

```text
on one axis, over an interval ending at a semantic instant,
the n-th moment of everything on that axis shall equal a target

flow compensation   target 0
velocity encoding   target set by venc
```

And they have **different owners**, because the two halves differ in what they act on. The
Handbook makes this structural rather than a judgement call: it puts them in **different
chapters**, by the same author.

```text
Ch. 9  Motion-Sensitizing Gradients   9.2  Flow-Encoding Gradients      -- puts information IN
Ch. 10 Correction Gradients           10.4 Gradient Moment Nulling      -- takes an artefact OUT
```

That the requirement model is shared is also the source's own: §9.2 and §10.4 contain the *same
sentence* about per-axis independence, word for word.

## 2. Velocity encoding is an ordinary encoding leaf

A bipolar pair, "two lobes of equal area and opposite polarity", net area zero so stationary spins
accumulate no phase. The public quantity is **Δm1, the change in first moment between the two
toggles**, through `VENC = π/(γ|Δm1|)` — a caller who sets m1 rather than Δm1 is out by two.
Closed form, per axis, needs nothing from outside, correctness checkable on its own events. It
belongs beside `PhaseEncode`, in a folder whose contract — "gradients, no ADC, imposing a phase
you intend to sample" — describes it exactly.

## 3. Flow compensation is not a Module, and the source says so in its first sentence

> "Gradient moment nulling … is the process of **modifying a gradient waveform** … The
> modifications can include changing the amplitude, duration, shape, or **number of gradient
> lobes** that make up the waveform." (§10.4, p. 331)

The waveform it modifies belongs to `CartesianLine`, `Excitation` or `PhaseEncode`. A leaf may not
rewrite another leaf's events, so a separate insertable "flow compensation" block is not the
abstraction — a **reshaped** readout or slice-selection waveform is. Hence `EXTEND_EXISTING`: a
compensation option on the leaf that owns the waveform, which can already compute its own moments.

The cross-leaf case — velocity-compensating a slice-selection waveform, which needs the first
moment the excitation already accumulated **from the RF isodelay point** — is kernel work by the
ownership table's own definition.

What a leaf should expose for that is **not a cached moment.** A moment depends on the placement
time, the semantic origin, the endpoint and the order, so `m1` is not a property a leaf can
publish truthfully. The faithful form is that a leaf exposes its *semantic instants* and its
*events* — the RF centre or isodelay instant, the echo instant, the gradients it intends — and the
joint owner integrates what it needs over the interval it cares about. Moment analysis is an
internal helper, not an interface every leaf must implement.

## 4. The correction that matters most — it is closed form

**Phase B's inference was too strong, and this is the main result of C3.**

Phase B saw Open4DFlow state a residual moment requirement and delegate to `gropt`, and concluded
SeqCraft needs a way to express requirements before realisation. The Handbook reduces the hard
canonical case — the velocity-compensated slice-selection waveform, accounting for `ms`, the slice
lobe's first moment from the RF isodelay point — to

> *"a quadratic equation for the unknown lobe area"* with an explicit physically significant root
> (§10.4, p. 343, Eqs. 10.70–10.73)

Target, neighbour's contribution, hardware limits in; lobe areas and widths out. **Algebra, not
search.**

The correction has to be stated precisely, because it is easy to over-shoot:

```text
WRONG       Open4DFlow uses GrOpt, therefore SeqCraft needs a general solver layer
CORRECT     the canonical cases do not justify introducing a general solver layer
ALSO WRONG  coupled gradient design never needs an optimizer
```

The Handbook gives closed forms for the canonical constructions and **does not establish them for
the general problem** — several moment orders at once, several semantic windows, a minimum-TE
search, fixed waveform segments, non-zero endpoints, arbitrary pre-existing waveforms, or
eddy-current, concomitant-field and PNS constraints. A numerical optimizer remains a legitimate
*realization backend* there. What the evidence rules out is founding SeqCraft's public semantics
on one.

So the questions the earlier directive asked, answered:

| | |
|---|---|
| **known before the solve** | the target moment and order, the hardware limits on the axis, the time available |
| **depends on the surrounding sequence** | `ms` from the excitation; the readout's own moment to the echo; where the echo is |
| **what would the candidate output?** | **events** — for the canonical cases. Not moments, not requirements. A requirement interface is needed only for the general case, which is not proposed |
| **where v0 stops being expressive** | see [`gap_note.md`](gap_note.md) — and it is narrower than "needs a solver" |

One simplification worth keeping: once m0 = 0, the first moment is **independent of the chosen
time origin** (p. 343, citing Simonetti 1991). So a moment requirement need not agree a clock,
provided the zeroth moment is already nulled.

## 5. The convergent finding

The two families **do** compose, and the efficient composition is the thing SeqCraft cannot
express. §9.2.3, p. 290: appending the bipolar to an already-compensated imaging waveform gives
five lobes, and *"the resulting five lobes can be combined into three"*, which *"always reduces the
minimum TE"*.

Merging means the efficient waveform spans territory currently split between `Excitation` or
`CartesianLine` and the encoding. **That is the same limitation C2 measured**, from the opposite
direction:

```text
C2  the rewinder of repetition n and the prephaser of n+1 are one lobe    220 us / axis / TR
C3  the encoding lobe and the compensated imaging lobe are one lobe       ~0.4-0.9 ms
```

Two fine scans, two sequence families, neither looking for it, same cause:

> **The most efficient physical waveform may span boundaries that SeqCraft currently assigns to
> separate leaf/module instances.**

Put that way rather than "one abstraction cannot own a waveform several leaves emit", because the
preferred resolution is that those leaves never emit the contested pieces — a joint owner holds
their coupled design **before events are materialised**, which keeps `Module → LogicBlock → tree →
compiler` intact and adds nothing between the layers.

The two cases are not identical, and Phase E should not merge them. C2's may be a case where two
already-decided events could be fused with **no change of physical semantics** — a candidate for a
compiler optimisation. C3's is not: its merged shape depends on VENC and a moment target, so it
must be designed with the physics in view. Recovering C2's 220 µs by teaching the compiler about
bSSFP would be the wrong fix.

## 6. v0

**Handled well.** The ownership test produced the asymmetric answer cleanly and without a
tiebreaker — it is what said flow compensation has nothing to own. `EXTEND_EXISTING` was already in
the status vocabulary and is exactly right. The licence discipline held: Open4DFlow stayed
discovery-only, the physics was re-derived from the Handbook, and nothing was adapted.

**Reasoned outside v0.** Three things, all in [`gap_note.md`](gap_note.md): one record holds one
status and this candidate was two; the ownership model has no vocabulary for a kernel that
*jointly realizes* a waveform from its leaves' physical facts rather than merely *placing* them;
and the acceptance claim that matters here is a sum across several modules' events, which rules E
and F, both per-module, do not reach — and which the validator must integrate for itself rather
than ask a module to report.

## 7. Did domain evidence change the outcome again?

**Yes — third for three, and this time it argued against building something.**

```text
C1  supplied a family that one code witness could not establish
C2  corrected a contract that three code witnesses agreed on
C3  corrected an architectural inference drawn from the only code witness there is
```

C3 is the most consequential. A code-only reading said "a real implementation calls a solver,
therefore we need a requirement abstraction". The domain source said the canonical cases are a
quadratic. Had Phase B's reading gone straight to design, SeqCraft would have acquired a
constraint-solver layer for a problem that has a closed-form answer — and would still have had the
merged-waveform problem, which is what actually blocks both C2 and C3.

That is the strongest case yet for `evidence_class: domain-reference`, and it is still recorded by
hand rather than in the schema.

## 8. What implementing `VelocityEncode` added, 2026-09-24

The scan's physics held. Three things came out of building it that the record could not have.

**One acceptance bullet was false, and the spec-consistency pass caught it.** The record said the
two toggles are identical "in duration, zeroth moment, and moments above the first". They are not:
one toggle is the other's negation, so the second moment negates with it. The only moment the two
share is the zeroth, and they share it because it is zero. The bullet is corrected in place, and
the test now asserts both the true invariant *and* that `m2` is not equal, so the false version
cannot return unnoticed.

**The open API question was for a human, and it was answered as one.** A single instance owns
`venc_m_s`, the waveform design and the pair invariant; `build(polarity=+1|-1)` selects the
encoding state, matching `PhaseEncode`'s design-once/build-per-acquisition idiom. The name is
`polarity` rather than `sign`, and it means **the sign of the emitted first moment** — not of a
velocity, and not of a reconstructed phase.

That definition earned itself immediately. In Layer 3 the measured phase difference comes back
*negative* where the textbook relation is positive, which is either a sign error in the module or
the simulator's phase convention — and the flow measurement alone cannot distinguish them. A
**static control**, a stationary spin at a known offset under a deliberately un-nulled zeroth
moment, shows the same `−1` ratio. So it is the simulator's convention throughout. Had `polarity`
been defined on velocity or on reconstructed phase, the module would have had to take a position
in that argument; defined on the first moment, it never does.

**Layer 3 was also a check on the instrument.** Every Layer 1 assertion here is made with
`sc.moments(..., 1)`, whose first-moment path had little independent exercise before this
candidate. A simulation that moves spins reaches `delta_m1` by a completely different route, and
the two agree to 1e-5 — so the analyser is corroborated along with the module. That is the rule
the `b_value` correction produced, applied prospectively for the first time rather than after a
defect.

Unchanged: the closed form is algebra, exactly as §4 argued. No moment-requirement IR and no
numerical optimiser were needed or added, and `flow_compensation.yaml` is still the record that
carries the hard part.
