# What generalised, what stayed family-specific, what still needs a human

Skill v0's summary judgement after three pilots. The useful content is the middle column: the
things that did **not** generalise are what stops this becoming a sequence-wrapper generator.

---

## 1. What generalised cleanly

**The pipeline order.** Every pilot ran corpus → discovery → evidence → contract → ownership →
claim → adapter → extraction → validation → decision, and in each case the ordering did real work.
The one ordering constraint that carried the most weight: **no extraction before the acceptance
claim is written down.** Radial's rotation-equivariance test was defined against the official
reference and passed there *before the module existed*, so it cannot be encoding the module's
behaviour. That property is only available if the order is respected.

**The ownership test.** *If a component cannot determine the correct value using only its own
physical contract and parameters, it should not own that event.* This separated every disputed
case across three families and never needed a tiebreaker. It is sharper than "avoid duplication"
and occasionally disagrees with it — which is the point.

**Provenance and licence recording.** Mechanical, and it never once caused an argument.

**The four rules, and the reason each is kept with its failure.** A rule stated alone gets
optimised away by the next person who finds it inconvenient; a rule stated with the sequence that
compiled legally and selected nothing does not.

**Stop outcomes.** `NO_NEW_MODULE` needed to be a first-class result before the first candidate,
not after a disappointing one. Two things were correctly not created — `RadialGRETR` and
`seqcraft/trajectory/` — and neither felt like failure at the time because the vocabulary had room
for them.

**The thirteen-field core.** Derived by intersection, not design. `target.untouched` surviving all
three is the mildly surprising one: in practice the load-bearing sentence is as often what the
candidate does *not* own as what it does. `Refocusing` appears in all three untouched lists.

---

## 2. What stayed family-specific, and must

**The acceptance criterion.** Three pilots, three criteria, no overlap:

| | | why nothing else would do |
|---|---|---|
| TSE | event identity | it was promoting an existing implementation; anything weaker would not have proved the extraction changed nothing |
| Radial | trajectory agreement | two legitimate constructions of the same trajectory share **no events** |
| GRE3D | lattice + signed z physics + timing | the reference has no DC sample where `CartesianLine` puts one, and the selective path has no executable reference at all |

Freezing any one of these into the tooling would have broken the other two. This is why
`SequenceFingerprint` is not an equality rule: generic measurement primitives → family-specific
check set → structured `DifferenceReport`.

**The interior of `validation`.** The strongest single piece of evidence in the whole exercise:
`validation` appears in all three records and **not one of its sub-keys does**. Three acceptance
claims needed three disjoint shapes of evidence. The schema therefore requires the section to
exist, requires `layers` inside it, and says nothing else.

**What counts as a witness.** Independence is a judgement about a corpus, not a property of a
file. Two ported implementations of one design are one witness. A workflow can *require* the
judgement; it cannot make it.

**Whether Layer 3 applies.** For GRE3D it was cheap and it caught three real mistakes. For Radial
it requires infrastructure that does not exist and must not be built to satisfy a checklist. Same
ladder, opposite correct answers.

**Tolerances.** Every one in the record was chosen against a measured floor for that family.
There is no portable number.

---

## 3. What still needs human judgement

These are not gaps to close. They are the reasons this is a **supervised** workflow, and an
attempt to automate any of them would be the thing that makes it dangerous.

**Is this one contract, or two, or none?** The skill can force the question, marshal the evidence
and refuse to proceed without an answer. Nothing in it answers "is a radial GRE a kernel." That
took a human looking at four lines of composition and deciding they were not a shared physical
solve.

**Is the reference wrong?** Three references disagreeing may mean one is wrong, or that all three
encode a historical constraint nobody remembers. `REFERENCES_DISAGREE` records the disagreement
and stops; picking a winner to make progress is exactly the failure it prevents.

**Is the comparator wrong?** Four times out of four, when a measurement disagreed with a
candidate, the measurement was at fault. The workflow can insist on known-answer tests; it cannot
generate the suspicion that an instrument is lying. This remains the sharpest open weakness and is
recorded as such in `conformance_report.md` §5.

**Does the abstraction earn its name?** `WRAPPER_ONLY` has a test — does extraction shorten the
notebook, and is there a value the wrapper determines that its caller could not — but applying it
requires taste about what a caller should have to get right twice.

**When is a compiler change real?** The escalation makes the claim expensive and orders the
evidence. The final step — characterising the numerical floor from something other than the
candidate that wants the tolerance widened — is where judgement concentrates, and it is about to
matter for Spiral.

---

## 4. Readiness

Skill v0 is ready to be **reviewed**, not yet to be trusted on a new candidate.

What supports it: the workflow recovered every boundary decision, every acceptance claim and every
safety boundary across three pilots; the schema carries a real pilot with nothing lost; and it
produces three artifacts the pilots did not have — a dependency map with the docstring-mention trap
removed, per-layer validation status, and structured independence.

What qualifies it: the schema has been validated against three candidates from two sequence
families, all Cartesian-or-radial, all GREEN, none of which exercised `NO_NEW_MODULE`, `RED`, or a
compiler-change claim. **Every stop outcome is currently untested.** The first genuinely new
candidate is as much a test of the skill as of the candidate, and should be treated that way.

Recommended first target is a candidate with a plausible chance of stopping. `SaturationPrep` is a
good choice partly *because* it may well turn out to be `EXTEND_EXISTING` or `NO_NEW_MODULE`, and
watching the workflow decline to produce a Module would be worth more than another GREEN.
