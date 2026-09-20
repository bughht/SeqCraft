# Outcomes, including all the ways to stop

The workflow must be able to end without a Module. A corpus that contains no reusable physical
contract is a normal result, and reaching `NO_NEW_MODULE` at the contract stage is far cheaper
than reaching `RED - WRAPPER_ONLY` after an extraction.

Every outcome here is a **recommendation to a human reviewer**. The skill never promotes, never
merges, and never resolves its own YELLOW.

---

## Coverage classification — the coarse scan, before any decision

A coarse scan asks a different question from a fine scan. Not *should this be a Module* but
*what already covers it*, and the two claims below are not the same:

```text
"SeqCraft can already build it"
"SeqCraft already ships the right reusable abstraction"
```

Collapsing them into the word "duplicate" hides the gap between them, which is where
discoverability problems and missing abstractions live. Use:

| classification | meaning |
|---|---|
| `DIRECT_SHIPPED_DUPLICATE` | an equivalent **public shipped** abstraction already exists |
| `DEGENERATE_CASE_OF_EXISTING_MODULE` | the family is a parameter limit of a shipped abstraction |
| `COMPOSITION_COVERED` | existing modules express the acquisition, but no equivalent public abstraction exists |
| `NOTEBOOK_ONLY_EXISTING` | the abstraction exists only in an example or reference notebook |
| `PRIMITIVE_COMPOSITION_COVERED` | the sequence is small enough that plain composition **is** the intended API |

Two of these are not verdicts and must not be read as closing a family:

- `DEGENERATE_CASE_OF_EXISTING_MODULE` may travel with `ARCHITECTURE_REVISIT_CANDIDATE` — the
  parameter limit works, and the naming or the layering may still be wrong.
- `NOTEBOOK_ONLY_EXISTING` says where the code is, not whether it should be promoted. A notebook
  class that no other notebook needs is doing its job.

`PRIMITIVE_COMPOSITION_COVERED` is the one that keeps a batch scan honest: **a sequence can be
fundamental and common without deserving its own Module.** The test is whether there is any
cross-component solve for a class to own — an excitation, a delay and an ADC composed in order is
not one. A scan that turns every textbook sequence name into a class has stopped measuring
anything.

Which families currently sit here is **not** recorded in this file. See the coarse scan.

---

## GREEN

All of:

- an executable reference exists, or the acceptance claim explains why none is possible and what
  carries the proof instead;
- the candidate boundary is clear and the ownership question is answered per event;
- the acceptance claim was written **before** extraction, and both halves of rule C are stated;
- the family-specific check set passes at the stated tolerance;
- the parameter sweep passes;
- the relevant metamorphic checks pass;
- emitted-sequence inspection has been done for every promoted physical mode (rule B);
- shared-leaf dependency impact is mapped and every consumer classified, with no
  `UNKNOWN / NEEDS_REVIEW` remaining (rule D);
- the SeqCraft package suite passes;
- no compiler changes;
- provenance and licences recorded;
- per-layer validation status recorded, with a reason for any deferral.

```text
-> candidate may enter PR review
```

GREEN is scoped to the layers it covers. A Layer-1+2 GREEN is not an end-to-end claim and must not
be summarised as one.

---

## APPROVED_FOR_IMPLEMENTATION

The design is settled and the module does not exist yet.

The boundary, the ownership of every parameter, any mode contract rule A requires, and any policy
that departs from a reference have all been reviewed and accepted. What remains is writing it.

```text
-> implementation may begin, against the contract as recorded
```

**GREEN cannot express this**, and the distinction is not bureaucratic. GREEN means the checks
passed; here there is nothing to check, and a record claiming it would have to either fake an
emitted-sequence inspection or leave rule B unsatisfied. **YELLOW cannot express it either** —
YELLOW says stop and ask a human, and the human has already answered.

Two things this state must carry, or it is premature:

- every parameter's ownership decided, with any mode contract **written**, not planned;
- any **durable shape** the first implementation must respect from its first line, so that
  implementing a subset does not later change the meaning of what was built.

That second point is what the state is for. A design approved in the abstract and implemented
incrementally will encode the first increment's assumptions unless the full contract is on paper
before anyone starts.

## NO_NEW_MODULE

The corpus was searched and no reusable physical contract was found — the similarity is at the
level of names, protocols or house style rather than physics, or the shared part is already
expressed by existing modules plus caller policy.

```text
-> record the evidence and the reasoning; propose nothing
```

The record matters: the next scan over the same corpus should not repeat the work. Say what was
looked at and what the similarity turned out to be.

This is not failure. `RadialGRETR` and `seqcraft/trajectory/` were both considered and both
correctly not created; a radial GRE is `Excitation -> RadialReadout -> spoiler -> delay`, which is
composition, and the angle schedule is caller policy.

---

## YELLOW — needs human MRI or architecture review

Physics or boundary is unresolved in a way more measurement will not settle. Stop and ask.

### PHYSICAL_BOUNDARY_UNCLEAR
The physics is settled but ownership is not: two layers could each plausibly own an event, and the
"can it determine the value from its own parameters" test does not separate them. TSE sat here
while the physics was complete — the open question was whether the reusable unit included the
excitation.

### MODE_CONTRACT_INCOMPLETE
The boundary is settled and the physics is understood, but a candidate with named modes has not
published the table rule A requires. Distinct from `PHYSICAL_BOUNDARY_UNCLEAR`: nothing is
*unknown*, something is **unwritten**, and the two need different work — one needs analysis, the
other needs a document. Using the boundary code for a missing contract overstates the uncertainty
and hides that the remaining task is bounded.

The table is what proves named modes are modes of one thing rather than several implementations
behind a string flag, so it is owed **before** implementation and not after.

### REFERENCES_DISAGREE
Two references produce materially different sequences for the same nominal protocol, and it is not
yet known whether the difference is representational, historical, or a genuine physical
disagreement. Record the difference; do not pick a winner to make progress.

> **References disagreeing on protocol parameters does not imply disagreement on the abstraction
> boundary.** Escalate to this code only when the sources disagree about *what the module
> fundamentally owns*.

Two independent designs choosing different flip angles, durations or bandwidths, while agreeing on
the structure, the sign conventions and what must follow what, have not called the boundary into
question — they have shown you which values are **caller-owned protocol**. That is information, not
an obstacle. The boundary is in doubt when the sources disagree about which layer determines a
value, what the module must emit, or what may not be omitted.

### REFERENCE_NOT_INDEPENDENT
The references that agree turn out to be one witness — a port, a fork, or a shared house style
with identical arithmetic. Agreement among copies is not corroboration. This is why the `Evidence`
artifact requires an `independence` field whenever more than one reference is cited.

### VALIDATION_ORACLE_UNTRUSTED
The instrument that would produce the evidence has not itself been validated. An image through an
unreviewed reconstruction is two unvalidated things agreeing. The correct move is to defer that
layer with the reason recorded — not to lower the claim quietly and not to build the oracle as a
side quest.

Live instance: `RadialReadout` Layer 3, deferred until the non-Cartesian reconstruction path is
independently reviewed with **two** real consumers, because a reconstruction contract with one
consumer is that consumer's implementation under a more general name.

### SHARED_LEAF_IMPACT_UNRESOLVED
Rule D ran and at least one consumer is still `UNKNOWN / NEEDS_REVIEW`.

### COMPILER_CHANGE_CLAIM_NEEDS_MINIMAL_REPRODUCER
A compiler defect is claimed without a candidate-free reproducer proving a legal physical input is
rejected or altered. See the escalation in `rules.md`.

```text
-> human MRI / architecture review
```

---

## RED — do not promote

### WRAPPER_ONLY
The abstraction removes no meaningful duplication and owns no arithmetic anyone could get wrong.
Test: does extraction shorten the notebook, and is there a value the wrapper determines that its
caller could not? If the answer to either is no, it is a name for a sequence, not a module.

### DUPLICATES_EXISTING_MODULE
The contract is already expressed by an existing module, possibly under a different protocol name.
A multi-echo gradient echo is a spoiled GRE whose readout reads the line more than once; it got no
class for exactly this reason.

### EXTRACTION_CHANGES_THE_PHYSICS
The extracted module does not reproduce the reference's physics, and the difference is physical
rather than representational. A module that cannot be extracted without altering the sequence is
not a module.

Also RED: the reference cannot be reproduced at all; licensing or provenance is unresolved; the
candidate requires source-specific compiler behaviour.

```text
-> do not promote
```

---

## Recording an outcome

```yaml
status: NEW_KERNEL
traffic_light: YELLOW
reason_code: PHYSICAL_BOUNDARY_UNCLEAR
decided_by: human
```

`reason_code` is required for YELLOW and RED and must come from the catalogue above. A YELLOW
without a reason code cannot be triaged and cannot be closed — it is just a feeling about the
candidate.

When a YELLOW is later resolved, keep the superseded status in the record rather than overwriting
it. TSE's `candidate.yaml` keeps its `MORE_EVIDENCE_REQUIRED / YELLOW` note beside the final
GREEN, and that note is the more instructive half: it shows the physics was settled long before
the boundary was.
