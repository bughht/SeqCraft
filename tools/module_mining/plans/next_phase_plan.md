# SeqCraft Module Mining — Post-TSE Next Phase Plan

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../README.md) for which one to edit.

> **Superseded as a forward plan.** The current roadmap is
> [`post_v0_roadmap.md`](post_v0_roadmap.md). This file is retained for the reasoning it
> records, not for the plan it proposes; see [`README.md`](README.md) for the audit.

**Date:** 2026-09-18  
**Status:** Superseded — retained as the record of the Radial fine scan as it was planned  
**Current milestone:** TSE/FSE fine-scan pilot COMPLETE — GREEN  
**Next candidate:** `RadialReadout`

## 1. Current state

The TSE/FSE pilot is complete. `TSEShot`, `FSE2D`, the `CartesianLine` echo-geometry extension, reference adapters, layered comparator, parameter sweep, metamorphic checks, and boundary decision are already merged.

The final ownership rule established by that pilot is:

```text
leaf module     → intrinsic physics / waveform geometry
kernel          → coupling between multiple leaves
imaging module  → acquisition / scan policy
compiler        → Pulseq legality
```

Operationally:

> If a component cannot determine the correct value from its own physical contract and parameters, it should not own that event.

The TSE pilot should now be closed out in the planning docs as `COMPLETE — GREEN`. Historical YELLOW material should remain as history, but be explicitly marked as superseded by the later GREEN boundary decision after human review and extraction.

## 2. MATLAB scope

SeqCraft package-level Module implementations remain Python-only.

MATLAB Pulseq may be used as:

```text
physics reference
architecture evidence
optional .seq validation source
```

It is not a Module implementation target, frontend requirement, or reason to introduce cross-language runtime abstractions.

## 3. Next fine scan: RadialReadout

The next candidate should be `RadialReadout`.

Do not implement the class first. Repeat the evidence-first workflow:

```text
references
→ code archaeology
→ semantic decomposition
→ invariant table
→ thin reference adapter
→ comparator extension
→ candidate boundary
→ extraction
→ parameter sweep
→ GREEN / YELLOW / RED
```

Primary evidence should include:

```text
official PyPulseq radial GRE
official Pulseq MATLAB radial GRE
an independent radial implementation, preferably OpenMRF RAD or equivalent
```

The official PyPulseq implementation is the primary executable reference.

## 4. Why Radial is the right second candidate

TSE/FSE tested:

```text
kernel ownership
multiple-leaf coupling
cross-axis moments
echo-train timing
shot vs acquisition boundary
```

RadialReadout tests a different part of the architecture:

```text
readout leaf ownership
non-Cartesian k-space semantics
trajectory rotation
k-space centre semantics
ADC/gradient synchronization
full-spoke vs centre-out ownership
```

If the same fine-scan process works for both, the workflow has much stronger evidence of being general rather than TSE-specific.

`GRE3DTR` remains important, but it is structurally closer to existing `GRE2DTR` and therefore a weaker test of workflow generality.

## 5. Main Radial ownership questions

The fine scan should answer these before extraction:

1. What is the natural reusable unit?
   - full spoke;
   - centre-out spoke;
   - one generic radial readout.

2. Do full-spoke and UTE-style centre-out trajectories share one clean physical/API contract?
   - If yes, one `RadialReadout` may support both.
   - If not, prefer separate abstractions rather than forcing a mode flag.

3. Does the readout itself own angle/orientation?
   - Likely yes if it can determine the complete trajectory from its own parameters.

4. Does the readout itself own semantic `k=0` timing/location?
   - Likely yes. Higher layers should not reverse-engineer k-space centre from waveforms.

5. What remains acquisition policy?
   Likely:
   - spoke ordering;
   - golden-angle schedules;
   - number of spokes;
   - interleave/view-sharing policy.

Those should stay above the leaf.

## 6. Code archaeology

For each reference, trace:

```text
user parameters
→ gradient design
→ prephasing/start position
→ ADC timing
→ trajectory orientation
→ k-space centre placement
→ spoke ordering policy
```

Classify each item as:

```text
intrinsic readout physics
candidate parameter
caller/acquisition policy
Pulseq representation artifact
```

Do not let manual Pulseq block structure define the candidate API.

## 7. Radial invariant table

Write the invariant table before the final API.

At minimum measure:

```text
requested trajectory angle
actual k-space angle
k-space extent
Δk / sample spacing
ADC-gradient synchronization
semantic k=0 location
gradient limit
slew limit
full-spoke symmetry when applicable
```

Also record, where meaningful:

```text
prephaser moment
readout total moment
ADC centre sample index
centre timing relative to gradient
```

## 8. Radial metamorphic tests

A primary metamorphic test should be:

```text
RadialReadout(angle = φ)
≈ rotate(RadialReadout(angle = 0), φ)
```

Compare:

```text
Gx/Gy
k-space trajectory
ADC k-space samples
```

Additional useful checks:

```text
angle + π → same physical line with reversed orientation, if full-spoke
negating angle → reflected trajectory
FOV change → expected Δk scaling
matrix/sample-count change → expected extent/density change
```

These tests should complement, not replace, reference comparisons.

## 9. SequenceFingerprint evolution

Reuse the TSE structure:

```text
L0 metadata
L1 exact event digest when meaningful
L2 physical waveform timeline
L3 acquisition semantics
L4 family-specific physics
```

For Radial, L4 should emphasize:

```text
trajectory geometry
orientation
centre semantics
extent
rotation equivariance
```

Extend the comparator only when the Radial scan exposes a real missing concept.

Do not redesign the entire comparator in advance.

## 10. Do not add `seqcraft/trajectory/` yet

The coarse scan suggested a possible future layer:

```text
seqcraft/trajectory/
    radial.py
    spiral.py
    rewinder.py
```

Do not create it from Radial alone.

Decision rule:

```text
Radial only
    → keep design local if coherent

Radial + Spiral reveal the same numerical-design contract
    → consider a shared trajectory utility layer
```

A utility layer should appear only after at least two real consumers define a stable contract.

## 11. Radial GREEN criteria

`RadialReadout` may become GREEN when:

```text
candidate boundary is clear
at least one executable Python reference exists
provenance/license is recorded
waveform checks pass
trajectory invariants pass
rotation metamorphic tests pass
parameter sweep passes
package tests pass
no compiler change is required
```

The sweep should cover multiple:

```text
angles
matrices
FOVs
bandwidth/dwell settings
full-spoke / centre-out modes, if both are supported
```

If full-spoke and centre-out fail to share a clean contract, that is evidence for separate abstractions, not a failed fine scan.

## 12. What happens after Radial

If Radial reaches GREEN, compare the two pilots before doing broad conversion:

```text
TSE/FSE        → kernel-heavy Cartesian echo train
RadialReadout  → non-Cartesian readout leaf
```

Ask:

```text
Which ModuleCandidate fields mattered in both?
Which adapter concepts survived both?
Which comparator layers were reusable?
Which report fields helped human review?
Which ownership questions repeated?
```

Only then consider formalizing:

```text
ModuleCandidate schema
ReferenceSequence interface
SequenceFingerprint API
fine-scan report schema
module-miner skill workflow
```

The stable part should be workflow/evidence structure. Family-specific MRI physics checks must remain extensible.

## 13. Candidate order after Radial

Recommended order:

```text
1. RadialReadout
2. GRE3DTR
3. SaturationPrep
4. DiffusionPrep / DiffusionEncoding
5. SpiralReadout
```

Why:

- `GRE3DTR`: broadens 3D Cartesian coverage with relatively low architectural risk.
- `SaturationPrep`: tests preparation semantics and avoids abusing `Excitation`.
- `DiffusionPrep / Encoding`: major architecture stress test involving b-value/b-tensor, cross terms, and refocusing context; defer until the workflow is more stable.
- `SpiralReadout`: best point to decide whether non-Cartesian trajectory design deserves a shared utility layer.

## 14. Immediate execution plan

The next local session should:

```text
A. update the module-mining plan:
      TSE/FSE pilot = COMPLETE — GREEN

B. keep historical YELLOW evidence,
   but link it to the superseding GREEN boundary decision

C. create a new Radial fine-scan branch

D. collect the primary Radial references

E. perform code archaeology

F. write semantic decomposition

G. write the invariant table

H. decide whether full-spoke and centre-out share one contract

I. build the smallest Python reference adapter

J. extend the comparator only as required

K. add rotation/metamorphic tests

L. only then propose/extract `RadialReadout`

M. run a parameter sweep

N. record GREEN / YELLOW / RED

O. compare TSE and Radial workflows before formalizing tooling
```

## 15. Human-review checklist

For each future candidate ask:

```text
Does it represent one coherent MR concept?
Can it determine its own events from its own contract?
Is acquisition policy being pushed too low?
Is coupling between leaves being pushed too low?
Is the compiler being asked to understand MRI semantics?
Is a utility layer being created before multiple consumers define it?
Does the candidate shorten real sequence code without hiding important choices?
Can preservation be demonstrated independently of the new code's own tests?
```

## 16. Success condition

This phase succeeds when:

```text
TSE/FSE = GREEN
RadialReadout = GREEN
```

and the same evidence-first workflow has been shown to work for both a kernel-heavy Cartesian echo train and a non-Cartesian readout leaf.

Only then should the module-mining process begin to harden into a durable AI-assisted skill.
