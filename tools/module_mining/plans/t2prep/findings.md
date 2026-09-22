# T2 preparation — Phase C first fine scan

**Date:** 2026-09-22 (revised) · **Skill:** v0, frozen and unmodified
**Result:** `NEW_LEAF / APPROVED_FOR_IMPLEMENTATION`, for one mode

The conventional-control candidate, chosen to answer *"does the ordinary reusable-leaf path still
work cleanly?"* — and it does. Nothing has been implemented; the record is
[`candidate.yaml`](candidate.yaml) and the domain card is
[`../../domain_evidence/t2prep.md`](../../domain_evidence/t2prep.md).

## What changed, and why

The first pass of this scan reached `YELLOW / MODE_CONTRACT_INCOMPLETE`, and part of the reasoning
was wrong. It treated **one executable repository** as thin evidence for the family itself. That
conflates two different things:

```text
does this physical abstraction exist, with canonical semantics?   -> domain evidence
does OUR realisation match a particular implementation?           -> code witnesses
```

Reference coverage limits how strongly we may claim equivalence to an implementation. It does not
decide whether a mature abstraction exists. Under the corrected standard the boundary was never
the problem, and the two questions that genuinely blocked are both now answered.

## 1. Evidence, in three classes

**Domain reference (D1).** *Handbook of MRI Pulse Sequences*, §17.4 "Driven Equilibrium", pp.
888–895; §17.4.2 at pp. 893–895. The Handbook's own term is driven equilibrium and it says
plainly: *"This use of driven equilibrium is also called T2 preparation (T2-prep)"* (p. 888). It
gives the nonselective `90x–τ–180y–τ–90₋ₓ` structure, the weighting `Mz·exp(−Tprep/T2)`, multiple
180° pulses as the canonical way to lengthen `Tprep`, composite pulses as an optional B0/B1
robustness choice, the spoiler after the tip-up, and centric k-space ordering as *downstream*
policy. Full extraction with page-level provenance is in the card.

**Design witnesses.** Two implementations, different institutions, different languages:

| | OpenMRF `T2/` and `MLEV/` | `mrseq` `preparations/t2_prep.py` |
|---|---|---|
| licence | MIT | **Apache-2.0** |
| tip pair | adiabatic BIR-4 / AHP | **hard 90x; composite 270x + (−360x)** |
| train | 2 composites, or n×4 with a ±1 list | **MLEV-4, phase `+ + − −`** |
| public quantity | `prep_times`, edge-to-edge | **`echo_time`, centre-to-centre** |
| spoiler | 3 axes, each derated 1/√3, ramps staggered | 1 axis (z), optional |
| tests assert | — | block durations and a short-TE refusal |

`mrseq` cites Levitt & Freeman (1981) for the composite pulse and **Brittain et al. (1995)** for
the application — the same Brittain reference the Handbook names for cardiac blood–myocardial
contrast. Independent literature and independent code converging on one citation.

**Oracle.** Not used yet; `Mz·exp(−Tprep/T2)` is exactly a Bloch-simulation claim and is the
planned Layer 3.

### Measured, not read

I executed `mrseq` at seven protocols (TE 30–100 ms, 180-pulse 1 and 2 ms):

```text
realised echo time  ==  requested, to 0.00 us, at every protocol
refocusing points   ==  0.1250, 0.3750, 0.6250, 0.8750 of it   (MLEV-4, exact)
tip-up              ==  1.0000
15 RF pulses: 1 excitation + 4 composites x 3 + a 2-pulse composite tip-up
```

So the semantic timing definition is realisable **exactly**, which is the fact the whole timing
question turned on.

## 2. The timing semantics — resolved, and it was not a disagreement

The two implementations book the preparation time differently: OpenMRF tip-down-**end** to
tip-up-**start**, `mrseq` excitation-**centre** to tip-up-**centre**. They differ by one tip-pulse
duration — 2 to 10 ms against a 40 to 100 ms preparation, so it matters.

It is **two bookkeepings of one physical interval**, not two physics. The Handbook states the
weighting for *"the time between the DE 90° pulses"*, the period during which the magnetisation is
transverse. So SeqCraft takes the semantic definition and inherits neither implementation's block
arithmetic:

> **`prep_time_s` = effective tip-down rotation instant → effective tip-up rotation instant.**
> For the hard and composite pulses of the initial mode this reduces to RF centre to RF centre.

An adiabatic tip pair does *not* rotate at its pulse centre, so a future adiabatic mode must
define its own effective-rotation instant. That is written into the mode differential so adding
the mode cannot silently change what `prep_time_s` means.

## 3. Ownership — and two things the Handbook decided that code could not

The leaf owns the tip pair, the refocusing train, the timing solve and the spoiler. `caller` owns
placement, the readout, the collection of preparation times and any fitting. No kernel, no
compiler involvement. Two refinements came from the domain source:

**The input magnetisation is not the leaf's, and the Handbook says so explicitly** — `Mz` "must be
calculated by taking into account steady-state conditions; it is not simply the equilibrium
longitudinal magnetisation" (p. 893). So the leaf owns the *factor* it applies and cannot know the
quantity it multiplies. That is a clean ownership statement rather than a gap, and it sharpens the
claim: the module's claim is about `exp(−prep_time/T2)`, not about absolute signal.

**One canonical variant crosses the boundary.** The Handbook's B1-robust scheme (Fig. 17.36) puts a
dephasing lobe inside the preparation and its matching rephaser *inside the subsequent imaging
sequence*, "combined with the slice-selection rephaser". No single module can hold that
correctness condition — it is a cross-leaf coupling, so it is excluded from the initial mode and
flagged as kernel-layer work. Cross-referenced to C3, which is about exactly this class.

**`Refocusing` is still not reusable**, measured as before: it refuses `crush_cycles_per_voxel=0`,
and a crusher inside the preparation destroys what is being stored.

## 4. What the abstraction establishes, and what it does not

**Establishes**, all measurable without a scanner: interval symmetry on the emitted lattice; the
declared `prep_time_s` equal to the measured effective-rotation interval; the `+ + − −` phase
pattern read off the compiled file; zero net gradient area inside the preparation with the spoiler
moment entirely outside it; and the metamorphic relation that doubling `prep_time_s` doubles the
transverse interval while leaving pulse count, phases and spoiler unchanged.

**Does not establish:** `exp(−prep_time/T2)` weighting until Layer 3 runs — and that is the claim a
user actually cares about, so Layer 3 is not optional garnish here. **Nor B1 or B0 robustness.**
That claim belongs to the deferred adiabatic mode. The composite pulses in the initial mode do buy
some insensitivity — the Handbook says so — and no quantitative claim is made about it. Nor
absolute signal, per the ownership note above.

## 5. Recommendation

**Implement `mlev4-hard`, and only that.** Boundary, ownership, semantic timing, mode contract and
acceptance claim are all written and the mode has a permissively licensed witness that I executed
and measured. Deferred with triggers: the adiabatic modes, a refocusing-count API, the G1/G2
variant, MRF scheduling, and the three-axis crusher policy.

Corpus completeness is not a promotion condition, and four deferred modes do not block the first.

## 6. Skill observation — recorded, v0 unchanged

> Mature MRI physical abstractions may be established by authoritative domain evidence even when
> executable implementation witnesses are sparse. Code witnesses constrain implementation-specific
> claims; they do not determine the existence of the physical abstraction.

Concretely, v0 lacks an **evidence class**. Its `evidence[].role` vocabulary describes a
reference's relationship to *other references* and every value presumes the reference is an
implementation; `authoritative-external` is the nearest fit and means something else. This record
carries a hand-added `evidence_class: domain-reference` plus a `curated_card` pointer.

It changed the outcome here, which is why it is worth recording: the same evidence read through
v0's implementation-only lens gave YELLOW, and read with a domain class gives
`APPROVED_FOR_IMPLEMENTATION` with the code witnesses doing the job they are good for — fixing the
mode and the timing convention.

**Whether it earns formalisation depends on C2 and C3.** One case is not enough.
