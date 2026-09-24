# T2 preparation — Phase C first fine scan

**Date:** 2026-09-22 (revised) · **Skill:** v0, frozen and unmodified
**Result:** `NEW_LEAF / APPROVED_FOR_IMPLEMENTATION`, for one mode
**Implemented** 2026-09-23. A design correction followed the first implementation review; the scan
record below is left as it was written and [§7](#7-the-correction-2026-09-23--realisation-is-not-contract)
records what it got wrong.

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

## 7. The correction, 2026-09-23 — realisation is not contract

The scan above is accurate about what each witness *does*. It is wrong about what that makes
**ours**. The mode contract it wrote carried `mrseq`'s composite `270x / −360x` tip-up as the
mode's only form, and the first implementation then had to invent a rule — an amplitude-weighted
group centre — to say where `prep_time_s` ended inside an asymmetric composite. Both are the same
error, made once and then compounded:

```text
the implementation corpus tells us    "someone implemented it this way"
domain evidence tells us              "this is what the abstraction actually requires"
```

A special realisation from one corpus implementation must not be promoted into the reusable
physical contract merely because it is the permissively licensed executable witness. That is
exactly what happened here, and the reason it was easy to miss is that `mrseq` was the *only*
buildable witness, so its choices arrived looking like the only choices there were.

### The evidence, reclassified

| | what it is | what it therefore settles |
|---|---|---|
| **D1** Handbook §17.4 | domain reference | the **physical contract**: tip in, refocus, tip back, spoil — and nothing more |
| **E1** MRIquestions | educational synthesis *(new class)* | the **canonical picture** the field teaches: `90x → refocus → −90x → spoil` |
| **R1/R2** OpenMRF MLEV | design witness | the **refocusing structure**: MLEV-4, `+ + − −`, which says nothing about the tip-up |
| **R3** `mrseq` | design witness | an **established realisation variant** of the tip-up, and its own timing convention for it |

E1 is registered with `roles: [educational-synthesis]` and `explicitly_not: [design-witness,
oracle]`: it identifies the canonical form of a family, and **it must not override stronger domain
evidence**. It is cited here because it is the second independent statement that the simple `−90`
is the canonical tip-up, which is what makes the composite a realisation rather than the contract.

### What shipped instead

The default is the canonical `+90 → MLEV-4 → −90 → spoiler`, with `prep_time_s` from the `+90`'s
centre to the `−90`'s centre. The composite is `tip_up='composite_270_360'`, with `prep_time_s`
from the `+90`'s centre to the **270**'s centre — `mrseq`'s own convention for `mrseq`'s own
realisation, adopted rather than replaced by a rule of ours. The amplitude-weighted-centre rule is
removed: neither the domain evidence nor the witness establishes it, and a Layer 3 slope fit
cannot adjudicate it, because a constant offset does not change a slope.

Layer 3 then measured the two under B1 error and found **no consistent advantage for the
composite** — both become unreliable away from nominal B1, their errors are not identical, and at
B1 = 0.8 the composite is notably the worse of the two. That is enough for the one conclusion this
record needs: the composite must not be presented as a robustness feature. It is *not* a finding
about which part of the preparation dominates B1 sensitivity, and none was sought — the two
tip-ups are equivalent as ideal rotations, and the emitted sequence, with finite pulse durations
and real gaps between them, is not obliged to agree with that. Had the composite been promoted
into the contract as a robustness feature, the measurement would not have backed it.

### Where the lesson lives

In the fine-scan playbook, as *"Realisation is not contract"*, with the three questions to ask of
any detail before it becomes an invariant — not as a rule local to this candidate. `sources.yaml`
gained the `educational-synthesis` role and the `mriquestions` entry at the same time, and
`schema.py` gained the class.

The timing convention in §2 was reconciled *before* any code was written. This one was not caught
until the module existed and its docstring had to explain itself — and the corrective evidence was
in the record all along: §1 records D1's simple `−90` and then treats it as a base form the mode
had departed from, rather than as the contract the mode should have kept.
