# T2 preparation — Phase C first fine scan

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-21, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../README.md) for which one to edit.

**Date:** 2026-09-21 · **Skill:** v0, frozen · **Result:** `UNDECIDED / YELLOW / MODE_CONTRACT_INCOMPLETE`

The conventional-control candidate. The question it was chosen to answer was *"does the ordinary
reusable-leaf path still work cleanly after all the difficult calibration work?"* — and the answer
is **mostly yes, and the workflow stopped it for a reason that is about evidence rather than about
the workflow.**

Nothing was implemented. The record is [`candidate.yaml`](candidate.yaml).

---

## 1. What the evidence turned out to be

Searched all nine registered corpora. `pypulseq`, `pulseq-matlab`, `pulseq-cest`,
`PulseqDiffusion` and `Open4DFlow` have no T2 preparation at all. The MRzero playground appeared
to match on `MLEV` in three notebooks; the hits were case-insensitive fragments of base64 image
data and were discarded after checking.

**So: one independent witness, `OpenMRF`, and two designs inside it.**

| | `src_preparations/T2/` | `src_preparations/MLEV/` |
|---|---|---|
| tip-down / tip-up | adiabatic, BIR-4 default | adiabatic, AHP default |
| refocusing train | exactly **two** composite pulses | `n_mlev` × **four**, with an explicit ±1 cycling list |
| preparation duration | an **input**; spacing derived from it | an **output**; `n_composite × 1/f_SL` |
| inter-pulse gaps | `(prep_time − 2·rfc_dur) / 4`, four of them | minimal between pulses, longer at the ends |
| crusher | 3-axis, each derated to 1/√3 | identical helper, identical derating |

They are **not** two witnesses — same corpus, same author, shared crusher helper and shared
composite-pulse concept. They are two *designs*, which is evidence about the abstraction rather
than corroboration of either one's numbers.

## 2. The contract, stated from the numbers

`t_inter = (prep_time − 2·rfc_dur)/4`, used four times, with a composite pulse after the first and
the third. So the pulses sit at one quarter and three quarters of the preparation period, the gaps
either side of each are equal, and the two spin echoes coincide with the pulse centres. Summing
the block gives `tip_down + prep_time + tip_up`, which fixes the semantics precisely:

> **`prep_time` is the interval from the end of the tip-down pulse to the start of the tip-up
> pulse** — edge to edge, not centre to centre.

The composite refocusing pulse is constant amplitude with phase −y for the first quarter, +x for
the middle half, −y for the last quarter: the classic **90(−y) 180(+x) 90(−y)**. The two composites
carry opposite phase, `ref_phase ± π/2`, which is the alternation that cancels B1 error to first
order.

Four things a user could get wrong, each producing a legal sequence and a plausible image: the
tip-up phase sense, the interval symmetry (get it wrong and you have measured T2\* while calling it
T2), the refocusing phase alternation, and putting the crusher before the tip-up instead of after.

## 3. Ownership — settled

The refocusing train is **intrinsic to the leaf**. Both designs keep it entirely inside the
preparation and never expose it; its count and spacing follow from the preparation time and the
leaf's own pulse durations; nothing outside influences either. v0's ownership test answers in the
affirmative, and nothing in the evidence argues for a kernel or for a caller-assembled train.

**`Refocusing` is not reusable here, and that was measured rather than assumed.** It refuses
`crush_cycles_per_voxel=0` — "must be positive" — so it cannot emit a crusher-free pulse, and a
crusher inside a T2 preparation destroys the magnetisation being stored. It is also built to be
slice-selective with a balanced crusher pair, which is the opposite job. This is v0's corollary
exactly: *reuse of semantics does not oblige reuse of an emitted event.*

So the shape is a self-contained preparation leaf owning its tip pair, its train, its timing solve
and its crusher — the same shape as `IRPrep` and `SaturationPrep`, in a folder whose contract
already names the third seat.

## 4. What blocked it, and it is not the boundary

**The parameterisation is genuinely open.** One design takes the duration and derives the train;
the other takes the train and derives the duration. Rule A's mode table cannot be written honestly
from that, and picking one now would be picking the one that is easier to code.

**And the more serious one: both witnesses tip with an adiabatic 90, and SeqCraft cannot build
one.** `pypulseq.make_adiabatic_pulse` is documented "Make an adiabatic **inversion** pulse" and
offers only `hypsec` and `wurst`; those invert and cannot perform a 90° rotation. BIR-4 and AHP are
a different family, and OpenMRF generates BIR-4 through SigPy.

My first guess — that SeqCraft has no adiabatic support at all — was wrong, and checking it
mattered: `IRPrep` already routes `hypsec` and `wurst` to `make_adiabatic_pulse`. What is absent is
specifically the 90° adiabatic family.

The consequence is a claim-scope problem rather than a capability one. A `T2Prep` shipped today
would use a non-adiabatic tip pair — a variant **no witness implements**, and one that is *not* a
metamorphic transform of a witnessed variant either. That is the difference from the Spiral
precedent, where `in` really is `out` reversed and analytic plus metamorphic evidence carried the
unwitnessed modes. Here the unwitnessed thing is the pulse family, and B1 robustness is the whole
reason the witnesses chose adiabatic pulses.

**Cheapest unblock:** a permissively licensed non-adiabatic T2 preparation, which would witness
the buildable variant directly instead of by analogy. That is a registry question, not this scan's
to answer.

## 5. What this says about the workflow

**v0 handled this well, and that is the result the candidate was chosen to produce.**

- The ownership test settled the boundary without a tiebreaker, and settled the `Refocusing`
  question by making a claim that turned out to be measurable — the refusal is in the code.
- `MODE_CONTRACT_INCOMPLETE` was already in the reason-code vocabulary and fits exactly. I did not
  have to bend a label or invent one.
- Rule A's "no extraction before the mode table" is what stopped this, and it stopped it for the
  right reason: the evidence contains two parameterisations and does not choose between them.
- The evidence-tier rule — *reference coverage sets the evidence tier, not the abstraction
  boundary* — did real work. The boundary is settled on one witness; what one witness cannot do is
  license a claim about an unwitnessed pulse family.

One small thing v0 did not have a place for: the distinction between **one independent witness**
and **two designs by one author**. Both are `evidence[]` entries with an `independence` note, and
that is adequate — but `reference_situation.designs_examined` is a field I added by hand because
the schema has no way to say "the variation I am reasoning from is intra-corpus". Recorded as a
Phase E observation, not a schema change.

## 6. Recommendation

Do not implement yet. Two questions need a human, and neither is answered by more scanning:

1. Is the leaf's parameter the preparation time or the refocusing count?
2. Should a non-adiabatic `T2Prep` ship at all, given that no witness implements one — or should
   this wait for a witness of the buildable variant, or for a 90° adiabatic family?

If the answer to (2) is "ship the non-adiabatic variant with an explicit claim scope", this becomes
`APPROVED_FOR_IMPLEMENTATION` immediately — the boundary, the ownership and the acceptance claim
are all written and only the mode table is missing. If the answer is "wait for a witness", the
record stands as it is with a trigger.

**Either answer is cheap from here.** That is what a clean conventional case looks like when the
evidence is thin rather than the physics.
