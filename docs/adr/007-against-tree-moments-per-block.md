# ADR-007: The against-tree moments are compared per compiled block, in two stages

- Status: Proposed
- Date: 2026-10-26
- Amends: the m0 and m1 invariants of the always-on verification reviewed in
  [`refactor/phase6_verification_time_audit.md`](../refactor/phase6_verification_time_audit.md).
  The duration and label-address invariants are unchanged, and verification stays always on.

## Context

`verify_against_tree` compared each axis's m0 and m1 over the whole sequence, with m1 referenced
to `t = 0`, and read the sequence PyPulseq had registered. Registration is not exact. A
compressible shape is stored as differences rounded to `1e-7` of its peak, and an event whose
parameters print like an earlier one's under `%.6g` is stored as that one. A registration error
`δm0` in a shape placed at `t` moves the absolute m1 by `t·δm0`, so a legal waveform that compiled
at 0 s was refused at 1 s. A tolerance wide enough to absorb that at 10 s, 3e-4 s/m, would hide a
4.5 1/m blip played one raster late, 4.5e-5 s/m. One tolerance also had to cover two unrelated
error sources: SeqCraft's compile and PyPulseq's storage.

## Decision

Compare m0 and m1 **per compiled block and axis, each about the block's own start**, in two
stages.

1. **Compile.** The tree's placed events, restricted exactly to each block, against the ready
   blocks that emission hands PyPulseq, measured before registration. The tolerance is SeqCraft's
   own arithmetic, each term a stated bound: float summation, the float of superposing pieces, the
   endpoint snap, and a tick where a knot's time can land one from where the other side puts it.
2. **Registration.** Those ready blocks against what PyPulseq stored, read back block by block.
   The tolerance is PyPulseq's storage precision: compression rounding, bounded by its first-sample
   edge case at 1.5 quanta of the peak, and the six-digit deduplication key. That allowance belongs
   to this stage alone.

The whole-sequence m0 and m1 refusals are retired. For any partition, `M0 = Σ m0_k` and
`M1 = Σ (m1_k + a_k·m0_k)`. A whole-sequence comparison that tolerated every block's error could
therefore never fire, and a narrower one would bring the translated term `a_k·δm0_k` back.

A block-axis that legalization resampled onto the raster is reported as data
(`LegalizationResult.resampled`), never inferred from warning text. The compile stage compares it
with the tree resampled onto the gradient raster from the block's start. Merging a resampled
stretch's inner seams was considered and rejected: it is exact only where both ends of the stretch
are flat, since a sloped end keeps a boundary term of about `r²·g′/8`, and it coarsens the check
inside the stretch.

## Consequences

- A correct waveform compares the same wherever it sits. A refusal names its stage, block, start
  time, difference and tolerance.
- Area moved across a block boundary is refused, which the whole-sequence m0 could not see.
- PyPulseq's six-digit deduplication is accepted as its precision. An event that shares a key with
  an earlier one still plays the earlier one's values.
- m0 and m1 remain two linear functionals per block, not waveform fidelity.
- `superpose` interpolates in ticks from the block's start, in a commit of its own with the
  tolerance term it makes unnecessary.
- The endpoint snap, which zeroes a superposed block's edge value below `1e-6·max_grad`, stays and
  is modelled. Narrowing it to float residue was tried and reverted: PyPulseq's timing check refuses
  an end above 1e-9 Hz/m that the next block does not continue, and the snap is what keeps a summed
  end that cancels to nearly zero from being refused.
