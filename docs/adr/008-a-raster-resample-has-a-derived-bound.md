# ADR-008: A raster resample has a derived bound, and is refused beyond it

- Status: Proposed (2026-11-05)
- Date: 2026-11-05
- Amends: the resampled-block clause of
  [ADR-007](007-against-tree-moments-per-block.md). ADR-007's per-block, two-stage comparison is
  unchanged, including its comparison of a resampled block with the tree resampled the same way.

## Context

Legalization resamples a block-axis onto the gradient raster when no pulseq gradient can hold its
knots. The usual case is a trapezoid, with knots on raster edges, summed with a gradient sampled at
raster centres. ADR-007 compares the emitted block with the tree resampled the same way. That checks
that the compiler emits its approximation faithfully, but nothing bounds the approximation itself.
Its contract was the `resample` warning, which states only the largest amplitude change, as a
fraction of `max_grad`.

A review measured the approximation on `7f24c83`. The record is in the workspace, in
`docs/plans/module-mining/compiler_findings/resample_bound/`.

- **The amplitude figure does not measure the moments.** The same 1.87 % of `max_grad` went with no
  change in m0, for a raster-rate zigzag, and with a 7 Δk change, for a legal half-raster input.
- **Nothing bounded a half-raster input.** A gradient sampled every half raster
  (`oversampling=True`), summed with a trapezoid, moved k by 7 Δk in one block. That grows with
  its duration.
- **Errors add up across blocks.** Fifty resampled blocks, each within any per-block allowance,
  moved m0 by 1.75 Δk in total.
- **Some resamples were avoidable.** A raster-centre gradient that starts inside its block and is
  then cut is the same event with a delay, but it was resampled.

## Decision

**The bound.** Keeping only the raster-edge values replaces, in each raster interval, a
raster-centre sample's tent by its chord. Over a stretch of samples `c_i` the area changes
telescope:

    Δm0 = (r/8)·(D_end − D_start)
    Δm1 = (r/8)·(t_b·D_end − t_a·D_start) − (r²/8)·(g_b − g_a)     about the block's start

- `D` is the sample step across each end of the stretch: from the virtual sample that `first` or
  `last` implies, or across a cut.
- `t_a` and `t_b` are the stretch's ends, and `g_a` and `g_b` its values there.
- Edge knots are kept exactly, so a trapezoid in the sum adds nothing.

Within the slew limit `|D| ≤ max_slew·r`. A stretch therefore moves the running k by at most
`r²·max_slew/4` at every instant, whatever its duration.

**m1 needs no bound of its own, but the k bound must hold at every instant.** By parts, for a block
of length `L`:

    Δm1 = L·Δk(L) − ∫₀ᴸ Δk(t) dt,   so   |Δm1| ≤ 2·L·sup|Δk|.

The block's m0 alone is not enough. A half-raster input whose edges rise and then fall cancels its
m0 to 0.16 1/m, inside the bound, while k climbs 15.9 1/m and m1 moves by 0.016 s/m. The bound is
therefore placed on `sup|Δk|` over the block. That bounds both m0 and m1.

**Refused beyond it.** `_resampled` measures the resample's running k change exactly, on the union
of both knot sets and where the difference crosses zero. It raises `CompileError` when the change
exceeds `r²·max_slew/4` per non-edge piece in the block. The only tolerance is a relative 1e-9 for
float rounding, because a waveform built exactly at the slew limit reaches the bound.

- **Raster-centre inputs within the slew limit:** never refused. That follows from the
  derivation.
- **Half-raster inputs:** accepted only when measured within the bound, never assumed to be.
- **Inputs built past the slew limit:** refused when their resample shows it.

**Reported in moments.** Each resample's warning entry still gives the amplitude change first. It
now adds the largest k change against its bound, and the block's m0 and m1 change.

An axis resampled in more than one block gets a `resample_total` warning. It gives the net m0, which
is what a balanced design sees, and the summed magnitudes of the k and m1 changes, which
cancellation would otherwise hide.

There is no signed m1 total. Block m1 values share no origin, and about `t = 0` a total would bring
back the absolute-time lever arm that ADR-007 retired.

**Held exactly where possible.** A raster-centre run that starts at a raster edge inside its block,
with zero before it, or ends before the block does, with zero after it, is emitted as an arbitrary
gradient with a delay instead of being resampled.

PyPulseq stores a centre-sampled shape of four samples or fewer without its last half raster, and
ADR-007 refuses that at registration. A padded run that short therefore stays a bounded resample.
The run at the block's start keeps its existing behaviour.

## Alternatives considered

- **Warning-only, as written.** The stated figure does not bound k or m1, and a legal input moved k
  by 7 Δk under it.
- **Comparing the resampled block with the tree as given.** That refuses every resample, which
  ADR-007 already rejected, and it would still need a tolerance. The derived bound is that
  tolerance, applied where the resample happens.
- **A bound on block m0 only.** It misses a resample whose m0 cancels while k strays, and with it
  m1.
- **Emitting mixed edge and centre knots as a half-raster (oversampled) gradient.** It would be
  exact for every input SeqCraft accepts, but it depends on interpreter support for half-raster
  gradient time shapes, which is not established. It stays open.
- **Restoring a whole-sequence m1 check.** Its absolute-time false positives are what ADR-007
  removed.

## Consequences

- Every resample is within `r²·max_slew/4` of k per non-edge piece at every instant, and within
  twice its block's length times that in m1. A resample beyond that is refused, naming the block,
  axis, change and bound.
  - At a 10 µs raster and 150 T/m/s, the bound is 0.160 1/m: 0.035 Δk at a 220 mm FOV.
  - At 200 T/m/s and a 500 mm FOV, it is 0.106 Δk.
- The compile stage's comparison with the resampled tree is unchanged, so emission is still
  checked against the approximation it declared.
- **What changes for raster-centre gradients:**
  - one that starts inside its block and is cut, and one that ends before its block does, are now
    emitted exactly, so their sequences change, for the better;
  - one that `main` resampled, summed with a trapezoid, is emitted as before, with a longer warning.
- **Not bounded:** the accumulation across blocks. It is reported, not refused, because the
  compiler cannot know which repetition a residual matters for.
