# Repetition-level joint design — spike findings, before framework code

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-24, so that the pull request carries the direction behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../../README.md) for which one to edit.

**Status:** findings for review. No public API proposed.
**Date:** 2026-09-24
**Reads:** [`2026-09-24_repetition_level_joint_waveform_design.md`](2026-09-24_repetition_level_joint_waveform_design.md)

Measured on `spike/repetition-joint-design`, branched from PR #39's head. Every number comes from
building real events and integrating the emitted waveform.

---

## 1. The model

```text
MomentTarget(axis, order, origin, endpoint, value)
```

Both instants are required. Under a shift of origin `m1' = m1 - dt * m0`, so a first moment is
only a number once the origin is named. **`order = 0` does not read the origin**, which is why the
omission survived the first pass unnoticed; that is not a reason to leave it out of the semantics.

For the single-path GRE cases here:

```text
origin    the excitation's effective rotation instant, supplied by the repetition
endpoint  the echo -- or, in a train, the checkpoint being solved for
```

The origin is **supplied by the repetition**, not read off a particular class. For today's
`Excitation` it is `time_to_center()`, which is `pp.calc_rf_center` rather than the geometric
midpoint. A future RF module whose effective instant is defined differently should expose that
instant; nothing in the designer changes.

### How much the origin matters

At `ky = +16` over a 2.5 ms excitation-to-echo interval, `dt * m0 = 0.182 s/m` — **larger than an
entire VENC target of 0.167 s/m at 1.5 m/s**. The first pass of this spike measured `m1` about the
echo with `m0 != 0` and called it flow compensation; that quantity had no physical name and those
rows are withdrawn. The argument is worth keeping, the table is not.

The readout-axis case is the exception that hid it: `m0 = 0` at the echo there, so PR #39's `m1`
is origin-free. Measured about the echo, the block start and −5 ms it is `~1e-16` compensated and
`1.03e-01` uncompensated, all three.

## 2. The realisation: a basis, and a linear solve

```text
fixed physical contribution
+ adjustable waveform family in a window
        -> moment response matrix M
        -> solve  M a = target - fixed
```

Two fixed-duration lobes are two columns; the coefficients are their areas. Measured, over ten
`(ky, polarity)` states in one shared window of 2 x 630 us on `y`, origin at the excitation
centre:

```text
m0            exactly the encode area, every line
m1 difference 0.333333 against 0.333333 wanted, error <= 4.4e-16
m1 common     <= 1.4e-16
```

The limiting state is `(ky = -32, polarity = +1)` — neither an extreme of `ky` nor a fixed
polarity, which is `GRE3DTR`'s own "the limiting partition is a result of the coupling, never a
property of the index" reappearing in a larger state space.

**Representation and shape family are different things.** The adjacent two-lobe basis and the
decoupled base-plus-zero-area-bipolar basis span the **same** `(m0, m1)` set in the same window and
report identical utilisation -- they are two coordinate systems for one two-dimensional space, so
the decoupled one is better conditioned and more readable, not more capable. The first family that
adds genuine waveform freedom is the **raster-aligned split search**, which is 8 to 13 per cent
shorter at the same target (650 us against 710 at nominal limits, 910 against 1040 at the weakest
system tried).

## 3. Composition: augmentations resolve state, the designer sees absolutes

The C3 record makes velocity encoding's public quantity the **change** between states, with the
`+-delta/2` split following from the standalone toggles being exact waveform negations. So the two
augmentations constrain different components:

```text
VelocityEncode   m1(+) - m1(-)       = delta_m1
FlowComp         (m1(+) + m1(-)) / 2 = 0
```

and composing them **derives** `m1(state) = state * delta_m1 / 2` rather than assuming it.

Because they write different components they cannot contradict, so no precedence rule is needed;
a second claim on the same `(axis, order, component)` is a refusal. That resolution happens in the
augmentation layer. **The waveform solve receives one absolute target for one concrete state** and
never learns whether it came from flow compensation, velocity encoding, phase encoding or a future
acquisition policy.

## 4. Multi-echo is a recurrence, not one global matrix

`wave-gre-flow-comp` propagates moment state between echoes instead of solving the train at once:

```text
% If echo c-1 is already FC'ed, the maintained PE zeroth moment
% accumulates an additional M1_PE = -M0_PE*ESP at echo c.
% The zero-M0 inter-echo y FC lobe therefore targets +M0_PE*ESP.
```

That is `m1' = m1 - dt * m0` again — the same origin-shift law, applied as propagation. So:

```text
checkpoint n satisfied
    -> propagate the maintained state to checkpoint n+1   (m0 unchanged, m1 -= esp * m0)
    -> add the fixed contributions in the interval, about checkpoint n+1
    -> solve the adjustable window in that interval
    -> checkpoint n+1 satisfied
```

Each interval is its own small square solve. The model should therefore carry **checkpoints, the
windows between them, and a propagated moment state** — which is simpler than a train-wide matrix,
not harder.

## 5. The kernels already own the window decision

```text
GRE2DTR   winder_s = ceil_raster(max(ro.prephaser_duration_s,
                                     pe.min_duration_s,
                                     exc.rephaser_duration_s))
```

One shared window, sized from the limiting demand, across three axes — already. `GRE3DTR` does the
same across partitions, and its single reference lobe scaled per partition is the **1x1 case** of
the basis above, with `scale_grad` as the coefficient. Neither needs a new mechanism; the
requirement they size against needs an order and an origin.

| | GRE2DTR | GRE3DTR |
|---|---|---|
| window today | 520 us | 320 us |
| origin | 1.600 ms | 0.200 ms |
| echo | 6.218 ms | 3.198 ms |
| window with `m1` | 1440 us on `y` (venc 1.5 + flow comp) | 580 us on `z` (flow comp) |
| TE cost | +0.920 ms | +0.260 ms |

## 6. Slice selection needs no new leaf API

With the origin at the excitation's effective instant, the excitation's own contribution is just
another fixed moment. Measured on `GRE2DTR` from that instant to the echo:

```text
z   m0  +0.0000 1/m        the slice rephaser already nulls it
z   m1  -0.425956 s/m      the fine scan's `ms`, and 2.5x a VENC-1.5 target
```

The **kernel** computes it from events it already places, the same way it obtains the readout's
fixed moments. The fine scan's open question about what a leaf must expose is answered by needing
nothing new: no leaf publishes a moment.

## 7. What `wave-gre-flow-comp` classifies into

Registered as `stress-case`, MIT, at `0e1ec51`, admitted file by file rather than by admitting
HarmonizedMRI. Nothing is adapted from it.

| its logic | where it belongs |
|---|---|
| `makeM0PreservingM1CorrectedLobe` — base lobe + zero-area bipolar, solved for M1 | generic: the basis and the linear solve |
| `postGapToRef` / `tRef` — M1 referenced to a named instant | generic: origin/endpoint on the target |
| `designMinDurationM0M1LobeRefStart` — grow by `4*dt` until it fits | generic: the window search (bisection here) |
| `calcInterEchoCosExternalMoments` — neighbouring tails over an interval | generic: fixed contribution in a window |
| the `M1 = -M0*ESP` recurrence | generic: moment-state propagation |
| the `7point` shape-parameter search | realisation family, beneath the physical problem |
| cosine/sine wave gradients | wave-encoding-specific: a fixed contribution to this kernel |
| LIN/PAR tables, ESP, `t_adc_center` | GRE repetition timing and state |
| per-protocol flags and hard-coded targets | implementation convenience, and what a shared boundary removes |

**It does not contradict the boundary.** Everything it does that is generic is `(fixed moments +
origin/endpoint + adjustable basis/window + target)`; everything left over is either
wave-specific physics or timing policy that belongs to the kernel.

## 8. The cliff question is **not** closed, and the optimizer question is not either

An earlier draft of this note concluded "the window grows smoothly, there is no artificial cliff,
GrOpt is not needed". **That is withdrawn.** What the sweeps establish is narrower.

With the schedule free (AUTO), the minimum winder grows smoothly as the system is derated -- a
factor of 1.12 to 1.24 per step from 40 mT/m / 150 T/m/s down to 8/40, with no jump. Below that
the kernel refuses, and the refusal is **not** the joint design: the plain kernel refuses
identically at 6/30 with the same pypulseq amplitude violation, so that boundary belongs to the
base sequence.

With the schedule **pinned** at the nominal minimum TE and then derated, the same target degrades
smoothly too, and the diagnostic says by how much:

```text
 G      S        u = max(G_req/G_lim, S_req/S_lim)     fits
40    150                 0.32                         yes
20    100                 0.55                         yes
15     80                 0.72                         yes
10     60                 1.05                         no     <- misses by 5 %
 8     40                 1.39                         no
 2      8                 6.25                         no
```

So at 10/60 the fixed schedule misses by five per cent while AUTO finds a feasible design at a
24 per cent longer winder. That is the two-level design earning its keep, measured.

**What this does not establish.** It does not reproduce `wave-gre-flow-comp`'s pathology, which
needs fixed TE *and* fixed echo spacing, a multi-echo train, and the wave gradients occupying the
inter-echo windows -- none of which these single-echo GRE probes have. Until that is reproduced:

```text
optimizer required      NOT established
optimizer unnecessary   ALSO NOT established
```

No GrOpt comparison has been run. The right place for one is beside a reproduced cliff, not here.

## 9. The prototype

`src/seqcraft/modules/_joint.py` -- internal, provisional, 15 tests in
`tests/modules/test_joint_moment_design.py`. Neither kernel is modified; each supplies facts it
already computes, which is what makes this a test of the boundary rather than a rewiring.

```text
resolve()   claims -> ONE absolute target per repetition state
solve()     one state's absolute target -> areas, and never learns where it came from
```

`GRE2DTR` on `y` and `GRE3DTR` on `z` build the same `JointProblem` from their own facts and get
the same machinery. The 3D kernel's partition term enters as part of *its* base requirement, not
as anything an augmentation knows about. Flow compensation alone, velocity encoding alone, and
both together differ only in the claim list -- asserted, so a third branch would have to show up
here to pass.

Two things the prototype found rather than assumed:

* a difference claim whose partner state is not acquired raised `KeyError` from inside the
  resolver; it refuses now, because a phase difference needs two acquisitions to subtract;
* the state a difference relates is the encoding index, not the whole composite key, so `resolve`
  is called per encode index -- which is also how a kernel would loop.

`VelocityEncode.delta_m1_for(venc)` is new and is the formula without a waveform, so a joint
design reads the factor of two from the same place the standalone bipolar does. No second public
velocity abstraction, and no lazy-realisation refactor was needed.

## 10. Open

- The full TE/TR sweep needs the kernels actually rewired: today the prototype measures that the
  window grows and pins the kernels' existing `te_s`/`tr_s` refusals, but it cannot yet make a
  kernel build with the longer window and watch its own minimum move.
- Pathway-aware moments, deliberately deferred: pathway-aware `m0` first, `m1` later, neither now.
