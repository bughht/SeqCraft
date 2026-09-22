# Repetition augmentation — a direction, not a design

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-20, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../README.md) for which one to edit.

**Non-binding.** No API is proposed here and no implementation is authorised by it. It exists so
that a future contributor who needs flow compensation, diffusion encoding, an FID navigator or
velocity encoding finds the argument already made, rather than reaching for the first mechanism
that works.

## The gap, stated plainly

| | today |
|---|---|
| arbitrary overlap in `LogicBlock` | yes |
| semantic timing queries on modules (`time_to_center`, `time_to_echo`) | yes |
| manual composition around a kernel | yes |
| a generic repetition extension point | **no** |
| a standard physical context an augmentation could read | **no** |
| joint solving of base and augmentation requirements | **no** |

**None of that is currently a defect.** Every shipped module composes by overlap and semantic
timing, and nothing has yet needed to change the physical contract of a repetition that already
exists. The gap becomes relevant the first time something does.

## Two kinds of composition, which are not alike

**Independent insertion.** The added thing needs a semantic *window* and nothing else. It does not
change the correctness conditions of the gradients, the RF or the echo that were already there —
a navigator in dead time, a preparation before the train. Slot- or window-style composition is
sufficient, and `IRPrep` plus `GRE2D.time_to_center_line` is already an instance of it.

**Coupled augmentation.** The base repetition's requirements and the augmentation's requirements
are *the same solve*. Flow compensation, diffusion encoding and velocity encoding are all of this
kind: what the readout needs and what the augmentation needs are constraints on one set of
gradients, and satisfying them separately satisfies neither.

## The direction

> Future repetition augmentation should preferentially compose **physical requirements before
> waveform realization**, rather than patching an already-built repetition.

A flow-compensation module should not receive a built GRE repetition, measure its first moment,
and append a corrective lobe. That ordering makes the augmentation depend on decisions the kernel
has already made — which is the same failure the project rejected when it declined to nest
`GRE2DTR` inside `GRE3DTR`, and when it ruled that a mode may not be "another mode minus one
event".

```text
base requirements  +  augmentation requirements
                   ↓
            one physical solve
                   ↓
               LogicBlock
```

> When an augmentation needs information from the base sequence, expose **semantic physical state
> or requirements**, not raw implementation events.

Conceptually a flow compensation would want the moment state at a semantic echo, the M0/M1 it must
satisfy, the windows available either side of the echo, the semantic RF and echo times, and the
hardware limits — **not** a handle on a private gradient object. Diffusion would additionally want
the desired b-value or b-tensor, the refocusing timing, the imaging gradients' contribution, the
cross terms and the available windows.

This project already has the precedent for the shape of that: `Excitation.rephaser_area_per_m`
states a *requirement* while `build(rephase=False)` declines to realise it, which is what let
`GRE3DTR` fold the slab rephasing into its own z solve without either module reaching into the
other.

## Explicitly not proposed

```text
GRE3DTR(extra_modules=[...])        a generic plugin bus
public access to private gradients  semantic anchors inside LogicBlock
GradientMomentRequirement           DiffusionRequirement
a constraint-solver layer
```

Each may become right. None is justified now, and the rule that decides is unchanged:

> **One consumer: do not generalize. Two or more real consumers: inspect the shared physical
> contract.**

Revisit when two of these are real, not one: `FlowCompensation`, `DiffusionEncoding`,
`FIDNavigator`, velocity/motion encoding. The first to arrive should be built **candidate-local**,
and the second is what turns this note into a design.
