# Compiler constraint and support matrix

> **Historical.** This records the state at the time it was written.  The compiler was
> subsequently changed by [ADR-004](../adr/004-compile-returns-a-sequence.md): `compile` returns a
> bare `pypulseq.Sequence`, every legality failure raises, and `CompiledSequence`, `Report` and
> `Issue` no longer exist.  For the current shape see [`../compiler.md`](../compiler.md) and
> [`../architecture.md`](../architecture.md).

This matrix records the behavior accepted at the Phase 0 baseline. “Frozen” means later no-change
phases must preserve the behavior; it does not imply that every behavior is ideal.

## Event and extension support

| Event type | Placement/reservation | Boundary behavior | Emission/ownership | Status |
|---|---|---|---|---|
| `delay` | Occupies `[node_t, node_t + delay]` | Extends total duration; empty intervals become delay blocks | Represented by block duration rather than copied as a source event | Supported, frozen |
| `rf` | Own delay plus shape; reservation includes dead time and ringdown | Indivisible and exclusive | At most one RF per emitted block; delay is snapped to RF raster | Supported, frozen |
| `adc` | Own delay plus sample window; reservation includes trailing dead time | Indivisible and exclusive | At most one ADC per emitted block; delay is snapped to ADC raster | Supported, frozen |
| `trap` | Gradient waveform span on its axis | Cuttable when a required boundary demands it | Pass through, split, or superpose on the same axis | Supported, frozen |
| `grad` | Arbitrary/extended gradient span on its axis | Cuttable with exact PWL handling where representable | Preserves raster-centre arbitrary form; mixed-knot sums may be resampled with warning | Supported, frozen |
| `labelset` | Zero-duration point | Retimed to the first eligible ADC without crossing a barrier | Owned by the target ADC's block; orphan is emitted and warned | Supported, frozen |
| `labelinc` | Zero-duration point | Same targeting as `labelset` | Same-key increments commute; order-dependent set/inc groups fail | Supported, frozen |
| `trigger` | Pulse duration is indivisible | May share a block and overlap ADC/RF; cannot be split | PyPulseq trigger extension in containing block | Supported, frozen |
| `output` | Same reservation model as trigger | Same as trigger | PyPulseq digital output extension | Supported, frozen |
| `seqcraft_barrier` | Zero-duration pseudo-event | Forces a boundary unless it would split an indivisible event | Not emitted | Supported, frozen |
| `rot3D` | Recognized before scheduling | None | Rejected with build-time rotation guidance | Unsupported, explicit error |
| `soft_delay` | Recognized before scheduling | None | Rejected; plain delay is the current alternative | Unsupported, explicit error |
| `rf_shim` | Recognized before scheduling | None | Rejected | Unsupported, explicit error |
| Unknown `.type` | No silent fallback | None | Rejected with the handled-type list | Unsupported/version skew |

## Scheduling constraints

| ID | Constraint | Enforcement | Observable result |
|---|---|---|---|
| P-01 | Absolute leaf time is the sum of ancestor and local starts. | `flatten`, `placement.place_events` | Nested offsets and provenance paths are stable. |
| P-02 | A negative physical reservation start is illegal. | `compile_sequence` precheck | `CompileError` names the event, source, and likely TE cause. |
| P-03 | Gradient starts must lie on the gradient raster. | `compile_sequence` precheck | `CompileError` gives neighboring legal raster times. |
| P-04 | Same-time tree order is insertion order. | `LogicBlock.nodes`, `flatten` | Characterization test freezes the traversal input order. |
| P-05 | Compilation does not mutate source node times, identity, or numeric event content. | Characterization test | Reusing a tree produces the same source state. |
| L-01 | Two RF reservations, two ADC reservations, or RF/ADC reservations may not overlap. | `_check_exclusive` | Compilation fails before boundary selection. |
| L-02 | A boundary may not fall strictly inside RF, ADC, trigger, or output reservations. | `_Spans`, `_boundaries` | Indivisible hardware actions remain whole. |
| L-03 | Consecutive exclusive reservations have at least one boundary between them. | `_boundaries` | No emitted block contains multiple RF or ADC events. |
| L-04 | Gradient spans prefer not to be cut, but may be cut for a mandatory boundary or barrier. | `_boundaries`, `_axis_gradient` | Waveform is reconstructed across block seams. |
| L-05 | Different-axis gradient overlap is normal and silent. | Per-axis assembly | No merge warning. |
| L-06 | Same-axis gradients are summed and reported. | `_axis_gradient`, `_superpose` | One `grad_merge` warning per affected emitted block. |
| L-07 | A trap plus raster-centre arbitrary waveform may require raster resampling. | `_resampled` | `grad_resample` warning includes a measured maximum error. |
| L-08 | A barrier may not split an indivisible reservation. | `_barrier_conflict` | Clear `CompileError` names the barrier and blocking event. |
| L-09 | Point-only or empty trees do not create a zero-duration sequence. | Total-duration check | `CompileError: nothing to compile`. |
| E-01 | Event delay is relative to its emitted block and lies on the event's raster. | `_in_block_delay` | RF, ADC, and gradient delays are reproducible. |
| E-02 | Emitted block duration is a block-raster multiple and covers all events. | `_required_duration` | Short blocks fail before `Sequence.add_block`. |
| E-03 | Per-axis amplitude/slew is checked after merge. | `_limit_issues` | Violations are report errors, not compile exceptions. |
| E-04 | Simultaneous-axis vector norm is informational for hardware rotation risk. | `_limit_issues` | Norm violations are warnings. |
| E-05 | Origin is the longest common source path of a block's contributors. | `_common_path` | One stable provenance path per block. |
| V-01 | Emitted duration equals the tree horizon. | `CompiledSequence._verify` | A mismatch is a report error. |
| V-02 | Per-axis m0 and absolute-time m1 are preserved. | `_verify`, `moments` | Lost or shifted gradient pieces are report errors. |
| V-03 | Compiled labels match the tree-derived ADC states. | `_address_issues` | Mis-retimed labels are report errors. |
| V-04 | Imaging ADC addresses are unique unless marked noise/reference/navigation. | `_label_issues` | Duplicate addresses are report errors. |

## Boundary ownership

| Boundary source | Strength | May split gradient | May split RF/ADC/trigger/output | Current tie behavior |
|---|---|---:|---:|---|
| Sequence start/end | Hard | Yes | No interior split | Exact raster endpoints |
| Explicit barrier | Hard | Yes | No | Barrier time snapped to nearest block raster |
| Gap between exclusive reservations | Hard | Yes | No | Existing mark, then gradient edge, then raster midpoint fallback |
| Gradient/reservation edge | Natural candidate | No | No | Accepted only when outside all protected spans |
| Maximum block-duration subdivision | Hard storage constraint | Yes | Must remain outside protected spans by construction | Even raster subdivision |

The current compiler does not retain the candidate set, rejected reasons, or a boundary score. Phase
1–3 may expose those as private contracts, but they must reproduce the Phase 0 outputs before any
boundary-policy change is considered.
