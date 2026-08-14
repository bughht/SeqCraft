# Compiler current state

This document freezes the compiler's pre-refactor shape at commit
`afe576bafe17e9cce3edfc768b343ad958994560`. It describes what exists, not the intended Phase 1
design. The corresponding remote baseline is [GitHub Actions run #7](https://github.com/bughht/SeqCraft/actions/runs/31848625433), whose seven Linux/Windows, Python 3.11/3.12, lint, type, and example jobs passed.

## Public surface

The public compilation entry points are aliases of one function:

- `seqcraft.compile(root, system, *, geometry=None, name='', regime='default', definitions=None)`;
- `seqcraft.compile_sequence(...)`;
- `seqcraft.core.compile_sequence(...)`.

The public results are `CompiledSequence` and `WriteResult`. `CompiledSequence.seq` deliberately
exposes the underlying `pypulseq.Sequence`; `report`, `origins`, `definitions`, `duration_s`,
`n_blocks`, `check()`, `moments()`, `kspace()`, `pns()`, and `write()` are compatibility surface.

`compiler.py` already contains a private immutable `_Placed` dataclass. Phase 1 must adapt or rename
this real contract rather than introduce a second placement record. There is no ready-block record:
boundary selection, gradient legalization, PyPulseq event construction, block emission, provenance,
and limit reporting are interleaved inside `compile_sequence()`.

## Current control flow

```text
LogicBlock
  -> _place(root, opts) -> list[_Placed]
  -> negative/off-raster checks
  -> _check_exclusive
  -> _boundaries
  -> label target and order resolution
  -> per-boundary assembly loop
       -> _axis_gradient
            -> pass through | _superpose
            -> _as_arbitrary | _resampled
       -> _in_block_delay
       -> _required_duration
       -> _limit_issues
       -> pypulseq.Sequence.add_block
  -> CompiledSequence
       -> _verify(duration, m0, m1, labels)
       -> check(timing, labels, event sizes)
```

The direct private call graph is:

| Caller | Compiler-local callees |
|---|---|
| `_place` | `_Placed`, `_intrinsic_duration`, `_unsupported` |
| `_barrier_conflict`, `_gap_blocked` | `_covering` |
| `_boundaries` | `_Spans`, `_barrier_conflict`, `_gap_boundary`, `_gap_blocked` |
| `_axis_gradient` | `_in_block_delay`, `_superpose`, `_as_arbitrary`, `_resampled` |
| `compile_sequence` | `_place`, `_check_exclusive`, `_boundaries`, label helpers, `_adc_conflict`, `_axis_gradient`, `_in_block_delay`, `_required_duration`, `_limit_issues`, `_common_path`, `CompiledSequence` |
| `CompiledSequence._verify` | `moments`, `_address_issues` |
| `CompiledSequence.check` | PyPulseq timing check, `_label_issues`, `_event_size_issues` |

## Responsibility inventory

| Responsibility | Current symbols | Phase boundary implication |
|---|---|---|
| Placement | `_Placed`, `_intrinsic_duration`, `_unsupported`, `_place`, `_check_exclusive` | `_Placed` is real but private and compiler-coupled. |
| Label legalization | `_label_targets`, `_label_order_conflict`, `_expected_addresses`, `_orphan_label_issues` | Label ownership is resolved before block assembly but has no explicit stage result. |
| Boundary selection | `_Spans`, `_covering`, `_barrier_conflict`, `_gap_boundary`, `_gap_blocked`, `_boundaries` | Returns only times; reasons and candidate decisions are not retained. |
| Gradient legalization | `_superpose`, `_as_arbitrary`, `_axis_gradient`, `_resampled` | Operates while emitting a block, so there is no independently inspectable ready waveform. |
| Emission | `_in_block_delay`, `_required_duration`, `_adc_conflict`, the `compile_sequence` loop | Scheduling and PyPulseq construction share one loop. |
| Reporting | `_limit_issues`, `Report`, error builders embedded in the functions above | Stable messages are part of the characterized behavior. |
| Verification | `CompiledSequence._verify`, `_address_issues`, `moments`, `check`, `_event_size_issues`, `_label_issues` | Verification runs after emission and mutates only the result's cached report. |
| Output | `kspace`, `pns`, `write`, `_jsonable` | Not part of scheduling extraction; public compatibility still applies. |

## Test and sequence inventory

| Asset | Frozen behavior |
|---|---|
| `tests/logic/test_logic_block.py` | Recursive absolute offsets, insertion order, duration, barriers, and minimal tree model. |
| `tests/compiler/test_scheduling.py` | Axis overlap rules, RF/ADC exclusivity, barriers, rasters, limits, provenance, labels, input immutability, and zero-duration point-only trees. |
| `tests/compiler/test_boundaries.py` | Mandatory gaps, maximum block duration, and non-quadratic EPI boundary selection. |
| `tests/compiler/test_fidelity.py` + `fidelity.py` | Independent piecewise-linear waveform oracle for split, merge, arbitrary, EPI, spiral, and diffusion cases. |
| `tests/compiler/test_labels.py` | ADC-targeted label semantics, order conflicts, barriers, and orphan labels. |
| `tests/compiler/test_triggers.py` | Trigger/output indivisibility, overlap, barrier conflicts, and multiple extensions per block. |
| `tests/compiler/test_event_types.py` | Handled/unsupported PyPulseq event accounting and failure paths. |
| `tests/compiler/test_invariants.py` | Duration, exact m0/m1, address verification, and fault injection. |
| `tests/integration/test_sequences.py` | GRE, spin echo, spiral DTI, EPI DWI, physics, round trip, determinism, performance budget, and Phase 0 structural snapshots. |
| `examples/` | Getting started, spiral DTI build, and EPI DWI build notebooks in isolated smoke tests. |

The versioned artifact `tests/baselines/compiler_phase0.json` records the four integration recipes.
It includes block-duration, emitted-content and provenance digests; event and issue counts; moments;
placed-event, gradient-split and file-size observations; and local timing/memory samples. The stable
subset is asserted by the integration suite. Performance, file size, and the exact floating-point
content digest are observations tied to the recorded runner, not cross-platform exact assertions.

## Frozen behavior and known bugs

Frozen behavior includes public APIs, supported inputs, waveform/timing results, block structure,
provenance, warning/error kinds and messages, label targeting, and the intentional overlap rules in
the constraint matrix. A later intentional behavior change requires its own decision record and
baseline update.

The following are known bugs or debt, not compatibility promises:

- `NaN` and positive-infinite node starts reach integer timing conversion and leak `ValueError` or
  `OverflowError`; negative infinity happens to enter the normal negative-start `CompileError` path.
- PyPulseq emits repeated shape-restoration warnings while calculating some spiral trajectories.
  Current physics tests pass, but the warning surface needs an upstream/representation audit.
- Several RF fixtures intentionally trigger dead-time adjustment warnings; conflict tests also build
  deliberately over-limit pulses. These are expected test-construction warnings, not compiler output.
- CI constraints pin critical direct dependencies, not a complete transitive lock graph.
- Repository-wide formatting and full-source mypy remain explicit debt outside this refactor.

Phase 0 does not fix these items. Fixes must be isolated from no-change extraction phases and must
state whether they intentionally update the baseline.
