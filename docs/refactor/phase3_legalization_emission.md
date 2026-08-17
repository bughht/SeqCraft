# Phase 3 legalization and emission boundary

- Status: Complete; local and hosted gates passed
- Date: 2026-08-16
- Branch: `refactor/compiler-phases`
- Baseline: `45a404a`
- Policy change: None

## Scope delivered

- Added `LegalizationResult`, the frozen stage result containing the complete ready-block tuple and
  immutable transformation notes.
- Moved singleton scheduling, gradient interval assembly, ADC seam protection, limit and duration
  checks, provenance ownership and ready-block construction into `legalize_blocks`.
- Verify the complete ready-block sequence once before it leaves legalization.
- Reduced `emit_blocks` to ready-block iteration, explicit block duration, `Sequence.add_block` and
  contextual mapping of PyPulseq failures.
- Drain a private ready-block queue during emission so the immutable stage result does not increase
  peak memory while PyPulseq registers events.
- Added direct synthetic tests for both stage boundaries and a source-level dependency test that
  prevents emission from importing policy-bearing stages.
- Added legalization and emission to the strict mypy CI gate and updated the executable API
  reference and architecture documentation.

Phase 4 boundary-policy redesign remains deferred. The old Phase 5 mechanical-emission extraction
is complete as part of this phase, so there is no intermediate adapter or second compile path.

## Resulting contract

```text
PlacedEvent + fixed boundaries + label targets
        -> legalize_blocks
LegalizationResult(blocks, notes)
        -> emit_blocks
pypulseq.Sequence
```

Legalization is the last policy-bearing stage. Emission cannot see placed events, boundaries, label
targets, scanner limits or transformation-note storage. The public contract remains
`sc.compile(root, opts, ...) -> pypulseq.Sequence`.

## No-change evidence

A five-iteration before/after capture used the same Python environment and the `45a404a` source
snapshot. Both GRE and spin-echo recipes matched exactly on stable and observed fields: duration,
block counts and durations, event counts and content digest, moments, warnings, input and placed
event counts, gradient split counts, provenance digest and written file size.

| Recipe | Median wall-time ratio | Peak Python-memory ratio |
|---|---:|---:|
| GRE 2D | 0.974 | 1.000 |
| Spin echo 2D | 0.987 | 1.001 |

The first complete-tuple implementation retained ready events throughout emission and increased
peak Python memory by 69–77%. It was rejected. Draining a private queue after whole-result
verification restored the baseline without weakening the immutable stage contract.

## Local verification

| Gate | Result |
|---|---|
| Focused compiler suite | 136 passed |
| Non-heavy suite with coverage | 361 passed, 6 optional-dependency skips; 89% coverage |
| Source doctests | 44 passed |
| Ruff | passed |
| Strict mypy | 9 source modules passed, including legalization and emission |
| Executable API reference | 51 Python blocks; examples and full `__all__` index passed |
| Getting-started notebook smoke | passed |
| Same-runner semantic/structural baseline | exact match |
| Same-runner performance guardrail | passed |

The six skips require the optional `pulseq_systems` dependency and are unchanged from the baseline.

## Hosted verification

[GitHub Actions run 31989728932](https://github.com/bughht/SeqCraft/actions/runs/31989728932)
passed all seven jobs on PR #10:

- tests on Ubuntu and Windows with Python 3.11 and 3.12;
- Ruff linting;
- strict mypy type checking; and
- executable examples.

The hosted matrix and all local gates therefore satisfy the Phase 3 exit criteria.
