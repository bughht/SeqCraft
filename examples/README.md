# Examples

| | What it covers |
|---|---|
| [`01_getting_started.ipynb`](01_getting_started.ipynb) | Blocks, `Opts` and `compile`; the overlap rules, provenance, checking against physics, writing a file — and `sc.Module` at the end, once there is a reason for one. **Uses no modules**, which is the point: it is the demonstration that the compile path stands alone. |
| [`_parked/`](_parked/) | Two complete DTI acquisitions — spiral and EPI, one diffusion encoding between them. **They do not run against this version.** |
| [`lib/`](lib/) | Simulation and reconstruction helpers. Not part of the package; see [`lib/README.md`](lib/README.md). |

## Why `_parked/` exists

Both DTI folders were built on `seqcraft.modules`, which has been deleted rather than migrated
([ADR-003](../docs/adr/003-scanner-and-module-reform.md)). They are kept unmodified because they are the
**acceptance test for whatever module set is written next**: a library chosen in the abstract
acquires classes nobody needs and misses the ones everybody writes by hand, while one chosen by
making these two scans work again does not.

Read [`_parked/README.md`](_parked/README.md) for what each covers and where their physics went.

## Requirements

Building needs only `seqcraft`. Simulating and reconstructing need `MRzeroCore`, `torch` and
`sigpy`, which are **not** part of the package — see [`lib/README.md`](lib/README.md).

Outputs under any `seq/` are build products and are not tracked by git.
