# `matlab-frontend/` — the withdrawn MATLAB frontend, parked whole

**Frozen at `d27217d`. Not packaged, not imported, not linted, not tested.** Nothing in `src/`
refers to it and no CI job touches it. It is here so that a working implementation is one `git mv`
away if a MATLAB user ever turns up with the evidence that ADR-005 asked for and did not get.

This entry is a different kind from the rest of `salvage/`. The others are **physics lifted out**
of deleted classes and rewritten as plain functions. This is a **whole subsystem parked intact** —
MATLAB package, Python adapter, tests, fixtures and documentation — because its value is in the
wiring rather than in any one function, and a lifted fragment of it would be worth nothing.

Why it was withdrawn is [ADR-006](../../docs/adr/006-the-matlab-frontend-is-withdrawn.md). What it
was meant to be is [ADR-005](../../docs/adr/005-logicblock-tree-exchange.md), kept and marked
superseded. `matlab_interoperability.md` here is the user documentation as it stood.

## What is here, and where each path came from

| Here | Was | What it is |
|---|---|---|
| `matlab/+seqcraft/` | `matlab/+seqcraft/` | `LogicBlock` (a handle class), the abstract `Module`, `barrier`, `compile`, `writeLBTX`, and the private encoders |
| `matlab/tests/` | `matlab/tests/` | The MATLAB unit suite: handle semantics, shallow copy, the `Module` contract, and the writer against native `mr.make*` events |
| `examples/` | `examples/matlab/` | The same small GRE 2D twice — direct `LogicBlock` assembly, and through the example-only `GRE2DTR` |
| `python/exchange/` | `src/seqcraft/exchange/` | The LBTX JSON Schema, the generic reader, and the structural validator |
| `python/cli.py`, `python/__main__.py` | `src/seqcraft/cli.py`, `src/seqcraft/__main__.py` | The `compile-tree` process boundary MATLAB shelled out to |
| `tests/exchange/` | `tests/exchange/` | Codec, CLI and `crossval` MATLAB tests |
| `tests/fixtures/lbtx/` | `tests/fixtures/lbtx/` | `minimal.lbtx.json`, the handwritten document — see below |
| `matlab_interoperability.md` | `docs/matlab_interoperability.md` | Setup, the user contract, and how to run the cross-language tier |

Reviving it is those moves in reverse, plus the four things the removal took out of the build:
`jsonschema>=4.21` in `[project.dependencies]`, the `seqcraft = 'seqcraft.cli:main'` console script,
the wheel's `force-include` of `lbtx-1.schema.json`, and the two `.gitignore` negations that keep a
schema and a fixture from being swept up by `*.json`.

## What it does, in one line

MATLAB builds a tree of official `mr.make*` events, `writeLBTX` writes it as JSON with field names
normalized to snake_case, `python -m seqcraft compile-tree` compiles it with the ordinary Python
compiler, and `mr.Sequence.read` loads the `.seq` back. One direction only: MATLAB writes, Python
reads. The Python compile path never used it.

## What a revival has to answer first

The code works. These do not have answers, and finding three good ones is the actual cost:

1. **How is a module called?** In Python a module is a callable returning one block, and `__call__`
   is the one place the result is checked and named. MATLAB has no equivalent that is not a trap —
   emulating it means overloading `subsref`, which breaks object arrays — so this implementation
   split one lifecycle into a public `build` and a protected `buildImplicit`. That is a second name
   for one idea, in the class whose whole purpose is to have one. If it comes back, decide this
   before writing a second MATLAB module, not after. `buildImpl` at least matches what MATLAB users
   already know from `matlab.System`.
2. **How is the wire format tested on every pull request?** These MATLAB tests are marked
   `crossval` because public runners carry no MATLAB licence, so the only guard in the normal tier
   was `tests/fixtures/lbtx/minimal.lbtx.json` — **written by hand**, so it asserts what we believed
   MATLAB emits rather than what MATLAB emits. Two routes out: capture that fixture from a real
   MATLAB run with a checked-in script, or run MATLAB in CI with MathWorks' GitHub Action, which is
   free for public repositories. The second was never tried.
3. **Who maintains it?** A second frontend is a second implementation of the same concepts, with its
   own review and its own drift.

Two smaller faults to fix on the way back in:

- `writeLBTX` hardcodes `"pulseq": "1.5"` and the schema requires exactly that string, so the
  compatibility version is a claim and never a measurement. Read the installed Pulseq version.
- `examples/+seqcraft_examples/GRE2DTR.m` is named after the Python kernel and documented as its
  counterpart, but it designs every event inline instead of composing leaves, spoils two cycles per
  voxel where Python defaults to four, has no `te_s` and no RF-spoiling phase, and no test ever
  compared it against `sc.modules.GRE2DTR` — only against its sibling script. Either make it
  compose, or stop it claiming the name.

## Running it, if you do

```bash
git mv salvage/matlab-frontend/python/exchange src/seqcraft/exchange      # and the rest, per the table
export PULSEQ_MATLAB_PATH=/path/to/pulseq/matlab
pytest -m crossval tests/exchange/test_matlab.py
matlab -batch "addpath(getenv('PULSEQ_MATLAB_PATH')); addpath('matlab'); \
  results=runtests('matlab/tests'); assert(all([results.Passed]));"
```

Until then: read it, copy from it, port it. Do not import it from `seqcraft`.
