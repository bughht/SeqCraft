# ADR-006: The MATLAB frontend is withdrawn, and LBTX with it

- Status: Accepted
- Date: 2026-09-14
- Supersedes: [ADR-005](005-logicblock-tree-exchange.md)

## Context

[ADR-005](005-logicblock-tree-exchange.md) added a MATLAB frontend and a one-way exchange format:
`matlab/+seqcraft/`, `src/seqcraft/exchange/`, a `compile-tree` CLI, paired MATLAB GRE examples, and
their tests. It was deliberately the smallest thing that could work, and it did work -- MATLAB built
a tree of official `mr.make*` events, Python compiled it, and an `mr.Sequence` came back.

The reason to withdraw it is not that it was built badly. It is that **the concept it exists to
carry does not survive the language boundary**, and that the boundary cost is paid on every change
to the Python side.

**A module is a callable that returns one block, and MATLAB cannot say that.** In Python,
`readout(line=17)` reads correctly beside `add`, and `__call__` is the one place the framework
checks and names what came back. MATLAB has no equivalent that is not a trap: emulating it means
overloading `subsref`, which breaks object arrays and indexing. So the MATLAB port split the one
lifecycle into a public `build` and a protected `buildImplicit` -- a second name for one idea, in
the one class whose whole purpose is to have a single idea. Every later module would have paid that
tax, and every MATLAB user would have had to learn a contract that does not match the documented
one.

**Nothing MATLAB ran in CI, and nothing could.** The MATLAB tests are `crossval` because public
runners carry no MATLAB licence. The only thing guarding the wire format in the normal tier was
`tests/fixtures/lbtx/minimal.lbtx.json`, which was written by hand -- so it asserts what we believed
MATLAB emits, not what MATLAB emits. A field renamed in `mr.make*`, or a `camelToSnake` result that
collides, fails on a workstation weeks later or on a user's machine, not in the pull request that
caused it.

**The compatibility claim was a constant.** `writeLBTX` hardcoded `"pulseq": "1.5"` and the schema
required exactly that string. No code read the installed MATLAB Pulseq version, so an older Pulseq
declared 1.5 and was believed.

**The example did not demonstrate what it claimed.** `seqcraft_examples.GRE2DTR` was named after the
Python kernel and documented as corresponding to it, but it designed every event inline rather than
composing leaves, spoiled two cycles per voxel where Python defaults to four, had no `te_s` and no
RF-spoiling phase, and was never compared against `sc.modules.GRE2DTR` by any test -- only against
the sibling script. Two things with one name and different output is the failure mode the module
layer exists to prevent.

**And there is one maintainer.** A second frontend is a second implementation of the same concepts
with its own review, its own CI story and its own drift, and the evidence that MATLAB users need it
has not arrived.

## Decision

Take the MATLAB frontend and the exchange boundary out of the package, and park them whole in
[`salvage/matlab-frontend/`](../../salvage/matlab-frontend/): `matlab/`, `examples/matlab/`,
`src/seqcraft/exchange/`, `src/seqcraft/cli.py`, `src/seqcraft/__main__.py`, `tests/exchange/`,
`tests/fixtures/lbtx/` and `docs/matlab_interoperability.md`. Out of the build with them: the
`jsonschema` dependency, the `seqcraft` console script, the wheel's schema force-include, and the
two `.gitignore` negations that kept a schema and a fixture out of `*.json`.

Parked rather than deleted because the cost of reviving this is not the code — the code works. It
is the three questions below, and a deletion would have thrown away a working implementation
without making any of them easier to answer.

**Python is the only frontend.** `sc.compile` is unchanged, and was never on the LBTX path -- the
compiler, the `Module` contract, the warning policy and the return type are untouched by this
removal, which is the one good consequence of ADR-005 having kept the boundary thin.

A MATLAB user reaches Pulseq the way they did before SeqCraft existed: official MATLAB Pulseq
directly, or a `.seq` file written by Python and read with `mr.Sequence.read`. A `.seq` file is
already the portable artefact between the two languages, and it is the one both toolchains validate.

## Consequences

The exchange format's deferred list -- Python writer, MATLAB reader, round-trip, semantic hash,
provenance, diagnostics protocol -- is void rather than pending. `docs/serialization.md` returns to
saying that nothing there is implemented, which is true again.

If a MATLAB frontend is ever wanted, the implementation is in
[`salvage/matlab-frontend/`](../../salvage/matlab-frontend/), frozen at `d27217d`, with a README
recording where every path came from and what to restore in the build. This document says what a
revival would have to answer first: how a module is called in MATLAB without a second lifecycle,
how the wire format is tested on every pull request rather than on a workstation, and who maintains
it. Moving the code back is a `git mv`; those three answers are the actual cost.
