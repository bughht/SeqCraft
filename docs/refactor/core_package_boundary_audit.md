# Core package boundary audit charter — superseded

- Status: **Superseded.** Closed by the structure revision, which carried out the move this
  charter existed to authorise.
- Superseded by: [`docs/architecture.md`](../architecture.md) for the resulting layout, and
  `tests/test_layering.py` for the rule that keeps it.

## What it asked

Whether every module in `src/seqcraft/core/` had a coherent responsibility and dependency position
once the compiler had been turned into an explicit pipeline — and it deliberately withheld
implementation authority until that question had an evidence-backed answer.

## What the answer turned out to be

`core` did not. Its membership rule — *"what is required to get a logic block to a validated
`.seq`"* — is a property of the whole package rather than of any file in it, so it admitted the
scheduler, the unit table, the geometry, the report type and the exception hierarchy into one
4 700-line directory and gave no reason to keep any of them out. The audit's own framing names the
failure: it asked whether each module had "a coherent dependency position", which is a question a
directory called `core` cannot answer.

The replacement is four packages named for four questions — `scanner/`, `design/`, `compiler/`,
`result/` — ordered by the one direction the dependencies run. That rule is mechanical, so
`tests/test_layering.py` checks it per file rather than a document asserting it.

## What is worth keeping from it

The two mistakes it was written to prevent are both still real, and both still apply to whatever
moves next:

- retaining unrelated code in a package merely because it is already there;
- moving cohesive low-level code merely to make the directory look different.

The second is why `design/timing.py`, `design/units.py` and `design/logic.py` moved without a line
changing, and why `display.py` was **not** split into a package: three public functions with no
callers is not a subsystem yet.
