# Module mining -- experimental

Tooling for the AI-assisted module-mining project: it makes external sequence implementations
**measurable**, so that a candidate abstraction can be argued from numbers instead of from
source-code similarity.

Nothing here is part of SeqCraft's public API, nothing in `src/seqcraft` imports it, and the
layout, names and return shapes are free to change.

## Which copy to edit

[`plans/`](plans/) is a **mirror**, taken 2026-09-18, of the workspace planning repository's
`docs/plans/module-mining/`. It is here so that the pull request extracting `TSEShot` and `FSE2D`
carries the evidence that justified it, and so a reviewer does not need a second repository to
read why the boundary was drawn where it is.

The workspace copy remains the **working** one: that is where the plan is revised, where a new
candidate's evidence is written, and what the workspace docs index links to. Treat this copy as a
snapshot of the argument at the time of the PR, and re-mirror rather than edit in place — two
copies of a living document diverge, and the one in the code repository is the one nobody
remembers to update.

## Layout

```text
module_mining/
  plans/                          the plan and the evidence -- see "Which copy to edit" below
  reference.py                    one record per executed reference implementation
  fingerprint.py                  the comparator stack: L0 / L2 / L3 / L4
  candidates/tse/
    adapters/seqcraft_fse.py      R1  the notebook-local FSE2D, executed rather than copied
    adapters/pypulseq_tse.py      R2  PyPulseq write_tse.py, loaded by path and called
    adapters/matlab_seq.py        R3  a .seq written by Pulseq MATLAB, read back
    run_ringdown.py               E1  does the references' construction assume dead == ringdown?
    run_comparison.py             E0/E2/E3  both Python references at one protocol
```

## Rules the adapters follow

- **Execute the reference; never reimplement it.** The SeqCraft adapter runs the notebook's own
  cells. The PyPulseq adapter loads the example script from disk and calls its `main`.
- **Do not edit the source to make it comparable.** Where a reference hard-codes something the
  experiment needs to vary, the *harness* intervenes for the duration of one call --
  `pypulseq_tse.overridden_opts` patches the `Opts` factory and puts it back.
- **Preserve oddities as evidence.** A reference that warns, or that cannot be built at a given
  setting, is recorded as that outcome rather than worked around.

## What runs in CI, and what does not

Three different jobs, deliberately not conflated:

| | where |
|---|---|
| **is the comparator correct?** | `tests/module_mining/`, in the ordinary blocking test job |
| **does a promoted module still hold its accepted invariants?** | `tests/modules/`, beside every other module test |
| **do the external references still agree?** | **here**, run by hand or on a schedule — non-blocking |

The first is blocking because the comparator has been wrong three times and the candidates none.
Its tests use answers known independently of any candidate: an analytic gradient moment, a
synthetic rotation, a seeded error of chosen size, and the three historical bugs as regression
tests.

The third -- **external reference re-validation**: rerunning a candidate's references,
adapters and comparator end to end -- is not blocking because it depends on external checkouts, their versions, adapter
assumptions, large sweeps and candidate-specific tolerances. A red merge gate driven by a
measurement we are less sure of than the code is how correct production code gets "fixed" to
satisfy a broken instrument. If one comparison path proves stable across several candidates, it
can be promoted then.

**No candidate gets its own CI workflow.** Tests grow; CI *types* should not. A promoted module's
stable invariants become ordinary tests in `tests/modules/` — `test_radial_readout.py` carries the
rotation-equivariance check that was discovered here, and nothing about radial needed a new job.

## Running

Needs the `seqcraft-dev` environment, `nbformat`, and a PyPulseq checkout for R2
(`MODULE_MINING_PYPULSEQ`, default recorded in the reference inventory).

```bash
python tools/module_mining/candidates/tse/run_ringdown.py
python tools/module_mining/candidates/tse/run_comparison.py
```

Each writes a `*_result.json` beside itself. Those files are committed as evidence and are
regenerable; they are not inputs to anything.
