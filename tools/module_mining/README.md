# Module mining -- experimental

Tooling for the AI-assisted module-mining project: it makes external sequence implementations
**measurable**, so that a candidate abstraction can be argued from numbers instead of from
source-code similarity.

Nothing here is part of SeqCraft's public API, nothing in `src/seqcraft` imports it, and the
layout, names and return shapes are free to change.

## Which copy to edit

[`plans/`](plans/) is a **mirror**, taken 2026-09-19, of the workspace planning repository's
`docs/plans/module-mining/`. It is here so that the pull request promoting `TSEShot`, `FSE2D`,
`RadialReadout` and `GRE3DTR` carries the evidence that justified each boundary, and so a reviewer
does not need a second repository to read why the lines were drawn where they are.

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
  inspect_emitted.py              read a .seq back and report the physics it actually encodes
  candidates/tse/
    adapters/seqcraft_fse.py      R1  the notebook-local FSE2D, executed rather than copied
    adapters/pypulseq_tse.py      R2  PyPulseq write_tse.py, loaded by path and called
    adapters/matlab_seq.py        R3  a .seq written by Pulseq MATLAB, read back
    run_ringdown.py               does the references' construction assume dead == ringdown?
    run_comparison.py             both Python references at one protocol
    run_metamorphic.py            turbo factor as the transformation
    run_sweep.py                  the accepted invariants across the parameter grid
  candidates/radial/
    adapters/pypulseq_radial.py   the official write_radial_gre.py, loaded and called
    adapters/seqcraft_radial.py   RadialReadout at the same protocol
    checks.py                     the trajectory invariants, shared by the two scripts below
    run_archaeology.py            what the references agree on before any candidate exists
    run_validation.py             the candidate against the reference, sample for sample
  candidates/gre3d/
    run_validation.py             the three-axis delta-k comparison and the TE regression
```

`inspect_emitted.py` prints rather than asserts, and it is the newest piece here because of the
mistake it was written for. A `GRE3DTR` can satisfy every k-space and timing invariant, compile
legally, pass its whole suite, and still emit the **wrong pulse family** -- the first one played a
3 ms shaped sinc with the selection gradient removed and called that non-selective. No invariant
in `fingerprint.py` looks at pulse shape, because none of them had a reason to. So:

```bash
python tools/module_mining/inspect_emitted.py examples/gre_3d/seq/*.seq
```

reads the file back with no access to the Python that wrote it and reports, per RF event, whether
it is hard or shaped, how long it is, and which axes carry gradients while it plays. Run it on a
candidate's output before promoting anything whose API name makes a claim about the physics.

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

The first is blocking because **the comparator has been wrong four times and the candidates
none**. Its tests use answers known independently of any candidate: an analytic gradient moment, a
synthetic rotation, a seeded error of chosen size, and each of the four historical bugs as a
regression test -- the echo located by `argmin|kx|` on a reference with no DC sample, an equal-area
rule applied to the phase axis, bare spokes stacked so their trajectories accumulated, and a shot
index assumed to exist for a leaf with no excitation.

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
python tools/module_mining/candidates/tse/run_metamorphic.py
python tools/module_mining/candidates/tse/run_sweep.py
python tools/module_mining/candidates/radial/run_archaeology.py
python tools/module_mining/candidates/radial/run_validation.py
python tools/module_mining/candidates/gre3d/run_validation.py
```

Only the TSE scripts need the PyPulseq checkout's `write_tse.py`; `radial/` needs its
`write_radial_gre.py`, and `gre3d/` needs no external checkout at all.

Each writes a `*_result.json` beside itself. Those files are **not** tracked -- `.gitignore`
excludes `*.json` repo-wide, and scan output is exactly what that rule is for. They are
regenerable by re-running the script, and the reports quote the numbers that matter.
