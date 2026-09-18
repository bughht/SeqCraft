# Module mining -- experimental

Tooling for the AI-assisted module-mining project: it makes external sequence implementations
**measurable**, so that a candidate abstraction can be argued from numbers instead of from
source-code similarity.

Nothing here is part of SeqCraft's public API, nothing in `src/seqcraft` imports it, and the
layout, names and return shapes are free to change. The plan it serves lives in the workspace
repository under `docs/plans/module-mining/`.

## Layout

```text
module_mining/
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

## Running

Needs the `seqcraft-dev` environment, `nbformat`, and a PyPulseq checkout for R2
(`MODULE_MINING_PYPULSEQ`, default recorded in the reference inventory).

```bash
python tools/module_mining/candidates/tse/run_ringdown.py
python tools/module_mining/candidates/tse/run_comparison.py
```

Each writes a `*_result.json` beside itself. Those files are committed as evidence and are
regenerable; they are not inputs to anything.
