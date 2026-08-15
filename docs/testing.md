# Testing and continuous integration

SeqCraft tests real PyPulseq events and emitted sequences; compiler tests do not mock PyPulseq.
The CI environment constrains dependencies through `ci/constraints.txt` so an unchanged commit is
not made red by a new lint rule, type-checker release, or incompatible PyPulseq wheel.

SeqCraft currently requires Pulseq 1.5.1 behavior from the compatibility fork pinned in
`pyproject.toml` and `ci/constraints.txt`.  Its package metadata is the same `1.5.0.post1` used by
the incompatible PyPI wheel, so the immutable source URL is part of the dependency contract.

## Pull-request CI

Every push to `main` and every pull request runs four gates:

1. `lint`: the established Ruff correctness baseline on Python 3.11.
2. `types`: strict mypy checking of the pure-arithmetic core on Python 3.11.
3. `test`: pytest and source doctests on Linux and Windows with Python 3.11 and 3.12.
4. `examples`: isolated execution of the getting-started notebook on Linux with Python 3.11.
   The two DTI build notebooks are parked under `examples/_parked/` until a module library exists
   to rebuild them against, and are deliberately not executed.

The notebook runner copies `examples/` to a temporary directory before execution. Generated
sequence files therefore never modify the working tree.

To reproduce the main gates locally from the repository root:

```bash
python -m pip install --constraint ci/constraints.txt -e ".[dev,viz,rf]"
ruff check .
mypy src/seqcraft/design/timing.py src/seqcraft/design/units.py \
  src/seqcraft/design/module.py \
  src/seqcraft/scanner/opts.py \
  src/seqcraft/compiler/model.py src/seqcraft/compiler/placement.py \
  src/seqcraft/compiler/boundaries.py
pytest -n auto \
  -m "not slow and not bloch and not crossval and not hardware" \
  --cov=seqcraft --cov-report=term-missing
pytest --doctest-modules src/seqcraft
python tools/run_notebook_smoke.py
```

Ruff formatting remains a pre-commit hook.  It is not yet a repository-wide CI gate because the
existing tree has not had a dedicated format-only migration.  Likewise, mypy is intentionally
limited to the modules that the configuration describes as structurally typeable; expand that
list only after bringing an additional module to a clean baseline.  The README contains fenced
examples rather than doctests, so the former empty pytest invocation is not a gate.

## Tiers outside pull-request CI

Simulation and reconstruction notebooks require MRzeroCore, Torch, SigPy, and artifacts generated
by the build notebooks. Vendor hardware checks additionally require site-confidential `.asc` files
provided through `SEQCRAFT_ASC_DIR`. Those tiers run on a controlled lab or nightly environment,
not on public GitHub-hosted runners.

The public CI result therefore proves the compiler, the `Module` contract, documentation snippets,
file round-trip, and the getting-started notebook on the supported matrix. It does not prove vendor
hardware, full simulation, reconstruction, or external cross-validation.

## What the fixtures are made of

Every compiler and integration fixture is **raw pypulseq**. That is a standing rule, not an
accident of the current tree: the previous suite built its realistic trees out of module-library
classes, which made compiler coverage depend on whatever the library happened to contain and made
the library impossible to replace without also rewriting the compiler's tests. `tests/conftest.py`
supplies one `pp.Opts`, and the sequences are assembled from `pp.make_*` calls.

## Updating dependencies

Dependency updates are intentional maintenance changes. Update one related group in
`ci/constraints.txt`, run every gate above, and record any changed warning, error, waveform,
block-count, or notebook behavior. The constraints file pins direct CI-critical dependencies;
transitive packages remain resolver-managed until a compatibility issue justifies locking them.

## Compiler refactor baseline

The original Phase 0 capture is documented in
[`refactor/phase0_baseline.md`](refactor/phase0_baseline.md). Its machine-readable structural and
performance artifact lives at `tests/baselines/compiler_phase0.json`; the integration suite checks
the stable subset on every run.

**It was re-captured when the recipes moved to raw pypulseq**, so it no longer spans the Phase 0
boundary — the four module-built recipes it froze no longer exist. What it guards from here is that
the remaining refactor phases change block counts, boundaries and moments not at all. Regenerate it
only after reviewing an approved behavior change:

```bash
python tools/capture_compiler_baseline.py --iterations 3
```

Current compiler responsibilities and the event/boundary support matrix are recorded in
[`refactor/compiler_current_state.md`](refactor/compiler_current_state.md) and
[`refactor/compiler_constraint_matrix.md`](refactor/compiler_constraint_matrix.md).
