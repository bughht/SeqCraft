# Testing and continuous integration

SeqCraft tests real PyPulseq events and emitted sequences; compiler tests do not mock PyPulseq.
The CI environment constrains dependencies through `ci/constraints.txt` so an unchanged commit is
not made red by a new lint rule, type-checker release, or incompatible PyPulseq wheel.

SeqCraft currently requires Pulseq 1.5.1 behavior from the compatibility fork pinned in
`pyproject.toml` and `ci/constraints.txt`.  Its package metadata is the same `1.5.0.post1` used by
the incompatible PyPI wheel, so the immutable source URL is part of the dependency contract.

## Pull-request CI

Every push to `main` and every pull request runs five gates:

1. `lint`: the established Ruff correctness baseline on Python 3.14.
2. `types`: strict mypy checking of the pure-arithmetic core, hosted on Python 3.14.
3. `checks`: Linux with Python 3.14 -- pytest with coverage, the source doctests, the API
   reference and the prose attribution.
4. `test`: pytest on Linux with Python 3.10 and on Windows with Python 3.14.
5. `examples`: isolated execution of the **build-only** notebooks on Linux with Python 3.14 --
   the explicit allowlist `_NOTEBOOKS` in
   [`tools/run_notebook_smoke.py`](../tools/run_notebook_smoke.py), currently 22 notebooks.
   The allowlist is the source of truth and is not a filename pattern: it is every notebook
   that needs only seqcraft and matplotlib, which includes `01_getting_started.ipynb`, the
   `01_build.ipynb` of each example, and the `03` notebooks that are also build-only
   (`fse_2d/03_module_api`, and `03_flow_comp` under `gre_2d`, `gre_spiral_2d` and `megre_2d`).
   What is deliberately absent is the simulation tier -- every `02_simulate_and_reconstruct`
   and `bssfp_2d/03_flow_and_motion` -- which needs MRzeroCore, torch and sigpy plus a phantom
   download, and so belongs to the lab tier.

### What the matrix covers

The package suite runs in **three environments**:

```text
Linux   3.10   in `test`      minimum supported Python
Linux   3.14   in `checks`    current primary Python
Windows 3.14   in `test`      OS portability, at the primary Python
```

This is **representative compatibility coverage, not an exhaustive product** of the supported
dimensions. Each lane changes one primary compatibility dimension relative to the primary
environment, which makes a failure easier to localise:

```text
Linux 3.10  <-> Linux 3.14       Python differs, OS held    -- the supported version endpoints
Linux 3.14  <-> Windows 3.14     OS differs, Python held    -- OS portability
```

"Easier to localise" is the honest claim rather than "attributable". The 3.10 lane is a whole
**environment**, not only an interpreter: 3.10 also resolves older NumPy and matplotlib, so a
failure there narrows the search without proving the interpreter caused it.

**A combination SeqCraft supports is not required to appear here.** `requires-python` is `>=3.10`
with no upper bound; 3.11, 3.12 and 3.13 are supported and simply are not executed. Windows 3.11
used to run and was dropped on the same reasoning: once the minimum-Python Linux lane and the
primary-Python Windows lane are both present, it adds relatively little independent compatibility
information for the cost of another full-suite execution. That is a judgement about marginal
value, not a claim that it could detect nothing.

Python 3.10 is the `requires-python` floor, so **CI runs the minimum version we publicly claim to
support.** It is a compatibility lane like any other -- the same fast suite, without the
once-only validations. Being the minimum-version gate does not earn it coverage, doctest, API or
prose work.

Python 3.14 is the current stable endpoint, verified rather than assumed: the suite, the
notebooks and all four once-only validations were run there before it became the primary
environment.

### NumPy is pinned per interpreter

NumPy 2.3.0 raised its own support floor to Python 3.11 while ours stayed at 3.10, so a single
unconditional pin cannot install on the minimum Python we claim to support. `ci/constraints.txt`
uses environment markers:

```text
Python 3.10     ->  numpy 2.2.6    the last release supporting 3.10
Python >=3.11   ->  numpy 2.4.6    including 3.14, which resolves cp314 wheels
```

Keep the two markers mutually exclusive and exhaustive when bumping either line. Dropping the
marked row does **not** fail the 3.10 lane -- pip falls back to the `numpy>=1.24` floor in
`pyproject.toml`, and CI silently validates a version nobody chose.

### What runs once, and why

Coverage, the doctests, the API reference and the prose attribution run only in `checks`, on the
primary environment. That is a de-duplication decision rather than a claim that they could never
behave differently elsewhere -- the doctests and the API-reference blocks are executable code.
The package's OS and Python compatibility stays covered by the three suite environments above.

`[tool.mypy]` sets `python_version = '3.10'`, which keeps mypy's target-language semantics at the
declared minimum whatever version hosts the `types` job. The host still supplies the installed
stubs mypy reads, so it is not irrelevant -- it is simply not what fixes the target.

A pull request is also one concurrency group, so a new commit cancels the run for the commit it
supersedes.  Pushes to `main` group individually and cancel nothing.

The notebook runner copies `examples/` to a temporary directory before execution. Generated
sequence files therefore never modify the working tree.

Inside that copy it runs the allowlist **four notebooks at a time**. Almost half of a single
notebook's cost is fixed startup -- a kernel is 0.7 s and its imports another 1.7 s -- so the set
is dominated by per-notebook overhead rather than by any one notebook, and overlapping them is
what shortens the job. The bound is a ceiling rather than a target: each worker holds a kernel
subprocess with its own NumPy and matplotlib, so the limit is runner memory and CPU, not this
process.

Concurrency is safe here because the notebooks are independent, which is a property of the
notebooks and not an assumption: none reads another's output, no two write the same filename,
and the `seq/` directory they create is created with `exist_ok=True`.

**That independence is what the allowlist requires, and it is a precondition rather than a
preference.** A future notebook that depends on another notebook's output would violate it, so
it cannot simply be appended to the allowlist and left to run unordered. It would need explicit
dependency ordering, serialisation, or another ordering-aware validation path — which it may
well get while remaining a Layer-2 smoke notebook.

A failing notebook no longer stops the run. With several kernels already in flight, where a
serial run would have stopped is a property of the runner's CPU that minute, so every scheduled
notebook finishes and every failure is reported -- in allowlist order, with the failing cell and
the kernel's traceback.

To reproduce the main gates locally from the repository root:

```bash
python -m pip install --constraint ci/constraints.txt -e ".[dev,viz,rf]"
ruff check .
mypy src/seqcraft/design/timing.py src/seqcraft/design/units.py \
  src/seqcraft/design/module.py \
  src/seqcraft/scanner/opts.py \
  src/seqcraft/compiler/model.py src/seqcraft/compiler/placement.py \
  src/seqcraft/compiler/boundaries.py src/seqcraft/compiler/legalization.py \
  src/seqcraft/compiler/emission.py
pytest -n auto \
  -m "not slow and not bloch and not crossval and not hardware" \
  --cov=seqcraft --cov-report=term-missing --durations=20
pytest --doctest-modules src/seqcraft
python tools/check_api_reference.py
python tools/check_prose_attribution.py
python tools/run_notebook_smoke.py
```

That sequence is the `lint`, `types`, `checks` and `examples` gates together.  The `test` lanes
run the same pytest invocation without `--cov`, which is most of why they finish sooner.

`check_prose_attribution.py` blocks a short list of specific unsupported attributions and prints
its broader review words under `--review`, where they never fail.  See
[`writing_a_module.md`](writing_a_module.md).

`check_api_reference.py` executes every fenced `python` block in
[`api_reference.md`](api_reference.md) in one shared namespace, in document order, and checks its
index against every module's `__all__` in both directions -- a public name missing from the index,
or an indexed name the package does not export, both fail. Blocks that are illustrative rather than
runnable are skipped by name and each skip is printed, so the exemptions stay visible instead of
accumulating quietly. A reference nobody executes is a reference that is wrong within two commits.

Ruff formatting remains a pre-commit hook.  It is not yet a repository-wide CI gate because the
existing tree has not had a dedicated format-only migration.  Likewise, mypy is intentionally
limited to the modules that the configuration describes as structurally typeable; expand that
list only after bringing an additional module to a clean baseline.  The README contains fenced
examples rather than doctests, so the former empty pytest invocation is not a gate.

## Tiers outside pull-request CI

Simulation and reconstruction notebooks require MRzeroCore, Torch, SigPy, and artifacts generated
by the build notebooks. Vendor hardware checks additionally require a site's own `.asc` file,
which the caller passes explicitly and which is never committed here. Those tiers run on a controlled lab or nightly environment,
not on public GitHub-hosted runners.

The public CI result therefore proves the compiler, the `Module` contract, documentation snippets,
file round-trip, and the getting-started notebook on the three suite environments above. It does
not prove vendor hardware, full simulation, reconstruction, or external cross-validation.

## What the fixtures are made of

Every compiler and integration fixture is **raw pypulseq**. That is a standing rule, not an
accident of the current tree: the previous suite built its realistic trees out of module-library
classes, which made compiler coverage depend on whatever the library happened to contain and made
the library impossible to replace without also rewriting the compiler's tests. `tests/conftest.py`
supplies one `pp.Opts`, and the sequences are assembled from `pp.make_*` calls.

`tests/conftest.py` also holds the component assertions the package used to ship as
`seqcraft.testing`: `assert_deterministic`, `assert_pure`, `assert_output` and `assert_all`,
reachable through the `component_checks` fixture. They are here rather than in `src/` because a
package that ships assertions has to keep them working, and these are forty lines that only this
repository's own module tests use. Everything else that module asserted -- the raster, the limits,
that a block is well formed -- the compiler now checks, with a better message and on the *summed*
waveform.

## What a compiler test asserts

`sc.compile` raises on every legality failure and returns a bare `pypulseq.Sequence` otherwise, so
**a test that compiles at all has asserted most of what it used to assert explicitly**: limits on
the summed waveform, per-event sample counts, duplicate k-space addresses, pypulseq's own timing
audit, and the four against-the-tree invariants. What is left to write down is the number the test
is actually about.

Warnings are the other half. Use `pytest.warns(sc.SeqCraftWarning, match=...)` for the positive
form; for the negative -- *no* warning of a kind -- use `warnings.catch_warnings(record=True)` and
assert on the captured list, because `pytest.warns` has no clean negative.

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

The compiler's responsibilities and the event/boundary support matrix are in
[`compiler.md`](compiler.md).  The two `refactor/` documents that used to hold them --
[`compiler_current_state.md`](refactor/compiler_current_state.md) and
[`compiler_constraint_matrix.md`](refactor/compiler_constraint_matrix.md) -- are dated,
commit-pinned records of the pre-refactor shape and are marked historical; they are kept because
what changed is only legible against what was there.

The always-on verification cost, warning/error ownership and integer-tick boundary were reviewed
after the legalization/emission split in
[`refactor/phase6_verification_time_audit.md`](refactor/phase6_verification_time_audit.md).
The final public/private inventory, legacy-path audit and source dependency guards are recorded in
[`refactor/compiler_architecture_freeze.md`](refactor/compiler_architecture_freeze.md).
