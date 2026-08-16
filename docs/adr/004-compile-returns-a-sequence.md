# ADR-004: `compile` returns a `pypulseq.Sequence`, and every legality failure raises

- Status: Accepted
- Date: 2026-08-15
- Supersedes: the `CompiledSequence` / `Report` result model from ADR-003 and before
- Amends: [ADR-003](003-scanner-and-module-reform.md), whose Decision section describes
  `CompiledSequence` storing `opts` and a `pns(hardware)` method on it

## Context

`sc.compile` returned a `CompiledSequence`: the `pypulseq.Sequence`, a `Report`, per-block
provenance, the definitions, and the questions you ask afterwards — `check()`, `moments()`,
`kspace()`, `pns()`, `write()`. Findings were `Issue` values with a severity, and the contract was
"hard failures raise, soft findings report".

That split was inherited from a real problem. The reference implementation's `get_report()` did
this:

```python
ok, error_report = self.seq.check_timing()
if ok:
    print("Timing check passed successfully")
else:
    print("Timing check failed. Error listing follows:")
    [print(e) for e in error_report]
```

It printed and returned `None`, so a failing sequence could not be detected in code at all, and in
a notebook the message scrolled away. `Report` was the answer to that: a value, assertable,
renderable, writable into the provenance sidecar.

**The answer was one indirection short.** A `Report` a caller can decline to read fails the same
way `print` does. `out.check()` had to be *called*; `raise_if_failed()` had to be called after it.
Neither happens in a notebook cell that ends with `out`. The failure mode was unchanged: a `.seq`
gets written, the console refuses it an hour later, and the explanation is on an object nobody
looked at.

Three further observations decided it.

**The soft half turned out not to exist.** Going through the `Issue` kinds one by one: `grad_limit`,
`slew_limit`, `adc_samples_limit`, `rf_samples_limit`, `label`, `timing`, `duration`, `moment` and
`address` are all *errors* — they mean the sequence is not legal, and there is nothing to hand back.
`grad_merge`, `grad_resample`, `raster` and the two `_norm` kinds are not findings about legality at
all; they say what the compiler **did**. Nothing was left in the middle. The category "soft finding
worth reporting on a legal sequence" was empty.

**Python already has the mechanism for the second half.** A `SeqCraftWarning` goes through
`warnings`, which means `simplefilter('error')` makes it fatal, `-W` controls it from the command
line, `pytest.warns` asserts on it, and logging captures it — none of which a bespoke `Report` type
gets for free.

**The wrapper cost more than it carried.** `CompiledSequence` and `WriteResult` were 400 lines;
`report.py` was 190; the provenance sidecar was 168. Between them they justified a whole `result/`
package, a layering rule to keep `result/` from importing `compiler/`, and roughly 160 call sites
across the tests and examples that said `out.seq` where they meant the sequence.

## Decision

**`sc.compile(root, opts, *, name='', definitions=None) -> pypulseq.Sequence`.**

Every legality failure raises. `CompileError` for what the compiler cannot resolve by scheduling,
`HardwareLimitError` for what the machine cannot play, `DefinitionConflict` for two sources claiming
one key, `CompilerContractError` for a compiler bug. Each message names the offending number, the
tag path it came from, the time it happens, and two remedies with the values already computed.

Everything the compiler *did* is a `SeqCraftWarning`, **one aggregated warning per category**.
Aggregation is load-bearing rather than cosmetic: Python's default filter shows a warning once per
unique `(message, category, module, lineno)`, so one `warn` per merge would print the first and
silently swallow the other eleven — strictly worse than the count the report gave.

The definitions are applied to the sequence **during the compile**, which is what makes the
returned object self-sufficient: nothing has to survive until write time to put them there, so
there is nothing left for a wrapper to carry.

`result/`, `report.py`, `testing.py` and `design/sampling.py` are deleted. `analysis.py` takes the
four measurements — `sample`, `moments`, `kspace`, `pns` — under one entry shape: give it a tree,
get numbers back.

**An exception lives with the code that raises it.** Only what more than one package raises stays in
`errors.py`. Every one is re-exported at the root, and the layering test asserts the *identity*, so
a caller's `except sc.CompileError` is unchanged.

## Consequences

**The compiled bytes do not change.** This alters what is *returned* and how findings are
*delivered*, never what is emitted. `build_gre` and `build_se` write `.seq` files with the same
sha256 as before the change, and every structural field of the frozen baseline matches. That was
checked at each of the three gates and is asserted on every test run.

**Provenance is no longer returned.** `out.origin(i)` is gone. The tag paths are still computed and
are what every error message quotes — `from: tr.readout.prephaser` rather than `block 117` — but a
caller cannot enumerate them. Nothing was using them to do anything but assert on them.

**The provenance sidecar is a real regression**, and the only one. Nothing currently records the git
commit, dirty flag or package versions that produced a `.seq`. It is deferred deliberately rather
than overlooked: see [`serialization.md`](../serialization.md) for what a replacement has to decide,
and why designing the weak version first is how the strong one stops being written.

**`.check()` cannot come back.** Re-attaching findings to the returned object is the obvious way to
"improve" this API, and it would restore exactly the failure mode being removed.
`tests/test_layering.py::test_compile_returns_a_bare_pypulseq_sequence` asserts that the returned
object has neither a `report` nor a `check` attribute, so that regression is a failing test rather
than a review comment.

**Two defects surfaced while making the change**, both of which the old shape had hidden:
`make_arbitrary_grad` without `first=`/`last=` extrapolates from the end samples, so a test tree's
"lone spiral" started at −2680 Hz/m — a step from zero the amplifier cannot play, which compiled
silently until `check_timing` moved inside the compile. And the event-size check indexed a 0-based
`origins` list with pypulseq's 1-based block id, so it named the wrong block and could raise
`IndexError` on the last one.
