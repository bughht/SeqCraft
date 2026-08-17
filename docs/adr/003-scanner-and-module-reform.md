# ADR-003: `pp.Opts` as the scanner, and a module contract with no library

- Status: Accepted
- Date: 2026-08-15
- Supersedes: the `System`/`regime` scanner model and the `seqcraft.modules` library
- Partially superseded by: [ADR-004](004-compile-returns-a-sequence.md), which replaces only this
  ADR's `CompiledSequence`, report, provenance-result, and result-method decisions. The scanner,
  PNS hardware input, raw-PyPulseq fixture, and `Module` decisions remain authoritative.

## Context

Three concepts had accumulated that the compile path did not need, and that were coupled to each
other tightly enough that none could be removed alone.

**`core/system.py` (652 lines)** held `System`, `Limits`, scanner presets, named limit regimes,
rasters and PNS hardware loading. Measured against what the compiler actually reads, that is eight
`pp.Opts` fields — `max_grad`, `max_slew`, the four rasters, `adc_dead_time`, `rf_ringdown_time` —
plus the two per-event sample limits. `System` was touched at exactly three call sites outside its
own file. Everything it stored was already a field of the official `pypulseq.Opts`.

**`modules/` (15 files, 5 762 lines, 27 classes)** was a library of excitations, readouts,
encodings and preparations, written on a base class that held the scanner, resolved a limit regime,
ran a unit check behind the subclass's `__init__`, and scraped `__dict__` for provenance. The base
declared no abstract method, so it guaranteed nothing about what a subclass produced.

**The coupling that mattered most** was in the tests: `tests/compiler/test_fidelity.py` built its
realistic trees out of library classes. Compiler coverage therefore depended on whatever the library
happened to contain, and the library could not be replaced without also rewriting the compiler's
tests.

A declared `duration` property completed the trap. Because a caller placed the next thing by it, a
build argument was forbidden from changing a block's length — and the library needed a dedicated
assertion, `assert_duration_is_honest`, purely to catch that property lying.

## Decision

**The scanner is a `pypulseq.Opts`, and nothing else.** `compile_sequence(root, opts, ...)` replaces
`compile_sequence(root, system, regime=...)`. `CompiledSequence` stores `opts`. `pns(hardware)` takes
the response model explicitly, because PNS prediction is analysis rather than compilation. The
provenance sidecar records `vars(opts)`. `core/system.py` is deleted with no compatibility shim: a
`System` forwarding to `Opts` would preserve exactly the concept being removed.

Named regimes are replaced by a second `Opts`. `sc.opts.derate(opts, grad=0.85)` returns one, copying
every field it was not asked to change — which is what the multi-regime consistency check used to
guarantee, now unreachable rather than merely detected.

**`sc.Module` is a contract with four members** — `opts`, `tag`, `__call__`, and an abstract `build`
— living in `src/seqcraft/module.py`, beside the compile path rather than inside it or inside any
library. It declares no `duration`, runs no unit check, and walks no `__dict__`. `opts` is a required
keyword argument, because pypulseq's fallback is the process-global `Opts.default`.

**No concrete modules ship.** The 27 classes are deleted. The physics worth keeping — the exact
b-value integral and its solvers, the variable-density spiral trajectory, the EPI ramp-sampling
moment integral, the DTI direction tables — was lifted out as plain functions into `salvage/` before
the deletion, with no scanner object and no base class in their signatures.

**Every compiler and integration fixture is raw pypulseq**, and that is now a standing rule rather
than a state of the tree.

Two smaller removals follow from the same membership rule: `ordering.py` leaves the package (four of
its six functions had never had a caller, and ordering tables are sequence-programming choices
rather than physics), and `validate.py` loses `check_units`, `require_positive`, `require_int_in`,
`require_divides` and `suggest_field`, whose only consumer was the deleted base.

## Consequences

`sc.compile(tree, system)` and every `sc.modules.*` name are gone, with no deprecation window. That
is deliberate: the project is a 0.3.0 alpha, it has broken cleanly before, and a shim would keep
asserting the model being removed.

The DTI spiral and EPI notebooks stop working, and are parked under `examples/_parked/` rather than
deleted. They are the acceptance test for whatever module set is written next, so rewriting them is
what should *drive* that library rather than follow it.

The Phase 0 compiler baseline was re-captured, because the four module-built recipes it froze no
longer exist. It no longer spans the Phase 0 boundary; what it guards from here is that the
remaining refactor phases change block counts, boundaries and moments not at all.

What is genuinely lost, stated rather than discovered later:

1. **A named scanner catalogue.** `System.preset('cima_x')` is replaced by
   [PulseqSystems](https://github.com/nimpulseq/PulseqSystems) through `sc.opts.from_scanner`. That
   returns *specs* — `max_grad`, `max_slew`, `B0` — and not rasters, dead times, ringdown or sample
   limits, which is why `from_scanner` requires the four site constants as keyword arguments.
2. **Design-time regimes as names.** Now two `Opts` objects held by the caller. The compiler
   validated against one regime anyway, so nothing about the check changes.
3. **The cross-regime consistency guarantee.** Gone with the regimes; one `Opts` cannot disagree
   with itself.
4. **Automatic unit validation on module construction.** The bands and `check_fields` remain for
   dataclasses; calling anything is now the subclass's decision, taken where the plausible range is
   actually known.

## What was considered and rejected

**A `from_specs()` constructor** wrapping `pp.Opts` with the site constants as required arguments.
Rejected as redundant: `Opts(**specs)` already works, and `Opts` already raises `TypeError` on a
mistyped keyword. The requirement survives where it is not redundant — inside `from_scanner`, whose
vendor lookup *cannot* supply dead times, so returning an `Opts` with three zeros in it would be the
shortest available path to a file the console refuses.

**An `Interpreter` or `Site` value object** carrying the site constants. Rejected: a frozen
dataclass whose only job is to carry nine numbers into one function is the kind of thing this reform
is removing, and required keyword-only arguments give the same guarantee.

**Migrating the module library behind a shim.** Rejected: it was written against a base that no
longer exists, a `System` that is deleted, and a declared `duration` that is removed, and a shim
would keep asserting all three.

## Final interpretation

Read the references to `CompiledSequence` and its methods above as historical context. The current
forms are `sc.compile(tree, opts) -> pypulseq.Sequence` and
`sc.pns(tree, opts, hardware)`. [ADR-004](004-compile-returns-a-sequence.md) is authoritative for
the return type and diagnostic delivery; this ADR remains authoritative for scanner and module
ownership.
