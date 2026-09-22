# The rules, the failures behind them, and the compiler escalation

Each rule is here because something passed every check that existed and was still wrong. The
failure is kept with the rule, because a rule without its failure gets optimised away by the next
person who finds it inconvenient.

---

## A. Physical mode contract

### The failure

`GRE3DTR`'s non-selective mode inherited `Excitation`'s 3 ms shaped sinc and dropped the selection
gradient. The reasoning was "mode B is mode A minus Gz." The result:

- internally consistent;
- compiled legally;
- satisfied every k-space and timing invariant;
- passed the whole test suite;
- spent a soft pulse's duration selecting nothing.

It put 3 ms into every echo time for the privilege. Minimum TE went 4.377 -> **2.977 ms** when
fixed. No invariant in the comparator looks at pulse shape, because none had a reason to.

### The rule

Before implementing a candidate with meaningful modes or defaults, publish a table with one row
per mode and one column per axis of physical difference:

```text
RF family / shape
RF duration and source
selection-gradient presence/absence
intrinsic rephasing / moment requirement
encoding role
ADC / readout semantics
timing consequences
reference evidence
explicitly user-overridable choices
```

Then the differential, which is what catches "minus one event":

```text
mode A -> mode B

changes:
    ...

must remain invariant:
    ...
```

> **A mode may not be implemented as "another mode minus one event" unless reference evidence
> establishes that equivalence.**

Every official Pulseq 3D reference uses a block pulse for non-selective excitation. One look at
the evidence column would have settled it before any code existed. Worked example:
`tools/module_mining/plans/gre3d/mode_contract.md`.

---

## B. Emitted-sequence inspection

### The rule

At least one representative emitted `.seq` per promoted physical mode, read back without
reference to the Python that wrote it:

```bash
python tools/module_mining/inspect_emitted.py examples/gre_3d/seq/*.seq
```

Report, per mode: RF waveform/family, RF duration, selection gradients, rephaser realisation,
readout/ADC, and spoiler/rewinder where relevant.

> If the emitted `.seq` were inspected without reading the Python, would it tell the same physical
> story as the API name and the documentation?

This is a **review aid, not a gate** — it prints, it does not assert. Its value is that it reads
the artifact the scanner will actually play, through a different path than the one that built it.
A derived measurement computed from the same objects the module constructed cannot catch a wrong
pulse family; a file read back can.

---

## C. Claim-scoped reference validation

### The rule

Every reference comparison states both halves, adjacent, in the report and anywhere the result is
summarised:

```text
this comparison establishes:
    ...

this comparison does NOT establish:
    ...
```

"Validated against X" may never stand in for more than was measured.

### Why it needs enforcing

The gap is rarely a lie and usually a slide. `GRE3D`'s comparison against `gre3d.seq` establishes
that Δk agrees on three axes and that every probed line and partition lands where asked. It does
**not** establish that the two sequences are the same sequence, that the slab-selective path is
validated at all (no executable reference exists for it), or that any image has been reconstructed.

The same rule applied to a status table is what produced the layered `RadialReadout` claim: GREEN
for reference validation and acquisition build, and explicitly **not** "SeqCraft has an
independently validated end-to-end radial imaging and non-Cartesian reconstruction pipeline."

---

## D. Shared-leaf dependency impact

### The rule

A candidate that changes an existing leaf maps the blast radius before the change lands:

```text
changed leaf
    -> direct consumers
    -> transitive consumers
    -> examples / notebooks
```

and classifies every affected consumer as exactly one of:

```text
NO_BEHAVIOR_CHANGE
NEW_EARLY_REFUSAL_FOR_PREVIOUSLY_INVALID_INPUT
INTENTIONAL_BEHAVIOR_CHANGE
UNKNOWN / NEEDS_REVIEW
```

`UNKNOWN / NEEDS_REVIEW` is a legitimate answer and blocks promotion until resolved; a guess
dressed as `NO_BEHAVIOR_CHANGE` does not.

### The trap

**Grep for composition, not for the class name.** `CartesianLine` looked like it had five
consumers. `EPI2D` mentions the class only in a docstring at the top of the file and composes
nothing from it — the real consumers are `gre_2d_tr.py`, `gre_3d_tr.py`, `tse_shot.py` and
`radial_readout.py`. A name in prose is not a dependency, and treating it as one inflates the
review while hiding that nobody checked.

The distinction matters most for the second classification: when `CartesianLine` began refusing an
ADC sample count that is not a multiple of `adc_samples_divisor`, every consumer needed to be
checked for whether it had been *relying* on passing such a count through to the compiler. That is
`NEW_EARLY_REFUSAL_FOR_PREVIOUSLY_INVALID_INPUT` — new refusal, no previously-valid behaviour
changed — and it is only safe to write down after looking.

---

## E. Measure the module on the lattice it emits

> **A module that checks its output on a more convenient representation than the one it emits
> will certify something it cannot build.**

Derived quantities computed from a module's *internal* arrays are not the module's output. The
output is the event the factory accepts and the compiler plays, and its representation may differ
in ways that are invisible until something refuses it: sample positions, endpoint conventions,
interval widths at the edges, quantisation.

The check is mechanical: whatever property a module asserts about itself -- a limit respected, an
area delivered, a position reached -- compute it from the **emitted event's own knots**, and
prefer an independent measurement of the compiled sequence where one exists.

This is why `limits()`-style accessors are measurements rather than restatements of a solver's
constraints, and why a reported trajectory is integrated from emitted knots rather than carried
over from a design pass.

---

## F. Validate temporal occupancy after events are placed together

> **Individually correct events can become physically wrong once placed: one event can extend the
> compiled block past another waveform's designed support.**

Every event may satisfy its own contract and the assembly still be wrong, because occupancy is a
property of the *arrangement*. An event that reserves more time than its nominal content -- a
guard, a dead time, a raster round-up -- can push a block boundary past where another waveform was
designed to end, and the compiler then holds that block open and pads what is inside it.

The padding is not a compiler error. It is the correct response to a tree that asked for it.

So: after placement, check that every waveform's designed support **contains** the block it lands
in, not merely that each event is individually legal. Budgets computed by subtracting nominal
overheads are the usual way this is got wrong, because the true overhead includes the rounding.

---

## G. Composition-level claims are validated on the composition

### The failure

Not a historical one — this rule comes from a fine scan rather than from a defect, which is why it
is stated narrowly.

The flow-compensation candidate's central invariant is:

```text
the n-th moment of EVERY gradient on one axis,
from a semantic origin to the echo,
equals the target
```

`Excitation`, `PhaseEncode`, `CartesianLine` and the compensation all contribute. **Every one of
them can be individually correct and the sum wrong**, and no module-local check can see it. Rules
B, E and F are all module-local: inspect what *you* emit, measure on the lattice *you* emit,
validate occupancy after *your* events are placed.

### The rule

> When the property a candidate claims is a property of several modules together, validate it on
> the composition, not on the parts.

And:

> **The validator integrates; it does not ask.**

A module reporting its own `M1` while the validator checks that number against the module's own
target tests the module's arithmetic against itself. The derived quantity is computed from the
emitted events by something that shares no code with the module — which is the same argument rule
E makes about lattices, applied to quantities that span modules.

### What this is not

It is not a requirement that every candidate validate on a composition. Most claims are local and
rules B, E and F reach them. This applies when the *claim* spans modules — and a candidate whose
claim spans modules is usually a joint-realization kernel, which is the other thing the post-v0
scans found.

### A note on derived quantities

A moment is not an intrinsic property of a gradient event: it depends on the placement time, the
semantic origin, the integration endpoint and the order. So it is not something an event carries
or a leaf publishes — it is computed on demand from events and their placements. Keep moment
analysis an internal helper, and if it ever becomes a bottleneck, cache primitive integrals inside
that helper rather than changing the event API.

---

## Compiler-change escalation

The compiler owns Pulseq legality and knows nothing about any candidate. A candidate that will not
compile is, until proven otherwise, a candidate that is wrong.

```text
candidate exposes apparent compiler problem
        |
        v
remove candidate implementation from reproducer
        |
        v
minimal raw PyPulseq / LogicBlock case
        |
        v
prove a legal physical input is rejected or altered
        |
        v
characterise the numerical / error floor independently
        |
        v
only then consider a compiler change
```

Each step exists to defeat a specific way of being wrong:

| step | what it rules out |
|---|---|
| remove the candidate | the bug is in the candidate |
| minimal raw case | the bug is in the composition, not the compiler |
| prove a **legal physical** input is affected | the input was genuinely illegal and the refusal correct |
| characterise the error floor independently | the "failure" is a tolerance chosen to fit this candidate |

Without a candidate-free reproducer the outcome is
`YELLOW - COMPILER_CHANGE_CLAIM_NEEDS_MINIMAL_REPRODUCER`. A candidate that needs source-specific
compiler behaviour is `RED`.

The fourth step is the one most often skipped and it is the one that matters for tolerances: a
proposal to widen an m0 or m1 verification tolerance must characterise the numerical floor from
something other than the candidate that wants the tolerance widened. Otherwise the instrument is
being adjusted until it agrees with the thing it is supposed to be measuring.

This becomes live whenever a candidate arrives alongside an existing branch or PR that already
proposes compiler or tolerance changes. Treat such a branch as an **evidence corpus, not a
specification**, and run the escalation from the top regardless of how complete the proposed fix
looks. The roadmap records which candidates this is expected to apply to; this file does not.
