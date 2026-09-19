# The four rules, the failures behind them, and the compiler escalation

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

This will be live for Spiral: PR #23 proposes changes to m0/m1 verification tolerances. Treat that
PR as an evidence corpus, not a specification, and run the escalation from the top.
