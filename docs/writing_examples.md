# Writing example notebooks

**This is the authoring guide.** For the index of what examples exist, what each one teaches and
what they need installed, see [`examples/README.md`](../examples/README.md).

The audience here is whoever writes or edits a notebook under `examples/` -- contributor or agent.
It is the single source of truth for that; `AGENTS.md`, `CLAUDE.md` and the module-mining Skill
point here rather than restating any of it.

Its companion is [`writing_a_module.md`](writing_a_module.md), which covers reusable components
and the prose that documents them. The boundary the two of them are enforcing:

```text
tools/module_mining/     evidence, abstraction choice, validation argument
docs/writing_examples.md tutorials -- teach the MRI sequence or experiment
docs/writing_a_module.md components -- and public API/docstring prose
```

**Examples are tutorials and runnable demonstrations, not pull-request records, design reviews or
development postmortems.** A reader arrives knowing some MRI and no project history, wanting to
build something. Write for them.

The shape that works:

```text
concept  ->  physics  ->  SeqCraft API  ->  composition  ->  measurement
         ->  interpretation  ->  limitations  ->  references
```

**Sequence first, subtlety second.** An opening should normally establish, in this order: what
sequence or acquisition family this is; its basic physical structure — RF, gradients, encoding,
readout; what it is used to measure; and only then the particular question this notebook explores.
A reader who has seen only the **title and the first two paragraphs** should be able to answer the
first three. Keep the interesting sentence — put it where it lands on a reader who now knows what
it refers to, as a section heading or a conclusion.

**Do not assume the reader has completed a preceding example.** Cross-links may deepen context,
but an ordinary build or simulation notebook should establish its sequence or acquisition family,
its physical purpose and its main components **on its own**. A cross-link extends an explanation;
it does not supply the missing introduction.

Both rules are about ordinary **build** and **simulation/reconstruction** notebooks.
Special-purpose documents — API guides, preserved reference implementations, migration guides,
explicit validation notebooks — may lead with their document purpose instead, and should say what
that purpose is in the first paragraph. [`fse_2d/03_module_api.ipynb`](../examples/fse_2d/03_module_api.ipynb)
is the example: it opens by relating itself to the reference notebook beside it, which is exactly
what it is for.

| | |
|---|---|
| **start from the MRI concept and the user's goal** | title the notebook after the sequence or the problem, not after a conclusion about it. `# Off-resonance in spiral imaging`, not `# A spiral pays for its long readout in off-resonance` |
| **define domain terms before using them** | if b-value, VENC or turbo factor is the organising idea, say what it is first |
| **explain formulas directly** | write the equation and define every symbol; a citation supports the explanation, it is not the explanation |
| **keep the depth** | measurements, closed forms, tolerances and honest negative results are the value. Only the framing changes |
| **state limitations plainly** | what the example demonstrates, and what it does not |
| **end with `Summary` and `References`** | not "What this notebook established" |

Do **not** organise a notebook around:

```text
abstraction-selection history      why this is a kernel and not a leaf
rejected alternatives              what a different design would have cost
previous implementation bugs       what an earlier version of this module got wrong
reference-repository shortcomings  what some other project leaves out
licence and evidence arguments     which corpus witnessed what
solver / compiler avoidance        that no new layer was needed
```

All of that is real and worth recording — in `CHANGELOG.md`, in the pull request, or in
`tools/module_mining/plans/<candidate>/`. Page-level source provenance belongs in
`tools/module_mining/domain_evidence/`.

**Explain the MRI result, not why the check is trustworthy.** Measurements and controls belong in
the notebook; acceptance-criterion arguments, validator independence, why a module exists, what a
test would catch, and what evidence was sufficient to ship belong in the candidate/findings/PR. A
simulator or tool limitation belongs in the notebook only when the reader needs it to interpret or
reproduce the result.

> **The reader test.** If a sentence would not help a reader who knows MRI but knows nothing about
> SeqCraft's development history, move it out of the tutorial.

That is a test of the sentence, not a banned-word list — but these phrasings are reliable enough
signals to be worth stopping on when they appear:

```text
the claim / the contract              a physical statement, written as a promise the code makes
this module exists to                 motivation for the abstraction rather than for the physics
the witness / the corpus              evidence provenance
this establishes / proves             acceptance-criterion language
why this is a leaf / a kernel         layering rationale
the analyser was checked too          validator independence
measured rather than restated         an argument about the measurement instead of its result
```

Each can be rewritten into the physics it is standing in front of. "This is the claim the module
exists to make good on" is usually one sentence away from "the two toggles differ in first moment
by `delta_m1`"; the second is what the reader needed.

Avoid raising an imagined objection in order to answer it — "why this is not a defect", "the
check this notebook exists for", "three pieces, none of them new". State what the sequence does
and what the reader is about to build. The exception is a genuine **MRI** misconception, which is
worth confronting directly: `se_epi_2d/02` exists because most readers expect a spin echo to fix
EPI distortion, and it does not.

**Deliberately wrong cases earn their place when they teach physics or protocol design** — a
mismatched echo time that contaminates an ADC, a readout aligned to the wrong instant. Replaying a
historical *software* bug because it was instructive during development does not.

---

## What a simulation notebook is allowed to cost

MRzero's inner loop is dense complex linear algebra and it takes every core it is offered — sixteen
threads by default on a 32-core machine. A draft of [`fse_2d/02`](../examples/fse_2d/02_simulate_and_reconstruct.ipynb)
ran several 16-echo, 9-shot acquisitions in a single cell, minutes of saturated CPU with no output
until it finished, and took a workstation down.

So the spin-echo and EPI `02` notebooks each set `SIM_THREADS` **before importing torch**, print
their budget in the first cell, and keep anything expensive behind a named switch that is off by
default. `se_2d/02` is 38 s. `fse_2d/02` is 15 s of measurement plus 4 minutes of images, and the
images are the part behind the switch — every *number* in it is measured on a single spin, because a
ghost is a modulation of `ky` and a point object's k-space is that modulation.

`megre_2d/02` is three acquisitions and 24 SENSE reconstructions, and it has no switch because it
does not need one: it **picks up a GPU when there is one**, and the MRzero simulation — which is
essentially all of its cost — is a torch computation that moves there wholesale. Both `DEVICE` and
sigpy's `RECON_DEVICE` are chosen by *probing* rather than by asking, because
`torch.cuda.is_available()` reports the driver and comes back true on a machine whose devices are
all hidden, and because cupy needs a toolkit to compile against where torch needs only a driver to
talk to — so the two can disagree, and a machine can simulate on its GPU with no working cupy at
all. Either way the notebook is the same notebook, and every step prints what it cost.

The two EPI `02`s take the same bargain at a different scale: about a minute of one-voxel
measurement each, then `BRAIN_IMAGES` for the slab. `gre_epi_2d/02`'s switch also covers §7, which
is the one section that **cannot** be measured on a single spin — parallel imaging is a statement
about receive sensitivity, so it needs coils and an object or it needs nothing.

If you add a simulation to an example, price it first: the cost is voxels × repetitions × states, and
a 16-echo train has all three. Two of the three are not negotiable — the repetitions are the
sequence, and dropping `max_state_count` from 200 to 32 is 5 % wrong because a CPMG train's
stimulated pathways *are* the physics. The slab thickness was the one that turned out to be free, and
only because it was measured rather than assumed.

---

## Preserved reference notebooks may disagree with the package

A notebook that a shipped module was extracted from is sometimes kept as its **own** implementation
rather than replaced with an import -- [`fse_2d/01_build.ipynb`](../examples/fse_2d/01_build.ipynb)
is the case, with a test comparing the two event for event. That is deliberate: a reference which
is rewritten whenever the package changes cannot detect a regression in the package, so it has to
stay **able to disagree**. If you are editing such a notebook, do not "fix" it by importing the
module it is checking.
