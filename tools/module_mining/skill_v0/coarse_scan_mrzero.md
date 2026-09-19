# Coarse scan including MRzero-Core

**Corpus added:** `MRsources/MRzero-Core` @ `5732d86`, `documentation/playground_mr0/`.
**Registered as:** `discovery` + `oracle`, **explicitly not `design-witness`** — see §2.
**Date:** 2026-09-19. **No extraction performed, no Module proposed for implementation.**

This scan re-runs the coarse pass over the registry in `tools/module_mining/sources.yaml` with the
new corpus included. It reports families, duplicates and likely stop cases. It does not decide
anything.

---

## 1. What is in the corpus

26 notebooks; 18 build a sequence. A defined subset of 12 is exercised twice by the repository's
own CI — once for execution (`playground_test.yml`, one job per notebook) and once for numerical
regression (`tests/simulation_test/`, magnitude and phase NRMSE against stored reference files at
a tolerance of 0.01). That is a **maintained, executable** corpus, which is rarer than it sounds
and is the reason it is worth registering.

```text
FID                     mr0_FID_seq
SE / CPMG               mr0_SE_CPMG_seq
STE, 3 pulses 5 echoes  mr0_STE_3pulses_5echoes_seq
FLASH / GRE             mr0_FLASH_2D_seq, mr0_GRE_to_FLASH, mr00_FLASH_2D_ernstAngle_opt
EPI                     mr0_EPI_2D_seq
DWI SE-EPI              mr0_DWI_SE_EPI
diffusion-prep STEAM    mr0_diffusion_prep_STEAM_2D_seq
RARE / TSE              mr0_RARE_2D_seq, mr0_RARE_2D_seq_multi_shot, mr0_TSE_2D_multi_shot_seq
bSSFP                   mr0_bSSFP_2D_seq
DREAM                   mr0_DREAM_STE_seq, mr0_DREAM_STID_seq
IR-FLASH / T1 fitting   mr0_opt_FLASH_2D_IR_Fit_T1, mr0_opt_FLASH_2D_IR_voxelNN_T1
pTx / RF shim           pulseq_sim_pTx, pulseq_rf_shim  (out of SeqCraft's scope)
```

## 2. Independence — the finding that sets this corpus's role

**The playground is not independent of PyPulseq for sequence construction.** 21 of 26 notebooks
import `pypulseq`; all 18 sequence-building `mr0_*` notebooks do; and the test workflow pins
`pypulseq==1.4.2post1`. The sequences are **built with PyPulseq and simulated with MRzeroCore**.

So a playground notebook agreeing with a PyPulseq example is one design agreeing with itself, and
citing both for a candidate would over-count witnesses in exactly the way
`REFERENCE_NOT_INDEPENDENT` exists to prevent.

Where MRzero *is* independent is the other direction: the Bloch simulation is a measurement path
developed separately from Pulseq. Hence `oracle`, and hence `design-witness` explicitly excluded.

**Consequence worth flagging.** SeqCraft already depends on `MRzeroCore>=1.0` as the optional
`sim` extra, and it is the oracle behind every `02_simulate_and_reconstruct.ipynb`. So the corpus
being added for discovery is already the Layer-3 instrument. That is not circular for discovery —
finding a family and measuring a candidate are different jobs — but it does mean the simulation
layer has **one** oracle, and the Radial Layer-3 reasoning (a contract with one consumer, an
instrument with one implementation) applies to it too. Not actionable here; recorded.

## 3. Licence

**AGPL-3.0, plus a separate `EULA.txt`** covering the MRtwin code. This is the most restrictive
source in the registry — the other copyleft entries are GPL/AGPL without an additional agreement.

- Extract **physics, invariants and the existence of families only**.
- No implementation text may be copied or adapted into SeqCraft.
- SeqCraft's existing dependency is an optional extra, imported by examples, with nothing in
  `src/seqcraft` importing it. **Whether the AGPL + EULA combination constrains that usage is a
  licensing question, not a mining question**, and it is raised here rather than answered.

## 4. Families against current SeqCraft

Classified with the coverage vocabulary in `reference/outcomes.md`, because **"SeqCraft can
already build it" is not "SeqCraft already ships the right reusable abstraction"** and the first
draft of this scan collapsed both into the word "duplicate".

| family | coverage | note |
|---|---|---|
| FID | `PRIMITIVE_COMPOSITION_COVERED` | §4.1 |
| SE / CPMG | `DEGENERATE_CASE_OF_EXISTING_MODULE` + `ARCHITECTURE_REVISIT_CANDIDATE` | §4.2 |
| FLASH / spoiled GRE | `DIRECT_SHIPPED_DUPLICATE` | `GRE2DTR`, `GRE2D` |
| EPI | `DIRECT_SHIPPED_DUPLICATE` | `EPI2D` |
| RARE / TSE, multi-shot | `DIRECT_SHIPPED_DUPLICATE` | `TSEShot` + `FSE2D`; multi-shot is the segment table, caller policy |
| IR-FLASH / MPRAGE-like | `COMPOSITION_COVERED` (+ `NOTEBOOK_ONLY_EXISTING`) | §4.3 |
| DWI SE-EPI | `DEFER` | later architecture stress test; Wave 3, unchanged |
| pTx / RF shimming | `OUT_OF_SCOPE` | parallel transmit is not in SeqCraft's model |
| **bSSFP** | **`CANDIDATE`** | §5 |
| **stimulated echo / STEAM** | **`CANDIDATE — BOUNDARY OPEN`** | §6 |
| DREAM | **`YELLOW / DEFER`** | §7 |

Only three families are `DIRECT_SHIPPED_DUPLICATE`. The earlier "six duplicates" count was wrong
in a way that mattered: it treated "we can build this" as "we ship this", which is exactly the
conflation that hides missing abstractions and discoverability gaps.

### 4.1 FID — `PRIMITIVE_COMPOSITION_COVERED`

```text
Excitation -> delay / optional gradient -> ADC
```

No cross-component physical solve for a class to own. This is the canonical instance of the
classification, and the reason the classification exists: **a sequence can be fundamental and
common without deserving its own Module.** Without that category, a batch scan turns every
textbook sequence name into a class.

It is also a useful future calibration case — a prospective run on FID should reach
`NO_NEW_MODULE` or `RED - WRAPPER_ONLY`, and neither of those outcomes has ever been exercised.

### 4.2 SE / CPMG — a degenerate case, **not** a shipped duplicate

**SeqCraft ships no `SE2D`.** `sc.modules.__all__` is `CartesianLine`, `EPI2D`, `Excitation`,
`FSE2D`, `GRE2D`, `GRE2DTR`, `GRE3DTR`, `IRPrep`, `PhaseEncode`, `RadialReadout`, `Refocusing`,
`TSEShot`, `spoiler`. The only `SE2D` in the repository is notebook-local, in
`examples/se_2d/01_build.ipynb`.

What the package does have is `FSE2D(echoes=1)`, which is a conventional Cartesian 2D spin echo.
So:

```text
DEGENERATE_CASE_OF_EXISTING_MODULE   via FSE2D(echoes=1)
+ ARCHITECTURE_REVISIT_CANDIDATE
```

**A correction on the supporting evidence, because the claim scope is narrower than it looks.**
The proven event-for-event identity is `TSEShot` + `FSE2D` against **`examples/fse_2d/01`'s own
`FSE2D`**, at echoes 1, 8, 16 and 72 — `tests/modules/test_fse_notebook_matches_the_package.py`.
That is identity with the *FSE* notebook at `echoes=1`. It is **not** a measured comparison
against `examples/se_2d/01`'s `SE2D`, which is a separate class built at a different protocol and
has never been compared. The classification stands on the physics; the "event-for-event"
qualifier belongs to the `fse_2d` comparison and should not be transferred to `se_2d` without
running it. Running it would be cheap and is recorded here as an open item, not done.

**No `SE2D` should be added now.** It would add a name and some discoverability and **no new
physical contract**, which is `RED - WRAPPER_ONLY` by the skill's own test.

### 4.3 IR-FLASH / MPRAGE-like — composition, with a notebook class

`IRPrep` + `GRE2D` expresses it, and no equivalent public abstraction ships — so
`COMPOSITION_COVERED`. It is *also* `NOTEBOOK_ONLY_EXISTING`: `examples/mprage_2d/01_build.ipynb`
defines `MPRAGE2D` and `examples/mp2rage_2d/01_build.ipynb` defines `MP2RAGE2D`. Both
classifications are true and they answer different questions — the first says the physics is
covered, the second says where the convenience lives.

Unchanged unless a separate **timing contract** is later proven — an inversion-time placement
constraint that neither `IRPrep` nor `GRE2D` can determine alone. The MRzero notebooks here are
T1-fitting studies and supply no such evidence.

## 5. bSSFP — the strong candidate

`mr0_bSSFP_2D_seq.ipynb` builds a balanced SSFP repetition, and it has physics no current module
owns:

- **full gradient balance per TR on all three axes** — the slice rephaser `gzr1` is played before
  *and* after the readout, and the phase encode is rewound with `-phenc`; net zero moment each TR;
- **TR symmetric about the echo** — `minTR2 = minTR/2`, with the fill `TRd` split either side, so
  the echo sits at the centre of the repetition;
- **180° alternating RF phase cycling**, with the ADC phase following;
- **α/2 preparation at TR/2**, optional, to enter the steady state without oscillation.

None of these is `GRE2DTR` with spoiling disabled. Balance is a constraint across the whole
repetition that `GRE2DTR` does not express, and the symmetric-TR requirement couples the fill
delays to the echo position. Provisional layer: `kernel`, sibling of `GRE2DTR` and `GRE3DTR`.

**The API to avoid**, if the physical contract is genuinely different:

```python
GRE2DTR(balanced=True, spoiled=False, phase_cycle=...)
```

Flags that switch a module between physical contracts produce one oversized class whose invariants
depend on its arguments — and they are how application-specific branches (diffusion-sensitive
bSSFP, CEST with a bSSFP readout) end up inside a GRE class. A clean kernel boundary keeps those
as compositions.

Whether a `BSSFP2D` imaging module follows is a **separate** decision, made on usage evidence
after the kernel, exactly as `GRE3D` was declined for `GRE3DTR`.

**This candidate is not new.** The 2026-09-18 coarse scan already reached it from
`pulseq/pulseq` `writeTrufi.m` and OpenMRF's BSSFP readout, named it **`BalancedSSFPTR:
NEW_KERNEL`**, listed the same distinctive physics, and concluded: *"Do not add a `balanced=True`
option to `GRE2DTR`; the invariants are different enough to justify a separate kernel."*

So MRzero **corroborates an existing candidate** rather than producing one, and it does so from a
source that is not an independent design witness — the playground builds with PyPulseq. The
evidence that would carry a bSSFP fine scan is `writeTrufi.m` plus OpenMRF, both already in the
registry and one of them genuinely independent. MRzero's contribution here is a maintained,
executing, simulation-checked instance of the family, which is worth having as discovery and is
worth nothing as corroboration.

## 6. Stimulated echo / STEAM — the ambiguous candidate

**This is the recommended next supervised candidate**, and it was not chosen because it looks
likely to fail. It is chosen because the corpus contains the same three-pulse structure three
times for three different purposes, and the question of whether that is one contract or three uses
is genuinely open:

```text
mr0_STE_3pulses_5echoes    90 / 90 / 180, no preparation intent -- an observation of every
                           coherence pathway three pulses produce (5 echoes from 3 pulses)

mr0_diffusion_prep_STEAM   90 / 180(phi=90) / 90(phi=180) storage, then a 5 deg readout train
                           -- magnetisation tipped, dephased, STORED along z, then read out

mr0_DREAM_STE              two small-flip STEAM pulses (second at phi=180), then a sinc readout
                           train -- reads STE and FID together to compute a flip-angle map
```

Shared, and real: two pulses separated by τ, a mixing time TM during which the magnetisation is
stored longitudinally, and a dephasing/rephasing balance that must hold across TM. The fraction
stored and the gradient bookkeeping either side of TM are arithmetic somebody can get wrong.

Not shared, and the reason this is ambiguous: the middle pulse differs in flip *and* role across
the three; the first notebook is not a preparation at all; and DREAM's value comes from reading
**two** pathways rather than from preparing one.

The fine scan's job is to decide where the stable reusable unit sits — something nearer
`StimulatedEcho`, `STEAMKernel` or `STEAMPrep` — and, before that, to separate:

```text
shared coherence-pathway physics
    vs
application-specific timing / diffusion weighting / mapping policy
```

Plausible outcomes, all live, none predetermined:

```text
NEW_LEAF                  a preparation leaf, sibling of IRPrep
NEW_KERNEL                if the storage and the readout must be solved together
EXTEND_EXISTING           IRPrep generalises to a stored-magnetisation preparation
YELLOW - PHYSICAL_BOUNDARY_UNCLEAR   the three uses do not share one owner
NO_NEW_MODULE             three pulses and two delays is composition
RED - WRAPPER_ONLY        if extraction turns out to own no arithmetic
```

A candidate that can reach five different outcomes is worth more to a skill calibration than one
that can reach one. Note also that a STEAM preparation touches `IRPrep` — a shared leaf — so
**rule D fires**, which is the other thing worth exercising.

## 7. DREAM — `YELLOW / DEFER`, not a pre-judged stop case

**Reclassified.** The first draft of this scan called DREAM "likely `NO_NEW_MODULE`". That was a
coarse-scan verdict on a question the coarse scan cannot answer, and it is withdrawn.

DREAM is a **B1 / B0 / transmit-receive mapping** method built on stimulated-echo preparation and
a readout path. The architectural question is not *is DREAM important enough to deserve a class* —
that question is about the name, not the physics. It is:

> After the shared stimulated-echo / STEAM physics is extracted, does DREAM still contain a
> stable reusable physical contract of its own?

That cannot be answered before §6 is. If the answer decomposes cleanly into

```text
StimulatedEcho / STEAM  +  GRE-style readout  +  mapping-specific ordering and post-processing
```

then `NO_NEW_MODULE` is likely right. But if DREAM still owns irreducible cross-component timing
or echo-placement behaviour **reused across its variants** — the corpus has two,
`mr0_DREAM_STE_seq` and `mr0_DREAM_STID_seq`, which differ in which coherence pathway is read —
then a higher-level kernel may be justified.

```text
status: YELLOW / DEFER
revisit: after the StimulatedEcho / STEAM boundary is understood
```

### 7.1 BOLD is a different question and should not be grouped with DREAM

Worth stating because the two get filed together as "functional" methods and they are not alike.

```text
BOLD      a contrast / application, commonly acquired with GRE-EPI
DREAM     a mapping method built on stimulated-echo preparation
```

BOLD acquisition is expected to **compose existing physics** — a GRE or SE contrast mechanism, an
EPI readout, and time-series scheduling — rather than to need a `BOLD` module by name. Time-series
scheduling is acquisition policy, which is caller-owned. Nothing in this corpus changes that, and
no `BOLD` candidate is proposed.

## 7.2 A dependency note, recorded before anyone is tempted

Two structures this scan explicitly does **not** propose, kept here so the reasoning is on record
when the question returns.

**Do not make `FSE2D` depend on a future `SE2D`.** FSE is not a complete SE scan repeated many
times — the echo train's crusher window, echo spacing and midpoint placement are properties of the
train, not of a single spin echo iterated. A dependency in that direction would put the simple
case above the general one and force the train to reach inside decisions a single-echo class had
already made. This is the same argument that made `GRE3DTR` a **sibling** of `GRE2DTR` rather than
a wrapper around it.

If a reusable spin-echo contract exists, it is **lower** than either:

```text
            shared spin-echo physics
                (SpinEchoCore?)
                 /            \
          single SE          TSEShot
              |                 |
            SE2D              FSE2D
```

**`SpinEchoCore` is not proposed and must not be created yet.** The current evidence supports only
"it may exist". What would turn that into real evidence is comparing Cartesian SE, SE-EPI, DWI
SE-EPI and any other refocused readout, and finding they independently share

```text
excitation
180 degree refocusing
echo-centre semantics
TE symmetry
crusher / coherence handling
```

**while differing only in readout.** SeqCraft already has two of those four — `examples/se_2d` and
`examples/se_epi_2d` — and the corpus adds DWI SE-EPI, so the comparison is affordable. Until it
has been run, `FSE2D` and `TSEShot` stay exactly as they are.

## 8. Recommendation

```text
1. stimulated echo / STEAM   CANDIDATE - BOUNDARY OPEN; recommended first prospective run
2. bSSFP                     strong NEW_KERNEL candidate (BalancedSSFPTR), already identified
3. DREAM                     YELLOW / DEFER -- revisit after 1
4. spin-echo core comparison  evidence-gathering only; SE / SE-EPI / DWI SE-EPI
5. DWI / diffusion prep      unchanged: later architecture stress test
```

`SaturationPrep` and `Spiral` keep their existing roadmap positions. This scan inserts a
calibration run before them; it does not replace or reorder them.

**Why STEAM rather than bSSFP, on information gain rather than chance of GREEN.**

bSSFP is the safer candidate and that is the argument against it going first. It was already
identified in the 2026-09-18 scan as `BalancedSSFPTR: NEW_KERNEL`, its physics is well separated
from spoiled GRE, and it has independent evidence in `writeTrufi.m` and OpenMRF. A prospective run
would very likely reach GREEN — and would exercise the workflow on a case whose answer is already
believed. A fourth GREEN tells us little; every stop outcome would remain untested, which is the
gap `generalization_report.md` §4 names.

STEAM's outcome is genuinely unknown, and each possible result teaches something different:

| if it lands | what it exercises |
|---|---|
| `NEW_LEAF` / `NEW_KERNEL` | a preparation-layer candidate, a layer no pilot has covered |
| `EXTEND_EXISTING` | a status never exercised, and rule D against `IRPrep` |
| `YELLOW - PHYSICAL_BOUNDARY_UNCLEAR` | the first YELLOW the workflow has ever produced prospectively |
| `NO_NEW_MODULE` | the stop outcome the whole design rests on, never once reached |
| `RED - WRAPPER_ONLY` | the other untested stop |

It also puts pressure on parts of the skill nothing else has: the evidence is under **AGPL plus a
EULA**, so the licence and independence machinery is load-bearing rather than decorative; the three
uses come from one non-independent corpus, so `REFERENCE_NOT_INDEPENDENT` is live; and a
preparation candidate touches `IRPrep`, so rule D fires on a shared leaf.

The honest risk: STEAM may stall in a way that teaches us about MR rather than about the workflow,
and the AGPL corpus may not supply usable acceptance evidence at all — in which case the
Pulseq/OpenMRF STEAM and DREAM sources must carry it. If the fine scan cannot assemble acceptance
evidence from permissively licensed sources, that is itself a result worth having, and it argues
for bSSFP as the next run instead.

**No extraction is proposed here.** The next step, when authorised, is a prospective supervised
fine scan: reference evidence and acceptance claim first, outcome unknown at the start, and the
stop outcomes genuinely available.
