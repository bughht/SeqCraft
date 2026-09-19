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

| family | status | note |
|---|---|---|
| FID | `COVERED` | `Excitation` + an ADC. Pedagogical. |
| SE / CPMG | `COVERED` | `Refocusing`, `TSEShot`. |
| FLASH / GRE | `COVERED` | `GRE2DTR`, `GRE2D`. |
| EPI | `COVERED` | `EPI2D`. |
| RARE / TSE, multi-shot | `COVERED` | `TSEShot` + `FSE2D`; multi-shot is the segment table, which is caller policy. |
| IR-FLASH | `COVERED` | `IRPrep` + `GRE2D`; the notebooks are T1-fitting studies, not new physics. |
| DWI SE-EPI | `DEFER` | Diffusion is Wave 3 in the roadmap and stays there. |
| pTx / RF shimming | `OUT_OF_SCOPE` | Parallel transmit is not in SeqCraft's model. |
| **bSSFP** | **`CANDIDATE`** | §5 |
| **stimulated echo / STEAM** | **`CANDIDATE — AMBIGUOUS`** | §6 |
| DREAM | `LIKELY NO_NEW_MODULE` | §7 |

**Six of eleven families are duplicates of shipped modules.** That is the expected result for a
teaching corpus and is itself useful: it means the corpus's value here is the *tail*, not the head.

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

Plausible outcomes, all live:

```text
NEW_LEAF / preparation    a STEAMPrep sibling of IRPrep
EXTEND_EXISTING           IRPrep generalises to a stored-magnetisation preparation
YELLOW - PHYSICAL_BOUNDARY_UNCLEAR   the three uses do not share one owner
NO_NEW_MODULE             three pulses and two delays is composition
RED - WRAPPER_ONLY        if extraction turns out to own no arithmetic
```

A candidate that can reach five different outcomes is worth more to a skill calibration than one
that can reach one. Note also that a STEAM preparation touches `IRPrep` — a shared leaf — so
**rule D fires**, which is the other thing worth exercising.

## 7. DREAM — likely stop case, but not the calibration target

DREAM looks like a protocol rather than a physical contract: a STEAM pair, then a GRE readout
train, then a ratio of two measured signals. The ratio is reconstruction, the readout is
`GRE2DTR`, and the STEAM pair is §6's question. If §6 yields a preparation module, DREAM is
composition on top of it; if §6 yields nothing, DREAM has nothing to sit on.

It is listed as a likely `NO_NEW_MODULE` and is **not** proposed as the calibration candidate —
picking a family already judged likely to stop would be the hand-picking the brief rules out.

## 8. Recommendation

```text
1. bSSFP                strong CANDIDATE, kernel, evidence available without MRzero
2. stimulated echo      CANDIDATE - AMBIGUOUS, recommended next supervised run
3. DREAM                revisit only after 2
4. DWI SE-EPI           unchanged: Wave 3
```

`SaturationPrep` and `Spiral` stay in the roadmap in their existing positions. This scan adds a
calibration run before them; it does not replace or reorder them.

**No extraction is proposed.** The next step, when authorised, is the supervised fine scan for the
stimulated-echo candidate, run prospectively — reference evidence and acceptance claim first,
outcome unknown at the start.
