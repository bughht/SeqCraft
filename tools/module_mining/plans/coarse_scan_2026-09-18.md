# SeqCraft Open-Source Sequence Coarse Scan

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-18, so that the pull request extracting `TSEShot` and `FSE2D` carries the evidence
> that justified it.  The two copies are identical today and there is nothing keeping them
> that way; see [`../README.md`](../README.md) for which one to edit.

**Date:** 2026-09-18  
**Scope:** Official PyPulseq examples, official Pulseq MATLAB demos, selected research repositories, and the current SeqCraft module/example library.

## 1. Purpose

This is a **coarse architecture/mining scan**, not a conversion pass.

The goal is not to map every open-source sequence to a new SeqCraft class. The goal is to answer:

1. Which reusable MRI capabilities are already covered by `sc.modules`?
2. Which open-source sequences expose a real missing reusable primitive?
3. Which missing capability should extend an existing module instead of creating a sibling class?
4. Which complete sequence should remain composition-only?
5. Which notebook-local SeqCraft implementations now have enough external evidence to consider promotion?
6. Does the growing module set expose weaknesses in the current `rf / preparation / encoding / readout / kernel / imaging` taxonomy?

The central rule remains:

> A module is extracted from working sequence code; it is not designed in the abstract.

For this coarse scan, "candidate" means **a hypothesis to investigate in the fine scan**, not permission to implement it.

---

# 2. Corpus scanned

## 2.1 SeqCraft baseline

Current packaged modules:

| Layer | Existing capability |
|---|---|
| `rf/` | `Excitation`, `Refocusing` |
| `preparation/` | `IRPrep` |
| `encoding/` | `PhaseEncode` |
| `readout/` | `CartesianLine`, `EPI2D` |
| `kernel/` | `GRE2DTR` |
| `imaging/` | `GRE2D` |
| top-level function | `spoiler()` |

Important current capabilities that materially affect candidate mining:

- `Excitation` already supports sinc / SLR / Gaussian / block pulses and selective or non-selective use.
- `Refocusing` already handles crusher balance around the **effective RF centre**, not only nominal pulse midpoints.
- `IRPrep` already supports adiabatic inversion (`hypsec`, `wurst`) and an inversion crusher.
- `PhaseEncode` is axis-generic; a future 3D acquisition does **not** need a duplicate `PartitionEncode`.
- `CartesianLine` already supports:
  - exact prephaser arithmetic;
  - partial echo / readout partial Fourier;
  - `prephase=False` for spin-echo use;
  - multi-echo trains;
  - monopolar and bipolar echo trains;
  - explicit echo spacing;
  - semantic echo timing.
- `EPI2D` already encapsulates the difficult EPI readout mechanics, including alternating polarity, blips and advanced timing behaviour.
- `GRE2DTR` owns the repeating-unit timing and coupling of excitation, winders, readout, rewinder and spoiling.
- `GRE2D` owns acquisition ordering / line lists / dummy repetitions rather than baking sampling policy into lower modules.

This means the external corpus should be compared against a fairly capable base, not against raw PyPulseq.

## 2.2 SeqCraft notebook-local implementations

The current examples are important evidence because they already contain working implementations that have deliberately **not** all been promoted into the package:

| Example | Notebook-local composite |
|---|---|
| `examples/se_2d` | `SE2D` |
| `examples/fse_2d` | `FSE2D` |
| `examples/se_epi_2d` | `SEEPI2D` |
| `examples/gre_epi_2d` | `GREEPI2D` |
| `examples/mprage_2d` | `MPRAGE2D` |
| `examples/mp2rage_2d` | `MP2RAGE2D` |
| `examples/megre_2d` | no new class; uses existing GRE / Cartesian multi-echo capability |

These examples make it possible to distinguish:

- a genuinely missing physics primitive;
- a package-promotion question;
- a sequence that is already naturally expressible as composition.

## 2.3 Official PyPulseq

Scanned example families include:

- Cartesian GRE and labelled/soft-delay GRE;
- gradient-echo EPI and labelled EPI;
- spin-echo EPI and high-performance ramp-sampled SE-EPI;
- TSE;
- HASTE;
- 2D / 3D MPRAGE;
- radial GRE;
- UTE.

## 2.4 Official Pulseq MATLAB

The MATLAB demo corpus expands the coverage substantially:

- Cartesian GRE, 3D GRE, GRAPPA variants;
- multi-echo GRE;
- cine GRE;
- EPI / ramp-sampled EPI / spin-echo EPI;
- diffusion EPI;
- MPRAGE variants;
- TSE / HASTE;
- radial GRE;
- spiral;
- TrueFISP / bSSFP;
- UTE;
- ZTE / PETRA;
- PRESS;
- semi-LASER;
- selective RF and spin echo;
- labels, triggers and scanner-facing control examples.

## 2.5 Research repositories used as cross-validation

The coarse scan also sampled:

- `ritagnunes/PulseqDiffusion`
  - single-refocused DWI;
  - twice-refocused DWI;
  - b-value / direction handling.
- `usc-mrel/rtspiral_pypulseq`
  - 2D/3D real-time spiral;
  - explicit separation between spiral trajectory design, rewinders, RF kernels and sequence assembly.
- `BAMMri/Open4DFlow`
  - velocity encoding / first-moment gradient design;
  - cardiac triggering and 4D-flow acquisition.
- `OpenMRF/openmrf-core-matlab`
  - modular preparation library;
  - FAT, SAT, inversion, T2, spin-lock, CEST;
  - BSSFP, EPI, GRE, radial, spiral, TSE, UTE, PRESS and other readout families;
  - trajectory-design utilities separated from higher-level sequence assembly.
- `kherz/pulseq-cest`
  - CEST preparation protocols and reusable saturation-period concepts.

### License note

The scanner must record provenance and license before any conversion work.

At the time of this scan:

| Repository | Repository-level license signal |
|---|---|
| `pulseq/pypulseq` | MIT |
| `pulseq/pulseq` | MIT license file |
| `ritagnunes/PulseqDiffusion` | AGPL-3.0 |
| `usc-mrel/rtspiral_pypulseq` | GPL-3.0 |
| `BAMMri/Open4DFlow` | GPL-3.0 |
| `OpenMRF/openmrf-core-matlab` | MIT for original OpenMRF work, with third-party components covered separately |
| `kherz/pulseq-cest` | MIT |

The mining system should extract **physics and invariants**, not copy implementation text from incompatible codebases.

---

# 3. Status vocabulary for this scan

Every candidate is assigned one action:

| Status | Meaning |
|---|---|
| `COVERED` | Current packaged SeqCraft module already expresses the essential reusable physics. |
| `EXTEND_EXISTING` | A genuine missing capability exists, but it belongs naturally in an existing module. |
| `PROMOTE_NOTEBOOK` | SeqCraft already has a working notebook-local implementation; external evidence now strengthens the case for package promotion/refactoring. |
| `NEW_LEAF` | A new reusable physics/readout/preparation leaf is justified. |
| `NEW_KERNEL` | A new repeating acquisition unit is justified. |
| `NEW_IMAGING` | A complete-scan abstraction may be justified after its kernel is stable. |
| `COMPOSITION_ONLY` | The sequence is primarily an arrangement of existing reusable pieces and should not automatically become a packaged class. |
| `ARCHITECTURE_REVIEW` | The candidate exposes a taxonomy or non-module infrastructure issue that should be resolved before implementation. |
| `NO_MODULE` | Raw LogicBlock / pypulseq event / orchestration is already the appropriate abstraction. |

---

# 4. Official PyPulseq coarse scan

## 4.1 Cartesian GRE family

### Sources

- `write_gre.py`
- `write_gre_label.py`
- `write_gre_label_softdelay.py`

### Decomposition

The essential physics is:

```text
Excitation
  ↓
slice rephase + readout prephase + phase encode
  ↓
Cartesian readout + ADC
  ↓
phase rewind + gradient spoiler
  ↓
TR fill
  ↓
RF-spoiling schedule at acquisition level
```

### SeqCraft mapping

Already covered by:

- `Excitation`
- `PhaseEncode`
- `CartesianLine`
- `spoiler`
- `GRE2DTR`
- `GRE2D`

The labelled and soft-delay examples add scanner metadata/control but do not reveal a missing MRI physics primitive.

**Decision: `COVERED`.**

No `GRELabelled`, `GRESoftDelay`, or similar class should be introduced.

---

## 4.2 Multi-echo GRE

### Source

- Pulseq MATLAB `writeGRE_multiEcho_label.m`
- SeqCraft `examples/megre_2d`

### SeqCraft mapping

`CartesianLine` already owns:

- `echoes`;
- monopolar / bipolar polarity;
- flyback;
- echo spacing;
- true echo sample timing.

`GRE2DTR` and `GRE2D` simply inherit a longer readout.

**Decision: `COVERED`.**

This is strong confirmation that SeqCraft's existing decision **not** to create `MEGRE2D` is correct.

---

## 4.3 Gradient-echo EPI

### Sources

- `write_epi.py`
- `write_epi_label.py`
- Pulseq `writeEpi*.m`

### SeqCraft mapping

The difficult reusable portion is the EPI train itself, already represented by `EPI2D`.

The complete shot still needs an excitation and TE/TR/tail placement. SeqCraft already has a notebook-local `GREEPI2D`.

**Readout decision: `COVERED`.**

**Composite decision: `PROMOTE_NOTEBOOK` candidate**, but with a layer/naming change; see Section 9.

---

## 4.4 Spin-echo EPI

### Sources

- `write_epi_se.py`
- `write_epi_se_rs.py`
- Pulseq `writeEpiSpinEcho*.m`
- Pulseq high-performance ramp-sampled SE-EPI demos.

### SeqCraft mapping

Reusable pieces already exist:

- `Excitation`
- `Refocusing`
- `EPI2D(prephase=False)` or equivalent EPI spin-echo usage.

SeqCraft already has notebook-local `SEEPI2D`.

The high-performance demos mainly confirm that the readout abstraction must handle the ugly split-gradient / ramp-sampling details internally — which is already the reason `EPI2D` exists.

**Readout decision: `COVERED`.**

**Shot/kernel decision: `PROMOTE_NOTEBOOK`.**

No separate “ramp-sampled EPI” class should be created.

---

## 4.5 TSE / FSE / HASTE

### Sources

- PyPulseq `write_tse.py`
- PyPulseq `write_haste.py`
- Pulseq `writeTSE.m`
- Pulseq `writeTSEprop.m`
- Pulseq `writeHASTE.m`
- SeqCraft `examples/fse_2d`

### Repeated physics motif

```text
90 excitation
   ↓
[ crusher + 180 + crusher
  + phase encode
  + spin-echo readout
  + phase rewind ] × echo train
   ↓
tail/spoil
```

The reusable difficulty is not one gradient; it is the **CPMG/FSE echo-train timing and moment balance**.

SeqCraft's notebook-local `FSE2D` already solves:

- repeated refocusing;
- crusher balancing;
- spin-echo Cartesian readout;
- phase encoding per echo;
- echo spacing;
- echo-to-ky ordering;
- effective TE;
- dummy shots;
- partial-Fourier/HASTE-like operation.

### Candidate split

**Candidate A — `TSEShot` / `FSEShot`: `NEW_KERNEL` / extraction from notebook.**

This should own one excitation plus one echo train.

**Candidate B — `FSE2D`: `PROMOTE_NOTEBOOK` / `NEW_IMAGING`.**

This should own shot ordering/segmentation/dummies.

### HASTE decision

HASTE is a parameterization of the same physical family:

- one long echo train;
- partial Fourier;
- suitable ky ordering.

**HASTE: `COMPOSITION_ONLY` / configuration of FSE/TSE.**

Do not create a dedicated `HASTE` class unless the fine scan finds physics that cannot be represented by the FSE abstraction.

---

## 4.6 MPRAGE

### Sources

- PyPulseq `write_2Dt1_mprage.py`
- PyPulseq `write_3Dt1_mprage.py`
- PyPulseq `write_mprage.py`
- Pulseq `writeMPRAGE*.m`
- SeqCraft `examples/mprage_2d`
- SeqCraft `examples/mp2rage_2d`

### What is already covered

- inversion pulse + crusher: `IRPrep`;
- GRE repetition physics: `GRE2DTR`;
- 2D GRE train composition: existing examples.

### Real missing capability

The 3D implementations repeatedly need a second phase/partition axis.

`PhaseEncode` is already axis-generic, so the missing capability is **not** a `PartitionEncode` leaf.

The missing reusable level is the **3D GRE repetition / 3D acquisition composite**.

### Decisions

- `IRPrep`: `COVERED`.
- `PartitionEncode`: **do not create**; `PhaseEncode(axis='z')` is sufficient at leaf level.
- `GRE3DTR`: `NEW_KERNEL`.
- `GRE3D`: likely `NEW_IMAGING` after the kernel exists.
- `MPRAGE` / `MP2RAGE`: initially `COMPOSITION_ONLY` or later imaging-level promotion.

The main extraction target should be `GRE3DTR`, because both ordinary 3D GRE and inversion-prepared 3D GRE can reuse it.

---

## 4.7 Radial GRE

### Sources

- PyPulseq `write_radial_gre.py`
- Pulseq `writeRadialGradientEcho*.m`
- Pulseq `writeFastRadialGradientEcho*.m`
- OpenMRF radial readout.

### Repeated physics motif

A radial spoke repeatedly combines:

- a radial prephaser/start position;
- a readout gradient and ADC;
- an in-plane rotation angle;
- exact k-space centre/sample semantics.

The rotation angle changes per call; the spoke waveform itself is designed once.

### Candidate

**`RadialReadout` / `RadialSpoke`: `NEW_LEAF`, target `readout/`.**

Likely build arguments:

- angle / orientation;
- acquire;
- phase;
- possibly spoke index only as metadata, not waveform design.

The fine scan should determine whether full-spoke and center-out/UTE should be one class with a trajectory mode or two siblings.

---

## 4.8 UTE

### Sources

- PyPulseq `write_ute.py`
- Pulseq `writeUTE.m`
- Pulseq `writeUTE_rs.m`
- OpenMRF UTE readout.

UTE repeatedly introduces:

- asymmetric or center-out radial acquisition;
- very short RF-to-ADC timing;
- trajectory rotation;
- different prephasing from ordinary full-spoke radial GRE.

### Candidate decision

**`EXTEND_EXISTING` against the future radial candidate**, not immediately a standalone `UTE` class.

Fine scan question:

```text
Can RadialReadout naturally express both:
  full spoke through k=0
  center-out / half-spoke UTE
without a branchy API?
```

If yes, UTE is a kernel/imaging composition over `RadialReadout`.

If no, create a separate `CenterOutReadout` leaf.

Do **not** create a complete `UTE` class before resolving this leaf boundary.

---

# 5. Additional official Pulseq MATLAB candidates

## 5.1 Diffusion-weighted EPI

### Sources

- `writeEpiDiffusionRS.m`
- `writeEpiDiffusionRS_PMC.m`
- PulseqDiffusion single-refocused and twice-refocused DWI.

### Repeated physics

The reusable element is not “DWI-EPI”; EPI is only the readout.

The missing capability is diffusion weighting:

- one or more strong gradient lobes around refocusing pulses;
- direction vector;
- target b-value / b-tensor;
- exact timing;
- hardware limits;
- interaction with crusher/imaging gradients;
- single-refocused and twice-refocused variants.

### Candidate

**`DiffusionPrep` / `DiffusionEncoding`: `ARCHITECTURE_REVIEW` + high-priority new capability.**

Why architecture review is needed:

Current `encoding/` is described as gradients, no ADC, imposing a phase that will be sampled. That is a natural description of spatial phase encoding and also velocity encoding, but diffusion weighting is different: the ensemble dephasing is the contrast mechanism and the b-value depends on waveform timing and cross-terms.

Two viable designs should be compared in the fine scan:

1. **`preparation/DiffusionPrep`**
   - may contain refocusing and diffusion gradients;
   - exposes b-value/b-tensor and semantic timing;
   - directly reusable before EPI, spiral, radial, etc.

2. **`encoding/DiffusionEncoding`**
   - returns only gradient encoding;
   - caller places refocusing;
   - cleaner primitive, but exact b-value can become context-dependent.

The coarse scan does **not** freeze this API.

Strong recommendation: do not force diffusion into the current folder rule until a reference-preserving prototype exists.

---

## 5.2 Fat saturation / spectral saturation

### Sources

- high-performance EPI examples;
- spiral demos;
- diffusion repositories;
- OpenMRF `FAT` / `SAT`;
- spectroscopy water/fat suppression code.

A common motif is:

```text
frequency-selective RF saturation
        ↓
spoiler/crusher
```

`Excitation(thickness_mm=None)` can create a spectrally shaped RF waveform, but it intentionally marks RF use as `excitation`; that is the wrong semantics for saturation.

### Candidate

**`SaturationPrep`: `NEW_LEAF`, target `preparation/`.**

Recommended abstraction:

- frequency offset in Hz or ppm;
- pulse shape and duration;
- flip / B1;
- spoiler moment;
- semantic centre/tail timing.

`FatSat` should probably be an **example/preset policy**, not the fundamental module. A generic saturation prep can support fat saturation and become a building block for more advanced preparation trains.

---

## 5.3 Spiral

### Sources

- Pulseq `writeSpiral.m`
- `usc-mrel/rtspiral_pypulseq`
- OpenMRF spiral / spiral-TSE / rosette / cones infrastructure.

### Important observation

Spiral code naturally separates into two different things:

```text
trajectory design
k(t) / g(t), rewinder, constraints
        ↓
readout packaging
Gx/Gy + ADC + rotation + semantic k=0
        ↓
GRE / SE / FISP / fMRI / MRF kernel
```

The first object is **not** naturally a `Module`: it is numerical trajectory design and may not return a LogicBlock.

The research repositories independently separate spiral generation and rewinder design from sequence assembly.

### Candidates

**A. `SpiralReadout` / `SpiralInterleaf`: `NEW_LEAF`, target `readout/`.**

Owns:

- one designed interleaf;
- Gx/Gy;
- ADC timing/segmentation;
- k-space semantic information;
- per-call rotation/interleaf angle.

**B. Trajectory-design utilities: `ARCHITECTURE_REVIEW`, not a Module.**

Potential future location:

```text
seqcraft/
  trajectory/
    spiral.py
    radial.py
    rewinder.py
    utils.py
```

or under `seqcraft.design.trajectory`.

This layer should return numerical trajectories / gradient designs consumed by readout modules.

It should not become part of compiler/core.

---

## 5.4 bSSFP / TrueFISP

### Sources

- Pulseq `writeTrufi.m`
- OpenMRF BSSFP.

Distinctive reusable physics:

- fully balanced zeroth gradient moments per TR;
- RF/receiver phase cycling;
- readout and rewinds coupled across TR boundaries;
- fast minimum-TR timing.

This is not spoiled GRE with one flag changed.

### Candidate

**`BalancedSSFPTR`: `NEW_KERNEL`.**

Potential later complete scan:

**`BSSFP2D` / `BSSFP3D`: `NEW_IMAGING` only after kernel validation.**

Do not add a `balanced=True` option to `GRE2DTR`; the invariants are different enough to justify a separate kernel.

---

## 5.5 ZTE / PETRA

### Sources

- `writeZTE_Petra.m`
- sodium variant;
- OpenMRF center-out / 3D trajectory infrastructure.

The physics includes:

- gradient already present during RF/ADC;
- 3D radial direction set;
- gradient-to-gradient transitions;
- central k-space gap;
- PETRA single-point/shell filling.

This is too large to treat as one first-pass leaf module.

### Decision

- underlying center-out/3D radial readout: investigate after radial/UTE;
- direction-set/trajectory design: trajectory utility layer;
- PETRA: high-level imaging composition.

**Overall: `COMPOSITION_ONLY` until lower-level center-out primitives exist.**

---

## 5.6 PRESS / semi-LASER

### Sources

- `writePRESS.m`
- `writeSemiLaser.m`

They expose genuine reusable spectroscopy physics:

- orthogonal spatially selective RF pulses;
- crusher choreography;
- adiabatic refocusing;
- water suppression;
- long free-induction/spectroscopy ADC.

But the module boundaries differ substantially from imaging readout modules.

### Decision

For the current imaging-focused module expansion:

**`DEFER / ARCHITECTURE_REVIEW`.**

If spectroscopy becomes a first-class SeqCraft target, likely introduce a semantic application layer rather than calling `PRESS` a `readout/` merely because it contains an ADC.

A concrete nearer-term extension is:

**`Refocusing` adiabatic pulse support: `EXTEND_EXISTING` candidate.**

`IRPrep` already has adiabatic pulse plumbing; `Refocusing` currently targets sinc/SLR/Gaussian/block. The fine scan should test whether adiabatic refocusing fits cleanly in the same class or deserves a specialized sibling.

---

## 5.7 Cine GRE, triggers, gating and labels

### Sources

- `writeCineGradientEcho.m`
- EPI label examples;
- Open4DFlow triggering;
- PMC / trigger examples.

These repeatedly use:

- physiological input triggers;
- digital outputs;
- labels/counters;
- segmented ordering.

These are important for a complete scanner sequence, but they are generally **orchestration/control**, not reusable MR waveform physics.

### Decision

**`NO_MODULE` for now.**

Use raw pypulseq trigger/label events in LogicBlocks and let the imaging/kernel layer emit labels when it owns their semantics.

Do not create a `TriggerModule`, `LabelModule`, or `CardiacGRE` merely to wrap one event or one ordering policy.

Reconsider only if several future examples reveal repeated nontrivial timing arithmetic that genuinely shortens user code.

---

# 6. Research-repository cross-validation

## 6.1 PulseqDiffusion

The repository contains separate implementations for:

- SE-DWI-EPI;
- ramp-sampled SE-DWI-EPI;
- twice-refocused SE-DWI-EPI;
- b-value and b-vector utilities.

### What it validates

It strongly supports diffusion as an independent reusable capability rather than a property of EPI.

It also shows that the candidate must not be designed around only one Stejskal-Tanner pair. At minimum the fine scan should account for:

- single-refocused diffusion;
- twice-refocused diffusion;
- multiple directions;
- multiple b-values;
- TE optimization.

**Effect on candidate:** increases confidence in `DiffusionPrep/Encoding`.

---

## 6.2 rtspiral_pypulseq

The repository explicitly separates:

- spiral generation;
- gradient raster conversion;
- rewinder design;
- RF kernels;
- trigger/tagging/preparation kernels;
- final sequence assembly.

### What it validates

This is strong evidence for the proposed split:

```text
trajectory utility
    ↓
SpiralReadout
    ↓
contrast-specific kernel
```

A single `SpiralGRE` class would collapse too many independently reusable concerns.

**Effect on candidate:** high confidence for both `SpiralReadout` and a non-Module trajectory utility layer.

---

## 6.3 Open4DFlow

The sequence designs moment-targeted gradient waveforms across x/y/z, with explicit VENC-dependent first moments and cardiac/label orchestration.

### What it validates

Velocity encoding is a reusable capability independent of a particular k-space sampling family.

### Candidate

**`VelocityEncode`: `NEW_LEAF`, likely `encoding/`.**

Unlike diffusion, velocity encoding fits the current conceptual meaning of encoding more naturally:

> impose a controlled phase proportional to motion and sample that phase.

However, the exact moment design is coupled to surrounding imaging gradients. The fine scan must decide whether the module designs:

- only an incremental M1 encoding lobe, or
- a moment-compensated waveform given start/end moment constraints.

This should be solved numerically/physically, not hidden behind a thin wrapper.

---

## 6.4 OpenMRF

OpenMRF provides useful architecture evidence because it independently separates both preparation families and readout families.

Observed preparation families include:

- inversion;
- fat suppression;
- generic saturation;
- T2 preparation;
- spin lock;
- CEST;
- adiabatic spin lock.

Observed readout families include:

- GRE;
- EPI;
- BSSFP;
- radial;
- spiral;
- TSE;
- UTE;
- PRESS;
- spiral-TSE.

It also has trajectory design utilities for:

- radial;
- spiral;
- rosette;
- cones;
- Seiffert trajectories;
- rewinders.

### What it validates

1. SeqCraft's current split between preparation and readout is broadly useful.
2. The preparation category will become **composite**, not just a leaf-RF category.
3. Non-Cartesian readouts benefit from a trajectory utility layer.
4. T2 preparation deserves a real candidate.
5. BSSFP is independently reusable from spoiled GRE.

---

## 6.5 CEST repositories

CEST codebases emphasize reusable saturation preparation periods, often separately specified from the imaging readout.

### Candidate

**`CESTPrep` / `SaturationTrainPrep`: `NEW` but later priority.**

Do not force full CEST into the simple `SaturationPrep` candidate. A single fat-saturation pulse and a long pulse train with duty-cycle/B1/offset scheduling have different semantics and safety/timing needs.

A likely relationship is:

```text
SaturationPrep          # one saturation pulse + spoiler
SaturationTrainPrep     # repeated/off-resonant saturation schedule
CEST protocol           # composition/policy on top
```

The exact hierarchy should be decided during fine scan.

---

# 7. Consolidated ModuleCandidate map

## 7.1 Existing capabilities confirmed by the scan

| Capability | Current SeqCraft | Status | External evidence |
|---|---|---|---|
| selective excitation | `Excitation` | `COVERED` | GRE/EPI/TSE/MPRAGE/etc. |
| refocusing + crushers | `Refocusing` | `COVERED` | SE/TSE/HASTE |
| inversion prep | `IRPrep` | `COVERED` | MPRAGE/MP2RAGE |
| Cartesian phase encoding | `PhaseEncode` | `COVERED` | GRE/SE/TSE |
| Cartesian readout | `CartesianLine` | `COVERED` | GRE/SE |
| multi-echo Cartesian GRE readout | `CartesianLine` | `COVERED` | multi-echo GRE |
| EPI train | `EPI2D` | `COVERED` | GE-EPI/SE-EPI |
| gradient spoiler | `spoiler()` | `COVERED` | broad |
| spoiled GRE TR | `GRE2DTR` | `COVERED` | GRE/MPRAGE building block |
| complete 2D spoiled GRE | `GRE2D` | `COVERED` | GRE family |

## 7.2 Candidates

| Candidate | Proposed action | Proposed layer | Existing overlap | Evidence strength | Main question for fine scan |
|---|---|---|---|---|---|
| `SaturationPrep` | `NEW_LEAF` | `preparation/` | `Excitation`, `spoiler` | high | generic saturation vs FatSat-specific API |
| `TSEShot` / `FSEShot` | `NEW_KERNEL` / extract | `kernel/` | `Excitation`, `Refocusing`, `CartesianLine`, `PhaseEncode` | very high | exact reusable shot boundary/API |
| `FSE2D` | `PROMOTE_NOTEBOOK` | `imaging/` | notebook `FSE2D` | very high | extract shot first; preserve ordering freedom |
| `GREEPI2DTR` | `PROMOTE_NOTEBOOK` | `kernel/` | notebook `GREEPI2D`, `EPI2D` | high | naming and shot/TR contract |
| `SEEPI2DTR` | `PROMOTE_NOTEBOOK` | `kernel/` | notebook `SEEPI2D`, `EPI2D`, `Refocusing` | high | kernel vs full scan boundary |
| `GRE3DTR` | `NEW_KERNEL` | `kernel/` | `GRE2DTR`, `PhaseEncode` | high | reuse vs duplication with 2D |
| `GRE3D` | `NEW_IMAGING` | `imaging/` | future `GRE3DTR` | high | sampling API should be `(ky,kz)` data, not policy |
| `RadialReadout` | `NEW_LEAF` | `readout/` | `CartesianLine` arithmetic may inform it | very high | full-spoke + center-out in one abstraction? |
| UTE center-out mode | `EXTEND_EXISTING` candidate | `readout/` | future `RadialReadout` | high | extension vs separate `CenterOutReadout` |
| `SpiralReadout` | `NEW_LEAF` | `readout/` | none | very high | trajectory input/design contract; ADC segmentation |
| trajectory utilities | `ARCHITECTURE_REVIEW` | **not modules** | none | very high | public vs internal package location |
| diffusion weighting | `ARCHITECTURE_REVIEW` + new | `preparation/` or `encoding/` | `Refocusing` | very high | context-dependent b-value and refocusing ownership |
| `BalancedSSFPTR` | `NEW_KERNEL` | `kernel/` | `Excitation`, `CartesianLine`, `PhaseEncode` | high | balancing across TR boundaries |
| `VelocityEncode` | `NEW_LEAF` | `encoding/` | none | high | moment constraints and surrounding-gradient cross terms |
| `T2Prep` | `NEW` | `preparation/` | RF modules + spoilers | medium-high | composite preparation taxonomy |
| `SaturationTrainPrep` / `CESTPrep` | `NEW`, later | `preparation/` | future `SaturationPrep` | medium-high | generic train vs CEST policy |
| adiabatic refocusing | `EXTEND_EXISTING` | `rf/Refocusing` | `IRPrep` adiabatic support | medium | same class vs specialized sibling |
| ZTE/PETRA | `COMPOSITION_ONLY` initially | imaging | radial/trajectory primitives needed first | medium-high | wait for lower layers |
| MPRAGE / MP2RAGE | `COMPOSITION_ONLY` initially | imaging if promoted later | `IRPrep`, GRE kernels | high | avoid recipe proliferation |
| PRESS / semi-LASER | `DEFER / ARCHITECTURE_REVIEW` | possible spectroscopy domain | RF/refocusing/prep pieces | medium | whether MRS is first-class scope |
| cardiac trigger/gating | `NO_MODULE` | orchestration | raw Pulseq events | high | only revisit if nontrivial repeated arithmetic appears |
| labels/counters | `NO_MODULE` generically | kernel/imaging ownership | raw Pulseq labels | high | semantic labels stay with the layer that knows them |
| MRF schedules | `COMPOSITION_ONLY` | orchestration/imaging | multiple existing/future modules | medium | policy/schedule should not become a waveform leaf |

---

# 8. The strongest architecture finding: preparation cannot remain a "leaf only" category

Current documentation has an elegant simple taxonomy:

```text
rf/
preparation/
encoding/
readout/
kernel/
imaging/
```

but its current explanation mixes two different dimensions:

### Semantic role

- RF
- preparation
- encoding
- readout

### Composition depth

- kernel
- imaging

This is harmless with the current small library because most preparation modules are simple.

The coarse scan reveals the first real stress:

- `T2Prep` is semantically a preparation but contains several RF pulses and crushers.
- diffusion preparation may contain gradients plus one or more refocusing pulses.
- CEST preparation can contain a long train of saturation pulses and delays.
- WET / water suppression is a multi-pulse preparation.
- magnetization transfer / spin-lock preparations are also composite.

If the rule remains:

> "Anything composing modules from more than one leaf folder belongs in kernel/"

then `T2Prep` would incorrectly become a kernel, despite never containing an imaging readout or a repeating acquisition unit.

## Recommendation

Do **not** move existing modules now.

Instead, revise the conceptual rule when the first composite preparation is promoted:

```text
rf/           RF role primitives
encoding/     encoding primitives
readout/      ADC/readout primitives
preparation/  contrast/state preparation, simple OR composite
kernel/       repeating acquisition unit containing an imaging readout
imaging/      complete acquisition / scan composition
```

This keeps `T2Prep`, `DiffusionPrep`, CEST, etc. where users expect them.

The public imports remain flat, so this taxonomy change is cheap and does not need to affect user code.

---

# 9. Notebook-local promotion and layer-move recommendations

## 9.1 `FSE2D`

Current notebook-local `FSE2D` is a complete acquisition and therefore naturally belongs in `imaging/`.

However, the repeated physics inside it — one excitation plus N refocused echoes — is independently reusable.

### Recommendation

Extract:

```text
kernel/TSEShot.py
        ↓
imaging/FSE2D.py
```

or equivalent names.

HASTE then becomes an FSE configuration rather than a sibling class.

This is the strongest package-promotion candidate because:

- SeqCraft already has a working implementation;
- PyPulseq has independent TSE and HASTE examples;
- Pulseq MATLAB has independent TSE/HASTE examples;
- OpenMRF has a TSE family.

---

## 9.2 `GREEPI2D`

The notebook class name sounds like a complete 2D scan, but its reusable responsibility is a **single excitation + EPI shot with TE/TR/tail arithmetic**.

### Recommendation

If promoted, place the reusable unit in `kernel/` and rename around the unit, e.g.:

```text
GREEPI2DTR
GREEPIShot
```

Do not put the current shot object in `imaging/` merely because its k-space train is two-dimensional.

A later imaging class may stack segmented shots, references, dummies, slices and acquisition ordering.

---

## 9.3 `SEEPI2D`

Same conclusion as GREEPI:

- it is a spin-echo EPI shot/kernel;
- `EPI2D` stays the leaf readout;
- excitation/refocusing placement belongs in a reusable kernel.

Potential promotion:

```text
kernel/SEEPI2DTR
```

or a better shot-oriented name.

---

## 9.4 `MPRAGE2D` / `MP2RAGE2D`

These currently demonstrate that:

```text
IRPrep + GRE train composition
```

works without changing the compiler.

The coarse scan does not yet justify moving them into the package before 3D GRE primitives are addressed.

### Recommendation

Keep notebook-local for now.

Use the official MPRAGE corpus as evidence for:

1. `GRE3DTR`;
2. possibly `GRE3D`;
3. only then revisit whether MPRAGE deserves a packaged imaging abstraction.

---

# 10. Second architecture finding: non-Cartesian trajectories want a utility layer

Cartesian modules can design their gradients directly from a small number of closed-form parameters.

Spiral, rosette, cones, optimized radial and rewinders are different. Their core algorithm often produces:

```text
trajectory samples
    ↓
gradient waveform
    ↓
rewinder / moment correction
    ↓
ADC alignment
```

This numerical design is useful outside one Module and often has no natural `LogicBlock` output by itself.

## Recommendation

Add **no trajectory layer yet**, but reserve it explicitly as an architecture candidate:

```text
seqcraft/
  trajectory/
    spiral.py
    radial.py
    rewinder.py
```

Rules:

- pure numerical/design functions;
- no dependency on compiler;
- preferably no dependency on high-level modules;
- readout modules consume their output;
- validation can check gradient/slew and k-space independently.

This mirrors the successful separation observed in multiple spiral/reconstruction research codebases without copying their APIs.

---

# 11. Encoding taxonomy: diffusion and velocity should not be conflated

`PhaseEncode` encodes spatial position into phase/k-space.

`VelocityEncode` intentionally encodes velocity into phase and fits the broad notion of "encoding."

Diffusion weighting is different: its signal attenuation depends on the entire gradient history / b-tensor, not on a single phase value that remains for stationary spins.

## Recommendation

For now:

```text
encoding/
  phase_encoding.py
  velocity_encoding.py      # plausible future home

preparation/
  diffusion_prep.py         # provisional future home
```

Do not rename `encoding/` or add a new folder solely from the coarse scan.

The first reference-preserving diffusion prototype should be allowed to overturn this proposal.

---

# 12. Proposed near-term module tree if the first wave succeeds

This is **not** a commitment; it is the current coarse-scan target picture.

```text
modules/
  spoiler.py

  rf/
    excitation.py
    refocusing.py

  preparation/
    ir_prep.py
    saturation_prep.py
    diffusion_prep.py       # provisional
    t2_prep.py              # later

  encoding/
    phase_encoding.py
    velocity_encoding.py    # later

  readout/
    cartesian_line.py
    epi_2d.py
    radial_readout.py
    spiral_readout.py

  kernel/
    gre_2d_tr.py
    gre_3d_tr.py
    tse_shot.py
    gre_epi_2d_tr.py
    se_epi_2d_tr.py
    bssfp_tr.py

  imaging/
    gre_2d.py
    gre_3d.py
    fse_2d.py
```

Potential separate numerical layer:

```text
trajectory/
  spiral.py
  radial.py
  rewinder.py
```

Again: only candidates that survive fine-scan/reference validation should be added.

---

# 13. Priority recommendation

The priority should maximize **sequence families unlocked per new abstraction**, not maximize module count.

## Wave 0 — no new implementation

Freeze the coarse-scan result and fine-scan the highest-value boundaries.

## Wave 1 — strongest evidence and largest immediate coverage

### 1. TSE/FSE shot extraction

Why first:

- working SeqCraft implementation already exists;
- official PyPulseq TSE + HASTE;
- official Pulseq TSE + HASTE;
- research frameworks confirm the same family.

Unlocks:

- SE as turbo factor 1;
- FSE/TSE;
- HASTE;
- later variable-flip TSE.

### 2. `RadialReadout`

Unlocks:

- radial GRE;
- radial SE variants;
- UTE investigation;
- basis for some trajectory calibration and non-Cartesian examples.

### 3. `GRE3DTR`

Unlocks:

- ordinary 3D GRE;
- 3D MPRAGE;
- MP2RAGE-like trains;
- multiple mapping sequences.

### 4. `SaturationPrep`

Unlocks:

- fat saturation in EPI/GRE/spiral/DWI;
- common spectral prep;
- foundation for later saturation trains.

### 5. `DiffusionPrep` fine-scan prototype

Unlocks:

- SE-DWI-EPI;
- spiral DWI;
- other diffusion readouts.

This one should be prototyped early despite architecture uncertainty because the prototype is what will resolve the taxonomy.

## Wave 2 — high value but greater numerical/architectural work

- `SpiralReadout` + trajectory utility prototype;
- `BalancedSSFPTR`;
- `GREEPI2DTR` / `SEEPI2DTR` promotion;
- UTE/center-out radial refinement;
- `VelocityEncode`.

## Wave 3 — specialized preparation / domains

- `T2Prep`;
- CEST / saturation trains;
- spin-lock / T1rho;
- ASL preparation;
- ZTE/PETRA;
- spectroscopy;
- advanced MRF schedule families.

---

# 14. Revised provisional ModuleCandidate schema after the coarse scan

The scan suggests the original candidate schema needs a few explicit fields for overlap and architecture.

A **research YAML**, not runtime code:

```yaml
id: radial-readout
name: RadialReadout

capability:
  summary: "One in-plane radial spoke with prephasing/readout/ADC and per-call orientation."
  family: readout

status: NEW_LEAF

target:
  layer: readout
  proposed_path: modules/readout/radial_readout.py

existing_seqcraft:
  reused:
    - Excitation
    - spoiler
  overlaps_with:
    - CartesianLine
  replace_or_extend: null

evidence:
  - repo: pulseq/pypulseq
    path: examples/scripts/write_radial_gre.py
    role: authoritative-example
  - repo: pulseq/pulseq
    path: matlab/demoSeq/writeRadialGradientEcho.m
    role: authoritative-example
  - repo: OpenMRF/openmrf-core-matlab
    path: include_pulseq_toolbox/src_readouts/RAD/
    role: independent-reuse-evidence

physics:
  invariants:
    - k-space spoke orientation follows requested angle
    - k=0 sample has defined semantic position
    - ADC and readout gradient remain synchronized
  varying_parameters:
    - angle
    - fov
    - matrix
    - dwell_or_bandwidth
    - full_spoke_vs_center_out

architecture:
  needs_review: true
  note: "Fine scan must decide whether UTE center-out belongs in this class."

validation:
  reference_executable: unknown
  desired_checks:
    - waveform
    - k_adc
    - center_sample
    - rotation_equivariance
    - hardware_legality

priority:
  coverage_unlocked:
    - radial_GRE
    - radial_SE
    - UTE_candidate
  confidence: high
```

### New schema fields justified by this scan

The following fields now look important:

- `status`
- `target.layer`
- `existing_seqcraft.reused`
- `existing_seqcraft.overlaps_with`
- `existing_seqcraft.replace_or_extend`
- `architecture.needs_review`
- `architecture.note`
- `evidence[].role`
- `physics.invariants`
- `physics.varying_parameters`
- `priority.coverage_unlocked`

This is still intentionally loose YAML/Markdown.

No Pydantic/dataclass implementation is needed yet.

---

# 15. Suggested first fine-scan set

A fine scan should deliberately sample candidates with different failure modes.

| Reference | Why |
|---|---|
| SeqCraft `examples/fse_2d` + PyPulseq `write_tse.py` + Pulseq `writeTSE.m` | test package-promotion/mining from an already-working SeqCraft composite |
| PyPulseq `write_radial_gre.py` + Pulseq `writeRadialGradientEcho.m` | clean new readout leaf |
| PyPulseq/Pulseq UTE examples | test whether radial can be extended or needs a sibling |
| Pulseq `writeSpiral.m` + `rtspiral_pypulseq` | test trajectory utility vs Module boundary |
| Pulseq `writeEpiDiffusionRS.m` + PulseqDiffusion | test preparation/encoding taxonomy and context-dependent invariants |
| Pulseq `writeTrufi.m` + OpenMRF BSSFP | test genuinely distinct kernel |
| OpenMRF FAT + official fat-sat snippets | test small preparation extraction |
| PyPulseq `write_mprage.py` + Pulseq MPRAGE + SeqCraft MPRAGE notebook | test 3D GRE primitive vs high-level sequence wrapper |

This set is enough to stabilize `ModuleCandidate v1` before batch scanning much larger corpora.

---

# 16. What this scan says about SeqCraft's architecture

The strongest conclusion is **not** that SeqCraft needs many new layers.

The existing architecture remains sound:

```text
user / modules
    ↓
LogicBlock tree
    ↓
compiler
    ↓
Pulseq
```

The gaps are mostly in the reusable MRI library above LogicBlock.

Two small taxonomy refinements may become necessary as the library grows:

```text
1. preparation/ may contain composite preparations.
2. non-Cartesian trajectory design may live in a non-Module utility layer.
```

Neither requires changing `LogicBlock` or the compiler.

The current separation therefore survives the corpus scan well.

---

# 17. Main conclusion

A naive sequence-to-class conversion would produce many redundant classes:

```text
MEGRE
HASTE
MPRAGE
CineGRE
DWI-EPI
SpiralGRE
RadialGRE
UTE
...
```

The coarse scan instead points to a smaller reusable basis:

```text
Existing:
  Excitation
  Refocusing
  IRPrep
  PhaseEncode
  CartesianLine
  EPI2D
  GRE2DTR
  GRE2D

High-value additions:
  SaturationPrep
  TSE/FSE shot kernel
  GRE3DTR
  RadialReadout
  SpiralReadout
  DiffusionPrep / Encoding
  BalancedSSFPTR
  VelocityEncode

Infrastructure candidate:
  trajectory design utilities
```

From these, many sequence names become compositions rather than classes.

That is the direction most consistent with SeqCraft's current philosophy:

> **grow coverage by increasing compositional closure, not by mirroring a vendor sequence-name catalogue one class at a time.**
