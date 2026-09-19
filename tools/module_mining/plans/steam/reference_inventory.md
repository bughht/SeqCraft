# StimulatedEcho / STEAM — Reference Inventory

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** working name undecided (`StimulatedEcho` / `STEAMKernel` / `STEAMPrep` all proposed)
**Date:** 2026-09-19
**Stage:** prospective supervised run, evidence phase. **No extraction attempted.**

The terminal classification was **not** fixed before this inventory was assembled. The coarse scan
listed six possible outcomes and this document is allowed to reach any of them.

---

## 1. What exists, and under what licence

### R1 — MRzero playground, three notebooks

| | |
|---|---|
| repo | `MRsources/MRzero-Core` @ `5732d86` |
| paths | `documentation/playground_mr0/mr0_STE_3pulses_5echoes_seq.ipynb`, `mr0_diffusion_prep_STEAM_2D_seq.ipynb`, `mr0_DREAM_STE_seq.ipynb` |
| licence | **AGPL-3.0 + a separate `EULA.txt`** |
| executable | yes; 12-notebook subset runs in the repository's own CI with NRMSE 0.01 regression |
| registry role | `discovery` + `oracle`, **`design-witness` explicitly excluded** |

**Independence: these are one witness, and not a design witness at all.** All three build with
PyPulseq (`pypulseq==1.4.2post1` pinned in the test workflow) and simulate with MRzeroCore. They
come from one group, one house style, one repository. Agreement among them is one design agreeing
with itself.

**Licence consequence:** physics and invariants may be extracted; implementation text may not. So
R1 can establish *that* a family exists and *what* its physics is. It cannot serve as acceptance
evidence for a SeqCraft implementation, because the comparison would be against a source we may
not copy from and whose construction shares the PyPulseq arithmetic we would be testing.

### R2 — official PyPulseq: **nothing**

Searched `imr-framework/pypulseq` @ `f2c582b`. `examples/scripts/` contains GRE, EPI, SE-EPI, TSE,
HASTE, MPRAGE (2D and 3D), radial GRE, UTE, PNS and soft-delay demos. A repository-wide search for
`steam` and `stimulated` returns **no matches**.

### R3 — official Pulseq MATLAB: **nothing**

Searched `pulseq/pulseq` `matlab/`. `demoSeq/` contains 40+ demos including PRESS, semi-LASER,
TSE, TSEprop, Trufi, ZTE/PETRA, spiral, diffusion EPI and selective-RF. A search for `steam`,
`stimulated` and `dream` across `matlab/` returns **no matches**.

`writePRESS.m` and `writeSemiLaser.m` use three or more RF pulses but they are **90–180–180**
double spin echoes with voxel-selective gradients on three axes. That is a different coherence
pathway from a stimulated echo and is not a substitute witness.

### R4 — OpenMRF: **unverified**

`OpenMRF/openmrf-core-matlab` is in the registry with a modular preparation library (FAT, SAT,
inversion, T2, spin-lock, CEST). A T2-preparation is a 90–180–90 store-and-restore structure and
is therefore the most plausible permissively licensed witness to part of this candidate's physics.

**No local checkout exists and this has not been inspected.** It is recorded as an open item, not
as evidence. Nothing in this document depends on it.

---

## 2. The evidence position, stated plainly

```text
executable, permissively licensed reference for a stimulated echo:   NONE FOUND
executable reference under AGPL + EULA, non-independent, one group:  THREE NOTEBOOKS
plausible permissive partial witness, uninspected:                   OpenMRF T2-prep
```

This is a materially different starting position from all three completed pilots:

| | primary evidence | independent corroboration |
|---|---|---|
| TSE/FSE | SeqCraft's own notebook + official PyPulseq | MATLAB original (same design) |
| Radial | official PyPulseq, executed | a second Pulseq formulation **and** OpenMRF |
| GRE3D | shipped `gre3d.seq`, measured | `fmrifrey/lps`, independent, for the selective path |
| **STEAM** | **one AGPL corpus, three notebooks, one group** | **none identified** |

GRE3D is the closest precedent — its slab-selective path had a single witness — but that witness
was permissively licensed and *independent of the Pulseq house style*, and the path was carried by
analytic arithmetic that depended on no implementation at all. Neither condition holds here yet.
