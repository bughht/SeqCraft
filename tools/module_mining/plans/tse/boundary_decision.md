# TSE/FSE Fine Scan — Boundary Decision and Extraction

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-18, so that the pull request extracting `TSEShot` and `FSE2D` carries the evidence
> that justified it.  The two copies are identical today and there is nothing keeping them
> that way; see [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** TSE/FSE shot
**Date:** 2026-09-18
**Passes:** 7 — candidate boundary; then extraction
**Decided by:** human review (the four ownership questions raised in
[`sweep_report.md`](sweep_report.md) §3)

---

## 1. The decision

```text
NEW_KERNEL        kernel/TSEShot
PROMOTE_NOTEBOOK  imaging/FSE2D
EXTEND_EXISTING   readout/CartesianLine  -- intrinsic echo geometry only
```

| question | decision |
|---|---|
| Q1 `TSEShot` + `FSE2D`, or promote `FSE2D` whole? | **the split**, consistent with `GRE2DTR` -> `GRE2D` |
| Q2 where does the readout skew live? | the **geometry** in `CartesianLine`; the **compensation** in `TSEShot` |
| Q3 should `Refocusing` own the readout-axis lobe pair? | **no** |
| Q4 the three-axis coupled window solve? | kernel-owned; **no** generic abstraction yet |

### The rule behind them, which is the part that generalises

```text
leaf module     intrinsic physics / waveform geometry
kernel          coupling between multiple leaves
imaging module  acquisition / scan policy
compiler        Pulseq legality
```

and, operationally:

> If a component cannot determine the correct value using only its own physical contract and
> parameters, it should not own that event.

This is a sharper test than "avoid duplication", and it disagrees with it in at least one place
here — Q3, where the lobe pair straddles a refocusing pulse and looks like `Refocusing`'s work.
Its size depends on the *readout's* geometry, so a `Refocusing` that owned it would be carrying
downstream knowledge that is wrong for SE-EPI, a diffusion spin echo and spectroscopy.

### Project constraint recorded with it

> **Package-level SeqCraft Module implementations are Python-only for now.**

MATLAB Pulseq remains physics and architecture evidence and optional `.seq` validation; it is not
an implementation target. No MATLAB Module API, cross-language Module abstraction or MATLAB
implementation layer is to be designed. This is why the official PyPulseq `write_tse.py` was the
priority reference to close and `writeTSE.m` was not.

---

## 2. What was built

In the SeqCraft repository, branch `feat/module-mining-tse`:

| | |
|---|---|
| `src/seqcraft/modules/kernel/tse_shot.py` | `TSEShot` — one excitation and its echo train |
| `src/seqcraft/modules/imaging/fse_2d.py` | `FSE2D` — shots, and which lines each acquires |
| `src/seqcraft/modules/readout/cartesian_line.py` | `+ area_after_echo_per_m`, `+ echo_moment_imbalance_per_m` |
| `tests/modules/test_tse_shot.py` | 19 tests |
| `tests/modules/test_fse_2d.py` | 14 tests |
| `tests/modules/test_fse_notebook_matches_the_package.py` | 9 tests |
| `tests/modules/test_cartesian_line.py` | 3 added |
| `docs/api_reference.md`, `CHANGELOG.md` | the index, a worked example and the argument |

**`Refocusing` is untouched.** So is the compiler, `LogicBlock`, and every other module.

## 3. The proof that the extraction preserved behaviour

**Event for event, against the notebook it came from**, at four turbo factors:

| case | events | identical |
|---|---|---|
| `echoes=1`, 128 shots | 2176 | **yes** |
| `echoes=8`, +1 dummy | 1343 | **yes** |
| `echoes=16`, +1 dummy | 1335 | **yes** |
| `echoes=72`, `partial_fourier=0.625`, HASTE | 656 | **yes** |

Absolute times and content hashes of every leaf event, nesting and tags ignored — the package
nests a kernel where the notebook nests a method, and shape is not what has to match.

The fine-scan sweep was then re-run against the *extracted* modules rather than the notebook:
**8/8 cases pass every invariant, with numbers identical to the notebook's**. The extraction is
therefore checked against the evidence that justified it, not only against its own tests.

Repository checks: `ruff check .` clean; 902 passed / 19 skipped (was 857 before, so 45 new tests
and nothing broken); 72 doctests; strict mypy subset clean; `check_api_reference.py` clean with
the new example executing; all nine build notebooks execute.

## 4. Where the code pushed back, and what was done

Four places the guidance met the code and had to be resolved. None changes a decision; all four
are worth a reviewer's eye.

**4.1 The imbalance is the designed area minus the integrated one, not two integrations.**
`echo_moment_imbalance_per_m` returns `gx.area - 2*area_to_echo_per_m`. The more "measured"
spelling — integrating the lobe's own knots for the total as well — differs in the last few bits
(~6e-14 relative) and **moved every event**, because the lobes and the dephaser are designed from
it. `examples/se_2d/01` and `examples/fse_2d/01` both write the first spelling, and both are
pinned by CI. Preserving the pinned behaviour won; the docstring says which number it is.

**4.2 The `SEG` label: value from the acquisition, placement by the shot.** A naive reading of
"`FSE2D` owns segmentation" puts the label in `FSE2D`, which then has to compute where the first
readout starts — and moves the event. `TSEShot.build` takes `segment_index` and places it; the
shot cannot count its siblings, and the acquisition should not have to locate a readout.

**4.3 `TSEShot` owns `tr_s`, not `FSE2D`.** The guidance gave "shot duration" to `TSEShot` and
"complete scan assembly" to `FSE2D`, which leaves the repetition *period* ambiguous. `GRE2DTR`
owns `tr_s` and emits a block exactly `tr_s` long, with `GRE2D` stacking them; that precedent was
followed, since consistency with it was the stated reason for the split. `TSEShot` exposes
`shot_s` (what the events occupy), `min_tr_s` and `tr_s`. **Worth confirming.**

**4.4 The names carry the existing unit suffix.** `area_after_echo_per_m` and
`echo_moment_imbalance_per_m` rather than the conceptual `area_after_echo` / `echo_moment_imbalance`,
matching `area_to_echo_per_m` beside them.

## 5. Status

```text
GREEN
```

Against the playbook's list: an executable reference exists and runs; the candidate boundary is
decided; the comparator passes L0, L2, L3 and L4; the parameter sweep passes; the metamorphic
tests pass; the package test suite passes; **no compiler change was required**; provenance and
licences are recorded.

The MATLAB reference (E4) is still unrun. Under the Python-only constraint recorded in §1 it is
optional validation evidence rather than a gate, so it does not hold the candidate at YELLOW — but
it remains the one cheap piece of external evidence not yet collected.

What is left is ordinary review: a PR on `feat/module-mining-tse`, which has not been pushed.

## 6. What the pilot proved about the workflow

The stated goal was the workflow, not the class. Carried forward to the next candidate:

- **The comparator was wrong twice and measuring corrected it** — once about where the echo is
  (the reference has no DC sample), once about which axes an equal-area rule applies to (not the
  phase axis). A `True/False` comparator would have reported a false failure both times.
- **Reconciling the protocols produced findings before any comparison ran**: a B1 overrun the
  reference warns about and SeqCraft refuses, a readout time fixed in the script body, and a
  derived window assumed to be raster-legal.
- **The extraction was provable because the evidence came first.** The invariant table and the
  sweep existed before any code moved, so "did this change anything?" had an answer that was not
  the new code's own test suite.
- **The ownership rule did real work.** Asking *which layer can determine this value* settled Q2
  and Q3 in opposite directions — one piece of arithmetic moved down into a leaf, another was kept
  out of one — where "remove duplication" would have moved both.
