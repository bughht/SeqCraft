# StimulatedEcho / STEAM — Acceptance-Claim Options and Proposed Classification

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**For human review.** The skill proposes; it does not decide.

---

## 1. Acceptance-claim options, if the candidate were to proceed

Written as options because the boundary is not settled. Each is paired with the evidence it would
need, which is what separates them.

### Option A — event identity against a SeqCraft notebook

Not available. There is no SeqCraft notebook-local STEAM to promote. This is the TSE-shaped claim
and it requires an implementation already in the repository.

### Option B — trajectory / waveform agreement against an executable reference

The Radial-shaped claim. Requires an executable, permissively licensed reference at matched
protocol. **None exists** (`reference_inventory.md` §1). The only executable references are AGPL +
EULA and non-independent, so agreement with them would establish that we reproduced one group's
PyPulseq construction, not that the physics is right.

### Option C — analytic invariants, no executable reference

The GRE3D slab-selective-shaped claim. Would require stating the stimulated-echo contract as
signed moment and timing relations checkable without any external implementation. Plausible
targets:

```text
the dephasing moment before the store pulse equals the rephasing moment after the restore pulse
the stored state's phase is independent of position for zero applied moment
the echo appears at TE = tau after the third pulse, measured from RF effective centres
residual transverse magnetisation is spoiled to below a stated threshold
```

**This is the only option currently open**, and it is weaker here than it was for GRE3D. There,
the analytic path corroborated an *executable* reference on the other mode of the same kernel.
Here there would be no executable reference for anything.

### Option D — simulation agreement against MRzeroCore

Available in principle: MRzero is registered as an `oracle`, the playground notebooks run, and
SeqCraft already depends on MRzeroCore for Layer 3.

**It does not solve the problem.** The reference notebooks are *built with PyPulseq and simulated
with MRzeroCore*. A SeqCraft implementation simulated with the same oracle and compared to them
tests our construction against theirs through one simulator — one design and one instrument. That
is `VALIDATION_ORACLE_UNTRUSTED` and `REFERENCE_NOT_INDEPENDENT` at the same time.

---

## 2. Proposed terminal classification

```text
status:        UNDECIDED
traffic_light: YELLOW
reason_code:   REFERENCE_NOT_INDEPENDENT
```

Two findings, and the second is why the light is yellow rather than a stop:

**The broad candidate is not supported by the evidence.** The three uses do not share a module
contract — three pulse counts, three flip patterns, three gradient roles, two incompatible readout
couplings. What they share is a coherence pathway, which is spin physics rather than a boundary.
Considered on its own, "one STEAM module covering these three" would be `NO_NEW_MODULE`.

**But the narrow candidate cannot be assessed at all**, and that is the binding constraint. The
90/180/90 diffusion-prepared store of §1.2 may be a genuine contract. Deciding needs evidence that
does not exist in any permissively licensed source we hold: PyPulseq has no stimulated-echo demo,
Pulseq MATLAB has none, and the only executable references are AGPL + EULA, non-independent, and
from one group. `REFERENCE_NOT_INDEPENDENT` names that precisely.

### Why not the alternatives

| | why not |
|---|---|
| `NO_NEW_MODULE` | correct for the broad framing, but it would close the narrow one by silence. The narrow candidate was never assessed; recording it as "no module" would claim evidence we do not have. |
| `PHYSICAL_BOUNDARY_UNCLEAR` | the broad boundary is not unclear — it was analysed and came out negative. Using this code would overstate the uncertainty. |
| `RED - WRAPPER_ONLY` | would require an extraction to have been attempted and found empty. None was. |
| `GREEN` / `NEW_*` | no acceptance claim is available at all. |

### What would resolve it

```text
1. inspect OpenMRF's T2-preparation -- the most plausible permissive witness (uninspected, no checkout)
2. search for any permissively licensed stimulated-echo or STEAM implementation
3. if one is found, re-run this scan on the NARROW candidate only
4. if none is found, the narrow candidate is evidence-blocked, not rejected, and should be
   recorded as such rather than quietly dropped
```

---

## 3. What this run demonstrates about the skill

Recorded because the run was also a calibration.

**The outcome was not preloaded.** The coarse scan listed six possible results and named none as
expected; the prospective-run guardrail was committed before the evidence was assembled. The
result reached — a YELLOW on evidence grounds, plus a negative finding on the boundary the coarse
scan was most interested in — was not among the outcomes anyone had written down.

**It is the first stop the workflow has produced prospectively.** All three pilots were GREEN, and
no stop outcome had ever been exercised.

**Two gaps in the skill, found by using it:**

1. **The `status` vocabulary had no value for "stopped before deciding".** A YELLOW record still
   requires a `status`, and every available value asserts an action. `UNDECIDED` has been added,
   because forcing a provisional `NEW_LEAF` or `NO_NEW_MODULE` beside a YELLOW would put a
   decision in the record that nobody made. This is exactly the shape the historical TSE record
   needed when it read `MORE_EVIDENCE_REQUIRED / YELLOW`.

2. **Licence can block a candidate before physics does**, and nothing in the workflow said so.
   The pipeline orders evidence before contract, but the four prior candidates all had
   permissively licensed evidence, so the case never arose. It arises here, and it is generic.

**What it does not demonstrate:** that the workflow can refuse a candidate whose *physics* is
inadequate. This stop is about evidence availability. The `RED` outcomes remain untested, and a
run that reaches one on its own evidence is still outstanding.
