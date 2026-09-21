# PR #23 reconciliation

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-20, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../README.md) for which one to edit.

PR #23 was closed **superseded, not rejected**. This gives every meaningful capability in it
exactly one disposition, so that nothing useful is silently lost and nothing obsolete is silently
restored.

```text
PORTED                     in current main, independently re-derived or carried across
SIMPLIFIED / RE-DERIVED    the idea survives; the implementation does not
DEFERRED                   useful, not now, with a trigger
REJECTED / NOT ESTABLISHED do not port
```

It remains a **research corpus**, not a merge candidate. Nothing is cherry-picked wholesale.

---

## 1. Sequence side — settled by the Spiral calibration

| capability | disposition | note |
|---|---|---|
| four-mode Spiral family (`out`/`in`/`in-out`/`out-in`) | **PORTED** — independently re-derived | retained from the coupling test, not because PR23 had it |
| `g = 0` endpoint policy | **PORTED** — independently re-derived | and quantified before adoption: 0.20–0.78 % of readout duration |
| path / traversal separation | **PORTED** as an *internal* boundary | `writeSpiral.m` arrives at it independently; it is not two public modules |
| emitted-trajectory truth | **PORTED** | `k_per_m` integrated from emitted knots; agrees with `sc.kspace` to ~1e-12 1/m |
| ADC segmentation against the sample limit | **PORTED** | plus a defect PR23 did not have: the acquisition must finish *inside* its gradient |
| join correction on the assembled waveform | **PORTED** — re-derived as a rewinder sized from the emitted end | PR23's failure E, reproduced and fixed |
| `density` as sampled FOV multipliers | **PORTED**, with a recorded caveat | one witness; `writeSpiral.m` has no density parameter at all |
| labels (`LIN`/`ECO`/`SEG`/`SET`/`REF`/`IMA`) | **DEFERRED** | the current leaf emits none; trigger is a consumer that reads them |
| `Spiral2D` as a public class name | **REJECTED** | re-derived as `SpiralReadout`; the name was not the finding |
| Spiral `echoes > 1` | **DEFERRED** | an extension of the same leaf; contract written in `echoes_contract.md` |
| **m1 compiler defect** | **PORTED independently via PR #30** | reproduced with no Spiral module; fixed as a floor on peak \|g\| rather than PR23's larger relative coefficient |
| **m0 compiler claim** | **NOT ESTABLISHED — do not port** | the comment/code mismatch is real; the failure never reproduced. Revisit only with a candidate-free reproducer |
| `oblique_trapezoid` | **DEFERRED** | one consumer; kept private inside the module |
| `sample_quantum` | **DEFERRED** | one consumer |
| `segment_samples` | **DEFERRED** | one consumer; `check_event_sizes` already exists on main |
| `GRESpiral2D` / `SESpiral2D` public Modules | **REJECTED** | composition; the example defines no class |
| reconstruction inside `src/seqcraft` | **REJECTED architecturally** | would put an ESPIRiT dependency on the compile path — PR23's own argument, and it holds |

## 2. Reconstruction side — Phase A's actual work

`examples/noncartesian_recon.py`, 470 lines. Reviewed, not restored.

### Its four design arguments, each from a real failure

| | disposition |
|---|---|
| **the operator is composed from sigpy primitives**, so `A.H` is adjoint by construction — hand-writing both is where a forward and an adjoint come to disagree, and CG on a nearly-adjoint pair converges to the wrong image without raising | **PORTED** — the strongest argument in the file |
| **interleaves are solved together**; separate solves plus averaging throws away the joint conditioning that makes multi-shot work | **PORTED** |
| **DCF is a preconditioner, not a data weighting** — `sqrt(w)` on both sides is a *different estimator*, and on a variable-density spiral the weights span two orders of magnitude | **PORTED as the default**; the `weighting='data'` comparison path is **DEFERRED** |
| **the DCF knows about the density** — `\|k\|` and `\|dk/dt\|` are both wrong on a tapered spiral | **SIMPLIFIED**: Pipe–Menon generic by default, `\|k\|` radial weighting as a family option, the analytic spiral Jacobian kept **outside** the shared core |

### Feature-by-feature

| feature | disposition | why |
|---|---|---|
| trajectory/timing container | **PORTED, RENAMED** | it is trajectory-agnostic — Radial fits it unchanged. The historical name `SpiralReadout` must not survive on the reconstruction side; it now collides with the Module |
| explicit echo time, never inferred from `min(t_adc_s)` | **PORTED** | right for spiral-out for exactly one reason and wrong for `in`, `in-out`, multi-echo and EPI |
| refusal of sample times that look absolute | **PORTED** | the off-resonance term needs echo-relative time; absolute times give a structured, plausible artefact |
| `from_sidecar` `.npz` convention | **SIMPLIFIED** | keep the idea, let the notebooks own the file layout |
| SigPy NUFFT + coordinate adapter | **PORTED**, and **independently validated** | this is the part we own; see §3 |
| Pipe–Menon DCF | **PORTED** as the default | needs no analytic density and is what a real pipeline uses |
| `\|k\|` radial DCF | **PORTED** as a trajectory-family option | outside the generic interface |
| `\|dk/dt\| / FOV(\|k\|)` spiral Jacobian | **DEFERRED / notebook-local** | genuinely spiral-specific; does not belong in a shared core |
| CG solve | **PORTED**, simplified |  |
| coil sensitivities / ESPIRiT | **DEFERRED** | §5 of the brief starts single-coil; add only when a notebook needs it |
| off-resonance time segmentation | **PORTED**, and it is the point of the Spiral `02` |  |
| exact dense off-resonance operator | **DEFERRED** as a small numerical check | superseded as *validation* by the direct-DFT reference, which calls no SigPy |
| row-sum preconditioner | **DEFERRED** | not needed by a first B0-free reconstruction |
| spectral-norm estimate | **DEFERRED** | internal |
| field-map helper (two-echo) | **PORTED** | the Spiral `02` needs a field map to correct against |
| gradient-delay experiment (`delayed()`) | **DEFERRED** | trigger: a claim about measured trajectory fidelity. Hardware-sensitive and not needed by any current teaching claim |
| the eight-experiment spiral study | **DEFERRED** | replaced by one narrative a reader can follow |

## 3. What must be validated independently

SigPy is a backend, not an oracle. The adapter is ours:

```text
SeqCraft emitted trajectory -> physical k [1/m] + sample timing -> our adapter
    -> SigPy coordinate convention -> NUFFT
```

Checks that must not call SigPy:

```text
coordinate convention      scale, axis order, sign, orientation, pixel-centre
direct DFT reference       exp(-2 pi i k . r), explicit matrix, small matrix only
forward-model point test   one off-centre point; the complex PHASE must match analytically
adjoint test               <Ax, y> ~ <x, A^H y> on the composed operator
sample ordering            reconstruction consumes samples in acquisition order
timing semantics           block-relative vs echo-relative vs absolute, kept distinct
```

Both `RadialReadout` and `SpiralReadout` drive the same shared core. If it needs
`if radial: ... elif spiral: ...` for ordinary coordinate conversion or operator construction,
it is not shared — it is two implementations behind one name.

## 4. Anything else in PR #23

Also present, and dispositioned rather than ignored: the `.seq`-writing GRE and SE spiral build
notebooks (**REBUILD**, against the current architecture, not restored); `tests/modules/
test_spiral_2d.py` and `test_spiral_notebooks.py` (**RE-DERIVED** — the current suite tests family
relationships rather than restating each variant's numbers); the CHANGELOG rationale and README
entries (**SUPERSEDED**); and the `seqcraft[recon]` extra (**PORTED** — `sigpy` stays an optional
example dependency).
