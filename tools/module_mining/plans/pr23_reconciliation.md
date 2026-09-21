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

> **Closure pass, 2026-09-20.** Phase A2 finished, so every disposition below that read
> "REBUILD NOW" or "REVIEW + SIMPLIFY NOW" has been replaced by what actually happened. Only
> dispositions that A2 changed were touched; no archaeology was redone. **Every meaningful
> capability in PR #23 now has a final disposition**, and §5 lists what stays deferred so that a
> future reader does not mistake silence for rejection.

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
| `GRESpiral2D` / `SESpiral2D` public Modules | **REJECTED** | composition; `gre_spiral_2d/` and `se_spiral_2d/` both define no class, and `se_spiral_2d/01` shows what the composition actually owns — the placement arithmetic that puts an origin crossing on a spin echo |
| reconstruction inside `src/seqcraft` | **REJECTED architecturally** | would put an ESPIRiT dependency on the compile path — PR23's own argument, and it holds |

## 2. Reconstruction side — Phase A's actual work

`examples/noncartesian_recon.py`, 323 lines. **Reviewed and rebuilt, not restored** — and
smaller than PR #23's version because the parts with one consumer stayed out.

**Final disposition: SIMPLIFIED / PORTED AS EXAMPLE INFRASTRUCTURE.** It lives under `examples/`
beside `phantom.py`, has two real consumers and four notebooks, and was deliberately **not**
promoted into `src/seqcraft`. Shipping two consumers did not trigger promotion; that is the
decision, recorded in `spiral/candidate.yaml`'s `evidence_state` with the trigger that would
change it.

### Its four design arguments, each from a real failure

| | disposition |
|---|---|
| **the operator is composed from sigpy primitives**, so `A.H` is adjoint by construction — hand-writing both is where a forward and an adjoint come to disagree, and CG on a nearly-adjoint pair converges to the wrong image without raising | **PORTED** — the strongest argument in the file |
| **interleaves are solved together**; separate solves plus averaging throws away the joint conditioning that makes multi-shot work | **PORTED** |
| **DCF is a preconditioner, not a data weighting** — `sqrt(w)` on both sides is a *different estimator* | **REJECTED as stated; the behaviour is PORTED.** The distinction does not exist: both spellings assemble `A^H W A x = A^H W y`, since `(W^(1/2)A)^H (W^(1/2)A) = A^H W A`. Checked on the operator to 1e-16 relative. Neither is a preconditioner in the numerical sense — a preconditioner leaves the fixed point, a weight moves it. **Found by CI**, after a first test asserted the two converge differently and failed on floating-point ordering; the assertion was false, not flaky. The weighting *behaviour* is standard and stays the default; what was wrong was the name and the claim. Settled numerically: on noiseless consistent data weighted and unweighted agree to 0.0003 because both drive the residual to zero, and they separate once there is noise |
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
| coil sensitivities / ESPIRiT | **DEFERRED** | all four notebooks are single-coil and none needed it. Trigger unchanged: a notebook that does |
| off-resonance time segmentation | **PORTED**, and it is the point of `gre_spiral_2d/02` | 108 pixels back to 1 at eight segments; and four segments measured *worse* than none, which is in the notebook because it is the wrong lesson a reader would otherwise draw |
| exact dense off-resonance operator | **PORTED as a small numerical check**, in `gre_spiral_2d/02` §3 | kept small, as the brief asked: one table showing the segmented operator converging to the exact sum — 100 % error at one segment, 7e-3 at sixteen, with a floor there set by the interpolation between segment centres |
| row-sum preconditioner | **DEFERRED** | not needed by a first B0-free reconstruction |
| spectral-norm estimate | **DEFERRED** | internal |
| field-map helper (two-echo) | **REJECTED as unnecessary** | `gre_spiral_2d/02` corrects against a *known* imposed offset, which is the honest thing for a teaching notebook and is stated as a limit rather than hidden. A two-echo estimator would add a second unvalidated step to a demonstration of the first. Trigger: a claim that survives an *estimated* field map |
| gradient-delay experiment (`delayed()`) | **DEFERRED** | trigger: a claim about measured trajectory fidelity. Hardware-sensitive and not needed by any current teaching claim |
| the eight-experiment spiral study | **REJECTED as a structure** | replaced by three notebooks with one question each: how non-Cartesian reconstruction works, why a long spiral readout cares about off-resonance, and what refocusing changes. The individual experiments it contained are dispositioned in §5 |

### The notebooks PR #23 proposed

| | disposition |
|---|---|
| GRE spiral build notebook | **REBUILT** as `examples/gre_spiral_2d/01_build.ipynb` (merged in v0), extended in A2 with a single-shot file because the off-resonance claim needs a long readout |
| GRE spiral reconstruction notebook | **REBUILT** as `examples/gre_spiral_2d/02_simulate_and_reconstruct.ipynb`, against the current architecture and the shared contract |
| SE spiral build notebook | **REBUILT** as `examples/se_spiral_2d/01_build.ipynb`, from current modules — `Excitation`, `Refocusing`, `SpiralReadout` — and not from PR #23's structure. Its subject is semantic timing ownership, which PR #23's version did not isolate |
| SE spiral reconstruction notebook | **REBUILT** as `examples/se_spiral_2d/02_simulate_and_reconstruct.ipynb` |

Building the SE pair found a defect in `SpiralReadout` that the GRE pair could not:
`time_to_echo` was documented as block-relative and implemented arm-relative, which is invisible
for `variant='out'` and 300 µs wrong for `'in-out'`. Fixed with a regression before the notebook
was written; recorded in `spiral/findings.md` §10.

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

**Done, and it paid immediately.** `tests/examples/test_noncartesian_recon.py` implements every
line of that list against an explicit dense DFT that calls no SigPy, plus a source check that the
module contains no trajectory branching. The dense DFT caught a real convention error on first
run: SigPy's NUFFT is built on its *centred, unitary* FFT and divides by `sqrt(N_pixels)`, so the
adapter was wrong by exactly 15/16 at matrix 16 — a scale error invisible in any windowed
magnitude image and fatal the moment two reconstructions are compared.

Four notebooks later, the core still has no trajectory branch.

## 4. Anything else in PR #23

Also present, and dispositioned rather than ignored: the `.seq`-writing GRE and SE spiral build
notebooks (**REBUILD**, against the current architecture, not restored); `tests/modules/
test_spiral_2d.py` and `test_spiral_notebooks.py` (**RE-DERIVED** — the current suite tests family
relationships rather than restating each variant's numbers); the CHANGELOG rationale and README
entries (**SUPERSEDED**); and the `seqcraft[recon]` extra (**PORTED** — `sigpy` stays an optional
example dependency).

## 5. What stays deferred after Phase A

Listed together so that silence is not mistaken for rejection. Each already carries a trigger
above or in the relevant `evidence_state`; this is the index, not a second record.

```text
Spiral echoes > 1                 an extension of the same leaf.  Contract written in
                                  echoes_contract.md; the plural machinery is already exercised
                                  by out-in at echoes = 1

gradient-delay study              trigger: a claim about measured trajectory fidelity.  It is
                                  hardware-sensitive, it blurs a spiral in a way that looks like
                                  off-resonance, and no teaching claim here needs it

multi-echo / T2* study            trigger: a protocol with more than one echo per excitation.
                                  se_spiral_2d/02 explicitly does NOT quantify T2' -- it measures
                                  that it cannot, at one echo time

variable-density comparison       density is implemented and asserted; what is deferred is a
                                  second witness that the CONVENTION matches another
                                  implementation's

advanced reconstruction           coil sensitivities / ESPIRiT, row-sum preconditioner,
comparisons                       spectral-norm estimate, the analytic spiral Jacobian.  Each has
                                  no consumer among the four notebooks

two-echo field-map estimator      trigger: a claim that survives an estimated rather than a known
                                  field map

labels (LIN/ECO/SEG/SET/...)      trigger: a consumer that reads them

oblique_trapezoid,                one consumer each, still.  Four notebooks did not change that:
sample_quantum,                   they are examples, and an example is not a second consumer of a
segment_samples                   private helper

reconstruction promotion          examples/noncartesian_recon.py has two real consumers and stays
into src/seqcraft                 under examples/.  Not an omission -- see section 2
```

**One item was added rather than closed by Phase A**, and it is a **design-policy consistency
revisit rather than an open defect.** `SpiralReadout`'s arm uses a conservative vector-norm slew
ceiling on purpose, because that is what makes it rotationally invariant; its prephaser and
rewinder only need the package's per-axis hardware limits, and they meet them — per-axis slew
stays at or below 99.9 % of `max_slew` at every interleaf angle measured in `se_spiral_2d/01`.
The 138.9 % figure at 45° is the *vector norm* across two simultaneous ramps, which the package
contract treats as informational by design and not as a limit. So the open question is whether
the auxiliary lobes should adopt the arm's conservative policy for consistency or for extra
headroom — not whether anything is wrong. Recorded in `spiral/candidate.yaml`.
