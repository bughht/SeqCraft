# MEGRE flow compensation — what a multi-echo train actually requires

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-27, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../README.md`](../../README.md) for which one to edit.

**Status:** evidence pass complete, and **implemented**. The evidence below is what the
implementation was designed from; §12 records what was decided after it.
**Date:** 2026-09-27
**Measured against:** `main` at `cacdca7` (PR #41 merged)

---

## 1. The question

`PhysicalDesignScope` solves for **one adjustable region** ending at **one semantic endpoint**,
with `M0` and `M1` as the conditions. A multi-echo gradient echo has four echoes. So:

> What does first-moment compensation actually require across a multi-echo gradient-echo train,
> and which part of the current physical-design model is insufficient, if any?

This was deliberately run as a measurement, not as a design. Nothing was extended, no solver was
written, and the architecture boundary was left intact so that whatever it could not express would
show up as a number rather than as an opinion.

## 2. What MEGRE is, structurally

**There is no MEGRE module.** A multi-echo gradient echo is `GRE2DTR(echoes=N, polarity=...)` —
the spoiled GRE whose readout was told to read the line more than once. `examples/megre_2d/`
composes it through `GRE2D`. So the question is about `CartesianLine`'s train and about the scope
that designs the gradients before it, not about a new kernel.

```text
monopolar   4 echoes, ESP 2500 us, polarity_of = [ 1,  1,  1,  1], fly-back between lobes
bipolar     4 echoes, ESP 2080 us, polarity_of = [ 1, -1,  1, -1], lobes exactly contiguous
```

`ro.te_s` is the sanctioned source of echo instants; `echo_spacing_s` is the nominal period and
the two differ under `bipolar`, which turns out to matter (§6).

## 3. Measured, uncompensated

Protocol: 220 mm FOV, 64×64, 5 mm, 500 Hz/px, 4 echoes, **line 16 of 64** so that `M0y` is the
line's own non-zero `ky` rather than an accidental zero. Moments integrated from the RF effective
centre `t0` to each echo, referenced to `t0`.

| | echo 0 | echo 1 | echo 2 | echo 3 |
|---|---|---|---|---|
| **monopolar** `M1x` | 1.183505e-1 | 1.234496e-1 | 1.285487e-1 | 1.336477e-1 |
| **bipolar** `M1x` | 1.194443e-1 | −3.189596e-2 | 1.288989e-1 | −2.244142e-2 |

`M0x = 0` at every echo in both. `M0y = −72.7273 1/m`, `M1y = −1.316364e-1`, `M0z = 0`,
`M1z = −4.259556e-1` — **constant across every echo**, both polarities.

Monopolar climbs by an equal step. Bipolar alternates and does not return to where it started.

## 4. Measured, with `flow_comp` on x + y + z

`flow_comp=` has one echo in its design — the first.

| | echo 0 | echo 1 | echo 2 | echo 3 |
|---|---|---|---|---|
| **monopolar** `M1x` | 3.33e-16 | 5.099068e-3 | 1.019814e-2 | 1.529720e-2 |
| **bipolar** `M1x` | 0.0 | −1.513403e-1 | 9.454545e-3 | −1.418857e-1 |

**`M1y = 1.11e-16` and `M1z = 2.22e-15` at every echo, both polarities**, with `M0y` still holding
the line at −72.7273.

So: **`y` and `z` are compensated at every echo; `x` is compensated at the first echo only.**

### Why y and z come free

The event listing over `[TE_0, TE_last]`:

```text
monopolar   x: 7 events (4 readout lobes + 3 fly-backs)    y: 0    z: 0
bipolar     x: 4 events (contiguous lobes, no fly-back)    y: 0    z: 0
```

The phase-encode winder and the slice rephaser both finish before the first echo and nothing on
those axes plays again until after the last. A first moment about a **fixed** origin stops accruing
where the gradient stops — `∫ G(t)(t−t0) dt` has a zero integrand wherever `G = 0`, however large
the lever arm has grown. This is the same property found in Stage C; here it is what makes the
multi-echo result free on silent axes.

**Contract B (every echo) already holds on `y` and `z` with no new machinery.**

## 5. The recurrence

All moments about the same origin `t0`. Splitting the integral at the intermediate echo and
expanding `(t − t0) = (t − TE_i) + (TE_i − t0)`:

```text
M0(TE_{i+1}) = M0(TE_i) + a_i                       a_i = ∫ G dt   over [TE_i, TE_{i+1}]
M1(TE_{i+1}) = M1(TE_i) + mu_i + (TE_i − t0) * a_i  mu_i = ∫ G(t)(t − TE_i) dt   over the same
```

On the readout axis **every echo is a `k_x = 0` crossing**, so `a_i = M0(TE_{i+1}) − M0(TE_i) = 0`
and the lever-arm term vanishes:

```text
M1(TE_{i+1}) = M1(TE_i) + mu_i            M1(TE_n) = M1(TE_0) + sum_{i<n} mu_i
```

`mu_i` is the first moment of the inter-echo waveform **about its own start**. It references
neither `t0` nor `TE_0`. Measured:

```text
                      interval    net a_i     dM1 about t0    the interval about TE_i
monopolar             0->1        2.84e-13    0.005099        0.005099
                      1->2        7.96e-13    0.005099        0.005099
                      2->3        3.98e-13    0.005099        0.005099
bipolar               0->1       -2.27e-13   -0.151340       -0.151340
                      1->2       -1.71e-13    0.160795        0.160795
                      2->3        5.68e-14   -0.151340       -0.151340
```

Agreement to every printed digit. The stronger control: the compensated train's first echo is
**1.04 ms later** than the uncompensated train's, so every `(TE_i − t0)` lever arm differs — and
the increments are identical to nine digits.

```text
monopolar uncompensated   TE_0 3.114 ms   mu = [0.005099068  0.005099068  0.005099068]
monopolar flow_comp       TE_0 4.154 ms   mu = [0.005099068  0.005099068  0.005099068]
bipolar   uncompensated   TE_0 3.126 ms   mu = [-0.151340269 0.160794814 -0.151340269]
bipolar   flow_comp       TE_0 4.166 ms   mu = [-0.151340269 0.160794814 -0.151340269]
```

**One initial condition, then a propagation fixed entirely by the readout train.**

### Closed forms

```text
monopolar   mu_i = mu, one constant of the train
            M1(TE_n) = M1(TE_0) + n * mu

bipolar     mu_i alternates between mu_A and mu_B
            M1(TE_2j)   = M1(TE_0) + j * (mu_A + mu_B)
            M1(TE_2j+1) = M1(TE_0) + j * (mu_A + mu_B) + mu_A
```

Predicted from the first one (monopolar) or two (bipolar) increments, against all four measured
echoes: max error **7.6e-15** and **1.3e-15**.

## 6. Why the bipolar pair does not close — separately from monopolar

`mu_A + mu_B = +9.4545e-3`, not zero. The textbook statement that a bipolar train is compensated
on alternate echoes assumes it is.

The cause is already documented in `CartesianLine`: the echo is at sample `pre_echo_samples` on a
forward lobe and `num_samples − 1 − pre_echo_samples` on a reversed one, so **bipolar echo times
alternate by exactly two dwells**. Verified, including the sign:

```text
matrix  BW Hz/px  dwell us  fwd/rev sample   gap us   M1x at echo 2
    64       250    62.500        32 / 31    125.00      1.8364e-02
    64       500    31.500        32 / 31     63.00      9.4545e-03
    64      1000    16.000        32 / 31     32.00      5.1364e-03
    65       500    31.000        31 / 32    -62.00     -9.3182e-03
```

The gap is two dwells in every row; when the two sample indices swap, the gap and the residual
**both** change sign. The magnitude is what a strip of readout two dwells wide is worth,
`delta * k_max` with `k_max ≈ (N/2)/FOV = 145.45 1/m` — 9.16e-3 predicted against 9.45e-3 measured
for the 63 µs row, the remainder being ramp geometry.

So the even echoes of a bipolar train carry a small dwell-set residual **in addition to** the odd
echoes carrying the full `mu_A`. This is a refinement of the textbook claim, not a contradiction of
it, and it is a property of a real sample grid.

## 7. What one window can and cannot buy

Everything the pre-echo gradients do enters the recurrence **only** through `M1(TE_0)`. Nulling at
a different echo changes which `M1(TE_i)` is zero; the increments are untouched. Tested directly:

```text
                   uncompensated minus its own first        measured with flow_comp        diff
monopolar   [0, 5.099068e-3, 1.0198135e-2, 1.5297203e-2]   [3.3e-16, ...same...]        8.88e-15
bipolar     [0, -1.51340269e-1, 9.454545e-3, -1.41885723e-1]  [0, ...same...]           1.50e-15
```

The compensated train **is** the uncompensated train, rigidly shifted. And therefore:

```text
peak-to-peak M1x across the train
  monopolar   uncompensated 0.015297203    with flow_comp 0.015297203
  bipolar     uncompensated 0.160794814    with flow_comp 0.160794814
```

**Identical to nine digits.** The spread is `sum mu_i`, a property of the readout, and *no* choice
of pre-echo gradients moves it. A single adjustable window chooses an **offset**, nothing more.

### What that costs, in velocity phase

`phi = 2 pi M1 v`, at v = 0.10 m/s (slow venous / CSF-like), degrees, worst echo of the train:

```text
                none    null TE1    centred
monopolar       4.81        0.55       0.28
bipolar         4.64        5.45       2.89
```

**For bipolar, nulling the first echo makes the worst echo worse than doing nothing** — a rigid
shift cannot cancel an alternation, and moving echo 1 to zero moves echo 2 further away. This is
the single most protocol-relevant number in the pass.

## 8. What "MEGRE flow compensated" should mean

The contract is **per axis**, and it is **not the same for the two polarities**. Evidence and
practice both point away from a single global answer:

```text
y, z            contract B (every echo) -- and it already holds, for free, on both polarities
x monopolar     A and B are physically DISTINCT -- after nulling TE1 the later echoes carry
                5.099e-3, 1.020e-2, 1.530e-2 s/m, which is not zero.  Under THIS protocol that
                residual is smaller than the bipolar train's; that is a statement about this
                protocol, not a general equivalence of the two contracts
x bipolar       a pre-echo offset alone achieves B at exactly one echo; contract A is worse than
                nothing at the worst echo; contract C is the best a single window can do and
                still leaves 2.9 deg
```

That is **contract D — per-axis mixed** as the honest description of what the sequence supports
today, with the per-axis answer determined by whether the axis is silent across the train.

Full first-moment nulling at every echo of a multi-echo train is a published capability rather
than a default — Wu et al., *A fully flow-compensated multiecho susceptibility-weighted imaging
sequence*, Magn Reson Med 76(2):478-489, 2016 — which is consistent with the finding that it is
not reachable from the pre-echo gradients alone.

**Recommendation for the public surface:** `flow_comp=` as shipped means "M1 = 0 at the first
echo", which on silent axes happens to be every echo. Contract A and contract B are physically
distinct on `x` under both polarities, and the size of the difference is protocol dependent, so
neither should be described as an approximation of the other.

## 9. Is the one-window / one-endpoint scope sufficient?

**For `y` and `z`: yes, completely.** The scope already delivers contract B on any axis that is
quiet across the train, and it does so without knowing the train exists.

**For `x`: it is sufficient for contract A and contract C, and it cannot reach contract B.** Not
because the solver is weak, but because `sum mu_i` lies **outside the adjustable region** — it is
produced by gradients that play after the endpoint. No re-solve of a pre-echo window can change a
number that a later gradient produces.

Reaching contract B on `x` means `mu_i = 0` for every interval: a condition on the **inter-echo**
gradients. The condition is on the interval as a whole, not on any one event in it:

```text
over [TE_i, TE_{i+1}]        dM0 = 0        dM1 = 0

    fixed contribution       post-echo part of readout_i
                             + pre-echo part of readout_{i+1}

    adjustable contribution  the fly-back / transition waveform

    so the adjustable waveform must cancel the COMBINED fixed interval contribution --
    it is NOT the statement "the fly-back's own first moment about TE_i is zero"
```

`dM1` here is origin-independent **because** `dM0 = 0`, which is what keeps each interval a local,
decoupled two-constraint problem rather than a coupled family solve. What differs by polarity is
whether an adjustable contribution exists at all:

```text
monopolar   there IS an adjustable region between echoes -- the fly-back.  Its net interval area
            is fixed by the need to return to the same k_x; its shape is not.

bipolar     in the CURRENT contiguous waveform family there is no adjustable inter-echo region:
            the lobes abut, which is what makes the echo spacing shorter.  This is a statement
            about that waveform family, NOT about bipolar readouts in general -- a compensated
            bipolar family may insert or redesign inter-lobe transition waveforms, at the cost
            of a longer ESP.
```

So the answer to "one initial condition + known propagation + corrections at checkpoints, or a
genuinely larger family-level solve?" is: **the former.** The structure is a chain of independent
local problems, not a coupled optimisation, and the decoupling is a theorem, not a heuristic — it
follows from `a_i = 0`. Monopolar has the adjustable region that structure needs; bipolar would
have to grow one.

## 10. Classification

**CASE D, with a CASE B interior.**

- **CASE D** is the outer answer, on two independent grounds. The two polarities differ in
  *kind*, not degree: the monopolar train has an inter-echo degree of freedom and the current
  bipolar waveform family does not, so reaching contract B means different work in each case —
  reshaping an existing region versus growing one. And the right contract differs by axis —
  silent axes already satisfy contract B; the readout axis does not, under either polarity.
- **CASE B** describes what an extension would look like *if* monopolar full compensation is
  wanted: a repeated **local** correction with a proven origin-independent recurrence. It is
  explicitly **not CASE C** — no second semantic endpoint, no multi-checkpoint scope, no global
  optimiser. Each fly-back solves its own two-constraint problem.
- It is **not CASE A**: the current machinery cannot express contract B on `x`, and an adapter
  would not change that.

The practical consequence: **nothing in `PhysicalDesignScope` needs to change.** An all-echo
compensation capability is a condition on the inter-echo intervals of `CartesianLine`'s train,
which are the readout's own gradients, not the scope's.

## 11. Non-goals honoured

No multi-checkpoint scope, no second adjustable window, no checkpoint IR, no global optimiser, no
MEGRE-specific solver class, no compiler timing logic, no solver implementation of any kind. The
only code added is `examples/megre_2d/03_flow_comp.ipynb` plus its registration.

## 12. Roadmap principle established by this pass

> **Roadmap items are MRI capabilities first.** Module, kernel, function, example, or
> `PhysicalDesignScope` usage is decided only after the physical contract and reuse boundary are
> evidenced.

This pass is the worked example: "MEGRE flow compensation" looked like a solver feature and turned
out to be a per-axis, per-polarity contract question whose answer removes the need for the solver
feature on three of the four cases and relocates it to the readout on the fourth.

## 13. What was decided from this, and built

The contract was fixed as **every acquired echo**, for every requested axis, with no option to ask
for less: if meeting it costs echo spacing, that cost is physical design rather than a reason to
weaken the guarantee. That is a decision *against* the "contract D, per-axis mixed" description in
§8, which described what the sequence supported at the time rather than what it should promise.

The implementation follows §9's structure exactly and confirms its two predictions:

- **`PhysicalDesignScope` did not change.** The multi-echo problem is self-contained inside the
  leaf that owns both the echoes and the gradients between them, so it is local design in
  `CartesianLine`.
- **Each interval is a local, decoupled two-constraint problem**, as the `a_i = 0` argument said
  it would be. No second semantic endpoint, no multi-checkpoint IR, no global optimiser.

Two things the evidence pass got wrong, or stated too strongly, and the implementation corrected:

- The bipolar conclusion. §9 as first written said the contiguous lobes left "nothing to reshape".
  Narrowed during review to a statement about *that waveform family*, and the implementation then
  showed the family can simply grow a transition — a balanced pair carrying pure first moment,
  since the bipolar interval's net area is already zero. It costs 59 % of the echo spacing, which
  **reverses** the uncompensated ranking of the two readouts.
- The shape of the monopolar condition. Writing it as "the fly-back's first moment about `TE_i` is
  zero" is wrong; the condition is on the interval, and the adjustable part must cancel the
  *combined* fixed contribution of one lobe's post-echo tail and the next one's pre-echo head. A
  regression asserts the transition's own first moment is **not** zero, so the weaker statement
  cannot be reintroduced as a simplification.
