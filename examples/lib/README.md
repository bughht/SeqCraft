# `examples/lib`

Standalone helpers used by the example notebooks. **Not part of the seqcraft package.**

seqcraft builds sequences. Simulating and reconstructing them are separate jobs with their own
heavy dependencies, so they live here and are imported by path:

```python
import sys; sys.path.insert(0, '../lib')   # from within an example folder
from mr0_bridge import to_mr0, simulate, phantom_with_diffusion
from noncartesian_recon import Readout, reconstruct_shot
```

| File | What it does | Needs |
|---|---|---|
| `mr0_bridge.py` | Converts a `CompiledSequence` into an MRzeroCore sequence, with diffusion and B0. MRzero's own `.seq` reader is trapezoid-only, so a spiral cannot be imported through a file. | `MRzeroCore`, `torch` |
| `noncartesian_recon.py` | Exact phase-accumulation reconstruction with off-resonance correction, solved by conjugate gradient. Serves both readouts: a spiral's blur and an EPI's phase-encode shift are the same `exp(-2πi·df·t)` term. | `sigpy` |

## Two things that are per-readout, and are passed rather than guessed

Both were inferred from a spiral before the EPI example existed, and both were wrong for anything
whose echo is not its first sample.

| | spiral | ramp-sampled EPI |
|---|---|---|
| `Readout.t_echo_s` | the ADC delay — sampling begins at k=0 | ~⅓ of the way into the train, so `t_adc_s` is negative before it |
| `dcf=` in `reconstruct` | `'radial'`, the `\|k\|` Jacobian | `'speed'`, `\|dk/dt\|` — an EPI line is *sparsest* at k=0, so `'radial'` down-weights the signal |
| `k_reference=` in `to_mr0` | `'first_sample'` | `'echo'` — the readout prephaser is inside the same repetition, so rebasing on the first sample subtracts it away |

`Readout` raises if the sample times look absolute rather than readout-relative, which is the
failure that puts `2π·df·TE` of spatially varying phase into the operator.
