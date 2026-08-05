# `examples/lib`

Standalone helpers used by the example notebooks. **Not part of the seqcraft package.**

seqcraft builds sequences. Simulating and reconstructing them are separate jobs with their own
heavy dependencies, so they live here and are imported by path:

```python
import sys; sys.path.insert(0, 'lib')      # from within examples/
from mr0_bridge import to_mr0, simulate, phantom_with_diffusion
from spiral_recon import reconstruct
```

| File | What it does | Needs |
|---|---|---|
| `mr0_bridge.py` | Converts a `CompiledSequence` into an MRzeroCore sequence, with diffusion and B0. MRzero's own `.seq` reader is trapezoid-only, so a spiral cannot be imported through a file. | `MRzeroCore`, `torch` |
| `spiral_recon.py` | Exact phase-accumulation reconstruction with off-resonance correction, solved by conjugate gradient. | `sigpy` |
