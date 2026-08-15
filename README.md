# seqcraft

Composable, verifiable MRI pulse sequence programming on top of [pypulseq](https://github.com/imr-framework/pypulseq).

Three things, and one of them is pypulseq's.

| | |
|---|---|
| **`sc.LogicBlock`** | A tree of pulseq events, each with a start time. Two attributes, one method. Anything may overlap anything. |
| **`sc.compile`** | Turns the tree into legal pulseq blocks: finds boundaries, sums gradients that share an axis, and validates the result against the amplifier. |
| **`pp.Opts`** | The scanner. Not wrapped, not subclassed — the same object you pass to `pp.make_trapezoid`. |

Everything else is a way of producing logic blocks, and seqcraft imposes no structure on the code
that produces one beyond the block itself. `sc.Module` is the standard shape for a *reusable*
component you write; a plain function works just as well. See [*Writing your own*](#writing-your-own).

```python
import math
import pypulseq as pp
import seqcraft as sc

opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
               rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)

rf, gz, gz_reph = pp.make_sinc_pulse(flip_angle=math.radians(15), duration=1e-3,
                                     slice_thickness=5e-3, apodization=0.5, time_bw_product=4,
                                     delay=opts.rf_dead_time, use='excitation',
                                     system=opts, return_gz=True)
gx = pp.make_trapezoid('x', flat_area=64 * 4.0, flat_time=3.2e-3, system=opts)
adc = pp.make_adc(num_samples=64, duration=3.2e-3, delay=gx.rise_time, system=opts)
gx_pre = pp.make_trapezoid('x', area=-gx.area / 2, duration=1e-3, system=opts)

time_to_echo = gx.rise_time + gx.flat_time / 2          # where k = 0 sits inside the readout
t_winders = pp.calc_duration(gz)

seq = sc.LogicBlock('gre')
for i, line in enumerate(range(-32, 32)):
    t0 = i * 20e-3
    pe = pp.make_trapezoid('y', area=line * 4.0, duration=1e-3, system=opts)
    seq.add(t0, rf, gz)
    seq.add(t0 + t_winders, gz_reph, pe, gx_pre)        # z, y, x — one block, no coordination
    seq.add(t0 + 8e-3 - time_to_echo, gx, adc)          # this *is* the definition of TE

out = sc.compile(seq, opts)
out.check().raise_if_failed()
out.write('gre.seq')
```

The three winders land in one pulseq block because they coincide, and the compiler says nothing —
they are on different axes, so there is nothing wrong. That is the point of the design: you say
what you mean, and never think about block boundaries.

**Set the dead times.** pypulseq defaults `rf_dead_time`, `rf_ringdown_time` and `adc_dead_time` to
**zero**, which is wrong on every real scanner. A sequence built on those compiles cleanly, passes
every check, and is refused or silently mangled at the console. They belong to your installation,
not to the scanner model, so no preset and no vendor database can supply them.

---

## Install

```console
pip install -e ".[dev,viz]"
```

Requires the pinned `pypulseq` fork and `numpy>=1.24`. Optional extras: `viz` (matplotlib),
`systems` ([PulseqSystems](https://github.com/nimpulseq/PulseqSystems) vendor limits), `dev`
(pytest, ruff), `docs`.

The simulation and reconstruction helpers in `examples/lib/` need `MRzeroCore`, `torch` and `sigpy`.
They are **not** part of the package — seqcraft builds sequences, and simulating or reconstructing
them are downstream jobs with heavy dependencies of their own.

---

## What it does for you

**Overlap without hand-splitting blocks.** Place a slice rephaser, a phase blip and a readout
prephaser at the same time on three axes; the compiler puts them in one block. Two gradients on the
*same* axis are summed with a warning naming both sources. Two RF or two ADC events overlapping is
an error that names both and the overlap in microseconds — including when only their dead times
touch, which pypulseq would otherwise reject 40 000 blocks later.

**Limits checked where the truth is.** Two individually legal gradients on one axis can sum to an
illegal one: an area-100 and an area-200 trapezoid on a 40 mT/m, 150 T/m/s system reach 93 % of the
amplitude limit and **189 % of the slew limit**. No component can see that in isolation, so
amplitude and slew are measured on the compiled waveform.

**Errors that name what to change.** pypulseq says `Amplitude violation (117%)`. seqcraft says which
of the three parameters fixes it, and what value would work.

**Waveform fidelity you can check rather than trust.** What the compiler emits is compared against
what the tree said, exactly — including a boundary landing inside a waveform, two gradients summed,
an arbitrary waveform that must not be resampled, and long trains where float error accumulates. The
one case where pulseq's two gradient representations genuinely cannot both be held is *reported*
with a bound, rather than being inexact quietly.

**Files that cannot lie about themselves.** `write()` takes no geometry, matrix or FOV arguments —
everything written comes from what was compiled. A JSON sidecar records the versions, the git commit
and dirty flag, the definitions, the file's sha256, and every field of the `Opts` it was built
against.

---

## When to reach for it

For one sequence, once, write raw pypulseq — it is a fine tool for that.

seqcraft pays for itself when you have a *family* of sequences, a loop over more than two axes,
gradients that must overlap without you hand-splitting blocks, or files that have to be reproducible
six months later. And when it does not fit, `CompiledSequence.seq` is the pypulseq object itself.

---

## The module library

**There isn't one, on purpose.** seqcraft ships no concrete modules at all: no `SincExcitation`, no
`EPIReadout`, no `MonopolarDiffusion`. What it ships is the tree, the compiler, and the contract.

The previous library — 27 classes, 5 762 lines — was removed rather than migrated. A recipe is
somebody else's sequence choices baked into library code, and changing your own scan should never
mean editing a package. The physics worth keeping from it (the b-value solve, the variable-density
spiral trajectory, the EPI ramp-sampling moment integral) was lifted out as plain functions into
[`salvage/`](salvage/) before the deletion.

[ADR-003](docs/adr/003-scanner-and-module-reform.md) records what was decided and why. What replaces
the library, and on what principle, is deliberately still open — each primitive should be written
only when a real sequence needs it, with the raw-pypulseq path kept beside it so the module has to
earn its place.

---

## Writing your own

A component takes part in a sequence by returning a `LogicBlock`. That is the entire contract:
there is no method name to match and no registry to join.

```python
def spoiler(opts, *, area_per_m=800.0):                  # a function is a component
    g = pp.make_trapezoid('z', area=area_per_m, system=opts)
    return sc.LogicBlock('spoil').add(0.0, g)


class VelocityEncode:                                    # so is a class of any shape
    def __init__(self, opts, *, m1_s_per_m, axis='y'):
        self.lobe = pp.make_trapezoid(axis, area=..., system=opts)

    def pre(self):
        return sc.LogicBlock('venc_pre').add(0.0, self.lobe)

    def post(self):
        return sc.LogicBlock('venc_post').add(0.0, sc.events.derive(self.lobe, ...))
```

Two outputs named for what they are, rather than one `build(part=...)` — seqcraft has no opinion
either way.

`sc.testing.assert_output(component.pre, opts)` gives either of them the block-level contract
checks — a well-formed block, deterministic output, gradients on the raster, per-axis limits, and a
clean compile on its own.

`sc.Module` is the standard shape for a component you intend to *reuse*: parameters in, one
`LogicBlock` out. Four members, and no more — `opts`, `tag`, `__call__`, and the abstract `build`
you write. It declares no duration (the block measures itself), checks no units, and holds no
scanner wrapper.

```python
class PhaseEncode(sc.Module):
    def __init__(self, *, opts, fov_mm, matrix, axis='y', tag=None):
        super().__init__(opts=opts, tag=tag)
        self.dk = 1e3 / fov_mm
        self.g = pp.make_trapezoid(axis, area=self.dk * matrix / 2, system=opts)

    def build(self, *, line=0) -> sc.LogicBlock:
        scale = line * self.dk / float(self.g.area)
        return sc.LogicBlock().add(0.0, sc.events.derive(self.g, ...))
```

`sc.testing.assert_all(pe, line=17)` then runs the whole suite against it. See
[`docs/writing_a_module.md`](docs/writing_a_module.md).

---

## Examples

| | What it covers |
|---|---|
| [`01_getting_started.ipynb`](examples/01_getting_started.ipynb) | Blocks, `Opts` and `compile`; the overlap rules, provenance, the escape hatches, writing a file — and `sc.Module` at the end, once there is a reason for it. Uses no modules. |
| [`_parked/`](examples/_parked/) | Two complete DTI acquisitions, spiral and EPI. **They do not run against this version** — they were built on the deleted library, and are kept as the specification for what replaces it. |

---

## Documentation

- [`docs/architecture.md`](docs/architecture.md) — the layering, and what is deliberately absent.
- [`docs/compiler.md`](docs/compiler.md) — how block boundaries are chosen, how an event's own `delay` is handled, and what every warning means.
- [`docs/writing_a_module.md`](docs/writing_a_module.md) — the `Module` contract, components that inherit nothing, and what to assert.
- [`docs/testing.md`](docs/testing.md) — the test tiers and what each is for.

---

## Tests

```console
pytest tests --doctest-modules src/seqcraft
```

The compiler directory is the heart of it: one case per rule, plus the adversarial ones — a gradient
straddling an RF, a boundary that would fall inside an ADC window, a split mid-ramp, two RFs whose
dead times overlap. Every fixture is raw pypulseq, so compiler coverage never depends on whatever a
module library happens to contain. Every physics assertion is a number an independent calculation
gives, not a number the code happened to produce.

---

## A note on vendor data

Siemens `.asc` gradient descriptors carry proprietary PNS and acoustic-resonance coefficients. They
are **never** stored in, copied into or read from this repository:
`sc.hardware.load_hardware()` resolves them only through `$SEQCRAFT_ASC_DIR` and rejects anything
that looks like a path. `sc.hardware.synthetic_hardware()` provides a vendor-free stand-in so PNS
checks can run anywhere — it is not a real scanner and must never be used to clear a sequence for
human scanning.

## Licence

MIT.
