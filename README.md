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

seq_out = sc.compile(seq, opts)          # a pypulseq.Sequence, or an exception
seq_out.write('gre.seq')                 # pypulseq's own writer; the definitions are already set
```

`sc.compile` returns a `pypulseq.Sequence` and nothing else. If the tree cannot become a legal
sequence it **raises**; if the compiler had to change a waveform to make it legal it **warns**.
There is no report object to inspect and no result wrapper to unpack.

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

`sim` (MRzeroCore, torch) and `recon` (sigpy) are what `examples/gre_2d/02` needs. Simulating and
reconstructing are **not** part of the package — seqcraft builds sequences, and both are downstream
jobs with heavy dependencies of their own.

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

**Failures you cannot forget to check.** Every legality problem raises, with the offending
number, the tag path it came from, the time it happens, and two concrete remedies with the values
already computed. There is no report object, because an object carrying findings is an object whose
findings can go unread — and the way that fails is a `.seq` the console refuses an hour later.

> **No provenance sidecar, for now.** An earlier version wrote a JSON file beside the `.seq`
> recording versions, git state and every `Opts` field. It went with the result wrapper and nothing
> has replaced it yet, so a written `.seq` does not currently say what produced it.

---

## When to reach for it

For one sequence, once, write raw pypulseq — it is a fine tool for that.

seqcraft pays for itself when you have a *family* of sequences, a loop over more than two axes,
gradients that must overlap without you hand-splitting blocks, or files that have to be reproducible
six months later. And when it does not fit, what you are holding is already the pypulseq object.

---

## The module library

Six names, and every one of them was **extracted from a working sequence** rather than designed:

```python
gre = sc.modules.GRE2D(opts=opts, fov_mm=220.0, matrix=(64, 64), thickness_mm=5.0)
seq = sc.compile(gre(lines=range(64)), opts, name='gre_2d')
seq.write('gre_2d.seq')
```

| | |
|---|---|
| `Excitation` | an RF pulse and, when selective, its selection gradient and rephaser |
| `PhaseEncode` | one Cartesian phase-encode blip, designed once and scaled per line |
| `CartesianLine` | prephaser, readout gradient and ADC as one design |
| `spoiler` | *n* turns of phase across a voxel — a function, because it earns nothing more |
| `GRE2DTR` | one repetition of a spoiled 2D gradient echo |
| `GRE2D` | the complete scan |

The previous library — 27 classes, 5 762 lines — was removed rather than migrated, because a recipe
is somebody else's sequence choices baked into library code and changing your own scan should never
mean editing a package. [ADR-003](docs/adr/003-scanner-and-module-reform.md) records that decision.
What is here now is what came back under a stricter rule: **write it in the notebook first, in raw
pypulseq, simulate it until the image is right, and only then extract it with the compiled output
held fixed.** A module that cannot be extracted without altering the sequence is not a module; one
whose extraction does not shorten the notebook is a wrapper.

The line between a module and a recipe is who keeps the sequence choices. `CartesianLine` computes
a prephaser that cancels the readout's ramp — arithmetic nobody should have to get right twice.
`GRE2D` takes the **list of phase-encode lines to acquire** rather than an acceleration factor, and
ships no generator for it, because which lines to acquire is a sequence-programming choice and the
right answer depends on the coil array, the object and the reconstruction together.

[`examples/gre_2d/`](examples/gre_2d/) is where all six came from, and a test asserts that what the
notebook writes and what the package ships compile identically.

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

`sc.compile(sc.LogicBlock('probe').add(0.0, component.pre()), opts)` is the whole block-level
contract check: a component that only works when something else happens to be beside it is not
reusable, and the compile checks the raster, the limits and block legality with a better message
than a separate assertion could give.

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

The one check the compiler structurally **cannot** make is that a call leaves the module alone.
It validates a tree; it never sees the second call, so `self.g.amplitude = -self.g.amplitude` in a
per-call method compiles cleanly every TR and produces a plausible but wrong image. Three lines,
and [`docs/writing_a_module.md`](docs/writing_a_module.md) explains why:

```python
before = {k: sc.events.content_hash(v) for k, v in vars(pe).items() if hasattr(v, 'type')}
pe(line=17); pe(line=17)
assert {k: sc.events.content_hash(v) for k, v in vars(pe).items() if hasattr(v, 'type')} == before
```

---

## Examples

| | What it covers |
|---|---|
| [`docs/api_reference.md`](docs/api_reference.md) | **Every public name in the package**, by layer, with a runnable example for each. Executed by CI, so it cannot drift. |
| [`01_getting_started.ipynb`](examples/01_getting_started.ipynb) | Blocks, `Opts` and `compile`; the overlap rules, provenance, the escape hatches, writing a file — and `sc.Module` at the end, once there is a reason for it. Uses no modules. |
| [`gre_2d/01_build.ipynb`](examples/gre_2d/01_build.ipynb) | A spoiled 2D GRE three times over — raw pypulseq, the leaf modules composed inline, and the same composition written as a module. Every module in `sc.modules` came out of it. Needs nothing but seqcraft. |
| [`gre_2d/02_simulate_and_reconstruct.ipynb`](examples/gre_2d/02_simulate_and_reconstruct.ipynb) | The same three `.seq` files against a BrainWeb phantom — PD, T1, T2, T2′, D and a synthesised B0, all six simulated — with an eight-element receive ring, then one reconstruction across three samplings. |

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
