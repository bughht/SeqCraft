# seqcraft

Composable, verifiable MRI pulse sequence programming on top of [pypulseq](https://github.com/imr-framework/pypulseq).

Three concepts, and no more.

| | |
|---|---|
| **`LogicBlock`** | A tree of pulseq events, each with a start time. Two attributes, one method. Anything may overlap anything. |
| **`Module`** | A reusable sequence task. `__init__` designs, `build()` returns a logic block, and timing a caller needs is a plain property. |
| **`sc.compile`** | Turns the tree into legal pulseq blocks: finds boundaries, sums gradients that share an axis, and validates the result against the amplifier. |

```python
import pypulseq as pp
import seqcraft as sc

system = sc.System.preset('prisma')
exc = sc.modules.SincExcitation(system, flip_deg=15, duration_us=1000, slice_thickness_mm=5,
                                rephase=False)
ro = sc.modules.CartesianLine(system, fov_ro_mm=250, matrix_ro=64, readout_duration_us=3200,
                              prephase=False)
pe = sc.modules.PhaseEncode(system, fov_pe_mm=250, matrix_pe=64)

seq = sc.LogicBlock('gre')
for i, line in enumerate(range(-32, 32)):
    t0 = i * 20e-3
    seq.add(t0, exc.build(rf_phase_rad=sc.rf_spoil_phase(i)))
    seq.add(t0 + exc.duration, exc.rephaser())          # z  ─┐
    seq.add(t0 + exc.duration, pe.build(line=line))     # y   │ one block,
    seq.add(t0 + exc.duration, ro.prephaser_block())    # x  ─┘ no coordination
    seq.add(t0 + 8e-3 - ro.time_to_echo, ro.build())    # this *is* the definition of TE

out = sc.compile(seq, system)
out.check().raise_if_failed()
out.write('gre.seq')
```

The three winders land in one pulseq block because they coincide, and the compiler says nothing —
they are on different axes, so there is nothing wrong. That is the point of the design: you say
what you mean, and never think about block boundaries.

---

## Install

```console
pip install -e ".[dev,viz]"
```

Requires `pypulseq>=1.5.0` and `numpy>=1.24`. Optional extras: `viz` (matplotlib), `dev` (pytest,
ruff), `docs`.

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
amplitude limit and **189 % of the slew limit**. No module can see that in isolation, so amplitude
and slew are measured on the compiled waveform.

**Errors that name what to change.** pypulseq says `Amplitude violation (117%)`. seqcraft says which
of the three parameters fixes it, and what value would work.

**Physics you can check rather than trust.** The diffusion b-value is verified against numerical
integration of the built waveform to 0.5 %; `k` at the echo is checked on every axis, which is what
catches a missing slice rewinder; and the spiral generator solves the **discrete** amplitude and slew
constraints the hardware actually applies rather than a continuous-time model of them, so it runs at
100 % of peak slew and 99.8 % of it on average. That last one is the difference between a 19 ms
readout and a 29 ms one, and every limit check passed in both cases — which is the point: a sequence
being legal says nothing about it being efficient, and neither is visible without measuring.

**Timing properties that point where they say.** `EPIReadout.time_to_echo` lands 17.26 ms into a
49.92 ms train — the echo whose `ky` is zero, not the train's midpoint, which partial Fourier moves by
7.7 ms. What the ADC must cover is the *sampled* extent, so the readout lobe carries 540.53 1/m to
deliver 533.33 1/m and the prephaser cancels the difference; a prephaser of minus `k_max` would put
every shot 1.7 `dk` off-centre, which is a linear phase ramp across an image that looks perfectly
normal.

**Files that cannot lie about themselves.** `write()` takes no geometry, matrix or FOV arguments —
everything written comes from what was compiled. A JSON sidecar records the versions, the git commit
and dirty flag, the definitions and the file's sha256.

---

## When to reach for it

For one sequence, once, write raw pypulseq — it is a fine tool for that.

seqcraft pays for itself when you have a *family* of sequences, a loop over more than two axes,
gradients that must overlap without you hand-splitting blocks, or files that have to be reproducible
six months later. And when it does not fit, `RawEvents` wraps arbitrary pypulseq events as a module
and `CompiledSequence.seq` is the pypulseq object itself.

---

## Examples

Each folder is one complete scan — build, check, write, simulate, reconstruct, quantify. See
[`examples/README.md`](examples/README.md).

| | What it covers |
|---|---|
| [`01_getting_started.ipynb`](examples/01_getting_started.ipynb) | All three concepts, the overlap rules, the escape hatches, writing a file. |
| [`dti_spiral/`](examples/dti_spiral/) | The single-shot spin-echo **spiral** DTI at 1.88 mm, and the two-echo field map its 67 ms readout cannot do without. b-value against numerical integration, k at the echo on every axis, PNS against the site's own `.asc`. Then simulate and go through to an ADC map. |
| [`dti_epi/`](examples/dti_epi/) | The **same diffusion encoding** through a ramp-sampled **EPI** train — single-shot and two-shot, partial Fourier 0.75, no flat top at all. 20 Hz per pixel of phase-encode bandwidth, so 100 Hz of off-resonance displaces by five pixels — and the same operator that deblurs the spiral removes it. |

Two readouts, one encoding: the diffusivities have to agree, and where they do not, that is the
readout's error budget.

Each build notebook writes into its folder's `seq/`: the `.seq` files, a provenance `.seq.json` for
each, and one sidecar apiece. The `.traj.npz` carries the trajectory at **ADC sample times** with the
echo's own position recorded — neither readout's k-space is in the `.seq` in a form a reconstruction
can use — cross-checked against pypulseq's own calculation to 0.0000 1/m. The field map's `.meta.npz`
carries its echo times and line order.

Play the field map first: it takes five seconds and the diffusion reconstruction needs it.

---

## The module library

```
rf/         SincExcitation  SlabExcitation  HardExcitation  SincRefocusing  HardRefocusing
            SLRExcitation  SLRRefocusing  GaussSaturation  AdiabaticInversion
readout/    CartesianLine  EPIReadout  SpiralVDS  NoiseAcquisition
encoding/   PhaseEncode  PartitionEncode  Prephaser  Spoiler  Crusher
            MonopolarDiffusion  BipolarDiffusion  ArbitraryDiffusion  dti_directions
prep/       FatSat  InversionRecovery
control/    Delay  Trigger  Barrier  RawEvents
```

**There are no recipes.** A recipe is somebody else's sequence choices baked into library code, and
changing your own scan should never mean editing a package. The notebooks build their sequences from
modules, start to finish — the DTI one is about twenty lines of placement, and every number in it is
yours to change.

Writing your own module is one class with an `__init__` and a `build` — see
[`docs/writing_a_module.md`](docs/writing_a_module.md). Then
`sc.testing.assert_all(your_module)` gives it the same contract checks the built-in ones get.

---

## Documentation

- [`docs/architecture.md`](docs/architecture.md) — the three concepts, and what each core module is for.
- [`docs/compiler.md`](docs/compiler.md) — how block boundaries are chosen, how an event's own `delay` is handled, and what every warning means.
- [`docs/writing_a_module.md`](docs/writing_a_module.md) — the three conventions, and what to assert.

---

## Tests

```console
pytest tests --doctest-modules src/seqcraft
```

765 tests and doctests. The compiler directory is the heart of it: one case per rule, plus the
adversarial ones — a gradient straddling an RF, a boundary that would fall inside an ADC window, a
split mid-ramp, two RFs whose dead times overlap. Every physics assertion is a number an independent
calculation gives, not a number the code happened to produce.

---

## A note on vendor data

Siemens `.asc` gradient descriptors carry proprietary PNS and acoustic-resonance coefficients. They
are **never** stored in, copied into or read from this repository: `load_hardware()` resolves them
only through `$SEQCRAFT_ASC_DIR` and rejects anything that looks like a path.
`synthetic_hardware()` provides a vendor-free stand-in so PNS checks can run anywhere — it is not a
real scanner and must never be used to clear a sequence for human scanning.

## Licence

MIT.
