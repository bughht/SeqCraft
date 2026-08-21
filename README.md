# SeqCraft

**Composable, reusable MRI pulse sequence programming on top of [pypulseq](https://github.com/imr-framework/pypulseq).**

![Python](https://img.shields.io/badge/python-3.10%2B-blue) ![License](https://img.shields.io/badge/license-MIT-green) ![Built on](https://img.shields.io/badge/built%20on-pypulseq-orange)

## Pulseq is a list. SeqCraft makes it a tree.

A pulseq sequence is a **flat, sequential list of blocks**, and each block may hold at most **one RF,
one ADC, and one gradient per axis**. That is the right shape for hardware to execute, and the wrong
shape to write in: the moment two things you think of as separate want to happen at the same time — a
slice rephaser and a phase-encode blip, a preparation and the train it prepares, a diffusion lobe and
the refocusing pulse it straddles — they collide in one block. So you split waveforms at boundaries
by hand, keep the pieces in step yourself, and every component you write has to know what its
neighbours are doing.

`sc.LogicBlock` takes that off your side of the work:

- **Anything may overlap anything.** Two gradients on one axis, a gradient across an RF, an ADC while
  three axes play. You write what you mean at the time you mean it; `sc.compile` derives the legal
  pulseq blocks — splitting, summing and bounding — and validates them against the amplifier.
- **Blocks nest, to any depth.** A block holds events *or other blocks*, each with a start time
  relative to its parent. A scan holds shots, a shot holds an inversion and a train, a train holds
  repetitions, a repetition holds an excitation and a readout — all one kind of object. The flat list
  becomes a tree, which is the shape a pulse sequence diagram already has.
- **So a component can be written once and reused.** This is the part that follows from the other
  two: because nothing has to know where block boundaries will fall or what its neighbours are doing,
  a piece of a sequence can be handed around, nested, and retimed. That is what makes `sc.Module`
  possible — a piece of MR physics that takes parameters and returns a `LogicBlock`.

```
   what you write              the model              the compiler            the output
 ─────────────────        ─────────────────      ──────────────────      ─────────────────
  pypulseq events   ─►      sc.LogicBlock    ─►    block boundaries   ─►  pypulseq.Sequence
  your own modules        a tree of events and     same-axis sums          seq.write('x.seq')
  sc.modules.*            blocks, with relative    amplifier limits
                          start times
```

The three sections below are the whole tool: a tree by hand, a module, then a sequence of modules.
Print any of it with `.describe()`, draw it with `sc.plot_block`, measure it with `sc.moments`.

---

## Install

```console
git clone https://github.com/bughht/SeqCraft.git && cd SeqCraft
pip install -e ".[dev,viz]"
```

Python 3.10+, `numpy`, and the pinned `pypulseq` fork (installed automatically). Optional extras:
`viz` (matplotlib), `systems` (vendor limits), `rf` (SLR pulses), `sim` / `recon` (the simulation
notebooks).

> [!IMPORTANT]
> **Set the dead times.** pypulseq defaults `rf_dead_time`, `rf_ringdown_time` and `adc_dead_time`
> to **zero**, which is wrong on every real scanner: the sequence compiles cleanly, validates
> cleanly, and is refused or silently mangled at the console. They belong to your installation, so
> no preset can supply them.

---

## 1. Put pypulseq events on a timeline

Events are made by pypulseq, exactly as you already do it. seqcraft only says *when* each one
plays:

```python
import math
import pypulseq as pp
import seqcraft as sc

opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
               rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)

dk = 1e3 / 220.0                                     # k-space step for a 220 mm FOV, 1/m

rf, gz, gz_reph = pp.make_sinc_pulse(flip_angle=math.radians(15), duration=1e-3,
                                     slice_thickness=5e-3, delay=opts.rf_dead_time,
                                     use='excitation', system=opts, return_gz=True)
gx = pp.make_trapezoid('x', flat_area=64 * dk, flat_time=3.2e-3, system=opts)
adc = pp.make_adc(num_samples=64, duration=3.2e-3, delay=gx.rise_time, system=opts)
gx_pre = pp.make_trapezoid('x', area=-gx.area / 2, duration=1e-3, system=opts)
gy = pp.make_trapezoid('y', area=8 * dk, duration=1e-3, system=opts)     # phase-encode line 8

t_rf_center = rf.delay + pp.calc_rf_center(rf)[0]    # TE is measured from the pulse centre
echo_in_gx = gx.rise_time + gx.flat_time / 2         # and k = 0 sits here inside the readout

tr = sc.LogicBlock('tr')                             # a tree, tagged for error messages
tr.add(0.0, rf, gz)                                  # RF and slice-select: one instant
tr.add(pp.calc_duration(gz), gz_reph, gy, gx_pre)    # z, y, x together — overlap is free
tr.add(t_rf_center + 5e-3 - echo_in_gx, gx, adc)     # TE = 5 ms, by construction

seq = sc.compile(tr, opts)                           # a pypulseq.Sequence — 4 blocks, 7.25 ms
seq.write('tr.seq')                                  # pypulseq's own writer
```

`add` is the only method you need, and it has two shapes for one meaning. `add(t, *items)` — above —
puts these items `t` seconds into this block and returns the block, so calls chain. `add(rows)` takes
the whole schedule as a table of `[time, *items]` rows, which is what you want as soon as the times
are *computed* rather than written out:

```python
>>> table = sc.LogicBlock('tr').add([
...     [0.0,                             rf, gz],
...     [pp.calc_duration(gz),            gz_reph, gy, gx_pre],
...     [t_rf_center + 5e-3 - echo_in_gx, gx, adc],
... ])
>>> [n.start for n in table] == [n.start for n in tr]      # the same tree, built two ways
True
```

Rows go in in the order given and are never sorted, so `nodes` ends up identical either way. And
`tr.duration` is *measured* from the children rather than declared — a block cannot claim a length it
does not play.

Ask the tree what it looks like:

```python
>>> print(tr.describe())
tr  7.25 ms
  +0.0 us  rf
  +0.0 us  trap z
  +1260.0 us  trap z
  +1260.0 us  trap y
  +1260.0 us  trap x
  +4010.0 us  trap x
  +4010.0 us  adc
```

Seven events, three `add` calls, and **no blocks anywhere** — that is the point. The slice rephaser,
the phase blip and the readout prephaser all start at `+1260 µs`; the RF shares an instant with its
slice-select gradient; the ADC runs inside the readout gradient. Written as pulseq, deciding which of
those may share a block and where each waveform has to be cut is your problem. Here you named times,
and `sc.compile` turned the seven events into four legal blocks. `sc.plot_block(tr, opts)` draws the
same thing as a diagram.

---

## 2. Wrap that TR in a module

Section 1 is one repetition, for one phase-encode line. To get the other sixty-three — and to reuse
the whole thing inside a bigger sequence — make it a component: subclass `sc.Module`, design in
`__init__`, assemble in `build`, return one `LogicBlock`.

And do not write those events a second time. Three of the pieces are already modules —
`Excitation`, `PhaseEncode` and `CartesianLine` — and each carries the arithmetic you would
otherwise have to get right twice: the rephaser that follows a selective pulse, one blip designed
once and scaled per line, a prephaser that exactly cancels the readout's ramp. Composing them is
shorter than the events were:

```python
class GRETR(sc.Module):
    """One repetition of a 2D gradient echo, composed from the shipped leaf modules."""

    def __init__(self, *, opts, fov_mm=220.0, matrix=64, thickness_mm=5.0,
                 flip_deg=15.0, te_s=5e-3, bandwidth_hz_px=312.5, tag=None):
        super().__init__(opts=opts, tag=tag)
        self.exc = sc.modules.Excitation(opts=opts, flip_deg=flip_deg,
                                         thickness_mm=thickness_mm, duration_s=1e-3)

        # The blip and the readout prephaser play at the same instant, so the shorter is stretched
        # to match: every leaf reports its own minimum and accepts an override.
        blip = sc.modules.PhaseEncode(opts=opts, fov_mm=fov_mm, matrix=matrix, axis='y')
        read = sc.modules.CartesianLine(opts=opts, fov_mm=fov_mm, matrix=matrix,
                                        bandwidth_hz_px=bandwidth_hz_px)
        winder_s = max(blip.min_duration_s, read.prephaser_duration_s)
        self.pe = sc.modules.PhaseEncode(opts=opts, fov_mm=fov_mm, matrix=matrix, axis='y',
                                         duration_s=winder_s)
        self.ro = sc.modules.CartesianLine(opts=opts, fov_mm=fov_mm, matrix=matrix,
                                           bandwidth_hz_px=bandwidth_hz_px,
                                           prephaser_duration_s=winder_s)

        # Where the readout has to start for k = 0 to land at TE.  Quantised once: a computed start
        # time must sit on the gradient raster, and the compiler raises with this exact fix if not.
        self._read_start_s = sc.Raster(opts.grad_raster_time).ceil(
            self.exc.time_to_center() + te_s - self.ro.time_to_echo())

    def time_to_echo(self) -> float:
        """Seconds from the start of this block to k = 0 — two module times, added."""
        return self._read_start_s + self.ro.time_to_echo()

    def build(self, *, line: int = 0) -> sc.LogicBlock:
        return sc.LogicBlock().add([
            [0.0,                self.exc()],
            [self._read_start_s, self.pe(line=line), self.ro()],     # blip on y, readout on x
        ])
```

Two `add` calls, and no gradient areas, ramp times or dwell arithmetic anywhere: the leaves own
that. Calling the module runs `build`, tags the block with the class name — and the block it returns
is a **tree**, because each leaf contributed a block of its own:

```python
>>> gre_tr = GRETR(opts=opts)
>>> gre_tr(line=40)
LogicBlock(GRETR, 3 nodes, 7.23 ms)
>>> print(gre_tr(line=40).describe())
GRETR  7.23 ms
  +0.0 us  Excitation  1.80 ms
    +0.0 us  rf
    +0.0 us  trap z
    +1260.0 us  trap z
  +3670.0 us  PhaseEncode  0.32 ms
    +0.0 us  trap y
  +3670.0 us  CartesianLine  3.56 ms
    +0.0 us  trap x
    +320.0 us  trap x
    +320.0 us  adc
>>> sc.moments(gre_tr(line=32), order=0)['y']      # the centre line needs no blip
0.0
```

Three children where section 1 had seven events, and the phase blip and the readout prephaser still
land on the same instant — `+3670 µs` — one on y and one on x. Ask for `te_s=3e-3` instead and the
pair slides back to `+1670 µs`, alongside the slice rephaser still playing on z: three axes at once,
five pulseq blocks becoming three, and nothing in the module changed to allow it.

Three conventions make a module reusable, and all three are above:

- **`__init__` designs, `build` assembles.** Waveforms and timings are computed once; sixty-four
  lines are sixty-four cheap calls.
- **Calls are pure.** A call must not mutate the module or its events — `PhaseEncode` derives a
  scaled copy of its blip rather than rescaling one (`pp.scale_grad` plus `sc.events.derive`).
  Mutating in a per-line method still compiles, and makes line 64 differ from line 1 in a way no
  check can see.
- **Declare no duration, declare no position.** The block measures itself, and a module that must
  say *when* something happens inside it exposes a time instead: `time_to_echo()` is the one
  question the tree cannot answer, because a tree knows when its events play but not which instant
  among them is the echo.

Check any module in isolation with `sc.compile(sc.LogicBlock('probe').add(0.0, gre_tr(line=40)),
opts)`: a component that only works when something else happens to be beside it is not reusable.
`sc.modules.GRE2DTR` is this same composition with spoiler gradients, a phase-encode rewinder and an
RF-phase argument added — writing one out like this is how that one was written.

---

## 3. Build the sequence from modules

A block may hold blocks, which may hold blocks. So the same `add` builds every level: an inversion
and the train that follows it make a **shot**, and a handful of shots at different inversion times
make a T1-mapping **scan**.

```python
gre_tr = GRETR(opts=opts)                                                   # from section 2
inv = sc.modules.IRPrep(opts=opts, thickness_mm=None, spoil_voxel_mm=5.0)   # ships with seqcraft

tr_s = 12e-3
raster = sc.Raster(opts.grad_raster_time)
lines = sorted(range(64), key=lambda k: (abs(k - 32), k))    # centric: k = 0 acquired first

def shot(ti_s, seg=lines):
    """One inversion and the train recovering into it — itself just a LogicBlock."""
    # TI runs from the inversion's effective centre to the acquisition of k = 0.  Both ends are
    # module questions; the subtraction is the whole layout.
    t_train = raster.ceil(inv.time_to_center() + ti_s - gre_tr.time_to_echo())

    rows = [[0.0, inv()]]                                    # the inversion, then the train
    rows += [[t_train + i * tr_s, gre_tr(line=k)] for i, k in enumerate(seg)]
    return sc.LogicBlock('shot').add(rows)

tis = (100e-3, 300e-3, 700e-3, 1500e-3)
scan = sc.LogicBlock('ir_t1').add([[i * 4.0, shot(ti)] for i, ti in enumerate(tis)])   # 4 s apart

seq = sc.compile(scan, opts, name='ir_t1')   # 1547 blocks, 14.263 s
seq.write('ir_t1.seq')
```

`inv.time_to_center()` is 5.101 ms into its own block — a 10 ms hyperbolic secant inverts at its
centre, not at its start, and referencing TI to the block start would be a 5 ms error in the one
quantity the sequence exists to control. Neither that nor `gre_tr.time_to_echo()` is measurable from
the tree, which is exactly the division of labour: **modules know the physics, the tree knows the
times, the compiler knows pulseq.**

### The tree is the sequence diagram

Five levels, and every one of them is the same kind of object. Here is one shot cut down to two
lines so it fits on the page — the real one is the same shape, 64 repetitions wide:

```python
>>> print(sc.LogicBlock('ir_t1').add(0.0, shot(300e-3, lines[:2])).describe())
ir_t1  318.70 ms
  +0.0 us  shot  318.70 ms
    +0.0 us  IRPrep  11.34 ms
      +0.0 us  rf
      +10130.0 us  spoiler  1.21 ms
        +0.0 us  trap z
    +299470.0 us  GRETR  7.23 ms
      +0.0 us  Excitation  1.80 ms
        +0.0 us  rf
        +0.0 us  trap z
        +1260.0 us  trap z
      +3670.0 us  PhaseEncode  0.32 ms
        +0.0 us  trap y
      +3670.0 us  CartesianLine  3.56 ms
        +0.0 us  trap x
        +320.0 us  trap x
        +320.0 us  adc
    +311470.0 us  GRETR  7.23 ms
      +0.0 us  Excitation  1.80 ms
        +0.0 us  rf
        +0.0 us  trap z
        +1260.0 us  trap z
      +3670.0 us  PhaseEncode  0.32 ms
        +0.0 us  trap y
      +3670.0 us  CartesianLine  3.56 ms
        +0.0 us  trap x
        +320.0 us  trap x
        +320.0 us  adc
```

`trap x` sits inside `CartesianLine`, inside `GRETR`, inside `shot`, inside the scan — and every
offset is read against its own parent, never against the scan. That is what makes a subtree
portable: the shot does not know it is the second one, so moving it moves everything it contains.

```python
>>> len(scan), sum(1 for _ in sc.flatten(scan))   # four children; 1800 leaf events
(4, 1800)
>>> scan.nodes[2].start += 20e-3                  # 450 events later, from one number
>>> sorted({path for _, _, path in sc.flatten(scan)})
[('ir_t1', 'shot', 'GRETR', 'CartesianLine'), ('ir_t1', 'shot', 'GRETR', 'Excitation'), ('ir_t1', 'shot', 'GRETR', 'PhaseEncode'), ('ir_t1', 'shot', 'IRPrep'), ('ir_t1', 'shot', 'IRPrep', 'spoiler')]
```

Those paths are provenance, and nobody wrote them: they are the tags on the way down, and they are
what every warning and error message names. The rest of the tree is plain Python — `scan.nodes` is a
list, `sc.flatten` walks it, `lb.copy()` gives you a variant to retime, and adding one block object
at sixty-four times shares it rather than copying it.

Nothing in here coordinates. `inv` does not know a train follows it, `gre_tr` does not know what
preceded it, `shot` does not know it is one of four, and none of them knows where pulseq's block
boundaries will fall.

### The shipped modules

`sc.modules` has nine building blocks, each extracted from a working sequence rather than designed:
`Excitation`, `Refocusing`, `PhaseEncode`, `CartesianLine`, `EPI2D`, `spoiler`, `IRPrep`, `GRE2DTR`
and `GRE2D`. The last two are a whole repetition and a whole scan, so the GRE that section 2
composed also comes ready-made:

```python
gre = sc.modules.GRE2D(opts=opts, fov_mm=220.0, matrix=(64, 64), thickness_mm=5.0)
seq = sc.compile(gre(lines=range(64)), opts, name='gre_2d')      # 256 blocks, TE 4.62 ms
```

`GRE2D` takes the **list of phase-encode lines** rather than an acceleration factor: which lines to
acquire is a sequence-programming choice, and it stays yours. Segmenting the train across several
inversions — an MPRAGE — is the same tree one level deeper, and
[`examples/mprage_2d/`](examples/mprage_2d/) builds it.

`Refocusing` is the one that shows what a module is *for*. A refocusing pulse
conjugates k, so between consecutive refocusing centres every axis's gradient area before the echo
has to equal its area after it — measured to the RF's **effective centre**, not to the middle of the
block. Getting that wrong leaves a residual that alternates sign echo to echo and reads as a
hardware fault, and both pulseq reference implementations avoid it only by setting their transmit
dead time and ringdown to the same number:

```python
refoc = sc.modules.Refocusing(opts=opts, thickness_mm=6.25, crush_voxel_mm=5.0)

assert refoc.area_to_center_per_m == refoc.area_from_center_per_m    # 600.000000 both
assert refoc.time_to_center() == refoc().duration / 2                # exactly, to 0 ns
```

A spin echo and a sixteen-echo turbo spin echo are then the same composition with a longer list,
and [`examples/fse_2d/`](examples/fse_2d/) writes both.

`EPI2D` is the newest, and the one where the arithmetic is the module. It is the whole of k-space in
one shot — alternating readout lobes, blips on the zero crossings, one ADC per echo — and an EPI
ghosts because those lobes alternate. Two rules follow. The sampling window is **exactly** centred
in its lobe, because only then do the two polarities sample one k grid; write the lobe the way a
single *line* is correctly written and they sample two grids half a gradient raster apart, which is
an N/2 ghost the sequence made itself. And the blip sits on the readout's zero crossing with
`sc.barrier()` pinned there, because otherwise the compiler splits every readout lobe in the train
and the only report of it is a merge warning:

```python
epi = sc.modules.EPI2D(opts=opts, fov_mm=220.0, matrix=(128, 128), dwell_s=2.5e-6)

assert 2 * epi.guard_s + epi.num_samples * epi.dwell_s == epi.echo_spacing_s   # exactly
assert epi.k_read_per_m[epi.echo_sample(0)] == 0.0     # k = 0 is a sample, not a moment
```

Parallel imaging then needs nothing added: `blip_lines=R` is the acceleration and `lines` is the
table. The calibration band is more `build` calls too — but with a table of **length one**, which
is a Cartesian gradient echo on the EPI's own readout lobe rather than another EPI shot.

**Segmentation** needs exactly one thing: `phase_deg`, the receiver phase, because a spoiled EPI
needs its receiver locked to a transmitter that is advancing its carrier. Four shots that differ
from one another for *any* reason put a replica of the object at `Ny/4`, and it reads as an
under-sampling artefact.

---

## What the compiler checks

`sc.compile(tree, opts)` returns a `pypulseq.Sequence` and nothing else. If the tree cannot become a
legal sequence it **raises**; if it had to change a waveform to make it legal it **warns**. There is
no report object to unpack, so there is nothing to forget to read.

- **Block boundaries and overlap.** Different axes at one instant become one block. Two gradients on
  the *same* axis are summed, with a warning naming both sources.
- **Limits measured on the compiled waveform**, which is the only place the truth is: two
  individually legal gradients can sum to an illegal one, and no component can see that alone.
- **RF and ADC conflicts**, including when only their dead times touch — which pypulseq would
  otherwise reject 40 000 blocks later.
- **One informational warning you should expect.** Three axes ramping together exceed the
  *vector-norm* slew bound that per-axis limits imply, routinely and legally, so seqcraft reports it
  rather than raising.

Errors name the offending number, where it happened, the tag path it came from, and what to change:

```
HardwareLimitError: slew 189% of the 150 T/m/s limit on axis x.
  from   :  probe.a, probe.b
  at     :  0.000 ms (block 0)
  reached:  284.0 T/m/s
  fix
    lengthen the lobe, or lower the readout bandwidth
    or design that part against sc.opts.derate(opts, slew=0.52)
```

`sc.moments`, `sc.kspace`, `sc.sample` and `sc.pns` measure a tree directly, before any file is
written.

---

## Where to look next

Notebooks, each one a sequence that works rather than a feature tour:

- [`examples/01_getting_started.ipynb`](examples/01_getting_started.ipynb) — blocks, `Opts` and
  `compile`, the overlap rules and the escape hatches. Uses no modules.
- [`examples/gre_2d/`](examples/gre_2d/) — a spoiled 2D GRE three ways, then simulated and
  reconstructed. Six of the nine shipped modules came out of it.
- [`examples/mprage_2d/`](examples/mprage_2d/) — segmented and inversion-prepared, with the null
  point checked in simulation.
- [`examples/mp2rage_2d/`](examples/mp2rage_2d/) — two trains, the `SET` label that separates them,
  and the ratio that cancels the receive field.
- [`examples/se_2d/`](examples/se_2d/) — a spin echo, the area balance that makes it one, and the
  measurement that says the echo is a **T2** echo rather than a T2\* one. `Refocusing` came out of it.
- [`examples/fse_2d/`](examples/fse_2d/) — the same composition at sixteen echoes and then
  seventy-two: 4.3 minutes becomes 18 seconds, and what that costs is measured —
  including a ghost that every arithmetic check passes through.
- [`examples/gre_epi_2d/`](examples/gre_epi_2d/) — the whole of k-space in one shot, and the two
  rules that stop it ghosting on its own. Both are built the **wrong** way too, because both
  failures compile. `EPI2D` came out of it.
- [`examples/se_epi_2d/`](examples/se_epi_2d/) — the same readout after a refocusing pulse, and the
  measurement most readers expect to come out the other way: **a spin echo does not fix EPI
  distortion.**

All seven simulation notebooks share one phantom — [`examples/phantom.py`](examples/phantom.py) —
including the off-resonance map the EPI examples distort against, which is the phantom's own.

Documentation: [`api_reference.md`](docs/api_reference.md) (every public name, executed by CI),
[`architecture.md`](docs/architecture.md) (the layering, and what is deliberately absent),
[`compiler.md`](docs/compiler.md) (how boundaries are chosen, what every warning means),
[`writing_a_module.md`](docs/writing_a_module.md) (the `Module` contract in full).

Tests: `pytest tests --doctest-modules src/seqcraft`.

---

## Notes

**When not to reach for it.** For one sequence, once, raw pypulseq is a fine tool. seqcraft pays for
itself on a *family* of sequences, gradients that must overlap without hand-splitting blocks, or
files that have to be reproducible six months later — and when it does not fit, what you are holding
is already the pypulseq object.

**Vendor data stays out of this repository.** Siemens `.asc` descriptors carry proprietary
coefficients, so `sc.hardware.load_hardware()` reads them only through `$SEQCRAFT_ASC_DIR`.
`sc.hardware.synthetic_hardware()` is a vendor-free stand-in for PNS checks — not a real scanner,
and never to be used to clear a sequence for human scanning.

## Licence

MIT.
