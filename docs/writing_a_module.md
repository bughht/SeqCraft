# Writing a component

## What seqcraft actually requires

One thing: whatever takes part in a sequence returns a `LogicBlock`. There is no base class you must
inherit, no method name to match, no registry to join and no hook to override.

A function is a component:

```python
import pypulseq as pp
import seqcraft as sc


def spoiler(opts, *, twists=4, voxel_mm=5.0, axis='z'):
    """Dephase by `twists` cycles across a voxel."""
    area_per_m = twists / (voxel_mm / 1e3)
    g = pp.make_trapezoid(channel=axis, area=area_per_m, system=opts)
    return sc.LogicBlock('spoil').add(0.0, g)
```

So is a class of whatever shape the physics wants. Nothing here inherits anything, and the compiler
cannot tell:

```python
class VelocityEncode:
    """A bipolar pair straddling a refocusing pulse — two outputs, named for what they are."""

    def __init__(self, opts, *, venc_cm_s, axis='y'):
        self.m1_s_per_m = 1.0 / (2.0 * venc_cm_s / 100.0)
        self.lobe = pp.make_trapezoid(channel=axis, area=self._area(), system=opts)

    def pre(self):
        return sc.LogicBlock('venc_pre').add(0.0, self.lobe)

    def post(self):
        return sc.LogicBlock('venc_post').add(0.0, pp.scale_grad(self.lobe, -1.0))
```

```python
venc = VelocityEncode(opts, venc_cm_s=50)
first = venc.pre()
seq.add(t0, first)
seq.add(t0 + first.duration + refoc_duration, venc.post())

sc.compile(sc.LogicBlock('probe').add(0.0, venc.pre()), opts)   # it must compile on its own
```

Note the placement: `first.duration`, read off the block that was just built, rather than a
`lobe_duration` property. Building is cheap, and a measured number cannot disagree with what plays.

---

## `sc.Module` — the shape for something you will reuse

For a component you will use in several sequences, `sc.Module` is the standard shape: **parameters
in, one `LogicBlock` out.** Four members, and each earns its place.

```python
class Module(ABC):
    def __init__(self, *, opts: pp.Opts, tag: str | None = None) -> None
    def __call__(self, *args, **kwargs) -> LogicBlock       # the interface
    @abstractmethod
    def build(self, *args, **kwargs) -> LogicBlock          # what you write
```

**`opts` is required**, and it is the official `pypulseq.Opts`. It is the only scanner input a
module gets: rasters, dead times, ringdown, gradient and B1 limits, gamma, B0, sample limits. It is
required rather than defaulted because pypulseq's fallback is the *process-global* `Opts.default`,
which makes a sequence depend on import order. Pass it explicitly, always, including down to
submodules.

**`tag` is identity.** Tags become the provenance path `flatten` builds and the `from …` clause in
every compiler warning. Two instances of one class doing different jobs — two readouts, two spoilers
— are told apart by it. It defaults to the class name.

**`__call__` is the interface; `build` is what you write.** `module(...)` reads correctly at the
call site, next to `add`: `tr.add(t, readout(line=17))`. It is also the single place the framework
gets to check and name what came back, which is why `build` is not the public entry point.

```python
import pypulseq as pp
import seqcraft as sc


class VelocityEncode(sc.Module):
    """A bipolar gradient pair that encodes velocity along one axis."""

    def __init__(self, *, opts, venc_cm_s, axis='y', tag=None):
        super().__init__(opts=opts, tag=tag)
        # A phase of pi at the target velocity needs m1 = 1 / (2 v), in s/m.
        self.m1_s_per_m = 1.0 / (2.0 * float(venc_cm_s) / 100.0)
        self.lobe = pp.make_trapezoid(channel=axis, area=self._area(), system=opts)

    def build(self, *, sign=1.0) -> sc.LogicBlock:
        first = sc.events.derive(pp.scale_grad(self.lobe, sign))
        second = sc.events.derive(pp.scale_grad(self.lobe, -sign))
        return (sc.LogicBlock()
                .add(0.0, first)
                .add(float(pp.calc_duration(first)), second))
```

`pp.scale_grad` is the one to reach for, not a hand-written `derive(amplitude=…, area=…,
flat_area=…)`. A trapezoid carries three numbers that must move together and an arbitrary gradient
carries four (`waveform`, `first`, `last`, `area`); writing them out is that many chances to update
all but one, and the result still compiles — it just encodes a different line from the one it
reports. Wrapping it in `sc.events.derive` strips pypulseq's registration state from the copy, so
the stored design can be scaled again after a compile has registered one of its outputs.

```python
venc = VelocityEncode(opts=opts, venc_cm_s=50)
seq.add(t0, venc(sign=-1.0))
```

### The two failures the base catches

```python
class Broken(sc.Module): pass
Broken(opts=opts)
# TypeError: Can't instantiate abstract class Broken with abstract method build

class Wrong(sc.Module):
    def build(self): return pp.make_delay(1e-3)
Wrong(opts=opts)()
# TypeError: Wrong.build() returned SimpleNamespace, not a LogicBlock
```

The second is the one that matters. Without it, returning the wrong thing surfaces hundreds of lines
later inside `add`, which can only report that it was handed a `SimpleNamespace` — not which module
handed it over.

### What the base deliberately does not do

| Not there | Where the job goes |
|---|---|
| a `duration` property | `LogicBlock.duration`, measured from the block your call returned |
| a unit check | your `__init__`, where the plausible range is actually known |
| `params()` / provenance walking | a provenance writer that takes a mapping; nothing scrapes `__dict__` |
| `submodules()` | `vars(self)`, if anything ever needs it |
| a registry | nothing needs a string → class lookup |
| a scanner class, named limit regimes | `pp.Opts`. A part designed derated takes a second `Opts` from `sc.opts.derate` — an object, not a name to look up |
| a fixed timing vocabulary (`isodelay`, `time_to_echo`) | the individual module, where the physics is known |

Each is a real feature; none is required by *every* module, and a base class that carries what only
some subclasses need is how a small idea becomes a framework.

---

## The conventions that are worth following

None is enforced — the compiler never sees a module at all — but each exists because getting it
wrong produced a specific bug.

### `__init__` designs, `build` assembles

Waveforms are created in `__init__` and stored on `self`, the way `nn.Conv2d` allocates its weights
there. A call is cheap and returns a value, so a 30-direction diffusion encoding costs **one**
design and thirty cheap calls.

### Calls are pure

A call must not mutate the module or the events stored on it. It is called once per TR; a module
that mutates itself makes TR 500 differ from TR 1, and the difference is usually a sign flip that
produces a plausible but wrong image. Derive modified events with `sc.events.derive` — never assign
to them.

> `sc.events.derive` shallow-copies and strips pypulseq's registration state. Never `deepcopy` an
> event: `set_block` caches on it keyed by `id(sequence)`, and a stale key can silently
> mis-attribute the event to a different `Sequence` that lands at the same address.

### One call returns one block

Including when the parts are far apart in time. A diffusion encoding is one block with a hole in the
middle, and the caller drops the refocusing pulse into the hole. Overlap in the tree costs nothing,
so nothing has to be reserved and no artificial duration is added.

### A build argument *may* change the duration

This rule used to be the opposite, and the inversion is the point. There is no declared `duration`
any more, so there is no second number a build argument can invalidate. Build the block, then read
it:

```python
exc_block = exc()
tr.add(0.0, exc_block)
tr.add(exc_block.duration + gap, readout(line=k))
```

Building is cheap — the design happened in `__init__` — so nothing is lost by building before
placing. What is gained is that a call argument is now free to change the duration, which the
declared-duration design had to forbid: when `CartesianLine.build` accepted `prephase=False` it
silently invalidated `time_to_echo` and moved the echo 280 µs early.

### Semantic offsets are the module's, and are not standardised

`duration` is measurable and therefore belongs to the block. **Where k = 0 falls inside a readout is
not.** The block knows times, not meanings, and pinning the answer into it would stop the same block
being reused elsewhere — which is why `LogicBlock` has no marks or anchors.

So it belongs to whatever module knows the physics, expressed however that module's domain wants.
Where a build argument moves the echo, the query takes that argument too:

```python
ro.time_to_echo(shot=1)
```

The base standardises no names for these. What a name has to be honest about:

- **`isodelay`** — start of the block to the RF's effective centre. Read it from `rf.center`, not by
  recomputing a centre of mass: the two differ for asymmetric and minimum-phase pulses, and the
  designed value is the one timing needs.
- **`time_to_echo`** — start of the block to k = 0. **For a spiral that is the first sample**, not
  the middle of the window; conflating the two shifts the diffusion weighting rather than the image.
  For an EPI train it is neither: with partial Fourier 0.75 the `ky = 0` echo lands 17.2 ms into a
  49.9 ms train, so both the midpoint and the first sample are wrong, and by different amounts.

---

## Where an ADC goes, and on which raster

Three rasters constrain an ADC, and each one is a bug that compiles.

**Put the sampling offset in the ADC's own `delay`, and place the node where the gradient starts.**

```python
adc = pp.make_adc(num_samples=n, dwell=dwell, delay=offset, system=self.opts)   # right
out.add(lobe_start, adc)

adc = pp.make_adc(num_samples=n, dwell=dwell, delay=0.0, system=self.opts)      # wrong, twice
out.add(lobe_start + offset, adc)
```

The second form fails in two independent ways:

- **`pp.make_adc` silently raises a `delay` below `adc_dead_time` up to it**, and seqcraft
  *preserves* an event's own delay rather than folding it away. So `delay=0.0` becomes 10 µs, the
  node offset adds to it, and sampling begins 40 µs into the lobe where 30 was intended. A spiral's
  `time_to_echo` must therefore return `float(self.adc.delay)`, not `0.0`.
- **An RF or ADC node off the block raster is snapped**, with a `raster` warning rather than an
  error. The snap moves the event against the gradient by up to half a raster: at 2 MHz/m that is
  10 1/m of k, two and a half `dk`. Placing the node at the gradient's own start puts it on the
  raster by construction, and leaves the offset free.

The offset itself must land on the **RF** raster (1 µs), which is what pypulseq's own `check_timing`
requires of an ADC `delay` — not the 100 ns ADC raster, which only the *dwell* answers to. A 33.2 µs
offset produces one error per echo.

---

## Nesting needs no mechanism

A module that holds modules just holds them, passes `opts` down, and calls them:

```python
class FatSat(sc.Module):
    def __init__(self, *, opts, voxel_mm, tag=None):
        super().__init__(opts=opts, tag=tag)
        self.pulse = GaussSaturation(opts=opts, flip_deg=90, duration_us=8000)
        self.spoiler = Spoiler(opts=opts, twists=4, voxel_mm=voxel_mm)      # opts passed down

    def build(self) -> sc.LogicBlock:
        pulse = self.pulse()                                # place by what was just built
        return sc.LogicBlock().add(0.0, pulse).add(pulse.duration, self.spoiler())
```

`sc.flatten` then reports the whole path — `('gre', 'FatSat', 'Spoiler')` — so provenance is the tree
and needs no bookkeeping. **Not one tag string was written**, and every compiled block still traces
back to the module that produced it.

---

## Overlap is free

Do not arrange your module so its gradients avoid other modules'. Place them where they belong and
let the compiler deal with it:

```python
seq.add(t, gz_reph)                 # z
seq.add(t, pe(line=k))              # y   } one block, no coordination,
seq.add(t, prephaser())             # x   } and no warning
```

Different axes are silent. The same axis is summed with a warning naming both sources. Two RF or two
ADC events overlapping is an error. See [`compiler.md`](compiler.md).

---

## Units, and where the check went

The base no longer runs one, and there is no helper in the package to call either. The plausibility
bands keyed off the unit suffix — what turns `fov_mm=0.22` into *"0.22 looks like m. Did you mean
fov_mm=220?"* — left with `Geometry`, and are in [`salvage/geometry.py`](../salvage/geometry.py).
They are worth adopting when a module library grows the infrastructure to share them; copying them
into your own package is the intended use.

Until then: state the range where the range is actually known.

```python
if not 0.5 <= fov_mm <= 2000:
    raise ConfigurationError(...)
```

The convention that makes such a check possible is worth keeping either way: public parameters use
researcher-natural units (`fov_mm`, `te_ms`, `flip_deg`, `b_value_s_per_mm2`); anything derived is
strict SI with an SI suffix (`fov_m`, `te_s`, `flip_rad`). Convert exactly once, and never
auto-convert — a wrong unit is a wrong sequence, and silently fixing it hides the mistake.

---

## Errors should name what to change

When you refuse a configuration, say which parameter fixes it and what value would work:

```python
msg = format_error(
    f'a {flip:g} degree pulse in {duration:g} us needs {needed:.1f} uT, '
    f'above the {limit:.1f} uT limit.',
    {'flip_deg': flip, 'duration_us': duration, 'max_b1_uT': limit},
    [
        f'lengthen the pulse to at least {shortest:.0f} us',
        'reduce flip_deg',
        'use an adiabatic pulse, which inverts at lower peak B1',
    ],
)
raise ConfigurationError(msg)
```

pypulseq's own message for that case is `Amplitude violation (117%)`, which does not say which of the
three parameters to change. That difference is most of what a module is for.

---

## What to assert in your tests

**Compile it on its own.** That is one line, and it is most of the suite seqcraft used to ship:

```python
def test_the_lobe_compiles_alone(opts):
    sc.compile(sc.LogicBlock('probe').add(0.0, venc.pre()), opts)
```

A component that only works when something else happens to be beside it is not reusable, and the
compile checks everything a separate assertion could and more. An off-raster gradient start raises
`CompileError` naming the nearest raster point above and below; a lobe over the amplifier's slew
limit raises `HardwareLimitError` with the derating that would clear it — measured on the *summed*
waveform, which is the only place two individually legal gradients on one axis can be seen to reach
189 % of the limit together.

**Then check purity, because the compiler structurally cannot.** It validates a tree. It never sees
the second call. A component that designs once and assembles per TR is called once per TR, so
self-mutation, accumulation and nondeterminism are invisible to `sc.compile` by construction — each
individual call hands it a tree that is perfectly legal. The reference implementation's
`self.gx.amplitude = -self.gx.amplitude` inside a readout loop compiles cleanly every TR and
produces a plausible but wrong image.

```python
def test_my_module_is_pure():
    pe = PhaseEncode(opts=opts, fov_m=0.22, matrix=64, duration_s=600e-6)
    before = {k: sc.events.content_hash(v)
              for k, v in vars(pe).items() if hasattr(v, 'type')}
    pe(line=17); pe(line=17)
    assert {k: sc.events.content_hash(v)
            for k, v in vars(pe).items() if hasattr(v, 'type')} == before
```

Twice, not once. The canonical bug is an *involution*: comparing only before the first call and
after the second finds the module exactly where it started and reports nothing.

There is no duration check, and deliberately so: a module declares no duration, so there is no
second number that can disagree with the first, and nothing left to assert.

Add the **known values** yourself, because those are the ones that catch physics errors:

```python
def test_venc_encodes_the_right_velocity(opts):
    venc = VelocityEncode(opts=opts, venc_cm_s=50)
    raster = opts.grad_raster_time
    m1 = 0.0
    for node in venc():
        tt, wf = sc.events.waveform_of(node.item, raster)
        m1 += float(sc.events.trapz(wf * (tt + node.start), tt))
    assert m1 == pytest.approx(venc.m1_s_per_m, rel=1e-3)
```

**Measure from the waveform, not from the parameter.** The deleted diffusion module's b-value was
checked that way against numerical integration to 0.5 %, and it is what caught a published ramp
correction being used with the wrong `delta` convention — a 2–4 % b-value error that would have
biased every diffusivity by the same amount, invisibly. (That integral survives, as a plain
function, in [`salvage/bvalue.py`](../salvage/bvalue.py).)

Two more that are worth writing for any module with a refocusing pulse:

- **Where does k end up at the echo?** `out.kspace()['k_adc'][:, 0]` should be zero on every axis you
  meant to refocus. This is what catches a missing slice rewinder: a spin-echo spiral starts at k=0
  in x and y, which makes it tempting to think nothing needs rephasing on z — but the excitation's
  slice-select gradient leaves through-slice dephasing equal to its own tail, and the refocusing
  pulse does *not* undo it. Without the rewinder `k_z` is −525 1/m, which is 2.1 cycles across a 4 mm
  slice and about 95 % of the signal, and nothing else in the sequence looks wrong.
- **Does the scan take as long as it should?** TR is the time between exciting the *same* slice, so
  every slice belongs inside one TR period. Getting that wrong does not produce a wrong image — it
  produces a correct one that takes `n_slices` times longer.

Use `calculate_kspacePP` as the oracle for the first, rather than writing your own integrator: it
handles the refocusing conjugation itself.
