# Writing a component

## What seqcraft actually requires

One thing: whatever takes part in a sequence returns a `LogicBlock`. There is no base class to
inherit, no method name to match, no registry to join and no hook to override.

A function is a component:

```python
import pypulseq as pp
import seqcraft as sc


def spoiler(system, *, twists=4, voxel_mm=5.0, axis='z'):
    """Dephase by `twists` cycles across a voxel."""
    area_per_m = twists / (voxel_mm / 1e3)
    g = pp.make_trapezoid(channel=axis, area=area_per_m, system=system.default)
    return sc.LogicBlock('spoil').add(0.0, g)
```

So is a class of whatever shape the physics wants. Nothing here inherits anything, and the compiler
cannot tell:

```python
class VelocityEncode:
    """A bipolar pair straddling a refocusing pulse — two outputs, named for what they are."""

    def __init__(self, system, *, venc_cm_s, axis='y'):
        self.m1_s_per_m = 1.0 / (2.0 * venc_cm_s / 100.0)
        self.lobe = pp.make_trapezoid(channel=axis, area=self._area(), system=system.default)

    @property
    def lobe_duration(self):
        return float(pp.calc_duration(self.lobe))

    def pre(self):
        return sc.LogicBlock('venc_pre').add(0.0, self.lobe)

    def post(self):
        return sc.LogicBlock('venc_post').add(0.0, _scaled(self.lobe, -1.0))
```

```python
venc = VelocityEncode(system, venc_cm_s=50)
seq.add(t0, venc.pre())
seq.add(t0 + venc.lobe_duration + refoc.duration, venc.post())

sc.testing.assert_output(venc.pre, system)     # the block-level contract, for any callable
```

The old shape for this was one `build(part='pre')` with a keyword that meant two different
gradients. Where that reads better, keep it; where two methods read better, write two methods.

---

## The optional base

`sc.modules.Module` is a convenience base, and what the built-in library is written on. It holds the
scanner, resolves the limit regime to design against, checks units when your `__init__` returns, and
reports parameters for the provenance sidecar and `repr`. It declares no abstract method, so nothing
is required of a subclass — inheriting it buys you those four things and asks for nothing.

```python
import pypulseq as pp
import seqcraft as sc


class VelocityEncode(sc.modules.Module):
    """A bipolar gradient pair that encodes velocity along one axis."""

    def __init__(self, system, *, venc_cm_s, axis='y', regime='default'):
        super().__init__(system, regime=regime)
        self.venc_cm_s = float(venc_cm_s)
        self.axis = axis
        # A phase of pi at the target velocity needs m1 = 1 / (2 v), in s/m.
        self.m1_s_per_m = 1.0 / (2.0 * self.venc_cm_s / 100.0)
        self.lobe = pp.make_trapezoid(channel=axis, area=self._area(), system=self.opts)

    @property
    def duration(self):
        return 2.0 * float(pp.calc_duration(self.lobe))

    def build(self, *, sign=1.0):
        first = _scaled(self.lobe, sign)
        second = _scaled(self.lobe, -sign)
        out = sc.LogicBlock('venc')
        out.add(0.0, first)
        out.add(float(pp.calc_duration(first)), second)
        return out


def _scaled(grad, factor):
    return sc.events.derive(
        grad,
        amplitude=float(grad.amplitude) * factor,
        area=float(grad.area) * factor,
        flat_area=float(grad.flat_area) * factor,
    )
```

Then:

```python
venc = VelocityEncode(system, venc_cm_s=50)
sc.testing.assert_all(venc)      # the same contract checks the built-in modules get
```

---

## The three conventions the library follows

These are how the built-in modules are written and why. None of them is enforced — `Module` has no
abstract method and the compiler never sees a module at all — but each one exists because getting it
wrong produced a specific bug, so they are worth knowing before departing from them.

### `__init__` designs, `build` assembles

Waveforms are created in `__init__` and stored on `self`, the way `nn.Conv2d` allocates its weights
there. `build` is cheap and returns a value, so a 30-direction diffusion encoding costs **one**
design and thirty cheap builds.

`build` must not mutate the module or the events stored on it. It is called once per TR; a module
that mutates itself makes TR 500 differ from TR 1, and the difference is usually a sign flip that
produces a plausible but wrong image. Derive modified events with `sc.events.derive` — never assign
to them.

> `sc.events.derive` shallow-copies and strips pypulseq's registration state. Never `deepcopy` an
> event: `set_block` caches on it keyed by `id(sequence)`, and a stale key can silently
> mis-attribute the event to a different `Sequence` that lands at the same address.

### `build`'s arguments select the variant

`part='pre'`, `line=17`, `interleaf=3`, `slice_offset_m=z`. Ordinary keyword arguments with ordinary
defaults, and Python reports a typo as a `TypeError`.

**A build argument must not change the block's duration.** That is the one hard rule, and it exists
because callers place the *next* thing using the module's timing properties — which cannot see the
build argument. When `CartesianLine.build` accepted `prephase=False` it silently invalidated
`time_to_echo` and moved the echo 280 µs early. Anything that changes the duration belongs on
`__init__`:

```python
ro = sc.modules.CartesianLine(system, ..., prephase=False)   # right
ro.build(prephase=False)                                     # was wrong; no longer exists
```

The rule is about *timing properties*, not about the word `build`. A component whose outputs have no
shared duration to invalidate — `venc.pre()` and `venc.post()` above — does not have this problem,
and a separate method per output is often the clearer way to say so.

### Timing a caller needs is a property

`exc.isodelay`, `readout.time_to_echo`, `diff.lobe_duration`. The module has the domain knowledge, so
the module answers.

The block does not, and cannot: it does not know where it sits. Pinning the answer into it would also
stop the same block being reused elsewhere — which is why `LogicBlock` has no marks or anchors.

A property has to be honest about what it points at:

- **`duration`** — at least as long as the block `build` returns. It may exceed it (padding to the
  raster); it may never fall short, or whatever is placed next silently overlaps.
- **`isodelay`** — start of the block to the RF's effective centre. Read it from `rf.center`, not by
  recomputing a centre of mass: the two differ for asymmetric and minimum-phase pulses, and the
  designed value is the one timing needs.
- **`time_to_echo`** — start of the block to k=0. **For a spiral that is the first sample**, not the
  middle of the window; conflating the two shifts the diffusion weighting rather than the image. For
  an EPI train it is neither: with partial Fourier 0.75 the `ky = 0` echo lands 17.2 ms into a 49.9 ms
  train, so both the midpoint and the first sample are wrong, and by different amounts.

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
  *preserves* an event's own delay rather than folding it away. So `delay=0.0` becomes 10 µs, the node
  offset adds to it, and sampling begins 40 µs into the lobe where 30 was intended.
  `SpiralVDS.time_to_echo` returns `float(self.adc.delay)` rather than `0.0` for exactly this reason.
- **An RF or ADC node off the block raster is snapped**, with a `raster` warning rather than an error.
  The snap moves the event against the gradient by up to half a raster: at 2 MHz/m that is 10 1/m of
  k, two and a half `dk`. Placing the node at the gradient's own start puts it on the raster by
  construction, and leaves the offset free.

The offset itself must land on the **RF** raster (1 µs), which is what pypulseq's own `check_timing`
requires of an ADC `delay` — not the 100 ns ADC raster, which only the *dwell* answers to. A 33.2 µs
offset produces one error per echo.

---

## Nesting needs no mechanism

A module that contains modules just nests their blocks:

```python
class FatSat(sc.modules.Module):
    def __init__(self, system, *, voxel_mm):
        super().__init__(system)
        self.pulse = sc.modules.GaussSaturation(system, flip_deg=90, duration_us=8000, ...)
        self.spoiler = sc.modules.Spoiler(system, twists=4, voxel_mm=voxel_mm)

    @property
    def duration(self):
        return self.pulse.duration + self.spoiler.duration

    def build(self):
        pulse = self.pulse.build()
        return sc.LogicBlock('fatsat').add(0.0, pulse).add(pulse.duration, self.spoiler.build())
```

`sc.flatten` then reports the whole path — `('gre', 'fatsat', 'spoil')` — so provenance is the tree
and needs no bookkeeping.

---

## Overlap is free

Do not arrange your module so its gradients avoid other modules'. Place them where they belong and
let the compiler deal with it:

```python
seq.add(t, exc.rephaser())          # z
seq.add(t, pe.build(line=k))        # y   } one block, no coordination,
seq.add(t, ro.prephaser_block())    # x   } and no warning
```

Different axes are silent. The same axis is summed with a warning naming both sources. Two RF or two
ADC events overlapping is an error. See [`compiler.md`](compiler.md).

---

## Units are checked for you

Suffix a float attribute with its unit — `_mm`, `_us`, `_deg`, `_per_m`, `_s_per_mm2` — and
`sc.validate.check_units` rejects an implausible value when the constructor returns, with a hint.
`Module` runs it for you when a subclass's `__init__` returns; a component that inherits nothing
calls `sc.validate.check_units(self)` itself, since it takes any object.

```
UnitSanityError: M.slice_thickness_mm = 0.005 is outside the plausible range 0.01 .. 1000 mm.
  got  : 0.005
  hint : 0.005 looks like m. Did you mean slice_thickness_mm=5?
  note : seqcraft never auto-converts - a wrong unit is a wrong sequence.
```

The band comes from the suffix alone, so it is generous by design: `_mm` is checked against
0.01–1000 mm, which catches metres on a slice thickness but not on a field of view. Where a tighter
band is meaningful, state it — `require_in_range(self, 'fov_ro_mm', 0.5, 2000)`.

Public parameters use researcher-natural units (`fov_mm`, `te_ms`, `flip_deg`, `b_value_s_per_mm2`);
anything derived is strict SI with an SI suffix (`fov_m`, `te_s`, `flip_rad`). Convert exactly once,
and never auto-convert.

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

Two tiers, and neither asks what your component inherits.

**`sc.testing.assert_output(make, system)`** takes any callable returning a block, so it fits a
function, one method of several, or a module's `build`. It covers a well-formed block, deterministic
output, gradients on the raster, per-axis limits, and a clean compile on its own:

```python
sc.testing.assert_output(lambda: spoiler(system, twists=4), system)
sc.testing.assert_output(venc.pre, system)
sc.testing.assert_output(venc.post, system)
```

**`sc.testing.assert_all(module, **build_args)`** adds the checks that only mean something for the
`build()` convention — that the call mutates neither the module nor the events on it, that `duration`
is not optimistic, and that `isodelay` / `time_to_echo` fall inside the block — then runs
`assert_output` on the result. It reads `build`, `system`, `regime` and `duration` off the object and
never checks its type, so a class shaped that way passes whether or not it inherits `Module`.

```python
sc.testing.assert_all(venc)                 # the same contract checks the built-in modules get
```

Add the **known values** yourself, because those are the ones that catch physics errors:

```python
def test_venc_encodes_the_right_velocity(system):
    venc = VelocityEncode(system, venc_cm_s=50)
    raster = system.grad_raster.dt
    m1 = 0.0
    for node in venc.build():
        tt, wf = sc.events.waveform_of(node.item, raster)
        m1 += float(sc.events.trapz(wf * (tt + node.start), tt))
    assert m1 == pytest.approx(venc.m1_s_per_m, rel=1e-3)
```

**Measure from the waveform, not from the parameter.** The built-in diffusion module's b-value is
checked that way against numerical integration to 0.5 %, and it is what caught a published ramp
correction being used with the wrong `delta` convention — a 2–4 % b-value error that would have
biased every diffusivity by the same amount, invisibly.

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
handles the refocusing conjugation itself. Mine got the sign convention wrong twice before I gave up
and asked pypulseq.
