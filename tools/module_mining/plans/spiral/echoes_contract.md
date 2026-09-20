# SpiralReadout — the `echoes` contract

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-20, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

Separate from the variant table on purpose. `echoes` and `variant` are **orthogonal**, and writing
them as one table would produce a four-by-N Cartesian product describing one rule four times.

---

## 1. What one echo is

> **One echo is one traversal of the path that crosses the origin.** The echo instant is the
> crossing.

That definition is what makes `echoes` countable across variants that do not look alike. It also
explains the asymmetry the variant table found: `out`, `in` and `in-out` cross the origin once per
unit, `out-in` crosses it twice — so a single `out-in` unit **is** two echoes by this definition,
and reporting it as one would be the module lying about its own trajectory.

```text
out       1 crossing per unit    at the unit's start
in        1                      at the unit's end
in-out    1                      at the seam
out-in    2                      at both ends
```

## 2. The transition between echoes, and when it is free

After one unit the traversal sits at either the origin or `k_max`. The next unit must **start**
where the previous one **ended**, or something has to move k in between.

```text
one-arm variants (out, in)
    unit ends where the next must not start
    -> a FLYBACK is required: k_max -> 0 for `out`, 0 -> k_max for `in`
    -> it costs time and it acquires nothing

two-arm variants (in-out, out-in)
    unit ends where the next one starts, and both are at rest there
    -> consecutive units join CONTINUOUSLY, no connector, no flyback
    -> this is the trade `variant` carries, and it is the spiral-scale version of
       CartesianLine's monopolar/bipolar choice
```

**This is why `echoes` needs no `polarity` beside it.** The variant already says whether the train
pays for its transitions. A separate `polarity` would be a second name for the same fact.

## 3. Echo spacing

`echo_spacing_s` is **derived, not requested**: it is the time from one origin crossing to the
next, and it is fixed by the arm's traversal plus the transition. A caller cannot lengthen it
without changing `density`, `shots` or the scanner, and a caller cannot shorten it at all.

So the module **reports** it and refuses a requested value it cannot meet, naming which of the
three inputs would have to change. It does not pad to an asked-for spacing: padding between echoes
of one readout would put dead time where the trajectory expects continuity.

## 4. ADC and echo boundaries

The rule is the one the variant table already states, applied per train:

> **One ADC region per contiguous stretch of non-zero gradient**, segmented against the
> interpreter's sample limit, and every origin crossing falls inside a sampling window.

Consequences, which differ by variant and are the reason this section exists:

- **two-arm variants:** consecutive units are one continuous gradient, so **one ADC region can
  span several echoes**. The segmentation is then driven by the sample limit alone and has no
  relationship to the echo boundaries — a segment boundary may fall anywhere, **including on an
  echo**, and the contract requires only that the crossing be inside a window, not that it be
  inside a particular segment;
- **one-arm variants:** the flyback is a gap in acquisition, so each unit is its own region, and
  the segmentation restarts.

The failure this prevents is PR23's **D**: splitting per arm rather than per acquisition region
puts an in-out readout's `k = 0` into the dead-time guard between two ADC events instead of inside
a window.

## 5. What `echoes` does **not** change

`density`, `shots`, `fov_mm`, `matrix`, `dwell_s` and `angle_rad` are all independent of it. The
arm is designed once regardless of how many times it is traversed. That independence is what
keeps this a separate contract rather than a dimension of the variant table.

## 6. Incremental implementation is legitimate here

`echoes=1` is the whole contract with the train length set to one, and every rule above degenerates
correctly: no transition, one region, spacing undefined rather than wrong. So a first
implementation may ship `echoes=1` only, **provided** the reported echo instants are already a
sequence (mode contract §4) and `echo_spacing_s` is already derived rather than stored.

What must not happen is a first implementation in which "the echo time" is a scalar attribute,
because adding `echoes` — or `out-in` — then changes what that attribute *means*.
