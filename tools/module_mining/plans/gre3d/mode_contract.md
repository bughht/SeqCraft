# GRE3DTR — Mode Contract, Claim Scope and Dependency Impact

> **Mirrored from the workspace planning repository** (`docs/plans/module-mining/`) on
> 2026-09-19, so that the pull request carries the evidence behind the modules it adds.
> The two copies are identical today and there is nothing keeping them that way; see
> [`../../README.md`](../../README.md) for which one to edit.

**Candidate:** `GRE3DTR`
**Date:** 2026-09-19
**Why this exists:** written *after* the fact, because the mistake it would have prevented was
made. Every candidate with physical modes gets this table **before** implementation from here on;
see [`../fine_scan_playbook.md`](../fine_scan_playbook.md) §19.

---

## 1. The mistake this document is the answer to

The first implementation reasoned:

```text
slab-selective  = shaped RF + Gz
non-selective   = the same thing with Gz removed
```

So the non-selective mode emitted the slab path's **3 ms shaped sinc with no selection gradient**:
a soft pulse that selected nothing and spent three milliseconds of every echo time doing it. It
satisfied every k-space invariant, compiled legally, passed its notebook, and was not the
experiment its own mode name promised. Every official Pulseq 3D reference uses a short hard block
pulse.

Both modes are now derived from their own contract.

## 2. The mode contract

| | non-selective | slab-selective |
|---|---|---|
| selected by | `slab_thickness_mm=None` | `slab_thickness_mm=<mm>` |
| RF family | **hard `block`** | **shaped `sinc`** |
| RF duration, default source | **0.2 ms, this module** — `Excitation`'s 3 ms default is a *shaped* pulse's | 3 ms, `Excitation`'s own default |
| selection gradient | **absent** | present, on z, during the pulse |
| intrinsic rephasing requirement | none (`rephaser_area_per_m == 0`) | `-gz.area/2` from the RF effective centre, reported by `Excitation` and **realised by this kernel** |
| z encoding role | partition encode only | slab rephase **+** partition encode, one combined winder |
| ADC / readout semantics | identical | identical |
| timing consequence | min TE **2.977 ms** at the reference geometry | min TE **4.377 ms** — the pulse is the difference |
| reference evidence | `writeGradientEcho3D.m` (0.2 ms block), `writeMPRAGE*.m`, `write_3Dt1_mprage.py` — all block pulses | `fmrifrey/lps` `gre3d_write_seq.m` — 3 ms sinc, TBW 4, `slabfrac` of the FOV |
| deliberately overridable | `rf_pulse='sinc'` for a shaped but spatially non-selective pulse; `rf_duration_s`; flip | `rf_pulse`, `rf_duration_s`, `rf_time_bw_product` |

### Differential contract

```text
non-selective -> slab-selective

changes:
    RF family            block   -> sinc
    RF duration          0.2 ms  -> 3 ms
    selection gradient   absent  -> present on z
    z winder content     A_partition(p) -> A_slab + A_partition(p)
    minimum TE           2.977   -> 4.377 ms

must stay invariant:
    the kz reached at the echo, for a given partition
    dkz = 1/fov_z, and the centre-partition convention
    kz held through the readout
    ky and kx encoding, ADC timing, rewinding, spoiling
    one winder duration for every partition
```

The invariants are what the tests assert across both modes; the changes are what the mode contract
is *for*.

## 3. Emitted-sequence inspection

`python tools/module_mining/inspect_emitted.py examples/gre_3d/seq/*.seq`, read without reference
to the Python that wrote it:

| | `gre_3d_nonselective.seq` | `gre_3d_slab.seq` |
|---|---|---|
| first RF | 2 samples, **hard / block**, 0.200 ms | 3000 samples, **shaped / soft**, 2.999 ms |
| gradients with the RF | none → **non-selective** | `z` → **spatially selective** |
| acquisition | 256 ADC events, 1056 blocks, 1.46 s | 256 ADC events, 1056 blocks, 2.20 s |

Both files now tell the same physical story as their mode names. Run against the *first*
implementation, the left column would have read `3000 samples, shaped / soft, 2.999 ms` beside
`none → non-selective`, and the contradiction is visible without opening a single Python file.

## 4. Claim scope of the reference comparison

What `tools/module_mining/candidates/gre3d/run_validation.py` measured against `pulseq/pulseq`'s
shipped `gre3d.seq`:

```text
validated against the official GRE3D reference for:
    dk on all three axes           5.000 / 5.000 / 6.250 1/m, exact agreement
    the encoded index of every probed line and partition
    kz held through the readout
    the semantic centre conventions

NOT used to establish:
    the RF family or duration of either mode
    the slab-selective path at all -- the reference is non-selective
    the combined-z realisation, which the reference does not perform
    equal signal amplitude when TE differs
```

The third row is the one that matters. That comparison was correct, and it was never evidence
about the RF — which is precisely where the mistake was. A reference validates what it covers.

The slab-selective path's own evidence is `fmrifrey/lps` read for physics, the signed-moment
arithmetic, and the three degenerate limits, which depend on no external implementation.

## 5. Dependency impact of this phase's shared-leaf changes

Built from the imports, not from the names.

### `Excitation` — additive

```text
+ rephaser_area_per_m        a new property
+ build(rephase=False)       a new argument, default True
```

| consumer | impact |
|---|---|
| `GRE2DTR` | `NO_BEHAVIOR_CHANGE` — never passes `rephase`, so it gets the rephaser it always got |
| `TSEShot` | `NO_BEHAVIOR_CHANGE` — same |
| `GRE3DTR` | `INTENTIONAL_BEHAVIOR_CHANGE` — passes `rephase=False` when selective, and realises that moment itself |
| every example notebook | `NO_BEHAVIOR_CHANGE`, pinned by a content-hash test on the default build |

### `CartesianLine` — a tightened contract

```text
ADC sample counts must be a multiple of opts.adc_samples_divisor, refused at construction
```

| consumer | impact |
|---|---|
| `GRE2DTR` → `GRE2D` | `NEW_EARLY_REFUSAL_FOR_PREVIOUSLY_INVALID_INPUT` — a geometry that used to fail in the compiler now fails at construction, naming both arguments |
| `TSEShot` → `FSE2D` | same. This is the rule the HASTE example was already obeying in prose: 0.6 of 128 is 77 samples, which is why it runs at 0.625 |
| `RadialReadout` | same, inherited by composition rather than by a second copy of the rule |
| **`EPI2D`** | **`NO_BEHAVIOR_CHANGE` — it does not compose `CartesianLine`.** It owns its own train and ADC geometry and only *mentions* the class in a docstring comparison. A map built from names would have got this wrong |
| every example notebook | `NO_BEHAVIOR_CHANGE` — all pass; no shipped protocol was relying on an illegal count |

Nothing is classified `UNKNOWN`. The full suite passing is necessary and is not this map: the
suite says nothing broke, the map says who could have.
