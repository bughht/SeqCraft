"""
Phase A evidence spike: where does PNS belong in SeqCraft's physical-design loop?

**This is experimental evidence code, not production and not a supported tool.**  It exists so
that the tables in ``plans/current/2026-09-28_pns_aware_physical_design_spike.md`` can be
re-derived rather than trusted, and it reaches into private internals
(``physical_design._design``'s translation, ``design.scope.design_scope``'s ``admissible`` seam)
that no public caller should touch.  If any of those move, fix or delete this file; nothing
depends on it.

Public **synthetic** hardware only -- ``sc.hardware.synthetic_hardware()``, whose coefficients are
pypulseq's own illustrative ``safe_example_hw()`` values.  Every number here is development and
architecture evidence.  **None of it is a scanner-specific claim and none of it clears anything
for human scanning.**

Run with ``python tools/module_mining/spikes/pns_design_spike.py``.
"""

from __future__ import annotations

import warnings

import numpy as np
import pypulseq as pp

import seqcraft as sc
from seqcraft import physical_design as _pd
from seqcraft.design import _augment
from seqcraft.design import scope as _scope

OPTS = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s', B0=3.0,
               rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
               adc_samples_limit=8192)
HW = sc.hardware.synthetic_hardware()

FOV, N, THK, BW = 220.0, 64, 5.0, 500.0
FULL = sc.FlowCompensation(axis=('x', 'y', 'z'))


def gre(opts=OPTS, **kwargs):
    return sc.modules.GRE2DTR(opts=opts, fov_mm=FOV, matrix=(N, N), thickness_mm=THK,
                              bandwidth_hz_px=BW, **kwargs)


def stack(kernel, opts=OPTS):
    """One repetition per line, RF-spoiled, which is what the global evaluator sees."""
    out = sc.LogicBlock('scan')
    for index in range(N):
        out.add(index * kernel.tr_s,
                kernel(line=index, phase_deg=0.5 * 117.0 * index * (index + 1)))
    return out


def peak(tree, opts=OPTS):
    """Peak PNS as a fraction of the model's stimulation limit, and where it happened."""
    got = sc.pns(tree if isinstance(tree, sc.LogicBlock) else sc.LogicBlock('t').add(0.0, tree),
                 opts, HW)
    norm, times = np.asarray(got['norm']), np.asarray(got['t'])
    components = np.asarray(got['components'])
    if not norm.size:
        return 0.0, 0.0, np.zeros(3)
    at = int(np.argmax(norm))
    return float(norm[at]), float(times[at]), components[at]


def one(block):
    return sc.LogicBlock('one').add(0.0, block)


# ------------------------------------------------------------------ 1. the representative set
def families():
    """One repetition per waveform family, with the scan it assembles into."""
    out = {}
    kernel = gre()
    out['GRE 2D'] = (kernel(line=N // 4), stack(kernel))
    for polarity in ('monopolar', 'bipolar'):
        for tag, extra in (('plain', {}), ('all-echo FC', {'flow_comp': FULL})):
            kernel = gre(echoes=4, polarity=polarity, **extra)
            out[f'MEGRE {polarity} {tag}'] = (kernel(line=N // 4), stack(kernel))

    exc = sc.modules.Excitation(opts=OPTS, flip_deg=20.0, thickness_mm=THK, duration_s=1e-3)
    epi = sc.modules.EPI2D(opts=OPTS, fov_mm=FOV, matrix=(N, N), dwell_s=4e-6)
    spoil = sc.modules.spoiler(OPTS, cycles_per_voxel=4.0, voxel_mm=THK, axis='z')
    shot = sc.LogicBlock('gre_epi').add(0.0, exc()).add(float(exc().duration), epi(lines=range(N)))
    shot.add(shot.duration, spoil)
    out['EPI 2D single shot'] = (shot, sc.LogicBlock('scan').add(0.0, shot))

    prep = sc.modules.DiffusionSEPrep(opts=OPTS, thickness_mm=THK, b_s_per_mm2=1000.0, axis='x')
    dwi = sc.LogicBlock('dwi').add(0.0, prep()).add(float(prep().duration), epi(lines=range(N)))
    out['Diffusion SE-EPI b=1000'] = (dwi, sc.LogicBlock('scan').add(0.0, dwi))
    return out


def report_families():
    print('1. local (one repetition) against global (the assembled scan)\n')
    print(f'{"family":<27}{"local":>8}{"global":>8}{"g/l":>7}{"peak at":>10}{"dominant axis":>15}')
    for name, (repetition, scan) in families().items():
        local, at, components = peak(one(repetition))
        overall, _, _ = peak(scan)
        axis = 'xyz'[int(np.argmax(components))]
        print(f'{name:<27}{local:8.3f}{overall:8.3f}{overall / local:7.3f}'
              f'{at * 1e3:9.2f}ms{axis:>15}')
    print('\n   peaks are against the illustrative public example model, not a scanner.')
    print('   the dominant axis is x almost everywhere: the spoiler, at the end of TR.')


# ------------------------------------------------------------------ 2. which axis costs it
def report_axes():
    print('\n\n2. MEGRE monopolar, four echoes: which axis\'s compensation costs PNS\n')
    print(f'{"flow_comp axes":<16}{"local":>8}{"global":>8}{"peak at":>10}'
          f'{"x":>7}{"y":>7}{"z":>7}{"ESP/us":>9}')
    for axes in (None, 'x', 'y', 'z', ('x', 'y'), ('x', 'y', 'z')):
        extra = {} if axes is None else {'flow_comp': sc.FlowCompensation(axis=axes)}
        kernel = gre(echoes=4, polarity='monopolar', **extra)
        local, at, components = peak(one(kernel(line=N // 4)))
        overall, _, _ = peak(stack(kernel))
        print(f'{str(axes):<16}{local:8.3f}{overall:8.3f}{at * 1e3:9.2f}ms'
              f'{components[0]:7.3f}{components[1]:7.3f}{components[2]:7.3f}'
              f'{kernel.ro.echo_spacing_s * 1e6:9.1f}')
    print('\n   compensating x lowers the peak.  Compensating y leaves the LOCAL peak alone and')
    print('   raises the GLOBAL one, so its cost is accumulation across repetitions, not one TR.')


# ------------------------------------------------------------------ 3. what PNS costs in timing
def mark(value: float) -> str:
    """A peak, starred when it is inside the model's limit."""
    return f'{value:.3f}' + ('*' if value <= 1.0 else ' ')


def report_slew():
    print('\n\n3. holding max_grad at 40 mT/m, where does each become admissible?\n')
    print(f'{"max_slew":>9}{"plain local":>13}{"plain glob":>12}'
          f'{"FC local":>11}{"FC glob":>10}{"FC ESP/us":>11}{"FC TR/ms":>10}')
    for slew in (150.0, 100.0, 70.0, 50.0, 35.0, 25.0):
        opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=slew, slew_unit='T/m/s', B0=3.0,
                       rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
        plain = gre(opts, echoes=4, polarity='monopolar')
        done = gre(opts, echoes=4, polarity='monopolar', flow_comp=FULL)
        print(f'{slew:9.0f}{mark(peak(one(plain(line=N // 4)), opts)[0]):>13}'
              f'{mark(peak(stack(plain, opts), opts)[0]):>12}'
              f'{mark(peak(one(done(line=N // 4)), opts)[0]):>11}'
              f'{mark(peak(stack(done, opts), opts)[0]):>10}'
              f'{done.ro.echo_spacing_s * 1e6:11.1f}{done.tr_s * 1e3:10.2f}')
    print('\n   * passes.  At 150 T/m/s the compensated repetition passes and its scan does not,')
    print('   which is the local-pass / global-fail case -- here at the default slew rate.')


# ------------------------------------------------------------------ 4. the admissible seam
def spiral_scope(opts):
    exc = sc.modules.Excitation(opts=opts, flip_deg=15.0, thickness_mm=THK, duration_s=1e-3)
    arm = sc.modules.SpiralReadout(opts=opts, fov_mm=240.0, matrix=64, shots=8,
                                   dwell_s=4e-6, variant='in')
    angles = tuple(2.0 * np.pi * i / 8 for i in range(8))
    return arm, sc.PhysicalDesignScope(
        origin_s=exc.time_to_center(), before=exc(rephase=False),
        after=lambda angle: arm(angle_rad=angle, prephase=False),
        echo_in_after_s=arm.time_to_echo(0) - arm.prephaser_duration_s,
        axes=('x', 'y', 'z'), states=angles,
        design_states=(angles[0], angles[2]), min_window_s=arm.prephaser_duration_s,
    )


def substrate(scope, intents, opts):
    """What `physical_design._design` builds before it calls `design_scope`."""
    claims, encoding = _augment.claims_and_states(*intents)
    keys = tuple(_pd._Key(s, e) for s in scope.states for e in encoding)
    geometry = _pd._geometry(scope, opts)
    requirements = [_pd._requirement(scope, a, keys, claims, encoding, geometry)
                    for a in scope.axes]
    return geometry, requirements, keys


def candidate_block(scope, designed, key):
    """`before + designed + after` -- exactly what `PhysicalDesign.build` emits."""
    out = sc.LogicBlock('candidate')
    out.add(0.0, _pd._piece(scope.before, key.state))
    for axis in scope.axes:
        for at, event in designed.events_for(axis, key):
            out.add(at, event)
    out.add(designed.schedule.window_start_s + designed.window_s,
            _pd._piece(scope.after, key.state))
    return out


def region_only(scope, designed, key):
    """The designed region alone, so its own contribution can be separated from its neighbours."""
    out = sc.LogicBlock('region')
    for axis in scope.axes:
        for at, event in designed.events_for(axis, key):
            out.add(at - designed.schedule.window_start_s, event)
    return out


def report_seam():
    opts = pp.Opts(max_grad=38, grad_unit='mT/m', max_slew=140, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
                   adc_samples_limit=8192)
    arm, scope = spiral_scope(opts)
    print('\n\n4. the `admissible` seam: does widening the owned window reduce PNS?\n')
    for label, intents in (('no intent (m0 only)', (None, None)),
                           ('flow compensation', (FULL, None))):
        geometry, requirements, keys = substrate(scope, intents, opts)
        print(f'   {label}')
        print(f'{"window/us":>12}{"region alone":>14}{"whole candidate":>17}{"peak in arm":>13}')
        for steps in (145, 200, 300, 500, 900, 1500):
            window = steps * float(opts.grad_raster_time)
            got = _scope.realise_scope_at(geometry, requirements, opts, window_s=window)
            if got is None:
                print(f'{window * 1e6:12.0f}   not realisable')
                continue
            mine, _, _ = peak(one(region_only(scope, got, keys[0])), opts)
            whole, at, _ = peak(one(candidate_block(scope, got, keys[0])), opts)
            starts = got.schedule.window_start_s + got.window_s
            print(f'{window * 1e6:12.0f}{mine:14.3f}{whole:17.3f}{str(at >= starts):>13}')
        print()
    print('   the region\'s own contribution falls monotonically -- the lever converges.')
    print('   the whole candidate does not, because the peak is in the arm, which is `after`.')


# ------------------------------------------------------------------ 5. what a leaf can see
def report_leaf():
    print('\n\n5. what a leaf-local designer can see of its own repetition\n')
    print(f'{"case":<28}{"readout alone":>15}{"repetition":>12}{"leaf / rep":>12}')
    for tag, extra in (('MEGRE mono plain', dict(echoes=4, polarity='monopolar')),
                       ('MEGRE mono all-echo FC',
                        dict(echoes=4, polarity='monopolar', flow_comp=FULL)),
                       ('MEGRE bipolar all-echo FC',
                        dict(echoes=4, polarity='bipolar', flow_comp=FULL))):
        kernel = gre(**extra)
        leaf, _, _ = peak(one(kernel.ro()))
        whole, _, _ = peak(one(kernel(line=N // 4)))
        print(f'{tag:<28}{leaf:15.3f}{whole:12.3f}{leaf / whole:12.2f}')
    print('\n   a leaf sees somewhat over half to five sixths of what its repetition plays.')


def report_fill():
    """Widening the owned region against spending the same time as TE fill."""
    print('\n\n6. extra elapsed time as TE fill, rather than as a wider owned region\n')
    print(f'{"flow_comp":>16}{"TE/ms":>9}{"local":>9}{"global":>9}{"peak at":>11}{"dominant":>10}')
    for axes in ('y', ('x', 'y', 'z')):
        for te in (None, 5.5e-3, 9e-3, 16e-3):
            extra = {} if te is None else {'te_s': te}
            kernel = gre(echoes=4, polarity='monopolar',
                         flow_comp=sc.FlowCompensation(axis=axes), **extra)
            local, _, _ = peak(one(kernel(line=N // 4)))
            overall, at, components = peak(stack(kernel))
            print(f'{str(axes):>16}{kernel.te_s * 1e3:9.3f}{local:9.3f}{overall:9.3f}'
                  f'{at * 1e3:10.2f}ms{"xyz"[int(np.argmax(components))]:>10}')
        print()
    print('   fill moves the peak in BOTH directions depending on what is compensated -- worse on')
    print('   y alone, better on x+y+z.  Widening the owned region (section 4) only ever reduced')
    print('   what the owner contributes, so the two interventions are not interchangeable.')


def main():
    warnings.simplefilter('ignore', sc.SeqCraftWarning)
    print(__doc__.strip().splitlines()[0])
    print(f'hardware: {HW!r}\n')
    report_families()
    report_axes()
    report_slew()
    report_seam()
    report_leaf()
    report_fill()


if __name__ == '__main__':
    main()
