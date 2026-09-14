"""
Draw one GRE repetition as a sequence diagram, coloured by the block that wrote it.

A pulse sequence diagram is normally coloured by **event type** -- one colour for RF, one per
gradient axis, one for the ADC.  That picture is true of every Pulseq sequence ever written and
says nothing about this one.  This figure colours by **LogicBlock** instead, and draws a box around
each one, so what a reader sees is the thing seqcraft actually changes: the same waveforms, grouped
by the module that owns them.

Three facts then need no sentence to explain them.  ``Excitation``'s box spans RF and Gz, because a
slice-selective pulse is one idea on two lanes.  Three boxes overlap across one short window on
three axes -- the slice rephaser, the phase-encode blip and the readout prephaser, which is what
"anything may overlap anything" buys and exactly the window the compiler had to make one block of.
And the outer box is ``GRE2DTR``, which holds the others and emits the ``LIN`` label itself: a
module made of modules, which is the whole composition argument.  Along the bottom are the pulseq
blocks the compiler derived, with their boundaries drawn through everything above.

Everything here is measured.  The waveforms are the events the modules built, the boxes are the
tree's own nodes, and the blocks are what ``sc.compile`` returned; nothing is sketched, so the
figure cannot drift from the library without this script drawing a different picture.
``tests/test_readme_figure.py`` regenerates it and compares, which is what keeps that true.

Gradient lanes are scaled **per axis**, not against one shared maximum: a spoiler is fifteen times
a phase-encode blip, and one scale for all three draws the blip as a flat line.  This figure is
about which module owns which waveform, and the amplitudes are read off the module, not off it.

Run it after changing anything it draws::

    python tools/draw_sequence_figure.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pypulseq as pp

import seqcraft as sc
from seqcraft.design.events import knots_of

ROOT = Path(__file__).resolve().parents[1]
FIGURE = ROOT / 'docs' / 'figures' / 'gre-tr.svg'

#: The protocol the figure is drawn from.  Real enough to compile, chosen so the drawing is legible
#: with every axis still linear -- nothing here is compressed or faked to fit.
#:
#: ``bandwidth_hz_px`` is wide because a wider receive bandwidth is a shorter, taller readout: at
#: the 200 Hz/px default the readout is 6 % of the spoiler beside it and draws as a flat line.
#: ``rf_duration_s`` is half the default because a 3 ms pulse is half the repetition and crowds
#: everything else off the page.  It does not go lower: the slice-select gradient scales as
#: ``1 / duration``, and on a 2 mm slice a 1 ms pulse asks for 117 % of ``max_grad`` and is refused
#: -- which is the limit doing its job, not a drawing problem.  And ``thickness_mm`` matches the
#: in-plane voxel, 256 mm over 128, so both spoilers wind their cycles across the same distance and
#: come out the same length.  That is what per-axis cycles-per-voxel means, and it is easier to see
#: than to read.
PROTOCOL = dict(fov_mm=256.0, matrix=(128, 64), thickness_mm=2.0, bandwidth_hz_px=800.0,
                rf_duration_s=1.5e-3)
LINE = 40

#: The tag of the repetition itself.  Its own events are the kernel's, so they get no inner box.
KERNEL = 'GRE2DTR'

#: One colour per LogicBlock, never per event type.  That substitution is the whole figure.
COLOURS = {
    'Excitation': '#D2694A',
    'PhaseEncode': '#4F7DA8',
    'CartesianLine': '#5B8C5A',
    'spoiler': '#8A7FA8',
    KERNEL: '#B98A2E',
}

#: Lanes, top to bottom.  RF beside Gz so a slice-selective excitation is one contiguous box.
#:
#: Only the ones this repetition actually uses are drawn -- a spoiled GRE fires no trigger, and an
#: empty lane is a row of nothing asking to be read.  A sequence that does fire one gets the lane
#: back without anybody editing this list.
LANES = ('RF', 'Gz', 'Gy', 'Gx', 'ADC', 'Trigger', 'Label')
LANE_HEIGHT = {'RF': 60, 'Gz': 54, 'Gy': 54, 'Gx': 54, 'ADC': 32, 'Trigger': 24, 'Label': 26}

#: Wide enough that a box's label sits in the gap *above* it rather than over its own waveform.
LANE_GAP = 24

#: Narrow on purpose.  Time and amplitude share no units, so the aspect ratio is a drawing choice
#: rather than a fact -- and a wide one flattens every pulse into a plateau with a slope at each
#: end.  Squeezing the time axis is what makes an RF pulse look like a pulse.
GUTTER = 104
PLOT_X0, PLOT_X1 = GUTTER + 18, 752
TOP = 46
WIDTH = 790

INK = '#2B2926'
MUTED = '#7A736A'
FAINT = '#C3BBAD'
CARD = '#FCFBF9'
EDGE = '#E3DED4'
MONO = 'ui-monospace, SFMono-Regular, Menlo, monospace'
SANS = 'ui-sans-serif, -apple-system, Segoe UI, Helvetica, Arial, sans-serif'

#: Rough advance of the label font, for the halo behind a label and for deciding if one fits.
CHAR = 6.45


def lane_geometry(order):
    """Return ``{lane: (top, height)}`` and the y where the block strip begins."""
    tops, y = {}, TOP
    for lane in order:
        tops[lane] = (y, LANE_HEIGHT[lane])
        y += LANE_HEIGHT[lane] + LANE_GAP
    return tops, y + 14


def lanes_used(groups):
    """Return the lanes this repetition actually writes to, in canonical order."""
    used = {lane_of(event) for _, events in groups for _, event in events}
    return [lane for lane in LANES if lane in used]


def lane_of(event):
    """Return the lane an event is drawn on, or ``None`` if it is not drawn."""
    kind = getattr(event, 'type', None)
    if kind == 'rf':
        return 'RF'
    if kind in ('trap', 'grad'):
        return 'G' + str(event.channel)
    if kind == 'adc':
        return 'ADC'
    if kind in ('trigger', 'output'):
        return 'Trigger'
    if kind in ('labelset', 'labelinc'):
        return 'Label'
    return None


def shape_of(event, t0):
    """
    Return ``(times, values, span)`` for one event, in seconds and its own units.

    ``values`` is ``None`` for events with no waveform -- an ADC window or a label -- whose span is
    all the figure draws.
    """
    kind = getattr(event, 'type', None)
    delay = float(getattr(event, 'delay', 0.0) or 0.0)
    if kind in ('trap', 'grad'):
        t, g = knots_of(event, t0)
        return t, g, (float(t[0]), float(t[-1]))
    if kind == 'rf':
        t = np.asarray(event.t, dtype=float) + t0 + delay
        g = np.real(np.asarray(event.signal))
        return t, g, (float(t[0]), float(t[-1]))
    if kind == 'adc':
        start = t0 + delay
        end = start + float(event.num_samples) * float(event.dwell)
        return None, None, (start, end)
    return None, None, (t0 + delay, t0 + delay)


def instances(block):
    """
    Return one entry per child node: ``(tag, [(t, event), ...])``.

    The tree's own nodes are the grouping -- a module put its events in one block, so no clustering
    heuristic is needed to work out which events belong together.
    """
    out = []
    for node in block.nodes:
        item = node.item
        if isinstance(item, sc.LogicBlock):
            tag = item.tag or 'block'
            events = [(t, ev) for t, ev, _ in sc.flatten(item, node.start)]
        else:
            tag = KERNEL
            events = [(node.start, item)]
        drawn = [(t, ev) for t, ev in events if lane_of(ev) is not None]
        if drawn:
            out.append((tag, drawn))
    return out


def runs(indices):
    """Split sorted lane indices into contiguous runs, so a box never spans a lane it does not own."""
    out, run = [], [indices[0]]
    for index in indices[1:]:
        if index == run[-1] + 1:
            run.append(index)
        else:
            out.append(run)
            run = [index]
    out.append(run)
    return out


def build():
    """Return the repetition and the compiled block edges, in seconds."""
    opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    tr = sc.modules.GRE2DTR(opts=opts, **PROTOCOL)
    block = tr(line=LINE)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        seq = sc.compile(sc.LogicBlock('one_tr').add(0.0, block), opts, name='one_tr')
    edges, t = [0.0], 0.0
    for index in range(1, len(seq.block_events) + 1):
        t += seq.block_durations[index]
        edges.append(t)
    return block, edges


def peaks(groups):
    """Return ``{lane: peak}``, per axis rather than shared.  See this module's docstring."""
    out = {}
    for _, events in groups:
        for t0, event in events:
            lane = lane_of(event)
            _, values, _ = shape_of(event, t0)
            if values is not None:
                out[lane] = max(out.get(lane, 0.0), float(np.max(np.abs(values))))
    return out


def label(text, lx, ly, colour):
    """Return a label with a card-coloured halo, so it stays readable over a box or a waveform."""
    width = len(text) * CHAR + 8
    lx = min(lx, WIDTH - 10 - width)
    return (
        f'<rect x="{lx - 4:.2f}" y="{ly - 11:.2f}" width="{width:.2f}" height="15" rx="3" '
        f'fill="{CARD}" fill-opacity="0.9"/>'
        f'<text x="{lx:.2f}" y="{ly:.2f}" font-size="11.5" fill="{colour}" '
        f'font-family="{MONO}">{text}</text>'
    )


def render():
    """Return the whole figure as one SVG document."""
    block, edges = build()
    total = block.duration
    groups = instances(block)
    order = lanes_used(groups)
    tops, strip_y = lane_geometry(order)
    height = strip_y + 84

    def x(t):
        return PLOT_X0 + (t / total) * (PLOT_X1 - PLOT_X0)

    scale_of = peaks(groups)
    last_top, last_height = tops[order[-1]]

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {WIDTH} {height}" '
        f'width="{WIDTH}" height="{height}" font-family="{SANS}" role="img" '
        f'aria-labelledby="t d">',
        '<title id="t">One GRE repetition, coloured and boxed by the LogicBlock that wrote it</title>',
        '<desc id="d">Lanes for RF, Gz, Gy, Gx, ADC and Label carry the waveforms of one '
        'spoiled gradient-echo repetition. Each module that produced them is drawn as a labelled '
        'box in its own colour: Excitation spans the RF and Gz lanes, PhaseEncode the Gy lane, '
        'CartesianLine the Gx and ADC lanes, and two spoilers the Gx and Gz lanes. Three of those '
        'boxes overlap across one short window on three axes. An outer box marks the whole '
        'repetition, GRE2DTR, '
        'which emits the LIN label itself. Along the bottom are the four pulseq blocks the compiler '
        'derived, their boundaries drawn as dashed lines through everything above.</desc>',
        f'<rect x="0.5" y="0.5" width="{WIDTH - 1}" height="{height - 1}" rx="10" '
        f'fill="{CARD}" stroke="{EDGE}"/>',
    ]

    kernel_colour = COLOURS[KERNEL]
    out.append(f'<rect x="{x(0.0) - 16:.2f}" y="{TOP - 20:.2f}" '
               f'width="{x(total) - x(0.0) + 32:.2f}" '
               f'height="{last_top + last_height + 12 - (TOP - 20):.2f}" rx="7" fill="none" '
               f'stroke="{kernel_colour}" stroke-opacity="0.4" stroke-width="1" '
               f'stroke-dasharray="5 4"/>')
    out.append(f'<text x="{x(0.0) - 14:.2f}" y="{TOP - 26:.2f}" font-size="11.5" '
               f'fill="{kernel_colour}" font-family="{MONO}">{KERNEL}</text>')

    for lane in order:
        top, lane_height = tops[lane]
        mid = top + lane_height / 2
        out.append(f'<text x="{GUTTER - 12}" y="{mid + 5:.2f}" font-size="15" font-weight="700" '
                   f'fill="{MUTED}" text-anchor="end" font-family="{MONO}">{lane}</text>')
        out.append(f'<line x1="{x(0.0):.2f}" y1="{mid:.2f}" x2="{x(total):.2f}" y2="{mid:.2f}" '
                   f'stroke="{EDGE}" stroke-width="1"/>')

    for edge in edges[1:-1]:
        out.append(f'<line x1="{x(edge):.2f}" y1="{TOP - 20:.2f}" x2="{x(edge):.2f}" '
                   f'y2="{strip_y + 28:.2f}" stroke="{FAINT}" stroke-width="1" '
                   f'stroke-dasharray="3 4"/>')

    labels = []
    for tag, events in groups:
        if tag == KERNEL:
            continue
        colour = COLOURS.get(tag, MUTED)
        spans = [shape_of(event, t0)[2] for t0, event in events]
        lanes = sorted({order.index(lane_of(event)) for _, event in events})
        t_start = min(span[0] for span in spans)
        t_end = max(span[1] for span in spans)
        for run in runs(lanes):
            first, last = tops[order[run[0]]], tops[order[run[-1]]]
            y0, y1 = first[0] - 5, last[0] + last[1] + 5
            # Barely any horizontal padding: the readout ends exactly where the spoiler begins, and
            # a box drawn generously around each would overlap its neighbour and blur the one
            # boundary the figure is asking to be read.
            x0 = x(t_start) - 1.5
            x1 = max(x(t_end) + 1.5, x0 + 22)
            out.append(f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{x1 - x0:.2f}" '
                       f'height="{y1 - y0:.2f}" rx="5" fill="{colour}" fill-opacity="0.07" '
                       f'stroke="{colour}" stroke-opacity="0.55" stroke-width="1"/>')
            # Above the box, in the lane gap: a label inside it would sit on its own waveform, and
            # a phase-encode blip is exactly the shape a label is wide enough to hide.
            labels.append(label(tag, x0 + 2, y0 - 6, colour))

    for tag, events in groups:
        colour = COLOURS.get(tag, MUTED)
        for t0, event in events:
            lane = lane_of(event)
            top, lane_height = tops[lane]
            mid = top + lane_height / 2
            times, values, span = shape_of(event, t0)
            if lane in ('RF', 'Gx', 'Gy', 'Gz'):
                scale = (lane_height / 2 - 5) / scale_of[lane]
                points = ' '.join(f'{x(t):.2f},{mid - v * scale:.2f}' for t, v in zip(times, values))
                out.append(f'<polyline points="{points}" fill="none" stroke="{colour}" '
                           f'stroke-width="1.8" stroke-linejoin="round"/>')
            elif lane == 'ADC':
                x0, x1 = x(span[0]), x(span[1])
                out.append(f'<rect x="{x0:.2f}" y="{mid - 9:.2f}" width="{x1 - x0:.2f}" height="18" '
                           f'rx="2" fill="{colour}" fill-opacity="0.16" stroke="{colour}" '
                           f'stroke-width="1.2"/>')
                for k in range(1, 40):
                    tick = x0 + (x1 - x0) * k / 40
                    out.append(f'<line x1="{tick:.2f}" y1="{mid - 9:.2f}" x2="{tick:.2f}" '
                               f'y2="{mid + 9:.2f}" stroke="{colour}" stroke-opacity="0.4" '
                               f'stroke-width="0.6"/>')
                labels.append(label(f'{event.num_samples}', (x0 + x1) / 2 - 10, mid + 4, INK))
            elif lane == 'Label':
                px = x(span[0])
                out.append(f'<path d="M{px:.2f} {mid - 6:.2f} l6 6 l-6 6 l-6 -6 z" fill="{colour}"/>')
                labels.append(label(f'{event.label} {int(event.value)}', px + 11, mid + 4, MUTED))

    out.extend(labels)

    # "Pulseq / Blocks", not "blocks": a LogicBlock and a pulseq block are the two things this
    # figure exists to tell apart, and one of them may overlap its neighbours while the other is
    # the flat list that cannot.
    out.append(f'<text x="{GUTTER - 12}" y="{strip_y + 9:.2f}" font-size="15" font-weight="700" '
               f'fill="{MUTED}" text-anchor="end" font-family="{MONO}">Pulseq</text>')
    out.append(f'<text x="{GUTTER - 12}" y="{strip_y + 25:.2f}" font-size="15" font-weight="700" '
               f'fill="{MUTED}" text-anchor="end" font-family="{MONO}">Blocks</text>')
    for index in range(len(edges) - 1):
        x0, x1 = x(edges[index]), x(edges[index + 1])
        out.append(f'<rect x="{x0:.2f}" y="{strip_y:.2f}" width="{x1 - x0:.2f}" height="26" rx="4" '
                   f'fill="#EFEBE3" stroke="#D8D1C4"/>')
        out.append(f'<text x="{(x0 + x1) / 2:.2f}" y="{strip_y + 17:.2f}" font-size="11.5" '
                   f'fill="#4A443D" text-anchor="middle" font-family="{MONO}">{index + 1}</text>')

    axis_y = strip_y + 48
    out.append(f'<line x1="{x(0.0):.2f}" y1="{axis_y:.2f}" x2="{x(total):.2f}" y2="{axis_y:.2f}" '
               f'stroke="{FAINT}" stroke-width="1"/>')
    tick = 0.0
    while tick <= total + 1e-9:
        px = x(tick)
        out.append(f'<line x1="{px:.2f}" y1="{axis_y:.2f}" x2="{px:.2f}" y2="{axis_y + 5:.2f}" '
                   f'stroke="{FAINT}" stroke-width="1"/>')
        out.append(f'<text x="{px:.2f}" y="{axis_y + 18:.2f}" font-size="10.5" fill="{MUTED}" '
                   f'text-anchor="middle" font-family="{MONO}">{tick * 1e3:.0f}</text>')
        tick += 2e-3
    out.append(f'<text x="{x(total) + 14:.2f}" y="{axis_y + 18:.2f}" font-size="10.5" '
               f'fill="{MUTED}" text-anchor="middle" font-family="{MONO}">ms</text>')

    out.append('</svg>')
    return '\n'.join(out) + '\n'


def main():
    """Write the figure, and say whether it moved."""
    svg = render()
    previous = FIGURE.read_text(encoding='utf-8') if FIGURE.exists() else None
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    FIGURE.write_text(svg, encoding='utf-8')
    state = 'unchanged' if svg == previous else 'written'
    sys.stdout.write(f'{FIGURE.relative_to(ROOT).as_posix()}: {state}, {len(svg)} bytes\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
