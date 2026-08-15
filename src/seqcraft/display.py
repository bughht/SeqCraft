"""
Plotting helpers.

**The only module in seqcraft allowed to import matplotlib**, and it is imported lazily inside
each function so that ``import seqcraft`` stays cheap and side-effect free.  Every function
returns the figure rather than calling ``plt.show()``, so the caller decides what happens to it --
a notebook displays it, a script saves it, a test discards it.

Examples
--------
>>> import pypulseq as pp
>>> import seqcraft as sc
>>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s')
>>> block = sc.LogicBlock('spoiler').add(0.0, pp.make_trapezoid('z', area=500.0, system=opts))
>>> figure = sc.plot_block(block, opts)                 # doctest: +SKIP
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from .core import events as ev
from .core.errors import MissingExtraError, format_error
from .core.logic import BARRIER, LogicBlock, flatten

if TYPE_CHECKING:
    from collections.abc import Iterable

    from pypulseq.opts import Opts

    from .core.compiler import CompiledSequence

__all__ = ['plot_block', 'plot_kspace', 'plot_sequence', 'plot_trajectory']

_AXES = ('x', 'y', 'z')


def _pyplot() -> Any:
    """Import pyplot, with a message naming the extra to install."""
    try:
        import matplotlib.pyplot as plt  # noqa: PLC0415
    except ImportError as err:  # pragma: no cover - depends on the environment
        msg = format_error(
            'plotting needs matplotlib.',
            {'missing': 'matplotlib'},
            ['pip install "seqcraft[viz]"'],
        )
        raise MissingExtraError(msg) from err
    return plt


def _sample(root: LogicBlock, opts: Opts) -> tuple[np.ndarray, dict[str, np.ndarray], list[tuple]]:
    """
    Sample a tree onto a uniform grid: gradients per axis, plus the RF and ADC spans.

    Sampling the *tree* rather than the compiled sequence is deliberate -- it shows what the module
    meant, before the compiler decided where the block boundaries go, which is what you want when
    the question is "did I place this correctly".
    """
    raster = float(opts.grad_raster_time)
    placed = list(flatten(root))
    if not placed:
        return np.zeros(1), {}, []

    total = root.duration
    n = max(2, int(round(total / raster)) + 1)
    grid = np.arange(n) * raster
    grads: dict[str, np.ndarray] = {}
    marks: list[tuple] = []

    for start, event, path in placed:
        kind = getattr(event, 'type', None)
        delay = float(getattr(event, 'delay', 0.0) or 0.0)
        label = '.'.join(path) or kind or '?'
        if kind in ('trap', 'grad'):
            tt, wf = ev.waveform_of(event, raster)
            values = np.interp(grid - start, tt, wf, left=0.0, right=0.0)
            grads[event.channel] = grads.get(event.channel, np.zeros(n)) + values
        elif kind == 'rf':
            marks.append(('rf', start + delay, start + delay + float(event.shape_dur), label))
        elif kind == 'adc':
            span = float(event.num_samples) * float(event.dwell)
            marks.append(('adc', start + delay, start + delay + span, label))
        elif kind == BARRIER:
            marks.append(('barrier', start, start, label))
    return grid, grads, marks


def plot_block(root: LogicBlock, opts: Opts, *, title: str = '', figsize=(10.0, 4.5)) -> Any:
    """
    Plot one logic block: gradients per axis, with RF and ADC windows shaded.

    Parameters
    ----------
    root
        The block to draw, typically what a module call returned.
    opts
        Supplies the raster and gamma, so gradients can be shown in mT/m.
    title
        Figure title.  Defaults to the block's tag.
    figsize
        Matplotlib figure size.

    Returns
    -------
    matplotlib.figure.Figure

    Notes
    -----
    Draws the **tree**, not the compiled sequence.  When a block looks right here but the compiled
    sequence does not, the difference is the compiler's block boundaries -- and
    :func:`plot_sequence` shows those.
    """
    plt = _pyplot()
    grid, grads, marks = _sample(root, opts)
    gamma = float(opts.gamma)

    figure, axis = plt.subplots(figsize=figsize)
    for name in _AXES:
        if name in grads:
            axis.plot(grid * 1e3, grads[name] / gamma * 1e3, label=f'G{name}', linewidth=1.2)

    span = axis.get_ylim()
    for kind, start, end, label in marks:
        if kind == 'rf':
            axis.axvspan(start * 1e3, end * 1e3, color='tab:red', alpha=0.15, zorder=0)
            axis.text(start * 1e3, span[1], f' {label}', fontsize=7, va='top', color='tab:red')
        elif kind == 'adc':
            axis.axvspan(start * 1e3, end * 1e3, color='tab:green', alpha=0.15, zorder=0)
            axis.text(start * 1e3, span[0], f' {label}', fontsize=7, va='bottom', color='tab:green')
        else:
            axis.axvline(start * 1e3, color='0.4', linestyle=':', linewidth=1.0)

    axis.set_xlabel('time (ms)')
    axis.set_ylabel('gradient (mT/m)')
    axis.set_title(title or root.tag or 'logic block')
    axis.grid(alpha=0.25)
    if grads:
        axis.legend(loc='upper right', fontsize=8)
    figure.tight_layout()
    return figure


def plot_sequence(
    compiled: CompiledSequence,
    *,
    time_range: tuple[float, float] | None = None,
    figsize=(11.0, 6.0),
) -> Any:
    """
    Plot a compiled sequence, with the block boundaries the compiler chose marked.

    Parameters
    ----------
    compiled
        The result of :func:`~seqcraft.core.compiler.compile_sequence`.
    time_range
        ``(start, end)`` in seconds.  Defaults to the first 50 ms, because a whole acquisition is
        minutes long and drawing all of it says nothing.
    figsize
        Matplotlib figure size.

    Returns
    -------
    matplotlib.figure.Figure

    Notes
    -----
    The dotted vertical lines are block boundaries.  Seeing them is the point: they are what the
    compiler decided, and a surprising one is usually the explanation for a surprising merge
    warning.
    """
    plt = _pyplot()
    seq = compiled.seq
    lo, hi = time_range if time_range is not None else (0.0, min(50e-3, compiled.duration_s))
    gamma = float(compiled.opts.gamma)
    raster = float(compiled.opts.grad_raster_time)

    figure, (top, bottom) = plt.subplots(
        2, 1, figsize=figsize, sharex=True, height_ratios=(1, 2)
    )
    t = 0.0
    boundaries: list[float] = []
    traces: dict[str, list[tuple[np.ndarray, np.ndarray]]] = {a: [] for a in _AXES}
    rf_spans: list[tuple[float, float]] = []
    adc_spans: list[tuple[float, float]] = []

    for index in sorted(seq.block_events):
        duration = float(seq.block_durations[index])
        if t + duration >= lo and t <= hi:
            boundaries.append(t)
            block = seq.get_block(index)
            for name in _AXES:
                grad = getattr(block, f'g{name}', None)
                if grad is not None:
                    tt, wf = ev.waveform_of(grad, raster)
                    traces[name].append((tt + t, wf / gamma * 1e3))
            if getattr(block, 'rf', None) is not None:
                rf = block.rf
                rf_spans.append((t + float(rf.delay), t + float(rf.delay) + float(rf.shape_dur)))
            if getattr(block, 'adc', None) is not None:
                adc = block.adc
                start = t + float(adc.delay)
                adc_spans.append((start, start + float(adc.num_samples) * float(adc.dwell)))
        t += duration
        if t > hi:
            break
    boundaries.append(min(t, hi))

    for start, end in rf_spans:
        top.axvspan(start * 1e3, end * 1e3, color='tab:red', alpha=0.5)
    for start, end in adc_spans:
        top.axvspan(start * 1e3, end * 1e3, color='tab:green', alpha=0.5)
    top.set_yticks([])
    top.set_ylabel('RF / ADC')
    top.set_title(
        f'{compiled.definitions.get("Name", "sequence")}  --  {compiled.n_blocks} blocks, '
        f'{compiled.duration_s:.3f} s'
    )

    for name, colour in zip(_AXES, ('tab:blue', 'tab:orange', 'tab:purple')):
        first = True
        for times, values in traces[name]:
            bottom.plot(
                times * 1e3, values, color=colour, linewidth=1.1,
                label=f'G{name}' if first else None,
            )
            first = False
    for edge in boundaries:
        bottom.axvline(edge * 1e3, color='0.75', linestyle=':', linewidth=0.8, zorder=0)

    bottom.set_xlabel('time (ms)')
    bottom.set_ylabel('gradient (mT/m)')
    bottom.set_xlim(lo * 1e3, hi * 1e3)
    bottom.grid(alpha=0.25)
    bottom.legend(loc='upper right', fontsize=8)
    figure.tight_layout()
    return figure


def plot_kspace(compiled: CompiledSequence, *, figsize=(5.5, 5.5), max_points: int = 200_000) -> Any:
    """
    Plot the acquired k-space trajectory in the kx--ky plane.

    Parameters
    ----------
    compiled
        The compiled sequence.
    figsize
        Matplotlib figure size.
    max_points
        Decimate above this many samples, so a full acquisition still draws in a second.

    Returns
    -------
    matplotlib.figure.Figure
    """
    plt = _pyplot()
    k = compiled.kspace()['k_adc']
    step = max(1, k.shape[1] // max_points)
    kx, ky = k[0, ::step], k[1, ::step]

    figure, axis = plt.subplots(figsize=figsize)
    axis.plot(kx, ky, linewidth=0.4, alpha=0.8)
    axis.plot(kx, ky, ',', color='tab:red', alpha=0.4)
    axis.set_xlabel('$k_x$ (1/m)')
    axis.set_ylabel('$k_y$ (1/m)')
    axis.set_title(f'{compiled.definitions.get("Name", "sequence")} -- {k.shape[1]} samples')
    axis.set_aspect('equal')
    axis.grid(alpha=0.25)
    figure.tight_layout()
    return figure


def plot_trajectory(
    interleaves: Iterable[tuple[np.ndarray, np.ndarray]],
    *,
    title: str = '',
    figsize=(5.5, 5.5),
) -> Any:
    """
    Plot a set of k-space interleaves, one colour per shot.

    Parameters
    ----------
    interleaves
        An iterable of ``(kx, ky)`` array pairs in 1/m -- one pair per shot.  A spiral's shots, a
        radial set's spokes, one segment of an EPI train: this takes the numbers, so it works for
        any of them.
    title
        Figure title.  Defaults to the interleaf count.
    figsize
        Matplotlib figure size.

    Returns
    -------
    matplotlib.figure.Figure

    Notes
    -----
    Takes arrays rather than a readout object on purpose.  It previously took anything with
    ``.trajectory()``, ``.n_interleaves`` and ``.k_max_per_m`` -- a duck type only one class in one
    module library ever satisfied, which quietly made a plotting helper depend on that library.

    The picture to look at when deciding an interleaf count and a density: the radial gap between
    adjacent turns of the *combined* set is what has to stay within ``1/FOV``, and undersampling
    shows up here as visible white space at the edge long before it shows up in an image.

    Examples
    --------
    >>> import numpy as np, seqcraft as sc
    >>> theta = np.linspace(0, 8 * np.pi, 512)
    >>> shots = [(theta * np.cos(theta + p), theta * np.sin(theta + p)) for p in (0.0, 3.1)]
    >>> figure = sc.plot_trajectory(shots)              # doctest: +SKIP
    """
    plt = _pyplot()
    figure, axis = plt.subplots(figsize=figsize)
    limit = 0.0
    count = 0
    for kx, ky in interleaves:
        kx, ky = np.asarray(kx), np.asarray(ky)
        axis.plot(kx, ky, linewidth=0.7, alpha=0.85)
        limit = max(limit, float(np.max(np.abs(kx))), float(np.max(np.abs(ky))))
        count += 1
    if limit > 0.0:
        axis.set_xlim(-limit * 1.05, limit * 1.05)
        axis.set_ylim(-limit * 1.05, limit * 1.05)
    axis.set_xlabel('$k_x$ (1/m)')
    axis.set_ylabel('$k_y$ (1/m)')
    axis.set_title(title or f'{count} interleaves')
    axis.set_aspect('equal')
    axis.grid(alpha=0.25)
    figure.tight_layout()
    return figure
