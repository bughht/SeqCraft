"""
The one picture pypulseq cannot draw.

:func:`plot_block` shows a **tree**, before the compiler chose block boundaries -- which is what
you want when the question is "did I place this correctly".  For a compiled sequence use
``seq.plot()``: pypulseq's own plotter is better maintained and it is the picture everyone else in
the ecosystem reads.  The difference between the two is exactly the block structure the compiler
chose.

Two plotters used to live here and no longer do.  ``plot_sequence`` drew a compiled sequence,
which ``Sequence.plot()`` already does; ``plot_trajectory`` took bare ``(kx, ky)`` arrays, and
``sc.kspace(...)`` plus three lines of matplotlib draws them at the call site, where it is visible
which two of the three axes are being shown.

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

from .analysis import sample
from .design.events import AXES
from .errors import MissingExtraError, format_error

if TYPE_CHECKING:
    from pypulseq.opts import Opts

    from .design.logic import LogicBlock

__all__ = ['plot_block']


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
    sequence does not, the difference is the compiler's block boundaries -- and ``seq.plot()``
    shows those.
    """
    plt = _pyplot()
    grid, grads, marks = sample(root, opts)
    gamma = float(opts.gamma)

    figure, axis = plt.subplots(figsize=figsize)
    for name in AXES:
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
