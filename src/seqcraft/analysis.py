"""
Measuring a tree: one entry shape, four answers.

**Give it a tree, get numbers back.**  You should not have to know that PNS prediction needs a
compiled sequence while a moment does not, so none of these ask for one: they take the tree you
built and the scanner you built it against, and compile internally where compiling is what the
question requires.

This module sits *above* the compiler in the layering -- :func:`kspace` and :func:`pns` import it --
and below :mod:`seqcraft.display`, which draws what :func:`sample` returns.

Which of these is exact, and which is not
-----------------------------------------
The single most important thing to get right here.

=============  =====================================================  ========
Function       Basis                                                  Exact?
=============  =====================================================  ========
``sample``     uniform raster grid, **interpolated**                   **no**
``moments``    ``knots_of`` + ``pwl_moment`` over ``flatten(tree)``    **yes**
``kspace``     compiled, then ``calculate_kspacePP()``                 yes, at true ADC times
``pns``        compiled, then ``calculate_pns()``                      pypulseq's validated SAFE
=============  =====================================================  ========

:func:`moments` looks like it could be built on :func:`sample` now that they sit in one file.  It
must not be.  Sampling interpolates onto a uniform grid -- an arbitrary gradient's samples are at
raster *centres* and an extended trapezoid's knots are not uniformly spaced at all -- which is a
2.5 % amplitude error on a spiral, and the compiler's own self-check would then be comparing two
differently-wrong numbers.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from .compiler import compile_sequence
from .design import events as ev
from .design.logic import BARRIER, flatten

if TYPE_CHECKING:
    from types import SimpleNamespace

    from pypulseq.opts import Opts

    from .design.logic import LogicBlock

__all__ = ['kspace', 'moments', 'pns', 'sample']


def sample(
    root: LogicBlock,
    opts: Opts,
) -> tuple[np.ndarray, dict[str, np.ndarray], list[tuple]]:
    """
    Sample a tree onto a uniform grid: gradients per axis, plus the RF and ADC spans.

    Parameters
    ----------
    root
        The block to sample.  Nested blocks are flattened, so a component that nests is measured
        by its actual events rather than by its direct children.
    opts
        Supplies the gradient raster, which is the grid spacing.

    Returns
    -------
    grid
        Times in seconds from the start of `root`, spaced by ``opts.grad_raster_time``.
    grads
        ``channel -> amplitudes`` in Hz/m on `grid`, summed across everything on that axis.  Only
        the axes actually used appear.
    marks
        ``(kind, start, end, label)`` for each RF, ADC and barrier, in tree order.  ``kind`` is
        ``'rf'``, ``'adc'`` or ``'barrier'``; a barrier has ``start == end``.

    Notes
    -----
    **It samples the tree, not the compiled sequence** -- deliberately.  It shows what you *meant*,
    before the compiler chose block boundaries, which is what you want when the question is "did I
    place this correctly".  Use ``seq.plot()`` to see what the compiler made of it; the difference
    between the two pictures is exactly the block structure.

    **Lossy on purpose.**  The grid is uniform, which is convenient and inexact: an arbitrary
    gradient's samples sit at raster *centres* and an extended trapezoid's knots are not uniformly
    spaced at all, so both are **interpolated** onto the grid rather than reproduced on it.  On a
    spiral the difference is a 2.5 % amplitude error.  For a moment, a split or a sum use
    :func:`moments`, or :func:`~seqcraft.design.events.knots_of` with
    :func:`~seqcraft.design.events.pwl_moment` directly, which are exact.

    Examples
    --------
    >>> import numpy as np, pypulseq as pp, seqcraft as sc
    >>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s')
    >>> block = sc.LogicBlock('spoiler').add(0.0, pp.make_trapezoid('z', area=500.0, system=opts))
    >>> grid, grads, marks = sc.sample(block, opts)
    >>> sorted(grads), marks
    (['z'], [])
    >>> bool(np.allclose(np.diff(grid), opts.grad_raster_time))
    True
    >>> round(float(np.trapezoid(grads['z'], grid)), 1)          # the area, near enough to see
    500.0
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


def moments(tree: LogicBlock, order: int = 0) -> dict[str, float]:
    """
    Return the whole-tree gradient moment per axis, integrated from exact knots.

    Takes no ``Opts`` and does no compile: a moment is a property of the waveform, not of the
    scanner that plays it.

    Parameters
    ----------
    tree
        The block to measure.  Nested blocks are flattened, so a component that nests is measured
        by its actual events rather than by its direct children.
    order
        ``0`` for area -- k-space displacement, 1/m.  ``1`` for the first moment, s/m: flow and
        motion sensitivity, and whether a bipolar pair is truly nulled.  ``2`` for s^2/m,
        acceleration sensitivity.  Referenced to the start of `tree`.

    Returns
    -------
    dict
        ``axis -> moment``, over the axes actually used.

    Notes
    -----
    Integrated from :func:`~seqcraft.design.events.knots_of` and
    :func:`~seqcraft.design.events.pwl_moment`, never from :func:`sample`.  A moment is **linear
    in the waveform**, so summing per-event moments is exact whether the events overlap or not --
    no union of knot sets needs building.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s')
    >>> tree = sc.LogicBlock('spoil').add(0.0, pp.make_trapezoid('z', area=500.0, system=opts))
    >>> round(sc.moments(tree)['z'], 6)
    500.0
    """
    out: dict[str, float] = {}
    for start, event, _path in flatten(tree):
        if getattr(event, 'type', None) not in ev.GRADIENT_KINDS:
            continue
        axis = str(event.channel)
        out[axis] = out.get(axis, 0.0) + ev.pwl_moment(*ev.knots_of(event, start), order)
    return out


def kspace(tree: LogicBlock, opts: Opts) -> dict[str, np.ndarray]:
    """
    Return the k-space trajectory of `tree`, in 1/m.

    Compiles internally, then uses pypulseq's own calculation, so the sample times are the **true
    ADC sample times** rather than a raster approximation.

    Returns
    -------
    dict
        ``k_adc`` (3 x n_samples, at the ADC sample times), ``t_adc``, ``k`` (dense), ``t_k``,
        ``t_excitation``, ``t_refocusing``.

    Notes
    -----
    This exists rather than a direct call because ``calculate_kspacePP`` returns its tuple in a
    **different order** from ``calculate_kspace``, and getting that wrong silently swaps the
    trajectory for its timebase -- a wrong answer with no error.  The named return makes that
    unmistakable.
    """
    seq = compile_sequence(tree, opts)
    k_adc, t_adc, k, t_k, t_exc, t_refoc = seq.calculate_kspacePP()[:6]
    return {
        'k_adc': np.asarray(k_adc),
        't_adc': np.asarray(t_adc),
        'k': np.asarray(k),
        't_k': np.asarray(t_k),
        't_excitation': np.asarray(t_exc),
        't_refocusing': np.asarray(t_refoc),
    }


def pns(tree: LogicBlock, opts: Opts, hardware: SimpleNamespace) -> dict[str, Any]:
    """
    Predict peripheral nerve stimulation for `tree` against a gradient hardware model.

    Compiles internally, then delegates to pypulseq's SAFE model implementation.

    Parameters
    ----------
    tree, opts
        What to measure, and the scanner to compile it against.
    hardware
        The gradient *response* model, from :func:`seqcraft.hardware.load_hardware` or
        :func:`seqcraft.hardware.synthetic_hardware`.  **Required**: it describes how the body
        responds to being driven, not how hard the amplifier may be driven, so it has nothing to
        do with `opts` and is not carried on it.

    Returns
    -------
    dict
        ``ok``, ``peak`` (fraction of the stimulation limit), ``norm`` (the time-resolved curve),
        ``components`` (per-axis contributions) and ``t`` (their timebase).

    Notes
    -----
    **The full return matters.**  When ``ok`` is ``False``, ``peak`` says how much but ``norm``
    and ``t`` say *where*, which is what you need to fix it.

    This delegates rather than reimplementing.  ``dG/dt`` convolved with three exponentials is
    about sixty lines and looks reachable from :func:`sample`, but it is a *safety* calculation,
    pypulseq's implementation is validated against vendor behaviour, and a second one that can
    silently drift is the wrong thing to own.

    ``synthetic_hardware()`` is a conservative vendor-free stand-in for CI.  It is **not** a real
    scanner and must never be used to clear a human scan.
    """
    seq = compile_sequence(tree, opts)
    ok, pns_norm, components, t = seq.calculate_pns(hardware, do_plots=False)
    return {
        'ok': bool(ok),
        'peak': float(np.max(pns_norm)) if np.size(pns_norm) else 0.0,
        'norm': np.asarray(pns_norm),
        'components': np.asarray(components),
        't': np.asarray(t),
    }
