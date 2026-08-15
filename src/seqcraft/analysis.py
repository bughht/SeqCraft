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
from .design.logic import flatten
from .design.sampling import sample

if TYPE_CHECKING:
    from types import SimpleNamespace

    from pypulseq.opts import Opts

    from .design.logic import LogicBlock

__all__ = ['kspace', 'moments', 'pns', 'sample']


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
