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
import pypulseq as pp

from .compiler import compile_sequence
from .design import events as ev
from .design.logic import BARRIER, flatten

if TYPE_CHECKING:
    from types import SimpleNamespace

    from pypulseq.opts import Opts

    from .design.logic import LogicBlock

__all__ = ['b_value', 'kspace', 'moments', 'pns', 'sample']


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
    >>> round(float(ev.trapz(grads['z'], grid)), 1)              # the area, near enough to see
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


def b_value(tree: LogicBlock, opts: Opts, *, end_s: float | None = None) -> dict[str, float]:
    r"""
    Return the diffusion `b`-value per axis, in s/mm\ :sup:`2`, integrated from emitted gradients.

    .. math::

        b = (2\pi)^2 \int_0^{T} k(t)\cdot k(t)\, \mathrm{d}t
        \qquad
        k(t) = \int_0^{t} G(t')\, \mathrm{d}t'

    Bernstein, King & Zhou, *Handbook of MRI Pulse Sequences* (2004), §9.1, Eqs. 9.5--9.7.  In
    SeqCraft gradients are already in Hz/m, so :math:`\int G\,\mathrm{d}t` **is** the handbook's
    :math:`k(t)` in 1/m and no gyromagnetic ratio appears here.

    **This measures the tree, not any module's arithmetic.**  It shares no code with whatever
    designed the waveform, which is the point: a module that both computes its own `b` and is
    checked against its own number has tested nothing.  And `b` is a property of *everything* on
    the axis -- slice-select lobes, crushers and readout prephasers all contribute -- so this
    integrates the whole tree rather than one component.

    Parameters
    ----------
    tree
        The block to measure.  Nested blocks are flattened.
    opts
        Supplies the gradient raster, which is the integration grid.
    end_s
        Integrate to this time from the start of `tree`.  ``None`` integrates to the end.  For a
        spin echo the physically meaningful endpoint is the echo; integrating past it keeps
        accumulating and reports a number no experiment measures.

    Returns
    -------
    dict
        ``axis -> b`` over the axes actually used, plus ``'total'``, the sum -- which is the
        scalar `b` of a single-direction experiment.  A full b-matrix would need the cross terms
        and is deliberately not returned: no candidate has needed one, and the trace is what a
        b-value means.

    Notes
    -----
    **Refocusing pulses flip the sign of the accumulated k.**  That is what makes the spin-echo
    case work at all -- the handbook's §9.1 notes the two lobes "have the same polarity and are
    placed at either side of a refocusing RF pulse", which only integrates to a large `b` because
    the 180 conjugates what came before.  Every RF event with ``use='refocusing'`` flips it here,
    at its own effective centre.

    **The origin is the excitation, not the start of the tree.**  Integration runs from the first
    ``use='excitation'`` pulse's centre, with `k` measured from its value there, because that is
    when transverse magnetisation exists and where its phase starts.  It matters on the slice
    axis, whose rephaser undoes the slice-select area *after* the RF centre rather than half the
    lobe: measuring from ``t = 0`` leaves the other half as a constant offset, and a constant `k`
    integrates without bound -- reporting a `b` that grows with the echo time out of nothing.  A
    tree with no excitation in it is measured from its own start.

    Both of these were wrong in the first version of this function and neither was visible on the
    diffusion axis, where `k` is flat across the refocusing block and zero before the encoding.
    ``examples/dwi_se_epi_2d/02_simulate_and_reconstruct.ipynb`` is what found them: the Bloch
    simulator and this integral disagreed about a `b = 0` acquisition by a factor of nine.

    Examples
    --------
    >>> import numpy as np, pypulseq as pp, seqcraft as sc
    >>> opts = pp.Opts(max_grad=80, grad_unit='mT/m', max_slew=200, slew_unit='T/m/s',
    ...                rf_dead_time=100e-6, rf_ringdown_time=30e-6)
    >>> g = pp.make_trapezoid('x', amplitude=25e-3 * 42.576e6, rise_time=200e-6,
    ...                       flat_time=10e-3, system=opts)
    >>> tree = sc.LogicBlock('mono').add(0.0, g)
    >>> round(sc.b_value(tree, opts)['total'], 1)
    16.3

    One lobe on its own is a weak diffusion weighting; the pair straddling a refocusing pulse is
    what makes a useful `b`.  ``examples/dwi_se_epi_2d/`` shows the difference.
    """
    grid, grads, _spans = sample(tree, opts)
    if grid.size < 2:
        return {}
    # An RF pulse acts at its own effective centre, which is the pulse's delay -- transmit dead
    # time and any asymmetry of the envelope -- plus the centre within the waveform.  Omitting
    # the delay puts a refocusing pulse early, which is invisible on an axis whose k happens to
    # be flat across the block and badly wrong on one that is not.
    centres = [
        (start + float(event.delay) + float(pp.calc_rf_center(event)[0]),
         getattr(event, 'use', None))
        for start, event, _path in flatten(tree)
        if getattr(event, 'type', None) == 'rf'
    ]
    # Transverse magnetisation exists only after the excitation, and its phase origin is that
    # pulse's centre -- which is why a slice-select rephaser undoes the area *after* the centre
    # rather than half the lobe.  Measuring k from the start of the tree instead leaves that
    # half-lobe as a spurious constant offset, and a constant k integrates without bound.
    origin = min((t for t, use in centres if use == 'excitation'), default=float(grid[0]))
    # Where the accumulated k is conjugated, in increasing time order.
    flips = sorted(t for t, use in centres if use == 'refocusing')
    stop = float(grid[-1]) if end_s is None else float(end_s)
    raster = 0.5 * float(opts.grad_raster_time)
    keep = (grid >= origin - raster) & (grid <= stop + raster)

    out: dict[str, float] = {}
    for axis, waveform in grads.items():
        # The running integral on the raster, trapezoid because the sampled waveform is
        # piecewise linear between grid points.
        k = np.concatenate([[0.0], np.cumsum(
            0.5 * (waveform[1:] + waveform[:-1]) * np.diff(grid))])
        k = k - float(np.interp(origin, grid, k))
        # A refocusing pulse at `when` conjugates what has accumulated:
        #     k_after(t) = -k(when) + [k(t) - k(when)] = k(t) - 2 k(when)
        # applied in time order, so a second pulse conjugates the already-conjugated value.
        for when in flips:
            after = grid >= when
            if not after.any():
                continue
            k = np.where(after, k - 2.0 * float(np.interp(when, grid, k)), k)
        out[axis] = float((2.0 * np.pi) ** 2 * np.trapezoid(k[keep] ** 2, grid[keep]) / 1e6)
    if out:
        out['total'] = float(sum(out.values()))
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
