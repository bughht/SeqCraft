"""
What a compile returns.

:class:`CompiledSequence`
    The :class:`pypulseq.Sequence`, the compile report, and per-block provenance -- an *object*
    rather than a printed summary, so a failing check can be detected programmatically.  It
    answers the questions you ask after a build: how many blocks, how long, what came from where,
    what is the trajectory, does it pass every check, write it.  The reference implementation's
    ``get_report()`` printed ``check_timing()`` and returned ``None``, so a failing sequence could
    not be detected in code at all; this type is the answer to that.

:mod:`~seqcraft.result.provenance`
    The JSON sidecar written beside the ``.seq``, recording versions, git state, the definitions
    and the file hash, so a file six months old can still say what produced it.  Written by
    default -- :meth:`CompiledSequence.write` takes ``sidecar=False`` to suppress it, which is
    what the byte-comparison tests do, since a timestamp is exactly what they cannot compare.

:class:`~seqcraft.report.Issue` and :class:`~seqcraft.report.Report` are **not** here, though a
compile returns one.  They are the vocabulary the compiler writes findings in -- five of its stage
modules build ``Issue`` objects long before any ``CompiledSequence`` exists -- so they live at the
root beside :mod:`seqcraft.errors`, which is the other half of the same pair: hard failures raise,
soft findings report.  Keeping them here made ``compiler/`` import ``result/`` for a type, which
is backwards.

**Nothing here imports the compiler**, and a test asserts it.  That is not tidiness: it is what
lets a result be constructed, inspected and written without the transformation that made it, and
it is why :func:`~seqcraft.compiler.verification.verify_against_tree` -- which reads the private
IR -- is a compile stage rather than a method on :class:`CompiledSequence`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

from ..report import Issue, Report, ReportFailed
from .provenance import _jsonable, write_sidecar

if TYPE_CHECKING:
    from pathlib import Path
    from types import SimpleNamespace

    from pypulseq.opts import Opts

__all__ = ['CompiledSequence', 'Issue', 'Report', 'ReportFailed', 'WriteResult']


@dataclass(frozen=True)
class WriteResult:
    """What :meth:`CompiledSequence.write` produced."""

    path: Path
    sha256: str
    n_blocks: int
    duration_s: float
    sidecar: Path | None = None


@dataclass
class CompiledSequence:
    """
    A compiled sequence: the pypulseq object, the compile report, and provenance.

    Attributes
    ----------
    seq
        The :class:`pypulseq.Sequence`.  Yours to use directly; seqcraft never hides it.
    opts
        The scanner it was compiled against, kept so the post-compile checks and the provenance
        sidecar cannot disagree with what the boundaries were chosen for.
    origins
        One tag path per compiled block, so a block index traces back to the module that
        produced it.
    definitions
        The ``[DEFINITIONS]`` that will be written, already merged and collision-checked.  Plain
        pulseq keys: there is no geometry object here, because everything the file says about the
        scan came in as this mapping.
    report
        Empty, and last, and defaulted.  Everything it used to carry either raises or is a
        :class:`~seqcraft.errors.SeqCraftWarning` now, so the compiler no longer constructs one --
        which is what lets ``compiler/`` stop importing ``report`` at all.
    """

    seq: Any
    opts: Opts
    origins: tuple[tuple[str, ...], ...]
    definitions: dict[str, Any]
    tree_duration_s: float
    report: Report = field(default_factory=Report)
    _checked: Report | None = field(default=None, repr=False)

    # ------------------------------------------------------------------------ properties
    @property
    def n_blocks(self) -> int:
        """Number of pulseq blocks."""
        return len(self.seq.block_events)

    @property
    def duration_s(self) -> float:
        """Total duration, seconds.  ``Sequence.duration()`` returns a tuple; this does not."""
        return float(self.seq.duration()[0])

    def origin(self, block_index: int) -> tuple[str, ...]:
        """
        Return the tag path of the module that produced block `block_index`.

        Examples
        --------
        >>> import pypulseq as pp
        >>> import seqcraft as sc
        >>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
        ...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
        >>> inner = sc.LogicBlock('spoiler')
        >>> _ = inner.add(0.0, pp.make_trapezoid('z', area=500.0, system=opts))
        >>> out = sc.compile(sc.LogicBlock('tr').add(0.0, inner), opts)
        >>> out.origin(0)
        ('tr', 'spoiler')
        """
        return self.origins[block_index]

    def __repr__(self) -> str:
        """One line: name, block count, duration, error and warning counts."""
        return (
            f'CompiledSequence({self.definitions.get("Name", "?")}, {self.n_blocks} blocks, '
            f'{self.duration_s:.3f} s, {len(self.report.errors)} errors, '
            f'{len(self.report.warnings)} warnings)'
        )

    def check(self) -> Report:
        """
        Return the compile report, with the block count and duration attached.

        Everything this used to *find* now raises during the compile: limit violations, event
        sizes, duplicate k-space addresses and pypulseq's timing audit.  What is left is a view
        of what the compile *did*, so ``ok`` is always true and there is nothing to forget to
        check.  It is on its way out with the rest of this type.

        Examples
        --------
        >>> import pypulseq as pp
        >>> import seqcraft as sc
        >>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
        ...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
        >>> lb = sc.LogicBlock('t')
        >>> _ = lb.add(0.0, pp.make_trapezoid('x', area=100.0, system=opts))
        >>> sc.compile(lb, opts).check().ok
        True
        """
        if self._checked is not None:
            return self._checked
        out = Report(self.report.issues, subject=self.report.subject, values={
            'n_blocks': self.n_blocks,
            'duration_s': self.duration_s,
        })
        object.__setattr__(self, '_checked', out)
        return out

    # ---------------------------------------------------------------------------- output
    def kspace(self) -> dict[str, np.ndarray]:
        """
        Return the k-space trajectory, in 1/m.

        Returns
        -------
        dict
            ``k_adc`` (3 x n_samples, at the ADC sample times), ``t_adc``, ``k`` (dense),
            ``t_k``, ``t_excitation``, ``t_refocusing``.

        Notes
        -----
        ``calculate_kspacePP`` returns its tuple in a different order from
        ``calculate_kspace``; getting that wrong silently swaps the trajectory for its
        timebase.
        """
        k_adc, t_adc, k, t_k, t_exc, t_refoc = self.seq.calculate_kspacePP()[:6]
        return {
            'k_adc': np.asarray(k_adc),
            't_adc': np.asarray(t_adc),
            'k': np.asarray(k),
            't_k': np.asarray(t_k),
            't_excitation': np.asarray(t_exc),
            't_refocusing': np.asarray(t_refoc),
        }

    def pns(self, hardware: SimpleNamespace) -> Report:
        """
        Predict peripheral nerve stimulation against a gradient hardware model.

        Parameters
        ----------
        hardware
            The gradient response model, from :func:`seqcraft.hardware.load_hardware` or
            :func:`seqcraft.hardware.synthetic_hardware`.  **Required**: PNS prediction is
            analysis rather than compilation, and the model describes an amplifier's response
            rather than its limits, so it has nothing to do with the ``Opts`` this was compiled
            against and is not carried on it.  The synthetic model is a vendor-free stand-in for
            CI; it is **not** a real scanner and must never be used to clear a human scan.
        """
        model = hardware
        ok, pns_norm, _components, _t = self.seq.calculate_pns(model, do_plots=False)
        peak = float(np.max(pns_norm)) if np.size(pns_norm) else 0.0
        note = (
            ' (synthetic model -- not valid for clearing a human scan)'
            if getattr(model, 'is_synthetic', False)
            else ''
        )
        return Report(
            (
                Issue(
                    'pns',
                    'sequence',
                    f'peak PNS {peak * 100:.0f}% of the stimulation limit{note}',
                    'info' if ok else 'error',
                ),
            ),
            values={'peak_pns_fraction': peak},
        )

    def write(self, path: str | Path, *, sidecar: bool = True) -> WriteResult:
        """
        Write the ``.seq`` file, and by default a JSON provenance sidecar beside it.

        Takes no geometry, matrix or FOV arguments: everything written comes from what was
        compiled, so the file's metadata cannot disagree with what it plays.

        Parameters
        ----------
        path
            Destination ``.seq`` path.
        sidecar
            Also write ``<path>.json`` recording versions, git state, the definitions, the
            achieved duration and the file's sha256.
        """
        import hashlib  # noqa: PLC0415
        from pathlib import Path as _Path  # noqa: PLC0415

        target = _Path(path)

        for key, value in self.definitions.items():
            self.seq.set_definition(key, value)
        self.seq.set_definition('TotalDuration', self.duration_s)
        self.seq.write(str(target))
        digest = hashlib.sha256(target.read_bytes()).hexdigest()

        side: _Path | None = None
        if sidecar:
            side = write_sidecar(
                target,
                {
                    'definitions': {k: _jsonable(v) for k, v in self.definitions.items()},
                    # Every Opts field is a plain float, int, bool or None, so the scanner the
                    # sequence was built against records itself with no schema in between.
                    'opts': dict(vars(self.opts)),
                    'n_blocks': self.n_blocks,
                    'duration_s': self.duration_s,
                    'sha256': digest,
                    'issues': [
                        {'kind': i.kind, 'where': i.where, 'message': i.message,
                         'severity': i.severity}
                        for i in self.report.issues
                    ],
                },
            )
        return WriteResult(
            path=target,
            sha256=digest,
            n_blocks=self.n_blocks,
            duration_s=self.duration_s,
            sidecar=side,
        )
