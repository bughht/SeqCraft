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

from ..design import events as ev
from ..design.events import ADDRESS_KEYS, AXES
from ..report import Issue, Report, ReportFailed
from .provenance import _jsonable, write_sidecar

if TYPE_CHECKING:
    from collections.abc import Sequence
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
    report
        Everything the compile found: every same-axis merge, every limit violation, every
        snapped time.
    origins
        One tag path per compiled block, so a block index traces back to the module that
        produced it.
    definitions
        The ``[DEFINITIONS]`` that will be written, already merged and collision-checked.  Plain
        pulseq keys: there is no geometry object here, because everything the file says about the
        scan came in as this mapping.
    """

    seq: Any
    opts: Opts
    report: Report
    origins: tuple[tuple[str, ...], ...]
    definitions: dict[str, Any]
    tree_duration_s: float
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

    def moments(self, order: int = 0) -> dict[str, float]:
        """
        Return the whole-sequence gradient moment per axis, integrated from the compiled blocks.

        Parameters
        ----------
        order
            ``0`` for area in 1/m, ``1`` for s/m, ``2`` for s^2/m.  Referenced to the start of
            the sequence.

        Notes
        -----
        Integrated from each block's exact knots, not from raster samples.  The difference is
        not cosmetic: sampling is exact for ``order == 0`` and only then, and an arbitrary
        gradient's samples sit at raster *centres*, so a raster-sampled m0 quietly matched a
        raster-sampled tree even when the compiled waveform had moved.
        """
        out: dict[str, float] = dict.fromkeys(AXES, 0.0)
        t = 0.0
        # Block IDs are 1-based, and block_durations is a dict keyed by them, not a list.
        for index in sorted(self.seq.block_events):
            block = self.seq.get_block(index)
            for axis in AXES:
                grad = getattr(block, f'g{axis}', None)
                if grad is not None:
                    out[axis] += ev.pwl_moment(*ev.knots_of(grad, t), order)
            t += float(self.seq.block_durations[index])
        return out

    def check(self, *, allow_timing: Sequence[str] = ('TotalDuration',)) -> Report:
        """
        Run every post-compile check and return one report.

        Combines the compile report with ``Sequence.check_timing`` and label-address
        uniqueness.

        Parameters
        ----------
        allow_timing
            Substrings of ``check_timing`` messages to downgrade to information.  Defaults to
            the ``TotalDuration`` float-equality artifact, which pypulseq emits even on
            pulseq's own approved reference files.

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
        issues = list(self.report.issues)
        ok, errors = self.seq.check_timing()
        if not ok:
            for line in errors:
                text = str(line).strip()
                allowed = any(token in text for token in allow_timing)
                issues.append(Issue('timing', 'sequence', text, 'info' if allowed else 'error'))
        issues.extend(self._label_issues())
        issues.extend(self._event_size_issues())
        out = Report(tuple(issues), subject=self.report.subject, values={
            'n_blocks': self.n_blocks,
            'duration_s': self.duration_s,
        })
        object.__setattr__(self, '_checked', out)
        return out

    def _event_size_issues(self) -> list[Issue]:
        """
        Check every ADC and RF event against the interpreter's per-event sample limits.

        These limits live in ``Opts`` as ``adc_samples_limit`` and ``rf_samples_limit`` and default
        to ``0``, pypulseq's "no limit".  Nothing checked them until a 67 388-sample spiral readout
        reached a scanner, which refused the block with ``fRTEBFinish() failed for block type:
        ArbX ArbY ADC`` -- a message that names the block type and says nothing about samples.

        The limit is the vendor interpreter's, not the amplifier's, so it has to be set from the
        installation: :func:`seqcraft.scanner.opts.from_scanner` takes ``adc_samples_limit=`` for
        exactly this, and 8192 is the common Siemens value.  A readout longer than one event's worth has
        to be split into several ADCs, at the cost of ``adc_dead_time`` between them.
        """
        out: list[Issue] = []
        limits = self.opts
        for kind, attribute, count in (
            ('adc', 'adc_samples_limit', lambda block: int(block.adc.num_samples)),
            ('rf', 'rf_samples_limit', lambda block: int(np.size(block.rf.signal))),
        ):
            limit = int(getattr(limits, attribute, 0) or 0)
            if limit <= 0:
                continue
            worst = 0
            where = ''
            for index in sorted(self.seq.block_events):
                block = self.seq.get_block(index)
                if getattr(block, kind, None) is None:
                    continue
                samples = count(block)
                if samples > worst:
                    worst, where = samples, f'block {index} ({self.origin(index)})'
            if worst > limit:
                out.append(Issue(
                    f'{kind}_samples_limit',
                    where,
                    f'{worst} {kind.upper()} samples in one event, above the {limit} the '
                    f'interpreter accepts; split it into '
                    f'{-(-worst // limit)} events or lengthen the dwell',
                    'error',
                ))
        return out

    def _label_issues(self) -> list[Issue]:
        """
        Check that no two imaging ADCs write the same k-space address.

        The highest-value check available on a finished sequence: a duplicate address means two
        readouts landing in the same place, which catches a wrong slice order, an off-by-one
        partial-Fourier start and a mis-nested loop from one assertion.
        """
        try:
            labels = self.seq.evaluate_labels(evolution='adc')
        except (AttributeError, ValueError, IndexError):  # pragma: no cover - older pypulseq
            return []
        keys = [k for k in ADDRESS_KEYS if k in labels]
        if not keys:
            return []
        arrays = [np.atleast_1d(np.asarray(labels[k])) for k in keys]
        n = max(len(a) for a in arrays)
        arrays = [np.resize(a, n) for a in arrays]
        skip = np.zeros(n, dtype=bool)
        for flag in ('NOISE', 'REF', 'NAV'):
            if flag in labels:
                skip |= np.resize(np.atleast_1d(np.asarray(labels[flag])).astype(bool), n)
        addresses = [tuple(int(a[i]) for a in arrays) for i in range(n) if not skip[i]]
        duplicates = len(addresses) - len(set(addresses))
        if not duplicates:
            return []
        seen: set[tuple[int, ...]] = set()
        first: tuple[int, ...] = ()
        for address in addresses:
            if address in seen:
                first = address
                break
            seen.add(address)
        return [
            Issue(
                'label',
                'sequence',
                f'{duplicates} imaging ADC(s) repeat a k-space address; first repeat is '
                f'{dict(zip(keys, first))} -- two readouts are writing the same location',
                'error',
            )
        ]

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
