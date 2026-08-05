"""
Findings as data: :class:`Issue` and :class:`Report`.

Hard failures raise (see :mod:`seqcraft.core.errors`).  Soft findings -- a gradient at
98 % of the limit, a b-value 0.4 % below target after raster rounding -- become a
:class:`Report`, which is a value: it can be asserted on, rendered as a table, written into
the provenance sidecar, or turned into an exception on demand.

The behaviour this replaces is ``pSeq_Base.get_report()``::

    ok, error_report = self.seq.check_timing()
    if ok:
        print("Timing check passed successfully")
    else:
        print("Timing check failed. Error listing follows:")
        [print(e) for e in error_report]

That printed and returned ``None``, so a failing timing check could not be detected
programmatically -- and in a notebook the message scrolled away.

Constructing a :class:`Report` has **no side effects**: it never prints and never draws.

Examples
--------
>>> r = Report((Issue('raster', 'EPIReadout.blip', 'duration off raster'),), subject='epi')
>>> r.ok
False
>>> print(r)
epi: 1 issue
  [raster] EPIReadout.blip: duration off raster
>>> Report((), subject='gre').ok
True
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal, NamedTuple

from .errors import SeqCraftError

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping, Sequence

__all__ = ['Issue', 'Report', 'ReportFailed']

Severity = Literal['error', 'warning', 'info']


class ReportFailed(SeqCraftError):
    """Raised by :meth:`Report.raise_if_failed` when a report contains errors."""


class Issue(NamedTuple):
    """
    One finding.

    Parameters
    ----------
    kind
        Short machine-readable category: ``'raster'``, ``'grad_limit'``, ``'slew_limit'``,
        ``'rf'``, ``'adc'``, ``'anchor'``, ``'block'``, ``'label'``, ``'definition'``,
        ``'timing'``, ``'bvalue'``.
    where
        Location, normally ``'<module>.<block tag>'`` or a block index.
    message
        One sentence, including the offending number and its unit.
    severity
        ``'error'`` by default.  Warnings do not make :attr:`Report.ok` false.
    data
        Optional machine-readable detail, carried into the sidecar.
    """

    kind: str
    where: str
    message: str
    severity: Severity = 'error'
    data: Mapping[str, Any] | None = None


@dataclass(frozen=True)
class Report:
    """
    A collection of :class:`Issue` objects plus optional measured values.

    Parameters
    ----------
    issues
        The findings.
    subject
        What was checked, used in the rendered heading.
    values
        Measured quantities worth reporting even when nothing is wrong: achieved TE and
        TR, peak gradient per axis, b-value, block count.

    Notes
    -----
    :attr:`ok` is ``True`` when there are no ``'error'``-severity issues; warnings are
    reported but do not fail.  That split is what lets a build succeed while still
    recording, for example, that raster rounding moved the b-value by 0.4 %.
    """

    issues: tuple[Issue, ...] = ()
    subject: str = ''
    values: Mapping[str, Any] = field(default_factory=dict)

    # ------------------------------------------------------------------------ status
    @property
    def ok(self) -> bool:
        """``True`` when the report contains no errors (warnings are allowed)."""
        return not self.errors

    @property
    def errors(self) -> tuple[Issue, ...]:
        """Issues of severity ``'error'``."""
        return tuple(i for i in self.issues if i.severity == 'error')

    @property
    def warnings(self) -> tuple[Issue, ...]:
        """Issues of severity ``'warning'``."""
        return tuple(i for i in self.issues if i.severity == 'warning')

    def of_kind(self, kind: str) -> tuple[Issue, ...]:
        """Issues whose :attr:`Issue.kind` equals `kind`."""
        return tuple(i for i in self.issues if i.kind == kind)

    # --------------------------------------------------------------------- combining
    def merge(self, *others: Report, subject: str | None = None) -> Report:
        """Return a new report combining this one with `others`."""
        issues = list(self.issues)
        values = dict(self.values)
        for other in others:
            issues.extend(other.issues)
            values.update(other.values)
        return Report(tuple(issues), subject=subject or self.subject, values=values)

    @classmethod
    def combine(cls, reports: Iterable[Report], *, subject: str = '') -> Report:
        """Combine an iterable of reports into one."""
        reports = list(reports)
        if not reports:
            return cls((), subject=subject)
        return reports[0].merge(*reports[1:], subject=subject or reports[0].subject)

    # ---------------------------------------------------------------------- raising
    def raise_if_failed(self) -> Report:
        """
        Raise :class:`ReportFailed` if this report has errors; otherwise return ``self``.

        Returning ``self`` makes ``report = seq.check().raise_if_failed()`` read naturally.
        """
        if self.ok:
            return self
        raise ReportFailed(str(self))

    # -------------------------------------------------------------------- rendering
    def __str__(self) -> str:
        """Plain-text rendering, for logs and the CLI."""
        head = self.subject or 'report'
        n = len(self.issues)
        if not n:
            line = f'{head}: ok'
            return line + self._values_text()
        plural = 's' if n != 1 else ''
        lines = [f'{head}: {n} issue{plural}']
        for i in self.issues:
            mark = '' if i.severity == 'error' else f' ({i.severity})'
            lines.append(f'  [{i.kind}] {i.where}: {i.message}{mark}')
        return '\n'.join(lines) + self._values_text()

    def _values_text(self) -> str:
        if not self.values:
            return ''
        lines = [''] + [f'  {k} = {_fmt(v)}' for k, v in self.values.items()]
        return '\n'.join(lines)

    def _repr_html_(self) -> str:
        """HTML table, for Jupyter.  Never raises."""
        try:
            rows = ''.join(
                f'<tr><td style="padding:2px 8px"><code>{i.kind}</code></td>'
                f'<td style="padding:2px 8px"><code>{i.where}</code></td>'
                f'<td style="padding:2px 8px">{i.message}</td>'
                f'<td style="padding:2px 8px">{i.severity}</td></tr>'
                for i in self.issues
            )
            vals = ''.join(
                f'<tr><td style="padding:2px 8px"><b>{k}</b></td>'
                f'<td colspan="3" style="padding:2px 8px">{_fmt(v)}</td></tr>'
                for k, v in self.values.items()
            )
            status = 'ok' if self.ok else f'{len(self.errors)} error(s)'
            body = rows + vals or '<tr><td style="padding:2px 8px">nothing to report</td></tr>'
            return (
                f'<div><b>{self.subject or "report"}</b> &mdash; {status}'
                f'<table style="border-collapse:collapse;font-size:90%">{body}</table></div>'
            )
        except Exception as exc:  # noqa: BLE001 - a raising repr makes a notebook unusable
            return f'<pre>Report could not render: {exc!r}</pre>'

    def to_dict(self) -> dict[str, Any]:
        """JSON-safe representation, for the provenance sidecar."""
        return {
            'subject': self.subject,
            'ok': self.ok,
            'values': dict(self.values),
            'issues': [
                {
                    'kind': i.kind,
                    'where': i.where,
                    'message': i.message,
                    'severity': i.severity,
                    'data': dict(i.data) if i.data else None,
                }
                for i in self.issues
            ],
        }


def _fmt(value: Any) -> str:
    """Render a report value compactly."""
    if isinstance(value, float):
        return f'{value:.6g}'
    if isinstance(value, (list, tuple)) and value and isinstance(value[0], float):
        return '[' + ', '.join(f'{v:.6g}' for v in value) + ']'
    return str(value)
