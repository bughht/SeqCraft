r"""
:class:`FSE2D` -- a complete 2D turbo-spin-echo scan: shots, and which lines each one acquires.

``imaging/`` because it composes a kernel, and because everything it owns is **acquisition
policy**: how many shots there are, which ``ky`` each echo of each shot acquires, how many dummies
come first, and therefore which echo carries the centre of k-space and what the effective TE is.
The physics -- the echo spacing, the crusher window, the moment balance about every refocusing
pulse -- belongs to :class:`~seqcraft.modules.TSEShot` and is not repeated here.

``segments`` is data, not a pattern
-----------------------------------
A list of line lists, one per shot.  That single argument carries the turbo factor, the effective
TE and the ordering, and no generator ships for it: ``writeTSE.m``'s own reordering ``circshift``
is commented *"TSE echo time magic"*, which is a protocol choice rather than arithmetic.  Three
orderings out of one instance::

    interleaved   shot s takes lines s, s+shots, s+2*shots ...   each echo owns a band of ky
    linear        shot s takes a contiguous block                 each echo is a comb across ky
    centric       lines dealt outward from the centre             the earliest echoes near k = 0

The same 128 lines and the same waveforms in all three; what differs is which T2 weighting lands
where in k-space, and :meth:`echo_bands` measures it.

The failure mode a train has that a single shot does not
--------------------------------------------------------
Shot 0 starts from thermal equilibrium and every later shot starts from whatever the last TR
recovered to.  With interleaved segmentation each echo index owns a *comb* in k-space, so that
one difference is a periodic modulation -- and a periodic modulation of k-space is a **replica**
of the object, not a blur.  `dummy_shots` is the fix, and it is why the default is not obviously
zero for a multi-shot protocol.

HASTE is this class
-------------------
One shot, a long train and ``partial_fourier`` below 1 is what a vendor calls HASTE; ``echoes=1``
is a conventional spin echo.  Neither needs a class of its own, and the fine scan that produced
this module found no physics in either that this composition does not express.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

from ...design.logic import LogicBlock
from ...design.module import Module
from ...errors import ConfigurationError, SeqCraftWarning, format_error
from ..kernel.tse_shot import TSEShot

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from pypulseq.opts import Opts

__all__ = ['FSE2D']


class FSE2D(Module):
    """
    A complete 2D turbo spin echo: one shot per segment, stacked at ``TR``.

    Parameters
    ----------
    opts, fov_mm, matrix, thickness_mm, echoes, echo_spacing_s, tr_s, flip_deg,
    refocus_flip_deg, bandwidth_hz_px, partial_fourier, excitation_duration_s,
    refocus_duration_s, refocus_thickness_factor, crush_cycles_slice, crush_cycles_readout,
    spoil_cycles_per_voxel, spoil_axis
        Forwarded unchanged to :class:`~seqcraft.modules.TSEShot`, which owns all of them.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    shot : TSEShot
        The repeating unit.  Exposed so that ``fse.shot.te_s(3)`` or ``fse.shot.echo_spacing_s``
        can be read without rebuilding one.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s',
    ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> fse = FSE2D(opts=o, fov_mm=256.0, matrix=(128, 128), thickness_mm=5.0, echoes=16,
    ...             tr_s=2.0)
    >>> segments = [[s + n * 8 for n in range(16)] for s in range(8)]
    >>> fse.echo_of_center_line(segments)
    (0, 8)
    >>> round(fse.te_eff_s(segments) * 1e3, 3)
    96.302
    """

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: float | tuple[float, float],
        matrix: tuple[int, int],
        thickness_mm: float,
        echoes: int = 1,
        echo_spacing_s: float | None = None,
        tr_s: float | None = None,
        flip_deg: float = 90.0,
        refocus_flip_deg: float = 180.0,
        bandwidth_hz_px: float = 200.0,
        partial_fourier: float = 1.0,
        excitation_duration_s: float = 2.5e-3,
        refocus_duration_s: float = 4e-3,
        refocus_thickness_factor: float = 1.25,
        crush_cycles_slice: float = 3.0,
        crush_cycles_readout: float = 0.0,
        spoil_cycles_per_voxel: float = 4.0,
        spoil_axis: str | Iterable[str] = ('x', 'z'),
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.shot = TSEShot(
            opts=opts, fov_mm=fov_mm, matrix=matrix, thickness_mm=thickness_mm, echoes=echoes,
            echo_spacing_s=echo_spacing_s, tr_s=tr_s, flip_deg=flip_deg,
            refocus_flip_deg=refocus_flip_deg, bandwidth_hz_px=bandwidth_hz_px,
            partial_fourier=partial_fourier, excitation_duration_s=excitation_duration_s,
            refocus_duration_s=refocus_duration_s,
            refocus_thickness_factor=refocus_thickness_factor,
            crush_cycles_slice=crush_cycles_slice, crush_cycles_readout=crush_cycles_readout,
            spoil_cycles_per_voxel=spoil_cycles_per_voxel, spoil_axis=spoil_axis,
        )

    # ------------------------------------------------------------------ what it knows
    @property
    def echoes(self) -> int:
        """Echoes per shot -- the turbo factor.  Read from the shot."""
        return self.shot.echoes

    @property
    def center_line(self) -> int:
        """The phase-encode line that encodes ``k = 0``.  Read from the shot."""
        return self.shot.center_line

    # These take `segments` rather than reading `self`, mirroring `build`, because the answer
    # moves with the table -- the same reason `GRE2D.time_to_center_line` is a method.
    def echo_of_center_line(self, segments: Sequence[Sequence[int]]) -> tuple[int, int] | None:
        """``(shot, echo)`` of the readout that acquires :attr:`center_line`, or ``None``."""
        for shot, segment in enumerate(segments):
            lines = list(segment)
            if self.center_line in lines:
                return shot, lines.index(self.center_line)
        return None

    def te_eff_s(self, segments: Sequence[Sequence[int]]) -> float | None:
        """
        The effective echo time: the TE of whichever echo carries :attr:`center_line`.

        ``None`` when no segment acquires it, which is a legitimate table and an image with no DC
        term -- so this reports rather than raises, and :meth:`build` warns.
        """
        found = self.echo_of_center_line(segments)
        return None if found is None else self.shot.te_s(found[1])

    def echo_bands(self, segments: Sequence[Sequence[int]]) -> dict[int, int]:
        r"""
        ``echo -> span`` of ``|ky - center_line|`` across shots, in lines.

        Distance from the centre of k-space rather than signed ``ky``, and the difference is the
        whole usefulness of the number: a **centric** ordering deliberately puts its last echo at
        *both* edges, which is a span of the whole matrix in signed ``ky`` and a narrow band in
        distance.  Its point spread is symmetric and it blurs; a table whose one echo index is
        scattered over a *range of distances* modulates adjacent rows differently and ghosts.
        """
        table = [list(segment) for segment in segments]
        return {
            echo: max(abs(s[echo] - self.center_line) for s in table)
            - min(abs(s[echo] - self.center_line) for s in table)
            for echo in range(self.echoes)
        }

    # ------------------------------------------------------------------------- assembly
    def build(
        self,
        *,
        segments: Sequence[Sequence[int]],
        dummy_shots: int = 0,
        center_mm: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> LogicBlock:
        """
        Return the whole scan: `dummy_shots` that sample nothing, then one shot per segment.

        Parameters
        ----------
        segments
            One list of phase-encode indices per shot, each :attr:`echoes` long, in echo order.
        dummy_shots
            Shots played **before** the first acquired one, with ``acquire=False``, replaying the
            first segment's lines so that what they establish is the steady state the first real
            shot then samples.
        center_mm
            ``(x, y, z)`` centre of the imaging volume, millimetres, forwarded unchanged.
        """
        table = self._check(segments)
        count = self._require_dummies(dummy_shots)
        tr_s = self.shot.tr_s
        out = LogicBlock()
        shots = [table[0]] * count + list(table)
        for n, lines in enumerate(shots):
            # The segment index is this layer's to know -- a shot cannot count its siblings -- and
            # the shot's to place, because it knows where its first readout starts.
            out.add(n * tr_s, self.shot(lines=lines, acquire=n >= count,
                                        segment_index=max(0, n - count), center_mm=center_mm))
        return out

    # ---------------------------------------------------------------------- the refusals
    def _require_dummies(self, dummy_shots: int) -> int:
        """Return `dummy_shots` having checked it is not negative."""
        count = int(dummy_shots)
        if count < 0:
            msg = format_error(f'dummy_shots = {count} is negative.',
                               {'dummy_shots': dummy_shots}, [])
            raise ConfigurationError(msg)
        return count

    def _check(self, segments: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
        """
        What the compiler cannot see about a table.

        The band check **warns** rather than refuses, because a table that scatters one echo index
        across k-space is a legitimate experiment as well as a common mistake -- and the symptom,
        ghosting along the phase-encode direction, reads as motion rather than as a table problem.
        """
        table = tuple(tuple(int(line) for line in segment) for segment in segments)
        if not table or not all(table):
            msg = format_error('segments must be a non-empty list of non-empty shots.',
                               {'segments': table[:4]}, ['one list of lines per shot'])
            raise ConfigurationError(msg)
        # Length and range are the shot's own contract; ask it, so one refusal serves both.
        for segment in table:
            self.shot._check(segment)                                            # noqa: SLF001

        seen = [line for segment in table for line in segment]
        repeated = sorted({line for line in seen if seen.count(line) > 1})
        if repeated:
            msg = format_error(
                f'{len(repeated)} line(s) are acquired more than once: {repeated[:8]}.',
                {'shots': len(table), 'echoes': self.echoes},
                ['two readouts writing one k-space address is what the compiler rejects '
                 'downstream, with less to say about it than here'])
            raise ConfigurationError(msg)
        if self.center_line not in seen:
            warnings.warn(
                f'no segment acquires the centre line {self.center_line}, so the image has no DC '
                f'term and te_eff_s is undefined.', SeqCraftWarning, stacklevel=3)
        widest = max(1, self.matrix_y // self.echoes)
        wide = {echo: span for echo, span in self.echo_bands(table).items() if span > widest}
        if wide:
            worst = max(wide, key=lambda echo: wide[echo])
            warnings.warn(
                f'{len(wide)} echo index(es) fill widely separated ky: echo {worst} spans '
                f'{wide[worst]} lines of distance from the centre of k-space, against the '
                f'{widest} a contiguous band would.  Adjacent k-space rows then carry T2 '
                f'weightings many echoes apart, and the symptom is ghosting along the '
                f'phase-encode direction that reads as motion.  Not refused: a deliberate '
                f'experiment is a legitimate thing to build.', SeqCraftWarning, stacklevel=3)
        return table

    @property
    def matrix_y(self) -> int:
        """The phase-encode matrix size.  Read from the shot."""
        return self.shot.matrix[1]
