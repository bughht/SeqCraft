"""
:class:`GRE2D` -- a complete spoiled 2D gradient-echo scan.

``imaging/`` because it composes kernels, and ``imaging`` rather than ``sequence`` because
"sequence" already names ``pypulseq.Sequence`` and ``mr0.Sequence``, both of which appear in the
same notebooks; a third meaning of the word would be read as one of those.

The scan is one :class:`~seqcraft.modules.GRE2DTR` stacked at ``n * tr_s``.  What only a whole
acquisition can own lives here: the loop, the dummy repetitions, the RF-spoiling schedule, and
the validation of which lines are being acquired.

The sampling pattern is a build argument
----------------------------------------
``lines`` is the phase-encode indices in acquisition order -- a plain list, passed to ``build``
rather than to ``__init__``.  One instance therefore produces any number of patterns, which is
the *``__init__`` designs / ``build`` assembles* split doing its job: the waveforms are created
once and a three-way comparison is three cheap calls on one object::

    gre(lines=range(ny))                               # fully sampled
    gre(lines=sorted(acs | set(range(0, ny, 3))))      # R = 3 with a calibration block
    gre(lines=reversed(range(ny)))                     # ordering, with no new argument
    gre(lines=range(0, ny, 2))                         # shot 0 of a two-shot interleave

That single argument replaces acceleration, phase-encode partial Fourier, multi-shot and
ordering, because every one of them is a different list.  There is **no default**: fully sampled
is ``range(matrix[1])`` written out, so undersampling reads as a peer choice rather than as a
deviation.

``dummies`` is a build argument for the same reason.  It designs no waveform -- ``build`` plays
the first acquired line's gradients with ``acquire=False`` -- and it describes *what this
acquisition does*, exactly as ``lines`` does.  Three things follow: the timing query
:meth:`GRE2D.time_to_center_line` can mirror ``build`` rather than reading one of its two terms
off ``self``; one instance can be compared with and without them, which is a two-call experiment
rather than two objects; and a caller whose steady state is established somewhere else -- an
MPRAGE shot loop, where the equilibrium is shot-to-shot -- leaves the default alone instead of
re-passing zero to a constructor.

**No pattern generators ship.**  Which lines to acquire is a sequence-programming choice, which
ADR-003 keeps out of the library -- the reason ``salvage/ordering.py`` and
``salvage/geometry_pe.py`` are quarantined rather than imported.  A list of integers is data; a
function that produces one is policy, and shipping the second would give back exactly what
taking a list was meant to avoid.
"""

from __future__ import annotations

import warnings
from collections import Counter
from typing import TYPE_CHECKING

from ...design.logic import LogicBlock
from ...design.module import Module
from ...errors import ConfigurationError, SeqCraftWarning, format_error
from .._support import require_positive
from ..kernel.gre_2d_tr import GRE2DTR

if TYPE_CHECKING:
    from collections.abc import Iterable

    from pypulseq.opts import Opts

__all__ = ['GRE2D']

#: The quadratic phase increment that spoils best for a wide range of flip angles and TRs.
_DEFAULT_RF_SPOIL_DEG = 117.0


class GRE2D(Module):
    """
    A complete spoiled 2D gradient-echo acquisition.

    Parameters
    ----------
    opts, fov_mm, matrix, thickness_mm, flip_deg, te_s, tr_s, bandwidth_hz_px,
    partial_fourier, spoil_cycles_per_voxel, spoil_axis
        Forwarded unchanged to :class:`~seqcraft.modules.GRE2DTR`, which is the one
        :class:`~seqcraft.Module` this constructs and holds.  `partial_fourier` stays an
        ``__init__`` argument because it is the *readout* direction: it changes the prephaser
        area and the echo time inside the readout, which no list of phase-encode indices can
        express, and it is designed once.
    rf_spoil
        Quadratic phase increment on the excitation.  Default on.
    rf_spoil_deg
        The increment, in degrees.  ``None`` resolves to ``117.0`` when `rf_spoil` is on, and
        passing it with ``rf_spoil=False`` raises -- the same shape as `axis` alongside a
        non-selective :class:`~seqcraft.modules.Excitation`, so an argument that cannot take
        effect is reported rather than ignored.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    tr : GRE2DTR
        The repetition this scan stacks.  Readable, so ``gre.tr.te_s`` and
        ``gre.tr.ro.dwell_s`` answer without a second copy of any number.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
    ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> gre = GRE2D(opts=o, fov_mm=250.0, matrix=(64, 32), thickness_mm=5.0)
    >>> scan = gre(lines=range(32))
    >>> len(scan), abs(scan.duration - 32 * gre.tr.tr_s) < 1e-9
    (32, True)

    One instance, three patterns, three cheap calls:

    >>> acs = set(range(gre.center_line - 4, gre.center_line + 4))
    >>> len(gre(lines=sorted(acs | set(range(0, 32, 3)))))
    16

    And with and without dummies is one instance too, which is what makes showing what they fix
    an experiment rather than a second object:

    >>> len(gre(lines=range(32), dummies=8))
    40

    Notes
    -----
    **The spoiling phase counter runs across the dummies.**  Repetition *n* gets
    ``0.5 * rf_spoil_deg * n * (n + 1)`` with *n* counted from the first dummy, not from the
    first acquired line.  Restarting it at the acquisition is a real and quiet bug: the steady
    state the dummies established is a steady state *of a particular phase sequence*, and
    resetting the counter discards it at the moment it starts to matter.

    With ``rf_spoil=False`` every repetition gets ``phase_deg=0.0``.  Gradient spoiling is
    independent and unaffected -- the two are separate mechanisms, and turning off one while
    keeping the other is a legitimate thing to want to look at.

    **Validation is per call, and it is this module's whole responsibility over the pattern.**
    Out-of-range or duplicate indices raise.  A pattern **missing** :attr:`center_line` warns:
    not illegal, but an image reconstructed with no DC term reads as a windowing problem rather
    than as a sampling bug, and with no generators shipped that warning is the only thing
    standing between a caller and the failure ``salvage/geometry_pe.py``'s residue nudge existed
    to prevent.

    ``2D`` is in the name because ky-kz sampling is where 3D complexity lives -- elliptical
    shutter, CAIPIRINHA, radial or spiral ordering in the ky-kz plane.  ``GRE3D`` is written when
    MPRAGE needs it; whether the two merge is decided then, with both in front of us.
    """

    def __init__(
        self,
        *,
        opts: Opts,
        fov_mm: float | tuple[float, float],
        matrix: tuple[int, int],
        thickness_mm: float,
        flip_deg: float = 15.0,
        te_s: float | None = None,
        tr_s: float | None = None,
        bandwidth_hz_px: float = 200.0,
        partial_fourier: float = 1.0,
        rf_spoil: bool = True,
        rf_spoil_deg: float | None = None,
        spoil_cycles_per_voxel: float = 4.0,
        spoil_axis: str | Iterable[str] = ('x', 'z'),
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.tr = GRE2DTR(
            opts=opts, fov_mm=fov_mm, matrix=matrix, thickness_mm=thickness_mm,
            flip_deg=flip_deg, te_s=te_s, tr_s=tr_s, bandwidth_hz_px=bandwidth_hz_px,
            partial_fourier=partial_fourier, spoil_cycles_per_voxel=spoil_cycles_per_voxel,
            spoil_axis=spoil_axis,
        )
        self.rf_spoil = bool(rf_spoil)
        self.rf_spoil_deg = self._resolve_spoil_phase(rf_spoil_deg)

    # ------------------------------------------------------------------ what it knows
    @property
    def center_line(self) -> int:
        """The phase-encode line that encodes ``k = 0``.  Read from the repetition."""
        return self.tr.center_line

    def phase_deg(self, repetition: int) -> float:
        """
        Return the RF carrier phase for `repetition`, degrees, counted from the first dummy.

        Examples
        --------
        >>> import pypulseq as pp
        >>> from pypulseq.opts import Opts
        >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
        ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
        >>> gre = GRE2D(opts=o, fov_mm=250.0, matrix=(64, 32), thickness_mm=5.0)
        >>> [gre.phase_deg(n) for n in range(4)]
        [0.0, 117.0, 351.0, 702.0]
        """
        if not self.rf_spoil:
            return 0.0
        n = int(repetition)
        return 0.5 * self.rf_spoil_deg * n * (n + 1)

    def time_to_center_line(self, *, lines: Iterable[int], dummies: int = 0) -> float:
        """
        Seconds from the start of the block `build` returns to the acquisition of ``k = 0``.

        Parameters
        ----------
        lines, dummies
            **The same arguments** :meth:`build` takes, because the answer moves with both.

        Returns
        -------
        float
            ``(dummies + table.index(center_line)) * tr_s + tr.time_to_echo()``.

        Raises
        ------
        ConfigurationError
            On any table :meth:`build` refuses, and on one that omits :attr:`center_line` -- there
            is no answer then, and returning one would be worse than refusing.

        Examples
        --------
        >>> import pypulseq as pp
        >>> from pypulseq.opts import Opts
        >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
        ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
        >>> gre = GRE2D(opts=o, fov_mm=250.0, matrix=(64, 32), thickness_mm=5.0)
        >>> centric = [16, 15, 17, 14, 18]        # k = 0 first
        >>> linear = sorted(centric)              # k = 0 mid-train
        >>> round(gre.time_to_center_line(lines=centric) * 1e3, 3)
        6.218
        >>> round(gre.time_to_center_line(lines=linear) * 1e3, 3)
        25.358

        The same lines, two orderings, two answers -- which is why this is a method rather than a
        constant.  And the dummies move it by exactly their own duration:

        >>> gap = gre.time_to_center_line(lines=linear, dummies=20) - \\
        ...       gre.time_to_center_line(lines=linear)
        >>> abs(gap - 20 * gre.tr.tr_s) < 1e-12
        True

        Notes
        -----
        **This is what makes TI mean anything.**  An inversion time is the inversion's effective
        centre to the acquisition of ``k = 0``, so placing a train requires knowing where in it the
        centre falls -- which depends on the *ordering* and is therefore not a constant.  Centric
        puts ``k = 0`` first and linear puts it mid-train, and the same TI places those two trains
        hundreds of milliseconds apart.

        **The `dummies` term is not optional and is easy to drop.**  The block starts at the first
        dummy, not at the first acquired line, so omitting it under-reports by ``dummies * tr_s``
        -- with twenty dummies at a 20 ms TR that is 400 ms of TI error, silently, in the one
        number the sequence exists to control.  That is why this signature mirrors ``build``'s: the
        arguments that move the answer are the arguments that build the block.
        """
        table = self._check(lines)
        if self.center_line not in table:
            msg = format_error(
                f'lines omits the centre of k-space (line {self.center_line}), so there is no '
                f'k = 0 to measure to.',
                {'center_line': self.center_line, 'lines': table[:8]},
                [
                    'include the centre line in this train',
                    'under segmentation only the shot that carries it has a meaningful TI, and '
                    'this is the query that says which one that is',
                ],
            )
            raise ConfigurationError(msg)
        repetition = _require_dummies(dummies) + table.index(self.center_line)
        return repetition * self.tr.tr_s + self.tr.time_to_echo()

    # ----------------------------------------------------------------------- assembly
    def build(
        self,
        *,
        lines: Iterable[int],
        dummies: int = 0,
        center_mm: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> LogicBlock:
        """
        Return the whole scan: `dummies` repetitions, then one per entry of `lines`.

        Parameters
        ----------
        lines
            Phase-encode indices, **in acquisition order**.  Any iterable of ints; ``range`` and
            ``reversed`` both work, and the order is the acquisition order.
        dummies
            Repetitions played **before** the first acquired line, with ``acquire=False``, loading
            the first acquired line's gradients.  Describes this acquisition rather than any
            waveform, which is why it is here beside the pattern it precedes rather than on
            ``__init__``.
        center_mm
            ``(x, y, z)`` centre of the imaging volume, millimetres, forwarded unchanged.
        """
        table = self._check(lines)
        count = _require_dummies(dummies)
        tr_s = self.tr.tr_s
        out = LogicBlock()
        # The dummy repetitions load the first acquired line's gradients rather than a zero blip,
        # so what they establish is the steady state the first real repetition then samples.
        for n in range(count):
            out.add(n * tr_s, self.tr(line=table[0], phase_deg=self.phase_deg(n),
                                      acquire=False, center_mm=center_mm))
        for index, line in enumerate(table):
            n = count + index
            out.add(n * tr_s, self.tr(line=line, phase_deg=self.phase_deg(n),
                                      center_mm=center_mm))
        return out

    # -------------------------------------------------------------------- the refusals
    def _resolve_spoil_phase(self, rf_spoil_deg: float | None) -> float:
        """Return the phase increment, refusing one that cannot take effect."""
        if not self.rf_spoil:
            if rf_spoil_deg is not None:
                msg = format_error(
                    f'rf_spoil_deg={rf_spoil_deg:g} was passed with rf_spoil=False, so it '
                    f'cannot take effect.',
                    {'rf_spoil': False, 'rf_spoil_deg': rf_spoil_deg},
                    [
                        'drop rf_spoil_deg to keep RF spoiling off',
                        'or rf_spoil=True to use it',
                        'gradient spoiling is a separate mechanism and is unaffected either way',
                    ],
                )
                raise ConfigurationError(msg)
            return 0.0
        if rf_spoil_deg is None:
            return _DEFAULT_RF_SPOIL_DEG
        return require_positive(rf_spoil_deg, 'rf_spoil_deg')

    def _check(self, lines: Iterable[int]) -> tuple[int, ...]:
        """Return `lines` as a tuple, refusing what cannot be acquired and warning on DC."""
        table = tuple(int(line) for line in lines)
        ny = self.tr.matrix[1]
        if not table:
            msg = format_error(
                'lines is empty, so there is nothing to acquire.',
                {'matrix': self.tr.matrix},
                [f'lines=range({ny}) is the fully sampled table'],
            )
            raise ConfigurationError(msg)

        out_of_range = sorted({line for line in table if not 0 <= line < ny})
        if out_of_range:
            msg = format_error(
                f'{len(out_of_range)} line(s) fall outside 0 ... {ny - 1}.',
                {'offending': out_of_range[:8], 'matrix': self.tr.matrix},
                [
                    'lines are zero-based phase-encode indices, not k-space positions',
                    f'the centre of k-space is line {self.center_line}',
                ],
            )
            raise ConfigurationError(msg)

        counts = Counter(table)
        repeated = sorted(line for line, n in counts.items() if n > 1)
        if repeated:
            msg = format_error(
                f'{len(repeated)} line(s) appear more than once.',
                {'repeated': repeated[:8]},
                [
                    'two readouts writing one k-space address is what the compiler rejects '
                    'downstream, with less to say about it than here',
                    'averaging is a separate acquisition and a separate label',
                ],
            )
            raise ConfigurationError(msg)

        if self.center_line not in counts:
            warnings.warn(
                f'lines omits the centre of k-space (line {self.center_line}), so the image has '
                f'no DC term; that reads as a windowing problem rather than as a sampling one.',
                SeqCraftWarning,
                stacklevel=3,
            )
        return table


def _require_dummies(dummies: int) -> int:
    """
    Return `dummies` as a non-negative int.

    A free function because two call sites need it and both are build-time: ``build`` plays them
    and ``time_to_center_line`` counts them, and an answer that disagreed with the block it
    describes would be worse than either refusing.

    Examples
    --------
    >>> _require_dummies(0)
    0
    >>> _require_dummies(-1)
    Traceback (most recent call last):
        ...
    seqcraft.errors.ConfigurationError: dummies must not be negative, got -1.
    """
    count = int(dummies)
    if count < 0:
        msg = format_error(f'dummies must not be negative, got {count}.', {'dummies': count})
        raise ConfigurationError(msg)
    return count
