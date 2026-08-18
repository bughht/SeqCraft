"""
:class:`Excitation` -- an RF pulse and, when selective, its gradient and rephaser.

``rf/`` because the pulse's ``use`` is ``'excitation'``.  That is the folder's membership rule
rather than a description: a saturation or preparation pulse is a different ``use`` and belongs
somewhere else, and MRzero maps the same field onto its ``PulseUsage``.

What this module is, honestly
-----------------------------
**A coordinator, not a calculator.**  It was written expecting to compute its own rephaser,
because the classic error is to rewind half the selection gradient's *total* area rather than the
area after the pulse's effective centre -- which is the same number only for a symmetric pulse.
That expectation was tested rather than assumed, and this pypulseq computes it correctly for all
three shaped factories: ``make_sinc_pulse`` and ``make_gauss_pulse`` reference ``center_pos``,
and ``make_slr_pulse`` references ``calc_rf_center`` of the designed waveform, so a minimum-phase
SLR pulse whose peak sits three-quarters of the way through gets the right rephaser too.

So this module does not recompute it, and ``tests/modules/test_excitation.py`` asserts the
invariant against an asymmetric pulse rather than trusting the reading.  What is left for the
module is real but smaller: choosing among four factories with one vocabulary, refusing the
combinations that cannot mean anything, and holding the RF, the gradient and the rephaser in
fixed relationship so a caller places one block instead of three events.

One angular unit
----------------
``flip_deg``, ``phase_deg``, ``rf_spoil_deg`` -- every angle in this library is in degrees,
because that is how each is quoted in a protocol (15 degrees, 117 degrees).  Radians appear
exactly once, in the call to pypulseq, so no notebook carries a ``np.deg2rad`` and no signature
mixes conventions.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
import pypulseq as pp

from ...design.events import derive
from ...design.logic import LogicBlock
from ...design.module import Module
from ...errors import ConfigurationError, format_error
from .._support import require_axis, require_positive, shift_slice

if TYPE_CHECKING:
    from pypulseq.opts import Opts

    from ...design.events import Event

__all__ = ['Excitation']

#: The four pulse shapes, and the pypulseq factory each names.  Held by name and resolved on
#: use, so a pypulseq missing one fails through :func:`seqcraft._compat.require`'s single
#: complete list rather than through an import-time ``AttributeError`` in this file.
_FACTORIES: dict[str, str] = {
    'sinc': 'make_sinc_pulse',
    'slr': 'make_slr_pulse',
    'gauss': 'make_gauss_pulse',
    'block': 'make_block_pulse',
}

#: Arguments this module owns; `pulse_opts` may not reach past it to set one of them.
_RESERVED = frozenset({
    'delay', 'duration', 'flip_angle', 'freq_offset', 'phase_offset', 'return_delay',
    'return_gz', 'slice_thickness', 'system', 'time_bw_product', 'use',
})


class Excitation(Module):
    """
    An excitation pulse, with its selection gradient and rephaser when it is selective.

    Parameters
    ----------
    opts
        The scanner.
    flip_deg
        Flip angle, degrees.
    thickness_mm
        Slice or slab thickness, millimetres, or ``None`` for a non-selective pulse.  **No
        default**: there is none that is sensible for the module's most consequential property,
        so ``5.0`` or ``None`` is written out at every call site.
    duration_s
        Pulse duration.
    time_bw_product
        ``None`` defers to the factory's own default, which differs by shape -- 4.0 for a sinc,
        3.0 for a gauss.  Passing it with ``pulse='block'`` raises.
    pulse
        One of ``'sinc'``, ``'slr'``, ``'gauss'``, ``'block'``.  An unknown value raises, listing
        them.
    pulse_opts
        Forwarded to the chosen factory -- ``{'filter_type': 'pm', 'passband_ripple': 0.005}``.
        Tuning a stopband ripple is pulse *design*, so it lives behind an escape hatch rather
        than in the signature.  An unrecognised key raises, naming what this `pulse` accepts.
    axis
        Selection axis.  ``None`` means ``'z'`` when selective, and is the only legal value when
        not: with no selection gradient there is no axis to select along.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
    ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6)
    >>> exc = Excitation(opts=o, flip_deg=15.0, thickness_mm=5.0)
    >>> exc()                                   # rf, gz, gzr
    LogicBlock(Excitation, 3 nodes, 3.67 ms)
    >>> round(exc.time_to_center() * 1e6)       # rf_dead_time + half of a 3 ms sinc
    1600

    A hard pulse is the RF event alone:

    >>> hard = Excitation(opts=o, flip_deg=90.0, thickness_mm=None, pulse='block',
    ...                   duration_s=500e-6)
    >>> hard()
    LogicBlock(Excitation, 1 node, 0.63 ms)

    Notes
    -----
    **`thickness_mm=None` means no selection gradient and no rephaser.**  That is a hard pulse
    for ``'block'`` and a *spectrally* selective one for the shaped pulses -- a sinc or a gauss
    played with no gradient is how water excitation works.  ``pp.make_block_pulse`` has neither
    ``slice_thickness`` nor ``return_gz``, because a rectangular envelope's slice profile is a
    sinc with sidelobes nobody wants, so ``pulse='block'`` with a thickness raises and says so
    rather than surfacing an unexpected keyword.

    With ``thickness_mm=None``, passing `axis` raises and `position_mm` must be ``0.0``:
    offsetting a slice biases the RF against the selection gradient's amplitude, and with no
    gradient there is nothing to bias against.
    """

    def __init__(
        self,
        *,
        opts: Opts,
        flip_deg: float,
        thickness_mm: float | None,
        duration_s: float = 3e-3,
        time_bw_product: float | None = None,
        pulse: str = 'sinc',
        pulse_opts: dict[str, Any] | None = None,
        axis: str | None = None,
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.pulse = self._check_pulse(pulse)
        self.flip_deg = require_positive(flip_deg, 'flip_deg')
        self.duration_s = require_positive(duration_s, 'duration_s')
        self.selective = thickness_mm is not None
        self.thickness_mm = (
            require_positive(thickness_mm, 'thickness_mm') if self.selective else None
        )
        self.axis = self._check_axis(axis)
        self._check_combination(time_bw_product)

        kwargs: dict[str, Any] = {
            'flip_angle': np.deg2rad(self.flip_deg),
            'duration': self.duration_s,
            'system': opts,
            'use': 'excitation',
            'delay': float(opts.rf_dead_time),
        }
        if time_bw_product is not None:
            kwargs['time_bw_product'] = require_positive(time_bw_product, 'time_bw_product')
        if self.selective:
            kwargs['slice_thickness'] = self.thickness_mm / 1e3
            kwargs['return_gz'] = True
        kwargs.update(self._check_pulse_opts(pulse_opts))

        factory = getattr(pp, _FACTORIES[self.pulse])
        if self.selective:
            rf, gz, gzr = factory(**kwargs)
            if self.axis != 'z':
                gz, gzr = derive(gz, channel=self.axis), derive(gzr, channel=self.axis)
            self.rf, self.gz, self.gzr = rf, gz, gzr
        else:
            self.rf = factory(**kwargs)
            self.gz = self.gzr = None

    # ------------------------------------------------------------------ what it knows
    def time_to_center(self) -> float:
        """
        Seconds from the start of this module's block to the RF's effective centre.

        From ``pp.calc_rf_center``, plus the pulse's own delay, which carries the transmit dead
        time.  Spelled ``center`` rather than ``centre`` because it sits beside the pypulseq
        function it reads: identifiers follow the API they neighbour, prose follows the package.

        Not measurable from the tree -- an effective centre is a property of the waveform, and
        for a minimum-phase pulse it is nowhere near the midpoint.
        """
        return float(self.rf.delay) + float(pp.calc_rf_center(self.rf)[0])

    def time_to_rephaser(self) -> float:
        """
        Seconds from the start of this module's block to where the slice rephaser begins.

        **This is the earliest another axis may start doing work.**  The rephaser is on ``z`` and
        occupies the tail of this module's block, so a caller that waits for the whole block
        before starting its own gradients waits through it for nothing -- an encode on ``y`` and
        a prephaser on ``x`` play perfectly happily beside it, and every microsecond spent
        waiting is a microsecond added to TE.

        For a non-selective pulse there is no rephaser and the answer is the end of the block.

        See :attr:`rephaser_duration_s`, which is the other half of what a caller needs: the
        rephaser is the third participant in the winder coupling, not a phase before it.
        """
        if not self.selective:
            return float(pp.calc_duration(self.rf))
        return float(pp.calc_duration(self.gz))

    @property
    def rephaser_duration_s(self) -> float:
        """
        Seconds the slice rephaser occupies, or ``0.0`` when the pulse is non-selective.

        A caller overlapping its own gradients with the rephaser has to be at least this long
        before it starts anything that must follow the rephasing -- which is why this is reported
        rather than left to be measured off a block whose other events run past it.
        """
        if not self.selective:
            return 0.0
        return float(pp.calc_duration(self.gzr))

    # ----------------------------------------------------------------------- assembly
    def build(self, *, phase_deg: float = 0.0, position_mm: float = 0.0) -> LogicBlock:
        """
        Return the pulse, and its gradient pair when selective.

        Parameters
        ----------
        phase_deg
            Carrier phase for this repetition, degrees.  This is where an RF-spoiling schedule
            lands; the schedule itself belongs to whatever knows how many repetitions there are.
        position_mm
            Slice offset from isocentre along `axis`.  Requires a selective pulse.
        """
        rf = self.rf
        phase_rad = float(np.deg2rad(phase_deg))
        if phase_rad:
            rf = derive(rf, phase_offset=float(rf.phase_offset) + phase_rad)
        if position_mm:
            rf = shift_slice(rf, self._selection_gradient(), position_m=float(position_mm) / 1e3)

        out = LogicBlock().add(0.0, rf)
        if self.selective:
            out.add(0.0, self.gz).add(float(pp.calc_duration(self.gz)), self.gzr)
        return out

    def _selection_gradient(self) -> Event:
        """Return the selection gradient, or refuse an offset there is nothing to impose."""
        if not self.selective:
            msg = format_error(
                'position_mm needs a selection gradient, and this pulse is non-selective.',
                {'thickness_mm': None, 'pulse': self.pulse},
                [
                    'pass thickness_mm to make the pulse slice-selective',
                    'or leave position_mm at 0.0',
                    'a frequency-selective offset is a different thing: make_sinc_pulse takes '
                    'freq_ppm, which computes hertz from B0',
                ],
            )
            raise ConfigurationError(msg)
        return self.gz

    # -------------------------------------------------------------------- the refusals
    def _check_pulse(self, pulse: str) -> str:
        """Return `pulse` having checked it names one of the four factories."""
        if pulse not in _FACTORIES:
            listed = ', '.join(repr(name) for name in _FACTORIES)
            msg = format_error(
                f'pulse must be one of {listed}, got {pulse!r}.', {'pulse': pulse},
            )
            raise ConfigurationError(msg)
        return pulse

    def _check_axis(self, axis: str | None) -> str | None:
        """Return the selection axis: ``'z'`` by default when selective, ``None`` when not."""
        if not self.selective:
            if axis is not None:
                msg = format_error(
                    f'axis={axis!r} was passed with thickness_mm=None, so it cannot take effect.',
                    {'axis': axis, 'thickness_mm': None},
                    [
                        'pass thickness_mm to select a slice along axis',
                        'or drop axis: a non-selective pulse has no selection axis',
                    ],
                )
                raise ConfigurationError(msg)
            return None
        return require_axis('z' if axis is None else axis)

    def _check_combination(self, time_bw_product: float | None) -> None:
        """Refuse the two `pulse`-dependent combinations that cannot mean anything."""
        if self.pulse != 'block':
            return
        if self.selective:
            msg = format_error(
                "pulse='block' cannot be slice-selective.",
                {'pulse': 'block', 'thickness_mm': self.thickness_mm},
                [
                    "pass thickness_mm=None for a hard pulse",
                    "or pulse='sinc' / 'slr' / 'gauss' to select a slice",
                    "a rectangular envelope's slice profile is a sinc with sidelobes, so "
                    'pp.make_block_pulse offers neither slice_thickness nor return_gz',
                ],
            )
            raise ConfigurationError(msg)
        if time_bw_product is not None:
            msg = format_error(
                "time_bw_product has no meaning for pulse='block' at a fixed duration.",
                {'pulse': 'block', 'time_bw_product': time_bw_product},
                ['drop time_bw_product', "or pulse='sinc', whose bandwidth it sets"],
            )
            raise ConfigurationError(msg)

    def _check_pulse_opts(self, pulse_opts: dict[str, Any] | None) -> dict[str, Any]:
        """Return `pulse_opts` having checked every key against the chosen factory."""
        if not pulse_opts:
            return {}
        import inspect  # noqa: PLC0415  (only reached when an escape hatch is actually used)

        factory = getattr(pp, _FACTORIES[self.pulse])
        accepted = set(inspect.signature(factory).parameters) - _RESERVED
        unknown = sorted(set(pulse_opts) - accepted)
        if unknown:
            msg = format_error(
                f'pulse_opts key(s) {", ".join(repr(k) for k in unknown)} are not accepted by '
                f'pulse={self.pulse!r}.',
                {'pulse': self.pulse, 'accepted': ', '.join(sorted(accepted))},
                ['pulse_opts forwards pulse-design parameters only; this module owns the rest'],
            )
            raise ConfigurationError(msg)
        return dict(pulse_opts)
