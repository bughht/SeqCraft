"""
:class:`SaturationPrep` -- a spectrally selective pulse and the spoiler that destroys what it made.

``preparation/`` because the pulse's ``use`` is ``'saturation'``.  That is the folder's membership
rule, and this module is the sibling :class:`~seqcraft.modules.IRPrep` was always going to have:
``modules/__init__.py`` defines the folder as ``rf.use in {inversion, saturation, preparation}``
and only the first seat was taken.

Why this is not :class:`~seqcraft.modules.Excitation` with a different label
---------------------------------------------------------------------------
The two operations are opposites.  An excitation creates transverse magnetisation *to be read*,
and its selective form owns selection and rephasing so the signal survives to the echo.  A
saturation creates transverse magnetisation *to be thrown away*: nothing is rephased, a spoiler
follows immediately, and the point of the pulse is what is left **longitudinally** afterwards --
which is as close to nothing as the flip angle and the B1 field allow.

So ``Excitation(use='saturation')`` is not the missing API.  Adding a semantic-purpose flag to
``Excitation`` would turn a module with one physical contract into a generic RF pulse with a label,
and the label is the part that would then be carrying the meaning.  Reuse of waveform-design
machinery is not the same thing as reuse of the public physical abstraction.

The spectral offset is the physical contract, and its sign is load-bearing
--------------------------------------------------------------------------
`shift_ppm` is **signed, relative to water**, and fat is below water: about ``-3.4``.  The module
converts once,

.. math:: \\Delta f = \\text{shift\\_ppm} \\times 10^{-6} \\times B_0 \\times \\gamma

and emits that number as the pulse's ``freq_offset``.  Nothing else in the module depends on the
sign, which is exactly why it is dangerous: a sign error produces legal Pulseq, legal timing,
legal gradients and a correct-looking waveform, and saturates the wrong side of the spectrum.
The reference implementations reach the same number by different routes -- one carries the sign in
the ppm constant, the other applies it at the point of use -- and an implementation that mixed the
two conventions would saturate water and leave fat alone.

:attr:`SaturationPrep.offset_hz` reports the conversion so a caller and a test can check the chain
end to end, and ``tests/modules/test_saturation_prep.py`` traces it as far as the compiled
sequence rather than stopping at the constructor.

**Where :math:`B_0` comes from, and the trap in it.**  From ``opts.B0``, and nowhere else --
there is no field strength written into this module.  The chain is
``shift_ppm`` (portable, scanner-independent) → ``opts.B0`` (explicit, the scanner's) →
``freq_offset`` in hertz (what the scanner is told), and :attr:`SaturationPrep.b0_t` reports the
middle term so it is visible rather than implied.

That matters because **pypulseq's ``Opts`` defaults ``B0`` to 1.5 T when it is not given**, and
says nothing.  An ``Opts`` built for a 2.89 T system without passing ``B0=2.89`` puts fat at
-220 Hz instead of -424 Hz -- a legal sequence, a plausible number, and the wrong one.  This
module cannot tell an omitted ``B0`` from a deliberate 1.5 T one, so it reports what it used and
``tests/modules/test_saturation_prep.py`` pins the behaviour rather than leaving it to be
discovered.

**``freq_ppm`` was considered and rejected for the emitted event.**  pypulseq can carry a ppm
offset and let the interpreter resolve it against the scanner's own :math:`B_0`, which is in
principle more portable.  It also moves the conversion off the file, so the sign this module is
most likely to get wrong would no longer be visible in what we emit.  A number we can inspect is
worth more here than one the scanner computes.

What it does not do
-------------------
No selection gradient, and therefore no rephaser and no ``thickness_mm``: a shaped pulse played
with no gradient is spectrally selective, which is the whole mechanism.  A *spatial* saturation
slab plays a gradient and selects a region, and it is deliberately not folded in -- what it shares
with this is "prepare, then spoil", which is a waveform silhouette rather than a physical solve.
CEST saturation trains, with their duty-cycle and B1rms constraints, are a different problem again.

Protocol values are the caller's
--------------------------------
`flip_deg`, `duration_s` and `bandwidth_hz` have **no defaults**.  Two independent references
disagree on all three -- 110 degrees against 90, 8 ms against 12, a bandwidth derived from the
offset against a fixed 200 Hz -- and both work.  That is what a protocol parameter looks like, and
picking one lab's number to use as a default would encode a protocol as physics.

Examples
--------
>>> import pypulseq as pp
>>> import seqcraft as sc
>>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6, B0=2.89)
>>> fat = sc.modules.SaturationPrep(opts=opts, shift_ppm=-3.45, flip_deg=110.0,
...                                 duration_s=8e-3, bandwidth_hz=424.5, spoil_voxel_mm=0.1)
>>> round(fat.offset_hz, 2)
-424.5
>>> fat.rf.use
'saturation'
"""

from __future__ import annotations

import inspect
from typing import TYPE_CHECKING, Any

import numpy as np
import pypulseq as pp

from ...design.events import derive
from ...design.logic import LogicBlock
from ...design.module import Module
from ...errors import ConfigurationError, format_error
from .._support import require_axis, require_positive
from ..spoiler import spoiler

if TYPE_CHECKING:
    from pypulseq.opts import Opts

__all__ = ['SaturationPrep']

#: The shaped factories a spectrally selective pulse can come from.  ``'block'`` is absent: a
#: rectangular envelope's spectral profile is a sinc whose sidelobes reach water, which is the one
#: thing this module exists to avoid.  Held by name and resolved on use, as in ``rf/excitation.py``.
_FACTORIES: dict[str, str] = {
    'gauss': 'make_gauss_pulse',
    'sinc': 'make_sinc_pulse',
    'slr': 'make_slr_pulse',
}

#: Arguments this module owns; `pulse_opts` may not reach past it to set one of them.
_RESERVED = frozenset({
    'bandwidth', 'delay', 'duration', 'flip_angle', 'freq_offset', 'freq_ppm', 'phase_offset',
    'return_delay', 'return_gz', 'slice_thickness', 'system', 'time_bw_product', 'use',
})


class SaturationPrep(Module):
    """
    A spectrally selective saturation pulse and its spoiler, as one block.

    Parameters
    ----------
    opts
        The scanner.  ``opts.B0`` and ``opts.gamma`` set the ppm-to-hertz conversion, so a
        protocol written in ppm moves between field strengths without being rewritten.
    shift_ppm
        Chemical shift of the species to saturate, **signed, relative to water**.  Fat is below
        water: about ``-3.4``.  No default -- which species is being suppressed is the one thing
        the caller must say.
    flip_deg, duration_s, bandwidth_hz
        The pulse.  **No defaults**; see the module docstring.  `bandwidth_hz` is the excited
        band's full width, and the module refuses a combination whose band reaches water.
    pulse
        ``'gauss'`` (both references), ``'sinc'`` or ``'slr'``.  ``'block'`` is not offered.
    pulse_opts
        Forwarded to the chosen factory -- ``{'apodization': 0.42}`` for a gauss.  An unrecognised
        key raises, naming what this `pulse` accepts.
    spoil_cycles_per_voxel, spoil_axis, spoil_voxel_mm
        The spoiler.  `spoil_voxel_mm` is **required**: there is no length anywhere in this module
        to count cycles against, and guessing one would silently under- or over-spoil.  What
        matters physically is the product, so a caller wanting the references' "spoil to 0.1 mm"
        passes ``spoil_voxel_mm=0.1`` and leaves the cycles at 1.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    rf : Event
        The pulse, with ``use='saturation'`` and ``freq_offset`` in hertz.
    spoiler : LogicBlock
    offset_hz : float
    b0_t : float
        The field strength used for the conversion, read from `opts`.  Reported because
        ``pp.Opts`` defaults it to 1.5 T in silence, and a wrong one is invisible downstream.
    """

    def __init__(
        self,
        *,
        opts: Opts,
        shift_ppm: float,
        flip_deg: float,
        duration_s: float,
        bandwidth_hz: float,
        spoil_voxel_mm: float,
        pulse: str = 'gauss',
        pulse_opts: dict[str, Any] | None = None,
        spoil_cycles_per_voxel: float = 1.0,
        spoil_axis: str = 'z',
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.pulse = self._check_pulse(pulse)
        self.shift_ppm = self._check_shift(shift_ppm)
        self.flip_deg = require_positive(flip_deg, 'flip_deg')
        self.duration_s = require_positive(duration_s, 'duration_s')
        self.bandwidth_hz = require_positive(bandwidth_hz, 'bandwidth_hz')

        # The one conversion, done once, in one place.  Everything downstream reads `offset_hz`.
        # `b0_t` is kept rather than inlined so the field strength the conversion used is a fact
        # the caller can read back -- see the module docstring on Opts' silent 1.5 T default.
        self.b0_t = float(opts.B0)
        self.offset_hz = float(self.shift_ppm) * 1e-6 * self.b0_t * float(opts.gamma)
        self._check_clears_water()

        kwargs: dict[str, Any] = {
            'flip_angle': float(np.deg2rad(self.flip_deg)),
            'duration': self.duration_s,
            'freq_offset': self.offset_hz,
            'system': opts,
            'use': 'saturation',
            'delay': float(opts.rf_dead_time),
        }
        # One bandwidth, two ways of saying it.  `make_gauss_pulse` takes hertz directly; the
        # sinc and SLR factories take a time-bandwidth product, which is the same statement
        # multiplied by the duration -- and is how R3 expresses it, `timeBwProduct = BW * PW`.
        # The module owns the conversion so the caller states bandwidth once, in hertz.
        accepted = set(inspect.signature(getattr(pp, _FACTORIES[self.pulse])).parameters)
        if 'bandwidth' in accepted:
            kwargs['bandwidth'] = self.bandwidth_hz
        else:
            kwargs['time_bw_product'] = self.bandwidth_hz * self.duration_s
        kwargs.update(self._check_pulse_opts(pulse_opts))
        # No `return_gz`, so no selection gradient and no rephaser to drop.  The absence is the
        # contract, not an omission.
        self.rf = getattr(pp, _FACTORIES[self.pulse])(**kwargs)

        self.spoil_axis = require_axis(spoil_axis, 'spoil_axis')
        self.spoil_cycles_per_voxel = require_positive(
            spoil_cycles_per_voxel, 'spoil_cycles_per_voxel',
        )
        self.spoil_voxel_mm = require_positive(spoil_voxel_mm, 'spoil_voxel_mm')
        self.spoiler = spoiler(opts, cycles_per_voxel=self.spoil_cycles_per_voxel,
                               voxel_mm=self.spoil_voxel_mm, axis=self.spoil_axis)
        self._pulse_end_s = float(pp.calc_duration(self.rf))

    # ------------------------------------------------------------------ what it knows
    def time_to_center(self) -> float:
        """
        Seconds from the start of this module's block to the RF's effective centre.

        Same form and same origin as :meth:`IRPrep.time_to_center` and
        :meth:`Excitation.time_to_center`, deliberately: three modules answering one question the
        same way is what lets a caller add them along a timeline.

        Not measurable from the tree -- an effective centre is a property of the waveform.
        """
        return float(self.rf.delay) + float(pp.calc_rf_center(self.rf)[0])

    @property
    def band_edge_hz(self) -> float:
        """
        Hertz from water to the near edge of the excited band; negative means the band reaches it.

        ``abs(offset_hz) - bandwidth_hz / 2``.  This is the quantity
        :meth:`_check_clears_water` refuses on, exposed because a caller choosing a duration is
        choosing this number whether or not they know it.
        """
        return abs(self.offset_hz) - self.bandwidth_hz / 2.0

    # ----------------------------------------------------------------------- assembly
    def build(self, *, phase_deg: float = 0.0) -> LogicBlock:
        """
        Return the saturation pulse followed by its spoiler.

        Parameters
        ----------
        phase_deg
            Carrier phase, degrees.  Present because phase cycling a preparation across shots is a
            real thing to want; a single saturation keeps ``0.0``.

        Notes
        -----
        **The offset-induced phase is not compensated.**  One reference applies
        ``-2*pi*freq_offset*rf_center`` and the other does not, neither says why, and the
        transverse magnetisation this pulse makes is spoiled a few hundred microseconds later --
        so there is nothing left for the phase to matter to.  Recorded rather than silently
        chosen; if a use appears where it matters, that is evidence for a parameter, and until
        then a flag would be speculation.
        """
        rf = self.rf
        phase_rad = float(np.deg2rad(phase_deg))
        if phase_rad:
            rf = derive(rf, phase_offset=float(rf.phase_offset) + phase_rad)
        return LogicBlock().add(0.0, rf).add(self._pulse_end_s, self.spoiler)

    # -------------------------------------------------------------------- the refusals
    def _check_pulse(self, pulse: str) -> str:
        """Return `pulse` having checked it names a shaped factory."""
        if pulse not in _FACTORIES:
            listed = ', '.join(repr(name) for name in _FACTORIES)
            fixes = [f'pick one of {listed}']
            if pulse == 'block':
                fixes.append(
                    "'block' is not offered: a rectangular envelope's spectral profile is a sinc "
                    'whose sidelobes reach water, which is what this module exists to avoid'
                )
            msg = format_error(
                f'pulse must be one of {listed}, got {pulse!r}.', {'pulse': pulse}, fixes,
            )
            raise ConfigurationError(msg)
        return pulse

    def _check_shift(self, shift_ppm: float) -> float:
        """Return `shift_ppm` having refused the one value that cannot mean anything."""
        if shift_ppm == 0.0:
            msg = format_error(
                'shift_ppm=0.0 is water itself, so this pulse would saturate the signal.',
                {'shift_ppm': shift_ppm},
                [
                    'fat is about -3.4 ppm, signed and relative to water',
                    'a deliberate on-resonance preparation is a different module',
                ],
            )
            raise ConfigurationError(msg)
        return float(shift_ppm)

    def _check_clears_water(self) -> None:
        """Refuse a pulse whose excited band reaches the water line."""
        if self.band_edge_hz > 0.0:
            return
        msg = format_error(
            'the excited band reaches water, so this pulse would saturate the signal it is '
            'meant to preserve.',
            {
                'shift_ppm': self.shift_ppm,
                'b0_t': self.b0_t,
                'offset_hz': round(self.offset_hz, 2),
                'bandwidth_hz': self.bandwidth_hz,
                'band_edge_hz': round(self.band_edge_hz, 2),
            },
            [
                f'narrow the pulse: bandwidth_hz below {2 * abs(self.offset_hz):.1f} Hz',
                'or lengthen duration_s, which is the same statement in the time domain',
                'a larger |shift_ppm|, or a higher B0, moves the species further from water',
            ],
        )
        raise ConfigurationError(msg)

    def _check_pulse_opts(self, pulse_opts: dict[str, Any] | None) -> dict[str, Any]:
        """Return `pulse_opts` having checked every key against the chosen factory."""
        if not pulse_opts:
            return {}
        accepted = set(inspect.signature(getattr(pp, _FACTORIES[self.pulse])).parameters)
        accepted -= _RESERVED
        unknown = sorted(set(pulse_opts) - accepted)
        if unknown:
            named = ', '.join(repr(key) for key in unknown)
            msg = format_error(
                f'pulse_opts key(s) {named} are not accepted by pulse={self.pulse!r}.',
                {'pulse': self.pulse, 'accepted': ', '.join(sorted(accepted))},
                ['pulse_opts forwards pulse-design parameters only; this module owns the rest'],
            )
            raise ConfigurationError(msg)
        return dict(pulse_opts)
