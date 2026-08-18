"""
:class:`IRPrep` -- an inversion pulse and the spoiler that follows it.

``preparation/`` because the pulse's ``use`` is ``'inversion'``.  That is the folder's membership
rule rather than a description, and it is the same field ``rf/`` is read off; MRzero maps it onto
its own ``PulseUsage``, which is what makes the simulation in ``examples/mprage_2d/02`` a check on
the file rather than on an intention.

The spoiler is why this is a module rather than a wrapper
---------------------------------------------------------
An inversion is never perfect.  B1 inhomogeneity and off-resonance leave transverse magnetisation
behind, and without a crusher it survives a whole TI of recovery into the first readout, as a
stripe that *changes with TI* -- which reads as a contrast problem rather than as a spoiling one.
The pulse and its crusher are one design, and holding two events in fixed relationship is the job
a module exists for.

``spoil_cycles_per_voxel`` therefore defaults to ``8.0`` rather than the readout spoiler's ``4.0``:
there is a whole inversion time for anything left over to rephase in, and the recovery that follows
is long enough that a marginal crusher is not marginal by the time it is sampled.

Why the pulse is adiabatic by default
-------------------------------------
B1 insensitivity is the entire reason to use an inversion pulse of any sophistication.  A nominal
180 degrees that delivers 150 inverts nothing like as well -- and the error is *spatially varying*,
so it appears as shading that survives every uniformity correction downstream.  ``'hypsec'`` and
``'wurst'`` go to :func:`pypulseq.make_adiabatic_pulse`; the other three go to the same factories
:class:`~seqcraft.modules.Excitation` uses, for the cases where a short pulse matters more than
robustness.

Adiabatic pulses are long by construction -- the adiabatic condition is what sets the floor -- so
``duration_s`` defaults to ``10e-3``, matching ``make_adiabatic_pulse``'s own default.  That length
is exactly why :meth:`IRPrep.time_to_center` is load-bearing: TI is measured from the inversion's
*effective centre*, which a 10 ms hyperbolic secant puts 5 ms after the start of the block.
Referencing TI to the block start instead is a 5 ms error in the one quantity the sequence exists
to control.
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
from ..spoiler import spoiler

if TYPE_CHECKING:
    from pypulseq.opts import Opts

    from ...design.events import Event

__all__ = ['IRPrep']

#: The five pulse shapes, and the pypulseq factory each names.  The first two are adiabatic and
#: take no flip angle -- the sweep inverts or it does not.  Held by name and resolved on use, as
#: in ``rf/excitation.py``, so a pypulseq missing one fails through a single complete list.
_ADIABATIC: dict[str, str] = {'hypsec': 'hypsec', 'wurst': 'wurst'}
_FACTORIES: dict[str, str] = {
    'sinc': 'make_sinc_pulse',
    'slr': 'make_slr_pulse',
    'block': 'make_block_pulse',
}

#: Arguments this module owns; `pulse_opts` may not reach past it to set one of them.
_RESERVED = frozenset({
    'delay', 'duration', 'flip_angle', 'freq_offset', 'phase_offset', 'pulse_type',
    'return_delay', 'return_gz', 'slice_thickness', 'system', 'use',
})


class IRPrep(Module):
    """
    An inversion pulse and its spoiler, as one block.

    Parameters
    ----------
    opts
        The scanner.
    thickness_mm
        Slab thickness for a selective inversion, or ``None`` for a non-selective one.  **No
        default**, for the same reason as :class:`~seqcraft.modules.Excitation`: there is none
        that is sensible for the module's most consequential property.  ``None`` is the usual
        MPRAGE choice.
    pulse
        ``'hypsec'`` or ``'wurst'`` (adiabatic), or ``'sinc'``, ``'slr'``, ``'block'``.
    duration_s
        Pulse duration.  The default is ``make_adiabatic_pulse``'s own.
    pulse_opts
        Forwarded to the chosen factory -- ``{'beta': 1200.0, 'mu': 6.0}`` for a hyperbolic
        secant, ``{'time_bw_product': 6.0}`` for a sinc.  An unrecognised key raises, naming what
        this `pulse` accepts.
    spoil_cycles_per_voxel, spoil_axis, spoil_voxel_mm
        The crusher.  `spoil_voxel_mm` is the voxel dimension along `spoil_axis` that the cycles
        are counted across; it defaults to `thickness_mm` and is **required** when the inversion
        is non-selective, because there is then no length in this module to count against and
        guessing one would silently under- or over-spoil.
    axis
        Selection axis.  ``None`` means ``'z'`` when selective, and is the only legal value when
        not: with no selection gradient there is no axis to select along.
    tag
        Optional identity, as for any :class:`~seqcraft.Module`.

    Attributes
    ----------
    rf : Event
    gz : Event or None
        The pulse and, when selective, its selection gradient.  There is **no rephaser**: an
        inversion is not rephased, because no transverse magnetisation is meant to survive it --
        the spoiler's job is to make sure of that.

    Examples
    --------
    >>> import pypulseq as pp
    >>> from pypulseq.opts import Opts
    >>> o = Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
    ...          rf_dead_time=100e-6, rf_ringdown_time=30e-6)
    >>> inv = IRPrep(opts=o, thickness_mm=None, spoil_voxel_mm=5.0)
    >>> inv()                                       # the pulse, then the crusher
    LogicBlock(IRPrep, 2 nodes, 11.34 ms)
    >>> round(inv.time_to_center() * 1e3, 3)        # dead time plus half a 10 ms hypsec
    5.101

    The centre is **not** the middle of the block, because the spoiler is after the pulse:

    >>> round(inv().duration / 2 * 1e3, 3)
    5.67

    A slab-selective one adds the selection gradient, and no rephaser:

    >>> sel = IRPrep(opts=o, thickness_mm=8.0)
    >>> sel()
    LogicBlock(IRPrep, 3 nodes, 11.00 ms)

    Notes
    -----
    **`time_to_center` is load-bearing rather than bookkeeping.**  Every ``time_to_*`` in this
    package answers one question -- seconds from the start of the block *that module's* ``build``
    returns -- and they compose by addition only because they all share that origin.  A TI is laid
    out as ``t_shot + inv.time_to_center() + ti_s - gre.time_to_center_line(...)``, and every term
    in it is one of these.

    **`position_mm` stays `0.0` and `axis` stays `None`.**  Two arguments that both become
    inapplicable when ``thickness_mm=None``, and they get different treatments, matching
    ``Excitation``: passing **any** `axis` raises, while `position_mm` refuses only a **non-zero**
    value.  The rule is ``None`` for an argument with no identity value and the identity value
    itself when there is one -- there is no axis meaning "no axis", but zero is exactly the
    position a non-selective pulse delivers, so making it ``None`` would create two spellings of
    one thing and force every caller to branch on selectivity for a call that is valid either way.

    **A non-selective inversion inverts everything in the coil.**  A second slice acquired before
    the first has recovered therefore sees inverted magnetisation, which is fine for a single-slice
    example and is the reason multi-slice IR needs either a slab-selective inversion or a shot
    interval long enough to recover between slices.
    """

    def __init__(
        self,
        *,
        opts: Opts,
        thickness_mm: float | None,
        pulse: str = 'hypsec',
        duration_s: float = 10e-3,
        pulse_opts: dict[str, Any] | None = None,
        spoil_cycles_per_voxel: float = 8.0,
        spoil_axis: str = 'z',
        spoil_voxel_mm: float | None = None,
        axis: str | None = None,
        tag: str | None = None,
    ) -> None:
        super().__init__(opts=opts, tag=tag)
        self.pulse = self._check_pulse(pulse)
        self.duration_s = require_positive(duration_s, 'duration_s')
        self.selective = thickness_mm is not None
        self.thickness_mm = (
            require_positive(thickness_mm, 'thickness_mm') if self.selective else None
        )
        self.axis = self._check_axis(axis)
        self._check_combination()

        kwargs: dict[str, Any] = {
            'duration': self.duration_s,
            'system': opts,
            'use': 'inversion',
            'delay': float(opts.rf_dead_time),
        }
        if self.pulse in _ADIABATIC:
            kwargs['pulse_type'] = _ADIABATIC[self.pulse]
        else:
            # 180 degrees is not a parameter: a pulse that inverts by some other angle is a
            # different preparation with a different name, and the adiabatic pair have no flip
            # angle to give at all -- the sweep inverts or it does not.
            kwargs['flip_angle'] = float(np.pi)
        if self.selective:
            kwargs['slice_thickness'] = self.thickness_mm / 1e3
            kwargs['return_gz'] = True
        kwargs.update(self._check_pulse_opts(pulse_opts))

        factory = getattr(pp, 'make_adiabatic_pulse' if self.pulse in _ADIABATIC
                          else _FACTORIES[self.pulse])
        if self.selective:
            # The third return is a rephaser, and it is dropped on purpose: rephasing an
            # inversion would restore transverse magnetisation this module is about to crush.
            rf, gz, _rephaser = factory(**kwargs)
            if self.axis != 'z':
                gz = derive(gz, channel=self.axis)
            self.rf, self.gz = rf, gz
        else:
            self.rf, self.gz = factory(**kwargs), None

        self.spoil_axis = require_axis(spoil_axis, 'spoil_axis')
        self.spoil_cycles_per_voxel = require_positive(
            spoil_cycles_per_voxel, 'spoil_cycles_per_voxel',
        )
        self.spoil_voxel_mm = self._check_spoil_voxel(spoil_voxel_mm)
        self.spoiler = spoiler(opts, cycles_per_voxel=self.spoil_cycles_per_voxel,
                               voxel_mm=self.spoil_voxel_mm, axis=self.spoil_axis)
        self._pulse_end_s = float(pp.calc_duration(self.rf) if not self.selective
                                  else max(pp.calc_duration(self.rf), pp.calc_duration(self.gz)))

    # ------------------------------------------------------------------ what it knows
    def time_to_center(self) -> float:
        """
        Seconds from the start of this module's block to the RF's effective centre.

        From ``pp.calc_rf_center``, plus the pulse's own delay, which carries the transmit dead
        time.  Identical in form to :meth:`Excitation.time_to_center`, and deliberately so: two
        modules answering the same question with the same origin is what lets a timeline add them.

        Not measurable from the tree -- an effective centre is a property of the waveform.
        """
        return float(self.rf.delay) + float(pp.calc_rf_center(self.rf)[0])

    # ----------------------------------------------------------------------- assembly
    def build(self, *, phase_deg: float = 0.0, position_mm: float = 0.0) -> LogicBlock:
        """
        Return the inversion pulse, its gradient when selective, and the spoiler after it.

        Parameters
        ----------
        phase_deg
            Carrier phase, degrees.  An inversion normally keeps ``0.0``; it is here because
            phase cycling across shots is a real thing to want and costs nothing to allow.
        position_mm
            Slab offset from isocentre along `axis`.  Requires a selective pulse; ``0.0`` is
            valid either way, and is what a non-selective pulse delivers.
        """
        rf = self.rf
        phase_rad = float(np.deg2rad(phase_deg))
        if phase_rad:
            rf = derive(rf, phase_offset=float(rf.phase_offset) + phase_rad)
        if position_mm:
            rf = shift_slice(rf, self._selection_gradient(), position_m=float(position_mm) / 1e3)

        out = LogicBlock().add(0.0, rf)
        if self.selective:
            out.add(0.0, self.gz)
        return out.add(self._pulse_end_s, self.spoiler)

    def _selection_gradient(self) -> Event:
        """Return the selection gradient, or refuse an offset there is nothing to impose."""
        if not self.selective:
            msg = format_error(
                'position_mm needs a selection gradient, and this inversion is non-selective.',
                {'thickness_mm': None, 'pulse': self.pulse},
                [
                    'pass thickness_mm to make the inversion slab-selective',
                    'or leave position_mm at 0.0, which is what a non-selective pulse delivers',
                ],
            )
            raise ConfigurationError(msg)
        return self.gz

    # -------------------------------------------------------------------- the refusals
    def _check_pulse(self, pulse: str) -> str:
        """Return `pulse` having checked it names one of the five shapes."""
        if pulse not in _ADIABATIC and pulse not in _FACTORIES:
            listed = ', '.join(repr(name) for name in (*_ADIABATIC, *_FACTORIES))
            msg = format_error(
                f'pulse must be one of {listed}, got {pulse!r}.',
                {'pulse': pulse},
                ["'hypsec' and 'wurst' are adiabatic, which is why they are the default"],
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
                        'pass thickness_mm to select a slab along axis',
                        'or drop axis: a non-selective inversion has no selection axis',
                        'spoil_axis is a different argument and applies either way',
                    ],
                )
                raise ConfigurationError(msg)
            return None
        return require_axis('z' if axis is None else axis)

    def _check_combination(self) -> None:
        """Refuse the one `pulse`-dependent combination that cannot mean anything."""
        if self.pulse == 'block' and self.selective:
            msg = format_error(
                "pulse='block' cannot be slab-selective.",
                {'pulse': 'block', 'thickness_mm': self.thickness_mm},
                [
                    'pass thickness_mm=None for a hard inversion',
                    "or pulse='hypsec' / 'wurst' / 'sinc' / 'slr' to select a slab",
                    "a rectangular envelope's profile is a sinc with sidelobes, so "
                    'pp.make_block_pulse offers neither slice_thickness nor return_gz',
                ],
            )
            raise ConfigurationError(msg)

    def _check_spoil_voxel(self, spoil_voxel_mm: float | None) -> float:
        """
        Return the length the crusher's cycles are counted across, refusing to invent one.

        `thickness_mm` when the inversion is selective, and no fallback when it is not: a
        non-selective pulse carries no length at all, and the number that matters -- the imaging
        voxel along `spoil_axis` -- belongs to a sequence this module has never heard of.  A
        default would be a guess whose only symptom is faint residual banding that changes with
        TI, which is exactly the failure the crusher exists to prevent.
        """
        if spoil_voxel_mm is not None:
            return require_positive(spoil_voxel_mm, 'spoil_voxel_mm')
        if self.thickness_mm is not None:
            return self.thickness_mm
        msg = format_error(
            'spoil_voxel_mm is needed: a non-selective inversion has no thickness to count '
            'spoiler cycles across.',
            {'thickness_mm': None, 'spoil_axis': self.spoil_axis,
             'spoil_cycles_per_voxel': self.spoil_cycles_per_voxel},
            [
                'pass spoil_voxel_mm = the imaging voxel dimension along spoil_axis, which for '
                "the usual spoil_axis='z' is the slice thickness",
                'or pass thickness_mm for a slab-selective inversion, which then supplies it',
            ],
        )
        raise ConfigurationError(msg)

    def _check_pulse_opts(self, pulse_opts: dict[str, Any] | None) -> dict[str, Any]:
        """Return `pulse_opts` having checked every key against the chosen factory."""
        if not pulse_opts:
            return {}
        import inspect  # noqa: PLC0415  (only reached when an escape hatch is actually used)

        name = 'make_adiabatic_pulse' if self.pulse in _ADIABATIC else _FACTORIES[self.pulse]
        accepted = set(inspect.signature(getattr(pp, name)).parameters) - _RESERVED
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
