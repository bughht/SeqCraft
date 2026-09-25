"""
Physical intent for a repetition-level augmentation: what is wanted, not how to realise it.

These are **values**, not modules.  They take no ``opts``, build nothing, and designing a waveform
is not among the things they can do -- constructing one costs nothing and commits to nothing.

The split they exist to hold:

.. code-block:: text

    intent      says WHAT is wanted:  "no common-mode first moment at the echo, on y"
    the owner   says WHETHER and HOW: "this repetition owns no adjustable y interval"

So an intent validates only what is knowable without a sequence -- that an axis name is a logical
gradient channel, that a venc is positive -- and nothing about which axes a given repetition can
actually carry a moment claim on.  It cannot know that: the same
:class:`FlowCompensation` is honoured by :class:`~seqcraft.modules.GRE3DTR` on ``z`` and refused by
:class:`~seqcraft.modules.GRE2DTR`, whose ``z`` gradient is a slice rephaser that
:class:`~seqcraft.modules.Excitation` realises for itself.  Capability belongs to the repetition,
which is the only thing that knows what it will emit.

**This is a closed set.**  There is no base class, no protocol and no registration point, and a
third augmentation is a SeqCraft-owned type added here with its own keyword -- not something a
caller can supply.
"""

from __future__ import annotations

from dataclasses import dataclass

from .design.validation import require_axis, require_positive
from .errors import ConfigurationError, format_error

__all__ = ['FlowCompensation', 'VelocityEncoding']


@dataclass(frozen=True)
class FlowCompensation:
    """
    Ask that the **common-mode** first gradient moment be zero at the echo.

    Common mode, not "zero in every acquired state".  The two mean the same thing until something
    else constrains the same moment, and then they do not:

    .. code-block:: text

        alone, one acquired state    the common mode IS that state, so this reduces to m1 = 0.
                                     A spin moving at constant velocity along the axis arrives
                                     at the echo with the phase it would have had standing
                                     still -- which is what "flow compensated" means.  It is a
                                     statement about that phase term, not a promise about how
                                     much artefact a given acquisition will show

        with VelocityEncoding        this fixes the MEAN of the two states at zero while the
                                     velocity encoding fixes their separation, so the pair comes
                                     out at +delta/2 and -delta/2.  NEITHER state has m1 = 0,
                                     and neither should: a moving spin is meant to arrive with a
                                     different phase in each, and that difference is the
                                     velocity signal.  What the compensation removes here is the
                                     background phase common to both, which is what makes the
                                     subtraction a clean velocity map

    The symmetric half-delta is **derived**, not asked for.  Nothing computes it; it falls out of
    two independent constraints on the same moment, which is what lets the two compose without a
    precedence rule.

    Parameters
    ----------
    axis
        Logical gradient channel, or several: ``'y'`` or ``('y', 'z')``.  Compensating two axes
        is two independent requirements on two waveforms, not a coupled one, so asking for both
        here means the same as asking for each -- and it is spelled the way `spoil_axis` already
        is on these kernels.  Whether *this* repetition can carry the requirement on a given axis
        is the repetition's answer, not this object's.

    Examples
    --------
    >>> FlowCompensation(axis='y').axes
    ('y',)
    >>> FlowCompensation(axis=('y', 'z')).axes
    ('y', 'z')
    >>> FlowCompensation(axis='ky')
    Traceback (most recent call last):
        ...
    seqcraft.errors.ConfigurationError: axis must be one of 'x', 'y', 'z', got 'ky'.
    """

    axis: str | tuple[str, ...]

    def __post_init__(self) -> None:
        names = (self.axis,) if isinstance(self.axis, str) else tuple(self.axis)
        if not names:
            raise ConfigurationError(format_error(
                'flow compensation needs at least one axis.', {'axis': self.axis},
                ["pass axis='y', or a tuple such as axis=('y', 'z')"],
            ))
        axes = tuple(require_axis(name) for name in names)
        if len(set(axes)) != len(axes):
            raise ConfigurationError(format_error(
                f'flow compensation names an axis twice: {axes}.', {'axis': self.axis},
                ['one requirement per axis -- naming it twice asks for the same thing twice'],
            ))
        object.__setattr__(self, 'axis', self.axis if isinstance(self.axis, str) else axes)

    @property
    def axes(self) -> tuple[str, ...]:
        """The axes this asks for, always as a tuple, whether one was given or several."""
        return (self.axis,) if isinstance(self.axis, str) else tuple(self.axis)


@dataclass(frozen=True)
class VelocityEncoding:
    """
    Ask that the first gradient moment **differ** by a known amount between two acquired states.

    That difference is what makes velocity measurable: the two states give a moving spin phases
    that differ by ``2 pi delta_m1 v``, and their subtraction is a velocity map.  A stationary spin
    sees the same phase in both and subtracts to zero.

    The states are this object's to define and the caller's to schedule:

    .. code-block:: text

        this intent      the venc relation, and that there are two states, called +1 and -1
        the caller       whether both are acquired, and in what order
        the repetition   realises the ONE state it is handed, via build(encoding_state=...)

    Nothing here decides shot ordering, because which ordering is right is a protocol decision:
    state-innermost minimises the time between the two acquisitions that will be subtracted,
    state-outermost minimises the difference in their eddy-current history.

    Parameters
    ----------
    venc_m_s
        The velocity that earns a half-turn of phase difference between the states, in m/s.
        Velocities beyond it alias.  A smaller `venc_m_s` needs a larger first-moment difference,
        and so a longer waveform.
    axis
        Logical gradient channel.  Velocity is encoded along this axis only.

    Examples
    --------
    The relation is gamma-free -- it is a first moment in s/m, not an area:

    >>> venc = VelocityEncoding(venc_m_s=1.5, axis='z')
    >>> round(venc.delta_m1_s_per_m, 6)
    0.333333
    >>> venc.states
    (1, -1)

    Halving the venc doubles the first-moment difference it needs:

    >>> round(VelocityEncoding(venc_m_s=0.75, axis='z').delta_m1_s_per_m, 6)
    0.666667
    """

    venc_m_s: float
    axis: str

    def __post_init__(self) -> None:
        object.__setattr__(self, 'venc_m_s', require_positive(self.venc_m_s, 'venc_m_s'))
        object.__setattr__(self, 'axis', require_axis(self.axis))

    @property
    def delta_m1_s_per_m(self) -> float:
        """
        The first-moment difference this venc asks for, in s/m.

        Read from :meth:`~seqcraft.modules.VelocityEncode.delta_m1_for`, which stays the single
        source of the relation -- the standalone bipolar and this intent must not be able to
        disagree about what a venc means.
        """
        from .modules import VelocityEncode
        return float(VelocityEncode.delta_m1_for(self.venc_m_s))

    @property
    def axes(self) -> tuple[str, ...]:
        """
        The axis this encodes along, as a one-tuple, so routing reads one shape for both intents.

        Always one.  Encoding two directions is a different acquisition -- it needs more than two
        states -- rather than this intent applied twice, so it is not spelled here.
        """
        return (self.axis,)

    @property
    def states(self) -> tuple[int, int]:
        """
        The two encoding states, ``(+1, -1)``.

        Iterate these rather than writing the pair out, so that what a state *is* stays with the
        object that owns the physics:

        .. code-block:: python

            for state in venc.states:
                scan.add(at, tr(line=line, encoding_state=state))
        """
        return (1, -1)
