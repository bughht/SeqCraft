"""
Unit conversion and the two-tier naming rule.

seqcraft uses **two tiers** and converts between them exactly once, at the boundary
where a user-facing parameter enters a module:

===============================  ===============================================
tier                             convention
===============================  ===============================================
recipes / public parameters      researcher-natural units, unit **in the name**:
                                 ``fov_mm``, ``te_ms``, ``flip_deg``,
                                 ``max_grad_mT_m``, ``readout_duration_us``
internal / derived values        strict SI, SI suffix: ``fov_m``, ``te_s``,
                                 ``flip_rad``, ``max_grad_Hz_per_m``
===============================  ===============================================

There is no ``Quantity`` type and no ``pint``.  A unit type would poison numpy interop
for waveform arrays and would still not prevent the actual observed bug class, which is
a *plausible number in the wrong unit* -- ``fov=220`` meaning millimetres reaching code
that reads metres.  That is caught by range validation in
:mod:`seqcraft.core.validate`, not by a type.

Gradient and slew conversions depend on the gyromagnetic ratio, so they take `gamma`
explicitly rather than reading a global.  ``System`` supplies it.

Examples
--------
>>> mm(220)
0.22
>>> round(mT_per_m(40), 3)          # 40 mT/m at the proton gamma
1703040.0
>>> round(deg(90), 6)
1.570796
>>> s_per_mm2(1000)                 # b-value: s/mm^2 -> s/m^2
1000000000.0
"""

from __future__ import annotations

import math

__all__ = [
    'GAMMA_1H',
    'Hz_per_m_to_mT_per_m',
    'T_per_m_per_s',
    'cm',
    'deg',
    'kHz',
    'mT_per_m',
    'm_per_s2_to_Hz',
    'mm',
    'ms',
    'ppm',
    'rad_to_deg',
    's_per_mm2',
    's_per_m2_to_s_per_mm2',
    'uT',
    'um',
    'us',
]

#: Proton gyromagnetic ratio in Hz/T, matching ``pypulseq.Opts`` default ``gamma``.
GAMMA_1H = 42_576_000.0


# --------------------------------------------------------------------------- length
def mm(x: float) -> float:
    """Millimetres to metres."""
    return x / 1e3


def cm(x: float) -> float:
    """Centimetres to metres."""
    return x / 1e2


def um(x: float) -> float:
    """Micrometres to metres."""
    return x / 1e6


# ----------------------------------------------------------------------------- time
def ms(x: float) -> float:
    """Milliseconds to seconds."""
    return x / 1e3


def us(x: float) -> float:
    """Microseconds to seconds."""
    return x / 1e6


# ------------------------------------------------------------------------ frequency
def kHz(x: float) -> float:
    """Kilohertz to hertz."""
    return x * 1e3


def ppm(x: float) -> float:
    """Parts per million as a dimensionless fraction."""
    return x / 1e6


# ---------------------------------------------------------------------------- angle
def deg(x: float) -> float:
    """Degrees to radians."""
    return x * math.pi / 180.0


def rad_to_deg(x: float) -> float:
    """Radians to degrees."""
    return x * 180.0 / math.pi


# ------------------------------------------------------------------------- gradient
def mT_per_m(x: float, gamma: float = GAMMA_1H) -> float:
    """
    Gradient amplitude in mT/m to pulseq units (Hz/m).

    Examples
    --------
    >>> round(mT_per_m(80), 1)
    3406080.0
    """
    return x / 1e3 * gamma


def Hz_per_m_to_mT_per_m(x: float, gamma: float = GAMMA_1H) -> float:
    """Gradient amplitude in Hz/m back to mT/m, for reports and error messages."""
    return x / gamma * 1e3


def T_per_m_per_s(x: float, gamma: float = GAMMA_1H) -> float:
    """
    Slew rate in T/m/s to pulseq units (Hz/m/s).

    Examples
    --------
    >>> round(T_per_m_per_s(200), 1)
    8515200000.0
    """
    return x * gamma


def Hz_per_m_per_s_to_T_per_m_per_s(x: float, gamma: float = GAMMA_1H) -> float:
    """Slew rate in Hz/m/s back to T/m/s, for reports and error messages."""
    return x / gamma


def uT(x: float, gamma: float = GAMMA_1H) -> float:
    """B1 amplitude in microtesla to hertz."""
    return x / 1e6 * gamma


# -------------------------------------------------------------------------- b-value
def s_per_mm2(x: float) -> float:
    """
    b-value in s/mm^2 to SI (s/m^2).

    b-values are quoted in s/mm^2 universally, so this conversion appears at every
    diffusion call site.  Internally b is always s/m^2 and is printed back in s/mm^2.

    Examples
    --------
    >>> s_per_mm2(1000)
    1000000000.0
    """
    return x * 1e6


def s_per_m2_to_s_per_mm2(x: float) -> float:
    """b-value in SI (s/m^2) back to s/mm^2, for reports."""
    return x / 1e6


# ------------------------------------------------------------------------ mechanics
def m_per_s2_to_Hz(x: float, gamma: float = GAMMA_1H) -> float:
    """Convenience for first-moment work: (m/s^2) to Hz via gamma."""
    return x * gamma
