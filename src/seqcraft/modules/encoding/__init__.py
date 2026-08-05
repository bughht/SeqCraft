"""Encoding modules: phase/partition encodes, spoilers, diffusion."""

from __future__ import annotations

from .cartesian import Crusher, PartitionEncode, PhaseEncode, Prephaser, Spoiler
from .diffusion import (
    ArbitraryDiffusion,
    BipolarDiffusion,
    MonopolarDiffusion,
    direction_condition_number,
    dti_directions,
)

__all__ = [
    'ArbitraryDiffusion',
    'BipolarDiffusion',
    'Crusher',
    'MonopolarDiffusion',
    'PartitionEncode',
    'PhaseEncode',
    'Prephaser',
    'Spoiler',
    'direction_condition_number',
    'dti_directions',
]
