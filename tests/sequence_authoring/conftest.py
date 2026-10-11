"""
One built acquisition, shared by every check in this directory.

The `.seq` is generated into pytest's own temporary directory and never committed -- `*.seq` is
git-ignored, and a recorded file would stop testing the build that produced it.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pypulseq as pp
import pytest

import seqcraft as sc

sys.path.insert(0, str(Path(__file__).parent))

import reference_acquisition as ref  # noqa: E402 -- the fixture lives beside this file


@pytest.fixture(scope='session')
def built(tmp_path_factory) -> dict[str, object]:
    """Build once: the written path, what the script claims, and the warnings it raised."""
    path = tmp_path_factory.mktemp('sequence_authoring') / 'reference_acquisition.seq'
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        claims = ref.build(path)
    return {'path': path, 'warnings': list(caught), **claims}


@pytest.fixture(scope='session')
def emitted(built) -> pp.Sequence:
    """The file read back, which is the only thing these tests are allowed to measure."""
    seq = pp.Sequence()
    seq.read(str(built['path']))
    return seq


@pytest.fixture(scope='session')
def trajectory(emitted) -> dict[str, object]:
    """k-space and event times recomputed from the emitted waveforms, not from the authoring code."""
    k_adc, _k_full, t_excitation, _t_refocusing, t_adc = emitted.calculate_kspace()
    return {'k_adc': k_adc, 't_excitation': t_excitation, 't_adc': t_adc}


__all__ = ['built', 'emitted', 'sc', 'trajectory']
