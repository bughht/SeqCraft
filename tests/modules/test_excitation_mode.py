"""
The selective/non-selective excitation rule, where it now lives.

**Selective and non-selective are two physical excitation modes, not one mode with the selection
gradient switched off.**  Each has its own RF family: a slab wants a shaped pulse, no slab wants a
hard one.  That rule was established in ``GRE3DTR`` and is now in one place so a second 3D kernel
inherits it rather than re-deriving it.

What these tests pin is the *emitted pair* -- which RF, and which gradient -- because that pair is
what the mode is.  Asking only "is Gz present?" cannot tell a non-selective excitation from a
selective one whose gradient was deleted, and those are not the same design.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc
from seqcraft.errors import ConfigurationError
from seqcraft.modules import GRE3DTR
from seqcraft.modules._support import HARD_PULSE_S, resolve_excitation_mode

_SELECTION = ('slab_thickness_mm', 140.0)
_NO_SELECTION = ('slab_thickness_mm', None)


def _resolve(selective: bool, **kw) -> dict:
    kw.setdefault('rf_pulse', None)
    kw.setdefault('rf_duration_s', None)
    kw.setdefault('rf_time_bw_product', None)
    return resolve_excitation_mode(
        selective=selective, selection=_SELECTION if selective else _NO_SELECTION, **kw
    )


# ------------------------------------------------------------------ the mode chooses the family
def test_a_slab_defaults_to_a_shaped_pulse() -> None:
    assert _resolve(True)['pulse'] == 'sinc'


def test_no_slab_defaults_to_a_hard_pulse() -> None:
    """The point of the rule: the non-selective mode is not the selective one minus a gradient."""
    out = _resolve(False)

    assert out['pulse'] == 'block'
    assert out['duration_s'] == HARD_PULSE_S


def test_a_shaped_default_takes_no_duration_of_its_own() -> None:
    """
    ``Excitation`` already has a shaped default, and a second number here could drift from it.

    The hard pulse is the exception, and only because ``Excitation``'s default is a *shaped*
    pulse's duration -- which is not a duration a block pulse should inherit.
    """
    assert 'duration_s' not in _resolve(True)


# --------------------------------------------------------------- a deliberate choice is honoured
@pytest.mark.parametrize('selective', [True, False])
@pytest.mark.parametrize('pulse', ['sinc', 'slr', 'gauss', 'block'])
def test_an_explicit_pulse_is_never_overridden(selective: bool, pulse: str) -> None:
    """
    A caller who names a pulse family has made a design decision, in either mode.

    A shaped but spatially non-selective pulse is a real design -- a spectrally selective
    excitation, say.  It is an explicit alternative, not the default imaging mode, and the
    difference between those two is the whole subject of this file.
    """
    assert _resolve(selective, rf_pulse=pulse)['pulse'] == pulse


@pytest.mark.parametrize('selective', [True, False])
def test_an_explicit_duration_is_never_overridden(selective: bool) -> None:
    assert _resolve(selective, rf_duration_s=1.25e-3)['duration_s'] == 1.25e-3


def test_an_explicit_duration_beats_the_hard_pulse_default() -> None:
    """The mode supplies a duration only when the caller did not."""
    assert _resolve(False, rf_duration_s=0.5e-3)['duration_s'] == 0.5e-3
    assert _resolve(False)['duration_s'] == HARD_PULSE_S


# ------------------------------------------------------------------------------- the one refusal
@pytest.mark.parametrize('selective', [True, False])
def test_a_time_bandwidth_product_needs_a_shaped_pulse(selective: bool) -> None:
    """The parameter has no meaning for a block pulse, so it is refused rather than ignored."""
    with pytest.raises(ConfigurationError) as err:
        _resolve(selective, rf_pulse='block', rf_time_bw_product=4.0)

    assert 'block pulse' in str(err.value)
    assert 'slab_thickness_mm' in str(err.value), 'the refusal must name the selection argument'


def test_the_refusal_reaches_a_non_selective_default_too() -> None:
    """Nothing was named ``block`` here -- the *mode* chose it, and the refusal still applies."""
    with pytest.raises(ConfigurationError, match='block pulse'):
        _resolve(False, rf_time_bw_product=6.0)


@pytest.mark.parametrize('selective', [True, False])
def test_a_time_bandwidth_product_passes_through_to_a_shaped_pulse(selective: bool) -> None:
    assert _resolve(selective, rf_pulse='sinc', rf_time_bw_product=6.0)['time_bw_product'] == 6.0


# ----------------------------------------------------------- the emitted pair, through GRE3DTR
@pytest.fixture
def opts() -> pp.Opts:
    return pp.Opts(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6, B0=3.0)


def _emitted(tr: GRE3DTR) -> tuple[int, float, set[str]]:
    """``(RF sample count, RF duration, gradient channels)`` of the excitation it actually plays."""
    block = tr.exc(phase_deg=0.0)
    events = [ev for _, ev, _ in sc.flatten(block)]
    rf = next(ev for ev in events if getattr(ev, 'type', None) == 'rf')
    channels = {ev.channel for ev in events if getattr(ev, 'type', None) in {'grad', 'trap'}}
    return int(np.asarray(rf.signal).size), float(pp.calc_duration(rf)), channels


def test_the_two_modes_emit_different_rf_families(opts: pp.Opts) -> None:
    """
    The assertion the rule exists for, stated on the waveform rather than on an argument.

    A slab gets thousands of samples over milliseconds; no slab gets a two-sample block an order
    of magnitude shorter.  If the non-selective mode were the selective one with ``Gz`` removed,
    these two RF waveforms would be identical and only the channel set would differ.
    """
    base = {'opts': opts, 'fov_mm': (200.0, 200.0, 400.0), 'matrix': (64, 64, 16)}
    slab_n, slab_s, slab_ch = _emitted(GRE3DTR(**base, slab_thickness_mm=140.0))
    hard_n, hard_s, hard_ch = _emitted(GRE3DTR(**base, slab_thickness_mm=None))

    assert slab_ch == {'z'} and hard_ch == set(), 'only the slab mode plays a selection gradient'
    assert slab_n > 100 and hard_n <= 4, 'and the RF families differ, which is the actual mode'
    assert hard_s < slab_s / 5


def test_a_non_selective_slab_still_honours_an_explicit_shaped_pulse(opts: pp.Opts) -> None:
    """The expert path: a shaped pulse that selects no slab is a design, and stays available."""
    n, _, channels = _emitted(GRE3DTR(opts=opts, fov_mm=(200.0, 200.0, 400.0),
                                      matrix=(64, 64, 16), slab_thickness_mm=None,
                                      rf_pulse='sinc'))

    assert channels == set(), 'still no selection gradient'
    assert n > 100, 'but the caller asked for a shaped pulse and got one'
