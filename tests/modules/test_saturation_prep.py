"""
``SaturationPrep`` -- the sign of the offset first, because it is the failure that looks fine.

A wrong chemical-shift sign produces legal Pulseq, legal timing, legal gradients and a correct
waveform, and saturates water instead of fat.  Nothing downstream notices: the sequence compiles,
every k-space check passes, and the image comes back with the wrong tissue suppressed.  The two
reference implementations reach the same number by different routes -- one carries the sign in the
ppm constant, the other applies it at the point of use -- so an implementation that mixed the
conventions would be exactly this wrong.

So the tests trace the whole chain, ``shift_ppm -> offset_hz -> emitted rf.freq_offset``, as far as
the *compiled* sequence rather than stopping at the constructor.  The rest is the module's other
promises: that the pulse is labelled what it is, that no gradient and no rephaser are emitted, and
that the spoiler follows.
"""

from __future__ import annotations

import numpy as np
import pypulseq as pp
import pytest

import seqcraft as sc

#: R1's protocol: pypulseq ``write_epi_se_rs.py``, which is also ``writeEpiSpinEchoRS.m``'s.
FAT_PPM, FLIP_DEG, DURATION_S, BANDWIDTH_HZ = -3.45, 110.0, 8e-3, 424.5


@pytest.fixture
def opts() -> pp.Opts:
    """The reference protocol's scanner, at its 2.89 T."""
    return pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6, B0=2.89)


def prep(opts: pp.Opts, **overrides: object) -> sc.modules.SaturationPrep:
    """The reference protocol, with `overrides` applied."""
    kwargs: dict = {
        'opts': opts, 'shift_ppm': FAT_PPM, 'flip_deg': FLIP_DEG,
        'duration_s': DURATION_S, 'bandwidth_hz': BANDWIDTH_HZ, 'spoil_voxel_mm': 0.1,
    }
    kwargs.update(overrides)
    return sc.modules.SaturationPrep(**kwargs)


def only_rf(seq: pp.Sequence) -> object:
    """The one RF event in a compiled sequence."""
    events = [seq.get_block(n).rf for n in range(1, len(seq.block_events) + 1)
              if getattr(seq.get_block(n), 'rf', None) is not None]
    assert len(events) == 1, f'expected one RF event, found {len(events)}'
    return events[0]


# ------------------------------------------------------------------ the sign, end to end
def test_offset_hz_is_the_chemical_shift_converted_once(opts: pp.Opts) -> None:
    """``shift_ppm * 1e-6 * B0 * gamma``, and the module reports it."""
    expected = FAT_PPM * 1e-6 * opts.B0 * opts.gamma
    assert prep(opts).offset_hz == pytest.approx(expected)
    assert prep(opts).offset_hz == pytest.approx(-424.504008)


def test_fat_saturates_below_water(opts: pp.Opts) -> None:
    """A negative shift is a negative offset.  This is the test the whole module is for."""
    assert prep(opts).offset_hz < 0.0


def test_a_positive_shift_saturates_above_water(opts: pp.Opts) -> None:
    """The sign is carried, not assumed: the module does not force fat's side."""
    assert prep(opts, shift_ppm=+3.45).offset_hz == pytest.approx(+424.504008)


def test_the_compiled_sequence_carries_the_offset_we_computed(opts: pp.Opts) -> None:
    """The chain reaches the emitted event, not just the constructor."""
    module = prep(opts)
    seq = sc.compile(module(), opts, name='fat_sat')
    assert only_rf(seq).freq_offset == pytest.approx(module.offset_hz)
    assert only_rf(seq).freq_offset < 0.0


def test_the_offset_scales_with_field_and_nothing_else_does(opts: pp.Opts) -> None:
    """Metamorphic, and definable against the references before this module existed."""
    low = prep(opts)
    high = prep(pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
                        rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6,
                        B0=2 * opts.B0))
    assert high.offset_hz == pytest.approx(2 * low.offset_hz)
    assert high.time_to_center() == pytest.approx(low.time_to_center())
    assert np.allclose(np.abs(high.rf.signal), np.abs(low.rf.signal))


def test_the_offset_scales_with_the_shift_and_nothing_else_does(opts: pp.Opts) -> None:
    """The other half of the metamorphic pair."""
    one, two = prep(opts), prep(opts, shift_ppm=2 * FAT_PPM)
    assert two.offset_hz == pytest.approx(2 * one.offset_hz)
    assert np.allclose(np.abs(two.rf.signal), np.abs(one.rf.signal))


# ------------------------------------------------------------ the pulse says what it is
def test_the_emitted_pulse_is_labelled_a_saturation(opts: pp.Opts) -> None:
    """Rule B, as a test: the file must tell the same story as the API name."""
    seq = sc.compile(prep(opts)(), opts, name='fat_sat')
    assert only_rf(seq).use == 'saturation'


def test_an_excitation_built_the_same_way_would_not_be(opts: pp.Opts) -> None:
    """
    The gap this module exists to close, pinned.

    ``Excitation`` can already make a spectrally selective pulse -- it says so in its own
    docstring -- but it owns ``use`` and hard-codes ``'excitation'``, so a saturation built that
    way emits a file declaring the pulse an excitation.  If ``Excitation`` ever gains a purpose
    flag this test fails, and that is the conversation it exists to start.
    """
    other = sc.modules.Excitation(opts=opts, flip_deg=FLIP_DEG, thickness_mm=None, pulse='gauss',
                                  duration_s=DURATION_S, pulse_opts={'bandwidth': BANDWIDTH_HZ})
    assert only_rf(sc.compile(other(), opts, name='not_a_sat')).use == 'excitation'


# --------------------------------------------------------- what is deliberately absent
def test_no_gradient_is_played_during_the_pulse(opts: pp.Opts) -> None:
    """Spectral selectivity is the absence of a gradient, not a narrow one."""
    seq = sc.compile(prep(opts)(), opts, name='fat_sat')
    first = seq.get_block(1)
    assert all(getattr(first, axis, None) is None for axis in ('gx', 'gy', 'gz'))


def test_nothing_is_rephased(opts: pp.Opts) -> None:
    """
    Every gradient in the block is the spoiler, on one axis, with one sign.

    A rephaser would show up as a second z event of opposite sign.  There is no selection
    gradient to rephase, which is the point -- so the check is that the block's total z moment is
    the spoiler's and nothing cancels it.
    """
    module = prep(opts)
    events = [event for _, event, _ in sc.flatten(module()) if getattr(event, 'type', '') == 'grad'
              or getattr(event, 'type', '') == 'trap']
    assert len(events) == 1
    assert events[0].channel == 'z'


def test_the_spoiler_follows_the_pulse(opts: pp.Opts) -> None:
    """Order matters: a spoiler before the pulse spoils nothing."""
    module = prep(opts)
    times = {getattr(event, 'type', ''): t for t, event, _ in sc.flatten(module())}
    assert times['rf'] < times['trap']
    assert times['trap'] == pytest.approx(pp.calc_duration(module.rf))


def test_the_spoiler_moment_is_the_product_asked_for(opts: pp.Opts) -> None:
    """``cycles / voxel_m``, which is the references' ``area = 1/1e-4`` at 0.1 mm."""
    module = prep(opts, spoil_cycles_per_voxel=1.0, spoil_voxel_mm=0.1)
    trap = next(e for _, e, _ in sc.flatten(module()) if getattr(e, 'type', '') == 'trap')
    assert float(trap.area) == pytest.approx(1.0 / 1e-4)


# --------------------------------------------------------------------- timing it reports
def test_time_to_center_matches_the_emitted_waveform(opts: pp.Opts) -> None:
    """Same origin and same form as ``IRPrep`` and ``Excitation``, so a timeline can add them."""
    module = prep(opts)
    assert module.time_to_center() == pytest.approx(
        float(module.rf.delay) + float(pp.calc_rf_center(module.rf)[0]),
    )
    assert module.time_to_center() == pytest.approx(opts.rf_dead_time + DURATION_S / 2, abs=2e-5)


def test_the_block_holds_the_pulse_and_the_spoiler(opts: pp.Opts) -> None:
    """Two nodes, and the duration is the pulse plus the spoiler."""
    module = prep(opts)
    block = module()
    assert block.duration == pytest.approx(
        pp.calc_duration(module.rf) + module.spoiler.duration,
    )


# ------------------------------------------------------------------------- the refusals
def test_a_band_that_reaches_water_is_refused(opts: pp.Opts) -> None:
    """The one physical relationship the module owns, and it refuses rather than warns."""
    with pytest.raises(sc.ConfigurationError, match='reaches water'):
        prep(opts, bandwidth_hz=2 * abs(prep(opts).offset_hz) + 1.0)


def test_the_band_edge_is_reported_so_the_margin_is_visible(opts: pp.Opts) -> None:
    """A caller choosing a duration is choosing this number whether or not they know it."""
    module = prep(opts)
    assert module.band_edge_hz == pytest.approx(abs(module.offset_hz) - BANDWIDTH_HZ / 2)
    assert module.band_edge_hz > 0.0


def test_water_itself_is_refused(opts: pp.Opts) -> None:
    """``shift_ppm=0.0`` would saturate the signal the sequence exists to acquire."""
    with pytest.raises(sc.ConfigurationError, match='water itself'):
        prep(opts, shift_ppm=0.0)


def test_a_block_pulse_is_refused_with_the_reason(opts: pp.Opts) -> None:
    """A rectangular envelope's spectral profile has sidelobes that reach water."""
    with pytest.raises(sc.ConfigurationError, match='sidelobes reach water'):
        prep(opts, pulse='block')


@pytest.mark.parametrize('bad', ['flip_deg', 'duration_s', 'bandwidth_hz', 'spoil_voxel_mm'])
def test_non_positive_protocol_values_are_refused(opts: pp.Opts, bad: str) -> None:
    """Every protocol value is required and positive; none of them has a default to fall back on."""
    with pytest.raises(sc.ConfigurationError):
        prep(opts, **{bad: 0.0})


def test_pulse_opts_cannot_reach_what_the_module_owns(opts: pp.Opts) -> None:
    """``freq_offset`` is the module's; a caller setting it would break the sign chain."""
    with pytest.raises(sc.ConfigurationError, match='not accepted'):
        prep(opts, pulse_opts={'freq_offset': 0.0})


def test_pulse_opts_forwards_a_design_parameter(opts: pp.Opts) -> None:
    """``apodization`` is R3's, and it is a pulse-design choice rather than a module one."""
    plain, apodised = prep(opts), prep(opts, pulse_opts={'apodization': 0.42})
    assert not np.allclose(np.abs(plain.rf.signal), np.abs(apodised.rf.signal))


# ---------------------------------------------------------------- the other two shapes
@pytest.mark.parametrize('shape', ['gauss', 'sinc', 'slr'])
def test_every_offered_shape_builds_and_keeps_the_contract(opts: pp.Opts, shape: str) -> None:
    """Whatever the envelope, the label, the offset and the missing gradient hold."""
    seq = sc.compile(prep(opts, pulse=shape)(), opts, name=f'fat_sat_{shape}')
    rf = only_rf(seq)
    assert rf.use == 'saturation'
    assert rf.freq_offset < 0.0
    assert seq.get_block(1).gz is None


# ------------------------------------------------------- where B0 comes from, and its trap
def test_b0_is_read_from_opts_and_reported(opts: pp.Opts) -> None:
    """The middle term of the chain is a fact the caller can read back, not an implication."""
    module = prep(opts)
    assert module.b0_t == opts.B0
    assert module.offset_hz == pytest.approx(module.shift_ppm * 1e-6 * module.b0_t * opts.gamma)


def test_an_opts_without_b0_silently_means_1_5_tesla() -> None:
    """
    pypulseq's ``Opts`` defaults ``B0`` to 1.5 T and says nothing, so this module cannot tell an
    omitted field strength from a deliberate one.

    Pinned rather than worked around: the failure it causes -- a 2.89 T protocol whose author
    forgot ``B0=2.89``, putting fat at -220 Hz instead of -424 -- is legal, plausible and
    invisible downstream, and the defence is that ``b0_t`` is reported and this behaviour is a
    documented fact rather than a discovery.
    """
    silent = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
                     rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    assert silent.B0 == 1.5
    module = prep(silent, bandwidth_hz=200.0)
    assert module.b0_t == 1.5
    assert module.offset_hz == pytest.approx(-220.3, abs=0.1)
