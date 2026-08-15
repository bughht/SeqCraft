"""
The two operations ``pp.Opts`` does not do for itself.

``sc.opts`` deliberately does not wrap ``Opts``.  It holds :func:`derate`, because deriving one
``Opts`` from another by hand silently loses sixteen fields to the process-global default, and
:func:`from_scanner`, because the vendor amplitudes are data somebody else maintains -- and
because a lookup that returned an ``Opts`` with pypulseq's zero dead times in it would be the
shortest available path to a file the console refuses.
"""

from __future__ import annotations

import pytest
from pypulseq.opts import Opts

import seqcraft as sc
from seqcraft.core.errors import ConfigurationError, MissingExtraError, UnknownFieldError

# One installation's constants.  Not on any spec sheet, which is the whole point.
SITE = {
    'rf_dead_time': 100e-6,
    'rf_ringdown_time': 30e-6,
    'adc_dead_time': 10e-6,
    'max_b1': 851.52,
}


@pytest.fixture
def opts() -> Opts:
    """A 3 T scanner built the ordinary way -- there is no seqcraft constructor to use."""
    return Opts(
        max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s', B0=3.0,
        adc_samples_limit=8192, **SITE,
    )


# --------------------------------------------------------------------------------- derate
def test_derate_scales_what_it_is_asked_to(opts: Opts) -> None:
    epi = sc.opts.derate(opts, grad=0.9, slew=0.6)

    assert epi.max_grad == pytest.approx(opts.max_grad * 0.9)
    assert epi.max_slew == pytest.approx(opts.max_slew * 0.6)


def test_derate_copies_everything_else(opts: Opts) -> None:
    """
    The reason this function exists rather than a one-line ``Opts(max_grad=...)``.

    Also what the old multi-regime ``System`` needed a consistency check to guarantee: two regimes
    disagreeing on a raster or on gamma produced an unplayable file with no warning.  Deriving one
    ``Opts`` from another makes the disagreement unreachable instead of merely detected.
    """
    epi = sc.opts.derate(opts, grad=0.9, slew=0.6)

    untouched = set(vars(opts)) - {'max_grad', 'max_slew', 'rise_time'}
    assert {k: getattr(epi, k) for k in untouched} == {k: getattr(opts, k) for k in untouched}


def test_the_hand_written_version_loses_the_site_constants(opts: Opts) -> None:
    """
    The failure `derate` is the fix for, asserted so it cannot quietly stop being true.

    ``Opts`` fills every omitted argument from the process-global default, so the obvious way to
    write a derated scanner silently returns the dead times to zero.
    """
    by_hand = Opts(max_grad=opts.max_grad * 0.9)

    assert by_hand.rf_dead_time == 0
    assert by_hand.adc_dead_time == 0
    assert sc.opts.derate(opts, grad=0.9).rf_dead_time == 100e-6


def test_derate_never_mutates_its_argument(opts: Opts) -> None:
    before = dict(vars(opts))

    other = sc.opts.derate(opts, grad=0.5)

    assert vars(opts) == before
    assert other is not opts


def test_derate_returns_a_plain_opts(opts: Opts) -> None:
    """No subclass and no wrapper, so pypulseq's own ``make_*`` take it unchanged."""
    assert type(sc.opts.derate(opts)) is Opts


def test_derate_takes_absolute_overrides(opts: Opts) -> None:
    quiet = sc.opts.derate(opts, max_b1=sc.convert(12.0, 'uT', 'Hz'))

    assert sc.convert(quiet.max_b1, 'Hz', 'uT') == pytest.approx(12.0)


def test_derate_rejects_a_unit_hint(opts: Opts) -> None:
    """
    The one silent corruption ``Opts`` will not catch on its own.

    ``grad_unit='mT/m'`` with no ``max_grad`` reinterprets the carried-over Hz/m value as mT/m and
    derates the scanner by a factor of 42 576, without raising anything.
    """
    with pytest.raises(UnknownFieldError, match='grad_unit'):
        sc.opts.derate(opts, grad_unit='mT/m')


def test_derate_rejects_a_mistyped_override(opts: Opts) -> None:
    with pytest.raises(UnknownFieldError) as err:
        sc.opts.derate(opts, adc_deadtime=10e-6)

    assert 'adc_dead_time' in str(err.value)


def test_derate_drops_rise_time_so_the_slew_factor_holds() -> None:
    """
    ``rise_time`` is a rule for deriving ``max_slew``, not a limit of its own.

    Carried across a derate it would recompute the slew from the (unchanged) amplitude, so
    ``slew=0.5`` would silently do nothing at all.
    """
    opts = Opts(max_grad=40, grad_unit='mT/m', rise_time=200e-6, **SITE)
    assert opts.max_slew == pytest.approx(opts.max_grad / 200e-6)

    quiet = sc.opts.derate(opts, slew=0.5)

    assert quiet.max_slew == pytest.approx(opts.max_slew * 0.5)
    assert quiet.rise_time is None


# ------------------------------------------------------------------- the process-global default
@pytest.mark.parametrize('forbidden', ['set_as_default', 'reset_default'])
def test_the_global_default_cannot_be_reached_through_this_layer(
    opts: Opts, forbidden: str
) -> None:
    """
    The other thing ``Opts`` will not refuse: it honours ``set_as_default=True``.

    From then on every ``Opts()`` built anywhere in the process inherits this scanner, which is
    how a sequence comes to depend on import order.
    """
    with pytest.raises(ConfigurationError, match=forbidden):
        sc.opts.derate(opts, **{forbidden: True})


def test_deriving_an_opts_leaves_the_global_default_alone(opts: Opts) -> None:
    before = dict(vars(Opts.default))

    sc.opts.derate(opts, grad=0.5)

    assert vars(Opts.default) == before


# ----------------------------------------------------------------------------- from_scanner
def test_from_scanner_requires_the_site_constants() -> None:
    """
    A vendor database supplies amplitudes.  It cannot supply an installation's dead times.

    The lookup being a one-liner must not make them optional -- that would put pypulseq's
    zero-dead-time default back within reach by exactly the shortest path.
    """
    pytest.importorskip('pulseq_systems')

    with pytest.raises(TypeError, match='rf_dead_time'):
        sc.opts.from_scanner('Siemens Healthineers', 'MAGNETOM Prisma')


def test_from_scanner_builds_a_usable_opts() -> None:
    pytest.importorskip('pulseq_systems')

    opts = sc.opts.from_scanner('Siemens Healthineers', 'MAGNETOM Prisma',
                                adc_samples_limit=8192, **SITE)

    assert type(opts) is Opts
    assert opts.max_grad > 0
    assert opts.rf_dead_time == 100e-6
    assert opts.adc_samples_limit == 8192


def test_from_scanner_returns_per_axis_amplitudes() -> None:
    """
    The database divides a vector-norm figure by sqrt(3), which is the convention we want.

    The compiler treats a per-axis excess as an error and the vector norm as a warning, so a
    nominal 80 mT/m Prisma must arrive as 46.19 rather than 80 for the two to agree.
    """
    pytest.importorskip('pulseq_systems')

    opts = sc.opts.from_scanner('Siemens Healthineers', 'MAGNETOM Prisma', **SITE)

    assert sc.convert(opts.max_grad, 'Hz/m', 'mT/m') == pytest.approx(80.0 / 3**0.5, rel=1e-3)


@pytest.mark.parametrize(
    ('args', 'expected'),
    [
        (('Nonesuch Medical', 'MAGNETOM Prisma', None), 'manufacturer'),
        (('Siemens Healthineers', 'MAGNETOM Nonesuch', None), 'model'),
        (('Siemens Healthineers', 'MAGNETOM Cima X', 'Gemini'), 'gradient'),
    ],
)
def test_an_unknown_name_says_which_level_and_lists_the_alternatives(
    args: tuple[str, str, str | None], expected: str
) -> None:
    """
    ``get_pulseq_specs`` raises a bare ``KeyError`` naming the string, not the level.

    The third case is the one that costs an afternoon: the gradient names carry a suffix, so
    ``'Gemini'`` is wrong and ``'Gemini Gradients'`` is right, and nothing says so.
    """
    pytest.importorskip('pulseq_systems')

    with pytest.raises(ConfigurationError) as err:
        sc.opts.from_scanner(*args, **SITE)

    assert expected in str(err.value)
    assert 'available' in str(err.value)


def test_missing_package_names_the_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    """``pp.Opts`` still works without it, so the message says so rather than only "install"."""
    import builtins

    real_import = builtins.__import__

    def refuse(name: str, *args: object, **kwargs: object) -> object:
        if name == 'pulseq_systems':
            raise ImportError(name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', refuse)

    with pytest.raises(MissingExtraError) as err:
        sc.opts.from_scanner('Siemens Healthineers', 'MAGNETOM Prisma', **SITE)

    assert 'seqcraft[systems]' in str(err.value)
    assert 'pp.Opts' in str(err.value)
