"""
Two things worth doing to a ``pypulseq.Opts`` that ``Opts`` itself does not do.

seqcraft describes a scanner with the **official** :class:`pypulseq.Opts` and nothing else.  There
is no ``System``, no ``Limits``, and no named limit regime: a sequence that designs one part
against derated limits carries a second ``Opts``, which is one object more and one concept fewer.
Build the first one the ordinary way::

    opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s', B0=3.0,
                   rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)

There is deliberately no wrapper around that call.  What this module adds is the two operations
that constructor makes awkward or unsafe.

:func:`derate` -- **because deriving one ``Opts`` from another by hand loses the rest of it.**
``Opts`` fills every argument you omit from the *process-global* ``Opts.default``, so::

    pp.Opts(max_grad=opts.max_grad * 0.9).rf_dead_time      # -> 0, not yours

The derated scanner silently returns to zero dead times, pypulseq's rasters and a foreign gamma.
:func:`derate` copies all eighteen fields and changes only what it was asked to, which is what the
multi-regime ``System`` needed a cross-regime consistency check to guarantee.

:func:`from_scanner` -- **because the amplitudes are data somebody else maintains.**
`PulseqSystems <https://github.com/nimpulseq/PulseqSystems>`_ carries vendor ``max_grad``,
``max_slew`` and ``B0``, so those need not be copied off a spec sheet into your code.  It carries
nothing else, though: rasters, dead times, ringdown and sample limits belong to the *installation*
rather than to the magnet, and no vendor database has them.

That gap is the reason :func:`from_scanner` takes the four site constants as **required keyword
arguments**.  pypulseq defaults ``rf_dead_time``, ``rf_ringdown_time`` and ``adc_dead_time`` to
**zero** -- wrong on every real scanner -- and a sequence built on those defaults compiles cleanly,
validates cleanly, and is refused or silently mangled at the console.  A lookup that returned an
``Opts`` with three zeros in it would be the shortest path to that file.

Examples
--------
>>> import seqcraft as sc
>>> SITE = {                                   # this installation, not this scanner model
...     'rf_dead_time': 100e-6,
...     'rf_ringdown_time': 30e-6,
...     'adc_dead_time': 10e-6,
...     'max_b1': sc.convert(20, 'uT', 'Hz'),
...     'adc_samples_limit': 8192,
... }
>>> opts = sc.opts.from_scanner('Siemens Healthineers', 'MAGNETOM Cima X',
...                             'Gemini Gradients', **SITE)      # doctest: +SKIP
>>> epi = sc.opts.derate(opts, grad=0.85)      # what regime='epi' used to mean  # doctest: +SKIP

Notes
-----
Nothing here ever sets pypulseq's process-global ``Opts.default``.  ``set_as_default=`` and
``reset_default=`` are rejected rather than forwarded: a global scanner is how a sequence comes to
depend on import order, which is the bug this whole layer is arranged to prevent.

The `pulseq-systems <https://pypi.org/project/pulseq-systems/>`_ dependency is optional and used
only by :func:`from_scanner` -- ``pip install seqcraft[systems]``.  :func:`derate` needs nothing
beyond pypulseq.
"""

from __future__ import annotations

import difflib
import inspect
from typing import TYPE_CHECKING, Any

from pypulseq.opts import Opts

from ..core.errors import (
    ConfigurationError,
    MissingExtraError,
    UnknownFieldError,
    format_error,
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

__all__ = ['derate', 'from_scanner']

#: Every keyword ``Opts.__init__`` accepts, less the two that mutate the process-global default.
#: Read from the signature rather than listed here, so a pypulseq that grows a field does not need
#: an edit in this file.
_ACCEPTED: frozenset[str] = frozenset(
    inspect.signature(Opts.__init__).parameters
) - {'self', 'set_as_default', 'reset_default'}

#: The two ``Opts`` arguments this module refuses to forward, and why.
_FORBIDDEN = {
    'set_as_default': 'it would make every later Opts() inherit this scanner',
    'reset_default': 'it would discard a default another caller set',
}

#: Constructor arguments naming the unit the amplitude fields are quoted in.  Meaningful when
#: building from a spec sheet, and never when deriving from an ``Opts`` that already stores Hz.
_UNIT_HINTS = frozenset({'b1_unit', 'grad_unit', 'slew_unit'})


def _reject_forbidden(given: Iterable[str]) -> None:
    """
    Raise if a caller tried to reach pypulseq's process-global default through this layer.

    The one thing ``Opts`` will not refuse on its own: it honours ``set_as_default=True``, and
    from then on every ``Opts()`` built anywhere in the process inherits this scanner.
    """
    for name in given:
        if name in _FORBIDDEN:
            msg = format_error(
                f'seqcraft never passes {name}= to Opts.',
                {'reason': _FORBIDDEN[name]},
                [
                    'hold the Opts in a variable and pass it explicitly -- to pypulseq as '
                    'system=, to seqcraft as sc.compile(tree, opts)',
                ],
            )
            raise ConfigurationError(msg)


def _checked(given: Mapping[str, Any], allowed: frozenset[str], what: str) -> dict[str, Any]:
    """Return `given` unchanged, having verified every key names a settable ``Opts`` field."""
    _reject_forbidden(given)
    for name in given:
        if name not in allowed:
            fields: dict[str, object] = {'known': ', '.join(sorted(allowed))}
            close = difflib.get_close_matches(name, sorted(allowed), n=3)
            if close:
                fields['did you mean'] = ', '.join(close)
            raise UnknownFieldError(format_error(f'Opts has no {what} {name!r}.', fields))
    return dict(given)


def derate(opts: Opts, *, grad: float = 1.0, slew: float = 1.0, **overrides: Any) -> Opts:
    """
    Return a copy of `opts` with the limits scaled -- what a named regime used to be.

    Parameters
    ----------
    opts
        The scanner to derate.  Never mutated.
    grad, slew
        Multiplicative factors on ``max_grad`` and ``max_slew``.
    **overrides
        Absolute values for any ``Opts`` field, applied after the factors.  In **canonical units**
        -- Hz/m, Hz/m/s, Hz -- because that is what an existing ``Opts`` stores.  ``grad_unit`` and
        its siblings are therefore rejected: accepting ``grad_unit='mT/m'`` here would reinterpret
        the carried-over Hz/m amplitude as mT/m and derate the scanner by a factor of 42 576, which
        ``Opts`` would do silently.

    Returns
    -------
    Opts
        A new object.  Every field not scaled or overridden is copied, so a derated scanner cannot
        drift from the one it came from on rasters, dead times or gamma -- the disagreement the
        old multi-regime ``System`` needed a consistency check to catch, and which a hand-written
        ``pp.Opts(max_grad=...)`` reintroduces by falling back to the process-global default.

    Examples
    --------
    >>> import pypulseq as pp
    >>> import seqcraft as sc
    >>> opts = pp.Opts(max_grad=80, grad_unit='mT/m', max_slew=200, slew_unit='T/m/s',
    ...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    >>> epi = sc.opts.derate(opts, grad=0.9, slew=0.6)
    >>> round(epi.max_grad / opts.max_grad, 2), round(epi.max_slew / opts.max_slew, 2)
    (0.9, 0.6)
    >>> epi.rf_dead_time == opts.rf_dead_time      # everything else is carried across
    True
    >>> opts.max_slew == sc.opts.derate(opts).max_slew          # derate() never mutates
    True

    A peak-B1 ceiling is an override, not a special method:

    >>> round(sc.convert(sc.opts.derate(opts, max_b1=sc.convert(12, 'uT', 'Hz')).max_b1,
    ...                  'Hz', 'uT'), 1)
    12.0
    """
    fields = {key: value for key, value in vars(opts).items() if key in _ACCEPTED}
    fields['max_grad'] = float(opts.max_grad) * float(grad)
    fields['max_slew'] = float(opts.max_slew) * float(slew)
    # `rise_time` is a *rule* for deriving max_slew at construction, not a limit in its own right:
    # carried forward it would recompute the slew just set, so `slew=0.6` would silently do
    # nothing.  The number it implied is already in max_slew, and pypulseq reads it nowhere else.
    fields['rise_time'] = None
    fields.update(_checked(overrides, _ACCEPTED - _UNIT_HINTS, 'override'))
    return Opts(**fields)


def from_scanner(
    manufacturer: str,
    model: str,
    gradient: str | None = None,
    *,
    rf_dead_time: float,
    rf_ringdown_time: float,
    adc_dead_time: float,
    max_b1: float,
    **overrides: Any,
) -> Opts:
    """
    Look a scanner's amplitude limits up in PulseqSystems and build an ``Opts`` around them.

    Parameters
    ----------
    manufacturer, model, gradient
        As in ``pulseq_systems.get_pulseq_specs``.  `gradient` names a coil option and is needed
        only for models that have more than one.  An unknown name raises with the alternatives
        listed.
    rf_dead_time, rf_ringdown_time, adc_dead_time
        Seconds.  **Required, because the database cannot supply them and no correct default
        exists.**  pypulseq defaults all three to zero, which is wrong on every real scanner and
        produces a file the console refuses.
    max_b1
        Peak transmit amplitude to design against, in Hz -- ``sc.convert(20, 'uT', 'Hz')``.
        Required for the same reason: it belongs to the transmit chain and the loading, so no
        vendor database can know it.  Take it from the reference voltage the scanner reports and
        keep a margin; a 4 ms 180 degree SLR peaks at 15.4 uT and was rejected at 362 V.
    **overrides
        Any other ``Opts`` field -- ``adc_samples_limit`` (8192 is the common Siemens value), the
        rasters, ``gamma`` for another nucleus.

    Returns
    -------
    Opts
        A plain ``pypulseq.Opts``.  No wrapper, so pypulseq's own ``make_*`` take it unchanged.

    Raises
    ------
    MissingExtraError
        ``pulseq-systems`` is not installed.
    ConfigurationError
        The manufacturer, model or gradient is not in the database.  The message lists what is.

    Notes
    -----
    ``get_pulseq_specs`` returns **per-axis** amplitudes: where a vendor quotes a vector-norm
    figure the database divides it by sqrt(3), so a Prisma reads 46.19 mT/m rather than 80.  That
    is the convention the compiler wants -- it treats a per-axis excess as an error and the vector
    norm as a warning -- so the two agree without adjustment.

    Prefer :func:`derate` over the lookup's own ``scale_gradients`` / ``scale_slew_rate``: derating
    afterwards leaves the un-derated ``Opts`` in hand to validate the finished sequence against.
    """
    systems = _pulseq_systems()
    try:
        specs = systems.get_pulseq_specs(manufacturer, model, gradient)
    except KeyError as err:
        raise ConfigurationError(
            _lookup_message(systems, manufacturer, model, gradient)
        ) from err

    fields = _checked(specs, _ACCEPTED | _UNIT_HINTS, 'spec')
    fields.update(
        rf_dead_time=float(rf_dead_time),
        rf_ringdown_time=float(rf_ringdown_time),
        adc_dead_time=float(adc_dead_time),
        max_b1=float(max_b1),
    )
    fields.update(_checked(overrides, _ACCEPTED | _UNIT_HINTS, 'override'))
    return Opts(**fields)


def _pulseq_systems() -> Any:
    """Import ``pulseq_systems``, or raise naming the extra that provides it."""
    try:
        import pulseq_systems  # noqa: PLC0415
    except ImportError as err:
        msg = format_error(
            'looking a scanner up by name needs the pulseq-systems package.',
            {'provides': 'vendor max_grad, max_slew and B0 as data, maintained outside seqcraft'},
            [
                'pip install "seqcraft[systems]"',
                'or build the Opts directly: pp.Opts(max_grad=..., grad_unit=..., ...) -- the '
                'lookup only saves copying three numbers off a spec sheet',
            ],
        )
        raise MissingExtraError(msg) from err
    return pulseq_systems


def _lookup_message(systems: Any, manufacturer: str, model: str, gradient: str | None) -> str:
    """
    Build the error for a name PulseqSystems does not know, listing what it does.

    ``get_pulseq_specs`` raises a bare ``KeyError`` naming the string it failed on, which does not
    say *which* of the three levels was wrong -- and the gradient names carry a suffix
    (``'Gemini Gradients'``, not ``'Gemini'``) that nothing tells you about.
    """
    manufacturers = list(systems.list_manufacturers())
    if manufacturer not in manufacturers:
        return format_error(
            f'unknown manufacturer {manufacturer!r}.',
            {'available': ', '.join(manufacturers)},
        )
    models = list(systems.list_models(manufacturer))
    if model not in models:
        return format_error(
            f'unknown model {model!r} for {manufacturer!r}.',
            {'available': ', '.join(models)},
        )
    gradients = list(systems.list_gradients(manufacturer, model))
    return format_error(
        f'unknown gradient option {gradient!r} for {manufacturer!r} {model!r}.',
        {'available': ', '.join(gradients) or '(none; pass gradient=None)'},
        ['the names carry a suffix -- "Gemini Gradients", not "Gemini"'],
    )
