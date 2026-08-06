"""
Hardware description: :class:`Limits`, :class:`System`, and hardware-model loading.

A :class:`System` holds **named ``Opts`` regimes** rather than a single limit set.  Real
sequences need more than one: a diffusion encoding block may run at the full gradient
amplitude while the EPI readout must be derated for peripheral-nerve stimulation.  The
reference implementation modelled this as two unrelated attributes, ``self.system`` and
``self.epi_system``, each built from a separately copy-pasted literal dict -- so nothing
checked that the two agreed on raster times, and a mismatch would have produced an
unplayable file with no warning.  Here they are two regimes of one ``System`` and
agreement is asserted at construction.

Two further rules, both enforced:

* **seqcraft never passes ``system=None`` to pypulseq.**  Doing so silently falls back to
  the process-global ``Opts.default``, which makes a sequence depend on import order.
  Modules always resolve ``self.opts`` from their ``System``.
* **No vendor hardware file is ever read from inside this repository.**  Siemens ``.asc``
  gradient descriptors carry proprietary PNS/CNS response coefficients and forbidden
  acoustic-resonance bands.  :func:`load_hardware` resolves them through the
  ``SEQCRAFT_ASC_DIR`` environment variable only, and :func:`synthetic_hardware` provides
  a vendor-free stand-in so that PNS checks can run in CI.

Examples
--------
>>> sys_ = System.preset('generic_3t')
>>> sys_.gamma
42576000.0
>>> derated = sys_.derate('epi', grad=0.7, slew=0.55)
>>> round(derated.limits('epi').max_slew / sys_.limits('default').max_slew, 3)
0.55
>>> sys_ is derated                       # derate() never mutates
False
"""

from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass, field, replace
from pathlib import Path
from types import MappingProxyType, SimpleNamespace
from typing import TYPE_CHECKING, Any, Literal

from pypulseq.opts import Opts

from .errors import ConfigurationError, format_error
from .timing import Raster
from .units import GAMMA_1H, convert

if TYPE_CHECKING:
    from collections.abc import Mapping

__all__ = ['ASC_ENV_VAR', 'Limits', 'System', 'load_hardware', 'synthetic_hardware']

#: Environment variable naming a directory of vendor ``.asc`` files.  Never a path inside
#: this repository; ``.asc`` is gitignored.
ASC_ENV_VAR = 'SEQCRAFT_ASC_DIR'

#: Fields on which every regime of one System must agree, because blocks designed against
#: different values of these cannot legally coexist in a single ``.seq``.
_SHARED_FIELDS = (
    'grad_raster_time',
    'rf_raster_time',
    'adc_raster_time',
    'block_duration_raster',
    'adc_samples_divisor',
    'gamma',
    'B0',
)


@dataclass(frozen=True)
class Limits:
    """
    Gradient and RF limits in the units they are quoted in.

    Parameters
    ----------
    max_grad_mT_m
        Peak gradient amplitude, mT/m.
    max_slew_T_m_s
        Peak slew rate, T/m/s.
    rise_time_us
        If given, overrides `max_slew_T_m_s` with ``max_grad / rise_time``.  Some
        scanners are specified this way.

    Examples
    --------
    >>> Limits(max_grad_mT_m=80, max_slew_T_m_s=200).max_grad_mT_m
    80
    """

    max_grad_mT_m: float
    max_slew_T_m_s: float
    rise_time_us: float | None = None

    def to_opts(self, base: Opts) -> Opts:
        """Return a copy of `base` with these amplitude and slew limits applied."""
        gamma = base.gamma
        kwargs: dict[str, Any] = {
            'max_grad': convert(self.max_grad_mT_m, 'mT/m', 'Hz/m', gamma=gamma),
            'max_slew': convert(self.max_slew_T_m_s, 'T/m/s', 'Hz/m/s', gamma=gamma),
        }
        if self.rise_time_us is not None:
            kwargs['rise_time'] = convert(self.rise_time_us, 'us', 's')
        return _opts_with(base, **kwargs)


def _opts_with(base: Opts, **changes: Any) -> Opts:
    """Build a new ``Opts`` from `base`, overriding `changes`, without touching globals."""
    kwargs = {
        'max_grad': base.max_grad,
        'max_slew': base.max_slew,
        'max_b1': base.max_b1,
        'max_freq_offset': base.max_freq_offset,
        'rf_dead_time': base.rf_dead_time,
        'rf_ringdown_time': base.rf_ringdown_time,
        'adc_dead_time': base.adc_dead_time,
        'adc_raster_time': base.adc_raster_time,
        'rf_raster_time': base.rf_raster_time,
        'grad_raster_time': base.grad_raster_time,
        'block_duration_raster': base.block_duration_raster,
        'adc_samples_limit': base.adc_samples_limit,
        'rf_samples_limit': base.rf_samples_limit,
        'adc_samples_divisor': base.adc_samples_divisor,
        'gamma': base.gamma,
        'B0': base.B0,
    }
    kwargs.update(changes)
    # grad_unit/slew_unit default to Hz/m and Hz/m/s, which is what we pass.
    return Opts(**kwargs)


@dataclass(frozen=True, eq=False)
class System:
    """
    A scanner: one or more named limit regimes plus an optional hardware model.

    Parameters
    ----------
    regimes
        ``name -> Opts``.  Must contain ``'default'``.
    hardware
        Optional PNS/acoustic hardware model, as returned by
        :func:`load_hardware` or :func:`synthetic_hardware`.
    name
        Free-text scanner label, written into the provenance sidecar.
    source
        Where the hardware model came from.  For a vendor ``.asc`` this records the file
        *name* and a hash, never its contents.

    Notes
    -----
    ``__post_init__`` asserts that every regime agrees on the raster times, gamma, B0 and
    ADC sample divisor (see :data:`_SHARED_FIELDS`).  Nothing in pypulseq checks this, and
    a disagreement produces a file the scanner will reject.

    Examples
    --------
    >>> s = System.preset('generic_3t')
    >>> sorted(s.regime_names)
    ['default']
    >>> s.block_raster
    Raster(block, 10 us)
    """

    regimes: Mapping[str, Opts]
    hardware: SimpleNamespace | None = None
    name: str = 'unnamed'
    source: str | None = None
    _fingerprint: str = field(default='', repr=False, compare=False)

    # ------------------------------------------------------------------ construction
    def __post_init__(self) -> None:
        """Validate that 'default' exists and that all regimes agree on shared physics."""
        if 'default' not in self.regimes:
            msg = format_error(
                'System requires a regime named "default".',
                {'given': ', '.join(sorted(self.regimes)) or '(none)'},
            )
            raise ConfigurationError(msg)

        base = self.regimes['default']
        for regime, opts in self.regimes.items():
            if regime == 'default':
                continue
            for attr in _SHARED_FIELDS:
                a, b = getattr(base, attr), getattr(opts, attr)
                if a != b:
                    msg = format_error(
                        f'regimes "default" and {regime!r} disagree on {attr}.',
                        {'default': a, regime: b},
                        [
                            'build every regime with System.derate(), which copies the '
                            'shared fields',
                        ],
                    )
                    raise ConfigurationError(msg)

        object.__setattr__(self, 'regimes', MappingProxyType(dict(self.regimes)))
        object.__setattr__(self, '_fingerprint', self._compute_fingerprint())

    @classmethod
    def from_limits(
        cls,
        default: Limits,
        *,
        name: str = 'unnamed',
        b0_T: float = 3.0,
        gamma: float = GAMMA_1H,
        rf_dead_time_us: float = 100.0,
        rf_ringdown_time_us: float = 30.0,
        adc_dead_time_us: float = 10.0,
        grad_raster_us: float = 10.0,
        rf_raster_us: float = 1.0,
        adc_raster_ns: float = 100.0,
        block_raster_us: float = 10.0,
        adc_samples_divisor: int = 4,
        max_b1_uT: float = 20.0,
        adc_samples_limit: int = 0,
        rf_samples_limit: int = 0,
        hardware: SimpleNamespace | None = None,
        **named: Limits,
    ) -> System:
        """
        Build a :class:`System` from human-quoted limits.

        Parameters
        ----------
        default
            The ``'default'`` regime.
        b0_T
            Field strength in tesla, used for chemical-shift offsets (fat saturation).
        adc_samples_divisor
            Siemens requires the ADC sample count to be a multiple of 4.
        **named
            Additional regimes, e.g. ``epi=Limits(...)``.

        Examples
        --------
        >>> s = System.from_limits(Limits(80, 200), name='demo',
        ...                        epi=Limits(56, 110))
        >>> sorted(s.regime_names)
        ['default', 'epi']
        """
        base = Opts(
            max_grad=convert(default.max_grad_mT_m, 'mT/m', 'Hz/m', gamma=gamma),
            max_slew=convert(default.max_slew_T_m_s, 'T/m/s', 'Hz/m/s', gamma=gamma),
            max_b1=convert(max_b1_uT, 'uT', 'Hz', gamma=gamma),
            rf_dead_time=convert(rf_dead_time_us, 'us', 's'),
            rf_ringdown_time=convert(rf_ringdown_time_us, 'us', 's'),
            adc_dead_time=convert(adc_dead_time_us, 'us', 's'),
            grad_raster_time=convert(grad_raster_us, 'us', 's'),
            rf_raster_time=convert(rf_raster_us, 'us', 's'),
            adc_raster_time=convert(adc_raster_ns, 'ns', 's'),
            block_duration_raster=convert(block_raster_us, 'us', 's'),
            adc_samples_divisor=adc_samples_divisor,
            adc_samples_limit=adc_samples_limit,
            rf_samples_limit=rf_samples_limit,
            gamma=gamma,
            B0=b0_T,
        )
        if default.rise_time_us is not None:
            base = default.to_opts(base)
        regimes = {'default': base}
        regimes.update({key: lim.to_opts(base) for key, lim in named.items()})
        return cls(regimes=regimes, hardware=hardware, name=name)

    @classmethod
    def preset(cls, key: Literal['generic_3t', 'generic_1p5t', 'prisma', 'vida', 'cima_x']) -> System:
        """
        Return a named scanner preset.

        The numbers here are **published nominal specifications only** -- peak gradient,
        peak slew, field strength -- not vendor hardware descriptors.  No PNS model and no
        acoustic-resonance data is bundled with seqcraft; attach those at run time with
        :meth:`with_hardware` and :func:`load_hardware`.

        Examples
        --------
        >>> System.preset('prisma').limits('default').max_grad  # 80 mT/m in Hz/m
        3406080.0
        """
        presets: dict[str, tuple[Limits, float, str]] = {
            'generic_3t': (Limits(40.0, 150.0), 3.0, 'generic 3T'),
            'generic_1p5t': (Limits(33.0, 120.0), 1.5, 'generic 1.5T'),
            'prisma': (Limits(80.0, 200.0), 2.8936, 'Siemens Prisma'),
            'vida': (Limits(60.0, 200.0), 2.8936, 'Siemens Vida'),
            'cima_x': (Limits(200.0, 200.0), 2.8936, 'Siemens Cima.X'),
        }
        if key not in presets:
            msg = format_error(
                f'unknown System preset {key!r}.',
                {'available': ', '.join(sorted(presets))},
            )
            raise ConfigurationError(msg)
        limits, b0, label = presets[key]
        # Per-event sample limits are the *interpreter's*, not the amplifier's, and pypulseq
        # defaults them to 0 meaning "unlimited" -- which is true of the file format and false of
        # every scanner.  8192 ADC samples per event is the common Siemens value; a longer readout
        # has to be split into several ADCs.  Set them from your own installation if it differs:
        # too high and the sequence is refused at the console, too low and it is refused here.
        siemens = key in {'prisma', 'vida', 'cima_x'}
        return cls.from_limits(
            limits, name=label, b0_T=b0,
            adc_samples_limit=8192 if siemens else 0,
            rf_samples_limit=0,
        )

    # ------------------------------------------------------------------------ access
    @property
    def regime_names(self) -> tuple[str, ...]:
        """The names of every available regime."""
        return tuple(self.regimes)

    @property
    def default(self) -> Opts:
        """The ``'default'`` regime's ``Opts``."""
        return self.regimes['default']

    def limits(self, regime: str = 'default') -> Opts:
        """
        Return the ``Opts`` for `regime`.

        Raises
        ------
        ConfigurationError
            If `regime` is unknown.  Never silently falls back to ``'default'``, which
            would let a typo run the readout at full diffusion amplitude.
        """
        try:
            return self.regimes[regime]
        except KeyError:
            msg = format_error(
                f'unknown regime {regime!r}.',
                {'available': ', '.join(sorted(self.regimes))},
                ["add it with System.derate('name', grad=..., slew=...)"],
            )
            raise ConfigurationError(msg) from None

    def derate(
        self,
        name: str,
        *,
        grad: float = 1.0,
        slew: float = 1.0,
        max_grad_mT_m: float | None = None,
        max_slew_T_m_s: float | None = None,
    ) -> System:
        """
        Return a **new** System with an additional derated regime.

        Parameters
        ----------
        name
            Name of the new regime.
        grad, slew
            Multiplicative factors applied to the default regime's limits.
        max_grad_mT_m, max_slew_T_m_s
            Absolute overrides; take precedence over the factors.

        Examples
        --------
        >>> s = System.preset('cima_x').derate('epi', grad=0.9, slew=0.6)
        >>> round(s.limits('epi').max_slew / s.default.max_slew, 2)
        0.6
        """
        base = self.default
        max_grad = (
            self.convert(max_grad_mT_m, 'mT/m', 'Hz/m')
            if max_grad_mT_m is not None
            else base.max_grad * grad
        )
        max_slew = (
            self.convert(max_slew_T_m_s, 'T/m/s', 'Hz/m/s')
            if max_slew_T_m_s is not None
            else base.max_slew * slew
        )
        regimes = dict(self.regimes)
        regimes[name] = _opts_with(base, max_grad=max_grad, max_slew=max_slew)
        return replace(self, regimes=regimes, _fingerprint='')

    def with_hardware(self, hardware: SimpleNamespace, *, source: str | None = None) -> System:
        """Return a new System carrying a PNS/acoustic hardware model."""
        return replace(self, hardware=hardware, source=source, _fingerprint='')

    def with_max_b1(self, max_b1_uT: float) -> System:
        """
        Return a new System with a different peak-B1 limit on **every** regime.

        The limit belongs to the transmit chain and the loading, not to the gradient system, so no
        preset can know it: take it from the reference voltage the scanner reports and keep a margin.
        The default of 20 uT is generous enough to pass pulses a real amplifier refuses -- a 4 ms
        180 degree SLR peaks at 15.4 uT and was rejected at 362 V.

        Examples
        --------
        >>> import seqcraft as sc
        >>> system = sc.System.preset('cima_x').with_max_b1(12.0)
        >>> round(system.convert(system.default.max_b1, 'Hz', 'uT'), 1)
        12.0
        """
        limit = self.convert(float(max_b1_uT), 'uT', 'Hz')
        regimes = {
            name: _opts_with(opts, max_b1=limit) for name, opts in self.regimes.items()
        }
        return replace(self, regimes=regimes, _fingerprint='')

    # ------------------------------------------------- shared physics, single source
    @property
    def gamma(self) -> float:
        """Gyromagnetic ratio, Hz/T."""
        return float(self.default.gamma)

    @property
    def b0_T(self) -> float:
        """Field strength, tesla."""
        return float(self.default.B0)

    def convert(self, value: float, from_unit: str, to_unit: str | None = None) -> float:
        """
        Convert a value between units using **this scanner's** gamma and Larmor frequency.

        The same function as :func:`seqcraft.units.convert`, with `gamma` and `f0` filled in, so a
        chemical shift or a B1 amplitude cannot silently pick up the proton value on a system set up
        for another nucleus.

        Examples
        --------
        >>> import seqcraft as sc
        >>> system = sc.System.preset('prisma')
        >>> round(system.convert(system.default.max_grad, 'Hz/m', 'mT/m'), 1)
        80.0
        >>> round(system.convert(-3.4, 'ppm', 'Hz'))         # fat/water shift at 2.894 T
        -419
        """
        return convert(value, from_unit, to_unit, gamma=self.gamma, f0=self.gamma * self.b0_T)

    # ------------------------------------------------------------------------- rasters
    @property
    def grad_raster(self) -> Raster:
        """The gradient raster: gradient waveform samples land on it."""
        return Raster(float(self.default.grad_raster_time), 'gradient')

    @property
    def rf_raster(self) -> Raster:
        """The RF raster: RF waveform samples land on it."""
        return Raster(float(self.default.rf_raster_time), 'RF')

    @property
    def adc_raster(self) -> Raster:
        """The ADC raster: dwell times land on it."""
        return Raster(float(self.default.adc_raster_time), 'ADC')

    @property
    def block_raster(self) -> Raster:
        """The block-duration raster: every block duration must be a multiple."""
        return Raster(float(self.default.block_duration_raster), 'block')

    @property
    def adc_samples_divisor(self) -> int:
        """Required divisor of the ADC sample count (4 on Siemens)."""
        return int(self.default.adc_samples_divisor)

    # --------------------------------------------------------------------- reporting
    def describe(self) -> str:
        """Return a human-readable multi-line summary."""
        lines = [f'System {self.name!r}  B0 = {self.b0_T:.4g} T  gamma = {self.gamma:.6g} Hz/T']
        for regime, opts in self.regimes.items():
            lines.append(
                f'  {regime:<10} '
                f'max_grad = {self.convert(opts.max_grad, "Hz/m", "mT/m"):7.2f} mT/m'
                f'   max_slew = {self.convert(opts.max_slew, "Hz/m/s", "T/m/s"):7.2f} T/m/s'
            )
        lines.append(
            f'  rasters    grad {convert(self.grad_raster.dt, "s", "us"):.0f} us'
            f'   rf {convert(self.rf_raster.dt, "s", "us"):.0f} us'
            f'   adc {convert(self.adc_raster.dt, "s", "ns"):.0f} ns'
            f'   block {convert(self.block_raster.dt, "s", "us"):.0f} us'
        )
        if self.hardware is not None:
            lines.append(f'  hardware   {self.source or "attached"}')
        return '\n'.join(lines)

    def params(self) -> dict[str, Any]:
        """Return a JSON-safe description, for the provenance sidecar."""
        return {
            'name': self.name,
            'b0_T': self.b0_T,
            'gamma_Hz_per_T': self.gamma,
            'rasters_s': {
                'grad': self.grad_raster.dt,
                'rf': self.rf_raster.dt,
                'adc': self.adc_raster.dt,
                'block': self.block_raster.dt,
            },
            'adc_samples_divisor': self.adc_samples_divisor,
            'regimes': {
                regime: {
                    'max_grad_mT_m': self.convert(opts.max_grad, 'Hz/m', 'mT/m'),
                    'max_slew_T_m_s': self.convert(opts.max_slew, 'Hz/m/s', 'T/m/s'),
                    'rf_dead_time_us': convert(opts.rf_dead_time, 's', 'us'),
                    'rf_ringdown_time_us': convert(opts.rf_ringdown_time, 's', 'us'),
                    'adc_dead_time_us': convert(opts.adc_dead_time, 's', 'us'),
                }
                for regime, opts in self.regimes.items()
            },
            'hardware_source': self.source,
            'fingerprint': self.fingerprint,
        }

    @property
    def fingerprint(self) -> str:
        """Short stable hash of every regime and the hardware source."""
        return self._fingerprint or self._compute_fingerprint()

    def _compute_fingerprint(self) -> str:
        h = hashlib.sha256()
        for regime in sorted(self.regimes):
            opts = self.regimes[regime]
            h.update(regime.encode())
            for attr in sorted(vars(opts)):
                h.update(f'{attr}={getattr(opts, attr)!r};'.encode())
        h.update((self.source or '').encode())
        return h.hexdigest()[:12]


# --------------------------------------------------------------------------- hardware
def synthetic_hardware(name: str = 'synthetic_generic') -> SimpleNamespace:
    """
    Return a vendor-free PNS hardware model for tests and CI.

    Shaped like the output of ``pypulseq.utils.siemens.asc_to_hw`` so it can be handed
    straight to ``Sequence.calculate_pns``, but the coefficients are the illustrative
    values from pypulseq's own ``safe_pns_prediction.safe_example_hw()`` reference
    implementation, not measurements from any scanner.  **It is not a real scanner and
    must never be used to clear a sequence for human scanning** -- use
    :func:`load_hardware` with the site's own ``.asc`` for that.

    Examples
    --------
    >>> hw = synthetic_hardware()
    >>> hw.name
    'synthetic_generic'
    >>> hw.x.stim_limit
    30.0
    """

    def axis(stim_limit: float, stim_thresh: float, g_scale: float) -> SimpleNamespace:
        return SimpleNamespace(
            tau1=0.20,
            tau2=0.03,
            tau3=3.0,
            a1=0.4,
            a2=0.10,
            a3=0.50,
            stim_limit=stim_limit,
            stim_thresh=stim_thresh,
            g_scale=g_scale,
        )

    return SimpleNamespace(
        name=name,
        checkID=0,
        x=axis(30.0, 24.0, 0.4),
        y=axis(15.0, 12.0, 0.7),
        z=axis(25.0, 20.0, 0.3),
        acoustic_resonances=(),
        is_synthetic=True,
    )


def load_hardware(
    filename: str,
    *,
    cardiac_model: bool = False,
    directory: str | os.PathLike[str] | None = None,
) -> tuple[SimpleNamespace, tuple[dict[str, float], ...], str]:
    """
    Load a vendor ``.asc`` gradient descriptor from outside the repository.

    Parameters
    ----------
    filename
        Bare file name, e.g. ``'CimaX.asc'``.  A path containing directory separators is
        rejected: the whole point is that the location comes from the environment.
    cardiac_model
        Pass ``True`` for the CNS (cardiac) response model instead of PNS.
    directory
        Overrides ``$SEQCRAFT_ASC_DIR`` for this call.

    Returns
    -------
    hardware, resonances, source
        The hardware model, the forbidden acoustic-resonance bands, and a provenance
        string of the form ``'<filename> sha256:<12 hex>'``.  Only the file *name* and
        hash are recorded -- never the contents, which are vendor-confidential.

    Raises
    ------
    ConfigurationError
        If ``SEQCRAFT_ASC_DIR`` is unset, the file is missing, or `filename` is a path.

    Notes
    -----
    ``.asc`` files describe the gradient amplifier's PNS/CNS response and its forbidden
    acoustic-resonance frequencies.  They are site- and vendor-confidential, are excluded
    by ``.gitignore``, and are never bundled, copied into, or read from this repository.
    """
    if os.sep in filename or '/' in filename:
        msg = format_error(
            f'load_hardware() takes a bare file name, got a path: {filename!r}.',
            {'reason': 'the directory must come from $' + ASC_ENV_VAR},
        )
        raise ConfigurationError(msg)

    root = directory or os.environ.get(ASC_ENV_VAR)
    if not root:
        msg = format_error(
            f'no gradient hardware directory configured, cannot load {filename!r}.',
            {
                'expected': f'${ASC_ENV_VAR} pointing at a directory of vendor .asc files',
                'why': 'vendor .asc files carry proprietary PNS and acoustic-resonance data '
                'and are never stored in this repository',
            },
            [
                f'set {ASC_ENV_VAR} to the directory holding your .asc files',
                'or use seqcraft.System.preset(...) plus synthetic_hardware() for tests',
            ],
        )
        raise ConfigurationError(msg)

    path = Path(root) / filename
    if not path.is_file():
        msg = format_error(
            f'gradient hardware file not found: {filename!r}.',
            {'looked in': str(root)},
        )
        raise ConfigurationError(msg)

    from pypulseq.utils.siemens.asc_to_hw import (  # noqa: PLC0415
        asc_to_acoustic_resonances,
        asc_to_hw,
    )
    from pypulseq.utils.siemens.readasc import readasc  # noqa: PLC0415

    asc, _extra = readasc(str(path))
    hardware = asc_to_hw(asc, cardiac_model=cardiac_model)
    resonances = tuple(asc_to_acoustic_resonances(asc=asc))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()[:12]
    return hardware, resonances, f'{filename} sha256:{digest}'
