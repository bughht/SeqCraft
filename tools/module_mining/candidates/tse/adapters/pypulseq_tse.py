"""
R2 -- the PyPulseq ``write_tse.py`` example.

Loaded from disk and called, **never edited**.  The script's own ``main(...)`` is the builder, so
the harness only has to put the package it imports on the path and pass the protocol arguments it
already accepts.

What it does not accept as an argument is the scanner: ``system = pp.Opts(...)`` is written in the
body.  :func:`overridden_opts` therefore patches the ``Opts`` factory the module resolves, for the
duration of one call, and restores it afterwards.  That is a harness intervention rather than a
source change -- the script still designs everything itself -- and it is what makes the ringdown
experiment possible without forking the reference.

Set ``MODULE_MINING_PYPULSEQ`` to the checkout to use; the default is the local fork recorded in
``docs/plans/module-mining/tse/reference_inventory.md``, which is **not** ``pulseq/pypulseq``.
"""

from __future__ import annotations

import contextlib
import importlib.util
import os
import sys
import warnings
from pathlib import Path
from typing import Any

from ....reference import ReferenceSequence

CHECKOUT = Path(os.environ.get(
    'MODULE_MINING_PYPULSEQ', '/Users/yiyund/Code_mgh/pypulseq-matlab-like'))
SCRIPT = CHECKOUT / 'examples' / 'scripts' / 'write_tse.py'

PROVENANCE = {
    'repo': 'm-a-x-i-m-z/pypulseq-matlab-like',
    'path': 'examples/scripts/write_tse.py',
    'license': 'MIT',
    'role': 'supporting',
    'caveat': 'a fork; the official pulseq/pypulseq example has not been diffed against it',
}


def load() -> Any:
    """Import the reference script by path, with its own package importable."""
    src = CHECKOUT / 'src'
    if src.is_dir() and str(src) not in sys.path:
        sys.path.insert(0, str(src))
    spec = importlib.util.spec_from_file_location('module_mining_write_tse', SCRIPT)
    if spec is None or spec.loader is None:                           # pragma: no cover
        msg = f'cannot load the reference script at {SCRIPT}'
        raise FileNotFoundError(msg)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@contextlib.contextmanager
def overridden_opts(module: Any, **overrides: Any):
    """Force scanner options the script hard-codes, for one call, then put its factory back."""
    if not overrides:
        yield
        return
    factory = module.pp.Opts

    def patched(*args: Any, **kwargs: Any):
        return factory(*args, **{**kwargs, **overrides})

    module.pp.Opts = patched
    try:
        yield
    finally:
        module.pp.Opts = factory


def build(*, opts_overrides: dict[str, Any] | None = None, **protocol: Any) -> ReferenceSequence:
    """Run the reference at `protocol`, returning what it built."""
    module = load()
    arguments = {'plot': False, 'write_seq': False, 'test_report': False, **protocol}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with overridden_opts(module, **(opts_overrides or {})):
            sequence = module.main(**arguments)

    return ReferenceSequence(
        name='pypulseq_tse',
        sequence=sequence,
        parameters=arguments | {'opts_overrides': dict(opts_overrides or {})},
        definitions=dict(getattr(sequence, 'definitions', {}) or {}),
        provenance=PROVENANCE | {'checkout': str(CHECKOUT)},
        # The script states no echo spacing, no effective TE and no echo sample: `te` is a
        # requested first echo spacing, and everything else is implied by the construction.  That
        # absence is itself a finding, so it is recorded rather than filled in.
        semantic={'te_requested_s': arguments.get('te'), 'claims': 'none beyond te and tr'},
    )
