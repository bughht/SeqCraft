"""
The layout, asserted rather than described.

Five layers in one dependency order::

    errors  ->  design  ->  compiler  ->  analysis  ->  display

with :mod:`seqcraft.scanner` independent of all five, and the top-level facade allowed to reach
anywhere.

Both rules here were true at some point and quietly stopped being true, which is the whole
argument for testing them.  ``core`` once held the compiler *and* the geometry *and* the unit
table, so "the compile path is self-contained" was a claim about a directory rather than about
imports.  And ``CompiledSequence._verify`` took the compiler's private IR while living on the
result type, which made ``result -> compiler`` a real edge that no one had decided to add.

These are read off the **source**, not off ``sys.modules``: an import that only fires under
``TYPE_CHECKING`` still couples two modules for a reader, and an import that happens to be
transitively satisfied by import order is not evidence of anything.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

import seqcraft

ROOT = Path(seqcraft.__file__).parent

#: Later layers may import earlier ones, never the reverse.
#:
#: ``errors`` is a leaf and cross-cutting -- everything may raise -- so it comes first and may
#: import nothing else in the package.  ``analysis`` sits between the compiler and the display
#: because :func:`seqcraft.kspace` and :func:`seqcraft.pns` compile internally, and ``display``
#: draws what :func:`seqcraft.sample` returns.
ORDER = ['errors', 'design', 'compiler', 'analysis', 'display']

#: Nothing under ``compiler/`` may import these *at runtime*: they are beside the compile path.
OFF_THE_COMPILE_PATH = ['design.module', 'analysis', 'display', 'scanner']


def _imports(path: Path, *, runtime_only: bool = False) -> set[str]:
    """
    Return every intra-package module `path` imports, as dotted names from the package root.

    Relative imports are resolved against the importing module's own position, so
    ``from ..design.events import AXES`` inside ``compiler/emission.py`` comes back as
    ``design.events``.

    Parameters
    ----------
    runtime_only
        Skip imports inside ``if TYPE_CHECKING:``.  Those cost nothing at import time and exist to
        annotate a signature, which is a different claim from "this module runs that code".
    """
    tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    parts = path.relative_to(ROOT).with_suffix('').parts
    package = parts[:-1] if parts[-1] != '__init__' else parts[:-1]

    if runtime_only:
        for node in ast.walk(tree):
            for field, value in ast.iter_fields(node):
                if isinstance(value, list):
                    setattr(node, field, [
                        child for child in value
                        if not (isinstance(child, ast.If)
                                and ast.dump(child.test).find('TYPE_CHECKING') >= 0)
                    ])

    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.level:
            base = list(package[: len(package) - (node.level - 1)]) if node.level > 1 else list(package)
            target = [*base, *(node.module.split('.') if node.module else [])]
            out.add('.'.join(target))
            out.update('.'.join([*target, alias.name]) for alias in node.names)
        elif isinstance(node, ast.Import):
            out.update(
                name.name[len('seqcraft.'):]
                for name in node.names
                if name.name.startswith('seqcraft.')
            )
    return {name for name in out if name}


def _third_party_imports(path: Path) -> set[str]:
    """Return the top-level names of every absolute import in `path`."""
    tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out.update(alias.name.split('.')[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and not node.level and node.module:
            out.add(node.module.split('.')[0])
    return out


def _modules(package: str) -> list[Path]:
    """Every source file in one top-level package (or the single module of that name)."""
    directory = ROOT / package
    if directory.is_dir():
        return sorted(directory.rglob('*.py'))
    return [ROOT / f'{package}.py']


@pytest.mark.parametrize('index', range(len(ORDER)))
def test_no_package_imports_a_later_layer(index: int) -> None:
    """
    ``errors -> design -> compiler -> analysis -> display``, and never back up the chain.

    The edge this exists to forbid is ``compiler -> analysis``: the compiler must not reach for
    the measurements taken *of* what it produced.  It would be an easy one to add -- the
    self-check wants a moment per axis, and ``analysis.moments`` computes one -- and it would be
    wrong, because that moment walks the tree and the self-check needs the compiled side.  A
    self-check that shares code with the thing it checks compares a number with itself.
    """
    package = ORDER[index]
    later = set(ORDER[index + 1:])
    offenders = [
        (path.relative_to(ROOT).as_posix(), name)
        for path in _modules(package)
        for name in _imports(path)
        if name.split('.')[0] in later
    ]
    assert not offenders, (
        f'{package} imports a later layer: {offenders}. The order is {" -> ".join(ORDER)}.'
    )


def test_the_compile_path_imports_nothing_beside_it() -> None:
    """
    The compiler reads a tree and an ``Opts``, and nothing else it could do without.

    The component contract, the scanner package, the analysis toolbox and the display helpers are
    all either upstream or downstream of the transform.  A compile path that reached any of them
    could not be reused on its own, and the coupling would stay invisible until someone tried.

    One thing is deliberately **not** forbidden, and naming it is the point of writing the rule
    down rather than describing it: :func:`~seqcraft.compiler.legalization.check_limits` imports
    ``design.units`` to convert Hz/m to mT/m, so a limit violation is reported in the units the
    amplifier is specified in.  Banning it would either cost that message or duplicate the
    conversion factor, and a second copy of a factor is worse than an edge on the graph.

    ``design.geometry`` used to need an exemption here too, for the ``geometry=`` annotation.  It
    does not any more: the parameter is gone and the dataclass is in ``salvage/``.
    """
    offenders = [
        (path.relative_to(ROOT).as_posix(), name)
        for path in _modules('compiler')
        for name in _imports(path, runtime_only=True)
        if any(name == banned or name.startswith(f'{banned}.') for banned in OFF_THE_COMPILE_PATH)
    ]
    assert not offenders, f'the compile path reaches beside itself: {offenders}'


def test_display_is_the_only_module_that_imports_matplotlib() -> None:
    """
    ``import seqcraft`` stays cheap, and it stays cheap by there being one place to check.

    A fact about every file, rather than about whichever import happens to run first.
    """
    offenders = sorted(
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob('*.py')
        if path.name != 'display.py' and 'matplotlib' in _third_party_imports(path)
    )
    assert not offenders, f'matplotlib is imported outside display.py: {offenders}'


def test_importing_seqcraft_does_not_import_the_plotting_module() -> None:
    """
    The invariant the lazy ``_LAZY`` table exists to protect, checked in a fresh interpreter.

    Stated as "``seqcraft.display`` is not imported" rather than "matplotlib is not imported",
    because the second is not seqcraft's to promise: ``pypulseq/Sequence/calc_grad_spectrum.py``
    imports ``matplotlib`` at module level, so ``import pypulseq`` -- and therefore any import of
    seqcraft at all -- pulls it in regardless of what this package does.

    What *is* seqcraft's to promise is that nothing here reaches for it, which the previous test
    checks per file, and that the one heavyweight module stays behind ``__getattr__``.  A
    top-level ``from .display import plot_block`` added for convenience is exactly what would
    break this, and it would break it silently.

    ``seqcraft.testing`` used to be checked here too.  It is not deferred any more, it is gone --
    the two assertions it held are in ``tests/conftest.py``.
    """
    import subprocess
    import sys

    probe = (
        'import sys, seqcraft\n'
        'eager = sorted(m for m in ("seqcraft.display",) if m in sys.modules)\n'
        'print(eager)\n'
        'sys.exit(1 if eager else 0)\n'
    )
    result = subprocess.run(
        [sys.executable, '-c', probe], capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, (
        f'import seqcraft eagerly imported {result.stdout.strip()}; they belong in _LAZY'
    )
