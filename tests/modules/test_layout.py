"""
The folder taxonomy, asserted rather than described.

Two rules, and they are cheap to check while there are seven files -- which is exactly when a
taxonomy is worth pinning, because the first exception to it is always locally reasonable.

    The role folders hold ``sc.Module`` subclasses.  The top level holds what is not one.

and the flat re-export, which is what makes the taxonomy cheap to revise later: no caller writes
``seqcraft.modules.readout.cartesian_line``, so moving a file between folders breaks nothing.
"""

from __future__ import annotations

import importlib
import inspect
import pkgutil
from pathlib import Path

import pytest

import seqcraft as sc
import seqcraft.modules as modules

ROOT = Path(modules.__file__).parent

#: The folders whose files must each define a module subclass.
ROLE_FOLDERS = ('rf', 'encoding', 'readout', 'kernel', 'imaging')


def _defined_here(module) -> list[type]:
    """Classes a module defines itself, rather than imports."""
    return [
        value for value in vars(module).values()
        if inspect.isclass(value) and getattr(value, '__module__', '') == module.__name__
    ]


@pytest.mark.parametrize('folder', ROLE_FOLDERS)
def test_every_file_in_a_role_folder_defines_a_module(folder: str) -> None:
    """A role folder is for modules; anything else in one has been filed by habit."""
    for path in sorted((ROOT / folder).glob('*.py')):
        if path.name == '__init__.py':
            continue
        dotted = f'seqcraft.modules.{folder}.{path.stem}'
        defined = _defined_here(importlib.import_module(dotted))
        assert defined, f'{dotted} defines no class at all'
        assert all(issubclass(cls, sc.Module) for cls in defined), (
            f'{dotted} defines {[c.__name__ for c in defined]}, and {folder}/ is for '
            f'sc.Module subclasses'
        )


def test_nothing_at_the_top_of_modules_is_a_module_subclass() -> None:
    """
    The other half of the rule, and the one that caps what can accumulate up there.

    ``spoiler`` scores none of what a class buys -- it designs nothing per call, answers no
    timing question and holds no state -- so it is a function, and the rule needs no exception
    for it.
    """
    for path in sorted(ROOT.glob('*.py')):
        if path.name == '__init__.py':
            continue
        dotted = f'seqcraft.modules.{path.stem}'
        offenders = [
            cls.__name__ for cls in _defined_here(importlib.import_module(dotted))
            if issubclass(cls, sc.Module)
        ]
        assert not offenders, (
            f'{dotted} defines module subclass(es) {offenders}; those belong in a role folder'
        )


def test_every_public_name_is_reachable_flat() -> None:
    """
    Folders never appear in an import path.

    Every module subclass defined anywhere under ``modules/`` must be re-exported from
    ``sc.modules`` -- one added in a role folder and not re-exported is reachable only through
    the taxonomy, which is what the flat re-export exists to prevent.
    """
    exported = {getattr(modules, name) for name in modules.__all__}
    for info in pkgutil.walk_packages(modules.__path__, 'seqcraft.modules.'):
        if info.name.rsplit('.', 1)[-1].startswith('_'):
            continue
        for cls in _defined_here(importlib.import_module(info.name)):
            if issubclass(cls, sc.Module):
                assert cls in exported, f'{cls.__name__} is not re-exported from sc.modules'


def test_the_facade_reaches_it_without_a_second_import() -> None:
    """``import seqcraft`` is enough: ``modules`` is eager, unlike ``display``."""
    assert sc.modules is modules
    assert 'modules' in sc.__all__


def test_modules_do_not_reach_past_design() -> None:
    """
    A module needs pypulseq and the tree, and nothing downstream of them.

    The edge worth forbidding is ``modules -> compiler`` (or ``analysis``): a module that
    compiled internally to answer a question about itself would make every design call pay for a
    compile, and would make the library untestable against a tree it did not itself produce.
    """
    import ast

    banned = ('compiler', 'analysis', 'display', 'scanner')
    offenders = []
    for path in sorted(ROOT.rglob('*.py')):
        tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.level and node.module:
                head = node.module.split('.')[0]
                if head in banned:
                    offenders.append((path.relative_to(ROOT).as_posix(), node.module))
            elif isinstance(node, ast.Import):
                offenders.extend(
                    (path.relative_to(ROOT).as_posix(), alias.name)
                    for alias in node.names
                    if alias.name.startswith('seqcraft.')
                    and alias.name.split('.')[1] in banned
                )
    assert not offenders, f'modules/ reaches past design: {offenders}'
