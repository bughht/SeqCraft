"""
Execute every Python block in ``docs/api_reference.md``, and check its index against the package.

Two failure modes this exists to prevent, both of which are how a reference rots:

**An example that no longer runs.**  Every fenced ``python`` block is executed in one shared
namespace, in document order, exactly as a reader would follow it.  A block that is illustrative
rather than runnable is skipped by name (see ``_SKIP``) and each skip is printed, so the exemptions
stay visible rather than accumulating quietly.

**A name that moved, or one that was added and never written down.**  §9 is a table of every public
name and the module it lives in.  It is compared against ``__all__`` in both directions: a name in
the package but not the table is undocumented, and a name in the table but not the package is a
lie.
"""

from __future__ import annotations

import argparse
import importlib
import logging
import pathlib
import re
import sys

_ROOT = pathlib.Path(__file__).resolve().parents[1]
_DOC = _ROOT / 'docs' / 'api_reference.md'

#: Fragments marking a block as illustrative rather than runnable, with the reason.
_SKIP = {
    'class LogicBlock:': 'a schema, not a statement',
    'class LegalizationResult:': 'a schema, not a statement',
    'class Node:': 'a schema, not a statement',
    'class Module(ABC):': 'a schema, not a statement',
    'class PlacedEvent:': 'a schema, not a statement',
    'class PulseqReadyBlock:': 'a schema, not a statement',
    'convert(value, from_unit': 'a signature, not a call',
    'from_scanner(manufacturer': 'a signature, not a call',
    'sc.opts.from_scanner(': 'needs the `systems` extra and a real scanner name',
    "load_hardware('CimaX.asc')": 'needs a vendor .asc via $SEQCRAFT_ASC_DIR',
    'sc.plot_block(': 'needs matplotlib, and would open a figure',
    'seq.plot()': "pypulseq's own plotter; would open a figure",
    'lb.nodes[2].start += 40e-6': 'a table of idioms, not a script',
    'gre(lines=range(ny))': 'a table of sampling patterns, not a script',
    "print(f'peak stimulation": 'a fragment of the surrounding example',
}

#: Modules whose ``__all__`` the index must account for, and the label §9 uses for each.
_INDEXED = {
    'seqcraft.errors': 'errors',
    'seqcraft.analysis': 'analysis',
    'seqcraft.display': 'display',
    'seqcraft.design.logic': 'design.logic',
    'seqcraft.design.module': 'design.module',
    'seqcraft.design.events': 'design.events',
    'seqcraft.design.timing': 'design.timing',
    'seqcraft.design.units': 'design.units',
    'seqcraft.scanner.opts': 'scanner.opts',
    'seqcraft.scanner.hardware': 'scanner.hardware',
    'seqcraft.modules': 'modules',
    'seqcraft.compiler': 'compiler',
    'seqcraft.compiler.errors': 'compiler.errors',
    'seqcraft.compiler.model': 'compiler.model',
    'seqcraft.compiler.placement': 'compiler.placement',
    'seqcraft.compiler.boundaries': 'compiler.boundaries',
    'seqcraft.compiler.legalization': 'compiler.legalization',
    'seqcraft.compiler.emission': 'compiler.emission',
    'seqcraft.compiler.verification': 'compiler.verification',
    'seqcraft._compat': '_compat',
}

#: Names the index deliberately omits: package re-exports of a whole submodule, and the two
#: type aliases that are documented in prose but are not functions or classes.
_NOT_INDEXED = {'events', 'logic', 'module', 'timing', 'units', 'hardware', 'opts'}


def _blocks(text: str) -> list[tuple[int, str]]:
    """Return every fenced ``python`` block as ``(line number, source)``."""
    out = []
    for match in re.finditer(r'^```python\n(.*?)^```', text, re.S | re.M):
        line = text.count('\n', 0, match.start()) + 1
        out.append((line, match.group(1)))
    return out


def _run_examples(text: str) -> list[str]:
    """Execute every runnable block in one shared namespace; return the failures."""
    namespace: dict[str, object] = {'__name__': '__api_reference__'}
    failures: list[str] = []
    for line, source in _blocks(text):
        skip = next((why for frag, why in _SKIP.items() if frag in source), None)
        if skip is not None:
            logging.info('  skipped block at line %-4d  (%s)', line, skip)
            continue
        try:
            exec(compile(source, f'{_DOC.name}:{line}', 'exec'), namespace)  # noqa: S102
        except Exception as err:  # noqa: BLE001 - every failure is worth reporting, not the first
            failures.append(f'{_DOC.name}:{line}: {type(err).__name__}: {err}')
    return failures


def _index(text: str) -> dict[str, set[str]]:
    """Parse §9 into ``name -> {module, ...}``.  A name may legitimately appear twice."""
    section = text.split('# 9. Index', 1)[-1]
    out: dict[str, set[str]] = {}
    for name, where in re.findall(r'^\| `([^`]+)` \| `([^`]+)` \|', section, re.M):
        out.setdefault(name, set()).add(where)
    return out


def _check_index(text: str) -> list[str]:
    """Compare §9 against every indexed module's ``__all__``, in both directions."""
    documented = _index(text)
    problems: list[str] = []
    actual: dict[str, set[str]] = {}

    for dotted, label in _INDEXED.items():
        module = importlib.import_module(dotted)
        for name in getattr(module, '__all__', []):
            if name in _NOT_INDEXED:
                continue
            actual.setdefault(name, set()).add(label)

    for name, wheres in sorted(actual.items()):
        if name not in documented:
            problems.append(f'{name} ({"/".join(sorted(wheres))}) is public but not in the index')
        elif not (wheres & documented[name]):
            problems.append(
                f'{name} is indexed under {sorted(documented[name])} '
                f'but lives in {sorted(wheres)}'
            )

    for name in sorted(documented):
        if name not in actual and name not in {'BARRIER', 'Event', 'Item', 'trapz'}:
            problems.append(f'{name} is in the index but is not exported by any module')
    return problems


def main() -> None:
    """Execute the reference's examples and check its index; exit non-zero on any failure."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quiet', action='store_true')
    args = parser.parse_args()
    logging.basicConfig(level=logging.WARNING if args.quiet else logging.INFO,
                        format='%(message)s')

    text = _DOC.read_text(encoding='utf-8')
    total = len(_blocks(text))

    failures = _run_examples(text)
    problems = _check_index(text)

    logging.info('%d python blocks in %s', total, _DOC.relative_to(_ROOT).as_posix())
    if failures:
        logging.error('\n%d example(s) failed:', len(failures))
        for f in failures:
            logging.error('  %s', f)
    if problems:
        logging.error('\n%d index problem(s):', len(problems))
        for p in problems:
            logging.error('  %s', p)

    if failures or problems:
        sys.exit(1)
    logging.info('every example runs, and the index matches the package')


if __name__ == '__main__':
    main()
