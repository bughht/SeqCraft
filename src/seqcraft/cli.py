"""Command-line adapter for validating and compiling LBTX documents."""

from __future__ import annotations

import argparse
import json
import re
import sys
import warnings
from pathlib import Path
from typing import Any

from .compiler import compile_sequence
from .compiler.verification import CompilerContractError
from .errors import SeqCraftError, SeqCraftWarning
from .exchange import ExchangeError, read_logicblock

EXIT_OK = 0
EXIT_INPUT = 2
EXIT_COMPILE = 3
EXIT_ADAPTER = 4


def _warning_category(message: str) -> str:
    vocabulary = {
        'same-axis gradient merge': 'merge',
        'vector-norm limit': 'norm',
        'no ADC after': 'orphan_label',
        'resampled onto the raster': 'resample',
        'snapped to the block raster': 'snap',
    }
    return next((category for phrase, category in vocabulary.items() if phrase in message), 'compiler')


def _diagnostic(
    severity: str,
    category: str,
    error_type: str,
    message: str,
    source_path: str | None = None,
) -> dict[str, Any]:
    out: dict[str, Any] = {
        'severity': severity,
        'category': category,
        'error_type': error_type,
        'message': message,
    }
    if source_path:
        out['source_path'] = source_path
    return out


def _source_path(message: str) -> str | None:
    """Extract the compiler's stable provenance spelling without parsing a traceback."""
    explicit = re.search(r'^\s*from\s*:\s*([^\n]+)', message, flags=re.MULTILINE)
    if explicit:
        return explicit.group(1).strip()
    warning_site = re.search(r':\s*([A-Za-z0-9_.-]+)(?:[+, (]|$)', message)
    return warning_site.group(1) if warning_site else None


def _result(command: str, ok: bool, diagnostics: list[dict[str, Any]], **fields: Any) -> dict[str, Any]:
    return {
        'schema': 'seqcraft.diagnostic',
        'version': '0.1',
        'command': command,
        'ok': ok,
        'diagnostics': diagnostics,
        **fields,
    }


def _emit(result: dict[str, Any], style: str) -> None:
    if style == 'json':
        print(json.dumps(result, ensure_ascii=False, separators=(',', ':')))
        return
    if result['ok']:
        output = f" -> {result['output']}" if 'output' in result else ''
        print(f"{result['command']} succeeded{output}")
    for diagnostic in result['diagnostics']:
        stream = sys.stderr if diagnostic['severity'] in {'warning', 'error'} else sys.stdout
        location = f" ({diagnostic['source_path']})" if diagnostic.get('source_path') else ''
        print(f"{diagnostic['severity']}: {diagnostic['message']}{location}", file=stream)


def _run(args: argparse.Namespace) -> tuple[int, dict[str, Any]]:
    command = args.command
    try:
        root, opts, definitions = read_logicblock(args.input)
    except (ExchangeError, OSError) as err:
        category = 'exchange' if isinstance(err, ExchangeError) else 'io'
        path = err.path if isinstance(err, ExchangeError) else str(args.input)
        result = _result(command, False, [
            _diagnostic('error', category, type(err).__name__, str(err), path)
        ])
        return EXIT_INPUT, result

    if command == 'validate-tree':
        return EXIT_OK, _result(command, True, [], input=str(Path(args.input)))

    diagnostics: list[dict[str, Any]] = []
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always', SeqCraftWarning)
            sequence = compile_sequence(root, opts, definitions=definitions)
        diagnostics.extend(
            _diagnostic('warning', _warning_category(str(item.message)),
                        type(item.message).__name__, str(item.message),
                        _source_path(str(item.message)))
            for item in caught
            if isinstance(item.message, SeqCraftWarning)
        )
        sequence.write(str(args.output))
    except (SeqCraftError, CompilerContractError) as err:
        diagnostics.append(_diagnostic(
            'error', 'compiler', type(err).__name__, str(err), _source_path(str(err)),
        ))
        return EXIT_COMPILE, _result(command, False, diagnostics)
    except OSError as err:
        diagnostics.append(_diagnostic('error', 'io', type(err).__name__, str(err), str(args.output)))
        return EXIT_ADAPTER, _result(command, False, diagnostics)
    except Exception as err:  # pragma: no cover - last-resort stable adapter boundary
        diagnostics.append(_diagnostic('error', 'adapter', type(err).__name__, str(err)))
        return EXIT_ADAPTER, _result(command, False, diagnostics)

    return EXIT_OK, _result(
        command, True, diagnostics, input=str(Path(args.input)), output=str(Path(args.output)),
        duration_s=float(sequence.duration()[0]), block_count=len(sequence.block_events),
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog='seqcraft', description='SeqCraft command-line tools')
    subparsers = parser.add_subparsers(dest='command', required=True)
    for name in ('validate-tree', 'compile-tree'):
        command = subparsers.add_parser(name)
        command.add_argument('input', type=Path, help='LBTX .lb.json input')
        command.add_argument('--diagnostics', choices=('human', 'json'), default='human')
        if name == 'compile-tree':
            command.add_argument('--output', type=Path, required=True, help='Pulseq .seq output')
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the SeqCraft CLI and return its documented exit code."""
    args = _parser().parse_args(argv)
    code, result = _run(args)
    _emit(result, args.diagnostics)
    return code


if __name__ == '__main__':  # pragma: no cover
    raise SystemExit(main())
