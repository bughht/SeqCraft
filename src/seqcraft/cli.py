"""Small process boundary used by the MATLAB frontend."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .compiler import compile_sequence
from .compiler.verification import CompilerContractError
from .errors import SeqCraftError
from .exchange import ExchangeError, read_logicblock

EXIT_INPUT = 2
EXIT_COMPILE = 3
EXIT_ADAPTER = 4


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog='seqcraft')
    subparsers = parser.add_subparsers(dest='command', required=True)
    command = subparsers.add_parser('compile-tree')
    command.add_argument('input', type=Path, help='LBTX .lbtx.json input')
    command.add_argument('--output', type=Path, required=True, help='Pulseq .seq output')
    return parser


def main(argv: list[str] | None = None) -> int:
    """Compile one LBTX file, printing only failures and ordinary Python warnings."""
    args = _parser().parse_args(argv)
    try:
        root, opts, definitions = read_logicblock(args.input)
    except (ExchangeError, OSError) as err:
        print(err, file=sys.stderr)
        return EXIT_INPUT

    try:
        sequence = compile_sequence(root, opts, definitions=definitions)
        sequence.write(str(args.output))
    except (SeqCraftError, CompilerContractError) as err:
        print(err, file=sys.stderr)
        return EXIT_COMPILE
    except OSError as err:
        print(err, file=sys.stderr)
        return EXIT_ADAPTER
    return 0


if __name__ == '__main__':  # pragma: no cover
    raise SystemExit(main())
