"""Execute the build-only example notebooks in an isolated temporary directory."""

from __future__ import annotations

import logging
import shutil
import tempfile
from pathlib import Path

import nbformat
from nbclient import NotebookClient

_ROOT = Path(__file__).resolve().parents[1]
_NOTEBOOKS = (
    Path('01_getting_started.ipynb'),
    Path('dti_spiral/01_build.ipynb'),
    Path('dti_epi/01_build.ipynb'),
)


def _execute(path: Path) -> None:
    """Execute one notebook with its parent directory as the working directory."""
    notebook = nbformat.read(path, as_version=4)
    client = NotebookClient(
        notebook,
        timeout=900,
        kernel_name='python3',
        resources={'metadata': {'path': str(path.parent)}},
    )
    client.execute()


def main() -> None:
    """Copy examples to scratch space and fail on the first notebook cell error."""
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    source = _ROOT / 'examples'
    with tempfile.TemporaryDirectory(prefix='seqcraft-notebooks-') as scratch:
        examples = Path(scratch) / 'examples'
        shutil.copytree(source, examples)
        for relative in _NOTEBOOKS:
            logging.info('Executing %s', relative)
            _execute(examples / relative)


if __name__ == '__main__':
    main()
