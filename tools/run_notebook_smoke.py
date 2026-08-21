"""Execute the build-only example notebooks in an isolated temporary directory."""

from __future__ import annotations

import logging
import shutil
import tempfile
from pathlib import Path

import nbformat
from nbclient import NotebookClient

_ROOT = Path(__file__).resolve().parents[1]
# Build-only notebooks: seqcraft and matplotlib, and nothing else.  The simulation notebooks are
# deliberately absent -- they need MRzeroCore, torch and sigpy plus a phantom download, which is
# the lab/nightly tier rather than this one.  See docs/testing.md.
_NOTEBOOKS = (
    Path('01_getting_started.ipynb'),
    Path('gre_2d/01_build.ipynb'),
    Path('mprage_2d/01_build.ipynb'),
    Path('mp2rage_2d/01_build.ipynb'),
    Path('se_2d/01_build.ipynb'),
    Path('fse_2d/01_build.ipynb'),
    Path('gre_epi_2d/01_build.ipynb'),
    Path('se_epi_2d/01_build.ipynb'),
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
