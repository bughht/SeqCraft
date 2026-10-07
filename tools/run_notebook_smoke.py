"""Execute the build-only example notebooks in an isolated temporary directory."""

from __future__ import annotations

import logging
import os
import shutil
import tempfile
from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor
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
    Path('gre_3d/01_build.ipynb'),
    Path('mprage_2d/01_build.ipynb'),
    Path('mp2rage_2d/01_build.ipynb'),
    Path('se_2d/01_build.ipynb'),
    Path('fse_2d/01_build.ipynb'),
    Path('fse_2d/03_module_api.ipynb'),
    Path('gre_epi_2d/01_build.ipynb'),
    Path('se_epi_2d/01_build.ipynb'),
    Path('megre_2d/01_build.ipynb'),
    Path('bssfp_2d/01_build.ipynb'),
    Path('bssfp_3d/01_build.ipynb'),
    Path('gre_radial_2d/01_build.ipynb'),
    Path('fat_sat/01_build.ipynb'),
    Path('gre_spiral_2d/01_build.ipynb'),
    Path('se_spiral_2d/01_build.ipynb'),
    Path('dwi_se_epi_2d/01_build.ipynb'),
    Path('t2prep_gre_2d/01_build.ipynb'),
    Path('pc_gre_2d/01_build.ipynb'),
    Path('gre_2d/03_flow_comp.ipynb'),
    Path('gre_spiral_2d/03_flow_comp.ipynb'),
    Path('megre_2d/03_flow_comp.ipynb'),
)


#: Each worker owns a Jupyter kernel subprocess carrying its own NumPy and matplotlib, so the
#: ceiling here is runner CPU and memory rather than this process's orchestration -- which is also
#: why a thread pool is enough: every thread spends its time waiting on a kernel.  Measured
#: locally, 1/2/4/6 workers ran the set in 119/66/37/31 s, so four takes most of the available
#: gain while leaving headroom on a hosted runner.  Deliberately not `os.cpu_count()`.
_MAX_WORKERS = 4


def _worker_count(notebooks: Sequence[Path], cpus: int | None = None) -> int:
    """Workers to use: bounded by `_MAX_WORKERS`, the machine, and the work available."""
    available = cpus if cpus is not None else (os.cpu_count() or 1)
    return max(1, min(_MAX_WORKERS, available, len(notebooks)))


def run(
    notebooks: Sequence[Path],
    execute: Callable[[Path], None],
    max_workers: int,
) -> list[tuple[Path, Exception]]:
    """
    Execute every notebook, and return the ones that failed **in ``notebooks`` order**.

    Ordering the report by the allowlist rather than by completion is the whole reason this is a
    function: with several kernels running, the order failures arrive in is a property of the
    runner's CPU that week, and a gate whose output reshuffles is one nobody can diff.

    Unlike the serial version this does not stop at the first failure -- by the time one is
    raised the other workers are already mid-notebook, so stopping early would report a set
    that depends on timing.  Every scheduled notebook is allowed to finish and every failure is
    reported.
    """
    def one(relative: Path) -> tuple[Path, Exception | None]:
        logging.info('Executing %s', relative)
        try:
            execute(relative)
        except Exception as error:  # noqa: BLE001 - every failure is reported by main()
            # Deliberately not `logging.exception`: several kernels are running, so dumping a
            # traceback here interleaves it with other notebooks' progress lines.  The traceback
            # is printed once, in allowlist order, in the report at the end.
            logging.error('FAILED    %s', relative)  # noqa: TRY400
            return relative, error
        logging.info('Finished  %s', relative)
        return relative, None

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        results = dict(pool.map(one, notebooks))
    return [(n, error) for n in notebooks if (error := results.get(n)) is not None]


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
    """Copy examples to scratch space, execute the allowlist, and report every failure."""
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    source = _ROOT / 'examples'
    with tempfile.TemporaryDirectory(prefix='seqcraft-notebooks-') as scratch:
        examples = Path(scratch) / 'examples'
        # One shared copy, as before.  No notebook reads another's output and no two write the
        # same filename, so the only thing concurrent notebooks share is the `seq/` directory
        # their own cells create with `exist_ok=True`.
        shutil.copytree(source, examples)
        workers = _worker_count(_NOTEBOOKS)
        logging.info('Executing %d notebooks, %d at a time', len(_NOTEBOOKS), workers)
        failures = run(_NOTEBOOKS, lambda relative: _execute(examples / relative), workers)

    if failures:
        logging.error('')
        logging.error('%d of %d notebooks failed:', len(failures), len(_NOTEBOOKS))
        for relative, error in failures:
            # `str()` of an nbclient CellExecutionError is the failing cell and the kernel's own
            # traceback, which is the part a CI log is read for.
            logging.error('')
            logging.error('--- %s', relative)
            logging.error('%s', error)
        raise SystemExit(1)
    logging.info('All %d notebooks executed.', len(_NOTEBOOKS))


if __name__ == '__main__':
    main()
