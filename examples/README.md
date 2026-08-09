# Examples

Each folder is a complete scan, start to finish: build it, check its physics, write the files,
simulate them, reconstruct, quantify. Run a folder from inside it — the notebooks read and write
`./seq/` and import `../lib`.

| | What it covers |
|---|---|
| [`01_getting_started.ipynb`](01_getting_started.ipynb) | Blocks, modules and `compile`; the overlap rules, the escape hatches, writing a file. Not a scan. |
| [`dti_spiral/`](dti_spiral/) | Single-shot spin-echo **spiral** DTI at 1.88 mm, and the two-echo field map its 67 ms readout cannot do without. |
| [`dti_epi/`](dti_epi/) | The **same diffusion encoding** through a ramp-sampled **EPI** train — single-shot and two-shot, partial Fourier 0.75. |

Both DTI folders use one diffusion encoding and differ only in the readout, which is the point of
having both: the ADC maps have to agree, and where they do not, that is the readout's error budget.

| Folder | Build | Simulate and reconstruct |
|---|---|---|
| `dti_spiral/` | [`01_build.ipynb`](dti_spiral/01_build.ipynb) | [`02_simulate_and_reconstruct.ipynb`](dti_spiral/02_simulate_and_reconstruct.ipynb) |
| `dti_epi/` | [`01_build.ipynb`](dti_epi/01_build.ipynb) | [`02_simulate_and_reconstruct.ipynb`](dti_epi/02_simulate_and_reconstruct.ipynb) |

Run `01_build.ipynb` first; the reconstruction notebook reads the `.seq` files it writes. Within a
folder, play the field map before the diffusion scan — it takes a few seconds and the diffusion
reconstruction needs it.

## What is deliberately duplicated

**Each folder builds its own field map.** The dual-echo GRE cell is therefore in both build
notebooks, and the two copies will drift. That is a real cost, taken so that a folder is runnable on
its own with nothing above it, and so each map is parameterised to its own geometry — the EPI one is
matched to the EPI matrix, and a mismatched grid is a moiré generator rather than a smooth error.

## Requirements

Building needs only `seqcraft`. Simulating and reconstructing need `MRzeroCore`, `torch` and
`sigpy`, which are **not** part of the package — see [`lib/README.md`](lib/README.md).

Outputs under each `seq/` are build products and are not tracked by git.
