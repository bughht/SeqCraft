# Examples

| | What it covers |
|---|---|
| [`01_getting_started.ipynb`](01_getting_started.ipynb) | Blocks, `Opts` and `compile`; the overlap rules, provenance, checking against physics, writing a file — and `sc.Module` at the end, once there is a reason for one. **Uses no modules**, which is the point: it is the demonstration that the compile path stands alone. |
| [`gre_2d/`](gre_2d/) | A complete spoiled 2D gradient echo, built three ways and then simulated and reconstructed. This is where most of `sc.modules` came from. |
| [`mprage_2d/`](mprage_2d/) | An inversion-prepared segmented GRE: where `IRPrep` and `GRE2D.time_to_center_line` came from, and where an inversion time gets placed to the microsecond. Defines `MPRAGE2D` in its own notebook. |
| [`mp2rage_2d/`](mp2rage_2d/) | One inversion, two trains, and the `SET` label that separates them — plus the ratio that cancels the receive field. Defines `MP2RAGE2D` in its own notebook. |

## `gre_2d/`

| | |
|---|---|
| [`01_build.ipynb`](gre_2d/01_build.ipynb) | The sequence, three times over: raw pypulseq events, the four leaf modules composed inline, and the same composition written as a module. Every arithmetic check, three sampling patterns, three `.seq` files. **Needs nothing but `seqcraft`.** |
| [`02_simulate_and_reconstruct.ipynb`](gre_2d/02_simulate_and_reconstruct.ipynb) | Those three files against a BrainWeb phantom — PD, T1, T2, T2′, D and a synthesised B0, all six simulated — with an eight-element receive ring, then one reconstruction across three samplings. **Needs `seqcraft[sim,recon]`** and a one-off ~19 MB phantom download. |

`01` deliberately does **not** import `GRE2DTR` or `GRE2D`. It writes them, because the third pass
is where a reader sees a working composition become something reusable, and a pass that imported
the answer would teach nothing.
[`tests/modules/test_notebook_matches_the_package.py`](../tests/modules/test_notebook_matches_the_package.py)
asserts that what the notebook writes and what the package ships produce identical events, so the
tutorial cannot drift from the library without CI noticing.

Every module in `sc.modules` was extracted from this pair. That is the rule the library is built
on: a module that cannot be extracted without altering the sequence is not a module, and one whose
extraction does not shorten the notebook is a wrapper.

## `mprage_2d/` and `mp2rage_2d/`

| | |
|---|---|
| [`mprage_2d/01_build.ipynb`](mprage_2d/01_build.ipynb) | The timing, which is where the risk is: TI runs from `IRPrep.time_to_center()` to `GRE2D.time_to_center_line()`, and both ends are module methods because a tree of events cannot know either. Segmentation, dummy shots, the two minima that stop a train overlapping its own inversion, then `MPRAGE2D` written out. Two `.seq` files, one per line ordering. **Needs nothing but `seqcraft`.** |
| [`mprage_2d/02_simulate_and_reconstruct.ipynb`](mprage_2d/02_simulate_and_reconstruct.ipynb) | The null point at `TI = T1·ln2`, measured through a real shot — the one check arithmetic cannot make — and what the line ordering costs. **Needs `seqcraft[sim,recon]`.** |
| [`mp2rage_2d/01_build.ipynb`](mp2rage_2d/01_build.ipynb) | Two trains at two inversion times, and the stateful-label trap: `SET` has to be emitted before *every* train, because setting it once leaves every later shot's first train wearing the previous shot's value. The broken version is built, and the compiler refuses it. **Needs nothing but `seqcraft`.** |
| [`mp2rage_2d/02_simulate_and_reconstruct.ipynb`](mp2rage_2d/02_simulate_and_reconstruct.ipynb) | Split by `SET` read back out of the file, then `UNI = Re(S₁·conj(S₂))/(|S₁|²+|S₂|²)` — simulated against two different receive arrays to show that the magnitudes move and `UNI` does not. **Needs `seqcraft[sim]`.** |

`MPRAGE2D` and `MP2RAGE2D` are **defined in those notebooks and do not ship**: one consumer each,
and a module with one consumer belongs where that consumer is.
[`tests/modules/test_mprage_notebooks.py`](../tests/modules/test_mprage_notebooks.py) runs the two
build notebooks and asserts against what they defined, so a tutorial that stays a tutorial still
cannot drift silently.

## Requirements

Building needs only `seqcraft`. Simulating and reconstructing need `MRzeroCore`, `torch` and
`sigpy` — `pip install "seqcraft[sim,recon]"`.

## `data/`

Phantoms, shared by every example and downloaded on first use. `gre_2d/02` fetches BrainWeb
subject 05 (~19 MB) from the MRzero-Core repository into `examples/data/`; a second example wanting
a phantom should look there before fetching its own copy.

Not tracked by git — it is somebody else's data, it is the same for everyone, and one download is
cheaper than carrying it in the history. Outputs under any `seq/` are build products and are not
tracked either.

## What used to be here

Two DTI acquisitions — a spiral and an EPI — sat in `_parked/` as the acceptance test for whatever
module set was written next, along with the `lib/` helpers they needed. That job is discharged:
[`docs/adr/003`](../docs/adr/003-scanner-and-module-reform.md) asked for a module set chosen by
making real scans work, and `gre_2d/` is the one that was written against.

They are in git history, and the physics they depended on — the b-value solve, the
variable-density spiral, the EPI ramp-sampling moment integral — is in [`salvage/`](../salvage/) as
plain functions, which is where a future spiral or EPI module should start from. `mr0_bridge.py`
went with them: `02` uses `mr0.Sequence.import_file` on the written `.seq`, which is a stronger
check than converting the tree, because it tests the file a scanner would actually play.
