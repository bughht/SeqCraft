# MATLAB examples

These examples use official MATLAB Pulseq events and compile through the same Python SeqCraft
compiler as the Python frontend. Add official MATLAB Pulseq, SeqCraft's MATLAB frontend, and this
directory to the MATLAB path:

```matlab
addpath("/path/to/pulseq/matlab")
setenv("SEQCRAFT_PYTHON", "/path/to/environment/bin/python")
```

Run either script:

```matlab
run("/path/to/SeqCraft/examples/matlab/gre_2d_logicblock.m")
run("/path/to/SeqCraft/examples/matlab/gre_2d_module.m")
```

Both build the same small GRE 2D scan, return an official `mr.Sequence` in `seq`, and call
`seq.plot()` after the timing check. They write to separate files under `tempdir` by default. Set
`outputPath` before running a script to keep the `.seq` somewhere else; set `pythonExecutable` to
override `SEQCRAFT_PYTHON` for that run. For automated or headless execution, set
`plotSequence = false` before running the script.

Each script finds this checkout's `matlab/` and `examples/matlab/` directories relative to its own
file and adds them for the current MATLAB session. Official MATLAB Pulseq remains an external
dependency whose path the user supplies; the scripts raise a focused error if `mr.opts` is absent.
They never call `savepath`.

`gre_2d_logicblock.m` is the minimal path: native `mr.make*` events are placed directly in nested
`seqcraft.LogicBlock` objects. `gre_2d_module.m` uses the example-only
`seqcraft_examples.GRE2DTR` class. As in Python, `GRE2DTR.build(lineIndex)` returns exactly one
repetition; the script owns the acquisition loop. Its constructor designs the events once and its
protected `buildImplicit` method only assembles the requested line. The class lives in the example
package rather than the production `seqcraft` package, so it demonstrates the Module contract
without expanding the MATLAB frontend's public API.

The example repetition includes a phase-encode rewinder, and its readout prephaser is calculated at
the actual centre ADC sample. Cross-validation checks that every repetition has zero net `y`
moment, that the centre samples land on the requested Cartesian `(kx, ky, kz)`, and that `LIN`
labels describe those lines.

This small amount of duplication is intentional: the two scripts are a side-by-side proof that
the compiler path stands alone and that wrapping a reusable composition in a Module does not
change its sequence semantics.
