# MATLAB examples

These examples use official MATLAB Pulseq events and compile through the same Python SeqCraft
compiler as the Python frontend. Add official MATLAB Pulseq, SeqCraft's MATLAB frontend, and this
directory to the MATLAB path:

```matlab
addpath("/path/to/pulseq/matlab")
addpath("/path/to/SeqCraft/matlab")
addpath("/path/to/SeqCraft/examples/matlab")
setenv("SEQCRAFT_PYTHON", "/path/to/environment/bin/python")
```

Run either script:

```matlab
run("/path/to/SeqCraft/examples/matlab/gre_2d_logicblock.m")
run("/path/to/SeqCraft/examples/matlab/gre_2d_module.m")
```

Both build the same small GRE 2D scan and return an official `mr.Sequence` in `seq`. They write to
separate files under `tempdir` by default. Set `outputPath` before running a script to keep the
`.seq` somewhere else; set `pythonExecutable` to override `SEQCRAFT_PYTHON` for that run.

`gre_2d_logicblock.m` is the minimal path: native `mr.make*` events are placed directly in nested
`seqcraft.LogicBlock` objects. `gre_2d_module.m` uses the example-only
`seqcraft_examples.GRE2D` class. Its constructor designs the same events once and its protected
`buildImplicit` method only assembles the tree. The class lives in the example package rather than
the production `seqcraft` package so it demonstrates the Module contract without expanding the
MATLAB frontend's public API.

This small amount of duplication is intentional: the two scripts are a side-by-side proof that
the compiler path stands alone and that wrapping a reusable composition in a Module does not
change its sequence semantics.
