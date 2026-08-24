# MATLAB frontend

SeqCraft's MATLAB frontend mirrors the Python Module and LogicBlock concepts while using official
MATLAB Pulseq scanner options, event structs, and final Sequence objects. Python remains the only
SeqCraft compiler.

## Setup

Add both official MATLAB Pulseq and the SeqCraft MATLAB package to the MATLAB path:

```matlab
addpath("/path/to/pulseq/matlab")
addpath("/path/to/SeqCraft/matlab")
```

Choose a Python environment containing SeqCraft either globally or per compile:

```matlab
setenv("SEQCRAFT_PYTHON", "/path/to/environment/bin/python")
% or: seqcraft.compile(..., PythonExecutable="/path/to/python")
```

The executable must support `python -m seqcraft compile-tree ...`. On Windows pass the environment's
`python.exe`.

## Build and compile

Use the same official Pulseq objects a normal MATLAB sequence uses:

```matlab
opts = mr.opts( ...
    "MaxGrad", 40, "GradUnit", "mT/m", ...
    "MaxSlew", 150, "SlewUnit", "T/m/s", ...
    "rfDeadTime", 100e-6, ...
    "rfRingdownTime", 30e-6, ...
    "adcDeadTime", 10e-6, ...
    "B0", 3);

gx = mr.makeTrapezoid("x", opts, "Area", 80, "Duration", 1e-3);
adc = mr.makeAdc(64, opts, "Dwell", 10e-6, "Delay", gx.riseTime);

readout = seqcraft.LogicBlock("readout");
readout.add(0, gx, adc);

root = seqcraft.LogicBlock("gre");
root.add(0, readout);

definitions = struct("FOV", [0.22 0.22 0.005]);
seq = seqcraft.compile(root, opts, "gre.seq", Definitions=definitions);
```

`compile` performs four operations: writes a temporary LBTX file, invokes the existing Python
compiler, writes the requested `.seq`, then loads it with `mr.Sequence(opts).read(...)`. The returned
value is an official `mr.Sequence`, so the normal Pulseq workflow continues unchanged:

```matlab
seq.plot();
[ok, report] = seq.checkTiming();
[ktraj_adc, t_adc] = seq.calculateKspacePP();
```

Any write, Python compile, output-file, or MATLAB read failure raises an exception. There is no
`result.ok` value that must be checked separately.

## LogicBlock semantics

`LogicBlock` is a handle class because Python blocks are mutable references too. `add` changes the
same object and returns it only to permit chaining:

```matlab
block = seqcraft.LogicBlock("readout");
alias = block;
block.add(0, gx, adc);
assert(numel(alias.nodes) == 2)
```

Nodes stay in insertion order and their `start` values are seconds relative to the enclosing block.
Overlap is legal; the Python compiler legalizes it. `duration` is measured from child events with
`mr.calcDuration`. `copy()` copies the block and node container but deliberately shares nested
blocks, matching Python's shallow-copy contract.

A SeqCraft barrier is the one non-Pulseq item:

```matlab
root.add(2e-3, seqcraft.barrier("readout-boundary"));
```

## Writing a Module

As in Python, a module constructor designs events once and the build step only assembles them. A
MATLAB subclass implements protected `buildImplicit`; users call the base `build` method:

```matlab
classdef Readout < seqcraft.Module
    properties
        gx
        adc
    end

    methods
        function obj = Readout(opts)
            obj@seqcraft.Module(opts);
            obj.gx = mr.makeTrapezoid("x", opts, "Area", 80, "Duration", 1e-3);
            obj.adc = mr.makeAdc(64, opts, "Dwell", 10e-6, ...
                "Delay", obj.gx.riseTime);
        end
    end

    methods (Access = protected)
        function block = buildImplicit(obj, varargin)
            block = seqcraft.LogicBlock();
            block.add(0, obj.gx, obj.adc);
        end
    end
end
```

```matlab
readout = Readout(opts);
block = readout.build();
```

The base `build` verifies the return type and gives an unnamed block the module tag or class name.
`buildImplicit` is the MATLAB implementation hook corresponding to the Python subclass `build`;
there is still only one module-build lifecycle.

## Examples

[`examples/matlab`](../examples/matlab/) contains the same small GRE 2D built in two forms:

- `gre_2d_logicblock.m` places native `mr.make*` events directly in nested `LogicBlock` objects;
- `gre_2d_module.m` uses an example-only `seqcraft_examples.GRE2D` whose constructor designs events
  and whose protected `buildImplicit` method only assembles them.

Both compile to an official `mr.Sequence`, use identical parameters and definitions, and are checked
for matching duration, block, ADC, label, and timing semantics. The example class is deliberately
outside `+seqcraft`, so it does not expand the MATLAB frontend's production API.

## Explicit exchange files

Most users should call `compile`. To inspect, version, or transfer compiler input explicitly:

```matlab
seqcraft.writeLBTX(root, opts, "gre.lbtx.json", Definitions=definitions);
```

LBTX stores official event public fields, normalized to snake_case, rather than defining SeqCraft
event builders or payload types. Version 1 is intentionally one-way: MATLAB writes and Python reads.
It does not promise a Python writer, MATLAB reader/round-trip, semantic hash, provenance, or a full
diagnostics protocol.

## Tests

With `PULSEQ_MATLAB_PATH` pointing to official MATLAB Pulseq:

```bash
pytest -m crossval tests/exchange/test_matlab.py
matlab -batch "addpath(getenv('PULSEQ_MATLAB_PATH')); addpath('matlab'); \
  results=runtests('matlab/tests'); assert(all([results.Passed]));"
```

The normal Python test tier validates the small schema, reader, compile command, and handwritten
fixture. See [ADR-005](adr/005-logicblock-tree-exchange.md) for ownership and deferred features.
