# MATLAB interoperability preview

SeqCraft's MATLAB support is a frontend for LogicBlock Tree Exchange Format (LBTX) 0.1. MATLAB
builds or reads the same compiler-input tree as Python; the existing Python compiler remains the
only authority and writes the final Pulseq `.seq`.

## Setup

Install SeqCraft in a dedicated Python environment and add the repository MATLAB package to the
MATLAB path:

```matlab
addpath("/path/to/SeqCraft/matlab")
```

Choose the Python executable in one of two ways:

```matlab
setenv("SEQCRAFT_PYTHON", "/path/to/environment/bin/python")
% or pass PythonExecutable=... to validateTree/compileTree
```

On Windows use the environment's `python.exe`. Paths are shell-quoted by the adapter on Windows,
macOS, and Linux. The selected Python must be able to run `python -m seqcraft`; an editable install
or normal package installation both satisfy that.

## Build and compile a tree

Protocol scanner values use canonical SI/PyPulseq units: gradients in Hz/m, slew in Hz/m/s, RF in
Hz, and times in seconds.

```matlab
opts = seqcraft.scannerOpts(1703040, 6386400000, ...
    RFDeadTimeSeconds=100e-6, ...
    RFRingdownTimeSeconds=30e-6, ...
    ADCDeadTimeSeconds=10e-6, ...
    B0T=3);

gx = seqcraft.trapezoid("x", 100000, 100e-6, 800e-6, 100e-6);
adc = seqcraft.adc(64, 10e-6, DelaySeconds=100e-6, DeadTimeSeconds=10e-6);

readout = seqcraft.LogicBlock("readout");
readout = readout.add(0, gx, adc);       % value object: keep the returned value
root = seqcraft.LogicBlock("gre");
root = root.add(0, readout);

definitions = struct("FOV", [0.22 0.22 0.005]);
seqcraft.writeTree(root, opts, definitions, "gre.lb.json");
validation = seqcraft.validateTree("gre.lb.json");
compiled = seqcraft.compileTree("gre.lb.json", "gre.seq");
assert(validation.ok && compiled.ok)
```

`LogicBlock` preserves insertion order and relative starts; it never sorts nodes or legalizes
overlap. Assign the return from `add` because MATLAB classes in this frontend are value objects.

## Supported events

The MATLAB package provides `delay`, `trapezoid`, `arbitraryGradient`, `rf`, `adc`, `label`,
`trigger`, `output`, and `barrier`. Their arguments expose protocol physical fields directly.
Nested `LogicBlock` values cover the other item kind.

`readTree` returns `[root, opts, definitions, document]`. `writeTree` accepts the semantic first
three values and optional provenance. Unknown semantic fields and event types are rejected by the
Python validator; optional producer-specific data belongs in namespaced `extensions` at the JSON
protocol level.

## Diagnostics and failures

`validateTree` and `compileTree` always request JSON diagnostics. The returned struct contains:

- `ok` and `exit_code`;
- `diagnostics`, with severity, category, error type, message, and optional JSON source path;
- for a successful compile, `output`, `duration_s`, and `block_count`.

Exit code 2 means malformed/unsupported LBTX input, 3 means the compiler rejected a valid tree, and
4 means an adapter or file-system failure. Compiler warnings are returned as warning diagnostics;
they are not hidden or converted into success text.

If invocation fails before a JSON result is produced, MATLAB raises `seqcraft:CLIProtocol` and
includes the command status and captured output. Temporary files are caller-owned and retained on
failure so the exact compiler input remains available for diagnosis.

## Compatibility and testing

Version 0.1 rejects unknown semantic fields. Time values are decimal strings in seconds and survive
Python/MATLAB round trips without changing raster classification. Provenance does not affect the
semantic hash.

Run `matlab/tests/testLBTX.m` for MATLAB-local builder coverage and
`pytest -m crossval tests/exchange/test_matlab.py` for the bidirectional bridge. The normal Python
suite additionally checks CartesianLine echo/sample semantics and a Spiral arbitrary-gradient
stress case. See [ADR-005](adr/005-logicblock-tree-exchange.md) for the protocol decisions.
