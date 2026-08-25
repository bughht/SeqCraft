% GRE_2D_MODULE Build the same small GRE 2D scan through a SeqCraft Module.

exampleDir = fileparts(mfilename("fullpath"));
seqcraftRoot = fileparts(fileparts(exampleDir));
addpath(fullfile(seqcraftRoot, "matlab"), "-begin");
addpath(exampleDir, "-begin");
if isempty(which("mr.opts"))
    error("seqcraft_examples:MissingPulseq", ...
        "Add the official MATLAB Pulseq matlab directory to the MATLAB path.");
end

if ~exist("outputPath", "var") || strlength(string(outputPath)) == 0
    outputPath = fullfile(tempdir, "seqcraft_gre_2d_module.seq");
end
if ~exist("pythonExecutable", "var") || ...
        strlength(strtrim(string(pythonExecutable))) == 0
    pythonExecutable = string(getenv("SEQCRAFT_PYTHON"));
end
if strlength(strtrim(string(pythonExecutable))) == 0
    pythonExecutable = string(input( ...
        "Python executable with SeqCraft installed: ", "s"));
end
if strlength(strtrim(string(pythonExecutable))) == 0
    error("seqcraft_examples:MissingPython", ...
        "Enter a Python executable or set SEQCRAFT_PYTHON before running the example.");
end
if ~exist("plotSequence", "var")
    plotSequence = true;
end

matrix = [16 8];
fovM = [0.22 0.22];
sliceThicknessM = 5e-3;

opts = mr.opts( ...
    "MaxGrad", 40, "GradUnit", "mT/m", ...
    "MaxSlew", 150, "SlewUnit", "T/m/s", ...
    "rfDeadTime", 100e-6, ...
    "rfRingdownTime", 30e-6, ...
    "adcDeadTime", 10e-6, ...
    "B0", 3);

greTR = seqcraft_examples.GRE2DTR( ...
    opts, ...
    Matrix=matrix, ...
    FOVM=fovM, ...
    SliceThicknessM=sliceThicknessM, ...
    FlipAngleRad=deg2rad(10), ...
    TRSeconds=20e-3, ...
    DwellSeconds=20e-6);

root = seqcraft.LogicBlock("gre_2d");
for lineIndex = 0:(matrix(2) - 1)
    root.add(lineIndex * greTR.trSeconds, greTR.build(lineIndex));
end

definitions = struct("FOV", [fovM sliceThicknessM], "Matrix", matrix);
seq = seqcraft.compile( ...
    root, opts, string(outputPath), ...
    Definitions=definitions, ...
    PythonExecutable=string(pythonExecutable));
[timingOk, timingReport] = seq.checkTiming();
assert(timingOk, strjoin(timingReport, newline));
if plotSequence
    seq.plot();
end

fprintf("Wrote %s\n", string(outputPath));
