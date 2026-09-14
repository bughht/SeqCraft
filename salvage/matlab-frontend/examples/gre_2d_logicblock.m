%% GRE 2D from native Pulseq events and LogicBlocks
% This example designs a small Cartesian GRE acquisition with official
% MATLAB Pulseq events, assembles the scan directly as nested LogicBlocks,
% and compiles it with the Python SeqCraft compiler. The companion
% gre_2d_module.m example builds the same sequence through GRE2DTR.
%
% Optional workspace inputs:
%   outputPath       - destination .seq path
%   pythonExecutable - Python executable with SeqCraft installed
%   plotSequence     - whether to call seq.plot() after compilation

%% Path setup
exampleDir = fileparts(mfilename("fullpath"));
seqcraftRoot = fileparts(fileparts(exampleDir));
addpath(fullfile(seqcraftRoot, "matlab"), "-begin");
addpath(exampleDir, "-begin");
if isempty(which("mr.opts"))
    error("seqcraft_examples:MissingPulseq", ...
        "Add the official MATLAB Pulseq matlab directory to the MATLAB path.");
end

%% Runtime and output configuration
if ~exist("outputPath", "var") || strlength(string(outputPath)) == 0
    sequenceDir = fullfile(exampleDir, "seq");
    if ~isfolder(sequenceDir)
        mkdir(sequenceDir);
    end
    outputPath = fullfile(sequenceDir, "gre_2d_logicblock.seq");
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

%% Sequence and scanner parameters
matrix = [16 8];
fovM = [0.22 0.22];
sliceThicknessM = 5e-3;
flipAngleRad = deg2rad(10);
trSeconds = 20e-3;
dwellSeconds = 20e-6;

opts = mr.opts( ...
    "MaxGrad", 40, "GradUnit", "mT/m", ...
    "MaxSlew", 150, "SlewUnit", "T/m/s", ...
    "rfDeadTime", 100e-6, ...
    "rfRingdownTime", 30e-6, ...
    "adcDeadTime", 10e-6, ...
    "B0", 3);

%% Excitation and readout event design
[rf, gz, gzRephaseMinimum] = mr.makeSincPulse( ...
    flipAngleRad, opts, ...
    "Duration", 1e-3, ...
    "SliceThickness", sliceThicknessM, ...
    "use", "excitation");

readoutArea = matrix(1) / fovM(1);
gx = mr.makeTrapezoid( ...
    "x", opts, "FlatArea", readoutArea, "FlatTime", matrix(1) * dwellSeconds);
adc = mr.makeAdc( ...
    matrix(1), opts, "Dwell", dwellSeconds, "Delay", gx.riseTime);

%% Coupled winder design
% All three simultaneous gradients use one raster-aligned duration and
% per-axis limits scaled so their vector sum respects the scanner limits.
winderMaxGrad = opts.maxGrad / sqrt(3);
winderMaxSlew = opts.maxSlew / sqrt(3);
gzRephaseConstrainedMinimum = mr.makeTrapezoid( ...
    "z", opts, "Area", gzRephaseMinimum.area, ...
    "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);
preEchoSample = floor(matrix(1) / 2);
flatBeforeEcho = adc.delay - gx.riseTime + (preEchoSample + 0.5) * dwellSeconds;
prephaserArea = -gx.amplitude * (0.5 * gx.riseTime + flatBeforeEcho);
gxPreMinimum = mr.makeTrapezoid( ...
    "x", opts, "Area", prephaserArea, ...
    "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);
centerLine = floor(matrix(2) / 2);
gyMaximum = mr.makeTrapezoid( ...
    "y", opts, "Area", centerLine / fovM(2), ...
    "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);
winderSeconds = max([ ...
    mr.calcDuration(gzRephaseConstrainedMinimum), ...
    mr.calcDuration(gxPreMinimum), ...
    mr.calcDuration(gyMaximum)]);
winderSeconds = ceil(winderSeconds / opts.gradRasterTime) * opts.gradRasterTime;

gzRephase = mr.makeTrapezoid( ...
    "z", opts, "Area", gzRephaseMinimum.area, "Duration", winderSeconds, ...
    "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);
gxPre = mr.makeTrapezoid( ...
    "x", opts, "Area", prephaserArea, "Duration", winderSeconds, ...
    "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);

%% Phase rewind and spoiler design
% The y rewind cancels each line's phase-encoding moment. The simultaneous
% x/y/z tail again uses vector-safe per-axis limits and a shared duration.
tailMaxGrad = opts.maxGrad / sqrt(3);
tailMaxSlew = opts.maxSlew / sqrt(3);
gyRewindMaximum = mr.makeTrapezoid( ...
    "y", opts, "Area", centerLine / fovM(2), ...
    "maxGrad", tailMaxGrad, "maxSlew", tailMaxSlew);
gxSpoilMinimum = mr.makeTrapezoid( ...
    "x", opts, "Area", 2 * readoutArea, ...
    "maxGrad", tailMaxGrad, "maxSlew", tailMaxSlew);
gzSpoilMinimum = mr.makeTrapezoid( ...
    "z", opts, "Area", 2 / sliceThicknessM, ...
    "maxGrad", tailMaxGrad, "maxSlew", tailMaxSlew);
tailSeconds = max([ ...
    mr.calcDuration(gyRewindMaximum), ...
    mr.calcDuration(gxSpoilMinimum), ...
    mr.calcDuration(gzSpoilMinimum)]);
tailSeconds = ceil(tailSeconds / opts.gradRasterTime) * opts.gradRasterTime;

%% Per-line phase encoding and labels
% Line indices are zero-based to match Pulseq LIN labels and GRE2DTR.build.
phaseEncodes = cell(1, matrix(2));
phaseRewinds = cell(1, matrix(2));
lineLabels = cell(1, matrix(2));
for lineIndex = 0:(matrix(2) - 1)
    phaseArea = (lineIndex - centerLine) / fovM(2);
    phaseEncodes{lineIndex + 1} = mr.makeTrapezoid( ...
        "y", opts, "Area", phaseArea, "Duration", winderSeconds, ...
        "maxGrad", winderMaxGrad, "maxSlew", winderMaxSlew);
    phaseRewinds{lineIndex + 1} = mr.makeTrapezoid( ...
        "y", opts, "Area", -phaseArea, "Duration", tailSeconds, ...
        "maxGrad", tailMaxGrad, "maxSlew", tailMaxSlew);
    lineLabels{lineIndex + 1} = mr.makeLabel("SET", "LIN", lineIndex);
end

gxSpoil = mr.makeTrapezoid( ...
    "x", opts, "Area", 2 * readoutArea, "Duration", tailSeconds, ...
    "maxGrad", tailMaxGrad, "maxSlew", tailMaxSlew);
gzSpoil = mr.makeTrapezoid( ...
    "z", opts, "Area", 2 / sliceThicknessM, "Duration", tailSeconds, ...
    "maxGrad", tailMaxGrad, "maxSlew", tailMaxSlew);

%% Repetition timing
excitationSeconds = mr.calcDuration(rf, gz);
readoutSeconds = mr.calcDuration(gx, adc);
readoutStart = excitationSeconds + winderSeconds;
tailStart = readoutStart + readoutSeconds;
contentSeconds = tailStart + tailSeconds;
assert(contentSeconds < trSeconds, "GRE line does not fit in the requested TR.");
trDelay = mr.makeDelay(trSeconds - contentSeconds);

%% LogicBlock acquisition loop
% Each child block is one TR; the root owns the full phase-encoding loop.
root = seqcraft.LogicBlock("gre_2d");
for lineIndex = 0:(matrix(2) - 1)
    line = seqcraft.LogicBlock("line_" + lineIndex);
    line.add(0, rf, gz);
    line.add(excitationSeconds, gzRephase, gxPre, phaseEncodes{lineIndex + 1});
    line.add(readoutStart, gx, adc, lineLabels{lineIndex + 1});
    line.add(tailStart, phaseRewinds{lineIndex + 1}, gxSpoil, gzSpoil);
    line.add(contentSeconds, trDelay);
    root.add(lineIndex * trSeconds, line);
end

%% Compile, validate, and inspect
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
