function tests = testLBTX
%TESTLBTX MATLAB-side value object and JSON writer tests.
tests = functiontests(localfunctions);
end

function testBuilderPreservesInsertionOrder(testCase)
root = seqcraft.LogicBlock("root");
root = root.add(1e-3, seqcraft.delay(2e-3));
root = root.add(0, seqcraft.barrier("first-in-file-second-in-time"));
value = root.toStruct();

verifyEqual(testCase, numel(value.nodes), 2);
verifyEqual(testCase, string(value.nodes{1}.start_s), "0.001");
verifyEqual(testCase, string(value.nodes{2}.start_s), "0");
end

function testReadWriteKeepsSemanticDocumentValid(testCase)
matlabRoot = fileparts(fileparts(mfilename("fullpath")));
repositoryRoot = fileparts(matlabRoot);
fixture = fullfile(repositoryRoot, "tests", "fixtures", "lbtx", "minimal.lb.json");
output = string(tempname) + ".lb.json";
cleanup = onCleanup(@() deleteIfPresent(output));

[root, opts, definitions] = seqcraft.readTree(fixture);
seqcraft.writeTree(root, opts, definitions, output);
decoded = jsondecode(fileread(output));

verifyEqual(testCase, string(decoded.schema), "seqcraft.logicblock");
verifyEqual(testCase, string(decoded.version), "0.1");
verifyEqual(testCase, numel(decoded.root.nodes), 2);
end

function testEveryEventHelperProducesTheProtocolDiscriminator(testCase)
events = { ...
    seqcraft.delay(1e-3), ...
    seqcraft.trapezoid("x", 1000, 10e-6, 20e-6, 10e-6), ...
    seqcraft.arbitraryGradient("y", [0 100 0], [5e-6 15e-6 25e-6]), ...
    seqcraft.rf([1+2i 3+4i], [0 1e-6], ShapeDurationSeconds=2e-6, CenterSeconds=1e-6), ...
    seqcraft.adc(4, 1e-6), ...
    seqcraft.label("LIN", "SET", 2), ...
    seqcraft.label("REP", "INC", 1), ...
    seqcraft.trigger("physio1", 20e-6), ...
    seqcraft.output("osc0", 20e-6), ...
    seqcraft.barrier("cut")};
types = cellfun(@(item) string(item.event_type), events);

verifyEqual(testCase, types, ["delay", "trap", "grad", "rf", "adc", "labelset", ...
    "labelinc", "trigger", "output", "seqcraft_barrier"]);
verifyEqual(testCase, events{4}.payload.signal_hz.imag, [2 4]);
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
