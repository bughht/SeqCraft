function tests = testLBTX
%TESTLBTX MATLAB-side Module, LogicBlock, and native Pulseq writer tests.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testRoot = fileparts(mfilename("fullpath"));
addpath(fullfile(testRoot, "helpers"));
testCase.TestData.opts = mr.opts( ...
    "MaxGrad", 40, "GradUnit", "mT/m", ...
    "MaxSlew", 150, "SlewUnit", "T/m/s", ...
    "rfDeadTime", 100e-6, "rfRingdownTime", 30e-6, ...
    "adcDeadTime", 10e-6, "B0", 3);
end

function testAddUsesPythonLikeHandleSemantics(testCase)
root = seqcraft.LogicBlock("root");
alias = root;
root.add(1e-3, mr.makeDelay(2e-3));
root.add(0, seqcraft.barrier("first-in-file-second-in-time"));

verifyEqual(testCase, numel(alias.nodes), 2);
verifyEqual(testCase, alias.nodes{1}.start, 1e-3);
verifyEqual(testCase, alias.nodes{2}.start, 0);
verifyEqual(testCase, alias.duration, 3e-3, AbsTol=1e-15);
end

function testCopyOwnsNodesButSharesNestedBlocks(testCase)
nested = seqcraft.LogicBlock("nested");
original = seqcraft.LogicBlock("original");
original.add(0, nested);
copied = original.copy();

copied.add(1e-3, mr.makeDelay(1e-3));
verifyEqual(testCase, numel(original.nodes), 1);
verifyEqual(testCase, numel(copied.nodes), 2);

nested.add(0, mr.makeDelay(0.5e-3));
verifyEqual(testCase, original.nodes{1}.item.duration, 0.5e-3, AbsTol=1e-15);
verifyEqual(testCase, copied.nodes{1}.item.duration, 0.5e-3, AbsTol=1e-15);
end

function testModuleBuildImplicitIsFinalized(testCase)
module = TestReadoutModule(testCase.TestData.opts);
block = module.build();

verifyClass(testCase, block, "seqcraft.LogicBlock");
verifyEqual(testCase, block.tag, "TestReadoutModule");
verifyEqual(testCase, string(block.nodes{1}.item.type), "trap");
end

function testModuleRejectsNonLogicBlockOutput(testCase)
module = BadModule(testCase.TestData.opts);
verifyError(testCase, @() module.build(), "seqcraft:InvalidModuleOutput");
end

function testWriterAcceptsNativePulseqObjects(testCase)
opts = testCase.TestData.opts;
gx = mr.makeTrapezoid("x", opts, "Area", 80, "Duration", 1e-3);
adc = mr.makeAdc(64, opts, "Dwell", 10e-6, "Delay", gx.riseTime);
rf = mr.makeArbitraryRf([1+1i 1-1i], pi / 1200, opts, "use", "excitation");

readout = seqcraft.LogicBlock("readout");
readout.add(0, gx, adc);
root = seqcraft.LogicBlock("native-events");
root.add(0, readout);
root.add(2e-3, rf);
root.add(3e-3, seqcraft.barrier("end"));
root.add(3e-3, mr.makeDelay(10e-6));
root.add(3e-3, mr.makeLabel("SET", "LIN", 1));
root.add(3e-3, mr.makeTrigger("physio1", "Duration", 10e-6));
root.add(3e-3, mr.makeDigitalOutputPulse("osc0", "Duration", 10e-6));

output = string(tempname) + ".lbtx.json";
cleanup = onCleanup(@() deleteIfPresent(output));
seqcraft.writeLBTX(root, opts, output, Definitions=struct("FOV", [0.22 0.22 0.005]));
document = jsondecode(fileread(output));

verifyEqual(testCase, string(document.format), "seqcraft.logicblock");
verifyEqual(testCase, document.version, 1);
verifyEqual(testCase, document.opts.max_grad, opts.maxGrad);
verifyEqual(testCase, string(document.root.tag), "native-events");
encoded = fileread(output);
verifySubstring(testCase, encoded, '"rise_time"');
verifySubstring(testCase, encoded, '"signal"');
verifySubstring(testCase, encoded, '"real"');
verifySubstring(testCase, encoded, '"imag"');
verifySubstring(testCase, encoded, '"barrier"');
verifySubstring(testCase, encoded, '"labelset"');
verifySubstring(testCase, encoded, '"trigger"');
verifySubstring(testCase, encoded, '"output"');
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
