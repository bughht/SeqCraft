function sequence = compile(root, opts, outputPath, options)
%COMPILE Compile a LogicBlock tree and return an official mr.Sequence.
arguments
    root (1, 1) seqcraft.LogicBlock
    opts (1, 1) struct
    outputPath (1, 1) string
    options.Definitions (1, 1) struct = struct()
    options.PythonExecutable (1, 1) string = ""
end

exchangePath = string(tempname) + ".lbtx.json";
cleanup = onCleanup(@() deleteIfPresent(exchangePath));
seqcraft.writeLBTX(root, opts, exchangePath, Definitions=options.Definitions);

[status, output] = runCompiler(exchangePath, outputPath, options.PythonExecutable);
if status ~= 0
    error("seqcraft:CompileFailed", ...
        "Python SeqCraft compiler failed with exit code %d:\n%s", status, strtrim(output));
end
if ~isfile(outputPath)
    error("seqcraft:MissingOutput", ...
        "Python SeqCraft compiler succeeded without writing '%s'.", outputPath);
end
if strlength(strtrim(output)) > 0
    warning("seqcraft:CompilerWarning", "%s", strtrim(output));
end

sequence = mr.Sequence(opts);
sequence.read(outputPath);
end
