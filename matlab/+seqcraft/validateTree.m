function result = validateTree(inputPath, options)
%VALIDATETREE Validate one LBTX file through the authoritative Python reader.
arguments
    inputPath (1, 1) string
    options.PythonExecutable (1, 1) string = ""
end
result = seqcraft.runCli("validate-tree", inputPath, ...
    PythonExecutable=options.PythonExecutable);
end
