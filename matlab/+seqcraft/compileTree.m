function result = compileTree(inputPath, outputPath, options)
%COMPILETREE Compile one LBTX file with Python SeqCraft and return diagnostics.
arguments
    inputPath (1, 1) string
    outputPath (1, 1) string
    options.PythonExecutable (1, 1) string = ""
end
result = seqcraft.runCli("compile-tree", inputPath, OutputPath=outputPath, ...
    PythonExecutable=options.PythonExecutable);
end
