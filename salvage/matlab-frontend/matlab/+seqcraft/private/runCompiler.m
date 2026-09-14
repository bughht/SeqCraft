function [status, output] = runCompiler(inputPath, outputPath, pythonExecutable)
%RUNCOMPILER Invoke the Python-only compiler adapter.
python = pythonExecutable;
if strlength(python) == 0
    python = string(getenv("SEQCRAFT_PYTHON"));
end
if strlength(python) == 0
    python = "python";
end
parts = [shellQuote(python), "-m", "seqcraft", "compile-tree", ...
    shellQuote(inputPath), "--output", shellQuote(outputPath)];
[status, output] = system(strjoin(parts, " ") + " 2>&1");
end
