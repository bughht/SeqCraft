function result = runCli(command, inputPath, options)
%RUNCLI Invoke the Python LBTX adapter and decode its structured diagnostic.
arguments
    command (1, 1) string {mustBeMember(command, ["validate-tree", "compile-tree"])}
    inputPath (1, 1) string
    options.OutputPath (1, 1) string = ""
    options.PythonExecutable (1, 1) string = ""
end
python = options.PythonExecutable;
if strlength(python) == 0
    python = string(getenv("SEQCRAFT_PYTHON"));
end
if strlength(python) == 0
    python = "python";
end
parts = [seqcraft.shellQuote(python), "-m", "seqcraft", command, ...
    seqcraft.shellQuote(inputPath)];
if command == "compile-tree"
    if strlength(options.OutputPath) == 0
        error("seqcraft:MissingOutput", "OutputPath is required for compile-tree.");
    end
    parts = [parts, "--output", seqcraft.shellQuote(options.OutputPath)];
end
parts = [parts, "--diagnostics", "json"];
[status, output] = system(strjoin(parts, " "));
try
    result = jsondecode(strtrim(output));
catch exception
    error("seqcraft:CLIProtocol", ...
        "SeqCraft CLI returned non-JSON output (status %d): %s", status, output);
end
result.exit_code = status;
end
