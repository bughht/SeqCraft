function writeLBTX(root, opts, path, options)
%WRITELBTX Write a MATLAB LogicBlock tree for the Python SeqCraft compiler.
arguments
    root (1, 1) seqcraft.LogicBlock
    opts (1, 1) struct
    path (1, 1) string
    options.Definitions (1, 1) struct = struct()
end

document = struct( ...
    "format", "seqcraft.logicblock", ...
    "version", 1, ...
    "pulseq", "1.5", ...
    "opts", encodeValue(opts, true), ...
    "definitions", encodeValue(options.Definitions, false), ...
    "root", encodeBlock(root));

encoded = jsonencode(document, PrettyPrint=true);
file = fopen(path, "w", "n", "UTF-8");
if file < 0
    error("seqcraft:FileError", "Could not open '%s' for writing.", path);
end
cleanup = onCleanup(@() fclose(file));
fprintf(file, "%s\n", encoded);
end
