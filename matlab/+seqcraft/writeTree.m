function document = writeTree(root, opts, definitions, path, options)
%WRITETREE Write one MATLAB-built tree as an LBTX 0.1 JSON document.
arguments
    root (1, 1) seqcraft.LogicBlock
    opts (1, 1) struct
    definitions (1, 1) struct
    path (1, 1) string
    options.Provenance = []
end
document = struct( ...
    "schema", "seqcraft.logicblock", ...
    "version", "0.1", ...
    "units", "SI", ...
    "opts", opts, ...
    "definitions", definitions, ...
    "root", root.toStruct());
if ~isempty(options.Provenance)
    document.provenance = options.Provenance;
end
encoded = jsonencode(document, PrettyPrint=true);
file = fopen(path, "w", "n", "UTF-8");
if file < 0
    error("seqcraft:FileError", "Could not open '%s' for writing.", path);
end
cleanup = onCleanup(@() fclose(file));
fprintf(file, "%s\n", encoded);
end
