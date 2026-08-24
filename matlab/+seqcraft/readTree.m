function [root, opts, definitions, document] = readTree(path)
%READTREE Read an LBTX document into a MATLAB LogicBlock and semantic structs.
arguments
    path (1, 1) string
end
document = jsondecode(fileread(path));
if string(document.schema) ~= "seqcraft.logicblock" || string(document.version) ~= "0.1"
    error("seqcraft:UnsupportedVersion", "Expected seqcraft.logicblock version 0.1.");
end
root = seqcraft.LogicBlock.fromStruct(document.root);
opts = document.opts;
definitions = document.definitions;
end
