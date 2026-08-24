function value = encodeBlock(block)
%ENCODEBLOCK Convert one LogicBlock into the LBTX-owned tree structure.
nodes = cell(1, numel(block.nodes));
for index = 1:numel(block.nodes)
    source = block.nodes{index};
    item = source.item;
    if isa(item, "seqcraft.LogicBlock")
        node = struct("start_s", source.start, "block", encodeBlock(item));
    elseif string(item.type) == "seqcraft_barrier"
        node = struct("start_s", source.start, "barrier", char(string(item.tag)));
    else
        node = struct("start_s", source.start, "event", encodeValue(item, true));
    end
    nodes{index} = node;
end
value = struct("tag", char(block.tag), "nodes", {nodes});
end
