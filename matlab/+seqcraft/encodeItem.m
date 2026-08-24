function value = encodeItem(item)
%ENCODEITEM Convert a builder item to its JSON-ready discriminated structure.
if isa(item, "seqcraft.LogicBlock")
    value = struct("kind", "block", "block", item.toStruct());
    return
end
if ~isstruct(item) || ~isfield(item, "kind") || string(item.kind) ~= "event"
    error("seqcraft:InvalidItem", ...
        "A LogicBlock item must be a seqcraft.LogicBlock or an event from seqcraft.event().");
end
if string(item.event_type) == "rf" && isempty(item.payload.use)
    item.payload.use = string(missing);
end
value = item;
end
