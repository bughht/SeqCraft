function value = label(name, operation, labelValue)
%LABEL Construct an LBTX label set or increment event.
arguments
    name (1, 1) string
    operation (1, 1) string {mustBeMember(operation, ["SET", "INC"])}
    labelValue (1, 1) double {mustBeInteger}
end
eventType = "labelset";
if operation == "INC"
    eventType = "labelinc";
end
value = seqcraft.event(eventType, struct("label", name, "value", labelValue));
end
