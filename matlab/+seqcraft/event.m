function value = event(eventType, payload)
%EVENT Construct one explicitly-discriminated LBTX event.
arguments
    eventType (1, 1) string
    payload (1, 1) struct
end
supported = ["delay", "trap", "grad", "rf", "adc", "labelset", "labelinc", ...
    "trigger", "output", "seqcraft_barrier"];
if ~any(eventType == supported)
    error("seqcraft:UnsupportedEvent", "Unsupported LBTX event type '%s'.", eventType);
end
value = struct("kind", "event", "event_type", eventType, "payload", payload);
end
