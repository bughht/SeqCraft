function value = decimal(number)
%DECIMAL Encode one finite binary64 time as a round-trippable decimal string.
arguments
    number (1, 1) double
end
if ~isfinite(number)
    error("seqcraft:InvalidTime", "LBTX time values must be finite.");
end
value = sprintf("%.17g", number);
end
