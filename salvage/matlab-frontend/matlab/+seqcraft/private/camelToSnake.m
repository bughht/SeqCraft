function output = camelToSnake(input)
%CAMELTOSNAKE Normalize one public MATLAB Pulseq field name for LBTX.
output = lower(regexprep(char(input), "([a-z0-9])([A-Z])", "$1_$2"));
end
