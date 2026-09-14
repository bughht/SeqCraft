function output = encodeValue(input, normalizeNames)
%ENCODEVALUE Convert MATLAB values to finite, JSON-safe LBTX values.
if islogical(input)
    output = input;
elseif isnumeric(input)
    if any(~isfinite(real(input)), "all") || any(~isfinite(imag(input)), "all")
        error("seqcraft:InvalidNumber", "LBTX numeric values must be finite.");
    end
    if isreal(input)
        output = input;
    else
        output = struct("real", real(input), "imag", imag(input));
    end
elseif ischar(input)
    output = input;
elseif isstring(input)
    if any(ismissing(input), "all")
        error("seqcraft:InvalidString", "LBTX strings cannot contain missing values.");
    end
    if isscalar(input)
        output = char(input);
    else
        output = cellstr(input);
    end
elseif iscell(input)
    output = cell(size(input));
    for index = 1:numel(input)
        output{index} = encodeValue(input{index}, normalizeNames);
    end
elseif isstruct(input)
    if ~isscalar(input)
        output = cell(size(input));
        for index = 1:numel(input)
            output{index} = encodeValue(input(index), normalizeNames);
        end
        return
    end
    output = struct();
    names = fieldnames(input);
    for index = 1:numel(names)
        sourceName = names{index};
        targetName = sourceName;
        if normalizeNames
            targetName = camelToSnake(sourceName);
        end
        if isfield(output, targetName)
            error("seqcraft:FieldCollision", ...
                "Fields collide after LBTX name normalization: '%s'.", targetName);
        end
        output.(targetName) = encodeValue(input.(sourceName), normalizeNames);
    end
else
    error("seqcraft:UnsupportedValue", ...
        "LBTX cannot encode MATLAB value of type %s.", class(input));
end
end
