function value = shellQuote(text)
%SHELLQUOTE Quote one command argument for the current platform shell.
text = char(string(text));
if ispc
    quote = char(34);
    value = string([quote, strrep(text, quote, [quote, quote]), quote]);
else
    quote = char(39);
    replacement = [char(39), char(34), char(39), char(34), char(39)];
    value = string([quote, strrep(text, quote, replacement), quote]);
end
end
