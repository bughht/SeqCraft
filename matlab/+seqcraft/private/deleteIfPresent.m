function deleteIfPresent(path)
%DELETEIFPRESENT Delete one temporary file if it exists.
if isfile(path)
    delete(path);
end
end
