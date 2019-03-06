function ok = isarray(X)
    ok = isnumeric(X) || islogical(X) || isa(X, 'file_array');
end