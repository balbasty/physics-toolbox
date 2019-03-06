function ok = isboolean(X)
    ok = (isnumeric(X) || islogical(X)) && isscalar(X);
end