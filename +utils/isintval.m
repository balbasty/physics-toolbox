function ok = isintval(X)
    ok = isnumeric(X) && all(X == round(X));
end