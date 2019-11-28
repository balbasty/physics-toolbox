function ok = isintval(X)
% Check that a variable is a valid integer array
%
% FORMAT ok = isintval(x)
    ok = isnumeric(X) && all(X == round(X));
end