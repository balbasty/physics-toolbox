function ok = isboolean(X)
% Check that a variable is a valid boolean: scalar AND (numeric OR logical)
%
% FORMAT ok = isboolean(x)
    ok = (isnumeric(X) || islogical(X)) && isscalar(X);
end