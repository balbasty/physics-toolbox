function ok = isboundary(X)
% Check that a variable is a valid boundary condition:
% 0, 'c', 'circulant' OR 1, 'n', 'neumann'
%
% FORMAT ok = isboundary(x)
    ok = (isnumeric(X) && isscalar(X) && 0 <= X && X <= 1) || ...
         (ischar(X)    && any(strcmpi(X, {'c','circulant','n','neumann'})));
end