function [ind, n] = symIndices(k, side)
% FORMAT [ind, k] = utils.symIndices(n, ['n'])
% FORMAT [ind, n] = utils.symIndices(k, 'k')
% k - Length of the linearized and sparse symmetric matrix
% n - Side size of the corresponding square matrix
%
% Returns a converter between row/col and linear indices for sparse
% symmetric matrices.
%
% 1) If needed, finds n the size of the corresponding square matrix
% 2) Returns a n*n matrix which links subscripts (ex: (1,2) or (2,3)) to a
% corresponding linear index when the data is stored "sparsely" to save
% redunduncies.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

    if nargin == 2
        side = strcmpi(side, 'n');
    else
        side = true;
    end

    % 1) Compute n
    if ~side
        n = round(max(roots([1 1 -2*k])));
    else
        n = k;
        k = n*(n+1)/2;
    end
    
    % 2) Build the correspondance matrix
    ind = diag(1:n);
    l = n;
    for i1=1:n
        for i2=(i1+1):n
            l = l + 1;
            ind(i1, i2) = l;
            ind(i2, i1) = l;
        end
    end
    
    if side
        n = k;
    end
end