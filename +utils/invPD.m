function A = invPD(A)
% FORMAT iA = utils.invPD(A)
% A  - A positive-definite square matrix
% iA - Its inverse
%
% Stable inverse of a positive-definite matrix.
% Eigendecomposition is used to compute a more stable inverse.

% John Ashburner

    [V,D] = eig(A);
    if any(diag(D) <= 0)
        warning('Matrix has negative eigenvalues')
        D(D <= 0) = eps; % Threshold negative eigenvalues
    end
    D     = loaddiag(D);
    A     = real(V * (D \ V'));

end

function A = loaddiag(A)
% FORMAT A = loaddiag(A)
% A  - A square matrix
%
% Load A's diagonal until it is well conditioned for inversion.

    factor = 1e-7;
    while rcond(A) < 1e-5
        A = A + factor * max([diag(A); eps]) * eye(size(A));
        factor = 10 * factor;
    end

end