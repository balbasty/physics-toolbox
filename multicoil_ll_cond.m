function llm = multicoil_ll_cond(x, s, rho, A, mask)
% Compute conditional log-likelihood
%
% FORMAT ll = multicoil_ll_cond(x, s, rho, A, mask)
%
% x    - (File)Array [Nx Ny Nz Nc] - Complex coil images
% s    - (File)Array [Nx Ny Nz Nc] - Complex log-sensitivity profiles
% rho  - (File)Array [Nx Ny Nz]    - Complex mean image
% A    -       Array [Nc Nc]       - Noise precision matrix
% mask -       Array [Nx Ny]       - K-space sampling mask 
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 5
    mask = 1;
end

Nx = size(x,1);
Ny = size(x,2);
Nz = size(x,3);
Nc = size(x,4);
Nvox = Nx*Ny*Nz;

% -------------------------------------------------------------------------
% Compute log-likelihood (conditional)
llm = 0;
for z=1:Nz

    % -----------------------------------------------------------------
    % Load one slice of the complete coil dataset
    xz = loadarray(x(:,:,z,:), @single);
    xz = reshape(xz, [], Nc);

    % -----------------------------------------------------------------
    % Load one slice of the (previous) mean
    rhoz = loadarray(rho(:,:,z,:), @single);
    rhoz = reshape(rhoz, [], 1);

    % -----------------------------------------------------------------
    % Load one slice of the complete sensitivity dataset + correct
    sz   = loadarray(s(:,:,z,:), @single);
    sz   = reshape(sz, [], Nc);
    sz   = exp(sz);
    rhoz = bsxfun(@times, rhoz, sz);

    % -----------------------------------------------------------------
    % If incomplete sampling: push+pull coil-specific means
    prhoz = multicoil_pushpullwrap(reshape(rhoz, [Nx Ny 1 Nc]), mask);
    prhoz = reshape(prhoz, [], Nc);

    llm = llm - sum(double(real(dot(rhoz,(prhoz-2*xz)*A,1))));

end
llm = 0.5 * Nvox * llm;