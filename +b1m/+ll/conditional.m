function llm = conditional(coils, sens, mean, prec, mask)
% Compute conditional log-likelihood
%
% FORMAT ll = b1m.ll.conditional(coils, sens, mean, prep, [mask])
%
% coils - (File)Array [Nx Ny Nz Nc] - Complex coil images
% sens  - (File)Array [Nx Ny Nz Nc] - Complex log-sensitivity profiles
% mean  - (File)Array [Nx Ny Nz]    - Complex mean image
% prec  -       Array [Nc Nc]       - Noise precision matrix
% mask  -       Array [Nx Ny]       - K-space sampling mask 
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 5
    mask = 1;
end

Nx = size(coils,1);
Ny = size(coils,2);
Nz = size(coils,3);
Nc = size(coils,4);
Nvox = Nx*Ny*Nz;

% -------------------------------------------------------------------------
% Compute log-likelihood (conditional)
llm = 0;
for z=1:Nz

    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    xz = loadarray(coils(:,:,z,:), @single);
    xz = reshape(xz, [], Nc);

    % ---------------------------------------------------------------------
    % Load one slice of the (previous) mean
    rz = loadarray(mean(:,:,z,:), @single);
    rz = reshape(rz, [], 1);

    % ---------------------------------------------------------------------
    % Load one slice of the complete sensitivity dataset + correct
    sz = loadarray(sens(:,:,z,:), @single);
    sz = reshape(sz, [], Nc);
    sz = exp(sz);
    rz = bsxfun(@times, rz, sz);

    % ---------------------------------------------------------------------
    % If incomplete sampling: push+pull coil-specific means
    prz = b1m.adjoint_forward(reshape(rz, [Nx Ny 1 Nc]), mask);
    prz = reshape(prz, [], Nc);

    llm = llm - sum(double(real(dot(rz,(prz-2*xz)*prec,1))));

end
llm = 0.5 * Nvox * llm;