function llm = multicoil_ll_cond(x, s, rho, A)
% Compute conditional log-likelihood
%
% FORMAT ll = multicoil_ll_cond(x, s, rho, A)
%
% x   - (File)Array [Nx Ny Nz Nc (2)] - Complex coil images
% s   - (File)Array [Nx Ny Nz Nc (2)] - Complex log-sensitivity profiles
% rho - (File)Array [Nx Ny Nz  1 (2)] - Complex mean image
% A   -       Array [Nc Nc]           - Noise precision matrix
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Compute log-likelihood (conditional)
llm = 0;
for z=1:size(rho, 3) 

    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    if size(x, 5) == 2
        % Two real components
        x1 = reshape(single(x(:,:,z,:,:)), [], size(x,4), 2);
        x1 = x1(:,:,1) + 1i*x1(:,:,2);
    else
        % One complex volume
        x1 = reshape(single(x(:,:,z,:)), [], size(x,4));
    end
    if isa(A, 'gpuArray')
        x1 = gpuArray(x1);
    end

    % ---------------------------------------------------------------------
    % Load one slice of the mean
    if size(rho, 5) == 2
        % Two real components
        rho1 = reshape(single(rho(:,:,z,:,:)), [], 1, 2);
        rho1 = rho1(:,:,1) + 1i*rho1(:,:,2);
    else
        % One complex volume
        rho1 = reshape(single(rho(:,:,z,:,:)), [], 1);
    end
    if isa(A, 'gpuArray')
        rho1 = gpuArray(rho1);
    end

    % ---------------------------------------------------------------------
    % Load one slice of the complete sensitivity dataset + correct
    if size(s, 5) == 2
        % Two real components
        s1 = reshape(single(s(:,:,z,:,:)), [], size(x,4), 2);
        s1 = s1(:,:,1) + 1i*s1(:,:,2);
    else
        % One complex volume
        s1 = reshape(single(s(:,:,z,:,:)), [], size(x,4));
    end
    s1 = single(exp(double(s1)));
    if isa(A, 'gpuArray')
        s1 = gpuArray(s1);
    end
    rho1 = bsxfun(@times, rho1, s1);
    clear s1

    % ---------------------------------------------------------------------
    % Compute gradient

    llm = llm - 0.5 * sum(double(real(dot((rho1 - x1) * A, rho1 - x1, 2))));

end