function llm = multicoil_sense_ll_cond(x, s, rho, A, msk, dir)
% Compute conditional log-likelihood
%
% FORMAT ll = multicoil_ll_cond(x, s, rho, A, msk, dir)
%
% x   - (File)Array [Nx Ny Nz Nc (2)] - Complex coil images
% s   - (File)Array [Nx Ny Nz Nc (2)] - Complex log-sensitivity profiles
% rho - (File)Array [Nx Ny Nz  1 (2)] - Complex mean image
% A   -       Array [Nc Nc]           - Noise precision matrix
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Fudge factor to account for the true number of independent observations
Nm = sum(msk(:));
No = numel(msk);
ff = Nm/No;

% -------------------------------------------------------------------------
% Compute log-likelihood (conditional)
llm = 0;
for z=1:size(rho, 3) 

    % ---------------------------------------------------------------------
    % Load one slice of the mask
    if size(msk,3) > 1
        msk1 = msk(:,:,z);
    else
        msk1 = msk;
    end
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    if size(x, 5) == 2
        % Two real components
        x1 = single(x(:,:,z,:,:));
        x1 = x1(:,:,:,:,1) + 1i * x1(:,:,:,:,2);
    else
        % One complex volume
        x1 = single(x1(:,:,z,:,:));
    end
    if isa(A, 'gpuArray')
        x1 = gpuArray(x1);
    end

    % ---------------------------------------------------------------------
    % Load one slice of the mean
    if size(rho, 5) == 2
        % Two real components
        rho1 = single(rho(:,:,z,:,:));
        rho1 = rho1(:,:,:,:,1) + 1i*rho1(:,:,:,:,2);
    else
        % One complex volume
        rho1 = single(rho(:,:,z,:,:));
    end
    if isa(A, 'gpuArray')
        rho1 = gpuArray(rho1);
    end

    % ---------------------------------------------------------------------
    % Load one slice of the complete sensitivity dataset + correct
    if size(s, 5) == 2
        % Two real components
        s1 = double(s(:,:,z,:,:));
        s1 = single(exp(s1(:,:,:,:,1) + 1i*s1(:,:,:,:,2)));
    else
        % One complex volume
        s1 = single(exp(double(s(:,:,z,:,:))));
    end
    if isa(A, 'gpuArray')
        s1 = gpuArray(s1);
    end
    
    % ---------------------------------------------------------------------
    % Pull and preprocess
    rho1 = bsxfun(@times, rho1, s1); clear s1
    rho1 = multicoil_pullwrap(rho1, msk1, dir);
    rho1 = rho1 - x1; clear x1;
    
    % ---------------------------------------------------------------------
    % Compute gradient
    rho1 = reshape(rho1, [], size(x,4));
    llm = llm - 0.5 * ff * sum(double(real(dot(rho1 * A, rho1, 2))));

end