function [C,A] = noise(mean, coils, sens, gpu)
% Maximum likelihood covariance given a set of observed coil images, 
% log-sensitivity profiles and a mean image.
%
% FORMAT [C,A] = b1m.update.noise(mean, obs, sens, [gpu])
%
% mean  - (File)Array [Nx Ny Nz]    - Complex mean image
% coils - (File)Array [Nx Ny Nz Nc] - Complex coil images
% sens  - (File)Array [Nx Ny Nz Nc] - Complex log-sensitivities
% gpu                               - If 'gpu', compute on GPU
%
% cov   -       Array [Nc Nc]       - Noise covariance matrix
% prec  -       Array [Nc Nc]       - Noise precision matrix (= inv(cov))
%
% Nc = number of coils
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

warning('Must be updated to handle partial k-space')

if nargin < 4
    gpu = '';
end

% -------------------------------------------------------------------------
% Allocate output array
C = zeros(size(coils,4), 'double');
if strcmpi(gpu, 'gpu')
    C = gpuArray(C);
end

% -------------------------------------------------------------------------
% Process slice-wise to save memory
for z=1:size(mean, 3) 

    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    if size(coils, 5) == 2
        % Two real components
        x1 = reshape(single(coils(:,:,z,:,:)), [], size(coils,4), 2);
        x1 = x1(:,:,1) + 1i*x1(:,:,2);
    else
        % One complex volume
        x1 = reshape(single(coils(:,:,z,:)), [], size(coils,4));
    end
    if isa(C, 'gpuArray')
        x1 = gpuArray(x1);
    end

    % ---------------------------------------------------------------------
    % Load one slice of the mean
    if size(mean, 5) == 2
        % Two real components
        rho1 = reshape(single(mean(:,:,z,:,:)), [], 1, 2);
        rho1 = rho1(:,:,1) + 1i*rho1(:,:,2);
    else
        % One complex volume
        rho1 = reshape(single(mean(:,:,z,:,:)), [], 1);
    end
    if isa(C, 'gpuArray')
        rho1 = gpuArray(rho1);
    end
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete sensitivity dataset + correct
    if size(sens, 5) == 2
        % Two real components
        s1 = reshape(double(sens(:,:,z,:,:)), [], size(coils,4), 2);
        s1 = single(exp(s1(:,:,1) + 1i*s1(:,:,2)));
    else
        % One complex volume
        s1 = reshape(single(exp(double(sens(:,:,z,:,:)))), [], size(coils,4));
    end
    if isa(C, 'gpuArray')
        s1 = gpuArray(s1);
    end
    rho1 = bsxfun(@times, rho1, s1);
    clear s1

    % ---------------------------------------------------------------------
    % Compute and sum residuals
    rho1 = rho1 - x1;
    clear x1
    for n=1:size(coils,4)
        C(n,n) = C(n,n) + sum(double(real( rho1(:,n) .* conj(rho1(:,n)) )));
        for m=n+1:size(coils,4)
            tmp = sum(double(real( rho1(:,n) .* conj(rho1(:,m)) )));
            C(n,m) = C(n,m) + tmp;
            C(m,n) = C(m,n) + tmp;
            clear tmp
        end
    end
    clear rho1

end
C = C/(2*size(coils,1)*size(coils,2)*size(coils,3));
A = spm_matcomp('inv', C);
C = single(C);
A = single(A);
