function [C,A] = multicoil_cov(rho, x, s, gpu)
% Maximum likelihood covariance given a set of observed coil images, 
% log-sensitivity profiles and a mean image.
%
% FORMAT [C,A] = multicoil_cov(rho, x, s, (gpu))
%
% rho - (File)Array [Nx Ny Nz  1 (2)] - Complex mean image
% x   - (File)Array [Nx Ny Nz Nc (2)] - Complex coil images
% s   - (File)Array [Nx Ny Nz Nc (2)] - Complex log-sensitivity profiles
% gpu                                 - If 'gpu', compute on GPU
%
% C   -       Array [Nc Nc]           - Noise covariance matrix
% A   -       Array [Nc Nc]           - Noise precision matrix
%
% >> C = sum((b*rho-x)*(b*rho-x)')/(2J)     [where ' = conjugate transpose]
%
% Nc = number of coils
% Images can either be complex or have two real components that are then 
% assumed to be the real and imaginary parts.
% An output FileArray can be provided by using `rho` as an input. If not 
% provided, the output volume will have the same format as the input coil 
% volume.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 4
    gpu = '';
end

% -------------------------------------------------------------------------
% Allocate output array
C = zeros(size(x,4), 'double');
if strcmpi(gpu, 'gpu')
    C = gpuArray(C);
end

% -------------------------------------------------------------------------
% Process slice-wise to save memory
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
    if isa(C, 'gpuArray')
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
    if isa(C, 'gpuArray')
        rho1 = gpuArray(rho1);
    end
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete sensitivity dataset + correct
    if size(s, 5) == 2
        % Two real components
        s1 = reshape(double(s(:,:,z,:,:)), [], size(x,4), 2);
        s1 = single(exp(s1(:,:,1) + 1i*s1(:,:,2)));
    else
        % One complex volume
        s1 = reshape(single(exp(double(s(:,:,z,:,:)))), [], size(x,4));
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
    for n=1:size(x,4)
        C(n,n) = C(n,n) + sum(double(real( rho1(:,n) .* conj(rho1(:,n)) )));
        for m=n+1:size(x,4)
            tmp = sum(double(real( rho1(:,n) .* conj(rho1(:,m)) )));
            C(n,m) = C(n,m) + tmp;
            C(m,n) = C(m,n) + tmp;
            clear tmp
        end
    end
    clear rho1

end
C = C/(2*size(x,1)*size(x,2)*size(x,3));
A = spm_matcomp('inv', C);
C = single(C);
A = single(A);
