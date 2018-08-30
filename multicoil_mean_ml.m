function rho = multicoil_mean_ml(x, s, A, rho)
% Maximum likelihood mean given a set of observed coil images, 
% log-sensitivity profiles and a noise precision (= inverse covariance) 
% matrix.
%
% FORMAT rho = multicoil_mean_ml(x, s, A, (rho), ...)
%
% x   - (File)Array [Nx Ny Nz Nc (2)] - Complex coil images
% s   - (File)Array [Nx Ny Nz Nc (2)] - Complex log-sensitivity profiles
% A   -       Array [Nc Nc]           - Noise precision matrix
% rho - (File)Array [Nx Ny Nz  1 (2)] - Complex mean image [allocate]
%
% >> rho = (b'*A*b) ./ (b'*A*b)      [where ' = conjugate transpose]
%
% Nc = number of coils
% Images can either be complex or have two real components that are then 
% assumed to be the real and imaginary parts.
% An output FileArray can be provided by using `rho` as an input. If not 
% provided, the output volume will have the same format as the input coil 
% volume.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Allocate output volume
if nargin < 4
    rho = zeros(size(x,1), size(x,2), size(x,3), 1, size(x,5), 'single');
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
    if isa(A, 'gpuArray')
        x1 = gpuArray(x1);
    end

    % ---------------------------------------------------------------------
    % Load one slice of the complete bias dataset
    if size(s, 5) == 2
        % Two real components
        s1 = reshape(double(s(:,:,z,:,:)), [], size(x,4), 2);
        s1 = single(exp(s1(:,:,1) + 1i*s1(:,:,2)));
    else
        % One complex volume
        s1 = reshape(single(exp(double(s(:,:,z,:,:)))), [], size(x,4));
    end
    if isa(A, 'gpuArray')
        s1 = gpuArray(s1);
    end

    % ---------------------------------------------------------------------
    % Compute mean
    % rho = (s'*A*x) ./ (s'*A*s)
    Ab = double(s1 * A);
    rho1 = dot(Ab, double(x1), 2);
    rho1 = rho1 ./ single(dot(Ab, double(s1), 2));
    clear Ab
    rho1 = reshape(rho1, [size(rho,1) size(rho,2)]);
    
    % ---------------------------------------------------------------------
    % Write slice
    if size(rho, 5) == 2
        rho(:,:,z,1,1) = real(rho1);
        rho(:,:,z,1,2) = imag(rho1);
    else
        rho(:,:,z) = rho1;
    end

end