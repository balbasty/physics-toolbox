function rho = multicoil_mean_ml(x, s, A, rho, optim)
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

lat = [size(x,1), size(x,2), size(x,3)];
N   = size(x,4);
gpu_on = isa(A, 'gpuArray');
if gpu_on, loadarray = @loadarray_gpu;
else,      loadarray = @loadarray_cpu; end

% -------------------------------------------------------------------------
% Allocate output volume
if nargin < 4
    rho = zeros(lat, 'like', x);
end
if nargin < 5
    optim = [true true];
end
optim = logical(optim);
if ~any(optim)
    warning('Nothing to optimise');
    return
end

% -------------------------------------------------------------------------
% Process slice-wise to save memory
for z=1:lat(3)

    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    xz = loadarray(x(:,:,z,:), @single);
    xz = reshape(xz, [], N);

    % ---------------------------------------------------------------------
    % Load one slice of the complete sensitivity dataset + correct
    sz   = loadarray(s(:,:,z,:), @single);
    sz   = reshape(sz, [], N);
    sz   = single(exp(double(sz)));

    % ---------------------------------------------------------------------
    % Compute mean
    % rho = (s'*A*x) ./ (s'*A*s)
    Ab = double(sz * A);
    rhoz = dot(Ab, double(xz), 2);
    rhoz = rhoz ./ single(dot(Ab, double(sz), 2));
    clear Ab
    
    if ~all(optim)
         rho0 = loadarray(rho(:,:,z,:), @single);
         rho0 = reshape(rho0, [], 1);
        if optim(1)
            phi  = angle(rho0); clear rh0
            mag  = real(exp(-1i*phi).*rhoz);
            rhoz = mag .* exp(1i*phi);
        else
            mag  = abs(rho0); clear rho0
            phi  = real(-(1/1i)*log(mag./rhoz));
            rhoz = mag .* exp(1i*phi);
        end
        clear mag phi
    end
    
    % ---------------------------------------------------------------------
    % Write slice
    rho(:,:,z) = reshape(rhoz, lat(1), lat(2));

end