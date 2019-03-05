function mean = mean_fullysampled(coils, sens, prec, mean, optim)
% Maximum likelihood mean given a set of observed coil images, 
% log-sensitivity profiles and a noise precision (= inverse covariance) 
% matrix.
%
% FORMAT mean = multicoil_mean_ml(coils, sens, prec, [mean], ...)
%
% coils - (File)Array [Nx Ny Nz Nc] - Complex coil images
% sens  - (File)Array [Nx Ny Nz Nc] - Complex log-sensitivity profiles
% prec  -       Array [Nc Nc]       - Noise precision matrix
% mean  - (File)Array [Nx Ny Nz]    - Complex mean image [allocate]
%
% >> r = (s'*A*x) ./ (s'*A*s)      [where ' = conjugate transpose]
%
% Nc = number of coils
% An output FileArray can be provided by using `rho` as an input. If not 
% provided, the output volume will have the same format as the input coil 
% volume.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

lat = [size(coils,1), size(coils,2), size(coils,3)];
NC  = size(coils,4);
gpu_on = isa(prec, 'gpuArray');
if gpu_on, loadarray = @utils.loadarray_gpu;
else,      loadarray = @utils.loadarray_cpu; end

% -------------------------------------------------------------------------
% Allocate output volume
if nargin < 4 || isempty(mean)
    mean = zeros(lat, 'like', coils);
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
    xz = loadarray(coils(:,:,z,:), @single);
    xz = reshape(xz, [], NC);

    % ---------------------------------------------------------------------
    % Load one slice of the complete sensitivity dataset + correct
    sz = loadarray(sens(:,:,z,:), @single);
    sz = reshape(sz, [], NC);
    sz = exp(sz);

    % ---------------------------------------------------------------------
    % Compute mean
    % rho = (s'*A*x) ./ (s'*A*s)
    Asz = double(sz * prec);
    rz  = dot(Asz, double(xz), 2);
    rz  = rz ./ single(dot(Asz, double(sz), 2));
    clear Asz
    
    if ~all(optim)
         rz0 = loadarray(mean(:,:,z), @single);
         rz0 = reshape(rz0, [], 1);
        if optim(1)
            phi  = angle(rz0); clear rz0
            mag  = real(exp(-1i*phi).*rz);
            rz = mag .* exp(1i*phi);
        else
            mag  = abs(rz0); clear rz0
            phi  = real(-(1/1i)*log(mag./rz));
            rz = mag .* exp(1i*phi);
        end
        clear mag phi
    end
    
    % ---------------------------------------------------------------------
    % Write slice
    mean(:,:,z) = reshape(rz, lat(1), lat(2));

end