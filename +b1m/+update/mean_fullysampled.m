function meanim = mean_fullysampled(coils, sens, prec, meanim, optim)
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
% Number of observed voxels (handle missing data)
Nvoxsub = zeros(Nc,1);
for z=1:Nz
    xz = loadarray(coils(:,:,z,:), @single);
    xz = reshape(xz, [], Nc);
    Nvoxsub = Nvoxsub + reshape(sum(isfinite(xz),1), Nc, 1);
end

% -------------------------------------------------------------------------
% Allocate output volume
if nargin < 4 || isempty(meanim)
    meanim = zeros(lat, 'like', coils);
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
    % Compute map of missing data
    cz = utils.gmm.lib('obs2code', xz);
    code_list = unique(cz)';
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete sensitivity dataset + correct
    sz = loadarray(sens(:,:,z,:), @single);
    sz = reshape(sz, [], NC);
    sz = exp(sz);

    % ---------------------------------------------------------------------
    % Compute mean
    % rho = (s'*A*x) ./ (s'*A*s)
    rz  = zeros(size(xz,1), 'single');
    for code=code_list
        mask = cz == code;
        bin  = utils.gmm.lib('code2bin', code, Nc);
        if ~any(bin)
            continue
        end
        A = utils.invPD(prec);
        A = A(bin,bin);
        A = utils.invPD(A);
        
        Asz = double(sz(mask,bin) * A);
        rz(mask)  = dot(Asz, double(xz(mask,bin)), 2);
        rz(mask)  = rz(mask) ./ single(dot(Asz, double(sz(mask,bin)), 2));
        clear Asz
    end

    if ~all(optim)
         rz0 = loadarray(meanim(:,:,z), @single);
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
    meanim(:,:,z) = reshape(rz, lat(1), lat(2));

end