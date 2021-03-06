function sens = phase(mean, coils, sens, senslog)
% FORMAT sens = b1m.init.phase(mean, coils, [sens], [senslog])
%
% mean  - [Nx Ny Nz]    Complex mean image
% coils - [Nx Ny Nz Nc] Complex coil images
% sens  - [Nx Ny Nz Nc] Complex log-sensitivities
% senslog -               Log-encoding of sensitivities [false]
%
% Initialise sensitivities with the Von Mises mean of angle(rho - x).
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 4
    senslog = false;
end

lat = [size(coils,1) size(coils,2) size(coils,3)];
Nc  = size(coils,4);
gpu_on = isa(mean, 'gpuArray');
if gpu_on, loadarray = @utils.loadarray_gpu;
else,      loadarray = @utils.loadarray_cpu; end

% -------------------------------------------------------------------------
% Allocate output volume
if nargin < 3
    sens = [];
end
if isempty(sens)
    sens = zeros(size(coils), 'like', coils);
end

% -------------------------------------------------------------------------
% Process slice-wise to save memory
sumsin = zeros(1,Nc);
sumcos = zeros(1,Nc);
sumabs = zeros(1,Nc);
for z=1:lat(3) 
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    xz = loadarray(coils(:,:,z,:), @single);
    xz = reshape(xz, [], Nc);
    mask = isfinite(xz);
    mz = abs(xz);
    xz = angle(xz);

    % ---------------------------------------------------------------------
    % Load one slice of the mean
    rz = loadarray(mean(:,:,z,:), @single);
    rz = reshape(rz, [], 1);
    rz = angle(rz);

    % ---------------------------------------------------------------------
    % Accumulate sin(diff) and cos(diff)
    sz = bsxfun(@minus, xz, rz);
    sinsz        = sin(sz);
    sinsz(~mask) = 0;
    cossz        = cos(sz);
    cossz(~mask) = 0;
    sumsin = sumsin + sum(mz.*sinsz,1);
    sumcos = sumcos + sum(mz.*cossz,1);
    sumabs = sumabs + sum(mz,1);
end
sumsin = sumsin./sumabs;
sumcos = sumcos./sumabs;

% -------------------------------------------------------------------------
% Compute mean shift (Von Mises maximum likelihood)
shift = atan2(sumsin,sumcos);

% -------------------------------------------------------------------------
% Write
for n=1:size(sens,4)
    sens1 = sens(:,:,:,n);
    if ~senslog
        sens1 = log(sens1);
    end
    sens1(~isfinite(sens1)) = 1;
    sens1 = sens1 - 1i*imag(sens1) + 1i*shift(n);
    if ~senslog
        sens1 = exp(sens1);
    end
    sens(:,:,:,n) = sens1;
    clear sens1
end
