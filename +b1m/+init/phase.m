function sens = phase(mean, coils, sens)
% FORMAT sens = b1m.init.phase(mean, coils, [sens])
%
% mean  - [Nx Ny Nz]    Complex mean image
% coils - [Nx Ny Nz Nc] Complex coil images
% sens  - [Nx Ny Nz Nc] Complex log-sensitivities
%
% Initialise sensitivities with the Von Mises mean of angle(rho - x).
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging


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
for z=1:lat(3) 
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    xz = loadarray(coils(:,:,z,:), @single);
    xz = reshape(xz, [], Nc);
    mask = isfinite(xz);
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
    sumsin = sumsin + sum(sinsz,1);
    sumcos = sumcos + sum(cossz,1);

end

% -------------------------------------------------------------------------
% Compute mean shift (Von Mises maximum likelihood)
shift = atan2(sumsin,sumcos);

% -------------------------------------------------------------------------
% Write
for n=1:size(sens,4)
    sens(:,:,:,n) = sens(:,:,:,n) - 1i*imag(sens(:,:,:,n)) + 1i*shift(n);
end