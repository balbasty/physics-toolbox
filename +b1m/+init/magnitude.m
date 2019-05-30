function sens = magnitude(mean, coils, sens)
% FORMAT sens = b1m.init.magnitude(mean, coils, [sens])
%
% mean  - [Nx Ny Nz]    Complex mean image
% coils - [Nx Ny Nz Nc] Complex coil images
% sens  - [Nx Ny Nz Nc] Complex log-sensitivities [zeros]
%
% Initialise sensitivities with the mean of log(rho - x).
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
sumlog = zeros(1,Nc);
Nvox   = zeros(1,Nc);
for z=1:lat(3) 
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    xz = loadarray(coils(:,:,z,:), @single);
    xz = reshape(xz, [], Nc);
    mask = isfinite(xz);
    xz = log(abs(xz)+eps);

    % ---------------------------------------------------------------------
    % Load one slice of the mean
    rz = loadarray(mean(:,:,z,:), @single);
    rz = reshape(rz, [], 1);
    rz = log(abs(rz)+eps);

    % ---------------------------------------------------------------------
    % Accumulate
    sz = bsxfun(@minus, xz, rz);
    sz(~mask) = 0;
    sumlog = sumlog + sum(sz,1);
    Nvox   = Nvox   + sum(mask,1);
end

% -------------------------------------------------------------------------
% Compute mean shift
sumlog = sumlog./Nvox;

% -------------------------------------------------------------------------
% Write
for n=1:size(sens,4)
    sens(:,:,:,n) = sens(:,:,:,n) - real(sens(:,:,:,n)) + sumlog(n);
end