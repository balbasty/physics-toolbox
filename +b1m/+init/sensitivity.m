function sens = sensitivity(mean, coils, sens)
% FORMAT sens = b1m.init.sensitivity(mean, coils, [sens])
%
% mean    - [Nx Ny Nz]    Complex mean image
% coils   - [Nx Ny Nz Nc] Complex coil images
% sens    - [Nx Ny Nz Nc] Complex log-sensitivities [zeros]
%
% Initialise sensitivities with a maximum-likelihood constant value.
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
sumrr = zeros(1,Nc);
sumrx = zeros(1,Nc);
for z=1:lat(3) 
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    xz = loadarray(coils(:,:,z,:), @single);
    xz = reshape(xz, [], Nc);
    mask = isfinite(xz);

    % ---------------------------------------------------------------------
    % Load one slice of the mean
    rz = loadarray(mean(:,:,z,:), @single);
    rz = reshape(rz, [], 1);

    % ---------------------------------------------------------------------
    % Accumulate
    rx = bsxfun(@times, conj(rz), xz);
    rx(~mask) = 0;
    rx = reshape(rx, [], Nc);
    sumrx = sumrx + sum(rx, 1);
    clear rx
    
    rr = real(conj(rz).*rz);
    rr = bsxfun(@times, rr, mask);
    rr = reshape(rr, [], Nc);
    sumrr = sumrr + sum(rr, 1);
    clear rr
end

% -------------------------------------------------------------------------
% Compute ML constant sensitivity
sumrx = sumrx./sumrr;

% -------------------------------------------------------------------------
% Write 
for n=1:size(sens,4)
    sens(:,:,:,n) = sumrx(n);
end