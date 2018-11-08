function s = multicoil_init_phase(rho, x, s)
% FORMAT s = multicoil_init_phase(rho, x, (s))
%
% rho - [Nx Ny Nz 1] Complex mean image
% x   - [Nx Ny Nz Nc] Complex coil images
% s   - [Nx Ny Nz Nc] Complex log-sensitivities
%
% Initialise sensitivities with the Von Mises mean of angle(rho - x).

lat = [size(x,1) size(x,2) size(x,3)];
N   = size(x,4);
gpu_on = isa(rho, 'gpuArray');
if gpu_on, loadarray = @loadarray_gpu;
else,      loadarray = @loadarray_cpu; end

% -------------------------------------------------------------------------
% Allocate output volume
if nargin < 3
    s = [];
end
if isempty(s)
    s = zeros(size(x), 'like', x);
end

% -------------------------------------------------------------------------
% Process slice-wise to save memory
sumsin = zeros(1,N);
sumcos = zeros(1,N);
for z=1:lat(3) 
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    xz = loadarray(x(:,:,z,:), @single);
    xz = reshape(xz, [], N);
    xz = angle(xz);

    % ---------------------------------------------------------------------
    % Load one slice of the mean
    rhoz = loadarray(rho(:,:,z,:), @single);
    rhoz = reshape(rhoz, [], 1);
    rhoz = angle(rhoz);

    % ---------------------------------------------------------------------
    % Accumulate sin(diff) and cos(diff)
    sz = bsxfun(@minus, xz, rhoz);
    sumsin = sumsin + sum(sin(sz),1);
    sumcos = sumcos + sum(cos(sz),1);

end

% -------------------------------------------------------------------------
% Compute mean shift (Von Mises maximum likelihood)
shift = atan2(sumsin,sumcos);

% ---------------------------------------------------------------------
% Write
for n=1:size(s,4)
    s(:,:,:,n) = s(:,:,:,n) - 1i*imag(s(:,:,:,n)) + 1i*shift(n);
end