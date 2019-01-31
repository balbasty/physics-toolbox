function s = multicoil_init_magnitude(rho, x, s)
% FORMAT s = multicoil_init_magnitude(rho, x, (s))
%
% rho - [Nx Ny Nz 1] Complex mean image
% x   - [Nx Ny Nz Nc] Complex coil images
% s   - [Nx Ny Nz Nc] Complex log-sensitivities
%
% Initialise sensitivities with the mean of log(rho - x).

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
sumlog = zeros(1,N);
for z=1:lat(3) 
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    xz = loadarray(x(:,:,z,:), @double);
    xz = reshape(xz, [], N);
    xz = log(abs(xz)+eps);

    % ---------------------------------------------------------------------
    % Load one slice of the mean
    rhoz = loadarray(rho(:,:,z,:), @double);
    rhoz = reshape(rhoz, [], 1);
    rhoz = log(abs(rhoz)+eps);

    % ---------------------------------------------------------------------
    % Accumulate sin(diff) and cos(diff)
    sz = bsxfun(@minus, xz, rhoz);
    sumlog = sumlog + sum(sz,1);

end

% -------------------------------------------------------------------------
% Compute mean shift
sumlog = sumlog/prod(lat);

% ---------------------------------------------------------------------
% Write
for n=1:size(s,4)
    s(:,:,:,n) = s(:,:,:,n) - real(s(:,:,:,n)) + sumlog(n);
end