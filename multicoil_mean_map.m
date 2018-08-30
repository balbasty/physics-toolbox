function rho = multicoil_mean_map(x, s, A, rho, prm, vs)
% Maximum a posteriori mean given a set of observed coil images, 
% log-sensitivity profiles and a noise precision (= inverse covariance) 
% matrix.
%
% FORMAT rho = multicoil_mean_map(x, s, A, (rho), ...)
%
% x   - (File)Array [Nx Ny Nz Nc (2)] - Complex coil images
% s   - (File)Array [Nx Ny Nz Nc (2)] - Complex log-sensitivity profiles
% A   -       Array [Nc Nc]           - Noise precision matrix
% rho - (File)Array [Nx Ny Nz  1 (2)] - Complex mean image
% prm -
% vs  -
%
% Nc = number of coils
% Images can either be complex or have two real components that are then 
% assumed to be the real and imaginary parts.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 6
    vs = [1 1 1];
end

% -------------------------------------------------------------------------
% Allocate gradient and Hessian
g = zeros(size(rho,1),size(rho,2),size(rho,3),2,'single');
H = zeros(size(rho,1),size(rho,2),size(rho,3),1,'single');

% -------------------------------------------------------------------------
% Compute conditional part (slice-wise to save memory)
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
    % Load one slice of the mean
    if size(rho, 5) == 2
        % Two real components
        rho1 = reshape(single(rho(:,:,z,:,:)), [], 2);
        rho1 = rho1(:,1) + 1i*rho1(:,2);
    else
        % One complex volume
        rho1 = reshape(single(rho(:,:,z,:)), [], 1);
    end
    if isa(A, 'gpuArray')
        rho1 = gpuArray(rho1);
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
    % Compute Hessian
    tmp = real(dot(s1 * A, s1, 2));
    H(:,:,z) = reshape(tmp, size(H,1), size(H,2));
    
    % ---------------------------------------------------------------------
    % Compute gradient
    tmp =  rho1 .* tmp - dot(s1 * A, x1, 2);
    g(:,:,z,1) = reshape(real(tmp), size(g,1), size(g,2));
    g(:,:,z,2) = reshape(imag(tmp), size(g,1), size(g,2));

end

% -------------------------------------------------------------------------
% Compute prior part
if size(rho, 5) == 2
    % Two real components
    rho0 = single(rho(:,:,:,:,:));
else
    % One complex volume
    rho0 = cat(5, single(real(rho(:,:,:,:,:))), single(imag(rho(:,:,:,:,:))));
end
g(:,:,:,1) = g(:,:,:,1) + spm_field('vel2mom', rho0(:,:,:,1,1), [vs prm]);
g(:,:,:,2) = g(:,:,:,2) + spm_field('vel2mom', rho0(:,:,:,1,2), [vs prm]);

% -------------------------------------------------------------------------
% Gauss-Newton
rho0(:,:,:,1,1) = rho0(:,:,:,1,1) - spm_field('fmg', H, g(:,:,:,1), [vs prm 2 2]);
rho0(:,:,:,1,2) = rho0(:,:,:,1,2) - spm_field('fmg', H, g(:,:,:,2), [vs prm 2 2]);


% -------------------------------------------------------------------------
% Write on disk
if size(rho, 5) == 2
    rho(:,:,:,:,:) = rho0;
else
    rho(:,:,:) = rho0(:,:,:,1,1) + 1i * rho0(:,:,:,1,2);
end