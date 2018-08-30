function [rho,llm,llp,ok] = multicoil_sense_mean_map(x, s, A, rho, msk, dir, prm, vs, llm, llp)
% Maximum a posteriori mean given a set of observed coil images, 
% log-sensitivity profiles, a noise precision (= inverse covariance) 
% matrix and a K-space sampling scheme.
%
% /!\ Warning: the 3rd dimension should *NEVER* be accelerated.
%
% FORMAT [rho,llm,llp] = multicoil_sense_mean_map(x, s, A, (rho), ...)
%
% x   - (File)Array [Nx Ny Nz Nc (2)] - Complex coil images
% s   - (File)Array [Nx Ny Nz Nc (2)] - Complex log-sensitivity profiles
% A   -       Array [Nc Nc]           - Noise precision matrix
% rho - (File)Array [Nx Ny Nz  1 (2)] - Complex mean image
% msk - (File)Array [Nx Ny Nz]        - K-space sampling scheme
% dir -                               - List of accelerated directions
% prm -
% vs  -
%
% Nc = number of coils
% Images can either be complex or have two real components that are then 
% assumed to be the real and imaginary parts.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 7
    vs = [1 1 1];
end
if nargin < 9 || isnan(llm)
    llm = multicoil_sense_ll_cond(x,s,rho,A,msk,dir);
end
if nargin < 9 || isnan(llp)
    llp = multicoil_ll_mean_prior(rho,prm,vs);
end

fprintf('Update Mean (MAP):');

% -------------------------------------------------------------------------
% Fudge factor to account for the true number of independent observations
Nm = sum(msk(:));
No = numel(msk);
ff = Nm/No;

% -------------------------------------------------------------------------
% Allocate gradient and Hessian
g = zeros(size(rho,1),size(rho,2),size(rho,3),2,'single');
H = zeros(size(rho,1),size(rho,2),size(rho,3),1,'single');

% -------------------------------------------------------------------------
% Compute conditional part (slice-wise to save memory)
for z=1:size(rho, 3) 

    % ---------------------------------------------------------------------
    % Load one slice of the mask
    if size(msk,3) > 1
        msk1 = msk(:,:,z);
    else
        msk1 = msk;
    end
    
    % ---------------------------------------------------------------------
    % Load one slice of the *pushed* coil dataset
    if size(x, 5) == 2
        % Two real components
        x1 = single(x(:,:,z,:,:));
        x1 = x1(:,:,:,:,1) + 1i * x1(:,:,:,:,2);
    else
        % One complex volume
        x1 = single(x(:,:,z,:));
    end
    if isa(A, 'gpuArray')
        x1 = gpuArray(x1);
    end
    x1 = multicoil_pushwrap(x1, msk1, dir);
    x1 = reshape(x1, [], size(x,4));
    
    % ---------------------------------------------------------------------
    % Load one slice of the *native* bias dataset
    if size(s, 5) == 2
        % Two real components
        s1 = double(s(:,:,z,:,:));
        s1 = single(exp(s1(:,:,:,:,1) + 1i * s1(:,:,:,:,2)));
    else
        % One complex volume
        s1 = single(exp(double(s(:,:,z,:))));
    end
    if isa(A, 'gpuArray')
        s1 = gpuArray(s1);
    end
    
    % ---------------------------------------------------------------------
    % Load one slice of the *pulled* fit dataset (= mean x bias)
    if size(rho, 5) == 2
        % Two real components
        rho1 = single(rho(:,:,z,:,:));
        rho1 = rho1(:,:,:,:,1) + 1i * rho1(:,:,:,:,2);
    else
        % One complex volume
        rho1 = single(rho(:,:,z,:));
    end
    if isa(A, 'gpuArray')
        rho1 = gpuArray(rho1);
    end
    rho1 = bsxfun(@times, rho1, s1);
    rho1 = multicoil_pullwrap(rho1, msk1, dir);
    rho1 = reshape(rho1, [], size(s,4));
    s1   = reshape(s1, [], size(s,4));
    
    % ---------------------------------------------------------------------
    % Residuals / precompute
    rho1 = rho1 - x1; clear x1

    % ---------------------------------------------------------------------
    % Compute Hessian
    tmp = real(dot(s1 * A, s1, 2));
    H(:,:,z) = ff^2 * reshape(tmp, size(H,1), size(H,2));
    
    % ---------------------------------------------------------------------
    % Compute gradient
    tmp = dot(s1 * A, rho1, 2); clear s1 rho1
    g(:,:,z,1) = ff * reshape(real(tmp), size(g,1), size(g,2));
    g(:,:,z,2) = ff * reshape(imag(tmp), size(g,1), size(g,2));

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
drhor = spm_field('fmg', H, g(:,:,:,1), [vs prm 2 2]);
drhoi = spm_field('fmg', H, g(:,:,:,2), [vs prm 2 2]);

% -------------------------------------------------------------------------
% Constant parts of the log-likelihood
Ldrhor = spm_field('vel2mom', drhor, [vs prm]);
part1r  = reshape(rho0(:,:,:,:,1), 1, []) * reshape(Ldrhor, [], 1);
part2r  = reshape(drhor, 1, []) * reshape(Ldrhor, [], 1);
clear Ldrhor
Ldrhoi = spm_field('vel2mom', drhoi, [vs prm]);
part1i  = reshape(rho0(:,:,:,:,2), 1, []) * reshape(Ldrhoi, [], 1);
part2i  = reshape(drhoi, 1, []) * reshape(Ldrhoi, [], 1);
clear Ldrhoi

drho = cat(5,drhor, drhoi);
clear drhor drhoi

% -------------------------------------------------------------------------
% Line-Search
llp0   = llp;
llm0   = llm;
armijo = 1;
ok     = false;
for ls=1:6
    
    % ---------------------------------------------------------------------
    % Prior term
    llp = llp0 + armijo * (part1r + part1i) ...
               - 0.5 * armijo^2 * (part2r + part2i);
           
    % ---------------------------------------------------------------------
    % Conditional term
    llm = multicoil_sense_ll_cond(x, s, rho0 - armijo * drho, A, msk, dir);
    
    % ---------------------------------------------------------------------
    % Check progress
    if (llm + llp) > (llm0 + llp0)
        ok = true;
        break;
    else
        armijo = armijo/2;
    end
    
end

% -------------------------------------------------------------------------
% Write on disk
if ok
    fprintf(' :D (%d)', ls);
    rho0 = rho0 - armijo * drho;
    if size(rho, 5) == 2
        rho(:,:,:,:,:) = rho0;
    else
        rho(:,:,:) = rho0(:,:,:,1,1) + 1i * rho0(:,:,:,1,2);
    end
else
    fprintf(' :(');
    llm = llm0;
    llp = llp0;
end
fprintf('\n');