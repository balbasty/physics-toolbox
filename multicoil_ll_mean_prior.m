function llp = multicoil_ll_mean_prior(rho, prm, vs, part)
% Compute prior log-likelihood
%
% FORMAT ll = multicoil_ll_mean_prior(rho, prm, vs)
%
% rho  - (File)Array [Nx Ny Nz 1 (2)] - Complex mean image
% prm  -       Array [1 3]            - Regularisation [a m b]
% vs   -       Array [1 3]            - Voxel size [1 1 1]
% part -       Array [1 2]            - Factor for magnitude/phase [1 1] 
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 4
    part = [1 1];
end
if nargin < 3
    vs = [1 1 1];
end

% -------------------------------------------------------------------------
% Compute log-likelihood (prior)

rho1 = single(rho(:,:,:,:,:));
rho1 = cat(5, real(rho1), imag(rho1));

llpr = spm_field('vel2mom', rho1(:,:,:,:,1), [vs part(1)*prm]);
llpr = -0.5 * double(reshape(llpr, 1, [])) * double(reshape(rho1(:,:,:,:,1), [], 1));
llpi = spm_field('vel2mom', rho1(:,:,:,:,2), [vs part(2)*prm]);
llpi = -0.5 * double(reshape(llpi, 1, [])) * double(reshape(rho1(:,:,:,:,2), [], 1));
llp = llpr + llpi;