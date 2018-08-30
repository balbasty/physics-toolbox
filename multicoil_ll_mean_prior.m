function llp = multicoil_ll_mean_prior(rho, prm, vs)
% Compute prior log-likelihood
%
% FORMAT ll = multicoil_ll_mean_prior(rho, prm, vs)
%
% rho   - (File)Array [Nx Ny Nz 1 (2)] - Complex mean image
% prm -         Array [1 3]            - Regularisation [a m b]
% vs  -         Array [1 3]            - Voxel size [1 1 1]
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 3
    vs = [1 1 1];
end

% -------------------------------------------------------------------------
% Compute log-likelihood (prior)

if size(rho, 5) == 2
    % Two real components
    rho1 = single(rho(:,:,:,:,:));
else
    % One complex volume
    rho1 = single(rho(:,:,:,:,:));
    rho1 = cat(5, real(rho1), imag(rho1));
end
llpr = spm_field('vel2mom', rho1(:,:,:,:,1), [vs prm]);
llpr = -0.5 * double(reshape(llpr, 1, [])) * double(reshape(rho1(:,:,:,:,1), [], 1));
llpi = spm_field('vel2mom', rho1(:,:,:,:,2), [vs prm]);
llpi = -0.5 * double(reshape(llpi, 1, [])) * double(reshape(rho1(:,:,:,:,2), [], 1));
llp = llpr + llpi;