function llp = multicoil_ll_prior(s, prm, vs)
% Compute prior log-likelihood
%
% FORMAT ll = multicoil_ll_prior(s, prm, vs)
%
% s   - (File)Array [Nx Ny Nz Nc (2)] - Complex log-sensitivity profiles
% prm -       Array [1 3] or [Nc 3]   - Regularisation (/ coil) [a m b]
% vs  -       Array [1 3]             - Voxel size [1 1 1]
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 3
    vs = [1 1 1];
end
if size(prm, 1) == 1
    prm = repmat(prm, [size(s,4) 1]);
end

% -------------------------------------------------------------------------
% Compute log-likelihood (prior)
llp = 0;
for n=1:N
    if size(s, 5) == 2
        % Two real components
        s1 = single(s(:,:,:,n,:));
    else
        % One complex volume
        s1 = single(s(:,:,:,n,:));
        s1 = cat(5, real(s1), imag(s1));
    end
    llpr = spm_field('vel2mom', s1(:,:,:,:,1), [vs prm(n,:)]);
    llpr = -0.5 * double(reshape(llpr, 1, [])) * double(reshape(s1(:,:,:,:,1), [], 1));
    llpi = spm_field('vel2mom', s1(:,:,:,:,2), [vs prm(n,:)]);
    llpi = -0.5 * double(reshape(llpi, 1, [])) * double(reshape(s1(:,:,:,:,2), [], 1));
    llp = llp + llpr + llpi;
end