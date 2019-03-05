function llp = mean(mean, prm, vs, part)
% Compute prior log-likelihood of the mean
%
% FORMAT ll = b1m.ll.mean(mean, prm, [vs], [part])
%
% mean - (File)Array [Nx Ny Nz] - Complex mean image
% prm  -       Array [1 3]      - Regularisation [a m b]
% vs   -       Array [1 3]      - Voxel size [1 1 1]
% part -       Array [1 2]      - Factor for magnitude/phase [1 1] 
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

r1 = single(mean(:,:,:,:,:));
r1 = cat(5, real(r1), imag(r1));

llpr = spm_field('vel2mom', r1(:,:,:,:,1), [vs part(1)*prm]);
llpr = -0.5 * double(reshape(llpr, 1, [])) * double(reshape(r1(:,:,:,:,1), [], 1));
llpi = spm_field('vel2mom', r1(:,:,:,:,2), [vs part(2)*prm]);
llpi = -0.5 * double(reshape(llpi, 1, [])) * double(reshape(r1(:,:,:,:,2), [], 1));
llp = llpr + llpi;