function llp = multicoil_ll_prior(s, prm, gamma, alpha, bnd, optim, vs)
% Compute prior log-likelihood
%
% FORMAT ll = multicoil_ll_prior(s, prm, vs)
%
% s   - (File)Array [Nx Ny Nz Nc (2)] - Complex log-sensitivity profiles
% RegCoilFactor - Vector [Nc]      - Reg modulator / coil          [1]
% RegCompFactor - [Mag Phase]      - Reg modulator / component     [1 1]
% vs  -       Array [1 3]             - Voxel size [1 1 1]
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

spm_field('boundary', bnd); 
optim = logical(optim);

% -------------------------------------------------------------------------
% Compute log-likelihood (prior)
llp = 0;
for n=1:size(s,4)
    if size(s, 5) == 2
        % Two real components
        s1 = single(s(:,:,:,n,optim));
    else
        % One complex volume
        if all(optim)
            s1 = single(s(:,:,:,n));
            s1 = cat(5, real(s1), imag(s1));
        elseif optim(1)
            s1 = real(single(s(:,:,:,n)));
        elseif optim(2)
            s1 = imag(single(s(:,:,:,n)));
        end
    end
    s1 = reshape(s1, size(s1,1),size(s1,2),size(s1,3),size(s1,5));
    llp1 = spm_field('vel2mom', s1, [vs alpha(n) * prm], gamma(optim));
    llp1 = -0.5 * double(reshape(llp1, 1, [])) * double(reshape(s1, [], 1));
    llp  = llp + llp1;
end