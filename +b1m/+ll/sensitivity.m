function ll = sensitivity(sens, prm, partfactor, coilfactor, bnd, optim, vs)
% Compute prior log-likelihood of the sensitivities
%
% FORMAT ll = b1m.ll.sensitivity(SensMaps, RegStructure, RegPartFactor, 
%                                RegCoilFactor, RegBoundary, SensOptim, 
%                                VoxelSize)
%
% SensMaps      - Complex log-sensitivity profiles          (File)Array [Nx Ny Nz Nc]
% RegStructure  - Regularisation Structure (abs memb bend)  [0 0 1]
% RegPartFactor - Regularisation factor per Re/Im component [1E6]
% RegCoilFactor - Regularisation factor per coil            [1/Nc]
% RegBoundary   - Boundary conditions for sensitivities     [1=neumann]
% SensOptim     - Optimize real and/or imaginary parts      [true true]
% VoxelSize     - Voxel size                                [1]
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 7
    vs = [1 1 1];
    if nargin < 6
        optim = [true true];
        if nargin < 5
            bnd = 1;
            if nargin < 4
                coilfactor = 1/size(sens,4);
                if nargin < 3
                    partfactor = [1E6 1E6];
                    if nargin < 2
                        prm = [0 0 1];
                    end
                end
            end
        end
    end
end

spm_field('boundary', bnd); 
optim = logical(optim);

% -------------------------------------------------------------------------
% Compute log-likelihood (prior)
ll = zeros(1,size(sens,4));
for n=1:size(sens,4)
    if size(sens, 5) == 2
        % Two real components
        s1 = single(sens(:,:,:,n,optim));
    else
        % One complex volume
        if all(optim)
            s1 = single(sens(:,:,:,n));
            s1 = cat(5, real(s1), imag(s1));
        elseif optim(1)
            s1 = real(single(sens(:,:,:,n)));
        elseif optim(2)
            s1 = imag(single(sens(:,:,:,n)));
        end
    end
    s1 = reshape(s1, size(s1,1),size(s1,2),size(s1,3),size(s1,5));
    llp1 = spm_field('vel2mom', gather(s1), [vs coilfactor(n) * prm], partfactor(optim));
    llp1 = -0.5 * double(reshape(llp1, 1, [])) * double(reshape(s1, [], 1));
    ll(n)= llp1;
end