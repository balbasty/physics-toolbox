function Kf = forward(f,vs,bnd)
% FORMAT Kf = b1m.bending.forward(f, [vs], [bnd])
% f   - A 3D array Nx*Ny*Nz
% vs  - Voxel size [1 1 1]
% bnd - Boundary condition (0=circulant, [1]=neumann)
% Kf  - A 4D array Nx*Ny*Nz*3
%
% Compute second spatial derivatives along x/y/z in each voxel.
% This is the forward form of the bending energy

dim = [size(f) 1];
dim = dim(1:3);

vs = vs.^2; % pre-square

if nargin < 3, bnd = 1; end
switch bnd
    case 0
        bnd = 'circulant';
    case 1
        bnd = 'symmetric';
    otherwise
        error('Boundary condition must be 0 or 1')
end
f = utils.pad(f, [1 1 1], bnd, 'both');
% f = f(1:(end-1),1:(end-1),1:(end-1));

S0 = struct('type', '()', 'subs', {{':',':',':'}});
Kf = zeros([dim 3], 'like', f);
for d=1:3
    S = S0;
    for dd=1:3
        if dd~=d
            S.subs{dd} = 2:(size(f,dd)-1);
        end
    end
    Kf(:,:,:,d) = diff(subsref(f,S),2,d)/vs(d);
end
