function f = adjoint(Kf,vs,bnd)
% FORMAT f = b1m.bending.adjoint(Kf, [vs], [bnd])
% Kf  - A 4D array Nx*Ny*Nz*3
% vs  - Voxel size [1 1 1]
% bnd - Boundary condition (0=circulant, [1]=neumann)
% f   - A 3D array Nx*Ny*Nz
%
% This is the adjoint form of the bending energy

dim = [size(Kf) 1 1];
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
Kf = utils.pad(Kf, [1 1 1 0], bnd, 'both');

S0 = struct('type', '()', 'subs', {{':',':',':',':'}});
f = zeros(dim, 'like', Kf);
for d=1:3
    S = S0;
    S.subs{4} = d;
    for dd=1:3
        S.subs{dd} = 2:(size(Kf,dd)-1);
    end
    f = f - 2 * subsref(Kf, S) / vs(d);
    
    S = S0;
    S.subs{4} = d;
    for dd=1:3
        if dd~=d
            S.subs{dd} = 2:(size(Kf,dd)-1);
        else
            S.subs{dd} = 3:size(Kf,dd);
        end
    end
    f = f + subsref(Kf, S) / vs(d);
    
    S = S0;
    S.subs{4} = d;
    for dd=1:3
        if dd~=d
            S.subs{dd} = 2:(size(Kf,dd)-1);
        else
            S.subs{dd} = 1:(size(Kf,dd)-2);
        end
    end
    f = f + subsref(Kf, S) / vs(d);
end
