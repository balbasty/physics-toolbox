function pre = prepare_precomputed(fB1t, fB1r)

pre = struct;

% -------------------------------------------------------------------------
% Transmit B1 field
% -------------------------------------------------------------------------
if ~isempty(fB1t)
    V = numel(fB1t);
    pre.B1t.dat  = cell(1,V);
    pre.B1t.vol  = zeros(1,V);        % Volume index
    pre.B1t.mat0 = zeros(4,4,V);      % Voxel-to-world matrix (native)
    pre.B1t.mat  = zeros(4,4,V);      % Voxel-to-world matrix (registered)
    pre.B1t.dim  = zeros(3,V);      % Dimensions
    pre.B1t.nii  = cell(1,V);         % Nifti struct of individual volume
    [pre.B1t.nii{:}] = deal(nifti);
    for v=1:V
        if V == 1
            pre.B1t.vol(v) = Inf;
        else
            pre.B1t.vol(v) = v;
        end
        pre.B1t.nii{v}      = nifti(fB1t{v});
        pre.B1t.dat{v}      = pre.B1t.nii{v}.dat;
        pre.B1t.mat0(:,:,v) = pre.B1t.nii{v}.mat;
        pre.B1t.mat(:,:,v)  = pre.B1t.nii{v}.mat;
        dim                 = [pre.B1t.nii{v}.dat.dim 1];
        pre.B1t.dim(:,v)    = dim(1:3)';
    end
end

% -------------------------------------------------------------------------
% Receive B1 field
% -------------------------------------------------------------------------
if ~isempty(fB1r)
    V = numel(fB1r);
    pre.B1r.dat  = cell(1,V);
    pre.B1r.vol  = zeros(1,V);        % Volume index
    pre.B1r.mat0 = zeros(4,4,V);      % Voxel-to-world matrix (native)
    pre.B1r.mat  = zeros(4,4,V);      % Voxel-to-world matrix (registered)
    pre.B1r.dim  = zeros(3,V);      % Dimensions
    pre.B1r.nii  = cell(1,V);         % Nifti struct of individual volume
    [pre.B1r.nii{:}] = deal(nifti);
    for v=1:V
        if V == 1
            pre.B1r.vol(v) = Inf;
        else
            pre.B1r.vol(v) = v;
        end
        pre.B1r.nii{v}      = nifti(fB1r{v});
        pre.B1r.dat{v}      = pre.B1r.nii{v}.dat;
        pre.B1r.mat0(:,:,v) = pre.B1r.nii{v}.mat;
        pre.B1r.mat(:,:,v)  = pre.B1r.nii{v}.mat;
        dim                 = [pre.B1r.nii{v}.dat.dim 1];
        pre.B1r.dim(:,v)    = dim(1:3)';
    end
end

end