function [g,mat] = upsample(f, dim)
% Trilinear upsampling of volume.
%
% FORMAT [f,mat] = utils.upsample(f, dim)
% f   - Input volume (First 3 dimensions must be spatial)
% dim - Output dimensions
% mat - Voxel-to-voxel transformation matrix
%
% Upsample only applied to the first three dimensions, which are considered
% to be spatial dimensions.

dim0 = [size(f) 1];
dimf = dim0(4:end);     % Feature dimensions
dim0 = dim0(1:3);       % Original spatial dimensions
dim  = [dim(:)' 1 1 1];
dim  = dim(1:3);        % Output spatial dimensions

% --- Create transformation matrix
scale = dim0 ./ dim;
mat = [diag(scale) 0.5*(1-scale(:));
       zeros(1,3) 1];

% --- Build transformation field
t = zeros([dim 3], 'single');
[t(:,:,:,1),t(:,:,:,2),t(:,:,:,3)] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
t = reshape(t, [], 3);
t = bsxfun(@plus, t * mat(1:3,1:3)', mat(1:3,4)');
t = reshape(t, [dim 3]);

% --- Upsample
f = reshape(f, [dim0 prod(dimf)]);
g = zeros([dim prod(dimf)], 'like', f);
spm_diffeo('boundary', 1);
for k=1:size(f,4)
    g(:,:,:,k) = spm_diffeo('bsplins', single(f(:,:,:,k)), t, [1 1 1  1 1 1]);
end
g = reshape(g, [dim dimf]);

