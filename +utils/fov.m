function [dim,mat] = fov(dims,mats,vs,fov,mat)
% Compute common field-of-view (dimensions + voxel-to-world mapping) from
% individual FOVs.
%
% FORMAT [dim,mat] = utils.fov(dims, mats, [vs], [fov])
% dims - [3xN]   Dimensions of each volume
% mats - [4x4xN] V2W matrix of each volume
% vs   - [3x1]   Voxel size of the recon maps [default: first/specified volume]
% fov  -         Index of FOV to use ([0] = bounding box)
% dim  - [3x1]   Output dimension
% mat  - [4x4]   Output V2W matrix

N  = size(dims,2);
if nargin < 3 || isempty(vs)
    vs = [NaN NaN NaN];
end
if nargin < 4
    fov = 0;
end
if nargin < 5
    mat = 1;
end
if size(mats, 3) ~= N 
    error('Number of volumes not consistant');
end

% Use orientation of the specified/first volume
if isscalar(mat)
    mat = mats(:,:,mat);
end

% Ensure target voxel size
vs0 = sqrt(sum(mat(1:3,1:3).^2));
vs(~isfinite(vs)) = vs0(~isfinite(vs));
mat = mat / diag([vs0 1]) * diag([vs 1]);

% Specified volume: adapt to voxel size
if ~isscalar(fov)
    dim = fov;
    return
end
if fov > 0
    dim = dims(:,fov)';
    dim = ceil(dim .* vs0 ./ vs);
    return
end

% No specified volume: compute embedding field of view
minpos =  [Inf Inf Inf];
maxpos = -[Inf Inf Inf];
for n=1:N
    mat1 = mat\mats(:,:,n);
    corners = [1         1         1; ...
               1         1         dims(3,n); ...
               1         dims(2,n) dims(3,n); ...
               1         dims(2,n) 1; ...
               dims(1,n) dims(2,n) dims(3,n); ...
               dims(1,n) dims(2,n) 1; ...
               dims(1,n) 1         dims(3,n); ...
               dims(1,n) 1         1];
    corners = bsxfun(@plus, corners * mat1(1:3,1:3)', mat1(1:3,4)');
    minpos  = min([minpos; corners], [], 1);
    maxpos  = max([maxpos; corners], [], 1);
end
shift = eye(4);
shift(1:3,4) = (minpos-1)';
mat = mat * shift;
minpos = [minpos 1];
maxpos = [maxpos 1];
dim    = ceil(max(maxpos,minpos) - min(maxpos,minpos)) + 1;
dim    = dim(1:3)';
    