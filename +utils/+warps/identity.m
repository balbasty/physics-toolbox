function t = identity(dm)
% Create an identity warp (in voxel space)
%
% FORMAT t = utils.warps.identity(dim)
% dim - 3D dimensions of the lattice.
% t   - 3D transformation field (type: single)

dm = [dm 1 1 1];
dm = dm(1:3);

t = zeros(dm, 'single');
[t(:,:,:,1),t(:,:,:,2),t(:,:,:,3)] = ndgrid(1:dm(1), 1:dm(2), 1:dm(3));