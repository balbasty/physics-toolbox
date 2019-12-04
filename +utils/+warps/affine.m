function t = affine(M,dm)
% Create an identity warp (in voxel space)
%
% FORMAT t = utils.warps.affine(M,dim)
% M   - Voxel-to-voxel affine matrix
% dim - Dimensions (3D) of the lattice.
% t   - Transformation field (3D)

t = utils.warps.identity(dm);
t = utils.warps.compose(t,M);