function scales = compute_scales(nbscales,dim,mat,sub)
% Create a structure of parameters for processing data in a multi-scale 
% manner.
%
% FORMAT scales = mpm.compute_scales(nbscales, dim, mat, sub)
% nbscales - Number of scale levels
% dim      - Volume dimension at the finer level
% mat      - Voxel-to-world matrix at the finer level
% sub      - Subsampling to apply at the finer level

dim = dim(:)';

vs = sqrt(sum(mat(1:3,1:3).^2));
scales = struct('dim', dim, 'mat', mat, 'vs', vs, 'sub', sub, 'scl', 1, 'ff', 1);

for i=2:nbscales
    scales(i).dim = ceil(scales(i-1).dim/2);
    scales(i).scl = dim./scales(i).dim;
    M = [diag(scales(i).scl) 0.5*(1-scales(i).scl(:));
         zeros(1,3) 1];
    scales(i).mat = mat * M;
    scales(i).vs  = dim./scales(i).dim.*vs;
    if isfinite(sub)
        scales(i).sub = scales(i-1).sub/2;
    else
        scales(i).sub = mean(scales(i).vs)/2;
    end
    scales(i).ff = 1; % scales(i-1).ff * 2;
end

end
