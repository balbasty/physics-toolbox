function y = vel2mom(y, prec, vs)
% Compute gradient of L2 regularisation operator.
%
% FORMAT Lf = mpm.l1.vel2mom(f, prec, vs)
% f    - {nx ny nz nf}         - Input image (single channel)
% prec - {1}                   - Regularisation precision for each channel
% vs   - {3}                   - Voxel size
% Lf   - {nx ny nz nf}         - Gradient of regulariser


    spm_field('boundary', 1);
    y = spm_field('vel2mom', single(y), [vs 0 prec 0]);

end