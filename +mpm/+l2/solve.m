function dy = solve(H, g, mode, prec, vs)
% Solver for L2 spatial regularisation.
%
% FORMAT d = solve(H, g, mode, prec, vs)
% H    - {nx ny nz nf(nf+1)/2} - Field of (sparse) Hessian matrices
% g    - {nx ny nz nf}         - Field of gradients
% mode - {nf}                  - Regularisation mode for each channel (0|1|2)
% prec - {nf}                  - Regularisation precision for each channel
% vs   - {3}                   - Voxel size
% d    - {nx ny nz nf}         - Step: d = H\g

    fmg = [2 2];
    spm_field('boundary', 1);
    if numel(mode) == 1, mode = mode * ones(size(prec)); end
    prec = prec .* (mode > 0);
    dy = spm_field(single(H), single(g), [vs 0 1 0 fmg], prec(:)');

end