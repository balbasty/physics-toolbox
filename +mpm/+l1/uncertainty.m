function u = uncertainty(H, lam, vs, w)
% Diagonal uncertainty of reconstruction maps with MTV regularisation
%
% FORMAT u = mpm.l1.uncertainty(H, prec, vs, w)
% H    - {nx ny nz nf(nf+1)/2} - Field of (sparse) Hessian matrices
% prec - {nf}                  - Regularisation precision for each channel
% vs   - {3}                   - Voxel size
% w    - {nx ny nz}            - Weights (reweighted Gaussian approximation)
% u    - {nx ny nz nf}         - Posterior uncertainty (= variance)
%
% /!\ Warning: this function assumes that all channels have L1
%     regularisation (mixed L1/L2 not handled yet)

    vs  = vs(:)';
    lam = lam(:)';

    spm_field('boundary', 1);
    spm_diffeo('boundary', 1);
    u   = spm_field('diaginv1', H, w, [vs 1], lam);
    scl = spm_diffeo('kernel', [3 3 3], [vs 0 1 0 0 0]);
    scl = double(abs(scl(1,1,1)));
    u   = bsxfun(@times, u, reshape(scl*lam, 1, 1, 1, []));

end