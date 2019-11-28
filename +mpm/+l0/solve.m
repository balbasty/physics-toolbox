function dy = solve(H,g)
% Solver for L0 spatial regularisation.
%
% FORMAT d = solve(H,g)
% H - {nx ny nz nf(nf+1)/2} - Field of (sparse) Hessian matrices
% g - {nx ny nz nf}         - Field of gradients
% d - {nx ny nz nf}         - Step: d = H\g

    spm_field('boundary', 1);
    dy = spm_field(H, g);

end