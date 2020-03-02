function [ll,g] = prior(m,prec,vs)
% Compute prior-part of the gradient w.r.t. the mean.
%
% FORMAT [ll,g] = b1p.bs.mean.prior(m,[prec],[vs])
% m    | {r} {Nx Ny Nz} | Mean image
% prec | {r} {3}        | Prior precision [0 1 0]
% vs   | {r} {1 3}      | Voxel size [1 1 1]
% ll   | {r} {1}        | Log-likelihood of prior term
% g    | {r} {Nx Ny Nz} | Gradient of the prior term
%
% {r|c}    - Real or complex data
% Nx/Ny/Nz - Image dimensions

if nargin < 2, prec = [0 1 0]; end
if nargin < 3, vs   = [1 1 1]; end

if any(prec > 0)
    m  = single(m());
    g  = spm_field('vel2mom', m, [vs prec]);
    ll = -0.5 * m(:)'*g(:);
else
    g  = single(0);
    ll = 0;
end
