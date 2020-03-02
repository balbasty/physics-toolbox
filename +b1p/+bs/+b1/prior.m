function [ll,g] = prior(b1,prec,vs)
% Compute prior-part of the gradient w.r.t. the B1 field.
%
% FORMAT [ll,g] = b1p.bs.b1.prior(b1,[prec],[vs])
% b1   | {r} {Nx Ny Nz} | B1 field
% prec | {r} {1}        | Prior precision [1]
% vs   | {r} {1 3}      | Voxel size [1 1 1]
% ll   | {r} {1}        | Log-likelihood of prior term
% g    | {r} {Nx Ny Nz} | Gradient of the prior term
%
% {r|c}    - Real or complex data
% Nx/Ny/Nz - Image dimensions

if nargin < 2, prec = [0 0 1]; end
if nargin < 3, vs   = [1 1 1]; end

b1 = single(b1());
g  = spm_field('vel2mom', b1, [vs prec]);
ll = -0.5 * b1(:)'*g(:);
