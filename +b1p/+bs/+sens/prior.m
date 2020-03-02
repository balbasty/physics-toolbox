function [ll,g] = prior(s,prec,vs)
% Compute prior-part of the gradient w.r.t. the receive field.
%
% FORMAT [ll,g] = b1p.bs.b1.prior(b1,[prec],[vs])
% s    | {c} {Nx Ny Nz} | Receive field
% prec | {r} {1}        | Prior precision [1]
% vs   | {r} {1 3}      | Voxel size [1 1 1]
% ll   | {r} {1}        | Log-likelihood of prior term
% g    | {r} {Nx Ny Nz} | Gradient of the prior term
%
% {r|c}    - Real or complex data
% Nx/Ny/Nz - Image dimensions

if nargin < 2, prec = 1;       end
if nargin < 3, vs   = [1 1 1]; end

s  = single(s());
g  = spm_field('vel2mom', s, [vs prec]);
ll = -0.5 * s(:)'*g(:);
