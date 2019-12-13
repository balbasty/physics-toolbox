function [model,llx,llb] = update(dat,model,opt)
% Update the B1 field.
%
% FORMAT [model,llx,llb] = b1p.bs.b1.update(dat,model,opt)
% dat   - Input data structure (holding acquired images and parameters)
% model - Model structure (holding inferred quantities)
% opt   - Option structure
% ll    - {r} {1} Log-likelihood of data term
% g     - {r} {Nx Ny Nz} Gradient 
% H     - {r} {Nx Ny Nz} Hessian 
%
% `dat` has fields:
% . coils       - {c} {Nx Ny Nz Nc Np} - Observed coil images with +/- pulse
% . prec        - {r} {Nc Np}          - Noise precision
% . pulse.sign  - {r} {1 Np}           - Sign of Bloch-Siegert pulse
% . pulse.const - {r} {1}              - Constant of Bloch-Siegert pulse
% `model` has fields
% . meanim      - {c} {Nx Ny Nz}       - Object signal (mean image)
% . sens        - {c} {Nx Ny Nz Nc}    - Receive sensitivities
% . b1          - {r} {Nx Ny Nz}       - Phase shift induced by B1+ (K*B1^2)
% `opt` has fields:
% . vs          - {r} {1 3}            - Voxel size [1 1 1]
% . prec        - {r} {1}              - Prior precision [1]
% . log         - {l} {1}              - Log-encoding [true]
% . verbose     - {r} {1}              - Verbosity level [1]
%
% {r|c}    - Real or complex data
% Nx/Ny/Nz - Image dimensions
% Nc       - Number of receive coils
% Np       - Number of pulses

if nargin < 3, opt = struct; end
opt = utils.setdefault(opt, 'vs',          [1 1 1]);
opt = utils.setdefault(opt, 'b1.log',      true);
opt = utils.setdefault(opt, 'b1.prec',     1);
opt = utils.setdefault(opt, 'mean.prec',   1);
opt = utils.setdefault(opt, 'verbose',     1);

if opt.verbose > 0, fprintf('B1 field\n'); end

b1 = single(model.b1());

if opt.verbose > 0, fprintf('\tGradient: data '); end
[llx,gx,Hx] = b1p.bs.b1.gradient(dat, model, opt);
if opt.verbose > 0, fprintf('\n'); end

if opt.verbose > 0, fprintf('\tGradient: prior '); end
[llb,gb]    = b1p.bs.b1.prior(b1, opt.b1.prec, opt.vs);
if opt.verbose > 0, fprintf('.\n'); end

if opt.verbose > 0, fprintf('\tUpdate: '); end
model.b1(:,:,:) = b1 - spm_field(Hx, gx+gb, [opt.vs 0 0 opt.b1.prec]);
if opt.verbose > 0, fprintf('.\n'); end