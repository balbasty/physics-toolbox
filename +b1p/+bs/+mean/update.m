function [model,llx,llm] = update(dat,model,opt)
% Update the mean image.
%
% FORMAT [model,llx,llb] = b1p.bs.mean.update(dat,model,opt)
% dat   - Input data structure (holding acquired images and parameters)
% model - Model structure (holding inferred quantities)
% opt   - Option structure
% ll    - {r} {1} Log-likelihood of data term
% g     - {r} {Nx Ny Nz} Gradient 
% H     - {r} {Nx Ny Nz} Hessian 
%
% `dat` has fields:
% . coils        - {c} {Nx Ny Nz Nc Np} - Observed coil images with +/- pulse
% . prec         - {r} {Nc Np}          - Noise precision
% . pulse.sign   - {r} {1 Np}           - Sign of Bloch-Siegert pulse
% . pulse.factor - {r} {1}              - Constant of Bloch-Siegert pulse
% `model` has fields
% . mean        - {c} {Nx Ny Nz}       - Object signal (mean image)
% . sens        - {c} {Nx Ny Nz Nc}    - Receive sensitivities
% . b1          - {r} {Nx Ny Nz}       - Phase shift induced by B1+ (K*B1^2)
% `opt` has fields:
% . vs          - {r} {1 3}            - Voxel size [1 1 1]
% . b1.prec     - {r} {1}              - Prior precision [1]
% . b1.log      - {l} {1}              - Log-encoding [true]
% . meanim.prec - {r} {1}              - Prior precision [1]
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

if opt.verbose > 0, fprintf('Mean image\n'); end

meanim = single(model.mean());
meanim = cat(4, real(meanim), imag(meanim));
prec   = [opt.mean.prec.abs opt.mean.prec.mem opt.mean.prec.ben];

if opt.verbose > 0, fprintf('\tGradient: data '); end
[llx,gx,Hx] = b1p.bs.mean.gradient(dat, model, opt);
if opt.verbose > 0, fprintf('\n'); end

if opt.verbose > 0, fprintf('\tGradient: prior '); end
[llm,gm] = b1p.bs.mean.prior(meanim, prec, dat.vs);
gx       = gx + gm; clear gm
if opt.verbose > 0, fprintf('.\n'); end


if opt.verbose > 0, fprintf('\tUpdate: '); end
meanim(:,:,:,1) = meanim(:,:,:,1) - spm_field(Hx, gx(:,:,:,1), [dat.vs prec]);
meanim(:,:,:,2) = meanim(:,:,:,2) - spm_field(Hx, gx(:,:,:,2), [dat.vs prec]);
model.mean(:,:,:) = complex(meanim(:,:,:,1), meanim(:,:,:,2));
if opt.verbose > 0, fprintf('.\n'); end