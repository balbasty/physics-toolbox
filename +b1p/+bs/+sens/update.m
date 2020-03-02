function [model,llx,lls] = update(dat,model,opt)
% Update the receive fields.
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
% . mean        - {c} {Nx Ny Nz}       - Object signal (mean image)
% . sens        - {c} {Nx Ny Nz Nc}    - Receive sensitivities
% . b1          - {r} {Nx Ny Nz}       - Log B1 field or induced phase shift (B1^2)
% `opt` has fields:
% . vs          - {r} {1 3}            - Voxel size [1 1 1]
% . b1.log      - {l} {1}              - Log-encoding [true]
% . b1.prec     - {r} {1}              - Prior precision [1]
% . mean.prec   - {r} {1}              - Prior precision [1]
% . sens.prec   - {r} {1}              - Prior precision [1]
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
opt = utils.setdefault(opt, 'sens.prec',   1);
opt = utils.setdefault(opt, 'verbose',     1);

if opt.verbose > 0, fprintf('Receive field\n'); end

prec = [opt.sens.prec.abs opt.sens.prec.mem opt.sens.prec.ben];

for c=1:size(dat.coils,4)

    if opt.verbose > 0, fprintf('#%d\n', c); end
    
    sens = single(model.sens(:,:,:,c,:));
    sens = cat(4, real(sens), imag(sens));
    
    if opt.verbose > 0, fprintf('\tGradient: data '); end
    [llx,gx,Hx] = b1p.bs.sens.gradient(c, dat, model, opt);
    if opt.verbose > 0, fprintf('\n'); end

    if opt.verbose > 0, fprintf('\tGradient: prior '); end
    [lls,gs] = b1p.bs.sens.prior(sens, prec, opt.vs);
    gx       = gx + gs; clear gs
    if opt.verbose > 0, fprintf('.\n'); end

    if opt.verbose > 0, fprintf('\tUpdate: '); end
    sens(:,:,:,1) = sens(:,:,:,1) - spm_field(Hx, gx(:,:,:,1), [opt.vs prec]);
    sens(:,:,:,2) = sens(:,:,:,2) - spm_field(Hx, gx(:,:,:,2), [opt.vs prec]);
    model.sens(:,:,:,c) = complex(sens(:,:,:,1), sens(:,:,:,2));
    if opt.verbose > 0, fprintf('.\n'); end
    
end