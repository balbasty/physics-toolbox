function [ll,g,H] = gradient(c,dat,model,opt)
% Compute data-part of the gradient/Hessian w.r.t. the sensitivity.
%
% FORMAT [ll,g,H] = b1p.bs.sens.gradient(c,dat,model,opt)
% c     - Index of coil to update
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
% . b1          - {r} {Nx Ny Nz}       - Phase shift induced by B1+ (K*B1^2)
% `opt` has fields:
% . b1.log      - {l} {1}              - Log-encoding of b1 [true]
% . verbose     - {r} {1}              - Verbosity level [1]
%
% {r|c}    - Real or complex data
% Nx/Ny/Nz - Image dimensions
% Nc       - Number of receive coils
% Np       - Number of pulses

if nargin < 3, opt = struct; end
opt = utils.setdefault(opt, 'b1.log',      true);
opt = utils.setdefault(opt, 'verbose',     1);

Np = size(dat.coils,5); % Number of pulses

ll = 0;
g = single(0);
H = single(0);

m = single(model.mean());
% Loop over pulses
for p=1:Np
    w  = dat.pulse.sign(p);                         % sign of BS pulse
    k  = dat.pulse.const;                           % const of BS pulse
    b1 = single(model.b1());                        % square of b1
    if opt.b1.log, b1 = exp(2*b1); end
    b1 = exp(1i*w*k*b1);                            % complex phase shift
    y0 = m .* b1;                                   % shifted object
    
    lam = dat.prec(c,p);                            % precision
    x   = single(dat.coils(:,:,:,c,o));             % observed
    msk = ~isfinite(x);                             % mask of missing data
    x(msk)  = 0;
    y0(msk) = 0;
    s   = single(model.sens(:,:,:,c));              % sensitivity
    r   = y0 .* s - x;                              % fit
    ll  = ll - 0.5 * lam * sum(abs(r).^2, 'double'); 
    if nargout > 1
        g = g + lam * conj(r) .* y0;
        g = cat(4, real(g), imag(g));
    end
    if nargout > 2
        H = H + lam * real(conj(y0) .* y0);
    end
    if opt.verbose > 0, fprintf('.'); end
end

% =========================================================================
%   SYMBOLIC DERIVATIONS
% =========================================================================
% 
% m     = sym('m','real');               % True object signal * b1 shift
% mc    = sym('mc','real');
% s     = sym('s','real');               % Receive sensitivity 
% sc    = sym('sc','real');
% x     = sym('x','real');               % Observed complex image
% xc    = sym('xc','real');
% lam   = sym('lam','positive');         % Noise precision
% 
% ms  = m .* s;
% msc = mc .* sc;
% r   = ms - x;
% rc  = msc - xc;
% 
% L = 0.5 * lam * (rc*r);                 % Sum-of-squares (negative log-likelihood)
% 
% g = diff(L,s);
% H = diff(g,sc);
% 
% g_check  = 0.5 * lam * rc*m;
% simplify(g-g_check)
% 
% H_check = 0.5 * lam * mc*m;
% simplify(H-H_check)