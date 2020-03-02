function [ll,g,H] = gradient(dat,model,opt)
% Compute data-part of the gradient/Hessian w.r.t. the mean signal.
%
% FORMAT [ll,g,H] = b1p.bs.mean.gradient(dat,model,opt)
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

Nc = size(dat.coils,4); % Number of coils
Np = size(dat.coils,5); % Number of pulses

ll = 0;
g = single(0);
H = single(0);

meanim = single(model.mean());
% Loop over pulses
for p=1:Np
    w  = dat.pulse.sign(p);                         % sign of BS pulse
    k  = dat.pulse.factor;                          % sign of BS pulse
    b1 = single(model.b1());                        % square of b1
    if opt.b1.log, b1 = exp(2*b1); end
    b1 = exp(1i*w*k*b1);
    % Loop over coils
    for c=1:Nc
        lam = dat.prec(c,p);                        % precision
        s   = single(model.sens(:,:,:,c)) .* b1;    % sensitivity + shift
        y   = meanim .* s;                          % fit
        x   = single(dat.coils(:,:,:,c,p));         % observed
        msk = ~isfinite(x);                         % mask of missing data
        x(msk) = 0;
        y(msk) = 0;
        r = y - x;
        ll  = ll - 0.5 * lam * sum(abs(r(:)).^2, 'double'); 
        if nargout > 1
            g = g + lam * conj(r) .* s;
        end
        if nargout > 2
            H = H + lam * real(conj(s).*s);
        end
        if opt.verbose > 0, fprintf('.'); end
    end
end

if nargout > 1
    g = cat(4, real(g), imag(g));
end

% =========================================================================
%   SYMBOLIC DERIVATIONS
% =========================================================================
% 
% m     = sym('m','real');               % True object signal
% mc    = sym('mc','real');
% s     = sym('s','real');               % Receive sensitivity * b1 shift
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
% g = diff(L,m);
% H = diff(g,mc);
% 
% g_check  = 0.5 * lam * rc*s;
% simplify(g-g_check)
% 
% H_check = 0.5 * lam * sc*s;
% simplify(H-H_check)