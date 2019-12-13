function [ll,g,H] = gradient(dat,model,opt)
% Compute data-part of the gradient/Hessian w.r.t. the B1 field.
%
% FORMAT [ll,g,H] = b1p.bs.b1.gradient(dat,model,opt)
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
% . log         - {l} {1}              - Log-encoding [true]
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
    k  = dat.pulse.const;                           % sign of BS pulse
    b1 = single(model.b1());                        % square of b1
    if opt.b1.log, b1 = exp(2*b1); end
    y0 = meanim .* exp(1i*w*k*b1);
    % Loop over coils
    for c=1:Nc
        lam = dat.prec(c,p);                        % precision
        y   = single(model.sens(:,:,:,c)) .* y0;    % fit
        x   = single(dat.coils(:,:,:,c,o));         % observed
        msk = ~isfinite(x);                         % mask of missing data
        x(msk) = 0;
        y(msk) = 0;
        ll  = ll - 0.5 * lam * sum(abs(y-x).^2, 'double'); 
        if nargout > 1
            g = g - lam * w * k * imag(conj(y).*x);
            if opt.b1.log, g = 2 * g * b1; end
        end
        if nargout > 2
            H = H + lam * abs(y).^2;
            if opt.b1.log, H = 2 * H * b1.^2; end
        end
        if opt.verbose > 0, fprintf('.'); end
    end
end

% =========================================================================
%   SYMBOLIC DERIVATIONS
% =========================================================================
% 
% -------------------------------------------------------------------------
%   Parameterised by b1^2
% -------------------------------------------------------------------------
%
% m     = sym('m');                       % True object signal * receive modulation
% b1    = sym('b1','positive');           % Squared B1+ field
% x     = sym('x');                       % Observed complex image
% lam   = sym('lam','positive');          % Noise precision
% pulse = sym('pulse','real');            % Sign of Bloch-Siegert pulse
% 
% mp = m * exp(1i*b1*pulse);              % Modulated object signal + phase shift
% r  = mp - x;                            % Complex residuals
% 
% L = 0.5 * lam * (r'*r);                 % Sum-of-squares (negative log-likelihood)
% 
% g = diff(L,b1);
% H = diff(g,b1);
% 
% sym_imag = @(x) (x-conj(x))/(2*1i);     % Symbolic-friendly definition of imag
% sym_real = @(x) (x+conj(x))/2;          % Symbolic-friendly definition of real
% 
% g_check  =  lam * pulse * sym_imag(mp'*r);
% g_check2 = -lam * pulse * sym_imag(mp'*x);
% simplify(g-g_check)
% simplify(g-g_check2)
% 
% H_check = lam * pulse^2 * sym_real(mp'*x);
% simplify(H-H_check)
% % ^ note that we always have: pulse^2 == 1
%
% -------------------------------------------------------------------------
%   Parameterised by log(b1)
% -------------------------------------------------------------------------
%
% m     = sym('m');                       % True object signal * receive modulation
% p     = sym('b1','real');               % Squared B1+ field
% x     = sym('x');                       % Observed complex image
% lam   = sym('lam','positive');          % Noise precision
% pulse = sym('pulse','real');            % Sign of Bloch-Siegert pulse
% 
% b1 = exp(2*p);
% mp = m * exp(1i*b1*pulse);              % Modulated object signal + phase shift
% r  = mp - x;                            % Complex residuals
% 
% L = 0.5 * lam * (r'*r);                 % Sum-of-squares (negative log-likelihood)
% 
% g = diff(L,p);
% H = diff(g,p);
% 
% sym_imag = @(x) (x-conj(x))/(2*1i);     % Symbolic-friendly definition of imag
% sym_real = @(x) (x+conj(x))/2;          % Symbolic-friendly definition of real
% 
% g_check  =  lam * pulse * 2 * b1 * sym_imag(mp'*r);
% g_check2 = -lam * pulse * 2 * b1 * sym_imag(mp'*x);
% simplify(g-g_check)
% simplify(g-g_check2)
% 
% H_check = lam * pulse^2 * 4 * b1^2 * sym_real(mp'*x) ...
%          -lam * pulse   * 4 * b1   * sym_imag(mp'*x);
% simplify(H-H_check)
% H_opt_check = lam * pulse^2 * 4 * b1^2 * sym_real(mp'*x);
% % At optimum, sym_imag(mp'*x) == 0
% % Note that we always have: pulse^2 == 1
