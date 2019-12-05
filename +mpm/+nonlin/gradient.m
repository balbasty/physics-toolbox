function [llx,g,H] = gradient(prm, in, out, pre, opt)
% FORMAT [llx,g,H] = mpm.nonlin.gradient(prm, in, out, [pre], [opt])
%
% INPUT
% -----
% prm - Parameters to optimise (e.g., {'A' 'R1' 'R2s' 'MT'})
% in  - Structure of input (observed) echos (FLASH [+ MT pulse])
% out - Structure of output (fitted) parameters (A, R1, R2, MT)
% pre - Structure of pre-computed parameters (B1t, B1r)
% opt - Structure of options
%
% OUTPUT
% ------
% llx - Data log-likelihood
% g   - Gradient of parameters
% H   - Hessian of parameters

% TODO: some of this should be done slice-wise. I don't know yet if loading
% the current log-maps should be done inside the slice loop, or outside to
% save some i/o.
% If done outside the loop, it means holding up to 6+6+10=22 volumes in
% memory at once (4 param maps + 2 b1 maps + gradient + Hessian).
%
% Currently (nothing slice-wise) I hold max 22 (12+10) volumes as well

if nargin < 4 || isempty(pre), pre = struct; end
if nargin < 5 || isempty(opt), opt = struct; end
if ~isfield(opt, 'subsample'), opt.subsample = Inf;   end
if ~isfield(opt, 'verbose'),   opt.verbose   = 0;     end

% -------------------------------------------------------------------------
% Read volume info
% -------------------------------------------------------------------------
ydim = out.dim;                 % Model dimensions
ymat = out.mat;                 % Model orientation matrix
xdim = in.dim;                  % Observed dimensions
xmat = in.mat;                  % Observed orientation matrix

ct  = in.type;                  % Contrast
mtw = strcmpi(in.SO, 'MT');     % Has MT pulse
TR  = in.TR;                    % Repetition time
FA  = in.FA;                    % Nominal flip angle
lam = 1./in.var;                % Noise precision

% -------------------------------------------------------------------------
% Undersampling
% -------------------------------------------------------------------------
vs     = sqrt(sum(in.mat(1:3,1:3).^2)); % Voxel size
skip   = opt.subsample * ones(1,3);
skip(~isfinite(skip)) = vs(~isfinite(skip));
skip   = round(skip./vs);
skip   = max(skip, 1);
Mskip  = [diag(skip) (ones(3,1)-skip(:)); zeros(1,3) 1];
xmat   = xmat * Mskip;
xdim0  = xdim;
xdim   = ceil(xdim0./skip);
factor = prod(xdim0./xdim);

% -------------------------------------------------------------------------
% Signal fit (TE=0)
% -------------------------------------------------------------------------

% --- Create sampling field
t = utils.warps.affine(ymat\xmat, xdim);

% --- Proton density (x receive sensitivity)
y0 = utils.pull(single(out.logA.dat()), t);
b1m = 1;
if isfield(pre, 'B1m')
    if     isfield(pre.B1m, ct),    B1m = pre.B1m.(ct);
    elseif isfield(pre.B1m, 'all'), B1m = pre.B1m.all;
    else,                           B1m = pre.B1m; end
    if strcmpi(B1m.unit, 'p.u.'), norm = 100; else, norm = 1; end
    b1m = pull(single(B1m.map()), xdim, B1m.mat\xmat); % /norm
    % b1m = log(max(b1m/norm, eps('single')));
end
y0 = exp(y0) .* b1m; clear b1m

% --- True flip angle (x transmit sensitivity)
b1p = 1;
if isfield(pre, 'B1p')
    if     isfield(pre.B1p, ct),    B1p = pre.B1p.(ct);
    elseif isfield(pre.B1p, 'all'), B1p = pre.B1p.all;
    else,                           B1p = pre.B1p; end
    if strcmpi(B1p.unit, 'p.u.'), norm = 100; else, norm = 1; end
    b1p = pull(single(B1p.map()), xdim, B1p.mat\xmat)/norm;
    % b1p = log(max(b1p/norm, eps('single')));
end
a    = FA .* b1p; clear b1p
cosa = cos(a);
% cosa = 1 - 0.5 * sina.^2;
y0   = y0 .* sin(a); clear a

% --- T1 exponential decay: exp(-R1*TR)
TRxR1      = utils.pull(single(out.logR1.dat()), t);
TRxR1      = exp(TRxR1) .* TR;
exp_mTRxR1 = exp(-TRxR1);

% --- MT ratio
if mtw
    exp_mLogMT = exp(-utils.pull(single(out.logMT.dat()), t));
    mt = 1./(1 + exp_mLogMT);                  % MTratio
    one_minus_mt = 1 - mt;              % 1 - MTratio
    d_one_minus_mt = -exp_mLogMT./(1 + exp_mLogMT).^2; % first derivative
    clear emt
end

% --- Signal fit (TE=0)
if mtw
    y0 = y0 .* one_minus_mt.* (1 - exp_mTRxR1) ./ (1 - cosa .* one_minus_mt .* exp_mTRxR1);
else
    y0 = y0 .* (1 - exp_mTRxR1) ./ (1 - cosa .* exp_mTRxR1);
end

% --- R2 decay
r2 = exp(utils.pull(single(out.logR2s.dat()), t));

% -------------------------------------------------------------------------
% Allocate gradient/Hessian
llx = 0;
F = numel(prm); % Number of features
if nargout > 1
    g = zeros([xdim F], 'single');
    [ind,K] = utils.symIndices(F); % Number of Hessian components
    if nargout > 2
        H = zeros([xdim K], 'single');
    end
    % Decide indices
    iA  = find(strcmpi('logA',   prm));
    iR1 = find(strcmpi('logR1',  prm));
    iR2 = find(strcmpi('logR2s', prm));
    iMT = find(strcmpi('logMT',  prm));
end
for e=1:numel(in.echoes)
    opt.verbose > 0 && fprintf('.');
    % ---------------------------------------------------------------------
    % Compute residuals
    dat     = in.echoes{e}.dat;
    y       = y0 .* exp(-r2 * in.echoes{e}.TE);       % Echo fit
    x       = single(dat(1:skip(1):end, ...
                         1:skip(2):end, ...
                         1:skip(3):end));             % Oberved echo
    msk     = isfinite(y) & isfinite(x) & (x > 0);    % Mask of observed voxels
    y(~msk) = 0;
    x(~msk) = 0;
    clear msk
    r       = y-x;                                    % Residuals
    clear x

    % ---------------------------------------------------------------------
    % Compute log-likelihood
    llx = llx - 0.5 * lam * sum(r(:).^2, 'double');

    % ---------------------------------------------------------------------
    % Compute gradient and Hessian in observed space
    if nargout > 1
        % -----------------------------------------------------------------
        % Proton density (A)
        if ~isempty(iA)
            g(:,:,:,iA) = y;
        end
        % -----------------------------------------------------------------
        % T2 decay (R2)
        if ~isempty(iR2)
            g(:,:,:,iR2) = - in.echoes{e}.TE .* r2 .* y;
        end
        % -----------------------------------------------------------------
        % T1 decay (R1)
        if ~isempty(iR1)
            if mtw
                g(:,:,:,iR1) = one_minus_mt .* cosa ./ (1 - one_minus_mt .* cosa .* exp_mTRxR1) - 1 ./ (1 - exp_mTRxR1);
            else
                g(:,:,:,iR1) = cosa ./ (1 - cosa .* exp_mTRxR1) - 1 ./ (1 - exp_mTRxR1);
            end
            g(:,:,:,iR1) = -TRxR1 .* exp_mTRxR1 .* y .* g(:,:,:,iR1);
        end
        % -----------------------------------------------------------------
        % MT ratio (MT)
        if ~isempty(iMT) && mtw
            g(:,:,:,iMT) = 1 ./ (one_minus_mt .* (1 - one_minus_mt .* cosa .* exp_mTRxR1));
            g(:,:,:,iMT) = d_one_minus_mt .* y .* g(:,:,:,iMT);
        end
    end
    if nargout > 2
        % -----------------------------------------------------------------
        % The expectation of the Hessian is obtained from the gradient
        % (squared), up to the residuals.
        for k=1:F
            for kk=k:F
                H(:,:,:,ind(k,kk)) = lam * g(:,:,:,k) .* g(:,:,:,kk);
            end
        end
    end
    if nargout > 1
        % -----------------------------------------------------------------
        % Multiply gradient with residuals
        g = lam * bsxfun(@times, g, r);
    end
    clear r
end

% -------------------------------------------------------------------------
% Push gradient and Hessian to model space
if nargout > 1
    g = utils.push(factor * g, t, ydim);
    if nargout > 2
        H = utils.push(factor * H, t, ydim);
    end
end
clear t

llx = factor * llx;

% -------------------------------------------------------------------------
% Check gradients for MPM nonlinear fit
% -------------------------------------------------------------------------
% 
% syms TE TR sina cosa 'positive'   % Fixed parameters
% syms lR1 lR2 lA ld 'real'         % Optimised parameters (log)
% syms X 'real'                     % Observation
% 
% % Precompute useful values
% A  = exp(lA);
% R1 = exp(lR1);
% R2 = exp(lR2);
% d  = 1 / (1 + exp(-ld));
% md = 1 - d;
% dd = - exp(-ld)./(1 + exp(-ld)).^2; % Derivative of (1-d) w.r.t. ld
% e1 = exp(-TR*R1);
% e2 = exp(-TE*R2);
% 
% % Signal fit and negative log-likelihood
% S = A * sina * md * (1 - e1) / (1 - md * cosa * e1) * e2;
% R = S - X;
% L = 0.5*R^2;
% 
% % Compute gradients automatically
% g = [diff(L, lA); ...
%      diff(L, lR2); ...
%      diff(L, lR1); ...
%      diff(L, ld)];
% H = [diff(L, lA,  lA) diff(L, lA,  lR2) diff(L, lA,  lR1) diff(L, lA,  ld); ...
%      diff(L, lR2, lA) diff(L, lR2, lR2) diff(L, lR2, lR1) diff(L, lR2, ld); ...
%      diff(L, lR1, lA) diff(L, lR1, lR2) diff(L, lR1, lR1) diff(L, lR1, ld); ...
%      diff(L, ld,  lA) diff(L, ld,  lR2) diff(L, ld,  lR1) diff(L, ld,  ld)];
% H = subs(H, X, S);
% 
% % Check that our gradients are correct
% gg = [S; ...
%       -TE * R2 * S; ...
%       -TR * R1 * e1 * S * (md * cosa / (1 - md * cosa .* e1) - 1 ./ (1 - e1)); ...
%       dd * S / (md * (1 - md * cosa * e1))];
% HH = gg * gg.';
% gg = gg * R;
% 
% simplify(g-gg, 100)
% simplify(H-HH, 100)