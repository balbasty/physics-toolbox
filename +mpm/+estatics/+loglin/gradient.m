function [llx,g,H] = gradient(in, out, opt)
% Log-linear ESTATICS: log-likelihood, gradient, Hessian of the data term.
%
% FORMAT [llx,g,H] = mpm.estatics.loglin.gradient(in, out, [opt])
% in  - Input data structure
% out - Model data structure 
% opt - Structure of parameters with [optional] fields:
%   . scaled    - Scale noise precision by squared signal [false]
%   . subsample - Subsampling distance in mm [Inf=no]
%   . verbose   - Verbosity level [0]
%
% The 'scaled' option takes into account that 
%   Var[log x] \approx Var[x]/(E[x]^2)
% See: https://stats.stackexchange.com/questions/57715/
% It does not change anything when no regularisation is used.

% -------------------------------------------------------------------------
% Options
if nargin < 3, opt = struct; end
if ~isfield(opt, 'subsample'), opt.subsample = Inf;   end
if ~isfield(opt, 'verbose'),   opt.verbose   = 0;     end
if ~isfield(opt, 'scaled'),    opt.scaled    = false; end

% -------------------------------------------------------------------------
% Read volume info
ydim = [size(out.R2s.dat) 1];   % Model dimensions
ydim = ydim(1:3);
ymat = out.R2s.mat;             % Model orientation matrix
xmat = in.mat;                  % Observed orientation matrix

ct  = in.type;                  % Contrast name
mat = ymat\xmat;                % Observed-to-Model matrix
dim = in.dim;                   % Observed dimensions
lam = 1./in.var;                % Noise precision

% -------------------------------------------------------------------------
% Undersampling
vs    = sqrt(sum(in.mat(1:3,1:3).^2)); % Voxel size
skip  = opt.subsample * ones(1,3);
skip(~isfinite(skip)) = vs(~isfinite(skip));
skip  = round(skip./vs);
skip  = max(skip, 1);
Mskip = [diag(skip) (ones(3,1)-skip(:)); zeros(1,3) 1];
mat   = mat * Mskip;
dim0  = dim;
dim   = ceil(dim0./skip);
factor = prod(dim0./dim);

% -------------------------------------------------------------------------
% Load model data
y0  = utils.pull(single(out.(ct).dat()), mat, dim);
r2  = utils.pull(single(out.R2s.dat()),  mat, dim);

llx = 0;
if nargout > 1
    g = zeros([dim 2], 'single');
    if nargout > 2
        H = zeros([dim 3], 'single');
    end
end
for e=1:numel(in.echoes)
    opt.verbose > 0 && fprintf('.');
    % ---------------------------------------------------------------------
    % Compute residuals
    dat     = in.echoes{e}.dat;
    t       = in.echoes{e}.TE;                        % Echo time
    y       = y0 - t * r2;                            % Echo fit
    x       = log(single(dat(1:skip(1):end, ...
                             1:skip(2):end, ...
                             1:skip(3):end)));        % Oberved echo
    msk     = isfinite(y) & isfinite(x) & (x > 0);    % Mask of observed voxels
    y(~msk) = 0;
    x(~msk) = 0;
    clear msk
    r       = y-x;                                    % Residuals
    clear x
    if opt.scaled,  y = exp(2*y);
    else,           clear y; end

    % ---------------------------------------------------------------------
    % Compute log-likelihood
    if opt.scaled, llx = llx - 0.5 * lam * sum(y(:) .* r(:).^2, 'double');
    else,          llx = llx - 0.5 * lam * sum(r(:).^2, 'double'); end

    % ---------------------------------------------------------------------
    % Compute gradient and Hessian in observed space
    if nargout > 1
        g1 = lam * r;
        if opt.scaled, g1 = g1 .* y; end
        g(:,:,:,1) = g(:,:,:,1) + g1; 
        g(:,:,:,2) = g(:,:,:,2) + g1 * (-t); 
        clear g1
    end
    clear r
    if nargout > 2
        H1 = lam;
        if opt.scaled, H1 = H1 .* y; end
        H(:,:,:,1) = H(:,:,:,1) + H1;
        H(:,:,:,2) = H(:,:,:,2) + H1 * t^2;
        H(:,:,:,3) = H(:,:,:,3) + H1 * (-t);
        clear H1
    end
end

% -------------------------------------------------------------------------
% Push gradient and Hessian to model space
if nargout > 1
    g = utils.push(g, mat, ydim);
    if nargout > 2
        H = utils.push(H, mat, ydim);
    end
end

llx = factor * llx;
    