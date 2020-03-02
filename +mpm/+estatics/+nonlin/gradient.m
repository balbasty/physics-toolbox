function [llx,g,H] = gradient(in, out, opt)
% Non-linear ESTATICS: log-likelihood, gradient, Hessian of the data term.
%
% FORMAT [llx,g,H] = mpm.estatics.nonlin.gradient(in, out, [opt])
% in  - Input data structure
% out - Model data structure 
% opt - Structure of parameters with [optional] fields:
%       . subsample - Subsampling distance in mm [Inf=no]
%       . verbose   - Verbosity level [0]

% -------------------------------------------------------------------------
% Options
if nargin < 3, opt = struct; end
if ~isfield(opt, 'subsample'), opt.subsample = Inf;   end
if ~isfield(opt, 'verbose'),   opt.verbose   = 0;     end

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
    y       = exp(y0 - t * r2);                       % Echo fit
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
        g1 = lam * y .* r;
        g(:,:,:,1) = g(:,:,:,1) + g1; 
        g(:,:,:,2) = g(:,:,:,2) + g1 * (-t); 
        clear g1
    end
    clear r
    if nargout > 2
        H1 = lam * y.^2;
        H(:,:,:,1) = H(:,:,:,1) + H1;
        H(:,:,:,2) = H(:,:,:,2) + H1 * t^2;
        H(:,:,:,3) = H(:,:,:,3) + H1 * (-t);
        clear H1
    end
    clear y
end

% -------------------------------------------------------------------------
% Push gradient and Hessian to model space
if nargout > 1
    g = utils.push(factor*g, mat, ydim);
    if nargout > 2
        H = utils.push(factor*H, mat, ydim);
    end
end

llx = factor * llx;