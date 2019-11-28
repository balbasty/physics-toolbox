function varargout = residuals_estatics(n, in, out, subsample, which)
% Compute residuals of the ESTATICS model.
%
% FORMAT [...] = residuals_estatics(n, in, model, [subsample], [which])
% n         - Index of input volume to process
% in        - Structure of input file obtained using `prepare_input`
% model     - Structure of model files obtained using `estatics_nonlinfit`
%             or `estatics_loglinfit`.
% subsample - Subsample the input data? [false]
% which     - Specify variables to output (char or cell fo char).
%             Output variables can be:
%             ['r'] - Image of residuals
%             'r2'  - Image of squared residual
%             'll'  - Image of log-likelihood (assumed Gaussian, using 
%                     estimated noise variance)
%             's*'  - Prefixing the above options with an s returns their 
%                     sum across voxels.

    if nargin < 4 || isnan(subsample)
        subsample = Inf;
    end
    if nargin < 5
        which = {'r'};
    end
    if ~iscell(which)
        which = {which};
    end

    % ---------------------------------------------------------------------
    % Read individual info
    dat = in.dat{n};                % File array
    i   = in.idx(n);                % Parameter index (PD/T1/MT)
    v   = in.vol(n);                % Volume index
    t   = in.TE(n);                 % Echo time
    mat = out.mat\in.mat(:,:,v);    % Observed-to-Common matrix
    dim = in.dim(:,v)';             % Observed dimensions
    lam = 1./in.var(v);             % Noise precision
    
    % ---------------------------------------------------------------------
    % Undersampling
    vs    = sqrt(sum(in.mat(1:3,1:3,v).^2)); % Voxel size
    skip  = subsample * ones(1,3);
    skip(~isfinite(skip)) = vs(~isfinite(skip));
    skip  = round(skip./vs);
    skip  = max(skip, 1);
    Mskip = [diag(skip) (ones(3,1)-skip(:)); zeros(1,3) 1];
    mat   = mat * Mskip;
    dim0  = dim;
    dim   = ceil(dim0./skip);
    % factor = prod(dim0./dim);
    
    % ---------------------------------------------------------------------
    % Load data
    y   = pull(single(out.dat(:,:,:,i)), dim, mat);         % Intercept
    y   = y - t * pull(single(out.dat(:,:,:,4)), dim, mat); % Decay
    y   = exp(y);                                           % Echo fit
    x   = single(dat(1:skip(1):end, ...
                     1:skip(2):end, ...
                     1:skip(3):end));                       % Oberved echo
    msk = isfinite(x) & (x > 0);                            % Mask of observed voxels
    
    
    % ---------------------------------------------------------------------
    % Residuals
    r   = y-x;
    clear x y
    if cellfun(@(X) any(strcmpi(X, {'sr'})), which)
        sr = sum(r(msk), 'double');
    end
    
    % ---------------------------------------------------------------------
    % Squared residuals
    if cellfun(@(X) any(strcmpi(X, {'r2' 'sr2' 'll' 'sll'})), which)
        r2 = r.^2;
    end
    if ~cellfun(@(X) any(strcmpi(X, {'r'})), which)
        r = [];
    end
    if cellfun(@(X) any(strcmpi(X, {'sr2'})), which)
        sr2 = sum(r2(msk), 'double');
    end
    
    % ---------------------------------------------------------------------
    % Log-likelihood
    if cellfun(@(X) any(strcmpi(X, {'ll' 'sll'})), which)
        ll = - 0.5 * (lam * r2 - log(lam) + log(2*pi));
    end
    if ~cellfun(@(X) any(strcmpi(X, {'r2'})), which)
        r2 = [];
    end
    if cellfun(@(X) any(strcmpi(X, {'sll'})), which)
        sll = sum(ll(msk), 'double');
    end
    if ~cellfun(@(X) any(strcmpi(X, {'ll'})), which)
        ll = [];
    end
    
    clear msk
    
    % ---------------------------------------------------------------------
    % Fill varargout
    varargout = cell(1,numel(which));
    for i=1:numel(which)
        switch lower(which{i})
            case 'r'
                varargout{i} = r;
            case 'r2'
                varargout{i} = r2;
            case 'll'
                varargout{i} = ll;
            case 'sr'
                varargout{i} = sr;
            case 'sr2'
                varargout{i} = sr2;
            case 'sll'
                varargout{i} = sll;
        end
    end
    
end