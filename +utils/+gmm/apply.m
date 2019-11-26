function [Z,X] = apply(X, mean, prec, prop, varargin)
% Use a learned mixture to segment an image.
%
% FORMAT [Z, X] = spm_gmm('apply',X,MU,A,PI,...)         > Classical
% FORMAT [Z, X] = spm_gmm('apply',X,{MU,b},{V,n},a,...)  > Bayesian
%
% KEYWORD
% -------
% Missing    - Infer missing data: ['infer']/'remove'
% BinWidth   - 1x[P] Bin width (histogram mode: add bits of variance) [0]
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Parse inputs
p = inputParser;
p.FunctionName = 'utils.gmm.';
p.addRequired('X',                     @isnumeric);
p.addRequired('mean',                  @(X) isnumeric(X) || iscell(X));
p.addRequired('precision',             @(X) isnumeric(X) || iscell(X));
p.addRequired('prop',                  @(X) isnumeric(X) || iscell(X));
p.addParameter('Missing',        true, @islogical);
p.addParameter('BinWidth', 0,    @isnumeric);
p.addParameter('Template',       [],   @isnumeric);
p.parse(X, mean, prec, prop, varargin{:});
E          = p.Results.BinWidth;
Missing    = p.Results.Missing;
Template   = p.Results.Template;

MU = [];
b  = [];
A  = [];
V  = [];
N  = [];
PI = [];
a  = [];
logPI = [];

% -------------------------------------------------------------------------
% Read arguments
if ~iscell(mean)
    MU = mean;
elseif ~isempty(mean)
    MU = mean{1};
    if numel(mean) > 1
        b = mean{2};
    end
end
if ~iscell(prec)
    A = prec;
elseif ~isempty(prec)
    A = prec{1};
    if numel(prec) > 1
        n = prec{2};
        if sum(b) > 0
            V = A;
            A = bsxfun(@times, V, reshape(n,1,1,[]));
            prec = {V n};
        end
    end
end
if ~iscell(prop)
    PI = prop;
elseif ~isempty(prop)
    PI = prop{1};
end

% -------------------------------------------------------------------------
% Dimensions
P = size(MU,1);
K = size(MU,2);

% -------------------------------------------------------------------------
% Proportions/Dirichlet
if ~isempty(Template)
    % Compute logPI by combining Template and [1xK] proportions in PI, as
    % in:
    % Ashburner J & Friston KJ. "Unified segmentation".
    % NeuroImage 26(3):839-851 (2005).
    logPI = bsxfun(@times,Template,PI);    
    logPI = log(bsxfun(@times,logPI,1./sum(logPI,2)));    
elseif numel(PI) <= K
    PI = PI(:)';
    PI = utils.pad(PI, [0 K - numel(PI)], 'replicate', 'post');
    
    if abs(sum(PI)-1) > eps('single')
        % Dirichlet prior
        a     = PI;
        logPI = psi(a) - psi(sum(a));
    else
        logPI = log(max(PI,eps));
    end
else
    logPI = log(max(PI,eps));
end
clear PI

% -------------------------------------------------------------------------
% Reshape X (observations)
dimX = size(X);
if P == 1
    latX = dimX;
    if latX(2)==1
        latX = latX(1);
    end
else
    latX = dimX(1:end-1);
end
X = reshape(X, [], P);
N  = size(X, 1); % Number of observations
N0 = N;          % Original number of observations

% -------------------------------------------------------------------------
% Prepare missing data stuff (code image, mask, ...)
if ~any(any(isnan(X)))
    Missing = false;
end
if Missing
    % Deal with missing data
    code      = utils.gmm.lib('obs2code', X);     % Code image
    code_list = unique(code);                   % List of codes
    missmsk   = [];
else
    % Compute mask of removed rows
    missmsk   = any(isnan(X),2);
    code      = utils.gmm.lib('double2int', (2^P-1) * ones(sum(~missmsk),1));
    code_list = 2^P-1;
end
% Discard rows with missing values
if ~isempty(missmsk)
    X         = X(~missmsk,:);
    if size(logPI,1) > 1
        logPI = logPI(~missmsk,:);
    end
    N         = sum(~missmsk);
end
missmsk = find(missmsk); % saves a bit of memory

% -------------------------------------------------------------------------
% "Bin" variance
% > When we work with histograms, a bit of variance is lost due to the
% binning. Here, we assume some kind of uniform distribution inside the bin
% and consequently add the corresponding variance to the 2nd order moment.
if numel(E) < P
    E = utils.pad(E, [0 P - numel(E)], 'replicate', 'post');
end
E = (E.^2)/12;

% -------------------------------------------------------------------------
% Compute marginal log-likelihood
if Missing
    const = utils.gmm.lib('Const', mean, prec, code_list);
else
    const  = utils.gmm.lib('Const', mean, prec);
end
logpX = utils.gmm.lib('Marginal', X, [{MU} prec], const, {code,code_list}, E);

% -------------------------------------------------------------------------
% Compute responsibilities
Z = utils.gmm.lib('Responsibility', logpX, logPI);
clear logpX logPI
 
% -------------------------------------------------------------------------
% Infer missing values
if nargout >= 2 && Missing
    X = utils.gmm.lib('InferMissing', X, Z, {MU,A}, {code,code_list});
end
 
% -------------------------------------------------------------------------
% Replace discarded missing values
if sum(missmsk) > 0
    present = ones(N0, 1, 'logical');
    present(missmsk) = false;
    clear missing
    
    Z = expand(Z, present, N0, K, 0);
    if nargout >= 2
        X = expand(X, present, N0, P, NaN);
    end
end
    
% -------------------------------------------------------------------------
% Reshape everything (use input lattice)
Z = reshape(Z, [latX K]);
if nargout >= 2
    X = reshape(X, [latX P]);
end