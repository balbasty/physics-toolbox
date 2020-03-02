function [sens,meanim,prec] = fit_robust(varargin)
% Compute mode estimates (ML, MAP) of the parameters of a probabilistic 
% model of complex multicoil MR images.
%
%   This function performs a robust fit on complex datasets:
%       1. A fit using log-sensitivities is performed to get a good initial
%          estimate of the mean image. This fit is more robust because
%          log-sensitivities can be zero-centered, yielding a unique
%          optimum. However, smooth sensitivity phases are assumed, which
%          does not hold. Therefore, the output sensitivities cannot be
%          used as is.
%       2. The mean image is kept fixed and complex sensitivities are
%          fitted. Since the mean image is fixed, it is a convex problem
%          with a unique solution. A coarse to fine scheme is used to
%          accelerte convergence.
%       3. Finally, mean and sensitivity estimates are refined by
%          updating them in turn. The problem is not jointly convex, but
%          the starting estimates should be good enough to yield a good
%          optimum.
%
% FORMAT [sens,mean,prec] = b1m.fit_robust(coils, ...)
%
% REQUIRED
% --------
% coils - (File)Array [Nx Ny Nz Nch Nct] - Complex coil images
%
% KEYWORDS
% --------
% VoxelSize         - Voxel size                                [1]
% SamplingMask      - Mask of the k-space sampling scheme       [[]=fully sampled]
% Precision         - Noise precision matrix                    [NaN=estimate from fully-sampled]
% SensRegFactor     - Regularisation factor                     [1E5 1E4 1E4]
% SensRegDecFactor  - Start with more regularisation            [3 3 0] (log10)
% MeanRegFactor     - Regularisation factor                     [0 0 0]
% MeanRegDecFactor  - Start with more regularisation            [3 0 0] (log10)
% IterMax           - Total maximum number of iterations        [15 15 15]
% IterDec           - Number of decreasing reg iterations       [10 10 0]
% IterMin           - Total minimum number of iterations        [15 15 15]
% Tolerance         - Convergence threshold                     [1E-3 1E-3 1E-3]
% Verbose           - -1  = quiet
%                     [0] = moderate
%                      1  = verbose
%                      2  = plot mean image and log-likelihood
%                      3  = plot coil-wise fit (slow + disables parallel processing)
% Threads           - Number of parallel workers                [Inf=all]
%
% OUTPUT
% ------
% sens - Sensitivity maps   - (File)Array [Nx Ny Nz Nc]
% mean - Mean image         - (File)Array [Nx Ny Nz]
% prec - Noise precision    -       Array [Nc Nc]
%
% Nch = number of coils
% Nx  = Phase encode 1 
% Ny  = Phase encode 2 /or/ Slice
% Nz  = Frequency readout
% Nct = number of contrasts
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Parse input
% -------------------------------------------------------------------------
p  = inputParser;
p.FunctionName = 'b1m.fit';
p.addRequired('CoilImages',                           @utils.isarray);
p.addParameter('Precision',         NaN,              @isnumeric);
p.addParameter('SensRegFactor',     [1E5 1E4 1E4],    @(X) isnumeric(X) && numel(X) == 3);
p.addParameter('SensRegDecFactor',  [3 0 0],          @(X) isnumeric(X) && numel(X) == 3);
p.addParameter('MeanRegFactor',     [0 0 0],          @(X) isnumeric(X) && numel(X) == 3);
p.addParameter('MeanRegDecFactor',  [3 0 0],          @(X) isnumeric(X) && numel(X) == 3);
p.addParameter('VoxelSize',         [1 1 1],          @(X) isnumeric(X) && isrow(X) && numel(X) <= 3);
p.addParameter('Threads',           Inf,              @utils.isboolean);
p.addParameter('Tolerance',         [1E-3 1E-3 1E-3], @(X) isnumeric(X) && numel(X) == 3);
p.addParameter('IterMax',           [15 15 20],       @(X) utils.isintval(X) && numel(X) == 3);
p.addParameter('IterMin',           [15 15 20],       @(X) utils.isintval(X) && numel(X) == 3);
p.addParameter('IterDec',           [10 0 0],        @(X) utils.isintval(X) && numel(X) == 3);
p.addParameter('MaskBackground',    false,            @(X) utils.isarray(X) || utils.isboolean(X));
p.addParameter('Verbose',           0,                @utils.isboolean);
p.addParameter('SamplingMask',      [],               @utils.isarray);
p.parse(varargin{:});
coils           = p.Results.CoilImages;
prec            = p.Results.Precision;
sensfactor      = p.Results.SensRegFactor;
sensdecfactor   = p.Results.SensRegDecFactor;
meanfactor      = p.Results.MeanRegFactor;
meandecfactor   = p.Results.MeanRegDecFactor;
vs              = p.Results.VoxelSize;
tol             = p.Results.Tolerance;
itermax         = p.Results.IterMax;
itermin         = p.Results.IterMin;
iterdec         = p.Results.IterDec;
verbose         = p.Results.Verbose;
nbthreads       = p.Results.Threads;
mask            = p.Results.SamplingMask;
maskbg          = p.Results.MaskBackground;

% -------------------------------------------------------------------------
% Run
% -------------------------------------------------------------------------

% --- Estimate noise precision
if isempty(prec) || any(any(isnan(prec)))
    if verbose >= 1, fprintf([repmat('=', [1 80]) '\n']); end
    [prec, Nac] = b1m.init.noise(coils, mask);
    prec = prec / Nac;
end

% --- First pass: log-encoding
if verbose >= 1, fprintf([repmat('=', [1 80]) '\n']); end
if verbose >= 0, fprintf('#1: Fit mean + log-sensitivities\n'); end
[~,    meanim] = b1m.fit(coils, ...
                         'SensLog',             true, ...
                         'Precision',           prec, ...
                         'SamplingMask',        mask, ...
                         'MaskBackground',      maskbg, ...
                         'VoxelSize',           vs, ...
                         'SensRegFactor',       sensfactor(1), ...
                         'SensRegDecFactor',    sensdecfactor(1), ...
                         'MeanRegFactor',       meanfactor(1), ...
                         'MeanregDecFactor',    meandecfactor(1), ...
                         'Tolerance',           tol(1), ...
                         'IterMax',             itermax(1), ...
                         'IterMin',             itermin(1), ...
                         'IterDec',             iterdec(1), ...
                         'Verbose',             0, ...
                         'Threads',             nbthreads);
                         
% --- Second pass: fixed mean
if verbose >= 1, fprintf([repmat('=', [1 80]) '\n']); end
if verbose >= 0, fprintf('#2: Fit sensitivities\n'); end
[sens, meanim] = b1m.fit(coils, ...
                         'MeanImage',           meanim, ...
                         'MeanOptim',           false, ...
                         'SensLog',             false, ...
                         'Precision',           prec, ...
                         'SamplingMask',        mask, ...
                         'MaskBackground',      maskbg, ...
                         'VoxelSize',           vs, ...
                         'SensRegFactor',       sensfactor(2), ...
                         'SensRegDecFactor',    sensdecfactor(2), ...
                         'MeanRegFactor',       meanfactor(2), ...
                         'MeanregDecFactor',    meandecfactor(2), ...
                         'Tolerance',           tol(2), ...
                         'IterMax',             itermax(2), ...
                         'IterMin',             itermin(2), ...
                         'IterDec',             iterdec(2), ...
                         'Verbose',             2, ...
                         'Threads',             nbthreads);
                         
% --- Third pass: refine
if verbose >= 1, fprintf([repmat('=', [1 80]) '\n']); end
if verbose >= 0, fprintf('#3: Refine mean image + sensitivities\n'); end
[sens, meanim] = b1m.fit(coils, ...
                         'MeanImage',           meanim, ...
                         'SensMaps',            sens, ...
                         'SensLog',             false, ...
                         'Precision',           prec, ...
                         'SamplingMask',        mask, ...
                         'MaskBackground',      maskbg, ...
                         'VoxelSize',           vs, ...
                         'SensRegFactor',       sensfactor(3), ...
                         'SensRegDecFactor',    sensdecfactor(3), ...
                         'MeanRegFactor',       meanfactor(3), ...
                         'MeanregDecFactor',    meandecfactor(3), ...
                         'Tolerance',           tol(3), ...
                         'IterMax',             itermax(3), ...
                         'IterMin',             itermin(3), ...
                         'IterDec',             iterdec(3), ...
                         'Verbose',             2, ...
                         'Threads',             nbthreads);
                     
if verbose >= 1, fprintf([repmat('=', [1 80]) '\n']); end