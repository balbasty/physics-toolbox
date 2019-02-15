function [rho,llm,llp,ok,ls] = multicoil_mean_map(varargin)
% Maximum a posteriori mean given a set of observed coil images, 
% log-sensitivity profiles and a noise precision (= inverse covariance) 
% matrix.
%
% FORMAT [mean, ...] = multicoil_mean_map(coils, sens, (mean), ...)
%
% REQUIRED
% --------
% coils  - (File)Array [Nx Ny Nz Nc] - Complex coil images
% sens   - (File)Array [Nx Ny Nz Nc] - Complex log-sensitivity profiles
%
% OPTIONAL
% --------
% mean   - (File)Array [Nx Ny Nz]    - Complex mean image
%
% KEYWORD ARGUMENTS
% ----------------- 
% Precision     - Array [Nc Nc]    - Noise precision matrix        [eye(Nc)]
% RegStructure  - [Abs Mem Ben]    - Regularisation structure      [0 1 0]
% RegFactor     - Scalar           - Regularisation magnitude      [0=ML]
% RegPartFactor - [Mag Phase]      - Reg modulator magnitude/phase [1 1]
% RegBoundary   - Scalar | String  - Boundary condition            ['Neumann']
% VoxelSize     - Vector [3]       - Voxel size                    [1 1 1]
% LLPrior       - Scalar           - Previous prior log-likelihood [NaN]
% PartOptim     - [Mag Phase]      - Optimise magnitude/phase      [true true]
% SamplingMask  - Array [Nx Ny]    - Mask of the sampling scheme   [ones]
%
% OUTPUT
% ------
% mean - Updated mean
% llm  - Log-likelihood of the matching term
% llp  - Log-likelihood of the mean prior term
% ok   - Did we improve?
% ls   - Number of line-search steps used
%
% Nc = Number of coils
% Nx = Phase encode 1 
% Ny = Phase encode 2 /or/ Slice
% Nz = Frequency readout
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% =========================================================================
%
%                       PARSE AND PROCESS ARGUMENTS
%
% =========================================================================

% -------------------------------------------------------------------------
% Helper functions to check input arguments
% -------------------------------------------------------------------------
function ok = isarray(X)
    ok = isnumeric(X) || islogical(X) || isa(X, 'file_array');
end
function ok = isboundary(X)
    ok = (isnumeric(X) && isscalar(X) && 0 <= X && X <= 1) || ...
         (ischar(X)    && any(strcmpi(X, {'c','circulant','n','neumann'})));
end
function ok = isrealarray(X)
    function okk = isrealtype(T)
        okk = numel(T) > 7 || ~strcmpi(T(1:7),'complex');
    end
    if isa(X, 'file_array')
        ok = all(cellfun(@isrealtype, {X.dtype}));
    else
        ok = isreal(X);
    end
end

% -------------------------------------------------------------------------
% Parse input
% -------------------------------------------------------------------------
p = inputParser;
p.FunctionName = 'multicoil_mean_iter';
p.addRequired('CoilImages',                  @isarray);
p.addRequired('SensMaps',                    @isarray);
p.addOptional('MeanImage',      [],          @isarray);
p.addParameter('Precision',     1,           @isnumeric);
p.addParameter('RegStructure',  [0 1 0],     @(X) isnumeric(X) && numel(X) == 3);
p.addParameter('RegFactor',     0,           @(X) isnumeric(X) && isscalar(X));
p.addParameter('RegPartFactor', 1,           @(X) isnumeric(X) && numel(X) <= 2);
p.addParameter('RegBoundary',   1,           @isboundary);
p.addParameter('VoxelSize',     [1 1 1],     @(X) isnumeric(X) && numel(X) <= 3);
p.addParameter('LLPrior',       NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('PartOptim',     [true true], @(X) (isnumeric(X) || islogical(X)) && numel(X) == 2);
p.addParameter('SamplingMask',  1,           @isarray);
p.parse(varargin{:});
x           = p.Results.CoilImages;
s           = p.Results.SensMaps;
rho         = p.Results.MeanImage;
A           = p.Results.Precision;
prm         = p.Results.RegStructure;
regfactor   = p.Results.RegFactor;
regpart     = p.Results.RegPartFactor;
bnd         = p.Results.RegBoundary;
vs          = p.Results.VoxelSize;
llp         = p.Results.LLPrior;
optim       = p.Results.PartOptim;
mask        = p.Results.SamplingMask;

% -------------------------------------------------------------------------
% Post-process input
% -------------------------------------------------------------------------
Nc   = size(x,4);                        % Number of coils
lat  = [size(x,1) size(x,2) size(x,3)];  % Fully sampled lattice
Nvox = prod(lat);                        % Number of (fully sampled) voxels
% Precision: default = identity
if numel(A) == 1
    A = A * eye(Nc);
end
% Reg components: just change reg structure
regpart = padarray(regpart(:)', [0 max(0,2-numel(regpart))], 'replicate', 'post');
% Boundary: convert to scalar representation
switch bnd
    case {0, 'c', 'circulant'}
        bnd = 0;
    case {1, 'n', 'neumann'}
        bnd = 1;
    otherwise
        warning('Unknown boundary condition %s. Using Neumann instead', num2str(bnd))
        bnd = 1;
end
if regfactor > 0, spm_field('boundary', bnd); end
% Voxel size: ensure row vector + complete
vs = padarray(vs(:)', [0 max(0,3-numel(vs))], 'replicate', 'post');
% Optimisation: if observed images are real, optim = [true false]
% Reg components: just change reg structure
optim = padarray(optim(:)', [0 max(0,2-numel(optim))], 'replicate', 'post');
if isrealarray(x)
    optim(2) = false;
end
optim = logical(optim);
if all(~optim)
    warning('Nothing to update')
    return
end
% Correct regularisation based on the mean of the signal
if regfactor > 0
    regfactor = regfactor * double(mean(abs(x(:))))^2;
end

gpu_on = isa(A, 'gpuArray');
if gpu_on, loadarray = @loadarray_gpu;
else,      loadarray = @loadarray_cpu; end


% =========================================================================
%
%                           COMPUTE DERIVATIVES
%
% =========================================================================

% -------------------------------------------------------------------------
% Compute log-likelihood (prior)
% -------------------------------------------------------------------------
if isnan(llp) 
    if regfactor
        llp = multicoil_ll_mean_prior(rho, regfactor*prm, vs, regpart);
    else
        llp = 0;
    end
end

% -------------------------------------------------------------------------
% Allocate gradient and Hessian
% -------------------------------------------------------------------------
g   = zeros([lat sum(optim)],'like',loadarray(single(1)));
H   = zeros([lat 1],'like',loadarray(single(1)));
llm = 0;

if ~isempty(mask)
    propmask = sum(mask(:))/numel(mask);
else
    propmask = 1;
end

% -------------------------------------------------------------------------
% Compute conditional part (slice-wise to save memory)
% -------------------------------------------------------------------------
for z=1:lat(3)
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    xz = loadarray(x(:,:,z,:), @single);
    xz = reshape(xz, [], Nc);

    % ---------------------------------------------------------------------
    % Load one slice of the (previous) mean
    rhoz = loadarray(rho(:,:,z,:), @single);
    rhoz = reshape(rhoz, [], 1);

    % ---------------------------------------------------------------------
    % Load one slice of the complete sensitivity dataset + correct
    sz   = loadarray(s(:,:,z,:), @single);
    sz   = reshape(sz, [], Nc);
    sz   = exp(sz);
    rhoz = bsxfun(@times, rhoz, sz);

    % ---------------------------------------------------------------------
    % If incomplete sampling: push+pull coil-specific means
    prhoz = multicoil_pushpullwrap(reshape(rhoz, [lat(1:2) 1 Nc]),mask);
    prhoz = reshape(prhoz, [], Nc);
        
    llm  = llm - sum(double(real(dot(rhoz,(prhoz-2*xz)*A,1))));
    rhoz = [];
    
    % ---------------------------------------------------------------------
    % Compute Hessian
    tmp = propmask * Nvox * real(dot(sz, sz*A, 2));
    H(:,:,z) = reshape(tmp, lat(1:2));
    
    % ---------------------------------------------------------------------
    % Compute co-gradient
    tmp = Nvox * dot(sz, (prhoz - xz)*A, 2);
    gz  = zeros([size(tmp,1) sum(optim)], 'like' ,real(tmp(1)));
    i   = 1;
    if optim(1) % If optimise sensitivity magnitude
        gz(:,i) = real(tmp);
        i = i+1;
    end
    if optim(2) % If optimise sensitivity phase
        gz(:,i) = imag(tmp);
    end
    g(:,:,z,:)  = reshape(gz, lat(1), lat(2), 1, []);

    
    gz    = [];
    tmp   = [];
    prhoz = [];
end
llm = 0.5 * Nvox * llm;


if regfactor > 0
    g = gather(g);
    H = gather(H);
    
    % ---------------------------------------------------------------------
    % Load previous value
    % ---------------------------------------------------------------------
    if all(optim)
        rho0 = zeros([lat 2], 'single');
        rho0(:,:,:,1) = real(single(rho));
        rho0(:,:,:,2) = imag(single(rho));
    elseif optim(1)
        rho0 = real(single(rho));
    elseif optim(2)
        rho0 = imag(single(rho));
    end
    
    % ---------------------------------------------------------------------
    % Compute prior part
    % ---------------------------------------------------------------------
    g = g + spm_field('vel2mom', rho0, [vs regfactor*prm], regpart(optim));
end

% =========================================================================
%
%                           LINE SEARCH
%
% =========================================================================

% -------------------------------------------------------------------------
% Gauss-Newton
% -------------------------------------------------------------------------
if regfactor > 0
    drho = zeros(size(rho0), 'single');
    i = 1;
    if optim(1)
        drho(:,:,:,i) = spm_field(H, g(:,:,:,i), [vs regfactor*prm 2 2], regpart(1));
        i = i + 1;
    end
    if optim(2)
        drho(:,:,:,i) = spm_field(H, g(:,:,:,i), [vs regfactor*prm 2 2], regpart(2));
    end
else
    drho = bsxfun(@rdivide, g, H);
end
clear g H

% -------------------------------------------------------------------------
% Parts for log-likelihood (prior)
% -------------------------------------------------------------------------
if regfactor
    Ldrho = spm_field('vel2mom', drho, [vs regfactor*prm], regpart(optim));
    llp_part1 = double(reshape(rho0, 1, [])) * double(reshape(Ldrho, [], 1));
    llp_part2 = double(reshape(drho, 1, [])) * double(reshape(Ldrho, [], 1));
    clear Lds
else
    llp_part1 = 0;
    llp_part2 = 0;
end
clear rho0

if all(optim)
    drho = complex(drho(:,:,:,1),drho(:,:,:,2));
elseif optim(2)
    drho = complex(0,drho);
end
% drho = multicoil_pushpullwrap(drho,mask);

% -------------------------------------------------------------------------
% Line search
% -------------------------------------------------------------------------
llm0   = llm;
llp0   = llp;
armijo = 1;
ok     = false;
for ls=1:6
    
    % ---------------------------------------------------------------------
    % Compute log-likelihood (prior)
    llp = -0.5 * (armijo^2 * llp_part2 - 2 * armijo * llp_part1);
    llp = llp0 + llp;

    % ---------------------------------------------------------------------
    % Compute log-likelihood (conditional)
    llm = 0;
    for z=1:lat(3)
    
        % -----------------------------------------------------------------
        % Load one slice of the complete coil dataset
        xz = loadarray(x(:,:,z,:), @single);
        xz = reshape(xz, [], Nc);

        % -----------------------------------------------------------------
        % Load one slice of the (previous) mean
        rhoz = loadarray(rho(:,:,z,:), @single);
        rhoz = reshape(rhoz, [], 1);
        rhoz = rhoz - armijo*reshape(double(drho(:,:,z,:)), [], 1);

        % -----------------------------------------------------------------
        % Load one slice of the complete sensitivity dataset + correct
        sz   = loadarray(s(:,:,z,:), @single);
        sz   = reshape(sz, [], Nc);
        sz   = exp(sz);
        rhoz = bsxfun(@times, rhoz, sz);

        % -----------------------------------------------------------------
        % If incomplete sampling: push+pull coil-specific means
        prhoz = multicoil_pushpullwrap(reshape(rhoz, [lat(1:2) 1 Nc]), mask);
        prhoz = reshape(prhoz, [], Nc);

        llm = llm - sum(double(real(dot(rhoz,(prhoz-2*xz)*A,1))));
        rhoz  = [];
        prhoz = [];

    end
    llm = 0.5 * Nvox * llm;


    % ---------------------------------------------------------------------
    % Check progress
    if (llm+llp) >= (llm0+llp0)
        ok = true;
        break
    else
        armijo = armijo/2;
    end
    
end

% -------------------------------------------------------------------------
% If line-search failed: roll back
% -------------------------------------------------------------------------
if ~ok
    llm = llm0;
    llp = llp0;
end

% =========================================================================
%
%                                 SAVE
%
% =========================================================================

if ok
    if ~isa(rho, 'gpuArray')
        drho = gather(drho);
    end
    rho(:,:,:) = rho(:,:,:) - armijo * drho;
end

end