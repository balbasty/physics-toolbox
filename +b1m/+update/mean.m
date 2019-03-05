function [mean,llm,llp,ok,ls] = mean(varargin)
% Maximum a posteriori mean given a set of observed coil images, 
% log-sensitivity profiles and a noise precision (= inverse covariance) 
% matrix.
%
% FORMAT [mean, ...] = b1m.update.mean(coils, sens, (mean), ...)
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
p.FunctionName = 'b1m.update.mean';
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
coils       = p.Results.CoilImages;
sens        = p.Results.SensMaps;
mean        = p.Results.MeanImage;
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
Nc   = size(coils,4);                                     % Number of coils
lat  = [size(coils,1) size(coils,2) size(coils,3)]; % Fully sampled lattice
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
if isrealarray(coils)
    optim(2) = false;
end
optim = logical(optim);
if all(~optim)
    warning('Nothing to update')
    return
end
% Correct regularisation based on the mean of the signal
if regfactor > 0
    regfactor = regfactor / double(mean(abs(coils(:))))^2;
end

gpu_on = isa(A, 'gpuArray');
if gpu_on, loadarray = @utils.loadarray_gpu;
else,      loadarray = @utils.loadarray_cpu; end


% =========================================================================
%
%                           COMPUTE DERIVATIVES
%
% =========================================================================

% -------------------------------------------------------------------------
% Compute log-likelihood (mean prior)
% -------------------------------------------------------------------------
if isnan(llp) 
    if regfactor
        llp = b1m.ll.mean(mean, regfactor*prm, vs, regpart);
    else
        llp = 0;
    end
end

% -------------------------------------------------------------------------
% Allocate gradient and Hessian
% -------------------------------------------------------------------------
g   = zeros([lat sum(optim)],'like',loadarray(single(1)));       % Gradient
H   = zeros([lat 1],'like',loadarray(single(1)));                 % Hessian
llm = 0;                                       % Conditional log-likelihood

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
    xz = loadarray(coils(:,:,z,:), @single);
    xz = reshape(xz, [], Nc);

    % ---------------------------------------------------------------------
    % Load one slice of the (previous) mean
    rz = loadarray(mean(:,:,z,:), @single);
    rz = reshape(rz, [], 1);

    % ---------------------------------------------------------------------
    % Load one slice of the complete sensitivity dataset + correct
    sz = loadarray(sens(:,:,z,:), @single);
    sz = reshape(sz, [], Nc);
    sz = exp(sz);
    rz = bsxfun(@times, rz, sz);

    % ---------------------------------------------------------------------
    % If incomplete sampling: push+pull coil-specific means
    prz = b1m.adjoint_forward(reshape(rz, [lat(1:2) 1 Nc]), mask);
    prz = reshape(prz, [], Nc);
        
    llm = llm - sum(double(real(dot(rz,(prz-2*xz)*A,1))));
    clear rz
    
    % ---------------------------------------------------------------------
    % Compute Hessian
    Hz       = propmask * Nvox * real(dot(sz, sz*A, 2));
    H(:,:,z) = reshape(Hz, lat(1:2));
    clear Hz
    
    % ---------------------------------------------------------------------
    % Compute co-gradient
    tmp = Nvox * dot(sz, (prz - xz)*A, 2);
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

    clear gz tmp prz
end
llm = 0.5 * Nvox * llm;


if regfactor > 0
    g = gather(g);
    H = gather(H);
    
    % ---------------------------------------------------------------------
    % Load previous value
    % ---------------------------------------------------------------------
    if all(optim)
        mean0          = zeros([lat 2], 'single');
        mean0(:,:,:,1) = real(single(mean));
        mean0(:,:,:,2) = imag(single(mean));
    elseif optim(1)
        mean0 = real(single(mean));
    elseif optim(2)
        mean0 = imag(single(mean));
    end
    
    % ---------------------------------------------------------------------
    % Compute prior part
    % ---------------------------------------------------------------------
    g = g + spm_field('vel2mom', mean0, [vs regfactor*prm], regpart(optim));
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
    dmean = zeros(size(mean0), 'single');
    i = 1;
    if optim(1)
        dmean(:,:,:,i) = spm_field(H, g(:,:,:,i), [vs regfactor*prm 2 2], regpart(1));
        i = i + 1;
    end
    if optim(2)
        dmean(:,:,:,i) = spm_field(H, g(:,:,:,i), [vs regfactor*prm 2 2], regpart(2));
    end
else
    dmean = bsxfun(@rdivide, g, H);
end
clear g H

% -------------------------------------------------------------------------
% Parts for log-likelihood (prior)
% -------------------------------------------------------------------------
if regfactor
    Ldrho = spm_field('vel2mom', dmean, [vs regfactor*prm], regpart(optim));
    llp_part1 = double(reshape(mean0, 1, [])) * double(reshape(Ldrho, [], 1));
    llp_part2 = double(reshape(dmean, 1, [])) * double(reshape(Ldrho, [], 1));
    clear Lds
else
    llp_part1 = 0;
    llp_part2 = 0;
end
clear rho0

if all(optim)
    dmean = complex(dmean(:,:,:,1),dmean(:,:,:,2));
elseif optim(2)
    dmean = complex(0,dmean);
end
% Uncomment the line below to restrict the mean to lie inside the sampling
% pattern:
% dmean = b1m.adjoint_forward(dmean,mask);

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
        xz = loadarray(coils(:,:,z,:), @single);
        xz = reshape(xz, [], Nc);

        % -----------------------------------------------------------------
        % Load one slice of the (previous) mean
        rz = loadarray(mean(:,:,z,:), @single);
        rz = reshape(rz, [], 1);
        rz = rz - armijo*reshape(double(dmean(:,:,z,:)), [], 1);

        % -----------------------------------------------------------------
        % Load one slice of the complete sensitivity dataset + correct
        sz   = loadarray(sens(:,:,z,:), @single);
        sz   = reshape(sz, [], Nc);
        sz   = exp(sz);
        rz = bsxfun(@times, rz, sz);

        % -----------------------------------------------------------------
        % If incomplete sampling: push+pull coil-specific means
        prz = b1m.adjoint_forward(reshape(rz, [lat(1:2) 1 Nc]), mask);
        prz = reshape(prz, [], Nc);

        llm = llm - sum(double(real(dot(rz,(prz-2*xz)*A,1))));
        rz  = [];
        prz = [];

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
    if ~isa(mean, 'gpuArray')
        dmean = gather(dmean);
    end
    mean(:,:,:) = mean(:,:,:) - armijo * dmean;
end

end