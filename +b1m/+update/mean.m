function [meanim,llm,llp,ok,ls] = mean(varargin)
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
% Precision      - Array [Nc Nc]    - Noise precision matrix        [eye(Nc)]
% RegFactor      - Scalar           - Regularisation magnitude      [0=ML]
% VoxelSize      - Array [1 3]      - Voxel size                    [1]
% LLPrior        - Scalar           - Previous prior log-likelihood [NaN]
% SensLog        - Boolean          - Log-sensitivities             [false]
% SamplingMask   - Array [Nx Ny]    - Mask of the sampling scheme   [ones]
% BackgroundMask - Array [Nx Ny]    - Mask of the background        [zeros]
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
% Parse input
% -------------------------------------------------------------------------
p = inputParser;
p.FunctionName = 'b1m.update.mean';
p.addRequired('CoilImages',                   @utils.isarray);
p.addRequired('SensMaps',                     @utils.isarray);
p.addOptional('MeanImage',       [],          @utils.isarray);
p.addParameter('Precision',      1,           @isnumeric);
p.addParameter('RegFactor',      0,           @(X) isnumeric(X) && isscalar(X));
p.addParameter('VoxelSize',      [1 1 1],     @(X) isnumeric(X) && isrow(X) && numel(X) <= 3);
p.addParameter('LLPrior',        NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('SensLog',        false,       @utils.isboolean);
p.addParameter('SamplingMask',   [],          @utils.isarray);
p.addParameter('BackgroundMask', [],          @utils.isarray);
p.parse(varargin{:});
coils       = p.Results.CoilImages;
sens        = p.Results.SensMaps;
meanim      = p.Results.MeanImage;
prec        = p.Results.Precision;
reg         = p.Results.RegFactor;
vs          = p.Results.VoxelSize;
llp         = p.Results.LLPrior;
mask        = p.Results.SamplingMask;
senslog     = p.Results.SensLog;
bgmask      = p.Results.BackgroundMask;

% -------------------------------------------------------------------------
% Dimensions
% -------------------------------------------------------------------------
Nc   = size(coils,4);                                     % Number of coils
lat  = [size(coils,1) size(coils,2) size(coils,3)]; % Fully sampled lattice
Nvox = prod(lat);                        % Number of (fully sampled) voxels

% -------------------------------------------------------------------------
% Precision: default = identity

% -------------------------------------------------------------------------
if numel(prec) == 1
    prec = prec * eye(Nc);
end

% -------------------------------------------------------------------------
% Pad Voxel size
% -------------------------------------------------------------------------
vs = utils.pad(double(vs), [0 3-numel(vs)], 'replicate', 'post');

% -------------------------------------------------------------------------
% Regularisation structure
% -------------------------------------------------------------------------
if reg > 0
    prm = [0 1 0];              % Membrane energy
    spm_field('boundary', 1);   % Neumann boundary conditions
end

% =========================================================================
%
%                           COMPUTE DERIVATIVES
%
% =========================================================================

% -------------------------------------------------------------------------
% Compute log-likelihood (mean prior)
% -------------------------------------------------------------------------
if isnan(llp) 
    if reg
        llp = b1m.ll.mean(meanim, 'RegFactor', reg, 'VoxelSize', vs);
    else
        llp = 0;
    end
end

% -------------------------------------------------------------------------
% Allocate gradient and Hessian
% -------------------------------------------------------------------------
g   = zeros([lat 2],'single');                 % Gradient
H   = zeros([lat 1],'single');                 % Hessian
llm = 0;                                       % Conditional log-likelihood

% -------------------------------------------------------------------------
% Undersampling proportion
% -------------------------------------------------------------------------
if ~isempty(mask)
    hfactor = sum(mask(:))/numel(mask);
else
    hfactor = 1;
end

% -------------------------------------------------------------------------
% Missing data
% -------------------------------------------------------------------------
missing = any(isnan(coils(:)));

% -------------------------------------------------------------------------
% Compute conditional part (slice-wise to save memory)
% -------------------------------------------------------------------------
for z=1:lat(3)
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    xz = single(coils(:,:,z,:));
    xz = reshape(xz, [], Nc);
    
    % ---------------------------------------------------------------------
    % Load one slice of the (previous) mean
    rz = single(meanim(:,:,z,:));
    rz = reshape(rz, [], 1);

    % ---------------------------------------------------------------------
    % Load one slice of the complete sensitivity dataset + correct
    sz = single(sens(:,:,z,:));
    sz = reshape(sz, [], Nc);
    if senslog, sz = exp(sz); end
    mz = bsxfun(@times, rz, sz);
    
    % ---------------------------------------------------------------------
    % If incomplete sampling: push+pull coil-specific means
    pmz = b1m.adjoint_forward(reshape(mz, [lat(1:2) 1 Nc]), mask);
    pmz = reshape(pmz, [], Nc);
    
    if missing
    % ---------------------------------------------------------------------
    % Missing data in observed domain
        cz = utils.gmm.lib('obs2code', xz);
        code_list = unique(cz)';
        gz = zeros(size(xz,1), 1, 'single');
        Hz = zeros(size(xz,1), 1, 'single');
        for code=code_list
    
            sub = cz == code;
            bin  = utils.gmm.lib('code2bin', code, Nc);
            if ~any(bin)
                continue
            end
            % -------------------------------------------------------------
            % Extract sub + precompute
            A    = utils.invPD(prec);
            A    = A(bin,bin);
            A    = utils.invPD(A);
            mzs  = mz(sub,bin);
            xzs  = xz(sub,bin) * A;
            pmzs = pmz(sub,bin) * A;
            szs  = sz(sub,bin);
            % -------------------------------------------------------------
            % Compute Hessian
            Hz(sub) = hfactor * Nvox * real(dot(szs, szs*A, 2));
            % -------------------------------------------------------------
            % Compute co-gradient
            gz(sub) = Nvox * dot(szs, pmzs, 2);
            % -------------------------------------------------------------
            % Compute log-likelihood
            llm = llm + double(real(mzs(:)' * pmzs(:)) - 2*real(mzs(:)' * xzs(:)));
        end
    else
        % -----------------------------------------------------------------
        % Precompute
        xz  = xz * prec;
        pmz = pmz * prec;
        % -----------------------------------------------------------------
        % Compute Hessian
        Hz = hfactor * Nvox * real(dot(sz, sz*prec, 2));
        % -----------------------------------------------------------------
        % Compute co-gradient
        gz = Nvox * dot(sz, pmz - xz, 2);
        % -----------------------------------------------------------------
        % Compute log-likelihood
        llm = llm + double(real(mz(:)' * pmz(:)) - 2*real(mz(:)' * xz(:)));
    end
        
    % ---------------------------------------------------------------------
    % Background mask
    if ~isempty(bgmask)
        bgz = bgmask(:,:,z);
        bgz = reshape(bgz, [], 1);
        gz(bgz) = 0;
        Hz(bgz) = 0;
    end
    gz(~isfinite(gz)) = 0;
    Hz(~isfinite(Hz)) = 0;
    
    % ---------------------------------------------------------------------
    % Store slice
    gz         = cat(2, real(gz), imag(gz));
    g(:,:,z,:) = reshape(gz, lat(1), lat(2), 1, []);
    H(:,:,z)   = reshape(Hz, lat(1:2));
end
clear rz mz pmz sz gz Hz
llm = - 0.5 * Nvox * llm;


if reg > 0
    g = gather(g);
    H = gather(H);
    
    % ---------------------------------------------------------------------
    % Load previous value
    % ---------------------------------------------------------------------
    mean0          = zeros([lat 2], 'single');
    mean0(:,:,:,1) = real(single(meanim));
    mean0(:,:,:,2) = imag(single(meanim));
    
    % ---------------------------------------------------------------------
    % Compute prior part
    % ---------------------------------------------------------------------
    g = g + spm_field('vel2mom', mean0, [vs reg*prm]);
end

% =========================================================================
%
%                           LINE SEARCH
%
% =========================================================================

% -------------------------------------------------------------------------
% Gauss-Newton
% -------------------------------------------------------------------------
if reg > 0
    dmean = zeros(size(mean0), 'single');
    dmean(:,:,:,1) = spm_field(H, g(:,:,:,1), [vs reg*prm 2 2]);
    dmean(:,:,:,2) = spm_field(H, g(:,:,:,2), [vs reg*prm 2 2]);
else
    dmean = bsxfun(@rdivide, g, H);
end
clear g H

% -------------------------------------------------------------------------
% Parts for log-likelihood (prior)
% -------------------------------------------------------------------------
if reg
    Ldrho = spm_field('vel2mom', dmean, [vs reg*prm]);
    llp_part1 = double(mean0(:)' * Ldrho(:));
    llp_part2 = double(dmean(:)' * Ldrho(:));
    clear Lds
else
    llp_part1 = 0;
    llp_part2 = 0;
end
clear rho0

dmean = complex(dmean(:,:,:,1),dmean(:,:,:,2));

    
% -------------------------------------------------------------------------
% Background mask
% -------------------------------------------------------------------------
if ~isempty(bgmask)
    dmean(bgmask) = 0;
end

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
        xz = single(coils(:,:,z,:));
        xz = reshape(xz, [], Nc);
    
        % -----------------------------------------------------------------
        % Load one slice of the (previous) mean
        rz = single(meanim(:,:,z,:));
        rz = reshape(rz, [], 1);
        rz = rz - armijo*reshape(double(dmean(:,:,z,:)), [], 1);

        % -----------------------------------------------------------------
        % Load one slice of the complete sensitivity dataset + correct
        sz = single(sens(:,:,z,:));
        sz = reshape(sz, [], Nc);
        if senslog, sz = exp(sz); end
        mz = bsxfun(@times, rz, sz);
        
        % -----------------------------------------------------------------
        % If incomplete sampling: push+pull coil-specific means
        pmz = b1m.adjoint_forward(reshape(mz, [lat(1:2) 1 Nc]), mask);
        pmz = reshape(pmz, [], Nc);

        % -----------------------------------------------------------------
        % Missing data in observed domain
        if missing
            cz = utils.gmm.lib('obs2code', xz);
            code_list = unique(cz)';
            for code=code_list

                sub = cz == code;
                bin  = utils.gmm.lib('code2bin', code, Nc);
                if ~any(bin)
                    continue
                end
                % ---------------------------------------------------------
                % Extract sub + precompute
                A    = utils.invPD(prec);
                A    = A(bin,bin);
                A    = utils.invPD(A);
                mzs  = mz(sub,bin);
                xzs  = xz(sub,bin) * A;
                pmzs = pmz(sub,bin) * A;
                % ---------------------------------------------------------
                % Compute log-likelihood
                llm = llm + double(real(mzs(:)' * pmzs(:)) - 2*real(mzs(:)' * xzs(:)));
            end
        else
            % -------------------------------------------------------------
            % Precompute
            xz  = xz * prec;
            pmz = pmz * prec;
            % -------------------------------------------------------------
            % Compute log-likelihood
            llm = llm + double(real(mz(:)' * pmz(:)) - 2*real(mz(:)' * xz(:)));
        end

    end
    clear xz rz sz mz pmz
    llm = - 0.5 * Nvox * llm;


    % ---------------------------------------------------------------------
    % Check progress
    if (llm+llp) >= (llm0+llp0)
        ok = true;
        break
    else
        if armijo == hfactor
            % alpha = propmask should have improved the o.f.
            % if not, break out
            break
        end
        if ls < 5
            % no need to try alpha < hfactor
            armijo = max(armijo/2, hfactor);
        else
            % if last ls iteration: use alpha = hfactor
            armijo = hfactor;
        end
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
    meanim(:,:,:) = meanim(:,:,:) - armijo * dmean;
end