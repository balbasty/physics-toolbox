function [sens,llm,llp,ok,ls] = sensitivity(varargin)
% Maximum a posteriori sensitivity profiles given a set of observed coil 
% images, a mean image and a noise precision (= inverse covariance) 
% matrix.
%
% FORMAT [sens, ...] = multicoil_sensitivity(mean, coils, sens, ...)
%
% REQUIRED
% --------
% mean  - (File)Array [Nx Ny Nz]    - (Complex) mean image
% coils - (File)Array [Nx Ny Nz Nc] - (Complex) coil images
% sens  - (File)Array [Nx Ny Nz Nc] - (Complex) log-sensitivity maps
%
% Nc = number of coils
% Nx = Phase encode 1 
% Ny = Phase encode 2 /or/ Slice
% Nz = Frequency readout
%
% KEYWORDS
% --------
% Index         - Array         - Indices of coils to update    [1:Nc]
% Precision     - Array [Nc Nc] - Noise precision matrix        [eye(Nc)]
% RegFactor     - Array [Nc 2]  - Regularisation factor         [1]
% VoxelSize     - Array [1 3]   - Voxel size                    [1]
% Log           - Boolean       - Log-encoding                  [false]
% LLPrior       - Array [Nc 2]  - Previous prior log-likelihood [NaN]
% SamplingMask  - Array [Nx Ny] - Mask of the sampling scheme   [ones]
%
% OUTPUT
% ------
% sens                - Updated (complex) log-sensitivity maps
% llm  - Scalar       - Matching part of the log-likelihood (all coils)
% llp  - Array [Nc 2] - Prior part of the log-likelihood (per coil/part)
% ok   - Array [1 Nc] - True if a better value was found
% 
% The optimum is found numerically using complex MM-Newton optimisation.
% The inverse problem is real and is solved by full multigrid.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Regularisation structure
prm = [0 0 1];              % Bending energy
spm_field('boundary', 1);   % Neumann boundary conditions
fmg = [2 2];                % Full multigrid parameters

% -------------------------------------------------------------------------
% Number of coils
Nc = size(varargin{1}, 4);

% =========================================================================
%
%                       PARSE AND PROCESS ARGUMENTS
%
% =========================================================================

% -------------------------------------------------------------------------
% Parse input
p = inputParser;
p.FunctionName = 'b1m.update.sensitivity';
p.addRequired('MeanImage',                   @utils.isarray);
p.addRequired('CoilImages',                  @utils.isarray);
p.addRequired('SensMaps',                    @utils.isarray);
p.addParameter('Index',         [],          @(X) isnumeric(X) && isrow(X));
p.addParameter('Precision',     1,           @isnumeric);
p.addParameter('RegFactor',     1,           @(X) isnumeric(X) && size(X,1) <= Nc && size(X,2) <= 2);
p.addParameter('VoxelSize',     [1 1 1],     @(X) isnumeric(X) && isrow(X) && numel(X) <= 3);
p.addParameter('Log',           false,       @utils.isboolean);
p.addParameter('LLPrior',       NaN,         @(X) isnumeric(X) && size(X,1) <= Nc && size(X,2) <= 2);
p.addParameter('SamplingMask',  [],          @utils.isarray);
p.addParameter('BackgroundMask',[],          @utils.isarray);
p.parse(varargin{:});
meanim      = p.Results.MeanImage;
coils       = p.Results.CoilImages;
sens        = p.Results.SensMaps;
all_n       = p.Results.Index;
prec        = p.Results.Precision;
reg         = p.Results.RegFactor;
senslog     = p.Results.Log;
vs          = p.Results.VoxelSize;
llp         = p.Results.LLPrior;
mask        = p.Results.SamplingMask;
bgmask      = p.Results.BackgroundMask;

iscplx      = ~isreal(coils);

% -------------------------------------------------------------------------
% Store dimensions
dim  = [size(coils) 1 1 1];
Nc   = dim(4);              % Number of coils
Nct  = dim(5);              % Number of contrasts
lat  = dim(1:3);            % Fully sampled lattice
Nvox = prod(lat);           % Number of (fully sampled) voxels
Nx   = lat(1);              % Phase encode 1
Ny   = lat(2);              % Phase encode 2 /OR/ slice
Nz   = lat(3);              % Frequency readout

% -------------------------------------------------------------------------
% Coils to process: default = all + ensure row-vector
if isempty(all_n), all_n = 1:Nc; end

% -------------------------------------------------------------------------
% Precision: default = identity
if numel(prec) == 1, prec = prec * eye(Nc); end
diagprec = isdiag(prec);

% -------------------------------------------------------------------------
% Coil factor: copy value for all coils if needed
reg = utils.pad(reg, [Nc-size(reg,1) 2-size(reg,2)], 'replicate', 'post');
llp = utils.pad(llp, [Nc-size(llp,1) 0], 'replicate', 'post');
llp = utils.pad(llp, [0 2-size(llp,2)], 0, 'post');

% -------------------------------------------------------------------------
% Voxel size: ensure row vector + complete
vs = utils.pad(vs, [0 3-numel(vs)], 'replicate', 'post');

% -------------------------------------------------------------------------
% Hessian (hopefully) majorising factor: sampling proportion
if ~isempty(mask)
    hfactor = sum(mask(:))/numel(mask);
else
    hfactor = 1;
end

% =========================================================================
%
%                NESTED FUNCTION FOR FASTER LINE-SEARCH
%
% =========================================================================

% -------------------------------------------------------------------------
% Prepare stuff to save time in the loop
function llm = computellm(n,ds)
    if diagprec
        subload_n = n;
    else
        subload_n = 1:size(coils,4);
    end
    A11  = prec(subload_n,subload_n);
    Nx11 = size(coils,1);
    Ny11 = size(coils,2);
    Nz11 = size(coils,3);
    Nc11 = numel(subload_n);
    Nvox11 = Nx11*Ny11*Nz11;
    llm = 0;
    % ---------------------------------------------------------------------
    % Compute log-likelihood slice-wise to save memory
    for zz=1:Nz11

        % -----------------------------------------------------------------
        % Load one slice of the complete coil dataset
        xzz = single(coils(:,:,zz,subload_n));
        xzz = reshape(xzz, [], numel(subload_n));
        xzz = xzz * A11;

        % -----------------------------------------------------------------
        % Load one slice of the mean
        rzz = single(meanim(:,:,zz,:));
        rzz = reshape(rzz, [], 1);

        % -----------------------------------------------------------------
        % Load one slice of the delta sensitivity
        dszz = single(ds(:,:,zz,:));
        dszz = reshape(dszz, [], 1);
        
        % -----------------------------------------------------------------
        % Load one slice of the complete double dataset + correct
        szz = single(sens(:,:,zz,subload_n));
        szz = reshape(szz, [], Nc11);
        if diagprec, szz      = szz      - dszz;
        else,        szz(:,n) = szz(:,n) - dszz; end
        if senslog, szz = exp(szz); end
        mzz = bsxfun(@times, rzz, szz);

        % -----------------------------------------------------------------
        % If incomplete sampling: push+pull coil-specific means
        pmzz = b1m.adjoint_forward(reshape(mzz, [Nx11 Ny11 1 Nc11]),mask);
        pmzz = reshape(pmzz, [], Nc11);
        pmzz = pmzz * A11;
        
        % -----------------------------------------------------------------
        % Compute log-likelihood
        llm = llm + double(real(mzz(:)' * pmzz(:)) - 2*real(mzz(:)' * xzz(:)));
    end % < loop z
    llm = -0.5 * Nvox11 * llm;
end % < function computellm


% -------------------------------------------------------------------------
% Compute initial log-likelihood (prior)
llpnan = any(~isfinite(llp),2);
if any(llpnan)
    llp(llpnan,:) = b1m.ll.sensitivity(sens(:,:,:,llpnan), ...
                        'RegFactor', reg(llpnan,:), 'VoxelSize', vs);
end
    
% -------------------------------------------------------------------------
% For each coil
llm = zeros(1,Nc);
ok  = zeros(1,Nc);
for n=all_n
    
    % =====================================================================
    %
    %                      COMPUTE GRADIENT AND HESSIAN
    %
    % =====================================================================
    
    if diagprec, load_n = n;
    else,        load_n = 1:Nc; end
    A1  = prec(load_n,load_n);
    Nc1 = numel(load_n);
    
    % ---------------------------------------------------------------------
    % Allocate conjugate gradient and Hessian
    g   = zeros([lat (1+iscplx)], 'single');
    H   = zeros(lat, 'single');
    llm(n) = 0;
    
    % ---------------------------------------------------------------------
    % Compute gradient slice-wise to save memory
    for z=1:Nz
        
        % -----------------------------------------------------------------
        % Load one slice of the complete coil dataset
        xz = single(coils(:,:,z,load_n));
        xz = reshape(xz, [], numel(load_n));
        xz = xz * A1;

        % -----------------------------------------------------------------
        % Load one slice of the mean
        rz = single(meanim(:,:,z,:));
        rz = reshape(rz, [], 1);

        % -----------------------------------------------------------------
        % Load one slice of the complete sensitivity dataset + correct
        sz = single(sens(:,:,z,load_n));
        sz = reshape(sz, [], Nc1);
        if senslog, sz = exp(sz); end
        mz = bsxfun(@times, rz, sz);
        
        % -----------------------------------------------------------------
        % If incomplete sampling: push+pull coil-specific means
        pmz = b1m.adjoint_forward(reshape(mz, [Nx Ny 1 Nc1]), mask);
        pmz = reshape(pmz, [], Nc1);
        pmz = pmz * A1;
        
        % -----------------------------------------------------------------
        % Compute log-likelihood
        llm(n) = llm + double(real(mz(:)' * pmz(:)) - 2*real(mz(:)' * xz(:)));
        
        % -----------------------------------------------------------------
        % Compute gradient
        gz = pmz - xz;
        if senslog
            if diagprec, gz = Nvox * conj(mz) .* gz;
            else,        gz = Nvox * conj(mz(:,n)) .* gz(:,n); end
        else
            if diagprec, gz = Nvox * conj(rz) .* gz;
            else,        gz = Nvox * conj(rz) .* gz(:,n); end
        end
        
        % -----------------------------------------------------------------
        % Compute Hessian
        if senslog
            if diagprec, Hz  = Nvox * hfactor * A1 * real(conj(mz) .* mz);
            else,        Hz  = Nvox * hfactor * prec(n,n) * real(conj(mz(:,n)) .* mz(:,n)); end
        else
                         Hz  = Nvox * hfactor * prec(n,n) * real(conj(rz) .* rz);
        end
        
        % -----------------------------------------------------------------
        % Missing data
        gz(~isfinite(gz)) = 0;
        Hz(~isfinite(Hz)) = 0;
        
        % -----------------------------------------------------------------
        % Background
        if ~isempty(bgmask)
            bgz = bgmask(:,:,z);
            bgz = reshape(bgz, [], 1);
            gz(bgz) = 0;
            Hz(bgz) = 0;
        end
        
        % -----------------------------------------------------------------
        % Store slice
        if iscplx, gz = cat(2, real(gz), imag(gz)); end
        g(:,:,z,:) = reshape(gz, Nx, Ny, 1, []);
        H(:,:,z)   = reshape(Hz, Nx, Ny, 1);
        
    end % < loop z
    clear sz xz rz mz pmz gz Hz
    llm(n) = -0.5 * Nvox * llm(n);
   
    % =====================================================================
    %
    %                      COMPUTE GAUSS-NEWTON STEP
    %
    % =====================================================================

    % ---------------------------------------------------------------------
    % Gradient: add prior term
    s0  = single(gather(sens(:,:,:,n)));
    s0r = real(s0);
    s0i = imag(s0);
    clear sens0
    if reg(n,1) > 0
        g(:,:,:,1) = g(:,:,:,1) + spm_field('vel2mom', s0r, [vs prm], reg(n,1));
    end
    if iscplx && reg(n,2) > 0
        g(:,:,:,2) = g(:,:,:,2) + spm_field('vel2mom', s0i, [vs prm], reg(n,2));
    end
    
    % ---------------------------------------------------------------------
    % Gauss-Newton
    if reg(n,1) > 0
        ds = spm_field(H, g(:,:,:,1), [vs prm fmg], reg(n,1));
    else
        ds = g(:,:,:,1)./max(H,eps('single'));
    end
    if iscplx
        if reg(n,2) > 0
            dsi = spm_field(H, g(:,:,:,2), [vs prm fmg], reg(n,2));
        else
            dsi = g(:,:,:,2)./max(H,eps('single'));
        end
    end
    clear g H 

    % ---------------------------------------------------------------------
    % Parts for log-likelihood (prior)
    llpr_part1 = 0;
    llpr_part2 = 0;
    llpi_part1 = 0;
    llpi_part2 = 0;
    if reg(n,1) > 0
        Lds = spm_field('vel2mom', ds, [vs prm], reg(n,1));
        llpr_part1 = double(s0r(:)' * Lds(:));
        llpr_part2 = double(ds(:)' * Lds(:));
        clear s0r Lds
    end
    if iscplx && reg(n,2) > 0
        Lds = spm_field('vel2mom', dsi, [vs prm], reg(n,2));
        llpi_part1 = double(s0i(:)' * Lds(:));
        llpi_part2 = double(dsi(:)' * Lds(:));
        clear s0i Lds
    end
    
    % ---------------------------------------------------------------------
    % Convert to complex
    if iscplx, ds = complex(ds,dsi); clear dsi; end
    
    
    % =====================================================================
    %
    %                             LINE-SEARCH
    %
    % =====================================================================
    
    % ---------------------------------------------------------------------
    % Save previous values
    llm0   = llm(n);
    llp0   = llp(n,:);
    ok(n)  = false;
    armijo = 1;
    for ls=1:12

        % -----------------------------------------------------------------
        % Compute log-likelihood (prior)
        llp(n,1) = 0;
        if reg(n,1) > 0
            llp(n,1) = -0.5 * (armijo^2 * llpr_part2 - 2 * armijo * llpr_part1);
            llp(n,1) = llp0(1) + llp(n,1);
        end
        if iscplx
            llp(n,2) = 0;
            if reg(n,2)
                llp(n,2) = -0.5 * (armijo^2 * llpi_part2 - 2 * armijo * llpi_part1);
                llp(n,2) = llp0(2) + llp(n,2);
            end
        end

        % -----------------------------------------------------------------
        % Compute log-likelihood (conditional)
        llm(n) = computellm(n, armijo*ds);

        % -----------------------------------------------------------------
        % Check progress
        if (llm(n)+sum(llp(n,:))) >= (llm0+sum(llp0))
            ok(n) = true;
            break
        else
%             if armijo < 0.5*hfactor
%                 % alpha = hfactor should have improved the o.f.
%                 % if not, break out
%                 % (we leave a bit of room for precision issues, though)
%                 break
%             end
            if ls < 11
                % no need to try alpha < hfactor
                armijo = max(armijo/2, 0.5*hfactor);
            else
                % if last ls iteration: use alpha = hfactor
                armijo = 0.5*hfactor;
            end
        end

    end % < loop ls
    
            
    % =====================================================================
    %
    %                             WRITE OUTPUT
    %
    % =====================================================================
    
    % ---------------------------------------------------------------------
    % Save
    if ok(n)
        sens(:,:,:,n) = sens(:,:,:,n) - armijo * gather(ds);
    else
        llm(n)   = llm0;
        llp(n,:) = llp0;
    end
        
end % < loop n

if diagprec
    llm = sum(llm);
else
    llm = llm(end);
end

end % < function multicoil_sensitivity