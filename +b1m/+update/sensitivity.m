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
fmg = [2 2];

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
if isempty(all_n)
    all_n = 1:Nc;
end

% -------------------------------------------------------------------------
% Precision: default = identity
if numel(prec) == 1
    prec = prec * eye(Nc);
end
diagprec = isdiag(prec);

% -------------------------------------------------------------------------
% Coil factor: copy value for all coils if needed
reg = utils.pad(reg, [Nc-size(reg,1) 2-size(reg,2)], 'replicate', 'post');
llp = utils.pad(llp, [Nc-size(llp,1) 2-size(llp,2)], 'replicate', 'post');

% -------------------------------------------------------------------------
% Voxel size: ensure row vector + complete
vs = utils.pad(vs, [0 3-numel(vs)], 'replicate', 'post');

% -------------------------------------------------------------------------
% GPU
gpu_on = isa(prec, 'gpuArray');
if gpu_on, loadarray = @utils.loadarray_gpu;
else,      loadarray = @utils.loadarray_cpu; end

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
        load_n = n;
    else
        load_n = 1:Nc;
    end
    A1  = prec(load_n,load_n);
    Nc1 = numel(load_n);
    llm = 0;
    % ---------------------------------------------------------------------
    % Compute log-likelihood slice-wise to save memory
    for z=1:Nz

        % -----------------------------------------------------------------
        % Load one slice of the complete coil dataset
        xz = loadarray(coils(:,:,z,load_n), @single);
        xz = reshape(xz, [], numel(load_n));

        % -----------------------------------------------------------------
        % Load one slice of the mean
        rz = loadarray(meanim(:,:,z,:), @single);
        rz = reshape(rz, [], 1);

        % -----------------------------------------------------------------
        % Load one slice of the delta sensitivity
        dsz = loadarray(ds(:,:,z,:), @single);
        dsz = reshape(dsz, [], 1);
        
        % -----------------------------------------------------------------
        % Load one slice of the complete double dataset + correct
        sz = loadarray(sens(:,:,z,load_n), @single);
        sz = reshape(sz, [], Nc1);
        if diagprec
            sz = sz - dsz;
        else
            sz(:,n) = sz(:,n) - dsz;
        end
        dsz = [];
        if senslog
            sz = exp(sz);
        end
        mz  = bsxfun(@times, rz, sz);

        % -----------------------------------------------------------------
        % If incomplete sampling: push+pull coil-specific means
        pmz = b1m.adjoint_forward(reshape(mz, [Nx Ny 1 Nc1]),mask);
        pmz = reshape(pmz, [], Nc1);
        
        % -----------------------------------------------------------------
        % Compute log-likelihood
        llm = llm + sum(double(real(conj(reshape(mz,[],1)).*reshape(pmz*A1,[],1))), 'omitnan', 'double') ...
                  - 2*sum(double(real(conj(reshape(mz,[],1)).*reshape(xz*A1,[],1))), 'omitnan', 'double');
        
        mz  = [];
        pmz = [];
        xz  = [];
    end % < loop z
    llm = -0.5 * Nvox * llm;
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
    g   = zeros([lat 2], 'like', loadarray(single(1)));
    H   = zeros(lat, 'like', loadarray(single(1)));
    llm(n) = 0;
    
    % ---------------------------------------------------------------------
    % Compute gradient slice-wise to save memory
    for z=1:Nz
        
        % -----------------------------------------------------------------
        % Load one slice of the complete coil dataset
        xz = loadarray(coils(:,:,z,load_n), @single);
        xz = reshape(xz, [], numel(load_n));

        % -----------------------------------------------------------------
        % Load one slice of the mean
        rz = loadarray(meanim(:,:,z,:), @single);
        rz = reshape(rz, [], 1);

        % -----------------------------------------------------------------
        % Load one slice of the complete sensitivity dataset + correct
        sz = loadarray(sens(:,:,z,load_n), @single);
        sz = reshape(sz, [], Nc1);
        if senslog
            sz = exp(sz);
        end
        mz = bsxfun(@times, rz, sz);
        
        % -----------------------------------------------------------------
        % If incomplete sampling: push+pull coil-specific means
        pmz = b1m.adjoint_forward(reshape(mz, [Nx Ny 1 Nc1]), mask);
        pmz = reshape(pmz, [], Nc1);
        
        % -----------------------------------------------------------------
        % Compute log-likelihood
        llm(n) = llm(n) + sum(double(real(conj(reshape(mz,[],1)).*reshape(pmz*A1,[],1))), 'omitnan', 'double') ...
                        - 2*sum(double(real(conj(reshape(mz,[],1)).*reshape(xz*A1,[],1))), 'omitnan', 'double');
        
        % -----------------------------------------------------------------
        % Compute gradient
        gz = (pmz - xz) * A1;
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
            bgz = loadarray(bgmask(:,:,z), @logical);
            bgz = reshape(bgz, [], 1);
            gz(bgz) = 0;
            Hz(bgz) = 0;
        end
        
        % -----------------------------------------------------------------
        % Store slice
        gz = cat(2, real(gz), imag(gz));
        g(:,:,z,:)  = reshape(gz, Nx, Ny, 1, []);
        H(:,:,z)    = reshape(Hz, Nx, Ny, 1);
        
    end % < loop z
    clear sz xz rz mz pmz gz Hz
    llm(n) = -0.5 * Nvox * llm(n);
   
    % =====================================================================
    %
    %                      COMPUTE GAUSS-NEWTON STEP
    %
    % =====================================================================
    
    % ---------------------------------------------------------------------
    % Gather gradient & Hessian (if on GPU)
    g = gather(g);
    H = gather(H);

    % ---------------------------------------------------------------------
    % Gradient: add prior term
    s0  = single(gather(sens(:,:,:,n)));
    s0r = real(s0);
    s0i = imag(s0);
    clear sens0
    if reg(n,1) > 0
        g(:,:,:,1) = g(:,:,:,1) + spm_field('vel2mom', s0r, [vs prm], reg(n,1));
    end
    if reg(n,2) > 0
        g(:,:,:,2) = g(:,:,:,2) + spm_field('vel2mom', s0i, [vs prm], reg(n,2));
    end
    
    % ---------------------------------------------------------------------
    % Gauss-Newton
    if reg(n,1) > 0
        dsr = spm_field(H, g(:,:,:,1), [vs prm fmg], reg(n,1));
    else
        dsr = g(:,:,:,1)./max(H,eps('single'));
    end
    if reg(n,2) > 0
        dsi = spm_field(H, g(:,:,:,2), [vs prm fmg], reg(n,2));
    else
        dsi = g(:,:,:,2)./max(H,eps('single'));
    end
    clear g H 

    % ---------------------------------------------------------------------
    % Parts for log-likelihood (prior)
    llpr_part1 = 0;
    llpr_part2 = 0;
    llpi_part1 = 0;
    llpi_part2 = 0;
    if reg(n,1) > 0
        Lds = spm_field('vel2mom', dsr, [vs prm], reg(n,1));
        llpr_part1 = sum(s0r(:) .* Lds(:), 'double');
        llpr_part2 = sum(dsr(:) .* Lds(:), 'double');
        clear s0r Lds
    end
    if reg(n,2) > 0
        Lds = spm_field('vel2mom', dsi, [vs prm], reg(n,2));
        llpi_part1 = sum(reshape(s0i, 1, []) .* reshape(Lds, [], 1), 'double');
        llpi_part2 = sum(reshape(dsi, 1, []) .* reshape(Lds, [], 1), 'double');
        clear s0i Lds
    end
    
    % ---------------------------------------------------------------------
    % Convert to complex
    ds = complex(dsr,dsi);
    clear dsr
    
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
        llp(n,2) = 0;
        if reg(n,2)
            llp(n,2) = -0.5 * (armijo^2 * llpi_part2 - 2 * armijo * llpi_part1);
            llp(n,2) = llp0(2) + llp(n,2);
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