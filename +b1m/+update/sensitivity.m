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
% Index         - Array | Scalar   - Indices of coils to update    [1:Nc]
% Precision     - Array [Nc Nc]    - Noise precision matrix        [eye(Nc)]
% RegStructure  - [Abs Mem Ben]    - Regularisation structure      [0 0 1]
% RegCoilFactor - Vector [Nc]      - Reg modulator / coil          [1]
% RegPartFactor - [Mag Phase]      - Reg modulator / component     [1 1]
% RegBoundary   - Scalar | String  - Boundary condition            ['Neumann'] 
% VoxelSize     - Vector [3]       - Voxel size                    [1 1 1]
% LLPrior       - Scalar           - Previous prior log-likelihood [NaN]
% SensOptim     - [Mag Phase]      - Optimise magnitude/phase      [true true]
% SamplingMask  - Array [Nx Ny]    - Mask of the sampling scheme   [ones]
%
% Encoding      - String           - Sensitivity encoding   'freq'/['image']
% NbBasis       - [Mx My Mz]       - Number of DCT bases           [10 10 10]
% RegMatrix     - [NbBasis NbBasis]- Precomputed precision         []
%
% OUTPUT
% ------
% sens - Updated (complex) log-sensitivity maps
% llm  - Matching part of the log-likelihood (all coils)
% llp  - Prior part of the log-likelihood (all coils)
% ok   - True if a better value was found
% 
% The optimum is found numerically using complex Gauss-Newton optimisation.
% The inverse problem is real and is solved by full multigrid.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

%__________________________________________________________________________
% Development notes / Yael / 8 Nov 2018 
%
% This file is a bit complicated, as it tries to deal with various
% parameterisations of the sensitivity fields:
% - Sensitivity fields can be encoded directly in image space or in 
%   frequency space using discrete cosine basis functions.
% - It is possible to update only one of the (log)-field components, using
%   the `SensOptim` option. This also complicates stuff a bit.
%__________________________________________________________________________

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
p.addParameter('Index',         [],          @isnumeric);
p.addParameter('Precision',     1,           @isnumeric);
p.addParameter('RegStructure',  [0 0 1],     @(X) isnumeric(X) && numel(X) == 3);
p.addParameter('RegCoilFactor', 1,           @isnumeric);
p.addParameter('RegPartFactor', 1,           @(X) isnumeric(X) && numel(X) <= 2);
p.addParameter('RegBoundary',   1,           @utils.isboundary);
p.addParameter('VoxelSize',     [1 1 1],     @(X) isnumeric(X) && numel(X) <= 3);
p.addParameter('LLPrior',       NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('SensOptim',     [true true], @(X) (isnumeric(X) || islogical(X)) && numel(X) == 2);
p.addParameter('Encoding',      'image',     @ischar);
p.addParameter('NbBasis',       [10 10 10],  @(X) isnumeric(X) && numel(X) <= 3);
p.addParameter('RegMatrix',     [],          @(X) isnumeric(X));
p.addParameter('SamplingMask',  [],          @utils.isarray);
p.parse(varargin{:});
meanim      = p.Results.MeanImage;
coils       = p.Results.CoilImages;
sens        = p.Results.SensMaps;
all_n       = p.Results.Index;
prec        = p.Results.Precision;
prm         = p.Results.RegStructure;
coilfactor  = p.Results.RegCoilFactor;
partfactor  = p.Results.RegPartFactor;
bnd         = p.Results.RegBoundary;
vs          = p.Results.VoxelSize;
llp         = p.Results.LLPrior;
optim       = p.Results.SensOptim;
encoding    = p.Results.Encoding;
nbasis      = p.Results.NbBasis;
regmatrix   = p.Results.RegMatrix;
mask        = p.Results.SamplingMask;

% -------------------------------------------------------------------------
% Post-process input
Nc   = size(coils,4);                                     % Number of coils
lat  = [size(coils,1) size(coils,2) size(coils,3)]; % Fully sampled lattice
Nvox = prod(lat);                        % Number of (fully sampled) voxels
Nx   = lat(1);                                             % Phase encode 1
Ny   = lat(2);                                  % Phase encode 2 /OR/ slice
Nz   = lat(3);                                          % Frequency readout
% Coils to process: default = all + ensure row-vector
if isempty(all_n)
    all_n = 1:Nc;
end
all_n = all_n(:).';
% Precision: default = identity
if numel(prec) == 1
    prec = prec * eye(Nc);
end
diagprec = isdiag(prec);
% Reg components: just change reg structure
partfactor = padarray(partfactor(:)', [0 max(0,2-numel(partfactor))], 'replicate', 'post');
% Reg factor: ensure zero sum -> propagate their sum to reg components
coilfactor = padarray(coilfactor(:), [max(0,Nc-numel(coilfactor)) 0], 'replicate', 'post');
partfactor = partfactor * sum(coilfactor);
coilfactor = coilfactor/sum(coilfactor);
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
% Voxel size: ensure row vector + complete
vs = padarray(vs(:)', [0 max(0,3-numel(vs))], 'replicate', 'post');
if strcmpi(encoding, 'frequency')
    % Nb basis: ensure row vector + complete
    nbasis = padarray(nbasis(:)', [0 max(0,3-numel(nbasis))], 'replicate', 'post');
    % Create regularisation matrix
    if isempty(regmatrix)
        regmatrix = spm_bias_lib('regulariser', prm, lat, nbasis, vs);
    end
    % Create basis functions (B)
    [B1,B2,B3] = spm_bias_lib('dcbasis', lat, nbasis);
end
% Optimisation: if observed images are real, optim = [true false]
if utils.isrealarray(coils)
    optim(2) = false;
end
optim = logical(optim);
if all(~optim)
    warning('Nothing to update')
    return
end

% -------------------------------------------------------------------------
% GPU
gpu_on = isa(prec, 'gpuArray');
if gpu_on, loadarray = @utils.loadarray_gpu;
else,      loadarray = @utils.loadarray_cpu; end
 
% -------------------------------------------------------------------------
% Boundary condition (usually Neumann = null derivative)
spm_field('boundary', bnd); 

% -------------------------------------------------------------------------
% Undersampling proportion (for approximate Hessian)
if ~isempty(mask)
    propmask = sum(mask(:))/numel(mask);
else
    propmask = 1;
end

% =========================================================================
%
%                NESTED FUNCTION FOR FASTER LINE-SEARCH
%
% =========================================================================

% -------------------------------------------------------------------------
% Prepare stuff to save time in the loop
% --- Log-likelihood
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
        sz  = exp(sz);
        rz  = bsxfun(@times, rz, sz);
        clear sz

        % -----------------------------------------------------------------
        % If incomplete sampling: push+pull coil-specific means
        prz = b1m.adjoint_forward(reshape(rz, [Nx Ny 1 Nc1]),mask);
        prz = reshape(prz, [], Nc1);
        
        % -----------------------------------------------------------------
        % Compute log-likelihood
        llm = llm - sum(double(real(dot(rz,(prz-2*xz)*A1,1))), 'omitnan');
        
        rz  = [];
        prz = [];
        xz    = [];
    end % < loop z
    llm = 0.5 * Nvox * llm;
end % < function computellm


% -------------------------------------------------------------------------
% Compute log-likelihood (prior)
if isnan(llp)
    llp = b1m.ll.sensitivity(sens, prm, partfactor, coilfactor, bnd, optim, vs);
end
    
% -------------------------------------------------------------------------
% For each coil
for n=all_n
    
    % =====================================================================
    %
    %                      COMPUTE GRADIENT AND HESSIAN
    %
    % =====================================================================
    
    if diagprec
        load_n = n;
    else
        load_n = 1:Nc;
    end
    A1  = prec(load_n,load_n);
    Nc1 = numel(load_n);
    
    % ---------------------------------------------------------------------
    % Allocate conjugate gradient and Hessian
    switch lower(encoding)
        case 'image'
            g   = zeros([lat sum(optim)], 'like', loadarray(single(1)));
            H   = zeros(lat, 'like', loadarray(single(1)));
        case 'frequency'
            g   = zeros([nbasis sum(optim)], 'like', loadarray(single(1)));
            H   = zeros([nbasis nbasis], 'like', loadarray(single(1)));
    end
    llm = 0;
    
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
        sz = exp(sz);
        rz = bsxfun(@times, rz, sz);
        clear sz
        
        % -----------------------------------------------------------------
        % If incomplete sampling: push+pull coil-specific means
        prz = b1m.adjoint_forward(reshape(rz, [Nx Ny 1 Nc1]), mask);
        prz = reshape(prz, [], Nc1);
        
        % -----------------------------------------------------------------
        % Compute gradient
        llm = llm - sum(double(real(dot(rz,(prz-2*xz)*A1,1))), 'omitnan');
        
        tmp = (prz - xz) * A1;
        if diagprec
            tmp = Nvox * conj(rz) .* tmp;
        else
            tmp = Nvox * conj(rz(:,n)) .* tmp(:,n);
        end
        
        gz = zeros([size(tmp,1) sum(optim)], 'like', real(tmp(1)));
        i  = 1;
        if optim(1) % If optimise sensitivity magnitude
            gz(:,i) = real(tmp);
            i = i+1;
        end
        if optim(2) % If optimise sensitivity phase
            gz(:,i) = imag(tmp);
        end
        
        if diagprec
            Hz  = Nvox * propmask * A1 * real(conj(rz) .* rz);
        else
            Hz  = Nvox * propmask * prec(n,n) * real(conj(rz(:,n)) .* rz(:,n));
        end
        
        gz(~isfinite(gz)) = 0;
        Hz(~isfinite(gz)) = 0;
        
        switch lower(encoding)
            case 'image'
                g(:,:,z,:)  = reshape(gz, Nx, Ny, 1, []); clear gz
                H(:,:,z)    = reshape(Hz, Nx, Ny, 1);     clear Hz
            case 'frequency'
                b3z = B3(z,:);
                gz  = reshape(gz, Nx, Ny, 1, []);
                gz  = gather(gz); % dct only defined for CPU arrays
                gz  = dct(gz, [], 1);
                gz  = gz(1:nbasis(1),:,:,:);
                gz  = dct(gz, [], 2);
                gz  = gz(:,1:nbasis(2),:,:);
                gz  = loadarray(gz); % Back to GPU
                gz  = bsxfun(@times, gz, reshape(b3z, 1, 1, []));
                g   = g + gz; clear gz
                Hz  = reshape(Hz, Nx, Ny, 1);
                Hz  = gather(Hz); % spm_krutil only defined for CPU arrays
                Hz  = kron(b3z'*b3z,spm_krutil(double(Hz),B1,B2,1));
                H   = H + reshape(Hz, [nbasis nbasis]); clear Hz
        end
        
        clear tmp xz rz prz
    end % < loop z
    llm = 0.5 * Nvox * llm;
   
    % =====================================================================
    %
    %                      COMPUTE GAUSS-NEWTON STEP
    %
    % =====================================================================
    
    switch lower(encoding)
        case 'image'
            % -------------------------------------------------------------
            % Gather gradient & Hessian (if on GPU)
            g = gather(g);
            H = gather(H);
    
            % -------------------------------------------------------------
            % Gradient: add prior term
            if all(optim)
                sens0 = single(sens(:,:,:,n));
                sens0 = cat(4, real(sens0), imag(sens0));
            elseif optim(1)
                sens0 = real(single(sens(:,:,:,n)));
            elseif optim(2)
                sens0 = imag(single(sens(:,:,:,n)));
            end
            g  = g + spm_field('vel2mom', single(gather(sens0)), [vs coilfactor(n)*prm], partfactor(optim));

            % -------------------------------------------------------------
            % Gauss-Newton
            ds = zeros(size(sens0), 'single');
            i = 1;
            if optim(1)
                ds(:,:,:,i) = spm_field(H, g(:,:,:,i), [vs coilfactor(n)*prm 2 2], partfactor(1));
                i = i + 1;
            end
            if optim(2)
                ds(:,:,:,i) = spm_field(H, g(:,:,:,i), [vs coilfactor(n)*prm 2 2], partfactor(2));
            end
            clear g H

            % -------------------------------------------------------------
            % Parts for log-likelihood (prior)
            Lds = spm_field('vel2mom', ds, [vs coilfactor(n)*prm], partfactor(optim));
            llp_part1 = coilfactor(n) * double(reshape(sens0, 1, [])) * double(reshape(Lds, [], 1));
            llp_part2 = coilfactor(n) * double(reshape(ds, 1, [])) * double(reshape(Lds, [], 1));
            clear s0 Lds
    
        case 'frequency'
            
            % -------------------------------------------------------------
            % Convert to vector/matrix
            g = reshape(g, [], sum(optim));
            g = double(g); 
            
            H = reshape(H, prod(nbasis), prod(nbasis));
            H = double(H);
            H = H + 1e-7 * max(diag(H)) * eye(size(H));
            
            % -------------------------------------------------------------
            % Gradient: add prior term
            if all(optim)
                sens0 = single(sens(:,:,:,n));
                sens0 = cat(4, real(sens0), imag(sens0));
            elseif optim(1)
                sens0 = real(single(sens(:,:,:,n)));
            elseif optim(2)
                sens0 = imag(single(sens(:,:,:,n)));
            end
            sens0 = dct(sens0,[],1);
            sens0 = sens0(1:nbasis(1),:,:,:);
            sens0 = dct(sens0,[],2);
            sens0 = sens0(:,1:nbasis(2),:,:);
            sens0 = dct(sens0,[],3);
            sens0 = sens0(:,:,1:nbasis(3),:);
            sens0 = reshape(sens0, [], sum(optim));
            sens0 = loadarray(sens0, @double); % Load on GPU
            
            i = 1;
            if optim(1)
                g(:,i) = g(:,i) + coilfactor(n) * partfactor(1) * regmatrix * sens0(:,i);
                i = i + 1;
            end
            if optim(2)
                g(:,i) = g(:,i) + coilfactor(n) * partfactor(2) * regmatrix * sens0(:,i);
            end
            
            % -------------------------------------------------------------
            % Gauss-Newton
            ds = zeros(size(sens0), 'like', gpuArray(double(1)));
            i  = 1;
            if optim(1)
                ds(:,i) = (H + coilfactor(n) * partfactor(1) * regmatrix)\g(:,i);
                i = i + 1;
            end
            if optim(2)
                ds(:,i) = (H + coilfactor(n) * partfactor(2) * regmatrix)\g(:,i);
            end
            clear g H
            ds = reshape(ds, [nbasis sum(optim)]);
            
            % -------------------------------------------------------------
            % Parts for log-likelihood (prior)
            llp_part1 = 0;
            llp_part2 = 0;
            i = 1;
            if optim(1)
                Lds = coilfactor(n) * partfactor(1) * regmatrix * reshape(ds(:,:,:,i), [], 1);
                llp_part1 = llp_part1 + sens0(:,i)' * Lds;
                llp_part2 = llp_part2 + Lds' * Lds;
                i = i + 1;
            end
            if optim(2)
                Lds = coilfactor(n) * partfactor(2) * regmatrix * reshape(ds(:,:,:,i), [], 1);
                llp_part1 = llp_part1 + sens0(:,i)' * Lds;
                llp_part2 = llp_part2 + Lds' * Lds;
            end
            clear s0 Lds
            
            % -------------------------------------------------------------
            % Convert to image representation
            ds = single(ds);
            % ds = idct(idct(idct(ds,lat(1),1),lat(2),2),lat(3),3);
            ds = reshape(B1 * reshape(ds, nbasis(1), []), lat(1), nbasis(2), nbasis(3), []);
            ds = permute(ds, [2 3 1 4]);
            ds = reshape(B2 * reshape(ds, nbasis(2), []), lat(2), nbasis(3), lat(1), []);
            ds = permute(ds, [2 3 1 4]);
            ds = reshape(B2 * reshape(ds, nbasis(3), []), lat(3), lat(1), lat(2), []);
            ds = permute(ds, [2 3 1 4]);
    end
    
    if all(optim)
        ds = complex(ds(:,:,:,1),ds(:,:,:,2));
    elseif optim(2)
        ds = complex(0,ds);
    end
            
    % =====================================================================
    %
    %                             LINE-SEARCH
    %
    % =====================================================================
    
    % ---------------------------------------------------------------------
    % Save previous values
    llm0   = llm;
    llp0   = llp;
    ok     = false;
    armijo = 1;
    for ls=1:12

        % -----------------------------------------------------------------
        % Compute log-likelihood (prior)
        llp = -0.5 * (armijo^2 * llp_part2 - 2 * armijo * llp_part1);
        llp = llp0 + llp;

        % -----------------------------------------------------------------
        % Compute log-likelihood (conditional)
        llm = computellm(n, armijo*ds);

        % -----------------------------------------------------------------
        % Check progress
        if (llm+llp) >= (llm0+llp0)
            ok = true;
            break
        else
            armijo = armijo/2;
        end

    end % < loop ls
    
            
    % =====================================================================
    %
    %                             WRITE OUTPUT
    %
    % =====================================================================
    
    % ---------------------------------------------------------------------
    % Save
    if ok
        sens(:,:,:,n) = sens(:,:,:,n) - armijo * gather(ds);
    else
        llm = llm0;
        llp = llp0;
    end
        
end % < loop n

end % < function multicoil_sensitivity