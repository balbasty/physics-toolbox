function [s,llm,llp,ok,ls] = multicoil_sensitivity(varargin)
% Maximum a posteriori sensitivity profiles given a set of observed coil 
% images, a mean image and a noise precision (= inverse covariance) 
% matrix.
%
% FORMAT [sens, ...] = multicoil_sensitivity(mean, coils, sens, ...)
%
% REQUIRED
% --------
% mean  - (File)Array [Nx Ny Nz  1] - (Complex) mean image
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
% RegCompFactor - [Mag Phase]      - Reg modulator / component     [1 1]
% RegBoundary   - Scalar | String  - Boundary condition            ['Neumann'] 
% VoxelSize     - Vector [3]       - Voxel size                    [1 1 1]
% LLPrior       - Scalar           - Previous prior log-likelihood [NaN]
% SensOptim     - [Mag Phase]      - Optimise magnitude/phase      [true true]
% Parallel      - Scalar | Logical - Number of parallel workers    [false]
% SamplingMask  - Array [Nx Ny]    - Mask of the sampling scheme   [ones]
%
% Encoding      - String           - Sensitivity encoding  'image'/['frequency']
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
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

%__________________________________________________________________________
% Development notes / Yael / 8 Nov 2018 
%
% This file is a bit complicated, as it tries to deal with various
% parameterisations of the sensitivity fields:
% - It is possible to update only one of the (log)-field components, using
%   the `SensOptim` option. This also complicates stuff a bit.
%__________________________________________________________________________

% -------------------------------------------------------------------------
% Helper functions to check input arguments
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
p = inputParser;
p.FunctionName = 'multicoil_sensitivity';
p.addRequired('MeanImage',                   @isarray);
p.addRequired('CoilImages',                  @isarray);
p.addRequired('SensMaps',                    @isarray);
p.addParameter('Index',         [],          @isnumeric);
p.addParameter('Precision',     1,           @isnumeric);
p.addParameter('RegStructure',  [0 0 1],     @(X) isnumeric(X) && numel(X) == 3);
p.addParameter('RegCoilFactor', 1,           @isnumeric);
p.addParameter('RegCompFactor', 1,           @(X) isnumeric(X) && numel(X) <= 2);
p.addParameter('RegBoundary',   1,           @isboundary);
p.addParameter('VoxelSize',     [1 1 1],     @(X) isnumeric(X) && numel(X) <= 3);
p.addParameter('LLPrior',       NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('SensOptim',     [true true], @(X) (isnumeric(X) || islogical(X)) && numel(X) == 2);
p.addParameter('Parallel',      0,           @(X) (isnumeric(X) || islogical(X)) && isscalar(X));
p.addParameter('Encoding',      'image',     @ischar);
p.addParameter('NbBasis',       [10 10 10],  @(X) isnumeric(X) && numel(X) <= 3);
p.addParameter('RegMatrix',     [],          @(X) isnumeric(X));
p.addParameter('SamplingMask',  [],          @isarray);
p.parse(varargin{:});
rho         = p.Results.MeanImage;
x           = p.Results.CoilImages;
s           = p.Results.SensMaps;
all_n       = p.Results.Index;
A           = p.Results.Precision;
prm         = p.Results.RegStructure;
alpha       = p.Results.RegCoilFactor;
gamma       = p.Results.RegCompFactor;
bnd         = p.Results.RegBoundary;
vs          = p.Results.VoxelSize;
llp         = p.Results.LLPrior;
optim       = p.Results.SensOptim;
Nw          = p.Results.Parallel;
encoding    = p.Results.Encoding;
nbasis      = p.Results.NbBasis;
regmatrix   = p.Results.RegMatrix;
mask        = p.Results.SamplingMask;

% -------------------------------------------------------------------------
% Post-process input
N   = size(x,4);                        % Number of coils
lat = [size(x,1) size(x,2) size(x,3)];  % Fully sampled lattice
% Coils to process: default = all + ensure row-vector
if isempty(all_n)
    all_n = 1:N;
end
all_n = all_n(:).';
% Precision: default = identity
if numel(A) == 1
    A = A * eye(N);
end
diagprec = isdiag(A);
% Reg components: just change reg structure
gamma = padarray(gamma(:)', [0 max(0,2-numel(gamma))], 'replicate', 'post');
% Reg factor: ensure zero sum -> propagate their sum to reg components
alpha = padarray(alpha(:), [max(0,N-numel(alpha)) 0], 'replicate', 'post');
gamma = gamma * sum(alpha);
alpha = alpha/sum(alpha);
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
if isrealarray(x)
    optim(2) = false;
end
optim = logical(optim);
if all(~optim)
    warning('Nothing to update')
    return
end
% Parallel: convert to number of workers
if islogical(Nw)
    if Nw, Nw = inf;
    else,  Nw = 0;
    end
end
% if Nw > 0
%     warning('Parallel processing not implemented. Running sequential instead.')
%     Nw = 0;
% end

% -------------------------------------------------------------------------
% Boundary condition (usually Neumann = null derivative)
spm_field('boundary', bnd); 

% -------------------------------------------------------------------------
% Prepare stuff to save time in the loop
% --- GPU
gpu_on = isa(A, 'gpuArray');
if gpu_on, loadarray = @loadarray_gpu;
else,      loadarray = @loadarray_cpu; end
% --- Log-likelihood
function llm = computellm(n,ds)
    if diagprec
        load_n = n;
    else
        load_n = 1:N;
    end
    A1 = A(load_n,load_n);
    llm = 0;
    % ---------------------------------------------------------------------
    % Compute gradient slice-wise to save memory
    parfor(z=1:lat(3) , Nw) % < Uncomment for parallel processing
    % for z=1:lat(3)          % < Uncomment for sequential processing

        % -----------------------------------------------------------------
        % Enforce boundary condition -> needed with parfor
        spm_field('boundary', bnd); 

        % -----------------------------------------------------------------
        % Load one slice of the complete coil dataset
        xz = loadarray(x(:,:,z,load_n), @double);
        xz = reshape(xz, [], numel(load_n));

        % -----------------------------------------------------------------
        % Load one slice of the mean
        rhoz = loadarray(rho(:,:,z,:), @double);
        rhoz = reshape(rhoz, [], 1);

        % -----------------------------------------------------------------
        % Load one slice of the delta sensitivity
        dsz   = loadarray(ds(:,:,z,:), @double);
        dsz   = reshape(dsz, [], sum(optim));
        if all(optim)
            dsz = dsz(:,1) + 1i * dsz(:,2);
        elseif optim(2)
            dsz = 1i * dsz(:,2);
        end
        
        % -----------------------------------------------------------------
        % Load one slice of the complete double dataset + correct
        sz      = loadarray(s(:,:,z,load_n), @single);
        sz      = reshape(sz, [], numel(load_n));
        if diagprec
            sz = sz - dsz;
        else
            sz(:,n) = sz(:,n) - dsz;
        end
        dsz     = [];
        sz      = exp(sz);
        rhoz    = bsxfun(@times, rhoz, sz);
        sz      = []; % clear

        % -----------------------------------------------------------------
        % If incomplete sampling: push+pull coil-specific means
        prhoz = multicoil_pushpullwrap(reshape(rhoz, [lat(1:2) 1 numel(load_n)]),mask);
        prhoz = reshape(prhoz, [], numel(load_n));
        
        % -----------------------------------------------------------------
        % Compute log-likelihood
        llm = llm - 0.5 * prod(lat) * (...
                              sum(real(dot(rhoz,prhoz*A1,1))) ...
                            - 2*sum(real(dot(rhoz,xz*A1,1))));
        
        rhoz  = [];
        prhoz = [];
        xz    = [];
    end % < loop z
end % < function computellm
    
% -------------------------------------------------------------------------
% Compute log-likelihood (prior)
if isnan(llp)
    llp = multicoil_ll_prior(s, prm, gamma, alpha, bnd, optim, vs);
end
    
% -------------------------------------------------------------------------
% For each coil
for n=all_n
    
    if diagprec
        load_n = n;
    else
        load_n = 1:N;
    end
    A1 = A(load_n,load_n);
    
    % ---------------------------------------------------------------------
    % Allocate conjugate gradient and Hessian
    switch lower(encoding)
        case 'image'
            g   = zeros([lat sum(optim)], 'single');
            H   = zeros(lat, 'single');
        case 'frequency'
            g   = zeros([nbasis sum(optim)], 'single');
            H   = zeros([nbasis nbasis], 'single');
    end
    llm = 0;
    
    % ---------------------------------------------------------------------
    % Compute gradient slice-wise to save memory
    % parfor(z=1:lat(3) , Nw) % < Uncomment for parallel processing
    for z=1:lat(3)          % < Uncomment for sequential processing

        % -----------------------------------------------------------------
        % Enforce boundary condition -> needed with parfor
        spm_field('boundary', bnd); 
        
        % -----------------------------------------------------------------
        % Load one slice of the complete coil dataset
        xz = loadarray(x(:,:,z,load_n), @double);
        xz = reshape(xz, [], numel(load_n));

        % -----------------------------------------------------------------
        % Load one slice of the mean
        rhoz = loadarray(rho(:,:,z,:), @double);
        rhoz = reshape(rhoz, [], 1);

        % -----------------------------------------------------------------
        % Load one slice of the complete sensitivity dataset + correct
        sz   = loadarray(s(:,:,z,load_n), @double);
        sz   = reshape(sz, [], numel(load_n));
        sz   = exp(sz);
        rhoz = bsxfun(@times, rhoz, sz);
        sz   = []; % clear
        
        % -----------------------------------------------------------------
        % If incomplete sampling: push+pull coil-specific means
        prhoz = multicoil_pushpullwrap(reshape(rhoz, [lat(1:2) 1 numel(load_n)]),mask);
        prhoz = reshape(prhoz, [], numel(load_n));
        
        % -----------------------------------------------------------------
        % Compute gradient
        llm = llm - 0.5 * prod(lat) * (...
                              sum(real(dot(rhoz,prhoz*A1,1))) ...
                            - 2*sum(real(dot(rhoz,xz*A1,1))));
        
        tmp = (prhoz - xz) * A1;
        if diagprec
            tmp = prod(lat) * conj(rhoz) .* tmp;
        else
            tmp = prod(lat) * conj(rhoz(:,n)) .* tmp(:,n);
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
        
        % tmp = prhoz * A;
        % Hz  = prod(lat) * real(conj(rhoz(:,n)) .* tmp(:,n));
        if diagprec
            Hz  = prod(lat) * A1 * real(conj(rhoz) .* prhoz);
        else
            Hz  = prod(lat) * A(n,n) * real(conj(rhoz(:,n)) .* prhoz(:,n));
        end
        
        switch lower(encoding)
            case 'image'
                g(:,:,z,:)  = reshape(gz, lat(1), lat(2), 1, []);
                gz          = []; % clear
                H(:,:,z)    = reshape(Hz, lat(1), lat(2));
                Hz          = []; % clear
            case 'frequency'
                b3z = B3(z,:);
                gz  = reshape(gz, lat(1), lat(2), 1, []);
                gz  = dct(gz, [], 1);
                gz  = gz(1:nbasis(1),:,:,:);
                gz  = dct(gz, [], 2);
                gz  = gz(:,1:nbasis(2),:,:);
                gz  = bsxfun(@times, gz, reshape(b3z, 1, 1, []));
                g   = g + gz;
                gz  = []; % clear
                Hz  = reshape(Hz, lat(1), lat(2));
                Hz  = kron(b3z'*b3z,spm_krutil(double(Hz),B1,B2,1));
                H   = H + reshape(Hz, [nbasis nbasis]);
                Hz  = []; % clear
        end
        
        tmp   = []; % clear
        xz    = []; % clear
        rhoz  = []; % clear
        prhoz = []; % clear
    end % < loop z
    
    
    switch lower(encoding)
        case 'image'
            % -------------------------------------------------------------
            % Gather gradient & Hessian (if on GPU)
            g = gather(g);
            H = gather(H);
    
            % -------------------------------------------------------------
            % Gradient: add prior term
            if all(optim)
                s0 = zeros([lat 2], 'single');
                s1 = single(s(:,:,:,n));
                s0(:,:,:,1) = real(s1);
                s0(:,:,:,2) = imag(s1);
                clear s1
            elseif optim(1)
                s0 = real(single(s(:,:,:,n)));
            elseif optim(2)
                s0 = imag(single(s(:,:,:,n)));
            end
            g  = g + spm_field('vel2mom', s0, [vs alpha(n)*prm], gamma(optim));

            % -------------------------------------------------------------
            % Gauss-Newton
            ds = zeros(size(s0), 'single');
            i = 1;
            if optim(1)
                ds(:,:,:,i) = spm_field(H, g(:,:,:,i), [vs alpha(n)*prm 2 2], gamma(1));
                i = i + 1;
            end
            if optim(2)
                ds(:,:,:,i) = spm_field(H, g(:,:,:,i), [vs alpha(n)*prm 2 2], gamma(2));
            end
            clear g H

            % -------------------------------------------------------------
            % Parts for log-likelihood (prior)
            Lds = spm_field('vel2mom', ds, [vs alpha(n)*prm], gamma(optim));
            llp_part1 = alpha(n) * double(reshape(s0, 1, [])) * double(reshape(Lds, [], 1));
            llp_part2 = alpha(n) * double(reshape(ds, 1, [])) * double(reshape(Lds, [], 1));
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
                s0 = zeros([lat 2], 'single');
                s1 = single(s(:,:,:,n));
                s0(:,:,:,1) = real(s1);
                s0(:,:,:,2) = imag(s1);
                clear s1
            elseif optim(1)
                s0 = real(single(s(:,:,:,n)));
            elseif optim(2)
                s0 = imag(single(s(:,:,:,n)));
            end
            s0 = dct(s0,[],1);
            s0 = s0(1:nbasis(1),:,:,:);
            s0 = dct(s0,[],2);
            s0 = s0(:,1:nbasis(2),:,:);
            s0 = dct(s0,[],3);
            s0 = s0(:,:,1:nbasis(3),:);
            s0 = reshape(s0, [], sum(optim));
            s0 = double(s0);
            
            i = 1;
            if optim(1)
                g(:,i) = g(:,i) + alpha(n) * gamma(1) * regmatrix * s0(:,i);
                i = i + 1;
            end
            if optim(2)
                g(:,i) = g(:,i) + alpha(n) * gamma(2) * regmatrix * s0(:,i);
            end
            
            % -------------------------------------------------------------
            % Gauss-Newton
            ds = zeros(size(s0), 'double');
            i  = 1;
            if optim(1)
                ds(:,i) = (H + alpha(n) * gamma(1) * regmatrix)\g(:,i);
                i = i + 1;
            end
            if optim(2)
                ds(:,i) = (H + alpha(n) * gamma(2) * regmatrix)\g(:,i);
            end
            clear g H
            ds = reshape(ds, [nbasis sum(optim)]);
            
            % -------------------------------------------------------------
            % Parts for log-likelihood (prior)
            llp_part1 = 0;
            llp_part2 = 0;
            i = 1;
            if optim(1)
                Lds = alpha(n) * gamma(1) * regmatrix * reshape(ds(:,:,:,i), [], 1);
                llp_part1 = llp_part1 + s0(:,i)' * Lds;
                llp_part2 = llp_part2 + Lds' * Lds;
                i = i + 1;
            end
            if optim(2)
                Lds = alpha(n) * gamma(2) * regmatrix * reshape(ds(:,:,:,i), [], 1);
                llp_part1 = llp_part1 + s0(:,i)' * Lds;
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
            
    % ---------------------------------------------------------------------
    % Line-Search
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
        if (llm+llp) > (llm0+llp0)
            ok = true;
            break
        else
            armijo = armijo/2;
        end
        
    end % < loop ls
    
    % ---------------------------------------------------------------------
    % Save
    if ok
        if all(optim)
            ds = ds(:,:,:,1) + 1i * ds(:,:,:,2);
        elseif optim(2)
            ds = 1i * ds;
        end
        s(:,:,:,n) = s(:,:,:,n) - armijo * ds;
    else
        llm = llm0;
        llp = llp0;
    end
        
end % < loop n

end % < function multicoil_sensitivity