function [s,llm,llp,ok,ls] = multicoil_sensitivity(varargin)
% Maximum a posteriori sensitivity profiles given a set of observed coil 
% images, a mean image and a noise precision (= inverse covariance) 
% matrix.
%
% FORMAT [s,...] = multicoil_sensitivity(rho, x, s, ...)
%
% REQUIRED
% --------
% rho   - (File)Array [Nx Ny Nz  1 (2)] - (Complex) mean image
% x     - (File)Array [Nx Ny Nz Nc (2)] - (Complex) coil images
% s     - (File)Array [Nx Ny Nz Nc (2)] - (Complex) log-sensitivity maps
%
% Nc = number of coils
% Images can either be complex or have two real components that are then 
% assumed to be the real and imaginary parts.
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
% CentreFields  - Logical          - Enforce zeros-centered fields [false]
% SensOptim     - [Mag Phase]      - Optimise magnitude/phase      [true true]
% Parallel      - Scalar | Logical - Number of parallel workers    [false]
%
% OUTPUT
% ------
% s   - Updated (complex) log-sensitivity maps
% llm - Matching part of the log-likelihood (all coils)
% llp - Prior part of the log-likelihood (all coils)
% ok  - True if a better value was found
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
% 
% Note that I am thinking of adding yet another representation, where
% log-sensitivity fields are directly encoded by their discrete cosine
% components. This might help to deal better with small autocalibration
% regions, where the finite element approximation used in the
% regularization matrix becomes too poor.
%__________________________________________________________________________

% -------------------------------------------------------------------------
% Helper functions to check input arguments
function ok = isarray(X)
    ok = isnumeric(X) || isa(X, 'file_array');
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
p.addParameter('Encoding',      'frequency', @ischar);
p.addParameter('NbBasis',       [10 10 10],  @(X) isnumeric(X) && numel(X) <= 3);
p.addParameter('RegMatrix',     [],          @(X) isnumeric(X));
p.parse(varargin{:});
rho   = p.Results.MeanImage;
x     = p.Results.CoilImages;
s     = p.Results.SensMaps;
all_n = p.Results.Index;
A     = p.Results.Precision;
prm   = p.Results.RegStructure;
alpha = p.Results.RegCoilFactor;
gamma = p.Results.RegCompFactor;
bnd   = p.Results.RegBoundary;
vs    = p.Results.VoxelSize;
llp   = p.Results.LLPrior;
optim = p.Results.SensOptim;
Nw    = p.Results.Parallel;
encoding  = p.Results.Encoding;
nbasis    = p.Results.NbBasis;
regmatrix = p.Results.RegMatrix;

% -------------------------------------------------------------------------
% Post-process input
N   = size(x,4);
lat = [size(x,1) size(x,2) size(x,3)];

% Coils to process: default = all + ensure row-vector
if isempty(all_n)
    all_n = 1:N;
end
all_n = all_n(:).';
% Precision: default = identity
if numel(A) == 1
    A = A * eye(N);
end
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
% Nb basis: ensure row vector + complete
nbasis = padarray(nbasis(:)', [0 max(0,3-numel(nbasis))], 'replicate', 'post');
% Voxel size: ensure row vector + complete
vs = padarray(vs(:)', [0 max(0,3-numel(vs))], 'replicate', 'post');
% Create regularisation matrix
if isempty(regmatrix)
    regmatrix = spm_bias_lib('regulariser', prm, lat, nbasis, vs);
end
% Create basis functions (B)
[B1,B2,B3] = spm_bias_lib('dcbasis', lat, nbasis);
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
        xz = loadarray(x(:,:,z,:), @single);
        xz = reshape(xz, [], N);

        % -----------------------------------------------------------------
        % Load one slice of the mean
        rhoz = loadarray(rho(:,:,z,:), @single);
        rhoz = reshape(rhoz, [], 1);

        % -----------------------------------------------------------------
        % Load one slice of the delta sensitivity
        dsz   = loadarray(ds(:,:,z,:), @single);
        dsz   = reshape(dsz, [], size(ds,4));
        if all(optim)
            dsz = dsz(:,1) + 1i * dsz(:,2);
        elseif optim(2)
            dsz = 1i * dsz(:,2);
        end
        
        % -----------------------------------------------------------------
        % Load one slice of the complete sensitivity dataset + correct
        sz      = loadarray(s(:,:,z,:), @single);
        sz      = reshape(sz, [], N);
        sz(:,n) = sz(:,n) - dsz;
        dsz     = [];
        sz      = single(exp(double(sz)));
        rhoz    = bsxfun(@times, rhoz, sz);
        sz      = []; % clear

        % -----------------------------------------------------------------
        % Compute log-likelihood
        tmp = (rhoz - xz) * A;
        llm = llm - 0.5 * sum(double(real(dot(tmp, rhoz - xz, 2))));
        
        rhoz = [];
        xz   = [];
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
        xz = loadarray(x(:,:,z,:), @double);
        xz = reshape(xz, [], N);

        % -----------------------------------------------------------------
        % Load one slice of the mean
        rhoz = loadarray(rho(:,:,z,:), @double);
        rhoz = reshape(rhoz, [], 1);

        % -----------------------------------------------------------------
        % Load one slice of the complete sensitivity dataset + correct
        sz   = loadarray(s(:,:,z,:), @double);
        sz   = reshape(sz, [], N);
        sz   = single(exp(double(sz)));
        rhoz = bsxfun(@times, rhoz, sz);
        sz   = []; % clear
        
        % -----------------------------------------------------------------
        % Compute gradient
        
        tmp = (rhoz - xz) * A;
        
        llm = llm - 0.5 * sum(real(dot(tmp, rhoz - xz, 2)));
        
        tmp = rhoz(:,n) .* conj(tmp(:,n));
        
        gz = zeros([size(tmp,1) sum(optim)], 'like', real(tmp(1)));
        i  = 1;
        if optim(1) % If optimise sensitivity magnitude
            gz(:,i) = real(tmp);
            i = i+1;
        end
        if optim(2) % If optimise sensitivity phase
            gz(:,i) = -imag(tmp);
        end
        
        Hz = A(n,n) * real(conj(rhoz(:,n)) .* rhoz(:,n));
        
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
        
        tmp  = []; % clear
        xz   = []; % clear
        rhoz = []; % clear
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
    for ls=1:6
        
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