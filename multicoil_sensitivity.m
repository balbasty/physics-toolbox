function [s,llm,llp,ok,ls] = multicoil_sensitivity(varargin)
% Maximum a posteriori sensitivity profiles given a set of observed coil 
% images, a mean image and a noise precision (= inverse covariance) 
% matrix.
%
% FORMAT [s,...] = multicoil_sensitivity(rho, x, s, ...)
%
% MANDATORY
% ---------
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
    if size(X,5) == 2
        ok = false;
        return
    end
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
p.addParameter('CentreFields',  false,       @(X) (isnumeric(X) || islogical(X)) && isscalar(X));
p.addParameter('SensOptim',     [true true], @(X) (isnumeric(X) || islogical(X)) && numel(X) == 2);
p.addParameter('Parallel',      0,           @(X) (isnumeric(X) || islogical(X)) && isscalar(X));
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
centre_fields = p.Results.CentreFields;
Nw    = p.Results.Parallel;

% -------------------------------------------------------------------------
% Post-process input
N = size(x,4);

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
% Voxel size: ensure row vector
vs = double(vs(:)');
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
if Nw > 0
    warning('Parallel processing not implemented. Running sequential instead.')
    Nw = 0;
end
% Centre fields: ensure logical
centre_fields = logical(centre_fields);

% -------------------------------------------------------------------------
% Boundary condition (usually Neumann = null derivative)
spm_field('boundary', bnd); 

% -------------------------------------------------------------------------
% For each coil
for n=all_n
    
    % ---------------------------------------------------------------------
    % Compute log-likelihood (prior)
    if isnan(llp)
        llp = multicoil_ll_prior(s, prm, gamma, alpha, bnd, optim, vs);
    end
    
    % ---------------------------------------------------------------------
    % Prepare weights: beta(n,m) = [n == m] - alpha(n)
    beta    = repmat(-alpha(n), [1 N]);
    beta(n) = 1 + beta(n);
    
    % ---------------------------------------------------------------------
    % Allocate conjugate gradient and Hessian
    g = zeros(size(x,1), size(x,2), size(x,3), sum(optim), 'single');
    H = zeros(size(x,1), size(x,2), size(x,3), 1, 'single');
    if isa(A, 'gpuArray')
        g = gpuArray(g);
        H = gpuArray(H);
    end
    llm = 0;
    
    % ---------------------------------------------------------------------
    % Compute gradient slice-wise to save memory
    % parfor(z=1:size(rho, 3), Nw) % < Uncomment for parallel processing
    for z=1:size(rho,3)          % < Uncomment for sequential processing

        gz = g(:,:,z,:);
        Hz = H(:,:,z);
        
        % -----------------------------------------------------------------
        % Load one slice of the complete coil dataset
        if size(x, 5) == 2
            % Two real components
            x1 = reshape(single(x(:,:,z,:,:)), [], size(x,4), 2);
            x1 = x1(:,:,1) + 1i*x1(:,:,2);
        else
            % One complex volume
            x1 = reshape(single(x(:,:,z,:)), [], size(x,4));
        end
        if isa(A, 'gpuArray')
            x1 = gpuArray(x1);
        end

        % -----------------------------------------------------------------
        % Load one slice of the mean
        if size(rho, 5) == 2
            % Two real components
            rho1 = reshape(single(rho(:,:,z,:,:)), [], 1, 2);
            rho1 = rho1(:,:,1) + 1i*rho1(:,:,2);
        else
            % One complex volume
            rho1 = reshape(single(rho(:,:,z,:,:)), [], 1);
        end
        if isa(A, 'gpuArray')
            rho1 = gpuArray(rho1);
        end

        % -----------------------------------------------------------------
        % Load one slice of the complete sensitivity dataset + correct
        if size(s, 5) == 2
            % Two real components
            s1 = reshape(double(s(:,:,z,:,:)), [], size(x,4), 2);
            s1 = single(exp(s1(:,:,1) + 1i*s1(:,:,2)));
        else
            % One complex volume
            s1 = reshape(single(exp(double(s(:,:,z,:,:)))), [], size(x,4));
        end
        if isa(A, 'gpuArray')
            s1 = gpuArray(s1);
        end
        rho1 = bsxfun(@times, rho1, s1);
        s1   = []; % clear
        
        % -----------------------------------------------------------------
        % Compute gradient
        
        tmp = (rho1 - x1) * A;
        
        llm = llm - 0.5 * sum(double(real(dot(tmp, rho1 - x1, 2))));
        
        if centre_fields
            rho1 = bsxfun(@times, rho1, beta);
            tmp  = dot(tmp, rho1, 2);
        else
            tmp = rho1(:,n) .* conj(tmp(:,n));
        end
        
        i = 1;
        if optim(1) % If optimise sensitivity magnitude
            gz(:,:,1,i) = gz(:,:,1,i) + reshape( real(tmp), size(gz,1), size(gz,2));
            i = i+1;
        end
        if optim(2) % If optimise sensitivity phase
            gz(:,:,1,i) = gz(:,:,1,i) + reshape(-imag(tmp), size(gz,1), size(gz,2));
        end
        g(:,:,z,:)  = gz;
        gz          = []; % clear
        
        if centre_fields
            tmp = real(dot(rho1, rho1 * A, 2));
        else
            tmp = A(n,n) * real(conj(rho1(:,n)) .* rho1(:,n));
        end
        H(:,:,z) = reshape(tmp, size(Hz,1), size(Hz,2));
        
        tmp  = []; % clear
        x1   = []; % clear
        rho1 = []; % clear
    end
    
    % ---------------------------------------------------------------------
    % Gather gradient & Hessian (if on GPU)
    g = gather(g);
    H = gather(H);
    
    % ---------------------------------------------------------------------
    % Gradient: add prior term
    if size(s, 5) == 2
        % Two real components
        if all(optim)
            s0 = single(s(:,:,:,n,:));
        elseif optim(1)
            s0 = single(s(:,:,:,n,1));
        elseif optim(2)
            s0 = single(s(:,:,:,n,2));
        end
    else
        % One complex volume
        if all(optim)
            s0 = single(s(:,:,:,n));
            s0 = cat(5, real(s0), imag(s0));
        elseif optim(1)
            s0 = real(single(s(:,:,:,n)));
        elseif optim(2)
            s0 = imag(single(s(:,:,:,n)));
        end
    end
    s0 = reshape(s0, [size(s0,1) size(s0,2) size(s0,3) size(s0,5)]);
    g  = g + spm_field('vel2mom', s0, [vs alpha(n) * prm], gamma(optim));
    clear s0

    % ---------------------------------------------------------------------
    % Gauss-Newton
    regH = alpha(n) * prm;
    if centre_fields
        regH = regH * beta(n);
    end
    ds = [];
    i  = 1;
    if optim(1)
        ds = cat(4,ds,spm_field(H, g(:,:,:,i), [vs regH 2 2], gamma(1)));
        i = i+1;
    end
    if optim(2)
        ds = cat(4,ds,spm_field(H, g(:,:,:,i), [vs regH 2 2], gamma(2)));
    end
    clear g H
    
    % ---------------------------------------------------------------------
    % Parts for log-likelihood (prior)
    Lds = spm_field('vel2mom', ds, [vs alpha(n) * prm], gamma(optim));
    if centre_fields
        sums = 0;
        sumb = sum(alpha .* beta(:)'.^2);
        for m=1:N
            if size(s, 5) == 2
                % Two real components
                if all(optim)
                    s1 = single(s(:,:,:,m,:));
                elseif optim(1)
                    s1 = single(s(:,:,:,m,1));
                elseif optim(2)
                    s1 = single(s(:,:,:,m,2));
                end
            else
                % One complex volume
                if all(optim)
                    s1 = single(s(:,:,:,m));
                    s1 = cat(4, real(s1), imag(s1));
                elseif optim(1)
                    s1 = real(single(s(:,:,:,m)));
                elseif optim(2)
                    s1 = imag(single(s(:,:,:,m)));
                end
            end
            sums = sums + alpha(m) * beta(m) * s1;
            clear s1
        end
        % part1 = (sum alpha*beta) * (ds)'L(ds)
        % part2 = -2 * (sum alpha*beta*s)'L(ds)
        part1 = sumb * double(reshape(ds, 1, [])) * double(reshape(Lds, [], 1));
        part2 = -2 * double(reshape(Lds, 1, [])) * double(reshape(sums, [], 1));
        clear Lds sums sumb
    else
        part = alpha(n) * double(reshape(ds, 1, [])) * double(reshape(Lds, [], 1));
        clear Lds
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
        if centre_fields
            llp = -0.5 * (armijo^2 * part1 + armijo * part2);
        else
            llp = -0.5 * armijo^2 * part;
        end
        llp = llp0 + llp;
        
        % -----------------------------------------------------------------
        % Compute log-likelihood (conditional)
        llm = computellm(A, rho, x, s, armijo * ds, n, beta, optim);
        
        % -----------------------------------------------------------------
        % Check progress
        if (llm+llp) > (llm0+llp0)
            ok = true;
            break
        else
            armijo = armijo/2;
        end
        
    end
    
    % ---------------------------------------------------------------------
    % Save
    if ok
        ds = reshape(ds, [size(ds,1) size(ds,2) size(ds,3) 1 size(ds,4)]);
        if size(s,5) == 1
            if all(optim)
                ds = ds(:,:,:,:,1) + 1i * ds(:,:,:,:,2);
            elseif optim(2)
                ds = 1i * ds;
            end
            msk = true;
        else
            msk = optim;
        end
        if centre_fields
            for m=1:size(x,4)
                s(:,:,:,m,msk) = s(:,:,:,m,msk) - beta(m) * armijo * ds;
            end
        else
            s(:,:,:,n,msk) = s(:,:,:,n,msk) - armijo * ds;
        end
    else
        llm = llm0;
        llp = llp0;
    end
        
end

end

function llm = computellm(A, rho, x, s, ds, n, beta, optim)

centre_fields = size(ds,4) == size(s,4);
if centre_fields
    beta = reshape(beta, [1 1 size(x,4)]);
end

llm = 0;
for z=1:size(rho, 3) 

    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    if size(x, 5) == 2
        % Two real components
        x1 = reshape(single(x(:,:,z,:,:)), [], size(x,4), 2);
        x1 = x1(:,:,1) + 1i*x1(:,:,2);
    else
        % One complex volume
        x1 = reshape(single(x(:,:,z,:)), [], size(x,4));
    end
    if isa(A, 'gpuArray')
        x1 = gpuArray(x1);
    end

    % ---------------------------------------------------------------------
    % Load one slice of the mean
    if size(rho, 5) == 2
        % Two real components
        rho1 = reshape(single(rho(:,:,z,:,:)), [], 1, 2);
        rho1 = rho1(:,:,1) + 1i*rho1(:,:,2);
    else
        % One complex volume
        rho1 = reshape(single(rho(:,:,z,:,:)), [], 1);
    end
    if isa(A, 'gpuArray')
        rho1 = gpuArray(rho1);
    end

    % ---------------------------------------------------------------------
    % Load one slice of the search direction + make it complex
    if all(optim)
        ds1 = reshape(single(ds(:,:,z,:)), [], 2);
        ds1 = ds1(:,1) + 1i * ds1(:,2);
    elseif optim(1)
        ds1 = reshape(single(ds(:,:,z)), [], 1);
    elseif optim(2)
        ds1 = reshape(single(ds(:,:,z)) * 1i, [], 1);
    end
    if isa(A, 'gpuArray')
        ds1 = gpuArray(ds1);
    end
    
    % ---------------------------------------------------------------------
    % Load one slice of the complete sensitivity dataset + correct
    if size(s, 5) == 2
        % Two real components
        s1 = reshape(single(s(:,:,z,:,:)), [], size(x,4), 2);
        s1 = s1(:,:,1) + 1i*s1(:,:,2);
    else
        % One complex volume
        s1 = reshape(single(s(:,:,z,:,:)), [], size(x,4));
    end
    if isa(A, 'gpuArray')
        s1 = gpuArray(s1);
    end
    if centre_fields
        s1 = s1 - bsxfun(@times, beta, ds1);
    else
        s1(:,n) = s1(:,n) - ds1;
    end
    s1   = single(exp(double(s1)));
    rho1 = bsxfun(@times, rho1, s1);
    s1   = []; % clear

    % ---------------------------------------------------------------------
    % Compute log-likelihood
    llm = llm - 0.5 * sum(double(real(dot((rho1 - x1) * A, rho1 - x1, 2))));

end

end % < function computellm