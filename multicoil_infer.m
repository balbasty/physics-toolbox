function [s,rho,A,ll,llm,llp] = multicoil_infer(varargin)
% Compute mode estimates (ML, MAP) of the parameters of a probabilistic 
% model of complex multicoil MR images.
%
% FORMAT [sens,mean,prec,ll] = multicoil_infer(coils, ...)
%
% REQUIRED
% --------
% coils - (File)Array [Nx Ny Nz Nc] - Complex coil images
%
% KEYWORDS
% --------
% Precision     - Noise precision matrix                    [NaN=estimate]
% RegCoilFactor - Regularisation factor per coil            [1/Nc]
% RegCoilComp   - Regularisation factor per Re/Im component [1E6]
% RegDecFactor  - Start with more regularisation            [3] (log10)
% VoxelSize     - Voxel size                                [1]
% SensOptim     - Optimize real and/or imaginary parts      [true true]
% CovOptim      - Optimize noise covariance                 [false]
% Tolerance     - Convergence threshold                     [1E-3]
% IterMax       - Total maximum number of iterations        [15]
% IterMin       - Total minimum number of iterations        [1]
% Verbose       - (-1=quiet,0=moderate,1=verbose,2=plot)    [0]
%
% ADVANCED
% --------
% SamplingMask  - Mask of the k-space sampling scheme       [[]=fully sampled]
% SensMaps      - Initial complex sensitivity profiles      (File)Array [Nx Ny Nz Nc]
% MeanImage     - Initial complex mean image                (File)Array [Nx Ny Nz]
% RegStructure  - Regularisation Structure (abs memb bend)  [0 0 1]
% RegBoundary   - Boundary conditions for sensitivities     ['neumann']
% Parallel      - Number of parallel workers                [0]
% LLCond        - Previous conditional log-likelihood       [NaN=compute]
% LLPrior       - Previous prior log-likelihood             [NaN=compute]
% LLPrev        - Log-likelihood of previous iterations     []
%
% OUTPUT
% ------
% sens - (Log)-Sensitivity maps - (File)Array [Nx Ny Nz Nc]
% mean - Mean image             - (File)Array [Nx Ny Nz Nc]
% prec - Noise precision        -       Array [Nc Nc]
% ll   - Log-likelihood
%
% An output FileArray can be provided by using `MeanImage` as an input. If  
% not provided, the output volume will have the same format as the input  
% coil volume.
%
% Nc = number of coils
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
        okk = numel(T) > 7 || strcmpi(T(1:7),'complex');
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
Nc = size(varargin{1},4);
p  = inputParser;
p.FunctionName = 'multicoil_infer';
p.addRequired('CoilImages',                  @isarray);
p.addParameter('SensMaps',      [],          @isarray);
p.addParameter('MeanImage',     [],          @isarray);
p.addParameter('Precision',     NaN,         @isnumeric);
p.addParameter('RegStructure',  [0 0 1],     @(X) isnumeric(X) && numel(X) == 3);
p.addParameter('RegCoilFactor', 1/Nc,        @isnumeric);
p.addParameter('RegCompFactor', 1E6,         @(X) isnumeric(X) && numel(X) <= 2);
p.addParameter('RegDecFactor',  3,           @(X) isnumeric(X) && isscalar(X));
p.addParameter('RegBoundary',   1,           @isboundary);
p.addParameter('VoxelSize',     [1 1 1],     @(X) isnumeric(X) && numel(X) <= 3);
p.addParameter('SensOptim',     [true true], @(X) (isnumeric(X) || islogical(X)) && numel(X) == 2);
p.addParameter('CovOptim',      false,       @(X) (isnumeric(X) || islogical(X)) && isscalar(X));
p.addParameter('Parallel',      0,           @(X) (isnumeric(X) || islogical(X)) && isscalar(X));
p.addParameter('Tolerance',     1E-3,        @(X) isnumeric(X) && isscalar(X));
p.addParameter('IterMax',       15,          @(X) isnumeric(X) && isscalar(X));
p.addParameter('IterMin',       1,           @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLCond',        NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLPrior',       NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLPrev',        [],          @(X) isnumeric(X));
p.addParameter('Verbose',       0,           @(X) (isnumeric(X) || islogical(X)) && isscalar(X));
p.addParameter('SamplingMask',  [],          @isarray);
p.parse(varargin{:});
rho           = p.Results.MeanImage;
x             = p.Results.CoilImages;
s             = p.Results.SensMaps;
A             = p.Results.Precision;
reg           = p.Results.RegStructure;
alpha         = p.Results.RegCoilFactor;
gamma         = p.Results.RegCompFactor;
decfactor     = p.Results.RegDecFactor;
bnd           = p.Results.RegBoundary;
vs            = p.Results.VoxelSize;
optim_cov     = p.Results.CovOptim;
optim         = p.Results.SensOptim;
tol           = p.Results.Tolerance;
itermax       = p.Results.IterMax;
itermin       = p.Results.IterMin;
verbose       = p.Results.Verbose;
llm           = p.Results.LLCond;
llp           = p.Results.LLPrior;
ll            = p.Results.LLPrev;
Nw            = p.Results.Parallel;
mask          = p.Results.SamplingMask;

if optim_cov
    warning('Covariance inference is deactivated for now.')
end

% -------------------------------------------------------------------------
% Time execution
% -------------------------------------------------------------------------
if verbose > -1
    fprintf('Processing started\n');
    start = tic;
end

% -------------------------------------------------------------------------
% Post-process input
% -------------------------------------------------------------------------
% Store a few useful values
Nx   = size(x,1);
Ny   = size(x,2);
Nz   = size(x,3);
Nc   = size(x,4);
Nvox = Nx*Ny*Nz;
if ~isempty(mask)
    propmask = sum(mask(:))/numel(mask);
else
    propmask = 1;
end
% Precision: default = estimate from magnitude
if numel(A) == 1
    A = A * eye(Nc);
end
% Pad voxel size
vs = padarray(vs(:)', [0 max(0,3-numel(vs))], 'replicate', 'post');
% Reg components: pad reg structure
gamma = padarray(gamma(:)', [0 max(0,2-numel(gamma))], 'replicate', 'post');
% Reg coil factor: ensure zero sum -> propagate their sum to reg components
alpha = padarray(alpha(:), [max(0,Nc-numel(alpha)) 0], 'replicate', 'post');
gamma = gamma * sum(alpha);
alpha = alpha/sum(alpha);
% Decreasing regularisation
decfactor0 = decfactor;
decfactor  = linspace(decfactor0, 0, itermax)';
decfactor  = 10.^decfactor;
gammafinal = gamma;
% Allocate mean
init_rho = isempty(rho);
if init_rho
    rho = zeros(Nx, Ny, Nz, 'like', x(1));
end
init_s = isempty(s);
if init_s
    s = zeros(Nx, Ny, Nz, Nc, 'like', x(1));
    if isnan(llp)
        llp = 0;
    end
end
% Parallel processing
if isnan(Nw),     Nw = 0; end
if ~isfinite(Nw), Nw = parcluster('local'); Nw = Nw.NumWorkers; end


% =========================================================================
%
%                           INITIAL ESTIMATES
%
% =========================================================================

if verbose > -1
    fprintf('Initialisation\n');
end

% -------------------------------------------------------------------------
% Estimate covariance precision
% -------------------------------------------------------------------------
if isempty(A) || any(any(isnan(A)))
    [~,A] = multicoil_init_cov(x);
end

% -------------------------------------------------------------------------
% Compute log-determinant of precision matrix
% -------------------------------------------------------------------------
C   = inv(A);
ldC = spm_matcomp('LogDet', C);

% -------------------------------------------------------------------------
% Initial estimate of the mean
% -------------------------------------------------------------------------
if verbose > 1
    multicoil_plot_mean(rho, C, ll, vs);
end
if init_rho
    if verbose > 0, fprintf('> Update mean: '); end
    [rho,llm,~,ok,ls] = multicoil_mean_map(...
        x, s, rho,          ...
        'Precision',    A,  ...
        'VoxelSize',    vs, ...
        'SamplingMask', mask);
    llm = llm - propmask*Nvox*ldC;
    if verbose > 0
        if ok, fprintf(' :D (%d)\n', ls);
        else,  fprintf(' :(\n');
        end
    end
    if verbose > 1
        multicoil_plot_mean(rho, C, ll, vs);
    end
end

% -------------------------------------------------------------------------
% Initial estimate of the sensitivity phase
% -------------------------------------------------------------------------
% Phase is periodic, and phase wraps should not be penalised. However,
% periodic domains are not handled by spm_field, which deals with
% regularization. To circumvent this issue, I try to never wrap the phase.
% A common problem, though, is that the same target value of pi/-pi might 
% be converged towards from two sides (negative and positive), which
% sometimes leads to bad local minima. The solution I found is to
% initialise the sensitivity phase using the Von Mises mean of pointwise
% differences between the mean and the coil image.
% See: * Bishop's PRML - chapter 2.3.8
%      * https://en.wikipedia.org/wiki/Von_Mises_distribution
% Let's do a few iterations of those to centre the mean image.
if verbose > 1
    multicoil_plot_mean(rho, C, ll, vs);
end
if init_s
    llp = 0;
    for i=1:3
        if verbose > 0
            fprintf('> Update sensitivity\n');
        end
        if optim(1)
            s   = multicoil_init_magnitude(rho, x, s);
        end
        if optim(2)
            s   = multicoil_init_phase(rho, x, s);
        end

        if init_rho
            if verbose > 0, fprintf('> Update mean: '); end
            [rho,llm,~,ok,ls] = multicoil_mean_map(...
                x, s, rho,          ...
                'Precision',    A,  ...
                'VoxelSize',    vs, ...
                'SamplingMask', mask);
            llm = llm - propmask*Nvox*ldC;
            if verbose > 0
                if ok, fprintf(' :D (%d)\n', ls);
                else,  fprintf(' :(\n');
                end
            end
            if verbose > 1
                multicoil_plot_mean(rho, C, ll, vs);
            end
        end
    end
end

% -------------------------------------------------------------------------
% Log-likelihood
% -------------------------------------------------------------------------

if isnan(llp)
    % > Initial log-likelihood (prior term)
    llp = multicoil_ll_prior(s, reg, gammafinal.*decfactor(1), alpha, bnd, optim, vs);
end
if isnan(llm)
    % > Initial log-likelihood (cond term)
    llm = multicoil_ll_cond(x,s,rho,A,mask) - propmask*Nvox*ldC;
end
ll = [ll (sum(llm)+sum(llp))];

if verbose > 1
    multicoil_plot_mean(rho, C, ll, vs);
end

% -------------------------------------------------------------------------
% Loop
% -------------------------------------------------------------------------
for it=1:itermax
    
    gamma = gammafinal .* decfactor(it);
    if it > 1
        llp = llp * decfactor(it) / decfactor(it-1);
    end
    
    if verbose > -1
        if verbose > 0
            fprintf([repmat('-', [1 78]) '\n']);
        end
        fprintf('Iteration %d\n', it);
        if decfactor0
            fprintf('* reg magnitude = %7.1e\n', gamma(1));
            fprintf('* reg phase     = %7.1e\n', gamma(2));
        end
    end
    
    % ---------------------------------------------------------------------
    % Update mean/covariance (closed-form)
    % ---------------------------------------------------------------------
%     if optim_cov
%         if verbose > 0
%             fprintf('> Update Covariance\n');
%         end
%         % rho     = multicoil_mean_ml(x, s, A, rho, optim);
%         % rho = multicoil_mean_map(x, s, A, rho, [0 1E-5 0], vs, mask);
%         [C,A]   = multicoil_cov(rho, x, s);
%         ldC = spm_matcomp('LogDet', C);
%         if verbose > 1
%             multicoil_plot_mean(rho, C, ll, vs);
%         end
%     end
    
    % ---------------------------------------------------------------------
    % Coil-wise sensitivity update (Gauss-Newton)
    % ---------------------------------------------------------------------
    if verbose > 0, fprintf('> Update Sensitivity:'); end
    if isdiag(A) && verbose < 2
        % If the precision matrix is diagonal, computation of the
        % log-likelihood and derivatives in implemented in a slightly
        % faster way (each sensitivity field only depends on one coil).
        %
        % * This means that we can distribute the processing of each coil!
        %
        % * This also means that we need to sum the log-likelihood of each 
        %   coil to get the complete model log-likelihood.
        llp0 = llp;
        if numel(llp0) == 1
            llp0 = llp0 * ones(1,Nc);
        end
        llp  = zeros(1,Nc);
        llm  = zeros(1,Nc);
        Ad   = diag(A);
        parfor(n=1:Nc, Nw)

            [s(:,:,:,n),llm(n),llp(n),ok,ls] = multicoil_sensitivity(...
                rho, x(:,:,:,n), s(:,:,:,n), ...
                'Precision',     Ad(n),     ...
                'RegStructure',  reg,       ...
                'RegCoilFactor', alpha(n),  ...
                'RegCompFactor', gamma,     ...
                'RegBoundary',   bnd,       ...
                'VoxelSize',     vs,        ...
                'SensOptim',     optim,     ...
                'LLPrior',       llp0(n),   ...
                'SamplingMask',  mask);
            if verbose > 0
                if ok, fprintf(' [%2d :D (%d)]', n, ls);
                else,  fprintf(' [%2d :(    ]', n);
                end
            end

        end
        llm = sum(llm);
    else
        if isdiag(A)
            % We do not parallelise, but still need to sum individual 
            % coil-wise log-likelihoods
            llm = 0;
        end
        for n=1:Nc

            if verbose > 0
                if verbose > 1
                    multicoil_plot_fit(n, x, s, rho, mask, vs);
                end
                if mod(n,4) == 1, fprintf('\n  | '); end
                fprintf('%2d', n);
            end
            [s,llm1,llp,ok,ls] = multicoil_sensitivity(...
                rho, x, s,              ...
                'Index',         n,     ...
                'Precision',     A,     ...
                'RegStructure',  reg,   ...
                'RegCoilFactor', alpha, ...
                'RegCompFactor', gamma, ...
                'RegBoundary',   bnd,   ...
                'VoxelSize',     vs,    ...
                'SensOptim',     optim, ...
                'LLPrior',       sum(llp),   ...
                'SamplingMask',  mask);
            if isdiag(A)
                llm = llm + llm1;
            else
                llm = llm1;
            end
            if verbose > 0
                if ok, fprintf(' :D (%d) | ', ls);
                else,  fprintf(' :(     | ');
                end
                if verbose > 1
                    multicoil_plot_fit(n, x, s, rho, mask, vs);
                end
            end

        end
    end
    if verbose > 0, fprintf('\n'); end
    
    % Add normalisation term
    llm = llm - propmask*Nvox*ldC;
    
    % ---------------------------------------------------------------------
    % Center sensitivity fields
    % ---------------------------------------------------------------------
    % We know that, at the optimum, sum{alpha_n*s_n} = 0
    % To converge faster, and to avoid bias towards the initial mean
    % estimate, we enforce this condition at each iteration.
    meansen = zeros(Nx, Ny, Nz, 'like', s(1));
    for n=1:Nc
        meansen = meansen + alpha(n) * single(s(:,:,:,n));
    end
    for n=1:Nc
        s(:,:,:,n) = s(:,:,:,n) - meansen;
    end
    clear meansen
    % Zero-centering the sensitivity changes the conditional and prior
    % terms of the og-likelihood.
    % However, in practice, the impact is small, and updating the 
    % log-likelihood is costly. Therefore, we keep it as is.
    % This means that the log-likelihood changes slightly and might drop
    % during the first iterations.
    
    % ---------------------------------------------------------------------
    % Update mean (Gauss-Newton)
    % ---------------------------------------------------------------------
    if verbose > 0, fprintf('> Update mean:\n '); end
    for i=1:4
        [rho,llm,~,ok,ls] = multicoil_mean_map(...
            x, s, rho,          ...
            'Precision',    A,  ...
            'VoxelSize',    vs, ...
            'SamplingMask', mask);
        if verbose > 0
            if ok, fprintf(' :D (%d) |', ls);
            else,  fprintf(' :(     |');
            end
        end
        if verbose > 1
            multicoil_plot_mean(rho, C, ll, vs);
        end
   end
   if verbose > 0, fprintf('\n'); end
   llm = llm - propmask*Nvox*ldC;
    
    % ---------------------------------------------------------------------
    % Update log-likelihood
    % ---------------------------------------------------------------------
    ll = [ll (sum(llm)+sum(llp))];
    if verbose > 1
        multicoil_plot_mean(rho, C, ll, vs);
    end
    
    % ---------------------------------------------------------------------
    % Check gain
    % ---------------------------------------------------------------------
    if it > 1
        llmin = min(ll, [], 'omitnan');
        llmax = max(ll, [], 'omitnan');
        gain  = (ll(end) - ll(end-1))/(llmax-llmin);
        if verbose > 0
            switch sign(gain)
                case  1,   sgn = '+';
                case -1,   sgn = '-';
                case  0,   sgn = '=';
                otherwise, sgn = '';
            end
            fprintf('> Gain: %20.10g (%s)\n', gain, sgn);
        end
        if it >= itermax
            if verbose > -1
                fprintf('Reached maximum number of iterations (%d)\n', it);
            end
            break
        end
        if (it >= itermin) && abs(gain) < tol
            if verbose > -1
                fprintf('Converged (%f < %f)\n', abs(gain), tol);
            end
            break
        end
    end
    
end


% -------------------------------------------------------------------------
% Time execution
% -------------------------------------------------------------------------
if verbose > -1
    stop = toc(start);
    fprintf('Processing finished: in %s\n', sec2ydhms(stop));
end

end % < function multicoil_infer